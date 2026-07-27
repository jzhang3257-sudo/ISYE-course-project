# model3_asymmetric.jl — ISyE 524 Course Project
# Model 3: Asymmetric Distance Sensitivity Analysis
#
# Replaces symmetric Haversine distances with perturbed asymmetric distances
# to simulate the effect of one-way streets and road geometry. Runs Model 1
# with each perturbed matrix and compares the optimal plan to the baseline.
#
# Key metrics:
#   (a) Plan divergence: Jaccard distance between active-arc sets
#   (b) Cost gap: % increase when evaluating the symmetric plan on asymmetric distances
#   (c) Integrality gap: does asymmetry affect the LP relaxation quality?
#
# Usage:
#   julia model3_asymmetric.jl              # full dataset
#   julia model3_asymmetric.jl subset_10    # 10-station subset
using CSV, DataFrames, JuMP, HiGHS, Random, Statistics

Random.seed!(524)

# ─── Parameters ───────────────────────────────────────────────────────────────
const SUBSET_NAME = length(ARGS) > 0 ? ARGS[1] : ""
const DATA_DIR = joinpath(@__DIR__, "..", "data", SUBSET_NAME)
const DELTAS = [0.05, 0.10, 0.20]  # perturbation magnitudes

# ─── Load data ───────────────────────────────────────────────────────────────
stations_df = CSV.read(joinpath(DATA_DIR, "stations.csv"), DataFrame)
distances_df = CSV.read(joinpath(DATA_DIR, "distances.csv"), DataFrame)

dist_sym = Matrix(distances_df[:, 2:end])  # symmetric baseline

surplus_idx = findall(stations_df.net_flow .> 0)
deficit_idx = findall(stations_df.net_flow .< 0)

S = length(surplus_idx)
D = length(deficit_idx)
supply = stations_df.net_flow[surplus_idx]
demand = -stations_df.net_flow[deficit_idx]

println("Model 3: Asymmetric Distance Sensitivity")
println("  Dataset: $(isempty(SUBSET_NAME) ? "full" : SUBSET_NAME)")
println("  Stations: $(S+D), Surplus: $S, Deficit: $D")
println("  Perturbation δ: $DELTAS")

# ─── Helper: solve Model 1 with a given distance matrix ──────────────────────
function solve_model1(dist_mat)
    m = Model(HiGHS.Optimizer)
    set_silent(m)

    @variable(m, y[1:S, 1:D] >= 0)

    total_sup = sum(supply)
    total_dem = sum(demand)
    shortfall = max(0, total_dem - total_sup)
    penalty = maximum(dist_mat) * 10

    if shortfall > 0
        @variable(m, y_dummy[1:D] >= 0)
        @objective(m, Min,
            sum(dist_mat[surplus_idx[i], deficit_idx[j]] * y[i, j] for i in 1:S, j in 1:D) +
            sum(penalty * y_dummy[j] for j in 1:D))
        @constraint(m, demand_con[j in 1:D],
            sum(y[i, j] for i in 1:S) + y_dummy[j] == demand[j])
    else
        @objective(m, Min,
            sum(dist_mat[surplus_idx[i], deficit_idx[j]] * y[i, j] for i in 1:S, j in 1:D))
        @constraint(m, demand_con[j in 1:D],
            sum(y[i, j] for i in 1:S) == demand[j])
    end

    @constraint(m, supply_con[i in 1:S],
        sum(y[i, j] for j in 1:D) <= supply[i])

    optimize!(m)

    y_val = value.(y)
    real_bikes = sum(y_val)
    real_obj = objective_value(m)
    if @isdefined(y_dummy)
        real_obj -= penalty * sum(value.(y_dummy))
    end
    nonzero = count(y_val[i, j] > 1e-6 for i in 1:S, j in 1:D)

    return (y_val = y_val, real_bikes = real_bikes, real_cost = real_obj,
            nonzero = nonzero, status = termination_status(m))
end

# ─── Baseline: symmetric solution ────────────────────────────────────────────
println("\n--- Baseline (symmetric) ---")
baseline = solve_model1(dist_sym)
println("  Cost: $(round(Int, baseline.real_cost)), Bikes: $(round(Int, baseline.real_bikes))")
println("  Nonzero arcs: $(baseline.nonzero)")

# ─── Test each perturbation level ─────────────────────────────────────────────
for δ in DELTAS
    println("\n--- δ = $δ ---")

    # Generate perturbed asymmetric distance matrix
    dist_asym = copy(dist_sym)
    perturbation = zeros(size(dist_sym))
    for i in 1:size(dist_sym, 1), j in 1:size(dist_sym, 2)
        if i != j
            # Nonnegative increases reflect that road distance is normally no
            # shorter than straight-line Haversine distance.
            eps_ij = δ * rand()  # 0 to δ
            dist_asym[i, j] = dist_sym[i, j] * (1.0 + eps_ij)
            perturbation[i, j] = eps_ij
        end
    end

    asymmetry = mean(abs.(dist_asym - dist_asym') ./ max.(dist_sym, 0.01))
    println("  Mean asymmetry: $(round(asymmetry*100; digits=1))%")

    # Solve with asymmetric distances
    asym_result = solve_model1(dist_asym)
    println("  Asym cost: $(round(Int, asym_result.real_cost)), Bikes: $(round(Int, asym_result.real_bikes))")
    println("  Asym nonzero arcs: $(asym_result.nonzero)")

    # (a) Plan divergence: which arcs differ?
    y_sym = baseline.y_val
    y_asym = asym_result.y_val
    sym_active = Set((i, j) for i in 1:S, j in 1:D if y_sym[i, j] > 1e-6)
    asym_active = Set((i, j) for i in 1:S, j in 1:D if y_asym[i, j] > 1e-6)

    # Use the symmetric difference divided by the union (Jaccard distance).
    # This gives a clear percentage between 0% and 100%.
    all_active = union(sym_active, asym_active)
    changed_arcs = symdiff(sym_active, asym_active)
    n_changed = length(changed_arcs)
    n_union = length(all_active)
    divergence = n_union > 0 ? n_changed / n_union * 100 : 0.0
    println("  Plan divergence: $(round(divergence; digits=1))% ($n_changed of $n_union union arcs differ)")

    # (b) Cost gap: evaluate symmetric plan on asymmetric distances
    sym_on_asym_cost = sum(dist_asym[surplus_idx[i], deficit_idx[j]] * y_sym[i, j]
                           for i in 1:S, j in 1:D)
    if asym_result.real_cost > 0
        cost_gap = (sym_on_asym_cost - asym_result.real_cost) / asym_result.real_cost * 100
    else
        cost_gap = 0.0
    end
    println("  Sym plan on asym distances: $(round(Int, sym_on_asym_cost))")
    println("  Cost gap: $(round(cost_gap; digits=1))%")

    # (c) Normalized L1 flow difference
    total_volume = sum(abs(y_sym[i, j] - y_asym[i, j]) for i in 1:S, j in 1:D)
    total_flow = sum(y_sym) + sum(y_asym)
    volume_change = total_flow > 0 ? total_volume / total_flow * 100 : 0.0
    println("  Normalized flow difference: $(round(volume_change; digits=1))%")
end

println("\nDone.")