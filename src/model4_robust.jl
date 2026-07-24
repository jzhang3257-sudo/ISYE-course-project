# model4_robust.jl — ISyE 524 Course Project
# Model 4: Robust Optimization under Supply/Demand Uncertainty
# Bertsimas–Sim budgeted uncertainty applied to Model 1 (Transportation LP).
#
# The key structural insight: per-station constraints contain only a single
# uncertain parameter each, so the budget Γ has no bite at the station level.
# Instead, we apply the budget at the aggregate level, where many uncertain
# parameters are coupled in a single constraint.
#
# Usage:
#   julia model4_robust.jl              # full dataset
#   julia model4_robust.jl subset_10    # 10-station subset
#   julia model4_robust.jl subset_25    # 25-station subset
#   julia model4_robust.jl subset_50    # 50-station subset
using CSV, DataFrames, JuMP, HiGHS, Statistics

# ─── Parameters ───────────────────────────────────────────────────────────────
const SUBSET_NAME = length(ARGS) > 0 ? ARGS[1] : ""
const DATA_DIR = joinpath(@__DIR__, "..", "data", SUBSET_NAME)

# Uncertainty: deviation magnitude as fraction of nominal
const ALPHA = 0.3  # ±30% of nominal flow

# ─── Load data ───────────────────────────────────────────────────────────────
stations_df = CSV.read(joinpath(DATA_DIR, "stations.csv"), DataFrame)
distances_df = CSV.read(joinpath(DATA_DIR, "distances.csv"), DataFrame)

dist_matrix = Matrix(distances_df[:, 2:end])

surplus_idx = findall(stations_df.net_flow .> 0)
deficit_idx = findall(stations_df.net_flow .< 0)

S = length(surplus_idx)
D = length(deficit_idx)

supply = stations_df.net_flow[surplus_idx]       # s̄_i > 0
demand = -stations_df.net_flow[deficit_idx]       # b̄_j > 0

# Deviation estimates: 30% of nominal
s_hat = ALPHA .* supply
b_hat = ALPHA .* demand

total_s = sum(supply)
total_d = sum(demand)
total_hat_s = sum(s_hat)
total_hat_d = sum(b_hat)

Γ_max = S + D  # full conservatism

println("Model 4: Robust Optimization (Bertsimas–Sim)")
println("  Dataset: $(isempty(SUBSET_NAME) ? "full" : SUBSET_NAME)")
println("  Surplus: $S, Deficit: $D, Total stations: $(S+D)")
println("  Nominal supply: $total_s, Nominal demand: $total_d")
println("  Deviation (α=$ALPHA): supply σ̂=$total_hat_s, demand σ̂=$total_hat_d")
println("  Γ range: 0 → $Γ_max")

# ─── Sweep over Γ ────────────────────────────────────────────────────────────
Γ_values = [0, 1, 2, 3, 5, 10, 20, Γ_max]
results = []
base_cost = 0.0  # initialize for summary

for Γ in Γ_values
    println("\n--- Γ = $Γ ---")

    m = Model(HiGHS.Optimizer)
    set_silent(m)

    # Decision variables: y[i,j] = bikes from surplus i to deficit j
    @variable(m, y[1:S, 1:D] >= 0)

    penalty = maximum(dist_matrix) * 10
    shortfall = max(0, total_d - total_s)
    if shortfall > 0
        @variable(m, y_dummy[1:D] >= 0)
        @objective(m, Min,
            sum(dist_matrix[surplus_idx[i], deficit_idx[j]] * y[i, j] for i in 1:S, j in 1:D) +
            penalty * sum(y_dummy[j] for j in 1:D)
        )
        @constraint(m, demand_con[j in 1:D],
            sum(y[i, j] for i in 1:S) + y_dummy[j] >= demand[j]
        )
    else
        @objective(m, Min,
            sum(dist_matrix[surplus_idx[i], deficit_idx[j]] * y[i, j] for i in 1:S, j in 1:D)
        )
        @constraint(m, demand_con[j in 1:D],
            sum(y[i, j] for i in 1:S) >= demand[j]
        )
    end

    # Per-station supply constraints (nominal)
    @constraint(m, supply_con[i in 1:S],
        sum(y[i, j] for j in 1:D) <= supply[i]
    )

    # ─── Robust aggregate constraints ─────────────────────────────────────
    if Γ > 0
        @variable(m, z_s >= 0)       # supply budget dual
        @variable(m, r_s[1:S] >= 0)  # supply box duals
        @variable(m, z_d >= 0)       # demand budget dual
        @variable(m, t_d[1:D] >= 0)  # demand box duals

        # Supply side: Σy + Γ·z_s + Σr_s ≤ Σs̄
        @constraint(m, robust_supply,
            sum(y[i, j] for i in 1:S, j in 1:D) + Γ * z_s + sum(r_s) <= total_s
        )
        @constraint(m, supply_dual[i in 1:S],
            z_s + r_s[i] >= s_hat[i]
        )

        # Demand side: Σy ≥ Σb̄ + Γ·z_d + Σt_d
        @constraint(m, robust_demand,
            sum(y[i, j] for i in 1:S, j in 1:D) >= total_d + Γ * z_d + sum(t_d)
        )
        @constraint(m, demand_dual[j in 1:D],
            z_d + t_d[j] >= b_hat[j]
        )
    end

    # ─── Solve ───────────────────────────────────────────────────────────────
    optimize!(m)

    status = termination_status(m)
    if status in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL]
        y_val = value.(y)
        real_bikes = sum(y_val)
        dummy_bikes = @isdefined(y_dummy) ? sum(value.(y_dummy)) : 0.0
        nonzero = count(y_val[i, j] > 1e-6 for i in 1:S, j in 1:D)

        real_obj = objective_value(m)
        if @isdefined(y_dummy)
            real_obj -= penalty * sum(value.(y_dummy))
        end
        avg_dist = real_bikes > 0 ? real_obj / real_bikes : 0.0
        is_int = all(y_val[i, j] ≈ round(y_val[i, j]) for i in 1:S, j in 1:D if y_val[i, j] > 1e-6)

        println("  Status: $status")
        println("  Real bikes: $(round(Int, real_bikes)), Dummy: $(round(Int, dummy_bikes))")
        println("  Nonzero arcs: $nonzero, All integer: $is_int")
        println("  Avg distance: $(round(avg_dist; digits=2)) km")
    else
        real_bikes = 0.0
        dummy_bikes = 0.0
        nonzero = 0
        real_obj = 0.0
        avg_dist = 0.0
        is_int = false
        println("  Status: $status (infeasible — supply-demand gap too large for Γ=$Γ)")
    end

    # Price of robustness
    if Γ == 0
        global base_cost = real_obj
        println("  Base cost (Γ=0): $(round(Int, real_obj))")
    elseif status in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL]
        price = real_obj - base_cost
        println("  Price of robustness: $(round(Int, price)) (+$(round(price/base_cost*100; digits=1))%)")
    else
        println("  Price of robustness: ∞ (infeasible)")
    end

    push!(results, (
        Γ = Γ,
        status = status,
        objective = objective_value(m),
        real_cost = real_obj,
        real_bikes = real_bikes,
        dummy_bikes = dummy_bikes,
        nonzero_arcs = nonzero,
        avg_distance = avg_dist,
        all_integer = is_int,
        z_s = Γ > 0 ? value(z_s) : 0.0,
        z_d = Γ > 0 ? value(z_d) : 0.0
    ))
end

# ─── Summary table ───────────────────────────────────────────────────────────
println("\n" * "="^80)
println("ROBUST OPTIMIZATION SUMMARY")
println("="^80)
println(rpad("Γ", 6), rpad("Status", 12), rpad("Real Cost", 12), rpad("Bikes", 8),
        rpad("Arcs", 6), rpad("Avg km", 8), rpad("Price", 10), rpad("%", 6))
println("-"^70)
base_cost_summary = results[1].real_cost
for r in results
    if r.status == MOI.OPTIMAL
        price = r.real_cost - base_cost_summary
        pct = base_cost_summary > 0 ? round(price/base_cost_summary*100; digits=1) : 0.0
        price_str = string(round(Int, price))
        pct_str = string(pct)
    else
        price_str = "∞"
        pct_str = "∞"
    end
    println(rpad(string(r.Γ), 6), rpad(string(r.status), 12),
            rpad(string(round(Int, r.real_cost)), 12), rpad(string(round(Int, r.real_bikes)), 8),
            rpad(string(r.nonzero_arcs), 6), rpad(string(round(r.avg_distance; digits=2)), 8),
            rpad(price_str, 10), rpad(pct_str, 6))
end

# Check feasibility limit
println("\nFeasibility check (Σb̄ + Γ·z_d + Σt_d ≤ Σs̄ - Γ·z_s - Σr_s):")
for r in results
    feasible = r.status == MOI.OPTIMAL
    println("  Γ=$(r.Γ): $(feasible ? "FEASIBLE" : "INFEASIBLE")")
end