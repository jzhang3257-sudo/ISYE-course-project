# model1_transportation.jl — ISyE 524 Course Project
# Model 1: Uncapacitated Transportation Problem (LP)
#
# Reads stations.csv and distances.csv, formulates the classic
# transportation LP, solves with HiGHS, and verifies integrality.
#
# Usage:
#   julia model1_transportation.jl              # uses data/ (full dataset)
#   julia model1_transportation.jl subset_10    # uses data/subset_10/
#   julia model1_transportation.jl subset_25    # uses data/subset_25/
#   julia model1_transportation.jl subset_50    # uses data/subset_50/
using CSV, DataFrames, JuMP, HiGHS

# ─── Load data ───────────────────────────────────────────────────────────────
const SUBSET = length(ARGS) > 0 ? ARGS[1] : ""
const DATA_DIR = joinpath(@__DIR__, "..", "data", SUBSET)

stations_df = CSV.read(joinpath(DATA_DIR, "stations.csv"), DataFrame)
distances_df = CSV.read(joinpath(DATA_DIR, "distances.csv"), DataFrame)

# Parse distance matrix: first column is station IDs, rest are distances
dist_ids = string.(names(distances_df)[2:end])  # skip empty first col name
dist_matrix = Matrix(distances_df[:, 2:end])     # n × n

# Separate surplus (supply) and deficit (demand) stations
surplus_idx = findall(stations_df.net_flow .> 0)
deficit_idx = findall(stations_df.net_flow .< 0)

S = length(surplus_idx)
D = length(deficit_idx)

supply = stations_df.net_flow[surplus_idx]       # s_i > 0
demand = -stations_df.net_flow[deficit_idx]       # b_j > 0

println("Model 1: Transportation LP")
println("  Surplus stations: $S, Deficit stations: $D")
println("  Total supply: $(sum(supply)), Total demand: $(sum(demand))")

# ─── Build JuMP model ────────────────────────────────────────────────────────
model = Model(HiGHS.Optimizer)
set_silent(model)

# Decision variables: y[i,j] = bikes moved from surplus i to deficit j
@variable(model, y[1:S, 1:D] >= 0)

# Handle supply shortage: add a dummy surplus node with high penalty cost.
# This ensures feasibility while keeping the real solution meaningful.
total_supply = sum(supply)
total_demand = sum(demand)
shortfall = max(0, total_demand - total_supply)
if shortfall > 0
    println("  Supply shortage: $shortfall bikes (supply=$total_supply, demand=$total_demand)")
    println("  Adding dummy source with penalty cost")
    penalty = maximum(dist_matrix) * 10
    @variable(model, y_dummy[1:D] >= 0)
    @objective(model, Min,
        sum(dist_matrix[surplus_idx[i], deficit_idx[j]] * y[i, j]
            for i in 1:S, j in 1:D) +
        sum(penalty * y_dummy[j] for j in 1:D)
    )
    @constraint(model, demand_con[j in 1:D],
        sum(y[i, j] for i in 1:S) + y_dummy[j] == demand[j]
    )
else
    @objective(model, Min,
        sum(dist_matrix[surplus_idx[i], deficit_idx[j]] * y[i, j]
            for i in 1:S, j in 1:D)
    )
    @constraint(model, demand_con[j in 1:D],
        sum(y[i, j] for i in 1:S) == demand[j]
    )
end

# Supply constraints: don't exceed available bikes at each surplus station
@constraint(model, supply_con[i in 1:S],
    sum(y[i, j] for j in 1:D) <= supply[i]
)

# ─── Solve ───────────────────────────────────────────────────────────────────
optimize!(model)

println("\n" * "="^60)
println("RESULTS")
println("="^60)
println("  Status: $(termination_status(model))")
println("  Objective (total bike-km): $(round(Int, objective_value(model)))")

# ─── Verify integrality ──────────────────────────────────────────────────────
y_val = value.(y)
is_integer_solution = all(y_val[i, j] ≈ round(y_val[i, j])
                          for i in 1:S, j in 1:D if y_val[i, j] > 1e-6)

println("  All y_ij integer? $is_integer_solution")
println("  (Expected: true — transportation constraint matrix is TU)")

# ─── Summary statistics ──────────────────────────────────────────────────────
nonzero_arcs = count(y_val[i, j] > 1e-6 for i in 1:S, j in 1:D)
total_bikes_moved = sum(y_val)
if @isdefined(y_dummy)
    dummy_bikes = sum(value.(y_dummy))
    total_bikes_moved += dummy_bikes
    println("  Dummy bikes (unmet demand): $(round(Int, dummy_bikes))")
end
println("  Nonzero arcs: $nonzero_arcs")
println("  Total bikes moved: $(round(Int, total_bikes_moved))")
if total_bikes_moved > 0
    real_obj = objective_value(model)
    if @isdefined(y_dummy)
        real_obj -= penalty * sum(value.(y_dummy))
    end
    real_bikes = total_bikes_moved - (@isdefined(y_dummy) ? sum(value.(y_dummy)) : 0)
    if real_bikes > 0
        println("  Avg distance per real bike: $(round(real_obj / real_bikes; digits=2)) km")
    end
end

# ─── Top 5 arcs by volume ────────────────────────────────────────────────────
arcs = [(i, j, y_val[i, j]) for i in 1:S, j in 1:D if y_val[i, j] > 1e-6]
sort!(arcs, by = x -> -x[3])
println("\n  Top 5 arcs by bike volume:")
for (i, j, vol) in arcs[1:min(5, end)]
    sid_s = stations_df.station_id[surplus_idx[i]]
    sid_d = stations_df.station_id[deficit_idx[j]]
    d = dist_matrix[surplus_idx[i], deficit_idx[j]]
    println("    $sid_s → $sid_d : $(round(Int, vol)) bikes, $(round(d; digits=2)) km")
end
