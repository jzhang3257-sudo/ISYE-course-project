# model2_routing.jl — ISyE 524 Course Project
# Model 2: Capacitated Vehicle Rebalancing with Routing (MILP)
#
# Uses arc sparsification (k-nearest neighbors) to keep the MILP tractable.
# Each station connects only to its M nearest neighbors + depot, reducing
# binary variables from O(N²·K) to O(N·M·K).  For the full 224-station
# dataset with M=15 this yields ~23k binaries instead of ~304k.
#
# Improvements over the dense version:
#   1. Arc sparsification — only M nearest neighbors per station
#   2. MTZ subtour elimination only on arcs that exist
#   3. Tighter big-M coefficients (min(supply, Q) instead of supply)
#   4. return_empty constraint — trucks must return to depot empty
#   5. Adaptive M_NEIGHBORS per dataset size
#
# Usage:
#   julia model2_routing.jl              # full dataset
#   julia model2_routing.jl subset_10    # 10-station subset
#   julia model2_routing.jl subset_25    # 25-station subset
#   julia model2_routing.jl subset_50    # 50-station subset
using CSV, DataFrames, JuMP, HiGHS, LinearAlgebra, Statistics

# ─── Parameters ───────────────────────────────────────────────────────────────
const SUBSET_NAME = length(ARGS) > 0 ? ARGS[1] : ""
const DATA_DIR = joinpath(@__DIR__, "..", "data", SUBSET_NAME)

# Truck parameters
const Q = 50          # truck capacity (bikes)
const MIP_GAP = 0.05  # 5% optimality gap tolerance
const TIME_LIMIT = 300 # seconds

# Adaptive fleet size based on problem scale
const K_BAR = isempty(SUBSET_NAME) ? 6 :
              SUBSET_NAME == "subset_10" ? 2 :
              SUBSET_NAME == "subset_25" ? 3 :
              SUBSET_NAME == "subset_50" ? 4 : 6

# Arc sparsification: each station connects to its M nearest neighbors + depot.
# This improves tractability but restricts the feasible route set, so the sparse
# solution is an approximation to the full dense routing model.
# Larger M → better solutions but slower solve times
const M_NEIGHBORS = isempty(SUBSET_NAME) ? 15 :
                    SUBSET_NAME == "subset_10" ? 5 :
                    SUBSET_NAME == "subset_25" ? 8 :
                    SUBSET_NAME == "subset_50" ? 10 : 15

# ─── Load data ───────────────────────────────────────────────────────────────
stations_df = CSV.read(joinpath(DATA_DIR, "stations.csv"), DataFrame)
distances_df = CSV.read(joinpath(DATA_DIR, "distances.csv"), DataFrame)

n_stations = nrow(stations_df)
dist_matrix = Matrix(distances_df[:, 2:end])  # n × n

# Separate station types
surplus_idx = findall(stations_df.net_flow .> 0)
deficit_idx = findall(stations_df.net_flow .< 0)
S_set = length(surplus_idx)
D_set = length(deficit_idx)

supply_amt = stations_df.net_flow[surplus_idx]       # s_i > 0
demand_amt = -stations_df.net_flow[deficit_idx]       # b_j > 0

println("Model 2: Capacitated Vehicle Routing MILP (sparse)")
println("  Dataset: $(isempty(SUBSET_NAME) ? "full" : SUBSET_NAME)")
println("  Stations: $n_stations (surplus=$S_set, deficit=$D_set)")
println("  Trucks: K̄=$K_BAR, capacity Q=$Q")
println("  Sparsification: M=$M_NEIGHBORS neighbors per station")
println("  Total supply: $(sum(supply_amt)), Total demand: $(sum(demand_amt))")

# ─── Build depot ─────────────────────────────────────────────────────────────
# Add a synthetic depot at the centroid of all stations
depot_lat = mean(stations_df.lat)
depot_lon = mean(stations_df.lon)
println("  Depot at centroid: ($(round(depot_lat; digits=4)), $(round(depot_lon; digits=4)))")

# Build extended distance matrix: N₀ = {0} ∪ N
function haversine_km(lat1, lon1, lat2, lon2)
    R = 6371.0
    dlat = deg2rad(lat2 - lat1)
    dlon = deg2rad(lon2 - lon1)
    a = sin(dlat/2)^2 + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dlon/2)^2
    return R * 2 * atan(sqrt(a), sqrt(1 - a))
end

N = n_stations
N0 = N + 1  # include depot
d = zeros(N0, N0)
d[2:end, 2:end] = dist_matrix
for i in 1:N
    d[1, i+1] = haversine_km(depot_lat, depot_lon, stations_df.lat[i], stations_df.lon[i])
    d[i+1, 1] = d[1, i+1]  # symmetric
end

# Build supply/demand arrays extended with depot (index 1 = depot, no supply/demand)
supply_all = zeros(N0)
demand_all = zeros(N0)
supply_all[surplus_idx .+ 1] .= supply_amt
demand_all[deficit_idx .+ 1] .= demand_amt

# Index sets for model
K = 1:K_BAR
all_stations = 2:N0        # 2..N+1 (skip depot)
S_idx = surplus_idx .+ 1   # surplus station indices in N₀
D_idx = deficit_idx .+ 1   # deficit station indices in N₀

# ─── Build sparse arc set ────────────────────────────────────────────────────
# For each station, find its M_NEIGHBORS nearest neighbors (excluding self)
# Arc set includes: depot→station, station→depot, and station→M-nearest
actual_M = min(M_NEIGHBORS, N - 1)

# Find nearest neighbors for each station (original 1:N indexing)
neighbors = Vector{Vector{Int}}(undef, N)
for i in 1:N
    sorted_idx = sortperm(dist_matrix[i, :])
    # Exclude self (distance 0) and take top M
    nbrs = [j for j in sorted_idx if j != i][1:actual_M]
    neighbors[i] = nbrs
end

# Build sparse arc list as Vector{Tuple{Int,Int}} (extended indexing: depot=1)
sparse_arcs = Tuple{Int,Int}[]

# Depot ↔ all stations
for j in all_stations
    push!(sparse_arcs, (1, j))    # depot → station
    push!(sparse_arcs, (j, 1))    # station → depot
end

# Station → M nearest stations (in extended indexing)
for i in 1:N
    for j in neighbors[i]
        push!(sparse_arcs, (i + 1, j + 1))
    end
end

# Precompute outgoing / incoming arc sets for each node (fast lookups)
out_arcs = [Int[] for _ in 1:N0]
in_arcs  = [Int[] for _ in 1:N0]
for (i, j) in sparse_arcs
    push!(out_arcs[i], j)
    push!(in_arcs[j], i)
end

# MTZ arcs: station-to-station arcs only (no depot)
mtz_arcs = Tuple{Int,Int}[]
for (i, j) in sparse_arcs
    if i != 1 && j != 1 && i != j
        push!(mtz_arcs, (i, j))
    end
end

# 2-cycle arcs: pairs where both directions exist
two_cycle_arcs = Tuple{Int,Int}[]
for (i, j) in sparse_arcs
    if i != 1 && j != 1 && i < j && (j, i) in sparse_arcs
        push!(two_cycle_arcs, (i, j))
    end
end

n_arcs = length(sparse_arcs)
n_mtz = length(mtz_arcs)
n_z_vars = n_arcs * K_BAR
println("  Sparse arcs: $n_arcs (station→station + depot links)")
println("  MTZ arcs: $n_mtz")
println("  Binary variables z: $n_z_vars (vs $(N0*N0*K_BAR) dense)")

# ─── Build JuMP model ────────────────────────────────────────────────────────
m = Model(HiGHS.Optimizer)
set_silent(m)
set_time_limit_sec(m, TIME_LIMIT)
set_attribute(m, "mip_rel_gap", MIP_GAP)

# Decision variables — sparse: only over arcs that exist
@variable(m, z[(i, j) in sparse_arcs, k in K], Bin)     # routing
@variable(m, l[(i, j) in sparse_arcs, k in K] >= 0)     # bike load on arc
@variable(m, p[i in S_idx, k in K] >= 0, Int)            # pickups at surplus stations
@variable(m, drop[i in D_idx, k in K] >= 0, Int)         # dropoffs at deficit stations
@variable(m, u[all_stations, k in K] >= 0)               # MTZ position

# Track unmet demand with soft penalty
total_demand = sum(demand_amt)
@variable(m, unmet[j in D_idx] >= 0)                     # bikes not delivered
penalty = maximum(d) * 10

# ─── Objective: minimize total distance + unmet demand penalty ───────────────
@objective(m, Min,
    sum(d[i, j] * z[(i, j), k] for (i, j) in sparse_arcs, k in K) +
    penalty * sum(unmet[j] for j in D_idx)
)

# ─── Constraints ─────────────────────────────────────────────────────────────

# (1) Supply: total pickups ≤ available
@constraint(m, supply_con[i in S_idx],
    sum(p[i, k] for k in K) <= supply_all[i]
)

# (2) Demand: total dropoffs + unmet = required
@constraint(m, demand_con[j in D_idx],
    sum(drop[j, k] for k in K) + unmet[j] == demand_all[j]
)

# (3) Flow conservation: in = out at each station (using sparse arcs)
@constraint(m, flow_con[i in all_stations, k in K],
    sum(z[(i, j), k] for j in out_arcs[i]) ==
    sum(z[(j, i), k] for j in in_arcs[i])
)

# (4) Each truck departs depot at most once
@constraint(m, depart_con[k in K],
    sum(z[(1, j), k] for j in out_arcs[1]) <= 1
)

# (5) Each truck returns to depot at most once
@constraint(m, return_con[k in K],
    sum(z[(i, 1), k] for i in in_arcs[1]) <= 1
)

# (6) Fleet size limit
@constraint(m, fleet_con,
    sum(z[(1, j), k] for j in out_arcs[1], k in K) <= K_BAR
)

# (7) Bike mass balance: net load change = pickups - dropoffs
@constraint(m, load_balance[i in all_stations, k in K],
    sum(l[(i, j), k] for j in out_arcs[i]) -
    sum(l[(j, i), k] for j in in_arcs[i]) ==
    (i in S_idx ? p[i, k] : 0) - (i in D_idx ? drop[i, k] : 0)
)

# (8) Capacity linking: load ≤ Q only when arc is used
@constraint(m, capacity_con[(i, j) in sparse_arcs, k in K],
    l[(i, j), k] <= Q * z[(i, j), k]
)

# (9) Start empty from depot
@constraint(m, start_empty[j in out_arcs[1], k in K],
    l[(1, j), k] == 0
)

# (10) Return empty to depot (FIX: was missing in dense version)
@constraint(m, return_empty[i in in_arcs[1], k in K],
    l[(i, 1), k] == 0
)

# (11) Pickup only if visited — tighter big-M: min(supply, Q)
@constraint(m, pickup_if_visited[i in S_idx, k in K],
    p[i, k] <= min(supply_all[i], Q) * sum(z[(j, i), k] for j in in_arcs[i])
)

# (12) Dropoff only if visited — tighter big-M: min(demand, Q)
@constraint(m, dropoff_if_visited[i in D_idx, k in K],
    drop[i, k] <= min(demand_all[i], Q) * sum(z[(i, j), k] for j in out_arcs[i])
)

# (13) MTZ subtour elimination — only on station→station arcs that exist
@constraint(m, mtz[(i, j) in mtz_arcs, k in K],
    u[i, k] - u[j, k] + N * z[(i, j), k] <= N - 1
)

# (14) MTZ bounds
@constraint(m, mtz_bound[i in all_stations, k in K],
    1 <= u[i, k] <= N
)

# (15) 2-cycle elimination — only for pairs where both directions exist
@constraint(m, two_cycle[(i, j) in two_cycle_arcs, k in K],
    z[(i, j), k] + z[(j, i), k] <= 1
)

# (16) Each station visited at most once (by at most one truck)
# Relaxed for stations with supply/demand exceeding truck capacity Q
@constraint(m, visit_once[i in all_stations],
    sum(z[(j, i), k] for j in in_arcs[i], k in K) <=
    max(1, ceil(max(supply_all[i], demand_all[i]) / Q))
)

# (17) Symmetry breaking: truck k only used if truck k-1 is used
@constraint(m, symmetry[k in 2:K_BAR],
    sum(z[(1, j), k] for j in out_arcs[1]) <=
    sum(z[(1, j), k-1] for j in out_arcs[1])
)

# ─── Solve ───────────────────────────────────────────────────────────────────
println("\nSolving...")
optimize!(m)

status = termination_status(m)
println("Status: $status")

if status in [MOI.OPTIMAL, MOI.TIME_LIMIT, MOI.ALMOST_OPTIMAL]
    println("Objective: $(round(Int, objective_value(m)))")
    println("Best bound: $(round(Int, objective_bound(m)))")
    gap = (objective_value(m) - objective_bound(m)) / objective_value(m) * 100
    println("Relative gap: $(round(gap; digits=1))%")
    println("Solve time: $(round(solve_time(m); digits=1)) s")

    # Extract results
    z_val = value.(z)
    l_val = value.(l)
    p_val = value.(p)
    d_val = value.(drop)
    unmet_val = value.(unmet)

    # Count used trucks
    used_trucks = sum(sum(z_val[(1, j), k] for j in out_arcs[1]) > 0.5 for k in K)
    println("Trucks used: $(round(Int, used_trucks)) / $K_BAR")

    # Total bikes handled
    total_picked = sum(p_val)
    total_dropped = sum(d_val)
    total_unmet = sum(unmet_val)
    println("Bikes picked up: $(round(Int, total_picked))")
    println("Bikes dropped off: $(round(Int, total_dropped))")
    println("Unmet demand: $(round(Int, total_unmet))")

    # Route summary
    println("\nTruck routes:")
    for k in K
        used = sum(z_val[(1, j), k] for j in out_arcs[1]) > 0.5
        if !used
            println("  Truck $k: not used")
            continue
        end
        # Trace route from depot using sparse arcs
        route = [1]
        current = 1
        visited = Set{Int}()
        max_steps = N + 2  # safety limit
        step = 0
        while step < max_steps
            step += 1
            next_stops = [j for j in out_arcs[current] if z_val[(current, j), k] > 0.5]
            if isempty(next_stops) || current in visited
                break
            end
            push!(visited, current)
            current = next_stops[1]
            push!(route, current)
            if current == 1  # back to depot
                break
            end
        end
        # Format route
        route_str = String[]
        for node in route
            if node == 1
                push!(route_str, "Depot")
            else
                sid = stations_df.station_id[node-1]
                if node in S_idx
                    pk = p_val[node, k]
                    if pk > 0.5
                        push!(route_str, "$sid(+$(round(Int, pk)))")
                    else
                        push!(route_str, sid)
                    end
                elseif node in D_idx
                    dp = d_val[node, k]
                    if dp > 0.5
                        push!(route_str, "$sid(-$(round(Int, dp)))")
                    else
                        push!(route_str, sid)
                    end
                else
                    push!(route_str, sid)
                end
            end
        end
        # Calculate route distance
        route_dist = sum(d[route[i], route[i+1]] for i in 1:length(route)-1)
        println("  Truck $k: $(join(route_str, " → "))")
        println("          Distance: $(round(route_dist; digits=2)) km")
    end

    # Total distance
    total_dist = sum(d[i, j] * z_val[(i, j), k] for (i, j) in sparse_arcs, k in K)
    println("\nTotal truck distance: $(round(Int, total_dist)) km")

    # Compare with Model 1 lower bound
    real_obj = objective_value(m) - penalty * total_unmet
    if total_dropped > 0
        println("Avg distance per delivered bike: $(round(real_obj / total_dropped; digits=2)) km")
    end
else
    println("Model did not reach optimality.")
    println("Termination status: $status")
end