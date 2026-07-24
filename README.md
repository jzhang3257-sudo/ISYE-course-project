# CitiBike Bike Rebalancing Optimization

ISyE 524 Course Project — Optimization models for CitiBike rebalancing in Jersey City (May 2026).

## File Structure

```
course project/
├── src/
│   ├── preprocess.py              # Data preprocessing: aggregate tripdata into net flows
│   ├── generate_subsets.py        # Generate 10/25/50-station subsets for progressive testing
│   ├── model1_transportation.jl   # Model 1: Transportation LP
│   ├── model2_routing.jl          # Model 2: Capacitated Vehicle Routing MILP
│   ├── model3_asymmetric.jl       # Model 3: Asymmetric Distance Sensitivity
│   └── model4_robust.jl           # Model 4: Robust Optimization (Bertsimas–Sim)
├── data/
│   ├── JC-202605-citibike-tripdata.csv   # Raw trip data (May 2026)
│   ├── stations.csv                      # Preprocessed station-level net flows
│   ├── distances.csv                     # Pairwise Haversine distances (km)
│   ├── subset_10/                        # 10-station subset (5 surplus + 5 deficit)
│   ├── subset_25/                        # 25-station subset
│   └── subset_50/                        # 50-station subset
├── report/
│   ├── Final Report.ipynb         # Main report with all 4 models & results
│   ├── Final Report Template.ipynb
│   ├── Project_Proposal.md
│   ├── Project_Proposal.pdf
│   └── Project_Proposal_fixed.docx
├── slides/
│   └── Project Overview.pptx
├── figures/
│   └── Project_Proposal/
├── reference/
│   └── sample project/            # Reference sample projects
├── Project.toml                   # Julia project dependencies
└── Manifest.toml                  # Julia dependency manifest
```

## Dependencies

- **Julia 1.x** with packages: `JuMP`, `HiGHS`, `CSV`, `DataFrames`, `Statistics`, `LinearAlgebra`, `Random`
- **Python 3.x** with `pandas` (for preprocessing only)

Install Julia dependencies:
```bash
cd "course project"
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Data Pipeline

1. `preprocess.py` — Reads raw CitiBike trip CSV, computes net flow per station (`outflows − inflows`), outputs `stations.csv` and pairwise `distances.csv`.
2. `generate_subsets.py` — Selects equal numbers of top surplus/deficit stations by `|net_flow|` magnitude to create `subset_10/25/50`.

## Models

All models use **JuMP + HiGHS** solver. Each script accepts an optional subset argument.

### Model 1: Transportation LP
**File:** `src/model1_transportation.jl`

Classic uncapacitated transportation problem. Minimizes total bike-km by moving bikes from surplus stations (supply) to deficit stations (demand) on a bipartite graph. The constraint matrix is **totally unimodular (TU)**, guaranteeing integer solutions from the LP relaxation.

- **Decision variables:** `yᵢⱼ` — bikes moved from surplus station `i` to deficit station `j`
- **Dummy penalty** for infeasibility when total supply < total demand
- **Key result:** 933 bikes, 200 nonzero arcs, avg 1.80 km/bike, solves instantly

### Model 2: Capacitated Vehicle Routing MILP
**File:** `src/model2_routing.jl`

Extends Model 1 with truck routing: a fleet of identical trucks (capacity Q=50) departs from a centroid depot, picks up bikes at surplus stations, drops off at deficit stations, and returns to the depot. Uses **MTZ subtour elimination** with strengthening constraints.

- **Decision variables:** `zᵢⱼₖ` (binary routing), `lᵢⱼₖ` (bike load on arc), `pᵢₖ` (pickups), `dropⱼₖ` (dropoffs), `uᵢₖ` (MTZ positions)
- **Strengthening:** 2-cycle elimination, visit-once, symmetry breaking between trucks
- **Key result:** Optimal at ≤25 stations (≤90s), 36% gap at 50 stations (300s limit), no feasible solution at full scale (224 stations) — documents MTZ scaling limitation honestly
- **Fleet size** adapts to instance: 2 (10 stations) → 6 (full)

### Model 3: Asymmetric Distance Sensitivity
**File:** `src/model3_asymmetric.jl`

Tests whether symmetric Haversine distances adequately represent real asymmetric road networks. Perturbs distances as `d_asym[i,j] = d_sym[i,j] · (1 + εᵢⱼ)` where `ε ~ Uniform(0, δ)`, then compares the optimal plan under symmetric vs. asymmetric distances.

- **Metrics:** plan divergence (% arcs changed), cost gap (% cost increase), volume-weighted change
- **δ tested:** 5%, 10%, 20%
- **Key result:** Even at δ=20%, cost gap ≤2% — Haversine approximation is adequate for operational use

### Model 4: Robust Optimization (Bertsimas–Sim)
**File:** `src/model4_robust.jl`

Applies Bertsimas–Sim budgeted uncertainty to Model 1. Supply and demand are uncertain within ±30% intervals, and the **budget Γ** controls conservatism at the aggregate level via dual auxiliary variables (box + budget constraints).

- **Budget Γ:** sweeps from 0 (nominal) to `|S|+|D|` (full conservatism)
- **Price of robustness:** cost increase relative to Γ=0
- **Key structural finding:** Feasible only at Γ=0; all Γ≥1 are infeasible — the system has a **structural supply-demand gap** that cannot be absorbed under any uncertainty budget >0

## Usage

```bash
# Preprocess data (run once)
python src/preprocess.py
python src/generate_subsets.py

# Run individual models
julia --project=. src/model1_transportation.jl           # full dataset
julia --project=. src/model1_transportation.jl subset_10 # 10-station subset
julia --project=. src/model2_routing.jl subset_25        # 25-station subset
julia --project=. src/model3_asymmetric.jl subset_50     # 50-station subset
julia --project=. src/model4_robust.jl                   # full dataset
```

## Key Findings

| Model | Core Finding |
|-------|-------------|
| Model 1 (LP) | 933 bikes redistributed across 200 arcs; instantaneous solve |
| Model 2 (MILP) | Optimal for ≤25 stations; MTZ fails to scale beyond 50 |
| Model 3 (Sensitivity) | Haversine distance error ≤2% of optimal cost |
| Model 4 (Robust) | System structurally infeasible under any uncertainty (Γ ≥ 1) |

See [Final Report.ipynb](report/Final Report.ipynb) for full results and discussion.
