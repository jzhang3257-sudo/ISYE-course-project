# Project Proposal — ISyE/CS/ECE 524 (Summer 2026)

## Bike Rebalancing in Jersey City: A Network Flow Approach to CitiBike Redistribution

**Team Members:** Xiaonan Meng (email), Junkai Zhang (email)

---

## 1. Problem Statement

Bike‑sharing systems face a persistent operational challenge: by the end of each day,
ridership patterns leave some stations **overflowing with bikes** while others sit
**nearly empty**. To prepare for the next morning's rush, an operator must deploy a fleet
of trucks to **reposition bikes from surplus stations to deficit stations**. The goal is
to decide, for each truck, which stations to visit and how many bikes to move along each
route, so that total **transportation cost is minimized** while **all (or most) demand
is satisfied**.

This is a classic **minimum‑cost network flow** problem with integer constraints, which
naturally connects to the linear programming, integer programming, and transportation
models covered in class. It is *modeling‑focused* (no algorithm tuning, no simulation, no
machine learning) and scales from a hand‑verifiable toy network to a realistic 157‑station
system using **real CitiBike trip data** from Jersey City.

We will solve the problem as a sequence of four progressively more sophisticated models,
mirroring the structure of high‑quality project reports from previous semesters:
(1) a pure LP transportation model with a total unimodularity proof,
(2) a capacitated MILP with truck routing and explicit MTZ subtour elimination,
(3) an asymmetric‑distance sensitivity analysis using real road‑network distances, and
(4) a robust optimization model with full Bertsimas–Sim dual derivation.
A cross‑cutting multi‑objective Pareto analysis of cost vs. unmet demand is applied to
Model 2 in the Results section.

---

## 2. Motivation & Background

Bike‑sharing has grown into one of the most visible urban transportation innovations of
the past decade. As of 2024, over 2,000 cities worldwide operate bike‑share programs, and
CitiBike — serving New York City and Jersey City — is the largest in North America with
over 30 million annual trips [1].

The **bike rebalancing problem (BRP)** was first formulated as a mathematical optimization
problem in the early 2010s as systems scaled beyond what manual redistribution could handle
[3]. It belongs to the family of **vehicle routing problems with pickup and delivery
(VRPPD)** and can be viewed as a generalization of the classical **transportation problem**
[2], with the added complexity of vehicle capacity constraints and (in its full form)
routing decisions.

A landmark 2013 paper by Raviv, Tzur, and Forma [3] formulated the static BRP as a
mixed‑integer program, showing that even moderate‑size instances (100+ stations) are
tractable with commercial solvers. Since then, researchers have explored dynamic
(time‑varying), stochastic (demand‑uncertain), and robust formulations [4].

**Why this project is interesting:**

- The data is **real, public, and immediately usable**: CitiBike publishes trip‑level CSV
  files each month with start/end stations, timestamps, and GPS coordinates.
- Rebalancing is a **genuinely hard operations problem** where small changes in truck
  capacity, fleet size, or demand patterns yield markedly different optimal solutions.
- The model family spans **LP → MILP → multi‑objective → robust optimization**, providing
  rich material for the Results & Discussion section.
- The station network has a natural **geographic layout**, enabling compelling map
  visualizations that make the results immediately interpretable.

**Simplifying assumptions (toy version):**

- Rebalancing occurs **overnight** (static single-period model). We do not model time
  windows or traffic.
- Each truck starts and ends at a depot, has a fixed capacity (number of bikes), and
  travels directly between stations. Routing (TSP) constraints are simplified to
  origin–destination arcs in the base model.
- Demand at each station is the **net bike flow** (arrivals minus departures) aggregated
  from trip data, treated as deterministic.
- All stations can be visited by any truck.

---

## 3. Data

### Source

We use the **CitiBike Jersey City trip data** for June 2024, publicly available at
[https://s3.amazonaws.com/tripdata/index.html](https://s3.amazonaws.com/tripdata/index.html).

### Raw Data Snapshot

| Field | Example |
|---|---|
| `ride_id` | `A17A516F45135535` |
| `started_at` | `2024-06-01 21:57:43` |
| `ended_at` | `2024-06-01 22:02:11` |
| `start_station_id` | `HB301` |
| `end_station_id` | `HB501` |
| `start_lat`, `start_lng` | `40.7423, -74.0351` |
| `end_lat`, `end_lng` | `40.7482, -74.0325` |
| `member_casual` | `casual` |

### Scale (June 2024, Jersey City)

| Metric | Value |
|---|---|
| Total trips | 111,115 |
| Unique stations | 157 |
| Data size (zipped) | 4.4 MB |
| Missing values (station ID) | 1 row |

### Data Processing Pipeline

1. **Download** the zipped CSV from the CitiBike S3 bucket.
2. **Aggregate** trips by station: for each station $i$, compute

   $$\text{net}_i = \text{arrivals}_i - \text{departures}_i$$

   Stations with $\text{net}_i > 0$ are **surplus** (supply nodes); stations with
   $\text{net}_i < 0$ are **deficit** (demand nodes).
3. **Build the distance matrix** $d_{ij}$ between every pair of stations using the
   Haversine formula applied to the provided latitude/longitude coordinates. We also
   include a **depot** node (e.g., a central maintenance facility) from which all trucks
   depart.

All processing is done with a short Python script; the resulting aggregated CSV (157 rows)
is imported into Julia for modeling.

---

## 4. Toy Example: Hand‑Verified 5‑Station Network

Before scaling to the full Jersey City dataset, we validate our approach on a miniature
5‑station network that can be solved by hand. This ensures our JuMP implementation is
correct and provides intuition before tackling 157 stations.

### Toy Network Data

| Station | Type | Net flow | Interpretation |
|---|---|---|---|
| A | Surplus | +8 | 8 bikes available to move |
| B | Surplus | +5 | 5 bikes available to move |
| C | Deficit | −6 | needs 6 bikes |
| D | Deficit | −4 | needs 4 bikes |
| E | Deficit | −3 | needs 3 bikes |

Total supply: $8+5=13$, total demand: $6+4+3=13$ — the system is balanced.

**Transportation cost matrix $d_{ij}$ (km, Haversine), surplus rows × deficit columns:**

| | C | D | E |
|---|---|---|---|
| **A** | 2.1 | 3.5 | 4.0 |
| **B** | 1.8 | 2.9 | 3.3 |

(Only surplus‑to‑deficit distances enter the transportation model; the full symmetric
$5 \times 5$ matrix is used later when truck routing between arbitrary station pairs is
introduced in Model 2.)

### Hand Solution (Model 1 — LP Transportation)

By inspection, the greedy minimum‑cost assignment is:
- A → C: 6 bikes @ 2.1 km = 12.6 bike‑km
- A → D: 2 bikes @ 3.5 km = 7.0 bike‑km
- B → D: 2 bikes @ 2.9 km = 5.8 bike‑km
- B → E: 3 bikes @ 3.3 km = 9.9 bike‑km

**Total:** 35.3 bike‑km. This is optimal (verified by enumeration of all feasible
integer assignments). We will confirm that Model 1's LP solver matches this exactly.

---

## 5. Mathematical Models

### Notation (shared across all models)

- $N$ = set of stations, indexed by $i, j$; $N_0 = N \cup \{0\}$ where $0$ is the depot.
- $S = \{i \in N : \text{net}_i > 0\}$ — surplus stations (supply).
- $D = \{i \in N : \text{net}_i < 0\}$ — deficit stations (demand).
- $s_i = \text{net}_i$ (bikes available at $i \in S$).
- $b_j = -\text{net}_j$ (bikes needed at $j \in D$).
- $d_{ij}$ = distance from station $i$ to station $j$ (km, via Haversine).
- $K$ = set of trucks, each with capacity $Q$ (bikes); $\bar{K}$ = fleet size limit.

---

### Model 1 — Uncapacitated Transportation Problem (LP)

*This is the simplest formulation: ignore trucks and treat each surplus station as a pure
source and each deficit station as a pure sink. It provides a lower bound for all
subsequent models and is solved by hand on the 5‑station toy network above.*

**Decision variable:** $y_{ij} \ge 0$ — number of bikes moved from surplus $i$ to deficit $j$.

$$
\begin{aligned}
\underset{y}{\text{minimize}}\quad & \sum_{i \in S}\sum_{j \in D} d_{ij}\, y_{ij} \\
\text{subject to}\quad
& \sum_{j \in D} y_{ij} \le s_i && \forall i \in S \quad\text{(supply capacity)}\\
& \sum_{i \in S} y_{ij} = b_j && \forall j \in D \quad\text{(demand satisfaction)}\\
& y_{ij} \ge 0 && \forall i\in S,\; j\in D.
\end{aligned}
$$

**Why the LP solution is automatically integral.** The transportation problem's
constraint matrix, after converting the supply inequalities to equalities via slack
variables, is (up to row signs) the **node–arc incidence matrix** of a directed bipartite
graph augmented with slack columns: each flow column $y_{ij}$ has exactly one $+1$ (at the
source row) and exactly one $-1$ (at the sink row), and each slack column has a single
$+1$. A classical sufficient condition for total unimodularity then applies: a matrix with
all entries in $\{0,\pm 1\}$ in which every column contains at most one $+1$ and at most
one $-1$ is **totally unimodular (TU)**. By the Hoffman–Kruskal theorem, a constraint
matrix is TU if and only if the polyhedron $\{x : Ax \le c,\, x \ge 0\}$ has integer
vertices for every integer right‑hand side $c$; hence with integer $s_i, b_j$, every basic
feasible solution is integer‑valued. Therefore, solving the LP with HiGHS automatically
yields integer $y_{ij}$ — no branching required.

Importantly, TU depends only on the **constraint matrix**, not the objective
coefficients. The property holds regardless of whether distances are symmetric or
asymmetric — a point we return to when discussing Model 3.

---

### Model 2 — Capacitated Vehicle Rebalancing with Routing (MILP)

*We now introduce trucks. Each truck has finite capacity $Q$, starts and ends at the depot,
picks up bikes at surplus stations, and drops them off at deficit stations along a
multi‑stop route. The critical insight is that **the route a truck travels and the
origin‑destination pairing of bikes are separate concepts**: a truck may travel
depot→A→B→C→depot, picking up at A, passing through B, and dropping off at C. We
therefore separate the routing variables ($z$) from the load‑tracking variables ($l$)
and the pickup/dropoff variables ($p$, $d$).*

**Decision variables:**

- $z_{ijk} \in \{0,1\}$ — truck $k$ traverses arc $(i,j)$ (routing).
- $l_{ijk} \in \mathbb{R}_{\ge 0}$ — bike load **on board** truck $k$ when traversing arc $(i,j)$.
- $p_{ik} \in \mathbb{Z}_{\ge 0}$ — bikes picked up by truck $k$ at station $i \in S$.
- $d_{jk} \in \mathbb{Z}_{\ge 0}$ — bikes dropped off by truck $k$ at station $j \in D$.
- $u_{ik} \in \mathbb{R}_{\ge 0}$ — auxiliary variable for MTZ subtour elimination.

**Complete formulation:**

$$
\begin{aligned}
\underset{z,\,l,\,p,\,d,\,u}{\text{minimize}}\quad & \sum_{k \in K}\sum_{i \in N_0}\sum_{j \in N_0} d_{ij}\, z_{ijk} \\
\text{subject to}\quad
& \sum_{k\in K} p_{ik} \le s_i && \forall i \in S &&\text{(1) supply}\\
& \sum_{k\in K} d_{jk} = b_j && \forall j \in D &&\text{(2) demand}\\
& \sum_{j\in N_0} z_{ijk} = \sum_{j\in N_0} z_{jik} && \forall i\in N,\; k\in K &&\text{(3) truck flow}\\
& \sum_{j\in N} z_{0jk} \le 1 && \forall k\in K &&\text{(4) depot departure}\\
& \sum_{i\in N} z_{i0k} \le 1 && \forall k\in K &&\text{(5) depot return}\\
& \sum_{k\in K}\sum_{j\in N} z_{0jk} \le \bar{K} && &&\text{(6) fleet size}\\
& \sum_{j\in N_0} l_{ijk} - \sum_{j\in N_0} l_{jik}
   = p_{ik} - d_{ik} && \forall i\in N,\; k\in K &&\text{(7) bike mass balance}\\
& l_{ijk} \le Q \cdot z_{ijk} && \forall i,j\in N_0,\; k\in K &&\text{(8) capacity linking}\\
& l_{0jk} = 0 && \forall j\in N,\; k\in K &&\text{(9) start empty}\\
& p_{ik} \le s_i \sum_{j\in N_0} z_{jik} && \forall i\in S,\; k\in K &&\text{(10) pickup if visited}\\
& d_{jk} \le b_j \sum_{i\in N_0} z_{ijk} && \forall j\in D,\; k\in K &&\text{(11) dropoff if visited}\\
& u_{ik} - u_{jk} + |N| \cdot z_{ijk} \le |N|-1 && \forall i,j\in N,\; i \neq j,\; k\in K &&\text{(12) MTZ subtour}\\
& z_{ijk} \in \{0,1\},\; l_{ijk} \ge 0,\; p_{ik},\,d_{jk} \in \mathbb{Z}_{\ge 0},\; u_{ik} \ge 0.
\end{aligned}
$$

**Why this formulation is correct.** Constraint (7) is the key addition: for each truck
$k$ at each station $i$, the net outflow of on‑board bikes equals the net
pickup/dropoff at that station. If station $i$ is a surplus station where truck $k$
picks up $p_{ik}$ bikes, the bike load on departing arcs exceeds the load on arriving
arcs by $p_{ik}$. If $i$ is a deficit station where truck $k$ drops off $d_{ik}$ bikes,
the departing load is *less* than the arriving load by $d_{ik}$. Constraints (10)–(11)
ensure pickups/dropoffs only happen at stations the truck actually visits. Constraint
(9) ensures trucks start empty from the depot, and constraint (5) combined with (7)
implies they return empty ($\sum_j l_{j0k}=0$).

Note the contrast with the naive approach that would set
$x_{ijk} =$ "bikes moved directly from $i$ to $j$" and link it to $z_{ijk}$.
With multi‑stop routing, the origin and destination of a bike are not necessarily
adjacent in the truck's route — the truck may carry bikes through intermediate stations.
Our $l_{ijk}$ formulation correctly tracks the *cumulative load* along the route.

**Why MTZ over DFJ.** The MTZ constraints (12) are polynomial: $O(|N|^2)$ per truck,
manageable for 157 stations. The alternative DFJ formulation requires an exponential
number of subset constraints and would demand a cutting‑plane (branch‑and‑cut) solver
callback — this constitutes *algorithm tuning*, which the course guidelines explicitly
discourage. MTZ's LP relaxation is weaker than DFJ's, but HiGHS' MILP solver handles
the resulting branch nodes efficiently at our scale. We document the integrality gap
(lower bound from LP relaxation vs. MILP optimum) in Section 8.2.

### Discussion: Integrality Gap

Unlike Model 1, the MILP with MTZ constraints will generally have a nonzero integrality
gap — the LP relaxation (with $z_{ijk} \in [0,1]$) produces a strictly smaller objective
than the MILP optimum. Quantifying this gap and the **root‑node gap** (gap after presolve
but before branching) is a key result in Section 8.2.

---

### Model 3 — Asymmetric Driving Distances: Sensitivity of the Rebalancing Plan to the Distance Metric (MILP)

*Model 2 uses Haversine (great‑circle) distances, which are symmetric ($d_{ij} = d_{ji}$).
In reality, one‑way streets, traffic patterns, and road geometry make actual travel
distances **asymmetric**: the shortest path from A to B may differ from B to A. This
raises a practical question: does the optimal rebalancing plan change when we use true
street‑network distances instead of straight‑line approximations?*

**What changes.** We replace the Haversine distance matrix $d_{ij}^{\text{sym}}$ with a
**driving distance matrix** $d_{ij}^{\text{asym}}$ obtained from **OpenStreetMap** via the
`osrm` or `openrouteservice` API (a single batch query for $157 \times 157$ pairs, run
once). The formulation of Models 1 and 2 is otherwise unchanged — same constraints,
same variables, different numbers in the objective.

**Important: the constraint matrix remains TU.** Changing the cost coefficients in the
objective function does **not** change total unimodularity — TU is a property of the
constraint matrix, not the objective. Model 1 with asymmetric distances will still
produce integer solutions. The structural change from Model 2 to Model 3 is therefore
**data‑driven**, not constraint‑driven; but the research question is genuinely important:
by how much does the optimal plan change when using real road networks?

We investigate three specific questions:

| Experiment | Question | Metric |
|---|---|---|
| (a) Plan divergence | How many station pairings differ between symmetric vs. asymmetric optimal plans? | % of arcs that change |
| (b) Cost gap | What is the percentage increase in total distance when evaluating the symmetric‑optimal plan on asymmetric distances? | $\frac{\text{asym-cost}(\text{sym-plan}) - \text{asym-cost}(\text{asym-plan})}{\text{asym-cost}(\text{asym-plan})}$ |
| (c) Integrality gap comparison | Does the MILP integrality gap (LP relaxation vs. MILP optimum) differ between symmetric and asymmetric distance matrices? | Gap % for each |

This mirrors a common industrial concern: can we approximate driving distances with
Haversine, or does it lead to materially different (and suboptimal) operational plans?

---

### Model 4 — Robust Optimization Under Supply/Demand Uncertainty

*The net bike flow $\text{net}_i$ computed from one month's data is a point estimate.
On any given night the true supply and demand may deviate. We apply **budgeted
uncertainty** (Bertsimas & Sim [5]) to protect against up to $\Gamma$ stations
simultaneously realizing their worst‑case deviations. We present the full derivation
of the robust counterpart, going from the uncertain constraint through the inner
maximization to the dualized linear formulation.*

**Step 1 — Uncertainty model.** Let the true supply at surplus station $i$ and true
demand at deficit station $j$ be:

$$\begin{aligned}
s_i(\zeta) &= \bar{s}_i + \hat{s}_i \,\zeta_i^s, \quad |\zeta_i^s| \le 1,\; \forall i\in S,\\
b_j(\zeta) &= \bar{b}_j + \hat{b}_j \,\zeta_j^d, \quad |\zeta_j^d| \le 1,\; \forall j\in D,
\end{aligned}$$

where $\bar{s}_i, \bar{b}_j$ are nominal values from June 2024 data and
$\hat{s}_i, \hat{b}_j$ are deviation magnitudes (e.g., standard deviation across days).
The perturbations are constrained by a **budget of uncertainty** $\Gamma$:

$$\mathcal{U}_\Gamma = \Big\{ (\zeta^s,\zeta^d) \;\Big|\; |\zeta_i^s|\le 1,\; |\zeta_j^d|\le 1,\;
   \sum_{i\in S} |\zeta_i^s| + \sum_{j\in D} |\zeta_j^d| \le \Gamma \Big\}.$$

$\Gamma = 0$ recovers the deterministic Model 2; $\Gamma = |S| + |D|$ is the fully
conservative (Soyster) solution in which every station simultaneously realizes its
worst‑case deviation.

**Step 2 — Why per‑station constraints degenerate, and where the budget bites.**
Each per‑station supply constraint $\sum_{k} p_{ik} \le s_i(\zeta)$ contains only a
*single* uncertain parameter, $\zeta_i^s$. Robust feasibility of an inequality is checked
constraint by constraint, so the worst case for station $i$ over $\mathcal{U}_\Gamma$ is
simply $\zeta_i^s = -\min(1,\Gamma)$, giving the counterpart

$$\sum_{k\in K} p_{ik} \;\le\; \bar{s}_i - \hat{s}_i \min(1,\Gamma) \qquad \forall i \in S,$$

and analogously $\sum_{k} d_{jk} \ge \bar{b}_j + \hat{b}_j \min(1,\Gamma)$ for demand.
For any $\Gamma \ge 1$ this collapses to the fully conservative Soyster solution: the
budget provides no interpolation, because a single‑parameter constraint cannot "share"
the budget with anyone. This is an important structural observation in its own right —
the celebrated Bertsimas–Sim machinery only has bite when **many uncertain parameters
appear in the same constraint**.

We therefore apply the budget where it genuinely couples the perturbations: at the
**system level**. We keep the per‑station constraints of Model 2 at their nominal
right‑hand sides, with one necessary modification: the demand equality (2) is relaxed to
$\sum_{k} d_{jk} \ge \bar{b}_j$. (This relaxation is required for consistency — an
equality would pin total dropoffs to exactly $\sum_j \bar{b}_j$, contradicting the
aggregate robust demand buffer introduced below.) In addition,
we require that the plan carry enough *aggregate* slack — spare pickup headroom on the
supply side and spare delivered bikes on the demand side — to absorb *any* realization
in $\mathcal{U}_\Gamma$, with the operator re‑allocating bikes locally on the night of
execution. The two aggregate uncertain constraints are:

$$\sum_{i\in S}\sum_{k\in K} p_{ik} \;\le\; \sum_{i\in S} s_i(\zeta)
\qquad\text{and}\qquad
\sum_{j\in D}\sum_{k\in K} d_{jk} \;\ge\; \sum_{j\in D} b_j(\zeta)
\qquad \forall \zeta \in \mathcal{U}_\Gamma.$$

Each of these couples $|S|$ (resp. $|D|$) uncertain parameters in a single constraint —
exactly the setting Bertsimas–Sim addresses.

**Step 3 — Inner maximization and its LP dual (supply side).** The aggregate supply
constraint must hold at the worst case, i.e.,

$$\sum_{i\in S}\sum_{k\in K} p_{ik} \;\le\; \sum_{i\in S} \bar{s}_i \;-\; \beta^s(\Gamma),
\qquad
\beta^s(\Gamma) = \max_{\eta}\Big\{ \sum_{i\in S} \hat{s}_i\,\eta_i \;\Big|\;
0 \le \eta_i \le 1 \;\forall i \in S,\;\; \sum_{i\in S} \eta_i \le \Gamma \Big\},$$

where $\eta_i = -\zeta_i^s$. (At the optimum the adversary spends no budget on the
demand‑side perturbations $\zeta^d$, since they do not appear in this constraint, so
they can be dropped from the inner problem without loss.) The inner maximization is a
bounded knapsack‑type LP; strong duality holds since it is feasible ($\eta = 0$) and
bounded. Its dual is:

$$
\begin{aligned}
\beta^s(\Gamma) = \min_{z,\,r}\quad & \Gamma \cdot z + \sum_{i\in S} r_{i} \\
\text{subject to}\quad & z + r_{i} \ge \hat{s}_{i} \quad \forall i\in S,\\
& z \ge 0,\; r_{i} \ge 0 \;\forall i \in S.
\end{aligned}
$$

Here $z$ is the dual variable of the budget constraint and $r_i$ is the dual variable of
the box constraint $\eta_i \le 1$.

**Step 4 — Robust counterpart (supply side).** Since any dual‑feasible $(z, r)$ gives an
upper bound on $\beta^s(\Gamma)$, embedding the dual *minimization* into our (outer)
minimization yields the exact **linear robust counterpart**:

$$\boxed{\;\sum_{i\in S}\sum_{k\in K} p_{ik} + \Gamma \cdot z + \sum_{i\in S} r_{i} \;\le\; \sum_{i\in S}\bar{s}_i\;}$$

with $z \ge 0,\; r_{i} \ge 0,\; z + r_{i} \ge \hat{s}_{i}\;\forall i\in S.$

**Step 5 — Robust counterpart (demand side).** For the aggregate demand constraint the
adversary pushes demand *up* ($\zeta_j^d = +1$ on the budgeted stations), so we need
$\sum_{j,k} d_{jk} \ge \sum_j \bar{b}_j + \beta^d(\Gamma)$ with
$\beta^d(\Gamma) = \max\{\sum_j \hat{b}_j \eta_j : 0 \le \eta \le 1, \sum_j \eta_j \le \Gamma\}$.
An identical dual derivation with dual variables $w$ (budget) and $t_j$ (boxes) yields:

$$\boxed{\;\sum_{j\in D}\sum_{k\in K} d_{jk} \;\ge\; \sum_{j\in D}\bar{b}_j + \Gamma \cdot w + \sum_{j'\in D} t_{j'}\;}$$

with $w \ge 0,\; t_{j'} \ge 0,\; w + t_{j'} \ge \hat{b}_{j'}\;\forall j'\in D.$

(One caveat, stated explicitly: since the per‑station demand constraints are kept at
their *nominal* right‑hand sides, the aggregate buffer guarantees system‑wide capacity to
absorb any $\Gamma$‑budget realization, but assumes bikes can be re‑allocated among
nearby deficit stations at execution time. A fully station‑level guarantee would require
the Soyster‑type counterparts of Step 2, which we also report as the conservative
benchmark at $\Gamma \ge 1$.)

**Step 6 — Full robust MILP.** We augment Model 2 (with (2) relaxed to $\ge$ as above)
with the two boxed aggregate robust
constraints and the auxiliary variables $z, r_i, w, t_{j'}$ (all continuous, $O(|S|+|D|)$
of them). Because total dropoffs equal total pickups (mass balance (7) with empty start
and return), overall feasibility requires
$\sum_j \bar{b}_j + \beta^d(\Gamma) \le \sum_i \bar{s}_i - \beta^s(\Gamma)$; we verify
this condition on the data for each $\Gamma$ and report the largest feasible budget.
The resulting formulation is a modestly larger MILP — solved with HiGHS for a
sweep of $\Gamma \in \{0, 1, 2, 5, 10, 20, 50, |S|+|D|\}$.

**Step 7 — Analysis plan.** For each $\Gamma$ we record:
- The **price of robustness**: $\text{cost}(\Gamma) - \text{cost}(0)$, the extra
  truck‑distance needed to immunize against $\Gamma$ worst‑case deviations.
- The optimal dual prices $z, w$ and which stations have $r_i > 0$ (resp. $t_j > 0$),
  i.e., whose deviation magnitude $\hat{s}_i$ exceeds the marginal value of budget —
  these are the stations that drive the protection cost.
- How does the routing plan change qualitatively — do trucks cluster near "safe"
  stations with small $\hat{s}_i$?

A key insight is that robust solutions should perform **better out‑of‑sample**
(Section 8.4): by hedging against demand fluctuations in June, they degrade less
when evaluated on July data. This is the classic bias–variance trade‑off.

---

## 6. Model Size & Solvability

| Model | Type | Variables | Solver | Est. solve time |
|---|---|---|---|---|
| Model 1 (Transportation LP) | LP | ~6,000 continuous | HiGHS | < 1 s |
| Model 2 (Routing MILP) | MILP | ~25,000 int + binary + continuous | HiGHS | ~10–60 s |
| Model 3 (Asymmetric MILP) | MILP | same as Model 2 | HiGHS | ~10–60 s |
| Model 4 (Robust MILP) | MILP | ~30,000 | HiGHS | ~30–120 s |

All models comfortably fit within HiGHS' capabilities on a modern laptop. The 157‑station
Jersey City network is large enough to be nontrivial but small enough for rapid
experimentation — ideal for a course project.

---

## 7. Implementation Plan (Julia + JuMP)

- **Data preprocessing:** a short Python script downloads the CitiBike zip, aggregates
  net flow per station ($\text{arrivals} - \text{departures}$), computes the Haversine
  distance matrix, and writes two small CSV files: `stations.csv` (name, lat, lon, net)
  and `distances.csv` ($157 \times 157$ matrix).
- **Julia modeling:** each model is implemented as a separate, well‑commented `JuMP`
  function that reads the same data files. We define a shared data structure
  (`struct BikeShareData`) so models differ only in objective and constraints.
- **Solver:** `HiGHS` (open‑source MILP) via `JuMP.set_optimizer`.
- **Visualization** with `Plots.jl`:
  - Jersey City map with stations colored by net flow (green = surplus, red = deficit)
  - Optimal bike movement arcs overlaid on the map
  - Pareto frontier curve (distance vs. unmet demand)
  - Robustness cost vs. $\Gamma$ (budget of uncertainty)

---

## 8. Planned Results & Discussion

### 8.1 In‑Sample Optimal Rebalancing Plan
Which stations are matched, how many bikes move along each arc, total truck‑distance.
The geographic map (Jersey City streets with colored station markers and flow‑width arcs)
is the centerpiece visualization.

### 8.2 LP vs. MILP Integrality Gap
Comparing Model 1's LP lower bound against Model 2's MILP optimum quantifies the cost
of indivisibility (truck routing, integer bikes). We also compare Model 1's LP relaxation
under symmetric vs. asymmetric distances (Model 3) to measure how asymmetry degrades the
LP's quality — does the integrality gap widen? Is it empirically worth it to use
asymmetric real‑road distances?

### 8.3 Multi‑Objective Pareto: Cost vs. Unmet Demand
*This is a **cross‑cutting analysis** applied to Model 2, not a separate model.* We
introduce a **soft‑deficit** variable $u_j \ge 0$ and optimize the weighted‑sum:

$$\min\; \sum_{k,i,j} d_{ij}\, z_{ijk} \;+\; \lambda \sum_{j\in D} u_j,$$

sweeping $\lambda$ from $0$ to a large value to trace the full **cost–coverage Pareto
frontier**. This answers a practical business question: if a budget only allows $T$
bike‑km of truck travel, what is the maximum achievable coverage? We discuss the
"diminishing returns" phenomenon — the last few stations are disproportionately expensive
to serve. This is analogous to the mean–variance efficient frontier in the Portfolio
Optimization sample project.

### 8.4 Out‑of‑Sample Validation
*Critical for demonstrating robustness of conclusions.* We:
- **Train** Models 2, 3, and 4 on June 2024 net flows.
- **Test** on July 2024 net flows: freeze the optimal rebalancing plan (which arcs are
  served, how many bikes) and evaluate its **realized cost** on July data.
- Compare **in‑sample optimal cost** (June) vs. **out‑of‑sample realized cost** (July)
  across all four models.
- Hypothesis: Model 4 (robust, $\Gamma>0$) has a higher in‑sample cost but the **smallest
  out‑of‑sample degradation**, because it hedges against demand fluctuations that
  materialize in the test month. This is the classic **bias–variance** trade‑off.

### 8.5 Sensitivity Analysis
- **Truck capacity $Q$:** how does the optimal plan change as trucks get larger?
- **Fleet size $\bar{K}$:** what is the marginal value of an additional truck?
- **Seasonal comparison:** compare June (summer) vs. January (winter) net flows to see
  how seasonal demand shifts affect the rebalancing plan.

### 8.6 Robust vs. Deterministic
How much extra distance must be budgeted to protect against $\Gamma$ stations deviating
from their mean net flow? Which stations are consistently "fragile" in the robust
solution? How does the **price of robustness** scale with problem size?

### 8.7 Limitations
Static single‑period model (no within‑day dynamics), simplified routing (no full TSP with
time windows), deterministic travel times. These are natural candidates for the "Future
Directions" section.

---

## 9. Conclusion & Possible Future Directions

We expect to deliver a comprehensive family of four bike rebalancing models — LP
transportation, capacitated MILP with MTZ routing, asymmetric‑distance sensitivity
analysis, and robust optimization (Bertsimas–Sim) — all solved on **real CitiBike data**
and visualized clearly. The model progression mirrors the structure of the highest‑graded
sample projects: each model introduces a **structural change** (integrality → routing
→ data‑driven sensitivity → uncertainty), not just a parameter variation.

Natural extensions include a **full Vehicle Routing Problem (VRP) formulation** with
explicit truck route sequencing (TSP constraints), a **dynamic (multi‑period) model**
that rebalances multiple times throughout the day, and a **two‑stage stochastic program**
that explicitly samples demand scenarios rather than using budgeted uncertainty.

---

## 10. Feasibility Checklist (vs. course guidelines)

| Requirement | How this project satisfies it |
|---|---|
| Modeling, not data collection | CitiBike data downloaded as‑is; 1‑page aggregation script; emphasis on 4 distinct formulations |
| Uses class techniques | LP, MILP (MTZ subtour elimination), total unimodularity theory, robust optimization (Bertsimas–Sim dual derivation), multi‑objective Pareto frontier |
| Toy → realistic progression | 5‑station hand‑verified transport → 157‑station Jersey City network with truck routing |
| Solvable in JuMP | `HiGHS` handles all models in seconds to low minutes |
| Avoids physics/deep learning, no algorithm tuning | Pure model formulation; off‑the‑shelf solver; MTZ chosen over DFJ precisely to avoid cutting‑plane/callback tuning |
| Rich results/discussion | LP vs MILP integrality gap + asymmetric vs symmetric comparison + Pareto frontier + out‑of‑sample validation + truck sensitivity + price of robustness vs Γ + seasonal comparison |
| Each model has a clear "why it's different" | TU→integral (M1), MILP with MTZ (M2), asymmetric sensitivity (M3), uncertainty → robust dual (M4) |

---

## References

[1] NABSA, *State of the Bike Share Industry Report*, North American Bikeshare &
Scootershare Association, 2024.

[2] F. L. Hitchcock, "The distribution of a product from several sources to numerous
localities," *Journal of Mathematics and Physics*, 20(1–4), pp. 224–230, 1941.

[3] T. Raviv, M. Tzur, and I. A. Forma, "Static repositioning in a bike‑sharing system:
models and solution approaches," *EURO Journal on Transportation and Logistics*, 2(3),
pp. 187–229, 2013.

[4] J. Brinkmann, M. W. Ulmer, and D. C. Mattfeld, "Short‑term strategies for
stochastic inventory routing in bike sharing systems," *Transportation Research
Procedia*, 10, pp. 293–302, 2015.

[5] D. Bertsimas and M. Sim, "The price of robustness," *Operations Research*, 52(1),
pp. 35–53, 2004.
