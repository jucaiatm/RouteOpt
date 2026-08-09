# VURP-M Branch-Price-and-Cut application

This directory contains an exact branch-price-and-cut implementation for the
**vessel--UAV routing problem with multiple visits in one UAV flight (VURP-M)**
under the total-completion-time objective.

The implementation follows the decomposition used in the accompanying VURP-M
MISOCP/LBBD manuscript:

- one vessel and one UAV;
- a UAV sortie visits one or more platforms;
- launch and recovery positions are continuous decisions in a rectangular
  feasible region `OMEGA`;
- the vessel may move while the UAV is airborne and may wait at recovery;
- UAV endurance is reset between sorties;
- the objective is the elapsed time from vessel departure until the vessel and
  UAV return to the depot.

The new application is isolated under `packages/application/vurp_m`; the
existing CVRP and VRPTW applications are not modified.

## Column definition

A generated column is one complete **ordered UAV sortie path**

```text
(i_1, i_2, ..., i_m)
```

and stores its visited set, first platform, last platform, and exact internal
UAV distance

```math
h_r = \sum_{j=1}^{m-1} d_{i_j,i_{j+1}}.
```

The restricted master contains a copy `lambda[k,r]` for sortie slot `k`. The
column is linked to the master variables for platform assignment, first/last
platform selection, and internal route length. Continuous launch/recovery
coordinates remain in the master and are not discretized.

## Exact algorithm

The executable implements the following loop at every branch-and-bound node:

1. Phase-I column generation with artificial variables.
2. Phase-II column generation for the total-time objective.
3. Exact subset dynamic-programming pricing for every sortie slot.
4. Dynamic outer approximation of all second-order-cone constraints.
5. Branching on `z[k]`, `y[i,k]`, `first[i,k]`, and `last[i,k]`.
6. An exact fixed-partition SOCP oracle using Gurobi convex QCP constraints.
7. Logic-based Benders value cuts or no-good cuts for integer partitions.

The internal-length linking row is oriented as

```math
\ell_k - \sum_r h_r\lambda_{kr} \ge 0.
```

This inequality is equivalent to equality at an optimum: increasing
`ell[k]` can only tighten the UAV-duration constraint and cannot improve the
minimization objective. The greater-than row has a nonnegative dual in the
minimization master. Consequently, for a fixed visited set, first platform,
and last platform, the shortest Hamilton path is the only path that can
minimize reduced cost. This establishes the correctness of the
Held--Karp-style exact pricing routine.

The static incompatibility rule from the VURP-M model is applied during both
master construction and pricing:

```math
2d_{ij} > (v_d+v_v)B
\quad\Longrightarrow\quad
\text{platforms }i,j\text{ cannot share a sortie}.
```

## Exactness scope

The solver is exact, up to the configured numerical tolerances, for the model
implemented here:

- `OMEGA` is an axis-aligned rectangle;
- platform service times and UAV payload constraints are absent, matching the
  referenced total-time formulation;
- exact subset pricing is used for every processed branch-and-bound node;
- an integer solution is accepted only after the fixed-partition convex SOCP
  has certified it.

Exact subset pricing is enabled up to 20 platforms by default. This is a
memory guard, not a relaxation. Raising `--max-exact-customers` preserves
exactness but requires exponentially more memory and time.

## Instance format

```text
NAME Example5
DIMENSION 5
VESSEL_SPEED 1.0
UAV_SPEED 3.0
ENDURANCE 6.0
DEPOT 0.0 0.0
OMEGA -2.0 12.0 -2.0 12.0
PLATFORMS
1 2.0 1.0
2 4.0 3.0
3 7.0 2.0
4 8.0 7.0
5 3.0 8.0
EOF
```

Platform identifiers must be consecutive integers `1,...,n`.

## Build

Configure the repository's `FindGUROBI.cmake` first. The root `build.py`
already performs this update. Then build the application independently:

```bash
cmake -S packages/application/vurp_m \
      -B packages/application/vurp_m/build \
      -DCMAKE_BUILD_TYPE=Release
cmake --build packages/application/vurp_m/build -j
```

The executable is written to:

```text
packages/application/vurp_m/bin/vurp_m_bpc
```

A convenience wrapper is also provided:

```bash
python3 packages/application/vurp_m/build.py
```

## Run

```bash
packages/application/vurp_m/bin/vurp_m_bpc \
    packages/application/vurp_m/instances/example5.txt \
    --time-limit 7200 \
    --threads 1
```

Useful options:

```text
--pricing-columns <count>
--oa-cuts <count>
--max-exact-customers <count>
--verbose-lp
--verbose-socp
```

## Validation protocol

Before using larger instances, compare small instances against the complete
MISOCP or the existing LBBD implementation. At minimum, verify:

1. identical optimal total completion time;
2. identical feasibility of every reported sortie partition;
3. equality between the final RMP reduced-cost calculation and explicit
   `c - pi^T a` recomputation;
4. no positive artificial variable when Phase I terminates;
5. all SOC residuals are within tolerance at the final solution.
