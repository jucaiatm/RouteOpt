# Example9_1 P0 performance hotfix

## Why this hotfix exists

The 7200-second Example9_1 run reported the baseline executable fields
`Status`, `Patterns`, and `Slot columns`. That executable is not the later
modular performance build. The baseline implementation also closes the SOC
outer approximation at every fractional node and globally retains all generated
tangents, which explains the observed OA-row explosion.

The patch in this directory applies to the modular performance source delivered
as `RouteOpt-vurpm-bpc-performance-source`. It implements the following P0
changes:

1. fractional-node OA tangents are local to the current branch;
2. root OA separation is capped at 6 rounds;
3. non-root OA separation is capped at 1 round;
4. at most 24 OA cuts are added per round by default;
5. full fractional SOC closure is no longer required before branching;
6. fixed sortie counts are processed from the incumbent's sortie count outward;
7. the initial heuristic reserves at least one SOCP candidate per sortie count;
8. near-parallel OA directions use a less pathological duplicate threshold.

Exactness is preserved: every node master remains a valid relaxation, and every
integer partition is evaluated by the exact fixed-partition SOCP oracle before
it can update the incumbent. A feasible SOCP result generates a value cut; an
infeasible result generates a no-good cut.

## Apply

From a checkout that already contains the modular performance source:

```bash
git apply packages/application/vurp_m/hotfixes/Example9_1_P0_hotfix.patch
```

New controls:

```text
--oa-cuts N
--root-oa-rounds N
--node-oa-rounds N
--ascending-sortie-order
```

Recommended first rerun:

```bash
vurp_m_bpc Example9_1.txt \
  --time-limit 7200 \
  --threads 1 \
  --pricing-threads 1 \
  --oa-cuts 24 \
  --root-oa-rounds 6 \
  --node-oa-rounds 1 \
  --solution Example9_1_p0.json
```

## Validation already performed

- C++20 syntax check of all application sources with GCC 14 and
  `-Wall -Wextra -Wpedantic -Werror`;
- core test passed;
- exactness test passed, including Held--Karp, enumeration pricing, exact
  labeling, heuristic pricing, and subset-row reduced-cost checks;
- performance primitive test passed.

A licensed Gurobi end-to-end run was not available in the patch-generation
environment, so no runtime claim is made before the Example9_1 rerun.
