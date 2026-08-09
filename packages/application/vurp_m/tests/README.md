# VURP-M regression tests

Run the solver-independent exact-pricing regression with:

```bash
python3 packages/application/vurp_m/tests/run_pricing_regression.py
```

The test compiles the production `ExactPricing` and Held--Karp implementation
against a minimal Gurobi C API link stub. It then compares the dynamic-programming
pricing value against exhaustive enumeration on a five-platform instance,
including branch-induced required/forbidden/first/last restrictions. It also
compares the fixed-endpoint Hamilton-path value against permutation enumeration.

This regression does **not** solve the RMP or the fixed-partition SOCP. A real
Gurobi installation and license are required for the production executable and
for end-to-end comparison against the MISOCP/LBBD implementation.
