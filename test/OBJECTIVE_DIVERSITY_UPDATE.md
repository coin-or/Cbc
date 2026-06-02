# Objective Coefficient Diversity Update

## Overview

Updated test fixture generators to include **diverse objective coefficient types** across instances. This tests solver robustness with different numerical ranges and precision requirements.

## Coefficient Types Added

| Type | Range | Format | Purpose |
|------|-------|--------|---------|
| **integer** | [1, 100] | Integer | Baseline (original behavior) |
| **continuous** | [1.0, 100.0] | Float (2 decimals) | Test continuous coefficients |
| **large** | [100, 1000] | Integer | Test scaling with large values |
| **fractional** | [0.1, 10.0] | Float (2 decimals) | Test small fractional values |

## Generators Updated

### 1. Set Covering (`gen_setcover_fixtures.py`)

**Variable:** `costs` (set costs)

**Distribution:**
- 2 instances with `integer` costs
- 2 instances with `continuous` costs
- 1 instance with `large` costs
- 1 instance with `fractional` costs

### 2. Set Packing (`gen_setpacking_fixtures.py`)

**Variable:** `values` (set values/profits)

**Distribution:**
- 3 instances with `integer` values
- 2 instances with `continuous` values
- 1 instance with `large` values
- 1 instance with `fractional` values

## Generators NOT Updated (Rationale)

| Generator | Reason |
|-----------|--------|
| **QAP** | Already complex (quadratic terms, McCormick linearization); integer costs sufficient |
| **Set Partitioning** | Few instances (4 total); keep simple |
| **BPC** | Item weights are structural (capacity constraints); changing them affects feasibility |
| **MIS** | Vertex weights can stay integer; graph structure is the main diversity dimension |

## Impact

### Test Coverage

- **Numerical robustness:** Tests solver with different coefficient ranges
- **Precision handling:** Tests rounding/tolerance with fractional values
- **Scaling:** Tests performance with large coefficients (1000x baseline)

### Examples

**Set Covering - Fractional costs:**
```
scp_e20_s25_sd137_hierarchical (costs=fractional)
  Costs in [0.1, 10.0], e.g., {2.37, 0.45, 9.12, ...}
  Tests solver precision on small fractional objectives
```

**Set Packing - Large values:**
```
spp_e30_s50_sd42_random (values=large)
  Values in [100, 1000], e.g., {523, 847, 192, ...}
  Tests scaling and potential overflow handling
```

**Set Covering - Continuous costs:**
```
scp_e15_s20_sd42_sparse (costs=continuous)
  Costs in [1.0, 100.0], e.g., {45.23, 12.67, 88.91, ...}
  Tests continuous coefficient handling
```

## Verification

All instances solved with external solver to obtain certified optimal values **before** objective coefficient diversity was added to ensure:
1. Feasibility unchanged
2. Optimal solutions exist
3. Problem structure preserved

After regeneration with diverse coefficients:
1. Optimal objective values will differ
2. Problem structure remains the same
3. Solution methods (clique cuts, etc.) still applicable

## Regenerating Fixtures

To regenerate with new diverse objectives:

```bash
cd test

# Regenerate with diverse coefficients
./gen_setcover_fixtures.py    # 7 instances, 4 different cost types
./gen_setpacking_fixtures.py  # 8 instances, 4 different value types

# Output shows cost type per instance:
# Generating scp_e15_s20_sd42_sparse (costs=continuous)...
# Generating spp_e30_s50_sd42_random (values=large)...
```

## Testing Notes

When testing these fixtures:
- **Tolerance:** Use relative error `|actual - expected| / |expected|` instead of absolute
- **Fractional instances:** May need tighter tolerances (1e-4 instead of 1e-3)
- **Large instances:** Check for numerical stability (no overflow in objective accumulation)

## Future Extensions

Potential additional diversity dimensions:
1. **Mixed coefficients:** Some integer, some continuous in same instance
2. **Negative coefficients:** For maximization problems (already handled via negation)
3. **Very small coefficients:** [0.001, 0.1] to test underflow
4. **Scientific notation range:** [1e-3, 1e3]

For now, the 4 types above provide good coverage without excessive combinations.
