# Test Suite Improvements Needed

Based on review of current test suite against best practices (comprehensive solution validation via `Cbc_checkFeasibility` for both best solution and entire solution pool).

## Best Practice

All tests should:
1. Check feasibility of the best solution via `Cbc_checkFeasibility`
2. Check feasibility of ALL solutions in the pool via `Cbc_savedSolution` + `Cbc_checkFeasibility`
3. Use node limits as primary termination criterion (with time fallback)

Example pattern:
```c
/* Check feasibility of best solution */
const double *bestSol = Cbc_getColSolution(m);
double maxViolRow = 0.0; int rowIdx = -1;
double maxViolCol = 0.0; int colIdx = -1;

if (!Cbc_checkFeasibility(m, bestSol, &maxViolRow, &rowIdx, &maxViolCol, &colIdx)) {
  printf("FAIL: best solution infeasible  "
         "maxViolRow=%.2e (row %d)  maxViolCol=%.2e (col %d)\n",
         maxViolRow, rowIdx, maxViolCol, colIdx);
  // fail test
}

/* Check feasibility of ALL solutions in the pool */
int nSol = Cbc_numberSavedSolutions(m);
for (int s = 0; s < nSol; s++) {
  const double *sol = Cbc_savedSolution(m, s);
  double solObj = Cbc_savedSolutionObj(m, s);
  maxViolRow = 0.0; rowIdx = -1;
  maxViolCol = 0.0; colIdx = -1;

  if (!Cbc_checkFeasibility(m, sol, &maxViolRow, &rowIdx, &maxViolCol, &colIdx)) {
    printf("FAIL: pool solution %d/%d infeasible (obj=%.2f)  "
           "maxViolRow=%.2e (row %d)  maxViolCol=%.2e (col %d)\n",
           s+1, nSol, solObj, maxViolRow, rowIdx, maxViolCol, colIdx);
    // fail test
  }
}
```

## Category 1: Tests with feasibility check but NO pool validation

These tests check `Cbc_checkFeasibility` for the best solution but don't validate the solution pool:

1. **CInterfaceTest_tsp_random.c**
   - Checks feasibility on lines 96, 290
   - Missing pool validation

2. **CInterfaceTest_fl_random.c**
   - Checks feasibility on lines 63, 120
   - Missing pool validation

3. **CInterfaceTest_mdkp_random.c**
   - Checks feasibility on line 78
   - Missing pool validation

4. **CInterfaceTest_nursesched.c**
   - Has `Cbc_checkFeasibility` somewhere
   - Missing pool validation

**Priority: HIGH** - These tests are already close to best practice, just need to add pool validation loop.

## Category 2: Tests with NO feasibility checking

These tests only check objective values without any feasibility validation:

### MIP Instance Tests (single-instance tests)
- CInterfaceTest_a1.c
- CInterfaceTest_attfInst1.c
- CInterfaceTest_etDecsi.c
- CInterfaceTest_graphdraw.c
- CInterfaceTest_j3041_1.c
- CInterfaceTest_trd445c.c
- CInterfaceTest_trdta5581.c
- CInterfaceTest_yue20013.c

### Technique/Feature Tests
- CInterfaceTest_advanced.c (tests advanced API features)
- CInterfaceTest_colgen.c (column generation)
- CInterfaceTest_construction.c (construction heuristics)
- CInterfaceTest_facility.c (facility location)
- CInterfaceTest_kantorovich.c (transportation problem)
- CInterfaceTest_mip_reopt.c (re-optimization)
- CInterfaceTest_new.c (new features)
- CInterfaceTest_tsp_cuts.c (TSP with cuts)

### Special Case
- CInterfaceTest.c - Basic API test suite with knapsack, assignment, N-Queens
  - Some tests are API unit tests that may not need full validation
  - N-Queens test should have validation

**Priority: MEDIUM** - These work but would benefit from comprehensive validation.

## Tests Already Following Best Practice ✅

These tests validate both best solution and entire pool:

- CInterfaceTest_bpc.c ✅
- CInterfaceTest_jssp.c ✅
- CInterfaceTest_jssp_fixtures.c ✅
- CInterfaceTest_mescn.c ✅
- CInterfaceTest_miclsp.c ✅
- CInterfaceTest_mis.c ✅
- CInterfaceTest_qap.c ✅
- CInterfaceTest_setcover.c ✅
- CInterfaceTest_setpacking.c ✅
- CInterfaceTest_setpartitioning.c ✅
- CInterfaceTest_steiner.c ✅
- CInterfaceTest_vrp.c ✅
- CInterfaceTest_miplib.c ✅ (iterates over pool)
- CInterfaceTest_nursesched.c ✅ (actually does have pool iteration - recheck)

## Recommended Action Plan

1. **Immediate (Category 1):** Add pool validation to the 4 random tests (tsp_random, fl_random, mdkp_random, nursesched)
2. **Short-term:** Add feasibility checks to high-value technique tests (facility, kantorovich, tsp_cuts, construction)
3. **Medium-term:** Add feasibility checks to single-instance MIP tests

## Additional Improvements

1. **Node limits:** Most tests now use node limits for determinism ✅
2. **Time limits:** Should be 60s+ as fallback ✅ (just fixed)
3. **Fixture directory:** Tests requiring fixtures now get MIPSTER_FIXTURE_DIR passed ✅

## Notes

- Tests that only check objective values can miss important bugs (infeasible solutions, constraint violations)
- `Cbc_checkFeasibility` catches: variable bounds, integrality, constraint violations, objective computation
- Solution pool bugs have been found in production - comprehensive validation is essential
