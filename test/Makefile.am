# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter           IBM    2006-04-13

AM_CPPFLAGS = -I$(srcdir)/../src -I$(srcdir)/../src/OsiCbc $(CBCUNITTEST_CFLAGS)

AM_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#              mip-sanity-data regression suite for Cbc                #
########################################################################
#
# Self-contained regression suite: solves every instance in the
# test/mip-sanity-data submodule (one Cbc thread per instance, suggested
# node/time limits from limits.tsv, up to nproc instances in parallel via
# GNU parallel) and validates each saved solution with cbc_validate_sol
# (feasibility, integrality, objective, optimality/infeasibility vs the
# best-known reference). See run-mip-sanity-tests for full documentation
# and options (`./run-mip-sanity-tests --help`). `make test` builds the
# validator and runs the whole suite; from cbc-workspace, `./test` is a
# convenience symlink to run-mip-sanity-tests.

test: cbc_validate_sol$(EXEEXT)
	$(srcdir)/run-mip-sanity-tests

.PHONY: test

########################################################################
#                    Manual C-interface solve utility                  #
########################################################################
#
# mip-c-solve: solves a single mip-sanity-data instance (by name, applying
# limits.tsv's suggested node/time limits) OR any arbitrary MPS/MPS.gz/LP
# file (given directly as `mip-c-solve <file> <timeLimit> <nodesLimit>`)
# in-process through Cbc's C interface (Cbc_C_Interface.h), then queries and
# validates the best solution and the whole solution pool (feasibility,
# objective consistency, and optimality/infeasibility claims vs bks.tsv when
# available). Unlike cbc_validate_sol, it does not rely on a `.sol` file
# produced by the `cbc` executable — it drives Cbc_solve() itself and leaves
# Cbc's own solve messages unsuppressed, so it doubles as a way to manually
# inspect the C interface end to end and confirm it behaves like the `cbc`
# command line for any instance. See `mip-c-solve --help` for usage.

noinst_PROGRAMS = cbc_validate_sol mip-c-solve mip-debug-cuts lp-relaxation-test \
	lazy-constraint-test mipstart-test coefstr-test rowred-test bkclique-bench \
	clqext-bench zerohalf-bench probing-bench

cbc_validate_sol_SOURCES = cbc_validate_sol.cpp

cbc_validate_sol_LDADD = ../src/libCbc.la

mip_c_solve_SOURCES = mip-c-solve.cpp

mip_c_solve_LDADD = ../src/libCbc.la

########################################################################
#           Row-cut-debugger diagnostic for wrong-answer bugs          #
########################################################################
#
# mip-debug-cuts: activates Osi's row-cut debugger (Cbc_activateRowCutDebugger)
# with a known-correct reference solution (e.g. a certified optimum from
# mip-sanity-data/sols) and solves the instance. Any cut, bound-fixing, or
# branching decision that would exclude the reference solution is flagged
# during the search -- useful for diagnosing WRONG_OBJ/INFEASIBLE_WRONG
# regressions. Requires a --debug (non-NDEBUG) Cbc build for the CoinAssert
# checks in CbcModel/CglXXX to actually raise. See `mip-debug-cuts --help`.

mip_debug_cuts_SOURCES = mip-debug-cuts.cpp

mip_debug_cuts_LDADD = ../src/libCbc.la

########################################################################
#              LP-relaxation solve/resolve tests (C interface)         #
########################################################################
#
# lp-relaxation-test: targeted tests for Cbc_solveLinearProgram()/
# Cbc_resolve() -- root LP relaxation solves under different LPMethod/
# DualPivot/LPReductions settings (including parallel LP racing and its
# >= 2 threads precondition), warm-start reoptimization (variable fixing and
# constraint add/remove), LP solution validation (bounds, reduced-cost
# complementary slackness), and the LP settings that must survive into a
# reoptimize -- notably INT_PARAM_PERT_VALUE, where Clp's default is
# "perturbation off" but Cbc's tuned default is the mild 50 that fast
# warm-started reoptimization relies on, and DBL_PARAM_PSI, where a Positive
# Edge pivot rule requested for the root solve must still be installed on
# every later resolve. Both are pure LP-performance regressions that no
# objective check would catch, since the LP still reaches the same optimum.
# Also covers the LP time limit: Clp keeps CPU and wall-clock deadlines as two
# independent parameters and checks both, so the limit must be armed on the one
# INT_PARAM_ELAPSED_TIME selects (it defaults to wall clock) and only that one,
# and "no limit" must disarm rather than arm a far-future deadline -- under LP
# racing, enforcing a wall-clock budget as CPU time would cut it short by
# roughly the thread count.
# Expected LP relaxation objectives were computed independently (two external
# solvers agreed exactly) and cross-checked against Cbc's own solve. Run as
# part of the regression suite by run-mip-sanity-tests (--no-lp-tests skips
# it). See `lp-relaxation-test --help`.

lp_relaxation_test_SOURCES = lp-relaxation-test.cpp

lp_relaxation_test_LDADD = ../src/libCbc.la

########################################################################
#               Lazy-constraint tests (C interface, TSP)               #
########################################################################
#
# lazy-constraint-test: targeted tests for the C interface's lazy-constraint
# machinery -- Cbc_addLazyConstraint() and Cbc_addCutCallback() with
# atSolution=1 -- driven by the TSP with subtour elimination, which is the
# canonical use of the feature (the SEC family is exponential, so it is stated
# lazily and only violated members ever reach a solver).
# The failure mode being guarded against is silent: a lazy constraint that is
# ignored still yields an optimal-looking answer, just to a relaxation. The
# assignment problem (degree constraints only) is exactly the TSP minus its
# SECs and is strictly weaker on every instance used, so the two are told apart
# by the objective alone. No fixture data is needed: instances come from a fixed
# LCG and each reference optimum is computed in-process three independent ways
# -- Held-Karp dynamic programming, the same MIP with SECs as ordinary rows, and
# the same MIP with them lazy -- which must agree while the SEC-free bound must
# come out strictly below. Published solutions are also walked as a successor
# array, since an objective check alone would pass on two subtours priced as if
# they were a tour.
# Covers enforcement, the LP relaxation staying untightened, repeated Cbc_solve()
# on one model (where the previous answer returns as a MIPStart, screened before
# the model is flagged as carrying lazy constraints), constraints registered
# between solves, Cbc_clone() deep-copying them, a SEC-violating MIPStart being
# rejected rather than trusted in both objective senses, the cut-callback and
# iterative-LP-cut-loop formulations, staged rows/columns reaching the solve at
# all, and the whole solution pool rather than just the best solution. Run as
# part of the regression suite by run-mip-sanity-tests (--no-lazy-tests skips
# it). See `lazy-constraint-test --help`.

lazy_constraint_test_SOURCES = lazy-constraint-test.cpp

lazy_constraint_test_LDADD = ../src/libCbc.la

########################################################################
#                  MIP-start tests (C interface)                       #
########################################################################
#
# mipstart-test: targeted tests for MIP starts through Cbc_setMIPStart() and
# Cbc_setMIPStartI(), which is the path python-mip uses. A MIP start is stated
# in the *original* problem while Cbc, by default, searches the *preprocessed*
# one, and callers typically supply only the nonzero integer values -- so
# CbcMipStart::computeCompleteSolution() has to fix what it recognizes, re-solve
# the LP to derive the rest, and fall back to a mini branch and bound when
# integers come back fractional. Each of those steps can fail quietly.
# The C interface is deliberately the vehicle rather than the -mipstart file:
# CbcMipStart::read() densifies before returning, building a vector over every
# column and overwriting the entries the file named, so a file listing three
# columns and a file listing all of them arrive identically and the input shapes
# a caller might send cannot be told apart on that path at all. Only the
# setMIPStart() path keeps them sparse.
# No fixture data is needed: the model is a small fixed-charge problem whose
# optimum is enumerated in-process over all open sets, with the continuous flows
# derived by an independent greedy fill, so a wrong answer is a wrong answer and
# not a stale expectation. Its continuous columns are undetermined by the
# integers until the LP runs, which is what makes an integers-only start a real
# test rather than a restatement of the answer.
# Covers the three input shapes -- nonzero integers only, the full solution, all
# nonzeros -- by name and by index, with preprocessing on and off, asserting for
# each that the start is accepted, the published point is integral and satisfies
# every row, the objective recomputed from the published columns matches the
# reported one, and columns the start never mentioned came back derived rather
# than left at zero. Also covers INT_PARAM_MIPSTART_FIX (the runtime mode that
# replaced the compile-time JUST_FIX_INTEGER macro, whose fix-continuous-too
# recovery path was unreachable as shipped), negative integer values -- where a
# positive-only filter on the solution auto-carried between solves would drop
# them silently -- an infeasible start being rejected rather than trusted, since
# accepting one installs a cutoff that can exclude the genuine optimum, starts
# naming columns the model does not have, and a suboptimal start being improved
# upon. Run as part of the regression suite by run-mip-sanity-tests
# (--no-mipstart-tests skips it). See `mipstart-test --help`.

mipstart_test_SOURCES = mipstart-test.cpp

mipstart_test_LDADD = ../src/libCbc.la

# coefstr-test: targeted tests for CbcCoefficientStrengthening, the coefficient
# tightening ("big-M strengthening") step of the pre-root-LP strengthening phase
# (step 3 of CbcSolver::preRootLPStrenghtening(), per -coefStrengthening).
#
# The pass shrinks an integer coefficient exceeding its row's slack down to that
# slack and adjusts the right-hand side to compensate, which tightens the LP
# relaxation while leaving the integer-feasible set alone. Both halves of that
# are asserted: the new coefficient, right-hand side and LP bound against values
# worked out by hand, and the feasible set by *enumerating* every integer point
# of the box before and after -- a re-solve reaching the same optimum would say
# nothing about the points the tightened system silently stopped admitting.
# Covers the textbook big-M row, a general integer column (the case Cgl's own
# binary-only `difference <= 1.0` gate skips, i.e. what this pass adds beyond the
# existing implementation), negative coefficients (right-hand side off the lower
# bound, not the upper), >= rows (the scale = -1 path), the eight rows and
# columns the rule must decline -- range and equality rows, infinite max
# activity, singletons, continuous, fixed and sub-unit-range columns,
# non-positive slack -- fractional slacks used as-is and near-integral ones
# snapped up rather than down, idempotence under a second pass, agreement between
# the row and column matrix views (a modifyCoefficient()-based implementation
# leaves the cached row copy stale), and the parameter end to end. Every fixture
# is small enough that Cbc reaches the right answer with the pass disabled, so
# each assertion was checked by mutation: all thirteen plausible wrong
# implementations tried are caught, including the three ways the deferred
# matrix copy (taken only once a change is known, to keep the no-change path
# off it) could fail to be applied. Run as part of the regression suite by
# run-mip-sanity-tests (--no-coefstr-tests skips it). See `coefstr-test --help`.

coefstr_test_SOURCES = coefstr-test.cpp

coefstr_test_LDADD = ../src/libCbc.la

# rowred-test: targeted tests for CbcRowReductions, the row-reduction step of the
# same pre-root-LP strengthening phase (step 4 of
# CbcSolver::preRootLPStrenghtening(), per -rowReductions).
#
# The pass deletes rows that cannot constrain the problem: rows all of whose
# columns are fixed, and rows that are duplicates or scalar multiples of another
# row -- the survivor inheriting the intersection of the two rows' bounds. Both
# halves are asserted, as in coefstr-test: the surviving rows and their merged
# bounds against values worked out by hand, and no integer-feasible point lost
# *or gained* by enumerating the whole box before and after. "Gained" matters as
# much as "lost" here, since dropping a row is only sound if the survivor still
# says everything the dropped row said, and a re-solve reaching the same optimum
# would show neither.
#
# The negative-scale cases carry most of the weight. Scaling a row by lambda < 0
# swaps its two bounds, and when one side is infinite the infinity has to swap
# with them rather than be scaled -- Cbc's +/-1e20 sentinel convention makes that
# an explicit test, unlike HiGHS's kHighsInf where the sign flip falls out of the
# arithmetic. lambda = -1 also exercises the hash's pivot tie-break: a row and
# its negation have equal-magnitude candidate pivots, so unless both resolve the
# tie the same way (largest magnitude, ties broken by smallest column index) they
# normalise to different vectors, never share a hash bucket and are never
# compared -- which makes every "1 row removed" assertion on a lambda = -1
# fixture a tie-break assertion too.
#
# Also covered: exact duplicates and duplicates whose indices arrive in a
# different order (the per-nonzero hashes are combined order-independently
# because CoinPackedMatrix makes no promise about index order within a vector), a
# scale factor of 1e6, a whole parallel class collapsing to one row,
# contradictory merged bounds (infeasible, and nothing deleted) versus bounds
# crossing by less than the feasibility tolerance (noise, snapped to an
# equation), duplicate and negated equations, independence from input row order,
# idempotence under a second pass, the counters adding up, and the eight cases
# the reduction must decline -- singletons, same support but not proportional,
# proportional only to 1e-7, disjoint and subset supports, partly-fixed rows, a
# model with nothing to find, and a model with no rows at all. The two end-to-end
# fixtures are picked so that the specific plausible bug moves the reported
# optimum: keeping the looser of two duplicate bounds reads 8 instead of 2, and
# forgetting the negative-lambda swap reads 6 instead of 8.
#
# Checked by mutation, as coefstr-test was: 18 of 26 plausible wrong
# implementations are caught, and the 8 that are not were each confirmed to be
# *equivalent* to the original rather than an assertion gap -- for every one, a
# sharpened variant that does break the underlying invariant is caught. Four are
# equivalent because OsiClp normalises each row to ascending column order, so the
# pivot tie-break, the smallest-vs-largest pivot choice, and the order-independent
# hash combination cannot be distinguished through this interface at all (they are
# there because CoinPackedMatrix promises no index order and the pass takes an
# OsiSolverInterface); two are re-checks the hash seed already implies (equal
# length, equal support); one is an early-out whose weighting is separately
# asserted; and one only stages a deletion that is never applied. The exercise
# was worth it for a different reason: it found a real bug, an all-fixed row
# violated by up to 1e-6 being deleted -- turning an infeasible model into a
# feasible-looking one -- because the removal test and the infeasibility test
# both used the feasibility tolerance and left no sliver between them.
# See `rowred-test --help`.

rowred_test_SOURCES = rowred-test.cpp

rowred_test_LDADD = ../src/libCbc.la

########################################################################
#          Clique separation / extension offline benchmarks            #
########################################################################
#
# bkclique-bench and clqext-bench: replay a single Bron-Kerbosch call captured
# from inside Cbc, so the BK algorithm can be worked on without paying for a
# full MIP solve per measurement. Both consume the fixtures written by
# CbcClqFixtureDump.hpp -- a Cbc built with -DCBC_DUMP_CLQSEP_FIXTURE emits, at
# each of the two BK call sites, the preprocessed problem, its optimal LP basis
# and solution, and the conflict graph serialized by
# CoinStaticConflictGraph::save().
#
# The graph is reloaded rather than rebuilt because rebuilding does not
# reproduce it: Cbc's graph is built before CglCliqueStrengthening deletes and
# adds rows, and OsiSolverInterface::checkCGraph() only rebuilds when the
# *column* count changes, so the graph a separator sees is routinely stale with
# respect to the matrix. That staleness is the behaviour under study, hence
# OsiSolverInterface::setCGraph() and the serialization. --rebuild-cgraph keeps
# the rebuild path available for direct comparison, and `bkclique-bench
# --self-test` asserts the round trip is byte-exact and field-for-field
# identical, degrees included (they are deliberately approximate, so recomputing
# them would give different values).
#
# bkclique-bench drives CglBKClique over N rounds with an LP re-solve between
# them, and reports the sum of cut violations as well as the cut count: a
# separator finding fewer but more violated cliques is doing better work, which
# a count alone cannot show. It builds a fresh generator each round, because
# CglBKClique::separateCliques() ends with `maxCallsBK_ = bk->numCalls()` -- a
# monotonic ratchet that otherwise shrinks the budget of every later round and
# makes the rounds incomparable.
#
# clqext-bench drives CglCliqueStrengthening instead, reporting rows extended
# and dominated, the change in matrix size, and the LP bound before and after.
# It refuses to run when the warm start misses optimality: getReducedCost()
# returns NULL unless the model is proven optimal, and extension methods 4 and 5
# then silently downgrade to 2.
#
# Neither is part of the regression suite; both are measurement tools.

bkclique_bench_SOURCES = bkclique-bench.cpp

bkclique_bench_LDADD = ../src/libCbc.la

clqext_bench_SOURCES = clqext-bench.cpp

clqext_bench_LDADD = ../src/libCbc.la

########################################################################
#        0-1/2 separation replay benchmark (measurement tool)          #
########################################################################
#
# zerohalf-bench replays one CglZeroHalf separation call captured from Cbc by
# CbcZeroHalfFixtureDump.hpp, so an experiment on the 0-1/2 algorithm costs
# milliseconds instead of a preprocess plus a root LP.
#
# It loads no auxiliary structure, and that is faithful rather than a shortcut:
# CglZeroHalf::refreshSolver() rebuilds every array it uses from the solver on
# each call, so reconstructing from the model is exactly what Cbc does. Contrast
# bkclique-bench, which must load a serialized conflict graph because Cbc caches
# a stale one.
#
# --max-pairs=0 is the control flag: it disables the O(cnt^3) pair-weakening
# stage, so the difference against the default prices that stage exactly.
#
# Not part of the regression suite; it is a measurement tool.

zerohalf_bench_SOURCES = zerohalf-bench.cpp

zerohalf_bench_LDADD = ../src/libCbc.la

########################################################################
#          probing replay benchmark (measurement tool)                 #
########################################################################
#
# probing-bench replays one CglProbing call captured from Cbc by
# CbcProbingFixtureDump.hpp. Same motive as the two benches above -- reaching the
# root cut pass costs a preprocess plus an LP -- but three things differ enough to
# be worth stating, because each one is a way to accidentally measure a call Cbc
# never makes.
#
# It reports row cuts and column cuts SEPARATELY. Probing is the only generator
# here that returns both, and the column cuts (bound tightenings) often carry the
# bound movement without adding a row at all. A row-cut count alone would score
# air05 -- 0.75s, 0 row cuts, 62 column cuts -- as pure waste.
#
# It calls generateCutsAndModify with a CglTreeProbingInfo, not generateCuts with
# a plain CglTreeInfo, because CbcCutGenerator.cpp:357-376 special-cases Probing
# that way. CglTreeProbingInfo::initializeFixing returns nonzero, which switches on
# the implication-saving path inside probe(); the plain call skips it and is
# therefore cheaper than Cbc's. --plain-info prices that difference rather than
# inheriting it silently.
#
# Its default --pass=0 is the expensive call, not merely the first: probe() sets
# justFix = -1 when pass == 0 at the root, which overrides the maxProbe==123
# candidate subsetting so every candidate column is probed.
#
# Control flag: --max-probe=0 disables probing outright (the gate is real, at
# CglProbing.cpp:2350), leaving tighten() and the column cuts it implies, so
# time(123) - time(0) prices probing and time(0) is the floor. Note --max-look is
# NOT a control flag despite appearances -- probe() overwrites maxStack with a
# hardcoded 80 -- and a positive --max-probe only permutes the candidate order
# rather than capping it. Both are documented at length in the source.
#
# Not part of the regression suite; it is a measurement tool.

probing_bench_SOURCES = probing-bench.cpp

probing_bench_LDADD = ../src/libCbc.la
