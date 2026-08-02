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
	lazy-constraint-test mipstart-test

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
