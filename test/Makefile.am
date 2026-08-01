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

noinst_PROGRAMS = cbc_validate_sol mip-c-solve mip-debug-cuts lp-relaxation-test

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
# warm-started reoptimization relies on, so a path that fails to forward it
# is a pure LP-performance regression no objective check would catch.
# Expected LP relaxation objectives were computed independently (two external
# solvers agreed exactly) and cross-checked against Cbc's own solve. Run as
# part of the regression suite by run-mip-sanity-tests (--no-lp-tests skips
# it). See `lp-relaxation-test --help`.

lp_relaxation_test_SOURCES = lp-relaxation-test.cpp

lp_relaxation_test_LDADD = ../src/libCbc.la
