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
# mip-c-solve: solves a single mip-sanity-data instance in-process through
# Cbc's C interface (Cbc_C_Interface.h), applying limits.tsv's suggested
# node/time limits, then queries and validates the best solution and the
# whole solution pool (feasibility, objective consistency, and
# optimality/infeasibility claims vs bks.tsv). Unlike cbc_validate_sol, it
# does not rely on a `.sol` file produced by the `cbc` executable — it drives
# Cbc_solve() itself and leaves Cbc's own solve messages unsuppressed, so it
# doubles as a way to manually inspect the C interface end to end. See
# `mip-c-solve --help` for usage.

noinst_PROGRAMS = cbc_validate_sol mip-c-solve lp-relaxation-test

cbc_validate_sol_SOURCES = cbc_validate_sol.cpp

cbc_validate_sol_LDADD = ../src/libCbc.la

mip_c_solve_SOURCES = mip-c-solve.cpp

mip_c_solve_LDADD = ../src/libCbc.la

########################################################################
#              LP-relaxation solve/resolve tests (C interface)         #
########################################################################
#
# lp-relaxation-test: targeted tests for Cbc_solveLinearProgram()/
# Cbc_resolve() -- root LP relaxation solves under different LPMethod/
# DualPivot/LPReductions settings, warm-start reoptimization (variable
# fixing and constraint add/remove), and LP solution validation (bounds,
# reduced-cost complementary slackness). Expected LP relaxation objectives
# were computed independently (two external solvers agreed exactly) and
# cross-checked against Cbc's own solve. See `lp-relaxation-test --help`.

lp_relaxation_test_SOURCES = lp-relaxation-test.cpp

lp_relaxation_test_LDADD = ../src/libCbc.la
