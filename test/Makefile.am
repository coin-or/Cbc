# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter           IBM    2006-04-13

AM_CPPFLAGS = -I$(srcdir)/../src -I$(srcdir)/../src/CoinUtils -I$(srcdir)/../src/Clp -I$(srcdir)/../src/Cgl $(CBCUNITTEST_CFLAGS)

AM_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#                   Install paths for fixtures / runner / desktop      #
########################################################################

# Fixtures install to $(datadir)/mipster/test/fixtures/
pkgtestfixturedir = $(datadir)/mipster/test/fixtures
dist_pkgtestfixture_DATA = \
  fixtures/A-1.mps.gz \
  fixtures/graphdraw-domain.mps.gz \
  fixtures/nursesched-sprint02.mps.gz \
  fixtures/trd445c.mps.gz \
  fixtures/trdta5581.mps.gz

# Runner script installs to $(bindir)/run-mipster-tests
dist_bin_SCRIPTS = run-mipster-tests

# Desktop entry installs to $(datadir)/applications/
appsdir = $(datadir)/applications
dist_apps_DATA = mipster-tests.desktop mipster-tests-debug.desktop

########################################################################
#                         unitTest for Cbc                             #
########################################################################

CBC_TEST_TGTS = test_cbc ctests

cbcflags =
if COIN_HAS_SAMPLE
#  cbcflags += -dirSample `$(CYGPATH_W) $(SAMPLE_DATA)`
endif
if COIN_HAS_MIPLIB3
  cbcflags += -dirMiplib `$(CYGPATH_W) $(MIPLIB3_DATA)` -unitTest
#  cbcflags += -dirMiplib `$(CYGPATH_W) $(MIPLIB3_DATA)` -miplib
endif

# We are using the solver executables to do the unit test

test: $(CBC_TEST_TGTS)

test_cbc: ../src/cbc$(EXEEXT)
	RUNNING_TEST="test_cbc" ../src/cbc$(EXEEXT) $(cbcflags)

if COIN_HAS_SAMPLE
test_cbc_generic: ../src/cbc-generic$(EXEEXT)
	@echo "A minimal test is all that cbc-generic can manage just yet."
	RUNNING_TEST="cbc-generic" ../src/cbc-generic$(EXEEXT) -directory `$(CYGPATH_W) $(SAMPLE_DATA)` -import p0033.mps -branch
endif


.PHONY: test

bin_PROGRAMS = CInterfaceTest CInterfaceTest_tsp_random CInterfaceTest_fl_random CInterfaceTest_mdkp_random CInterfaceTest_nursesched CInterfaceTest_a1 CInterfaceTest_graphdraw CInterfaceTest_trdta5581 CInterfaceTest_trd445c CbcSolverLpTest

########################################################################
#                         unitTest for CbcCInterface                   #
########################################################################

CInterfaceTest_SOURCES = CInterfaceTest.c
nodist_EXTRA_CInterfaceTest_SOURCES = dummy.cpp # force using C++ linker

CInterfaceTest_LDADD = ../src/libmipster.la

########################################################################
#                   TSP random-instance tests                           #
########################################################################

CInterfaceTest_tsp_random_SOURCES = CInterfaceTest_tsp_random.c
nodist_EXTRA_CInterfaceTest_tsp_random_SOURCES = dummy.cpp
CInterfaceTest_tsp_random_LDADD = ../src/libmipster.la

########################################################################
#                   FL random-instance tests (UFL + CFLP)              #
########################################################################

CInterfaceTest_fl_random_SOURCES = CInterfaceTest_fl_random.c
nodist_EXTRA_CInterfaceTest_fl_random_SOURCES = dummy.cpp
CInterfaceTest_fl_random_LDADD = ../src/libmipster.la

########################################################################
#                   MDKP random-instance tests                          #
########################################################################

CInterfaceTest_mdkp_random_SOURCES = CInterfaceTest_mdkp_random.c
nodist_EXTRA_CInterfaceTest_mdkp_random_SOURCES = dummy.cpp
CInterfaceTest_mdkp_random_LDADD = ../src/libmipster.la

########################################################################
#                   CbcSolver LP relaxation C++ tests                  #
########################################################################

CbcSolverLpTest_SOURCES = CbcSolverLpTest.cpp
CbcSolverLpTest_LDADD = ../src/libmipster.la

########################################################################
#                   Nursesched MPS-file tests (LP+BP, MIP)             #
########################################################################

# FIXTURE_DIR points to the directory holding the MPS fixture files so the
# compiled binary can locate them at runtime regardless of the working directory
# from which make runs.
FIXTURE_DIR_VAL = $(abs_srcdir)

CInterfaceTest_nursesched_SOURCES = CInterfaceTest_nursesched.c
nodist_EXTRA_CInterfaceTest_nursesched_SOURCES = dummy.cpp
CInterfaceTest_nursesched_CPPFLAGS = $(AM_CPPFLAGS) -DFIXTURE_DIR=\"$(FIXTURE_DIR_VAL)\"
CInterfaceTest_nursesched_LDADD = ../src/libmipster.la

########################################################################
#                   A-1 project-scheduling tests (LP reopt, MIP)       #
########################################################################

CInterfaceTest_a1_SOURCES = CInterfaceTest_a1.c
nodist_EXTRA_CInterfaceTest_a1_SOURCES = dummy.cpp
CInterfaceTest_a1_CPPFLAGS = $(AM_CPPFLAGS) -DFIXTURE_DIR=\"$(FIXTURE_DIR_VAL)\"
CInterfaceTest_a1_LDADD = ../src/libmipster.la

########################################################################
#                   graphdraw-domain tests (LP + node-limited MIP)     #
########################################################################

CInterfaceTest_graphdraw_SOURCES = CInterfaceTest_graphdraw.c
nodist_EXTRA_CInterfaceTest_graphdraw_SOURCES = dummy.cpp
CInterfaceTest_graphdraw_CPPFLAGS = $(AM_CPPFLAGS) -DFIXTURE_DIR=\"$(FIXTURE_DIR_VAL)\"
CInterfaceTest_graphdraw_LDADD = ../src/libmipster.la

########################################################################
#                   trdta5581 tests (integer infeasibility via BP)     #
########################################################################

CInterfaceTest_trdta5581_SOURCES = CInterfaceTest_trdta5581.c
nodist_EXTRA_CInterfaceTest_trdta5581_SOURCES = dummy.cpp
CInterfaceTest_trdta5581_CPPFLAGS = $(AM_CPPFLAGS) -DFIXTURE_DIR=\"$(FIXTURE_DIR_VAL)\"
CInterfaceTest_trdta5581_LDADD = ../src/libmipster.la

########################################################################
#                   trd445c tests (LP reopt, MIP optimal)             #
########################################################################

CInterfaceTest_trd445c_SOURCES = CInterfaceTest_trd445c.c
nodist_EXTRA_CInterfaceTest_trd445c_SOURCES = dummy.cpp
CInterfaceTest_trd445c_CPPFLAGS = $(AM_CPPFLAGS) -DFIXTURE_DIR=\"$(FIXTURE_DIR_VAL)\"
CInterfaceTest_trd445c_LDADD = ../src/libmipster.la

ctests: CInterfaceTest$(EXEEXT) CInterfaceTest_tsp_random$(EXEEXT) CInterfaceTest_fl_random$(EXEEXT) CInterfaceTest_mdkp_random$(EXEEXT) CInterfaceTest_nursesched$(EXEEXT) CInterfaceTest_a1$(EXEEXT) CInterfaceTest_graphdraw$(EXEEXT) CInterfaceTest_trdta5581$(EXEEXT) CInterfaceTest_trd445c$(EXEEXT) CbcSolverLpTest$(EXEEXT)
	MIPSTER_FIXTURE_DIR=$(abs_srcdir)/fixtures $(srcdir)/run-mipster-tests
