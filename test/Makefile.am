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
  fixtures/etDecsi.mps.gz \
  fixtures/graphdraw-domain.mps.gz \
  fixtures/j3041_1.mps.gz \
  fixtures/nursesched-sprint02.mps.gz \
  fixtures/trd445c.mps.gz \
  fixtures/trdta5581.mps.gz \
  fixtures/yue20013.1.150.mps.gz

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

bin_PROGRAMS = CInterfaceTest CInterfaceTest_tsp_random CInterfaceTest_fl_random CInterfaceTest_mdkp_random CInterfaceTest_nursesched CInterfaceTest_a1 CInterfaceTest_graphdraw CInterfaceTest_trdta5581 CInterfaceTest_trd445c CInterfaceTest_j3041_1 CInterfaceTest_yue20013 CInterfaceTest_etDecsi CbcSolverLpTest CInterfaceTest_advanced CInterfaceTest_colgen CInterfaceTest_construction CInterfaceTest_facility CInterfaceTest_jssp CInterfaceTest_kantorovich CInterfaceTest_mip_reopt CInterfaceTest_new CInterfaceTest_tsp_cuts

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

########################################################################
#                   j3041_1 tests (LP reopt, MIP optimal)             #
########################################################################

CInterfaceTest_j3041_1_SOURCES = CInterfaceTest_j3041_1.c
nodist_EXTRA_CInterfaceTest_j3041_1_SOURCES = dummy.cpp
CInterfaceTest_j3041_1_CPPFLAGS = $(AM_CPPFLAGS) -DFIXTURE_DIR=\"$(FIXTURE_DIR_VAL)\"
CInterfaceTest_j3041_1_LDADD = ../src/libmipster.la

########################################################################
#                   yue20013.1.150 tests (LP reopt, MIP optimal)      #
########################################################################

CInterfaceTest_yue20013_SOURCES = CInterfaceTest_yue20013.c
nodist_EXTRA_CInterfaceTest_yue20013_SOURCES = dummy.cpp
CInterfaceTest_yue20013_CPPFLAGS = $(AM_CPPFLAGS) -DFIXTURE_DIR=\"$(FIXTURE_DIR_VAL)\"
CInterfaceTest_yue20013_LDADD = ../src/libmipster.la

########################################################################
#                   etDecsi tests (MIP node-limited, all sol valid)   #
########################################################################

CInterfaceTest_etDecsi_SOURCES = CInterfaceTest_etDecsi.c
nodist_EXTRA_CInterfaceTest_etDecsi_SOURCES = dummy.cpp
CInterfaceTest_etDecsi_CPPFLAGS = $(AM_CPPFLAGS) -DFIXTURE_DIR=\"$(FIXTURE_DIR_VAL)\"
CInterfaceTest_etDecsi_LDADD = ../src/libmipster.la

ctests: CInterfaceTest$(EXEEXT) CInterfaceTest_tsp_random$(EXEEXT) CInterfaceTest_fl_random$(EXEEXT) CInterfaceTest_mdkp_random$(EXEEXT) CInterfaceTest_nursesched$(EXEEXT) CInterfaceTest_a1$(EXEEXT) CInterfaceTest_graphdraw$(EXEEXT) CInterfaceTest_trdta5581$(EXEEXT) CInterfaceTest_trd445c$(EXEEXT) CInterfaceTest_j3041_1$(EXEEXT) CInterfaceTest_yue20013$(EXEEXT) CInterfaceTest_etDecsi$(EXEEXT) CbcSolverLpTest$(EXEEXT) CInterfaceTest_advanced$(EXEEXT) CInterfaceTest_colgen$(EXEEXT) CInterfaceTest_construction$(EXEEXT) CInterfaceTest_facility$(EXEEXT) CInterfaceTest_jssp$(EXEEXT) CInterfaceTest_kantorovich$(EXEEXT) CInterfaceTest_mip_reopt$(EXEEXT) CInterfaceTest_new$(EXEEXT) CInterfaceTest_tsp_cuts$(EXEEXT)
	MIPSTER_FIXTURE_DIR=$(abs_srcdir)/fixtures $(srcdir)/run-mipster-tests

########################################################################
#                   Embedded-data tests (no fixtures needed)           #
########################################################################

CInterfaceTest_advanced_SOURCES = CInterfaceTest_advanced.c
nodist_EXTRA_CInterfaceTest_advanced_SOURCES = dummy.cpp
CInterfaceTest_advanced_LDADD = ../src/libmipster.la

CInterfaceTest_colgen_SOURCES = CInterfaceTest_colgen.c
nodist_EXTRA_CInterfaceTest_colgen_SOURCES = dummy.cpp
CInterfaceTest_colgen_LDADD = ../src/libmipster.la

CInterfaceTest_construction_SOURCES = CInterfaceTest_construction.c
nodist_EXTRA_CInterfaceTest_construction_SOURCES = dummy.cpp
CInterfaceTest_construction_LDADD = ../src/libmipster.la

CInterfaceTest_facility_SOURCES = CInterfaceTest_facility.c
nodist_EXTRA_CInterfaceTest_facility_SOURCES = dummy.cpp
CInterfaceTest_facility_LDADD = ../src/libmipster.la

CInterfaceTest_jssp_SOURCES = CInterfaceTest_jssp.c
nodist_EXTRA_CInterfaceTest_jssp_SOURCES = dummy.cpp
CInterfaceTest_jssp_LDADD = ../src/libmipster.la

CInterfaceTest_kantorovich_SOURCES = CInterfaceTest_kantorovich.c
nodist_EXTRA_CInterfaceTest_kantorovich_SOURCES = dummy.cpp
CInterfaceTest_kantorovich_LDADD = ../src/libmipster.la

CInterfaceTest_mip_reopt_SOURCES = CInterfaceTest_mip_reopt.c
nodist_EXTRA_CInterfaceTest_mip_reopt_SOURCES = dummy.cpp
CInterfaceTest_mip_reopt_LDADD = ../src/libmipster.la

CInterfaceTest_new_SOURCES = CInterfaceTest_new.c
nodist_EXTRA_CInterfaceTest_new_SOURCES = dummy.cpp
CInterfaceTest_new_LDADD = ../src/libmipster.la

CInterfaceTest_tsp_cuts_SOURCES = CInterfaceTest_tsp_cuts.c
nodist_EXTRA_CInterfaceTest_tsp_cuts_SOURCES = dummy.cpp
CInterfaceTest_tsp_cuts_LDADD = ../src/libmipster.la
