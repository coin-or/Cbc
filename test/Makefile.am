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
  fixtures/yue20013.1.150.mps.gz \
  fixtures/attfInst1.mps.gz \
  fixtures/cvrp_asymmetric.mps.gz \
  fixtures/cvrp_capacity_loose.mps.gz \
  fixtures/cvrp_capacity_medium.mps.gz \
  fixtures/cvrp_capacity_tight.mps.gz \
  fixtures/cvrp_demand_outliers.mps.gz \
  fixtures/cvrp_demand_skewed.mps.gz \
  fixtures/cvrp_demand_uniform.mps.gz \
  fixtures/cvrp_dense_uniform.mps.gz \
  fixtures/cvrp_few_large.mps.gz \
  fixtures/cvrp_geo_clustered.mps.gz \
  fixtures/cvrp_geo_ring.mps.gz \
  fixtures/cvrp_high_utilization.mps.gz \
  fixtures/cvrp_large_scale.mps.gz \
  fixtures/cvrp_loose.mps.gz \
  fixtures/cvrp_many_small.mps.gz \
  fixtures/cvrp_medium.mps.gz \
  fixtures/cvrp_medium_scale.mps.gz \
  fixtures/cvrp_size_large.mps.gz \
  fixtures/cvrp_size_small.mps.gz \
  fixtures/cvrp_small_scale.mps.gz \
  fixtures/cvrp_sparse_grid.mps.gz \
  fixtures/cvrp_tight.mps.gz \
  fixtures/vrppd_large_scale.mps.gz \
  fixtures/vrppd_large_uniform.mps.gz \
  fixtures/vrppd_many_requests.mps.gz \
  fixtures/vrppd_medium_clustered.mps.gz \
  fixtures/vrppd_small.mps.gz \
  fixtures/vrppd_small_scale.mps.gz \
  fixtures/vrppd_small_uniform.mps.gz \
  fixtures/vrppd_tight.mps.gz \
  fixtures/vrppd_tight_ring.mps.gz \
  fixtures/qap_n4_s42_uniform_euclidean.mps.gz \
  fixtures/qap_n5_s42_uniform_euclidean.mps.gz \
  fixtures/qap_n5_s137_sparse_manhattan.mps.gz \
  fixtures/qap_n6_s42_uniform_grid.mps.gz \
  fixtures/qap_n6_s137_hub_euclidean.mps.gz \
  fixtures/qap_n7_s42_uniform_euclidean.mps.gz \
  fixtures/qap_n7_s137_sparse_grid.mps.gz \
  fixtures/scp_e10_s15_sd42_uniform.mps.gz \
  fixtures/scp_e15_s20_sd42_sparse.mps.gz \
  fixtures/scp_e20_s30_sd42_dense.mps.gz \
  fixtures/scp_e20_s25_sd137_hierarchical.mps.gz \
  fixtures/scp_e30_s50_sd42_uniform.mps.gz \
  fixtures/scp_e40_s60_sd137_sparse.mps.gz \
  fixtures/scp_e50_s80_sd42_dense.mps.gz \
  fixtures/spp_e10_s15_sd42_random.mps.gz \
  fixtures/spp_e12_s18_sd42_k_partite.mps.gz \
  fixtures/spp_e15_s20_sd137_geometric.mps.gz \
  fixtures/spp_e20_s25_sd42_scheduling.mps.gz \
  fixtures/spp_e30_s50_sd42_random.mps.gz \
  fixtures/spp_e35_s60_sd137_k_partite.mps.gz \
  fixtures/spp_e40_s70_sd42_geometric.mps.gz \
  fixtures/spp_e50_s80_sd137_scheduling.mps.gz \
  fixtures/sprt_e10_s20_sd42_uniform.mps.gz \
  fixtures/sprt_e15_s30_sd42_overlapping.mps.gz \
  fixtures/sprt_e30_s60_sd42_uniform.mps.gz \
  fixtures/sprt_e40_s80_sd137_overlapping.mps.gz \
  fixtures/bpc_n10_c100_sd42_random_uniform.mps.gz \
  fixtures/bpc_n12_c150_sd42_clique_diverse.mps.gz \
  fixtures/bpc_n15_c100_sd137_bipartite_uniform.mps.gz \
  fixtures/bpc_n15_c120_sd42_geometric_tight.mps.gz \
  fixtures/bpc_n20_c150_sd42_random_diverse.mps.gz \
  fixtures/bpc_n25_c200_sd137_clique_uniform.mps.gz \
  fixtures/bpc_n30_c180_sd42_geometric_tight.mps.gz \
  fixtures/mis_n15_sd42_random_uniform.mps.gz \
  fixtures/mis_n20_sd42_geometric_unit.mps.gz \
  fixtures/mis_n16_sd137_grid_diverse.mps.gz \
  fixtures/mis_n18_sd42_k_partite_uniform.mps.gz \
  fixtures/mis_n20_sd137_tree_uniform.mps.gz \
  fixtures/mis_n30_sd42_random_uniform.mps.gz \
  fixtures/mis_n35_sd137_geometric_diverse.mps.gz \
  fixtures/mis_n36_sd42_grid_unit.mps.gz \
  fixtures/mis_n40_sd137_planar_uniform.mps.gz \
  fixtures/mis_n30_sd42_cycle_uniform.mps.gz \
  fixtures/jssp_ft06.mps.gz \
  fixtures/jssp_ft10.mps.gz \
  fixtures/jssp_la01.mps.gz \
  fixtures/jssp_la06.mps.gz \
  fixtures/jssp_la11.mps.gz \
  fixtures/jssp_orb01.mps.gz \
  fixtures/steiner_grid4x4_t3_corners.mps.gz \
  fixtures/steiner_grid5x5_t4_sparse.mps.gz \
  fixtures/steiner_geom15_r30_t3.mps.gz \
  fixtures/steiner_geom20_r35_t5.mps.gz \
  fixtures/steiner_complete10_unit_t4.mps.gz \
  fixtures/steiner_complete12_rand_t3.mps.gz \
  fixtures/steiner_sparse20_p03_t4.mps.gz \
  fixtures/steiner_sparse25_p025_t6.mps.gz

# Reference optimal solution files (for -debugCuts diagnostic pass)
pkgtestsoldir = $(datadir)/mipster/test/fixtures/solutions
dist_pkgtestsol_DATA = \
  fixtures/solutions/A-1.sol \
  fixtures/solutions/etDecsi.sol \
  fixtures/solutions/j3041_1.sol \
  fixtures/solutions/nursesched-sprint02.sol \
  fixtures/solutions/trd445c.sol \
  fixtures/solutions/yue20013.1.150.sol \
  fixtures/solutions/miclsp_i04_t10_s042_r110.sol \
  fixtures/solutions/miclsp_i04_t10_s042_r114.sol \
  fixtures/solutions/miclsp_i04_t10_s042_r120.sol \
  fixtures/solutions/miclsp_i04_t10_s042_r130.sol \
  fixtures/solutions/miclsp_i04_t10_s137_r110.sol \
  fixtures/solutions/miclsp_i04_t10_s137_r114.sol \
  fixtures/solutions/miclsp_i04_t10_s137_r120.sol \
  fixtures/solutions/miclsp_i04_t12_s042_r110.sol \
  fixtures/solutions/miclsp_i04_t12_s042_r114.sol \
  fixtures/solutions/miclsp_i04_t12_s042_r120.sol \
  fixtures/solutions/miclsp_i04_t12_s137_r114.sol

# Runner script installs to $(bindir)/run-mipster-tests
dist_bin_SCRIPTS = run-mipster-tests run-miplib3-tests fetch-miplib3

# Desktop entry installs to $(datadir)/applications/
appsdir = $(datadir)/applications
dist_apps_DATA = mipster-tests.desktop mipster-tests-debug.desktop

########################################################################
#                         unitTest for Cbc                             #
########################################################################

CBC_TEST_TGTS = miplib3_tests ctests

# We are using the solver executables to do the unit test

test: $(CBC_TEST_TGTS)

# MIPLIB3 regression tests: per-instance C-interface tests run in parallel.
# Uses fetch-miplib3 to locate or clone the instance repository.
miplib3_tests: CInterfaceTest_miplib$(EXEEXT)
	$(srcdir)/run-miplib3-tests

if COIN_HAS_SAMPLE
test_cbc_generic: ../src/cbc-generic$(EXEEXT)
	@echo "A minimal test is all that cbc-generic can manage just yet."
	RUNNING_TEST="cbc-generic" ../src/cbc-generic$(EXEEXT) -directory `$(CYGPATH_W) $(SAMPLE_DATA)` -import p0033.mps -branch
endif


.PHONY: test

bin_PROGRAMS = CInterfaceTest CInterfaceTest_tsp_random CInterfaceTest_fl_random CInterfaceTest_mdkp_random CInterfaceTest_miclsp CInterfaceTest_mescn CInterfaceTest_nursesched CInterfaceTest_a1 CInterfaceTest_graphdraw CInterfaceTest_trdta5581 CInterfaceTest_trd445c CInterfaceTest_j3041_1 CInterfaceTest_yue20013 CInterfaceTest_etDecsi CInterfaceTest_attfInst1 CbcSolverLpTest CInterfaceTest_advanced CInterfaceTest_colgen CInterfaceTest_construction CInterfaceTest_facility CInterfaceTest_jssp CInterfaceTest_kantorovich CInterfaceTest_mip_reopt CInterfaceTest_new CInterfaceTest_tsp_cuts CInterfaceTest_vrp CInterfaceTest_miplib CInterfaceTest_qap CInterfaceTest_setcover CInterfaceTest_setpacking CInterfaceTest_setpartitioning CInterfaceTest_bpc CInterfaceTest_mis CInterfaceTest_jssp_fixtures CInterfaceTest_steiner

noinst_PROGRAMS = gen_miclsp_solutions

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
#                   MICLSP tests (binary + integer + continuous + Big-M)#
########################################################################

CInterfaceTest_miclsp_SOURCES = CInterfaceTest_miclsp.c
nodist_EXTRA_CInterfaceTest_miclsp_SOURCES = dummy.cpp
CInterfaceTest_miclsp_LDADD = ../src/libmipster.la
CInterfaceTest_miclsp_CPPFLAGS = $(AM_CPPFLAGS) -DFIXTURE_DIR=\"$(FIXTURE_DIR_VAL)\"

########################################################################
#                   MESCN tests (supply chain: binary + integer + continuous)#
########################################################################

CInterfaceTest_mescn_SOURCES = CInterfaceTest_mescn.c
nodist_EXTRA_CInterfaceTest_mescn_SOURCES = dummy.cpp
CInterfaceTest_mescn_LDADD = ../src/libmipster.la

########################################################################
#                   VRP tests (CVRP + VRPPD routing problems)          #
########################################################################

CInterfaceTest_vrp_SOURCES = CInterfaceTest_vrp.c
nodist_EXTRA_CInterfaceTest_vrp_SOURCES = dummy.cpp
CInterfaceTest_vrp_LDADD = ../src/libmipster.la

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

########################################################################
#                   attfInst1 tests (node-limited MIP, all sol valid)  #
########################################################################

CInterfaceTest_attfInst1_SOURCES = CInterfaceTest_attfInst1.c
nodist_EXTRA_CInterfaceTest_attfInst1_SOURCES = dummy.cpp
CInterfaceTest_attfInst1_CPPFLAGS = $(AM_CPPFLAGS) -DFIXTURE_DIR=\"$(FIXTURE_DIR_VAL)\"
CInterfaceTest_attfInst1_LDADD = ../src/libmipster.la

ctests: CInterfaceTest$(EXEEXT) CInterfaceTest_tsp_random$(EXEEXT) CInterfaceTest_fl_random$(EXEEXT) CInterfaceTest_mdkp_random$(EXEEXT) CInterfaceTest_miclsp$(EXEEXT) CInterfaceTest_nursesched$(EXEEXT) CInterfaceTest_a1$(EXEEXT) CInterfaceTest_graphdraw$(EXEEXT) CInterfaceTest_trdta5581$(EXEEXT) CInterfaceTest_trd445c$(EXEEXT) CInterfaceTest_j3041_1$(EXEEXT) CInterfaceTest_yue20013$(EXEEXT) CInterfaceTest_etDecsi$(EXEEXT) CInterfaceTest_attfInst1$(EXEEXT) CbcSolverLpTest$(EXEEXT) CInterfaceTest_advanced$(EXEEXT) CInterfaceTest_colgen$(EXEEXT) CInterfaceTest_construction$(EXEEXT) CInterfaceTest_facility$(EXEEXT) CInterfaceTest_jssp$(EXEEXT) CInterfaceTest_kantorovich$(EXEEXT) CInterfaceTest_mip_reopt$(EXEEXT) CInterfaceTest_new$(EXEEXT) CInterfaceTest_tsp_cuts$(EXEEXT) CInterfaceTest_qap$(EXEEXT) CInterfaceTest_setcover$(EXEEXT) CInterfaceTest_setpacking$(EXEEXT) CInterfaceTest_setpartitioning$(EXEEXT) CInterfaceTest_bpc$(EXEEXT) CInterfaceTest_mis$(EXEEXT)
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

CInterfaceTest_miplib_SOURCES = CInterfaceTest_miplib.c
nodist_EXTRA_CInterfaceTest_miplib_SOURCES = dummy.cpp
CInterfaceTest_miplib_LDADD = ../src/libmipster.la

########################################################################
#                   Pure binary problem tests (clique cuts)            #
########################################################################

CInterfaceTest_qap_SOURCES = CInterfaceTest_qap.c
nodist_EXTRA_CInterfaceTest_qap_SOURCES = dummy.cpp
CInterfaceTest_qap_CPPFLAGS = $(AM_CPPFLAGS) -DFIXTURE_DIR=\"$(FIXTURE_DIR_VAL)\"
CInterfaceTest_qap_LDADD = ../src/libmipster.la

CInterfaceTest_setcover_SOURCES = CInterfaceTest_setcover.c
nodist_EXTRA_CInterfaceTest_setcover_SOURCES = dummy.cpp
CInterfaceTest_setcover_CPPFLAGS = $(AM_CPPFLAGS) -DFIXTURE_DIR=\"$(FIXTURE_DIR_VAL)\"
CInterfaceTest_setcover_LDADD = ../src/libmipster.la

CInterfaceTest_setpacking_SOURCES = CInterfaceTest_setpacking.c
nodist_EXTRA_CInterfaceTest_setpacking_SOURCES = dummy.cpp
CInterfaceTest_setpacking_CPPFLAGS = $(AM_CPPFLAGS) -DFIXTURE_DIR=\"$(FIXTURE_DIR_VAL)\"
CInterfaceTest_setpacking_LDADD = ../src/libmipster.la

CInterfaceTest_setpartitioning_SOURCES = CInterfaceTest_setpartitioning.c
nodist_EXTRA_CInterfaceTest_setpartitioning_SOURCES = dummy.cpp
CInterfaceTest_setpartitioning_CPPFLAGS = $(AM_CPPFLAGS) -DFIXTURE_DIR=\"$(FIXTURE_DIR_VAL)\"
CInterfaceTest_setpartitioning_LDADD = ../src/libmipster.la

CInterfaceTest_bpc_SOURCES = CInterfaceTest_bpc.c
nodist_EXTRA_CInterfaceTest_bpc_SOURCES = dummy.cpp
CInterfaceTest_bpc_CPPFLAGS = $(AM_CPPFLAGS) -DFIXTURE_DIR=\"$(FIXTURE_DIR_VAL)\"
CInterfaceTest_bpc_LDADD = ../src/libmipster.la

CInterfaceTest_mis_SOURCES = CInterfaceTest_mis.c
nodist_EXTRA_CInterfaceTest_mis_SOURCES = dummy.cpp
CInterfaceTest_mis_CPPFLAGS = $(AM_CPPFLAGS) -DFIXTURE_DIR=\"$(FIXTURE_DIR_VAL)\"
CInterfaceTest_mis_LDADD = ../src/libmipster.la

CInterfaceTest_jssp_fixtures_SOURCES = CInterfaceTest_jssp_fixtures.c
nodist_EXTRA_CInterfaceTest_jssp_fixtures_SOURCES = dummy.cpp
CInterfaceTest_jssp_fixtures_CPPFLAGS = $(AM_CPPFLAGS) -DFIXTURE_DIR=\"$(FIXTURE_DIR_VAL)\"
CInterfaceTest_jssp_fixtures_LDADD = ../src/libmipster.la

CInterfaceTest_steiner_SOURCES = CInterfaceTest_steiner.c
nodist_EXTRA_CInterfaceTest_steiner_SOURCES = dummy.cpp
CInterfaceTest_steiner_CPPFLAGS = $(AM_CPPFLAGS) -DFIXTURE_DIR=\"$(FIXTURE_DIR_VAL)\"
CInterfaceTest_steiner_LDADD = ../src/libmipster.la

########################################################################
#             gen_miclsp_solutions — reference solution generator      #
########################################################################

gen_miclsp_solutions_SOURCES = gen_miclsp_solutions.c
nodist_EXTRA_gen_miclsp_solutions_SOURCES = dummy.cpp
gen_miclsp_solutions_LDADD = ../src/libmipster.la
