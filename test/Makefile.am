# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter           IBM    2006-04-13

AM_CPPFLAGS = -I$(srcdir)/../src -I$(srcdir)/../src/CoinUtils -I$(srcdir)/../src/Clp -I$(srcdir)/../src/Cgl $(CBCUNITTEST_CFLAGS)

AM_LDFLAGS = $(LT_LDFLAGS)

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

bin_PROGRAMS = CInterfaceTest CInterfaceTest_tsp_random CInterfaceTest_fl_random

########################################################################
#                         unitTest for CbcCInterface                   #
########################################################################

CInterfaceTest_SOURCES = CInterfaceTest.c
nodist_EXTRA_CInterfaceTest_SOURCES = dummy.cpp # force using C++ linker

CInterfaceTest_LDADD = ../src/libCbc.la

########################################################################
#                   TSP random-instance tests                           #
########################################################################

CInterfaceTest_tsp_random_SOURCES = CInterfaceTest_tsp_random.c
nodist_EXTRA_CInterfaceTest_tsp_random_SOURCES = dummy.cpp
CInterfaceTest_tsp_random_LDADD = ../src/libCbc.la

########################################################################
#                   FL random-instance tests (UFL + CFLP)              #
########################################################################

CInterfaceTest_fl_random_SOURCES = CInterfaceTest_fl_random.c
nodist_EXTRA_CInterfaceTest_fl_random_SOURCES = dummy.cpp
CInterfaceTest_fl_random_LDADD = ../src/libCbc.la

ctests: CInterfaceTest$(EXEEXT) CInterfaceTest_tsp_random$(EXEEXT) CInterfaceTest_fl_random$(EXEEXT)
	export RUNNING_TEST="CInterfaceTest" ; ./CInterfaceTest$(EXEEXT) $(MIPLIB3_DATA) $(SAMPLE_DATA)
	export RUNNING_TEST="CInterfaceTest_tsp_random" ; ./CInterfaceTest_tsp_random$(EXEEXT)
	export RUNNING_TEST="CInterfaceTest_fl_random" ; ./CInterfaceTest_fl_random$(EXEEXT)
