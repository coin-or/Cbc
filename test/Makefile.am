# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter           IBM    2006-04-13

AM_CPPFLAGS = -I$(srcdir)/../src -I$(srcdir)/../src/OsiCbc $(CBCUNITTEST_CFLAGS)

AM_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#                         unitTest for Cbc                             #
########################################################################

CBC_TEST_TGTS = gamstests ositests

if COIN_HAS_CLP
  CBC_TEST_TGTS += test_cbc ctests
endif

if CBC_BUILD_CBC_GENERIC
if COIN_HAS_SAMPLE
  CBC_TEST_TGTS += test_cbc_generic
endif
endif

cbcflags =
ositestsflags =
if COIN_HAS_SAMPLE
#  cbcflags += -dirSample `$(CYGPATH_W) $(SAMPLE_DATA)`
  ositestsflags += -mpsDir=`$(CYGPATH_W) $(SAMPLE_DATA)`
endif
if COIN_HAS_MIPLIB3
  cbcflags += -dirMiplib `$(CYGPATH_W) $(MIPLIB3_DATA)` -unitTest
#  cbcflags += -dirMiplib `$(CYGPATH_W) $(MIPLIB3_DATA)` -miplib
endif
if COIN_HAS_NETLIB
  ositestsflags += -netlibDir=`$(CYGPATH_W) $(NETLIB_DATA)` -testOsiSolverInterface
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

gamstests: gamsTest$(EXEEXT)
	export RUNNING_TEST="gamsTest" ; ./gamsTest$(EXEEXT)

.PHONY: test

bin_PROGRAMS = gamsTest osiUnitTest CInterfaceTest

gamsTest_SOURCES = gamsTest.cpp
gamsTest_LDADD = ../src/libCbcSolver.la ../src/libCbc.la

########################################################################
#                         unitTest for OsiCbc                          #
########################################################################

osiUnitTest_SOURCES = osiUnitTest.cpp OsiCbcSolverInterfaceTest.cpp

osiUnitTest_LDADD = ../src/OsiCbc/libOsiCbc.la ../src/libCbc.la $(CBCUNITTEST_LFLAGS)

ositests: osiUnitTest$(EXEEXT)
	export RUNNING_TEST="osiUnitTest" ; ./osiUnitTest$(EXEEXT) $(ositestsflags)

########################################################################
#                         unitTest for CbcCInterface                   #
########################################################################

CInterfaceTest_SOURCES = CInterfaceTest.c
nodist_EXTRA_CInterfaceTest_SOURCES = dummy.cpp # force using C++ linker

CInterfaceTest_LDADD = ../src/libCbcSolver.la ../src/libCbc.la

ctests: CInterfaceTest$(EXEEXT)
	export RUNNING_TEST="CInterfaceTest" ; ./CInterfaceTest$(EXEEXT) $(MIPLIB3_DATA) $(SAMPLE_DATA)
