# Copyright (C) 2006, 2009 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                               libCbc                                 #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory
lib_LTLIBRARIES = libCbc.la libCbcSolver.la

# List all source files for this library, including headers
libCbc_la_SOURCES = \
	Cbc_C_Interface.cpp Cbc_C_Interface.h \
	CbcConfig.h \
	CbcBranchActual.cpp CbcBranchActual.hpp \
	CbcBranchBase.cpp CbcBranchBase.hpp \
	CbcBranchCut.cpp CbcBranchCut.hpp \
	CbcBranchDynamic.cpp CbcBranchDynamic.hpp \
	CbcBranchLotsize.cpp CbcBranchLotsize.hpp \
	CbcCompareActual.cpp CbcCompareActual.hpp \
	CbcCompareBase.hpp \
	CbcCountRowCut.cpp CbcCountRowCut.hpp \
	CbcCutGenerator.cpp CbcCutGenerator.hpp \
	CbcEventHandler.cpp CbcEventHandler.hpp \
	CbcFathom.cpp CbcFathom.hpp \
	CbcFathomDynamicProgramming.cpp CbcFathomDynamicProgramming.hpp \
	CbcFeasibilityBase.hpp \
	CbcHeuristic.cpp CbcHeuristic.hpp \
	CbcHeuristicDive.cpp CbcHeuristicDive.hpp \
	CbcHeuristicDiveCoefficient.cpp CbcHeuristicDiveCoefficient.hpp \
	CbcHeuristicDiveFractional.cpp CbcHeuristicDiveFractional.hpp \
	CbcHeuristicDiveGuided.cpp CbcHeuristicDiveGuided.hpp \
	CbcHeuristicDiveLineSearch.cpp CbcHeuristicDiveLineSearch.hpp \
	CbcHeuristicDivePseudoCost.cpp CbcHeuristicDivePseudoCost.hpp \
	CbcHeuristicDiveVectorLength.cpp CbcHeuristicDiveLength.hpp \
	CbcHeuristicFPump.cpp CbcHeuristicFPump.hpp \
	CbcHeuristicGreedy.cpp CbcHeuristicGreedy.hpp \
	CbcHeuristicLocal.cpp CbcHeuristicLocal.hpp \
	CbcHeuristicPivotAndFix.cpp CbcHeuristicPivotAndFix.hpp \
	CbcHeuristicRandRound.cpp CbcHeuristicRandRound.hpp \
	CbcHeuristicRINS.cpp CbcHeuristicRINS.hpp \
	CbcMessage.cpp CbcMessage.hpp \
	CbcModel.cpp CbcModel.hpp \
	CbcNode.cpp CbcNode.hpp \
	CbcStatistics.cpp CbcStatistics.hpp \
	CbcStrategy.cpp CbcStrategy.hpp \
	CbcTree.cpp CbcTree.hpp \
	CbcTreeLocal.cpp CbcTreeLocal.hpp

libCbcSolver_la_SOURCES = \
	CbcCbcParam.cpp \
	Cbc_ampl.cpp Cbc_ampl.h \
	ClpConstraintAmpl.hpp ClpAmplObjective.hpp ClpAmplStuff.cpp \
	CbcLinked.cpp CbcLinked.hpp \
	unitTestClp.cpp CbcSolver.cpp 

#Stefan: added unitTestClp.cpp because CbcSolver defines reference to it

## If we have the Ampl solver library, we need to add additional things
## also see lower down
## new idea is to add in source code anyway
##libCbc_la_SOURCES += Cbc_ampl.cpp Cbc_ampl.h


# This is for libtool (on Windows)

libCbc_la_LDFLAGS = $(LT_LDFLAGS)
libCbcSolver_la_LDFLAGS = $(LT_LDFLAGS)

# Some common setup for cbc and cbc-generic.

bin_PROGRAMS =
AM_CPPFLAGS = 

########################################################################
#                            cbc program                               #
########################################################################

if COIN_HAS_CLP

# Name of the executable compiled in this directory.  We want it to be
# installed in the 'bin' directory
bin_PROGRAMS += cbc

# List all source files for this executable, including headers
cbc_SOURCES = \
	CoinSolve.cpp

# Additional COIN libraries
cbc_LDADD = libCbcSolver.la libCbc.la

if COIN_HAS_CGL_PREINSTALLED
  LIBS += -L`$(CYGPATH_W) $(CGLINSTDIR)/lib` -lCgl
else
  cbc_LDADD += $(CGLOBJDIR)/src/libCgl.la
endif

if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiClp -lOsi
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp`
  cbc_LDADD += $(OSIOBJDIR)/src/OsiClp/libOsiClp.la \
	$(OSIOBJDIR)/src/libOsi.la
endif

if COIN_HAS_CLP_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CLPINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(CLPINSTDIR)/lib` -lClp \
	`cat $(CLPDOCDIR)/clp_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc`
  cbc_LDADD += $(CLPOBJDIR)/src/libClp.la
  LIBS  += `cat $(CLPOBJDIR)/clp_addlibs.txt`
endif

if COIN_HAS_COINUTILS_PREINSTALLED
  LIBS += -L`$(CYGPATH_W) $(COINUTILSINSTDIR)/lib` -lCoinUtils
else
  cbc_LDADD += $(COINUTILSOBJDIR)/src/libCoinUtils.la
endif

# Finally, the -rpath flag is used by libtool to make sure that the shared
# library is found (in the lib install directory) when we are using dynamic
# libraries.
cbc_LDFLAGS = -rpath $(libdir)

if COIN_HAS_CPX
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CPXINCDIR)`
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiCpx
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx`
  cbc_LDADD += $(OSIOBJDIR)/src/OsiCpx/libOsiCpx.la
endif
endif

# If we have the Ampl solver library, we need to add additional things
if COIN_HAS_ASL
#cbc_SOURCES += Cbc_ampl.cpp Cbc_ampl.h
  cbc_LDADD += $(ASLLIB)
  AM_CPPFLAGS += $(ASL_CPPFLAGS) -DCOIN_HAS_ASL
endif

if COIN_HAS_THIRDPARTY_BLAS
if COIN_HAS_THIRDPARTY_BLAS_PREINSTALLED
  LIBS += -L`$(CYGPATH_W) $(THIRDPARTY_BLASINSTDIR)/lib/ThirdParty` -lcoinblas
else
  cbc_LDADD += $(THIRDPARTY_BLASOBJDIR)/libcoinblas.la
endif
endif

if COIN_HAS_THIRDPARTY_LAPACK
if COIN_HAS_THIRDPARTY_LAPACK_PREINSTALLED
  LIBS += -L`$(CYGPATH_W) $(THIRDPARTY_LAPACKINSTDIR)/lib/ThirdParty` -lcoinlapack
else
  cbc_LDADD += $(THIRDPARTY_LAPACKOBJDIR)/libcoinlapack.la
endif
endif

endif


########################################################################
#                        cbc-generic program                           #
########################################################################

if CBC_BUILD_CBC_GENERIC

# Name of the executable compiled in this directory.  We want it to be
# installed in the 'bin' directory

bin_PROGRAMS += cbc-generic

# List all source files for this executable, including headers

cbc_generic_SOURCES = \
	CbcGenBaB.cpp \
	CbcGenCbcParam.cpp CbcGenCbcParam.hpp CbcGenCbcParamUtils.cpp \
	CbcGenCtlBlk.cpp CbcGenCtlBlk.hpp \
	CbcGeneric.cpp \
	CbcGenMessages.cpp CbcGenMessages.hpp \
	CbcGenOsiParam.cpp CbcGenOsiParam.hpp CbcGenOsiParamUtils.cpp \
	CbcGenParam.cpp CbcGenParam.hpp CbcGenParamUtils.cpp \
	CbcGenSolution.cpp \
	CbcGenSolvers.cpp

cbc_generic_LDADD = libCbcSolver.la libCbc.la

if COIN_HAS_CGL_PREINSTALLED
  LIBS += -L`$(CYGPATH_W) $(CGLINSTDIR)/lib` -lCgl
else
  cbc_generic_LDADD += $(CGLOBJDIR)/src/libCgl.la
endif

# List all additionally required solver and Osi libraries

if COIN_HAS_CLP
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiClp
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiClp/libOsiClp.la
endif
if COIN_HAS_CLP_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CLPINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(CLPINSTDIR)/lib` -lClp \
	`cat $(CLPDOCDIR)/clp_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc`
  cbc_generic_LDADD += $(CLPOBJDIR)/src/libClp.la
  LIBS += `cat $(CLPOBJDIR)/clp_addlibs.txt`
endif
endif

if COIN_HAS_CPX
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CPXINCDIR)`
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiCpx
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiCpx/libOsiCpx.la
endif
endif

if COIN_HAS_DYLP
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiDylp
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiDylp`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiDylp/libOsiDylp.la
endif
if COIN_HAS_DYLP_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(DYLPINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(DYLPINSTDIR)/lib` -lDylp \
	`cat $(DYLPDOCDIR)/dylp_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(DYLPSRCDIR)/src/Dylp` \
	-I`$(CYGPATH_W) $(DYLPSRCDIR)/src/DylpStdLib` \
	-I`$(CYGPATH_W) $(DYLPOBJDIR)/inc`
  cbc_generic_LDADD += $(DYLPOBJDIR)/src/Dylp/libDylp.la
  LIBS += `cat $(DYLPOBJDIR)/dylp_addlibs.txt`
endif
endif

if COIN_HAS_FMP
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiFmp
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiFmp`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiFmp/libOsiFmp.la
endif
endif

if COIN_HAS_THIRDPARTY_GLPK
if COIN_HAS_THIRDPARTY_GLPK_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(THIRDPARTY_GLPKINSTDIR)/include/coin/ThirdParty`
  LIBS += -L`$(CYGPATH_W) $(THIRDPARTY_GLPKINSTDIR)/lib/ThirdParty` -lcoinglpk
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(THIRDPARTY_GLPKSRCDIR)/glpk/include`
  cbc_generic_LDADD += $(THIRDPARTY_GLPKOBJDIR)/libcoinglpk.la
endif
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiGlpk
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiGlpk`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiGlpk/libOsiGlpk.la
endif
endif

if COIN_HAS_MSK
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(MSKINCDIR)`
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiMsk
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiMsk`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiMsk/libOsiMsk.la
endif
endif

if COIN_HAS_OSL
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSLINCDIR)`
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiOsl
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiOsl`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiOsl/libOsiOsl.la
endif
endif

if COIN_HAS_SPX
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SPXINCDIR)`
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiSpx
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiSpx`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiSpx/libOsiSpx.la
endif
endif

if COIN_HAS_SYMPHONY
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiSym
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiSym`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiSym/libOsiSym.la
endif
if COIN_HAS_SYMPHONY_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SYMPHONYINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(SYMPHONYINSTDIR)/lib` -lSym \
	`cat $(SYMPHONYDOCDIR)/sym_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SYMPHONYSRCDIR)/src`
  cbc_generic_LDADD += $(SYMPHONYOBJDIR)/src/libSym.la
  LIBS += `cat $(SYMPHONYOBJDIR)/sym_addlibs.txt`
endif
endif

if COIN_HAS_VOL
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiVol
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiVol`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiVol/libOsiVol.la
endif
if COIN_HAS_VOL_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(VOLINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(VOLINSTDIR)/lib` -lVol \
	`cat $(VOLDOCDIR)/vol_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(VOLSRCDIR)/src`
  cbc_generic_LDADD += $(VOLOBJDIR)/src/libVol.la
  LIBS += `cat $(VOLOBJDIR)/vol_addlibs.txt`
endif
endif

if COIN_HAS_XPR
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(XPRINCDIR)`
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiXpr
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiXpr`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiXpr/libOsiXpr.la
endif
endif

if COIN_HAS_ASL
  cbc_generic_LDADD += $(ASLLIB)
  AM_CPPFLAGS += $(ASL_CPPFLAGS) -DCOIN_HAS_ASL
endif

if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiClp -lOsi
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp`
  cbc_generic_LDADD += $(OSIOBJDIR)/src/OsiClp/libOsiClp.la \
	$(OSIOBJDIR)/src/libOsi.la
endif

if COIN_HAS_COINUTILS_PREINSTALLED
  LIBS += -L`$(CYGPATH_W) $(COINUTILSINSTDIR)/lib` -lCoinUtils
else
  cbc_generic_LDADD += $(COINUTILSOBJDIR)/src/libCoinUtils.la
endif

if COIN_HAS_THIRDPARTY_BLAS
if COIN_HAS_THIRDPARTY_BLAS_PREINSTALLED
  LIBS += -L`$(CYGPATH_W) $(THIRDPARTY_BLASINSTDIR)/lib/ThirdParty` -lcoinblas
else
  cbc_generic_LDADD += $(THIRDPARTY_BLASOBJDIR)/libcoinblas.la
endif
endif

if COIN_HAS_THIRDPARTY_LAPACK
if COIN_HAS_THIRDPARTY_LAPACK_PREINSTALLED
  LIBS += -L`$(CYGPATH_W) $(THIRDPARTY_LAPACKINSTDIR)/lib/ThirdParty` -lcoinlapack
else
  cbc_generic_LDADD += $(THIRDPARTY_LAPACKOBJDIR)/libcoinlapack.la
endif
endif

# Finally, the -rpath flag is used by libtool to make sure that the shared
# library is found (in the lib install directory) when we are using dynamic
# libraries.

cbc_generic_LDFLAGS = -rpath $(libdir)

endif

########################################################################
#                            Additional flags                          #
########################################################################

if COIN_HAS_CGL_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CGLINSTDIR)/include/coin`
  LIBS += `cat $(CGLDOCDIR)/cgl_addlibs.txt`
else
# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
# "top_srcdir" refers to the basic directory for the main package that is
# being compiled.
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglClique` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglDuplicateRow` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglFlowCover` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglGomory` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglKnapsackCover` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglMixedIntegerRounding` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglMixedIntegerRounding2` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglOddHole` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglPreProcess` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglProbing` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglRedSplit` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglTwomir` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglLandP` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglResidualCapacity` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglZeroHalf`
  LIBS += `cat $(CGLOBJDIR)/cgl_addlibs.txt`
endif

if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += `cat $(OSIDOCDIR)/osi_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSIOBJDIR)/inc`
  LIBS += `cat $(OSIOBJDIR)/osi_addlibs.txt`
endif

if COIN_HAS_COINUTILS_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(COINUTILSINSTDIR)/include/coin`
  LIBS += `cat $(COINUTILSDOCDIR)/coinutils_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc`
  LIBS += `cat $(COINUTILSOBJDIR)/coinutils_addlibs.txt`
endif

LIBS += $(ADDLIBS)

# List additional defines
AM_CPPFLAGS += -DCOIN_NO_CLP_MESSAGE -DUSE_CBCCONFIG

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
	../inc/config_cbc.h \
	Cbc_C_Interface.h \
	CbcBranchActual.hpp \
	CbcBranchBase.hpp \
	CbcBranchDynamic.hpp \
	CbcBranchLotsize.hpp \
	CbcBranchCut.hpp \
	CbcCompareActual.hpp \
	CbcCompareBase.hpp \
	CbcCutGenerator.hpp \
	CbcFathom.hpp \
	CbcEventHandler.hpp \
	CbcFeasibilityBase.hpp \
	CbcHeuristic.hpp \
	CbcHeuristicDive.hpp \
	CbcHeuristicDiveCoefficient.hpp \
	CbcHeuristicDiveFractional.hpp \
	CbcHeuristicDiveGuided.hpp \
	CbcHeuristicDiveLineSearch.hpp \
	CbcHeuristicDivePseudoCost.hpp \
	CbcHeuristicDiveVectorLength.hpp \
	CbcHeuristicFPump.hpp \
	CbcHeuristicGreedy.hpp \
	CbcHeuristicLocal.hpp \
	CbcHeuristicPivotAndFix.hpp \
	CbcHeuristicRandRound.hpp \
	CbcHeuristicRINS.hpp \
	CbcMessage.hpp \
	CbcModel.hpp \
	CbcNode.hpp \
	CbcStrategy.hpp \
	CbcSolver.hpp \
	CbcTree.hpp \
	CbcLinked.hpp \
	CbcTreeLocal.hpp \
	ClpConstraintAmpl.hpp \
	ClpAmplObjective.hpp 

#############################################################################
# Create the Config.h file that always defines HAVE_CONFIG_H and install it #
#############################################################################

# You only need to adapt the following line
ConfigHeader = CbcConfig.h

install-exec-local:
	echo "#ifndef HAVE_CONFIG_H" >bla
	echo "#define HAVE_CONFIG_H" >>bla
	echo "#endif" >> bla
	cat $(srcdir)/$(ConfigHeader) >> bla
	$(install_sh_DATA) bla $(DESTDIR)$(includecoindir)/$(ConfigHeader)
	rm -f bla

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/$(ConfigHeader)
