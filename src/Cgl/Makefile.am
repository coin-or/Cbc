# Makefile.am for Cgl — merged into Cbc build tree.
# All Cgl cut generators compiled as a single non-installed library
# that is linked into libCbc.la.

########################################################################
#                              libCgl                                  #
########################################################################

noinst_LTLIBRARIES = libCgl.la

libCgl_la_SOURCES = \
	CglConfig.h \
	Cgl012cut.cpp Cgl012cut.hpp \
	CglAllDifferent.cpp CglAllDifferent.hpp \
	CglBKClique.cpp CglBKClique.hpp \
	CglClique.cpp CglCliqueHelper.cpp CglClique.hpp \
	CglCliqueStrengthening.cpp CglCliqueStrengthening.hpp \
	CglCutGenerator.cpp CglCutGenerator.hpp \
	CglDuplicateRow.cpp CglDuplicateRow.hpp \
	CglFlowCover.cpp CglFlowCover.hpp \
	CglGMI.cpp CglGMI.hpp \
	CglGMIParam.cpp CglGMIParam.hpp \
	CglGomory.cpp CglGomory.hpp \
	CglKnapsackCover.cpp CglKnapsackCover.hpp \
	CglLandP.cpp CglLandP.hpp \
	CglLandPMessages.cpp CglLandPMessages.hpp \
	CglLandPSimplex.cpp CglLandPSimplex.hpp \
	CglLandPTabRow.cpp CglLandPTabRow.hpp \
	CglLandPUtils.cpp CglLandPUtils.hpp \
	CglLandPValidator.cpp CglLandPValidator.hpp \
	CglLiftAndProject.cpp CglLiftAndProject.hpp \
	CglMessage.cpp CglMessage.hpp \
	CglMixedIntegerRounding.cpp CglMixedIntegerRounding.hpp \
	CglMixedIntegerRounding2.cpp CglMixedIntegerRounding2.hpp \
	CglOddHole.cpp CglOddHole.hpp \
	CglOddWheel.cpp CglOddWheel.hpp \
	CglParam.cpp CglParam.hpp \
	CglPreProcess.cpp CglPreProcess.hpp \
	CglProbing.cpp CglProbing.hpp \
	CglRedSplit.cpp CglRedSplit.hpp \
	CglRedSplitParam.cpp CglRedSplitParam.hpp \
	CglRedSplit2.cpp CglRedSplit2.hpp \
	CglRedSplit2Param.cpp CglRedSplit2Param.hpp \
	CglResidualCapacity.cpp CglResidualCapacity.hpp \
	CglSimpleRounding.cpp CglSimpleRounding.hpp \
	CglStored.cpp CglStored.hpp \
	CglTreeInfo.cpp CglTreeInfo.hpp \
	CglTwomir.cpp CglTwomir.hpp \
	CglZeroHalf.cpp CglZeroHalf.hpp

# Use same flags as the rest of the Cbc build (includes CoinUtils plus
# merged Clp/Osi/OsiClp headers from src/Clp).
AM_CPPFLAGS = $(CBCLIB_CFLAGS) -I$(srcdir) -I$(top_srcdir)/src/CoinUtils -I$(top_srcdir)/src/Clp

AM_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#              Headers installed into $(includedir)/coin-or            #
########################################################################

includecoindir = $(includedir)/mipster
includecoin_HEADERS = \
	CglConfig.h \
	Cgl012cut.hpp \
	CglAllDifferent.hpp \
	CglBKClique.hpp \
	CglClique.hpp \
	CglCliqueStrengthening.hpp \
	CglCutGenerator.hpp \
	CglDuplicateRow.hpp \
	CglFlowCover.hpp \
	CglGMI.hpp \
	CglGMIParam.hpp \
	CglGomory.hpp \
	CglKnapsackCover.hpp \
	CglLandP.hpp \
	CglLandPValidator.hpp \
	CglLiftAndProject.hpp \
	CglMessage.hpp \
	CglMixedIntegerRounding.hpp \
	CglMixedIntegerRounding2.hpp \
	CglOddHole.hpp \
	CglOddWheel.hpp \
	CglParam.hpp \
	CglPreProcess.hpp \
	CglProbing.hpp \
	CglRedSplit.hpp \
	CglRedSplitParam.hpp \
	CglRedSplit2.hpp \
	CglRedSplit2Param.hpp \
	CglResidualCapacity.hpp \
	CglSimpleRounding.hpp \
	CglStored.hpp \
	CglTreeInfo.hpp \
	CglTwomir.hpp \
	CglZeroHalf.hpp

install-exec-local:
	$(install_sh_DATA) $(srcdir)/config_cgl.h $(DESTDIR)$(includecoindir)/CglConfig.h
