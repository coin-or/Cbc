# Makefile.am for Clp/Osi/OsiClp — merged into Cbc build tree.
# All LP solver sources are compiled as a single non-installed library
# that is linked into libCbc.la.

########################################################################
#                              libClp                                  #
########################################################################

noinst_LTLIBRARIES = libClp.la

libClp_la_SOURCES = \
	ClpConfig.h \
	OsiConfig.h \
	OsiClpConfig.h \
	CbcOrClpParam.cpp \
	CbcOrClpParam.hpp \
	ClpAbortHandler.hpp \
	ClpCholeskyBase.cpp \
	ClpCholeskyBase.hpp \
	ClpCholeskyDense.cpp \
	ClpCholeskyDense.hpp \
	ClpCholeskyUfl.cpp \
	ClpCholeskyUfl.hpp \
	ClpConstraint.cpp \
	ClpConstraint.hpp \
	ClpConstraintLinear.cpp \
	ClpConstraintLinear.hpp \
	ClpConstraintQuadratic.cpp \
	ClpConstraintQuadratic.hpp \
	ClpDualRowDantzig.cpp \
	ClpDualRowDantzig.hpp \
	ClpDualRowPivot.cpp \
	ClpDualRowPivot.hpp \
	ClpDualRowSteepest.cpp \
	ClpDualRowSteepest.hpp \
	ClpDummyMatrix.cpp \
	ClpDummyMatrix.hpp \
	ClpDynamicExampleMatrix.cpp \
	ClpDynamicExampleMatrix.hpp \
	ClpDynamicMatrix.cpp \
	ClpDynamicMatrix.hpp \
	ClpEventHandler.cpp \
	ClpEventHandler.hpp \
	ClpFactorization.cpp \
	ClpFactorization.hpp \
	ClpGubDynamicMatrix.cpp \
	ClpGubDynamicMatrix.hpp \
	ClpGubMatrix.cpp \
	ClpGubMatrix.hpp \
	ClpHelperFunctions.cpp \
	ClpHelperFunctions.hpp \
	ClpInterior.cpp \
	ClpInterior.hpp \
	ClpLinearObjective.cpp \
	ClpLinearObjective.hpp \
	ClpMatrixBase.cpp \
	ClpMatrixBase.hpp \
	ClpMessage.cpp \
	ClpMessage.hpp \
	ClpModel.cpp \
	ClpModel.hpp \
	ClpModelParameters.hpp \
	ClpNetworkBasis.cpp \
	ClpNetworkBasis.hpp \
	ClpNetworkMatrix.cpp \
	ClpNetworkMatrix.hpp \
	ClpNode.cpp \
	ClpNode.hpp \
	ClpNonLinearCost.cpp \
	ClpNonLinearCost.hpp \
	ClpObjective.cpp \
	ClpObjective.hpp \
	ClpOutput.cpp \
	ClpOutput.hpp \
	ClpPEDualRowDantzig.cpp \
	ClpPEDualRowDantzig.hpp \
	ClpPEDualRowSteepest.cpp \
	ClpPEDualRowSteepest.hpp \
	ClpPEPrimalColumnDantzig.cpp \
	ClpPEPrimalColumnDantzig.hpp \
	ClpPEPrimalColumnSteepest.cpp \
	ClpPEPrimalColumnSteepest.hpp \
	ClpPESimplex.cpp \
	ClpPESimplex.hpp \
	ClpPackedMatrix.cpp \
	ClpPackedMatrix.hpp \
	ClpParam.cpp \
	ClpParam.hpp \
	ClpParamUtils.cpp \
	ClpParamUtils.hpp \
	ClpParameters.cpp \
	ClpParameters.hpp \
	ClpPlusMinusOneMatrix.cpp \
	ClpPlusMinusOneMatrix.hpp \
	ClpPredictorCorrector.cpp \
	ClpPredictorCorrector.hpp \
	ClpPresolve.cpp \
	ClpPresolve.hpp \
	ClpPrimalColumnDantzig.cpp \
	ClpPrimalColumnDantzig.hpp \
	ClpPrimalColumnPivot.cpp \
	ClpPrimalColumnPivot.hpp \
	ClpPrimalColumnSteepest.cpp \
	ClpPrimalColumnSteepest.hpp \
	ClpPrimalQuadraticDantzig.hpp \
	ClpQuadraticObjective.cpp \
	ClpQuadraticObjective.hpp \
	ClpRacingSolver.cpp \
	ClpRacingSolver.hpp \
	ClpSimplex.cpp \
	ClpSimplex.hpp \
	ClpSimplexDual.cpp \
	ClpSimplexDual.hpp \
	ClpSimplexNonlinear.cpp \
	ClpSimplexNonlinear.hpp \
	ClpSimplexOther.cpp \
	ClpSimplexOther.hpp \
	ClpSimplexPrimal.cpp \
	ClpSimplexPrimal.hpp \
	ClpSolve.cpp \
	ClpSolve.hpp \
	ClpSolver.cpp \
	ClpSolver.hpp \
	Clp_C_Interface.cpp \
	Clp_C_Interface.h \
	IdiSolve.cpp \
	Idiot.cpp \
	Idiot.hpp \
	OsiAuxInfo.cpp \
	OsiAuxInfo.hpp \
	OsiBranchingObject.cpp \
	OsiBranchingObject.hpp \
	OsiChooseVariable.cpp \
	OsiChooseVariable.hpp \
	OsiClpSolverInterface.cpp \
	OsiClpSolverInterface.hpp \
	OsiColCut.cpp \
	OsiColCut.hpp \
	OsiCollections.hpp \
	OsiColumnType.hpp \
	OsiCut.cpp \
	OsiCut.hpp \
	OsiCuts.cpp \
	OsiCuts.hpp \
	OsiFeatures.cpp \
	OsiFeatures.hpp \
	OsiNames.cpp \
	OsiPresolve.cpp \
	OsiPresolve.hpp \
	OsiRowCut.cpp \
	OsiRowCut.hpp \
	OsiRowCutDebugger.cpp \
	OsiRowCutDebugger.hpp \
	OsiSolverBranch.cpp \
	OsiSolverBranch.hpp \
	OsiSolverInterface.cpp \
	OsiSolverInterface.hpp \
	OsiSolverParameters.hpp

AM_CPPFLAGS = $(CBCLIB_CFLAGS) -I$(srcdir) -I$(top_srcdir)/src -I$(top_builddir)/src -I$(top_srcdir)/src/Cgl

AM_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#              Headers installed into $(includedir)/coin-or            #
########################################################################

includecoindir = $(includedir)/coin-or
includecoin_HEADERS = \
	ClpAbortHandler.hpp \
	Clp_C_Interface.h \
	ClpCholeskyBase.hpp \
	ClpCholeskyDense.hpp \
	ClpCholeskyUfl.hpp \
	ClpConstraint.hpp \
	ClpConstraintLinear.hpp \
	ClpConstraintQuadratic.hpp \
	ClpDualRowDantzig.hpp \
	ClpDualRowPivot.hpp \
	ClpDualRowSteepest.hpp \
	ClpDummyMatrix.hpp \
	ClpDynamicExampleMatrix.hpp \
	ClpDynamicMatrix.hpp \
	ClpEventHandler.hpp \
	ClpFactorization.hpp \
	ClpGubDynamicMatrix.hpp \
	ClpGubMatrix.hpp \
	ClpInterior.hpp \
	ClpLinearObjective.hpp \
	ClpMatrixBase.hpp \
	ClpMessage.hpp \
	ClpModel.hpp \
	ClpModelParameters.hpp \
	ClpNetworkMatrix.hpp \
	ClpNonLinearCost.hpp \
	ClpNode.hpp \
	ClpObjective.hpp \
	ClpOutput.hpp \
	ClpPackedMatrix.hpp \
	ClpParameters.hpp \
	ClpParam.hpp \
	ClpParamUtils.hpp \
	ClpPlusMinusOneMatrix.hpp \
	ClpPresolve.hpp \
	ClpRacingSolver.hpp \
	ClpPrimalColumnDantzig.hpp \
	ClpPrimalColumnPivot.hpp \
	ClpPrimalColumnSteepest.hpp \
	ClpQuadraticObjective.hpp \
	ClpSimplex.hpp \
	ClpSimplexNonlinear.hpp \
	ClpSimplexOther.hpp \
	ClpSimplexDual.hpp \
	ClpSimplexPrimal.hpp \
	ClpSolve.hpp \
	Idiot.hpp \
	ClpPESimplex.hpp \
	ClpPEPrimalColumnDantzig.hpp \
	ClpPEPrimalColumnSteepest.hpp \
	ClpPEDualRowDantzig.hpp \
	ClpPEDualRowSteepest.hpp \
	CbcOrClpParam.hpp \
	CbcOrClpParam.cpp \
	OsiConfig.h \
	OsiAuxInfo.hpp \
	OsiBranchingObject.hpp \
	OsiChooseVariable.hpp \
	OsiColCut.hpp \
	OsiCollections.hpp \
	OsiColumnType.hpp \
	OsiCut.hpp \
	OsiCuts.hpp \
	OsiPresolve.hpp \
	OsiRowCut.hpp \
	OsiRowCutDebugger.hpp \
	OsiSolverBranch.hpp \
	OsiSolverInterface.hpp \
	OsiSolverParameters.hpp \
	OsiFeatures.hpp \
	OsiClpSolverInterface.hpp

install-exec-local:
	$(install_sh_DATA) config_clp.h $(DESTDIR)$(includecoindir)/ClpConfig.h
	$(install_sh_DATA) config_osiclp.h $(DESTDIR)$(includecoindir)/OsiClpConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/ClpConfig.h $(DESTDIR)$(includecoindir)/OsiClpConfig.h
