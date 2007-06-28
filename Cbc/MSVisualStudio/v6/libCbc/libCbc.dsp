# Microsoft Developer Studio Project File - Name="libCbc" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=libCbc - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "libCbc.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "libCbc.mak" CFG="libCbc - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "libCbc - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "libCbc - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "libCbc - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GR /GX /O2 /I "..\..\..\..\Cbc\src\\" /I "..\..\..\..\Cgl\src\CglTwomir" /I "..\..\..\..\Cgl\src\CglMixedIntegerRounding" /I "..\..\..\..\Cgl\src\CglMixedIntegerRounding2" /I "..\..\..\..\Cgl\src\CglFlowCover" /I "..\..\..\..\Cgl\src\CglClique" /I "..\..\..\..\Cgl\src\CglOddHole" /I "..\..\..\..\Cgl\src\CglKnapsackCover" /I "..\..\..\..\Cgl\src\CglGomory" /I "..\..\..\..\Cgl\src\CglPreProcess" /I "..\..\..\..\Cgl\src\CglProbing" /I "..\..\..\..\Cgl\src" /I "..\..\..\..\Clp\src" /I "..\..\..\..\Osi\src\OsiClp" /I "..\..\..\..\Osi\src" /I "..\..\..\..\CoinUtils\src" /I "..\..\..\..\BuildTools\headers" /D "NDEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "USE_CBCCONFIG" /D "COIN_NO_CLP_MESSAGE" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "libCbc - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GR /GX /ZI /Od /I "..\..\..\..\Cbc\src\\" /I "..\..\..\..\Cgl\src\CglTwomir" /I "..\..\..\..\Cgl\src\CglMixedIntegerRounding" /I "..\..\..\..\Cgl\src\CglMixedIntegerRounding2" /I "..\..\..\..\Cgl\src\CglFlowCover" /I "..\..\..\..\Cgl\src\CglClique" /I "..\..\..\..\Cgl\src\CglOddHole" /I "..\..\..\..\Cgl\src\CglKnapsackCover" /I "..\..\..\..\Cgl\src\CglGomory" /I "..\..\..\..\Cgl\src\CglPreProcess" /I "..\..\..\..\Cgl\src\CglProbing" /I "..\..\..\..\Cgl\src" /I "..\..\..\..\Clp\src" /I "..\..\..\..\Osi\src\OsiClp" /I "..\..\..\..\Osi\src" /I "..\..\..\..\CoinUtils\src" /I "..\..\..\..\BuildTools\headers" /D "_DEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "USE_CBCCONFIG" /D "COIN_NO_CLP_MESSAGE" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "libCbc - Win32 Release"
# Name "libCbc - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\Cbc_C_Interface.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcBranchActual.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcBranchBase.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcBranchCut.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcBranchDynamic.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcBranchLotsize.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcCbcParam.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcCompareActual.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcCountRowCut.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcCutGenerator.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcEventHandler.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcFathom.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcFathomDynamicProgramming.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcHeuristic.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcHeuristicFPump.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcHeuristicGreedy.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcHeuristicLocal.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcMessage.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcModel.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcNode.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcParam.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcStatistics.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcStrategy.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcTree.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcTreeLocal.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\Cbc_ampl.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\Cbc_C_Interface.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcBranchActual.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcBranchBase.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcBranchCut.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcBranchDynamic.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcBranchLotsize.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcCompareActual.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcCompareBase.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcConfig.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcCountRowCut.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcCutGenerator.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcEventHandler.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcFathom.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcFathomDynamicProgramming.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcFeasibilityBase.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcHeuristic.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcHeuristicFPump.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcHeuristicGreedy.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcHeuristicLocal.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcMessage.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcModel.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcNode.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\CbcOrClpParam.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcParam.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcStatistics.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcStrategy.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcTree.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cbc\src\CbcTreeLocal.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglClique\CglClique.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglCutGenerator.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglFlowCover\CglFlowCover.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglGomory\CglGomory.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglKnapsackCover\CglKnapsackCover.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglMixedIntegerRounding\CglMixedIntegerRounding.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglMixedIntegerRounding2\CglMixedIntegerRounding2.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglOddHole\CglOddHole.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglPreProcess\CglPreProcess.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglProbing\CglProbing.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Cgl\src\CglTwomir\CglTwomir.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpConfig.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpFactorization.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpMatrixBase.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpModel.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpObjective.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpParameters.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpPresolve.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpSimplex.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpSimplexOther.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpSolve.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\Coin_C_defines.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinDistance.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinError.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinFactorization.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinFileIO.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinFinite.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinFloatEqual.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinHelperFunctions.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinIndexedVector.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinMessage.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinMessageHandler.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinMpsIO.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPackedMatrix.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPackedVector.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPackedVectorBase.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPragma.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPresolveMatrix.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinShallowPackedVector.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinSort.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinTime.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinWarmStart.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinWarmStartBasis.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\BuildTools\headers\configall_system.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\BuildTools\headers\configall_system_msc.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Osi\src\OsiAuxInfo.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Osi\src\OsiClp\OsiClpSolverInterface.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Osi\src\OsiColCut.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Osi\src\OsiCollections.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Osi\src\OsiCut.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Osi\src\OsiCuts.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Osi\src\OsiPresolve.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Osi\src\OsiRowCut.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Osi\src\OsiRowCutDebugger.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Osi\src\OsiSolverBranch.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Osi\src\OsiSolverInterface.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Osi\src\OsiSolverParameters.hpp
# End Source File
# End Group
# End Target
# End Project
