/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/
/*
  This file is part of cbc-generic.
*/

/*
  This file contains routines related to handling Osi solvers. The technique
  is to maintain a map of OsiSolverInterface objects as prototypes of the
  available solvers.
*/

#include "CbcConfig.h"
#include "CoinPragma.hpp"

#include "OsiSolverInterface.hpp"

/*
  Include class definitions for the solvers that are available. If
  CBC_DEFAULT_SOLVER is not defined in CbcConfig.h, set it to the first
  available solver.

  NOTE: Processing of keyword parameters is made case-independent by forcing
	lower case comparison. Maps, on the other hand, are case-sensitive.
	The solver names given here must contain only lower case letters and
	must match the keywords used when defining the keywords for the
	solver parameter.
*/

#ifdef COIN_HAS_CLP
# include "OsiClpSolverInterface.hpp"
# ifndef CBC_DEFAULT_SOLVER
#   define CBC_DEFAULT_SOLVER "clp"
# endif
#endif

#ifdef COIN_HAS_CPX
# include "OsiCpxSolverInterface.hpp"
# ifndef CBC_DEFAULT_SOLVER
#   define CBC_DEFAULT_SOLVER "cpx"
# endif
#endif

#ifdef COIN_HAS_DYLP
# include "OsiDylpSolverInterface.hpp"
# ifndef CBC_DEFAULT_SOLVER
#   define CBC_DEFAULT_SOLVER "dylp"
# endif
#endif

#ifdef COIN_HAS_GLPK
# include "OsiGlpkSolverInterface.hpp"
# ifndef CBC_DEFAULT_SOLVER
#   define CBC_DEFAULT_SOLVER "glpk"
# endif
#endif

#ifdef COIN_HAS_MSK
# include "OsiMskSolverInterface.hpp"
# ifndef CBC_DEFAULT_SOLVER
#   define CBC_DEFAULT_SOLVER "msk"
# endif
#endif

#ifdef COIN_HAS_SPX
# include "OsiSpxSolverInterface.hpp"
# ifndef CBC_DEFAULT_SOLVER
#   define CBC_DEFAULT_SOLVER "spx"
# endif
#endif

#include "CoinParam.hpp"

#include "CbcModel.hpp"
#include "CbcGenCtlBlk.hpp"

#include "CbcGenParam.hpp"
#include "CbcGenCbcParam.hpp"
#include "CbcGenOsiParam.hpp"

namespace {

char svnid[] = "$Id: CbcGenSolvers.cpp 1173 2009-06-04 09:44:10Z forrest $" ;

}

/*
  Unnamed local namespace to hide the data structures used to maintain the
  vector of OsiSolverInterface objects.
*/

namespace {

/*
  Data types for a vector of OsiSolverInterface objects.
*/
typedef std::map<std::string, OsiSolverInterface*> solverMap_t ;
typedef solverMap_t::const_iterator solverMapIter_t ;

/*
  The solver map.
*/

solverMap_t solvers ;

} // end unnamed local namespace


namespace CbcGenSolvers {

/*
  Create a vector of solver prototypes and establish a default solver.

  Creating multiple solvers is moderately expensive; if you're not interested
  in experimenting with solvers other than clp, you're likely better off just
  working with the cbc main program (CoinSolve.cpp).

  The businesss with CBC_DEFAULT_SOLVER will select the first available
  solver as the default, unless overridden at compile time.

*/

OsiSolverInterface *setupSolvers ()

{
    /*
      Populate the vector of OsiSolverInterface objects.
    */
# ifdef COIN_HAS_CLP
    solvers["clp"] = new OsiClpSolverInterface ;
# endif
# ifdef COIN_HAS_CPX
    solvers["cpx"] = new OsiCpxSolverInterface ;
# endif
# ifdef COIN_HAS_DYLP
    solvers["dylp"] = new OsiDylpSolverInterface  ;
# endif
# ifdef COIN_HAS_GLPK
    solvers["glpk"] = new OsiGlpkSolverInterface  ;
# endif
# ifdef COIN_HAS_MSK
    solvers["msk"] = new OsiMskSolverInterface  ;
# endif
# ifdef COIN_HAS_SPX
    solvers["spx"] = new OsiSpxSolverInterface  ;
# endif
    /*
      Set the standard default values in each solver.
    */
    for (solverMapIter_t solverIter = solvers.begin() ;
            solverIter != solvers.end() ;
            solverIter++) {
        OsiSolverInterface *osi = solverIter->second ;
        osi->messageHandler()->setLogLevel(0) ;
        CbcOsiParamUtils::setOsiSolverInterfaceDefaults(osi) ;
    }
    /*
      If we don't have a default solver, we're deeply confused.
    */
    OsiSolverInterface *dflt_solver = solvers[CBC_DEFAULT_SOLVER] ;
    if (dflt_solver) {
        std::cout << "Default solver is " << CBC_DEFAULT_SOLVER << std::endl ;
    } else {
        std::cerr << "No solvers!" << std::endl ;
    }

    return (dflt_solver) ;
}


/*
  Cleanup routine to delete the vector of OsiSolverInterface objects.
*/

void deleteSolvers ()
{
    for (solverMapIter_t solverIter = solvers.begin() ;
            solverIter != solvers.end() ;
            solverIter++) {
        if (solverIter->second) delete solverIter->second ;
    }
}

/*
  The `push' routine for the solver parameter.

  The basic operation is to clone the requested solver and assign it to the
  current CbcModel object.
*/

int changeCbcSolver (CoinParam *param)

{
    assert (param != 0) ;
    CbcGenParam *genParam = dynamic_cast<CbcGenParam *>(param) ;
    assert (genParam != 0) ;
    CbcGenCtlBlk *ctlBlk = genParam->obj() ;
    assert (ctlBlk != 0) ;
    CoinMessageHandler *msghandler = ctlBlk->messageHandler() ;
    /*
      Setup to return nonfatal/fatal error (1/-1) by default.
    */
    int retval ;
    if (CoinParamUtils::isInteractive()) {
        retval = 1 ;
    } else {
        retval = -1 ;
    }
    /*
      Try to locate the solver specified by the user.
    */
    const std::string solverName = genParam->kwdVal() ;
    OsiSolverInterface *protoOsi = solvers[solverName] ;
    if (protoOsi == 0) {
        std::cerr
            << "Can't find solver \"" << solverName
            << "\" in the solvers vector." << std::endl ;
        return (retval) ;
    }
    ctlBlk->dfltSolver_ = protoOsi ;
    /*
      We have a solver.
    */
    ctlBlk->message(CBCGEN_NEW_SOLVER)
    << solverName << CoinMessageEol ;
    CbcModel *model = ctlBlk->model_ ;
    assert (model != 0) ;
    OsiSolverInterface *newOsi = protoOsi->clone() ;
    model->assignSolver(newOsi) ;

    return (0) ;
}


/*
  This routine sets up a solver parameter object. It doesn't initialise the
  object being acted upon (a CbcGenCtlBlk); that's done back in the calling
  routine where we're setting up the cbc-generic parameter vector.
*/

void setupSolverParam (CbcGenParam &solverParam)

{
    /*
      Basic setup: parameter type, name, parameter code.
    */
    solverParam.setType(CoinParam::coinParamKwd) ;
    solverParam.setName("solver") ;
    solverParam.setParamCode(CbcGenParam::SOLVER) ;
    /*
      Add the solvers and set the default value.
    */
    for (solverMapIter_t solverIter = solvers.begin() ;
            solverIter != solvers.end() ;
            solverIter++) {
        solverParam.appendKwd(solverIter->first) ;
    }
    solverParam.setKwdVal(CBC_DEFAULT_SOLVER) ;
    solverParam.setDisplay(true) ;
    solverParam.setPushFunc(changeCbcSolver) ;
    /*
      And add the help strings.
    */
    solverParam.setShortHelp("Specify underlying LP solver") ;
    solverParam.setLongHelp(
        "Select the underlying LP solver that will be used to solve the continuous relaxations of subproblems."
    ) ;
}

} // end namespace CbcGenSolvers

