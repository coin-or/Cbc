/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/
/*
  This file is part of cbc-generic.
*/

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <string>
#include <cassert>

#include "CoinFinite.hpp"
#include "CoinParam.hpp"

#include "CbcModel.hpp"

#include "CbcGenCtlBlk.hpp"
#include "CbcGenParam.hpp"
#include "CbcGenCbcParam.hpp"

/*! \file CbcGenParamUtils
    \brief Implementation functions for CbcGenParam parameters.
*/

namespace {

char svnid[] = "$Id: CbcGenCbcParamUtils.cpp 1173 2009-06-04 09:44:10Z forrest $" ;

}

namespace CbcCbcParamUtils {


/* Function to set up cbc (CbcModel) parameters.  */

void addCbcCbcParams (int &numberParameters, CoinParamVec &parameters,
                      CbcModel *model)

{
    CbcCbcParam *param ;

    param = new CbcCbcParam(CbcCbcParam::ALLOWABLEGAP,
                            "allow!ableGap",
                            "Stop when gap between best possible and incumbent is less than this",
                            0.0, 1.0e20) ;
    param->setDblVal(0.0) ;
    param->setPushFunc(pushCbcCbcDbl) ;
    param->setObj(model) ;
    param->setLongHelp(
        "If the gap between best solution and best possible solution is less than this then the search will be terminated. Also see ratioGap."
    ) ;
    parameters.push_back(param) ;

    param = new CbcCbcParam(CbcCbcParam::CUTOFF,
                            "cuto!ff", "All solutions must be better than this", -1.0e60, 1.0e60) ;
    param->setDblVal(1.0e50) ;
    param->setPushFunc(pushCbcCbcDbl) ;
    param->setObj(model) ;
    param->setLongHelp(
        "All solutions must be better than this value (in a minimization sense).  This is also set by cbc whenever it obtains a solution and is set to the value of the objective for the solution minus the cutoff increment."
    ) ;
    parameters.push_back(param) ;

    param = new CbcCbcParam(CbcCbcParam::DIRECTION,
                            "direction", "Minimize or maximize", "min!imize", 0) ;
    param->appendKwd("max!imize") ;
    param->appendKwd("zero") ;
    param->setObj(model) ;
    param->setLongHelp(
        "The default is minimize - use 'direction maximize' for maximization.\nYou can also use the parameters 'maximize' or 'minimize'."
    ) ;
    parameters.push_back(param) ;

    param = new CbcCbcParam(CbcCbcParam::INCREMENT,
                            "inc!rement",
                            "A new solution must be at least this much better than the incumbent",
                            -1.0e20, 1.0e20, model->getDblParam(CbcModel::CbcCutoffIncrement)) ;
    param->setPushFunc(pushCbcCbcDbl) ;
    param->setObj(model) ;
    param->setLongHelp(
        "Whenever a solution is found the bound on future solutions is set to the objective of the solution (in a minimization sense) plus the specified increment.  If this option is not specified, the code will try and work out an increment.  E.g., if all objective coefficients are multiples of 0.01 and only integer variables have entries in objective then the increment can be set to 0.01.  Be careful if you set this negative!"
    ) ;
    parameters.push_back(param) ;

    param = new CbcCbcParam(CbcCbcParam::INFEASIBILITYWEIGHT,
                            "inf!easibilityWeight",
                            "Each integer infeasibility is expected to cost this much",
                            0.0, 1.0e20, model->getDblParam(CbcModel::CbcInfeasibilityWeight)) ;
    param->setPushFunc(pushCbcCbcDbl) ;
    param->setObj(model) ;
    param->setLongHelp(
        "A primitive way of deciding which node to explore next.  Satisfying each integer infeasibility is expected to cost this much."
    ) ;
    parameters.push_back(param) ;

    param = new CbcCbcParam(CbcCbcParam::INTEGERTOLERANCE,
                            "integerT!olerance",
                            "For an optimal solution, no integer variable may be farther than this from an integer value",
                            1.0e-20, 0.5, model->getDblParam(CbcModel::CbcIntegerTolerance)) ;
    param->setPushFunc(pushCbcCbcDbl) ;
    param->setObj(model) ;
    param->setLongHelp(
        "When checking a solution for feasibility, if the difference between the value of a variable and the nearest integer is less than the integer tolerance, the value is considered to be integral. Beware of setting this smaller than the primal tolerance."
    ) ;
    parameters.push_back(param) ;

    param = new CbcCbcParam(CbcCbcParam::LOGLEVEL,
                            "bclog!Level", "Level of detail in Coin branch and Cut output",
                            -1, 63, model->messageHandler()->logLevel()) ;
    param->setPushFunc(pushCbcCbcInt) ;
    param->setObj(model) ;
    param->setLongHelp(
        "If set to 0 then there should be no output in normal circumstances. A value of 1 is probably the best value for most uses, while 2 and 3 give more information."
    ) ;
    parameters.push_back(param) ;

    param = new CbcCbcParam(CbcCbcParam::MAXIMIZE,
                            "max!imize", "Set optimization direction to maximize") ;
    param->setObj(model) ;
    param->setLongHelp(
        "The default is minimize - use 'maximize' for maximization.\n A synonym for 'direction maximize'."
    ) ;
    parameters.push_back(param) ;

    param = new CbcCbcParam(CbcCbcParam::MAXNODES,
                            "maxN!odes", "Maximum number of nodes to evaluate", 1, 2147483647) ;
    param->setObj(model) ;
    param->setLongHelp(
        "This is a repeatable way to limit search.  Normally using time is easier but then the results may not be repeatable."
    ) ;
    parameters.push_back(param) ;

    param = new CbcCbcParam(CbcCbcParam::MINIMIZE,
                            "min!imize", "Set optimization direction to minimize") ;
    param->setObj(model) ;
    param->setLongHelp(
        "The default is minimize - use 'maximize' for maximization.\nThis should only be necessary if you have previously set maximization. A synonym for 'direction minimize'."
    ) ;
    parameters.push_back(param) ;

    param = new CbcCbcParam(CbcCbcParam::MIPOPTIONS,
                            "mipO!ptions", "Dubious options for mip", 0, COIN_INT_MAX, 0, false) ;
    parameters.push_back(param) ;

    param = new CbcCbcParam(CbcCbcParam::MOREMIPOPTIONS,
                            "more!MipOptions", "More dubious options for mip", -1, COIN_INT_MAX, 0, false) ;
    parameters.push_back(param) ;

    param = new CbcCbcParam(CbcCbcParam::NUMBERMINI,
                            "miniT!ree", "Size of fast mini tree", 0, COIN_INT_MAX, 0, false) ;
    param->setObj(model) ;
    param->setLongHelp(
        "The idea is that I can do a small tree fast. This is a first try and will hopefully become more sophisticated."
    ) ;
    parameters.push_back(param) ;

    param = new CbcCbcParam(CbcCbcParam::NUMBERANALYZE,
                            "numberA!nalyze",
                            "Number of analysis iterations", -COIN_INT_MAX, COIN_INT_MAX, false) ;
    param->setObj(model) ;
    param->setLongHelp(
        "This says how many iterations to spend at the root node analyzing the problem.  This is a first try and will hopefully become more sophisticated."
    ) ;
    parameters.push_back(param) ;

    param = new CbcCbcParam(CbcCbcParam::CUTPASS,
                            "passC!uts", "Number of cut passes at root node",
                            -999999, 999999, model->getMaximumCutPassesAtRoot()) ;
    param->setObj(model) ;
    param->setLongHelp(
        "The default is 100 passes if less than 500 columns, 100 passes (but stop if the drop is small) if less than 5000 columns, 20 otherwise."
    ) ;
    parameters.push_back(param) ;

    param = new CbcCbcParam(CbcCbcParam::GAPRATIO,
                            "ratio!Gap",
                            "Stop when the gap between the best possible solution and the incumbent is less than this fraction of the larger of the two",
                            0.0, 1.0e20, model->getDblParam(CbcModel::CbcAllowableFractionGap)) ;
    param->setPushFunc(pushCbcCbcDbl) ;
    param->setObj(model) ;
    param->setLongHelp(
        "If the gap between the best solution and the best possible solution is less than this fraction of the objective value at the root node then the search will terminate.  See 'allowableGap' for a way of using absolute value rather than fraction."
    ) ;
    parameters.push_back(param) ;

    param = new CbcCbcParam(CbcCbcParam::TIMELIMIT_BAB,
                            "sec!onds", "Maximum seconds for branch and cut", -1.0, 1.0e12) ;
    param->setPushFunc(pushCbcCbcDbl) ;
    param->setObj(model) ;
    param->setLongHelp(
        "After this many seconds the program will act as if maximum nodes had been reached."
    ) ;
    parameters.push_back(param) ;

    param = new CbcCbcParam(CbcCbcParam::STRONGBRANCHING,
                            "strong!Branching",
                            "Number of variables to look at in strong branching", 0, 999999,
                            model->numberStrong()) ;
    param->setObj(model) ;
    param->setLongHelp(
        "In order to decide which variable to branch on, the code will choose up to this number of unsatisfied variables and try mini up and down branches.  The most effective one is chosen. If a variable is branched on many times then the previous average up and down costs may be used - see number before trust."
    ) ;
    parameters.push_back(param) ;

    param = new CbcCbcParam(CbcCbcParam::NUMBERBEFORE,
                            "trust!PseudoCosts", "Number of branches before we trust pseudocosts",
                            -1, 2000000, model->numberBeforeTrust()) ;
    param->setObj(model) ;
    param->setPushFunc(pushCbcCbcInt) ;
    param->setLongHelp(
        "Using strong branching computes pseudo-costs.  After this many times for a variable we just trust the pseudo costs and do not do any more strong branching."
    ) ;
    parameters.push_back(param) ;

    numberParameters = parameters.size() ;

    assert (numberParameters <= parameters.capacity()) ;

}

void loadCbcParamObj (const CoinParamVec paramVec, int first, int last,
                      CbcModel *obj)

{
    int i ;

    /*
      Load the CbcModel object into the parameters
    */
    for (i = first ; i <= last ; i++) {
        CbcCbcParam *cbcParam = dynamic_cast<CbcCbcParam *>(paramVec[i]) ;
        assert (cbcParam != 0) ;
        cbcParam->setObj(obj) ;
    }

    return ;
}

/*
  Set CbcModel defaults appropriate for cbc-generic.
*/

void setCbcModelDefaults (CbcModel *model)

{
    model->setIntParam(CbcModel::CbcMaxNumNode, (COIN_INT_MAX / 2)) ;
    model->setIntParam(CbcModel::CbcMaxNumSol, 999999) ;
    model->setIntParam(CbcModel::CbcFathomDiscipline, 0) ;

    model->setDblParam(CbcModel::CbcIntegerTolerance, 1.0e-6) ;
    model->setDblParam(CbcModel::CbcInfeasibilityWeight, 0.0) ;
    model->setDblParam(CbcModel::CbcCutoffIncrement, 1.0e-5) ;
    model->setDblParam(CbcModel::CbcAllowableGap, 1.0e-10) ;
    model->setDblParam(CbcModel::CbcAllowableFractionGap, 0.0) ;
    // One year is 60x60x24x365 = 31,536,000 seconds.
    model->setDblParam(CbcModel::CbcMaximumSeconds, 3.0e7) ;
    model->setDblParam(CbcModel::CbcCurrentCutoff, 1.0e100) ;
    model->setDblParam(CbcModel::CbcOptimizationDirection, 1.0) ;
    model->setDblParam(CbcModel::CbcCurrentObjectiveValue, 1.0e100) ;
    model->setDblParam(CbcModel::CbcCurrentMinimizationObjectiveValue, 1.0e100) ;
    model->setDblParam(CbcModel::CbcStartSeconds, 0.0) ;

    model->setNumberBeforeTrust(5) ;
    model->setNumberStrong(5) ;

    return ;
}


/*
  Function to push a double parameter.
*/

int pushCbcCbcDbl (CoinParam *param)

{
    assert (param != 0) ;

    CbcCbcParam *cbcParam = dynamic_cast<CbcCbcParam *>(param) ;
    assert (cbcParam != 0) ;

    CbcModel *model = cbcParam->obj() ;
    double val = cbcParam->dblVal() ;
    CbcCbcParam::CbcCbcParamCode code = cbcParam->paramCode() ;

    assert (model != 0) ;

    int retval = 0 ;
    /*
      Translate the parameter code from CbcCbcParamCode into the correct key for
      CbcDblParam.
    */
    CbcModel::CbcDblParam key ;
    switch (code) {
    case CbcCbcParam::INTEGERTOLERANCE: {
        key = CbcModel::CbcIntegerTolerance ;
        break ;
    }
    case CbcCbcParam::INFEASIBILITYWEIGHT: {
        key = CbcModel::CbcInfeasibilityWeight ;
        break ;
    }
    case CbcCbcParam::INCREMENT: {
        key = CbcModel::CbcCutoffIncrement ;
        break ;
    }
    case CbcCbcParam::ALLOWABLEGAP: {
        key = CbcModel::CbcAllowableGap ;
        break ;
    }
    case CbcCbcParam::GAPRATIO: {
        key = CbcModel::CbcAllowableFractionGap ;
        break ;
    }
    case CbcCbcParam::TIMELIMIT_BAB: {
        key = CbcModel::CbcMaximumSeconds ;
        break ;
    }
    case CbcCbcParam::CUTOFF: {
        key = CbcModel::CbcCurrentCutoff ;
        break ;
    }
    default: {
        std::cerr << "pushCbcCbcDbl: no equivalent CbcDblParam for "
                  << "parameter code `" << code << "'." << std::endl ;
        retval = -1 ;
        break ;
    }
    }

    bool result = model->setDblParam(key, val) ;
    if (result == false) {
        retval = -1 ;
    }

    return (retval) ;
}


/*
  Function to push an integer parameter.
*/

int pushCbcCbcInt (CoinParam *param)

{
    assert (param != 0) ;

    CbcCbcParam *cbcParam = dynamic_cast<CbcCbcParam *>(param) ;
    assert (cbcParam != 0) ;

    CbcModel *model = cbcParam->obj() ;
    int val = cbcParam->intVal() ;
    CbcCbcParam::CbcCbcParamCode code = cbcParam->paramCode() ;

    assert (model != 0) ;

    int retval = 0 ;
    /*
      Translate the parameter code from CbcCbcParamCode into the correct key for
      CbcIntParam, or call the appropriate method directly.
    */
    CbcModel::CbcIntParam key = CbcModel::CbcLastIntParam ;
    switch (code) {
    case CbcCbcParam::CUTPASS: {
        model->setMaximumCutPassesAtRoot(val) ;
        break ;
    }
    case CbcCbcParam::LOGLEVEL: {
        CoinMessageHandler *hndl = model->messageHandler() ;
        assert (hndl != 0) ;
        hndl->setLogLevel(val) ;
        break ;
    }
    case CbcCbcParam::NUMBERBEFORE: {
        model->setNumberBeforeTrust(val) ;
        break ;
    }
    default: {
        std::cerr << "pushCbcCbcInt: no equivalent CbcIntParam for "
                  << "parameter code `" << code << "'." << std::endl ;
        retval = -1 ;
        break ;
    }
    }

    if (key != CbcModel::CbcLastIntParam) {
        bool result = model->setIntParam(key, val) ;
        if (result == false) {
            retval = -1 ;
        }
    }

    return (retval) ;
}

} // end namespace CbcCbcParamUtils

