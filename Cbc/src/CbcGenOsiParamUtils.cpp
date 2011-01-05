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

#include "OsiSolverInterface.hpp"

#include "CbcModel.hpp"
#include "CbcGenCtlBlk.hpp"

#include "CbcGenParam.hpp"
#include "CbcGenCbcParam.hpp"
#include "CbcGenOsiParam.hpp"

/*! \file CbcOsiParamUtils
    \brief Implementation functions for CbcOsiParam parameters.
*/

namespace {

char svnid[] = "$Id: CbcGenOsiParamUtils.cpp 1173 2009-06-04 09:44:10Z forrest $" ;

}

namespace CbcOsiParamUtils {



/*
  Function to set up OSI parameters. Does not include solver-specific
  parameters.

  ALGORITHM is commented out in CoinSolve.
*/

void addCbcOsiParams (int &numberParameters, CoinParamVec &parameters,
                      OsiSolverInterface *osi)
{
    CbcOsiParam *param ;
    OsiHintParam key ;
    bool sense ;
    OsiHintStrength strength ;
    int ival ;
    double dval ;

    param = new CbcOsiParam(CbcOsiParam::KEEPNAMES,
                            "keepN!ames", "Whether to keep row and column names on import.",
                            "off", 1) ;
    param->appendKwd("on") ;
    param->setPushFunc(pushCbcOsiKwd) ;
    param->setObj(osi) ;
    param->setLongHelp(
        "Row and column names are human-friendly, but maintaining names takes up space and time. Specifying -keepnames off >before< importing a problem will discard any name information."
    ) ;
    parameters.push_back(param) ;


    (void) osi->getIntParam(OsiMaxNumIteration, ival) ;
    param = new CbcOsiParam(CbcOsiParam::MAXITERATION,
                            "maxIt!erations", "Iteration limit for OSI solver.",
                            0, COIN_INT_MAX, ival) ;
    param->setPushFunc(pushCbcOsiInt) ;
    param->setObj(osi) ;
    param->setLongHelp(
        "Limits the number of iterations the OSI solver can perform when solving a problem."
    ) ;
    parameters.push_back(param) ;


    (void) osi->getIntParam(OsiMaxNumIterationHotStart, ival) ;
    param = new CbcOsiParam(CbcOsiParam::MAXHOTITS,
                            "hot!StartMaxIts", "Iteration limit for OSI solver hot start.",
                            0, COIN_INT_MAX, ival) ;
    param->setPushFunc(pushCbcOsiInt) ;
    param->setObj(osi) ;
    param->setLongHelp(
        "Limits the number of iterations the OSI solver can perform when solving a problem from a hot start. In the context of cbc, this limits the number of iterations expended on each LP during strong branching."
    ) ;
    parameters.push_back(param) ;

    /*
      Simplified to on/off for OsiSolverInterface, where it goes in as a hint.
    */
    (void) osi->getHintParam(OsiDoPresolveInInitial, sense, strength) ;
    if (sense == true) {
        ival = 1 ;
    } else {
        ival = 0 ;
    }
    param = new CbcOsiParam(CbcOsiParam::PRESOLVE,
                            "presolve", "Whether to presolve problem", "off", ival) ;
    param->appendKwd("on") ;
    param->setPushFunc(pushCbcOsiHint) ;
    param->setObj(osi) ;
    param->setLongHelp(
        "Presolve analyzes the model to find such things as redundant constraints, constraints which fix some variables, constraints which can be transformed into bounds, etc.  For the initial solve of any problem this is worth doing unless you know that it will have no effect."
    ) ;
    parameters.push_back(param) ;


    param = new CbcOsiParam(CbcOsiParam::PRIMALTOLERANCE,
                            "primalT!olerance",
                            "For an optimal solution no primal infeasibility may exceed this value",
                            1.0e-20, 1.0e12) ;
    param->setPushFunc(pushCbcOsiDbl) ;
    param->setObj(osi) ;
    param ->setLongHelp(
        "Normally the default tolerance is fine, but you may want to increase it a bit if a primal run seems to be having a hard time"
    ) ;
    parameters.push_back(param) ;

    /*
      Simplified for OsiSolverInterface, which just takes a hint.
    */
    (void) osi->getHintParam(OsiDoScale, sense, strength) ;
    if (sense == true) {
        ival = 1 ;
    } else {
        ival = 0 ;
    }
    param = new CbcOsiParam(CbcOsiParam::SCALING,
                            "scal!ing", "Whether to scale problem", "off", ival) ;
    param ->appendKwd("on") ;
    param->setPushFunc(pushCbcOsiHint) ;
    param->setObj(osi) ;
    param ->setLongHelp(
        "Scaling can help in solving problems which might otherwise fail because of lack of accuracy.  It can also reduce the number of iterations.  It is not applied if the range of elements is small.  When unscaled it is possible that there may be small primal and/or infeasibilities."
    ) ;
    parameters.push_back(param) ;

    ival = osi->messageHandler()->logLevel() ;
    param = new CbcOsiParam(CbcOsiParam::SOLVERLOGLEVEL,
                            "slog!Level", "Level of detail in Solver output", -1, 63, ival) ;
    param->setPushFunc(pushCbcOsiLogLevel) ;
    param->setObj(osi) ;
    param ->setLongHelp(
        "If 0 then there should be no output in normal circumstances. 1 is probably the best value for most uses, while 2 and 3 give more information."
    ) ;
    parameters.push_back(param) ;

    numberParameters = parameters.size() ;
    assert (numberParameters <= parameters.capacity()) ;

}

void loadOsiParamObj (const CoinParamVec paramVec, int first, int last,
                      OsiSolverInterface *obj)

{
    int i ;
    /*
      Load the OsiSolverInterface object into the parameters
    */
    for (i = first ; i <= last ; i++) {
        CbcOsiParam *osiParam = dynamic_cast<CbcOsiParam *>(paramVec[i]) ;
        assert (osiParam != 0) ;
        osiParam->setObj(obj) ;
    }

    return ;
}


/*
  Function to set default values for solver appropriate for cbc-generic.
*/

void setOsiSolverInterfaceDefaults (OsiSolverInterface *osi)

{
    bool result ;

    /*
      OsiNameDiscipline isn't supported by all solvers, so check to see that it
      worked. If not, fall back to zero.
    */
    osi->setIntParam(OsiMaxNumIteration, 1000000) ;
    osi->setIntParam(OsiMaxNumIterationHotStart, 1000) ;
    result = osi->setIntParam(OsiNameDiscipline, 1) ;
    if (!result) {
        result = osi->setIntParam(OsiNameDiscipline, 0) ;
    }

    /*
      Primal and dual feasibility tolerances (OsiPrimalTolerance and
      OsiDualTolerance, respectively)  are left to the discretion of the solver.
    */
    osi->setDblParam(OsiDualObjectiveLimit, 1.0e100) ;
    osi->setDblParam(OsiPrimalObjectiveLimit, 1.0e100) ;
    osi->setDblParam(OsiObjOffset, 0.0) ;

    osi->setHintParam(OsiDoPresolveInInitial, true, OsiHintDo) ;
    osi->setHintParam(OsiDoDualInInitial, true, OsiHintIgnore) ;
    osi->setHintParam(OsiDoPresolveInResolve, false, OsiHintTry) ;
    osi->setHintParam(OsiDoDualInInitial, true, OsiHintTry) ;
    osi->setHintParam(OsiDoScale, true, OsiHintDo) ;
    osi->setHintParam(OsiDoCrash, true, OsiHintIgnore) ;
    osi->setHintParam(OsiDoReducePrint, true, OsiHintDo) ;
    osi->setHintParam(OsiDoInBranchAndCut, true, OsiHintTry) ;

    return ;
}


/*
  Function to push an integer parameter.
*/

int pushCbcOsiInt (CoinParam *param)

{
    assert (param != 0) ;

    CbcOsiParam *osiParam = dynamic_cast<CbcOsiParam *>(param) ;
    assert (osiParam != 0) ;

    OsiSolverInterface *osi = osiParam->obj() ;
    assert (osi != 0) ;
    int val = osiParam->intVal() ;
    CbcOsiParam::CbcOsiParamCode code = osiParam->paramCode() ;

    assert (osi != 0) ;
    /*
      Setup to return nonfatal/fatal error (1/-1) by default, so that all we need
      to do is correct to 0 (no error) if we're successful.
    */
    int retval ;
    if (CoinParamUtils::isInteractive()) {
        retval = 1 ;
    } else {
        retval = -1 ;
    }
    /*
      Translate the parameter code from CbcOsiParamCode into the correct key for
      CbcIntParam.
    */
    OsiIntParam key ;
    switch (code) {
    case CbcOsiParam::MAXITERATION: {
        key = OsiMaxNumIteration ;
        break ;
    }
    case CbcOsiParam::MAXHOTITS: {
        key = OsiMaxNumIterationHotStart ;
        break ;
    }
    default: {
        std::cerr << "pushCbcOsiIntParam: no equivalent OsiIntParam for "
                  << "parameter code `" << code << "'." << std::endl ;
        retval = -1 ;
        break ;
    }
    }

    bool setOK = osi->setIntParam(key, val) ;
    if (setOK == false) {
        retval = -1 ;
    }

    return (retval) ;
}
/*
  Function to push a double parameter.
*/

int pushCbcOsiDbl (CoinParam *param)

{
    assert (param != 0) ;

    CbcOsiParam *osiParam = dynamic_cast<CbcOsiParam *>(param) ;
    assert (osiParam != 0) ;

    OsiSolverInterface *osi = osiParam->obj() ;
    assert (osi != 0) ;
    double val = osiParam->dblVal() ;
    CbcOsiParam::CbcOsiParamCode code = osiParam->paramCode() ;

    assert (osi != 0) ;
    /*
      Setup to return nonfatal/fatal error (1/-1) by default, so that all we need
      to do is correct to 0 (no error) if we're successful.
    */
    int retval ;
    if (CoinParamUtils::isInteractive()) {
        retval = 1 ;
    } else {
        retval = -1 ;
    }
    /*
      Translate the parameter code from CbcOsiParamCode into the correct key for
      CbcDblParam.
    */
    OsiDblParam key ;
    switch (code) {
    case CbcOsiParam::PRIMALTOLERANCE: {
        key = OsiPrimalTolerance ;
        break ;
    }
    case CbcOsiParam::DUALTOLERANCE: {
        key = OsiDualTolerance ; ;
        break ;
    }
    case CbcOsiParam::DUALBOUND: {
        key = OsiDualObjectiveLimit ;
        break ;
    }
    default: {
        std::cerr << "pushCbcOsiDblParam: no equivalent OsiDblParam for "
                  << "parameter code `" << code << "'." << std::endl ;
        retval = -1 ;
        break ;
    }
    }

    bool setOK = osi->setDblParam(key, val) ;
    if (setOK == false) {
        retval = -1 ;
    }

    return (retval) ;
}


/*
  Function to push a keyword-valued parameter. This can translate into integer
  as well as string-valued parameters.
*/

int pushCbcOsiKwd (CoinParam *param)

{
    assert (param != 0) ;
    CbcOsiParam *osiParam = dynamic_cast<CbcOsiParam *>(param) ;
    assert (osiParam != 0) ;
    OsiSolverInterface *osi = osiParam->obj() ;
    assert (osi != 0) ;

    std::string str = osiParam->kwdVal() ;
    CbcOsiParam::CbcOsiParamCode code = osiParam->paramCode() ;

    int retval = 0 ;
    /*
      Figure out what we're doing and set the relevant field.
    */
    OsiIntParam key ;

    switch (code) {
    case CbcOsiParam::KEEPNAMES: {
        if (str == "on" || str == "off") {
            int discipline ;
            if (str == "on") {
                discipline = 1 ;
            } else {
                discipline = 0 ;
            }
            bool recog = osi->setIntParam(OsiNameDiscipline, discipline) ;
            if (recog == false) {
                std::cerr
                    << "pushCbcOsiKwdParam(KEEPNAMES): underlying solver does not "
                    << "recognise name discipline " << discipline << "."
                    << std::endl ;
                retval = +1 ;
            }
        } else {
            std::cerr
                << "pushCbcOsiKwdParam(KEEPNAMES): unrecognised keyword `"
                << str << "'." << std::endl ;
            retval = -1 ;
        }
        break ;
    }
    default: {
        std::cerr
            << "pushCbcGenKwdParam: unrecognised parameter code `"
            << code << "'." << std::endl ;
        retval = -1 ;
        break ;
    }
    }

    return (retval) ;
}


/*
  Function to set the solver's output level. To cover all the bases, adjust
  the message handler and set the hint. Nothing can go fatally wrong here,
  but we'll return non-fatal error if the solver rejects the hint. The
  implementor of an OSI has wide latitude with hints, and may elect to set a
  log level as part of handling the hint. Do that first and then explicitly
  set the message handler log level to be sure the new value isn't
  overridden.
*/

int pushCbcOsiLogLevel (CoinParam *param)

{
    assert (param != 0) ;
    CbcOsiParam *osiParam = dynamic_cast<CbcOsiParam *>(param) ;
    assert (osiParam != 0) ;
    OsiSolverInterface *osi = osiParam->obj() ;
    assert(osi != 0) ;

    int lvl = param->intVal() ;
    /*
      Setup to return nonfatal/fatal error (1/-1) by default, so that all we need
      to do is correct to 0 (no error) if we're successful.
    */
    int retval ;
    if (CoinParamUtils::isInteractive()) {
        retval = 1 ;
    } else {
        retval = -1 ;
    }
    /*
      Now try to do the right thing with a hint. Harder to say -- assume that log
      level 1 is `normal'.
    */
    OsiHintStrength strength ;
    bool sense ;
    if (lvl < 1) {
        strength = OsiHintDo ;
        sense = true ;
    } else if (lvl == 1) {
        strength = OsiHintIgnore ;
        sense = true ;
    } else if (lvl == 2) {
        strength = OsiHintTry ;
        sense = false ;
    } else {
        strength = OsiHintDo ;
        sense = false ;
    }

    bool setOK = osi->setHintParam(OsiDoReducePrint, sense, strength) ;

    /*
      Recover the message handler and set the log level directly.
    */
    CoinMessageHandler *hndl = osi->messageHandler() ;
    assert (hndl != 0) ;
    hndl->setLogLevel(lvl) ;

    if (setOK) {
        return (0) ;
    } else {
        return (retval) ;
    }
}


/*
  Function for parameters that are enabled/disabled with a hint.
*/
int pushCbcOsiHint (CoinParam *param)

{
    assert (param != 0) ;
    CbcOsiParam *osiParam = dynamic_cast<CbcOsiParam *>(param) ;
    assert (osiParam != 0) ;
    OsiSolverInterface *osi = osiParam->obj() ;
    assert(osi != 0) ;
    /*
      Setup to return nonfatal/fatal error (1/-1) by default, so that all we need
      to do is correct to 0 (no error) if we're successful.
    */
    int retval ;
    if (CoinParamUtils::isInteractive()) {
        retval = 1 ;
    } else {
        retval = -1 ;
    }
    /*
      Set the sense for the hint.
    */
    std::string kwd = param->kwdVal() ;
    bool sense ;
    if (kwd == "off") {
        sense = false ;
    } else {
        sense = true ;
    }
    /*
      Grab the parameter code and translate to an OSI parameter key.
    */
    CbcOsiParam::CbcOsiParamCode code = osiParam->paramCode() ;
    OsiHintParam key ;
    switch (code) {
    case CbcOsiParam::PRESOLVE: {
        key = OsiDoPresolveInInitial ;
        break ;
    }
    case CbcOsiParam::SCALING: {
        key = OsiDoScale ;
        break ;
    }
    default: {
        std::cerr << "pushCbcOsiHint: no equivalent OsiHintParam for "
                  << "parameter code `" << code << "'." << std::endl ;
        retval = -1 ;
        break ;
    }
    }

    bool setOK = osi->setHintParam(key, sense, OsiHintDo) ;

    if (setOK) {
        return (0) ;
    } else {
        return (retval) ;
    }
}

} // end namespace CbcOsiParamUtils

