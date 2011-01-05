/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/
/*
  This file is part of cbc-generic.
*/

#ifndef CbcOsiParam_H
#define CbcOsiParam_H

/* \file CbcGenOsiParam.hpp
   \brief Declarations for parameters that act on a OsiSolverInterface object.
*/

/*
  $Id: CbcGenOsiParam.hpp 1173 2009-06-04 09:44:10Z forrest $
*/

/*! \class CbcOsiParam
    \brief Class for control parameters that act on a OsiSolverInterface object.

    Adds parameter type codes and push/pull functions to the generic parameter
    object.
*/

class CbcOsiParam : public CoinParam {

public:

    /*! \name Subtypes */
//@{

    /*! \enum CbcOsiParamCode
        \brief Enumeration for parameters that control an OsiSolverInterface
           object

      These are parameters that control the operation of an OsiSolverInterface
      object. CBCOSI_FIRSTPARAM and CBCOSI_LASTPARAM are markers to allow
      convenient separation of parameter groups.
    */
    typedef enum { CBCOSI_FIRSTPARAM = CbcCbcParam::CBCCBC_LASTPARAM + 1,

                   ALGORITHM, ALLSLACK, AUTOSCALE, BARRIER, BARRIERSCALE,
                   BASISIN, BASISOUT, BIASLU, CHOLESKY, CRASH, CROSSOVER,
                   DUALBOUND, DUALPIVOT, DUALSIMPLEX, DUALTOLERANCE, FAKEBOUND,
                   GAMMA, IDIOT, KEEPNAMES, KKT, MAXITERATION, MAXHOTITS, NETLIB_BARRIER,
                   NETLIB_DUAL, NETLIB_PRIMAL, NETWORK, OBJSCALE, PERTURBATION,
                   PERTVALUE, PFI, PLUSMINUS, PRESOLVE, PRESOLVEOPTIONS, PRESOLVEPASS,
                   PRIMALPIVOT, PRIMALSIMPLEX, PRIMALTOLERANCE, REALLY_SCALE,
                   RESTORE, REVERSE, RHSSCALE, SAVE, SCALING, SLPVALUE, SOLVERLOGLEVEL,
                   SPARSEFACTOR, SPECIALOPTIONS, SPRINT, TIGHTEN,

                   CBCOSI_LASTPARAM

                 } CbcOsiParamCode ;

//@}

    /*! \name Constructors and Destructors

      Be careful how you specify parameters for the constructors! There's great
      potential for confusion.
    */
//@{
    /*! \brief Default constructor */

    CbcOsiParam() ;

    /*! \brief Constructor for a parameter with a double value

      The default value is 0.0. Be careful to clearly indicate that \p lower and
      \p upper are real (double) values to distinguish this constructor from the
      constructor for an integer parameter.
    */
    CbcOsiParam(CbcOsiParamCode code, std::string name, std::string help,
                double lower, double upper, double dflt = 0.0,
                bool display = true) ;

    /*! \brief Constructor for a parameter with an integer value

      The default value is 0.
    */
    CbcOsiParam(CbcOsiParamCode code, std::string name, std::string help,
                int lower, int upper, int dflt = 0,
                bool display = true) ;

    /*! \brief Constructor for a parameter with keyword values

      The string supplied as \p firstValue becomes the first keyword.
      Additional keywords can be added using appendKwd(). Keywords are numbered
      from zero. It's necessary to specify both the first keyword (\p
      firstValue) and the default keyword index (\p dflt) in order to
      distinguish this constructor from the string and action parameter
      constructors.
    */
    CbcOsiParam(CbcOsiParamCode code, std::string name, std::string help,
                std::string firstValue, int dflt, bool display = true) ;

    /*! \brief Constructor for a string parameter

      The default string value must be specified explicitly to distinguish
      a string constructor from an action parameter constructor.
    */

    CbcOsiParam(CbcOsiParamCode code, std::string name, std::string help,
                std::string dflt, bool display = true) ;

    /*! \brief Constructor for an action parameter */

    CbcOsiParam(CbcOsiParamCode code, std::string name, std::string help,
                bool display = true) ;

    /*! \brief Copy constructor */

    CbcOsiParam(const CbcOsiParam &orig) ;

    /*! \brief Clone */

    CbcOsiParam *clone() ;

    /*! \brief Assignment */

    CbcOsiParam &operator=(const CbcOsiParam &rhs) ;

    /*! \brief  Destructor */

    ~CbcOsiParam() ;

//@}

    /*! \name Methods to query and manipulate a parameter object */
//@{

    /*! \brief Get the parameter code  */

    inline CbcOsiParamCode paramCode() const {
        return (paramCode_) ;
    }

    /*! \brief Set the parameter code */

    inline void setParamCode(CbcOsiParamCode code) {
        paramCode_ = code ;
    }

    /*! \brief Get the underlying OsiSolverInterface object */

    inline OsiSolverInterface *obj() const {
        return (obj_) ;
    }

    /*! \brief Set the underlying OsiSolverInterace object */

    inline void setObj(OsiSolverInterface *obj) {
        obj_ = obj ;
    }

//@}


private:

    /*! \name Data */
//@{

    /// Parameter code
    CbcOsiParamCode paramCode_ ;

    /// OsiSolverInterface object
    OsiSolverInterface *obj_ ;

//@}

} ;



/*
  Declare the utility functions.
*/

namespace CbcOsiParamUtils {
void addCbcOsiParams(int &numParams, CoinParamVec &paramVec,
                     OsiSolverInterface *osi) ;
void loadOsiParamObj(const CoinParamVec paramVec,
                     CbcGenCtlBlk *ctlBlk) ;
void setOsiSolverInterfaceDefaults(OsiSolverInterface *osi) ;

int pushCbcOsiLogLevel(CoinParam *param) ;
int pushCbcOsiInt(CoinParam *param) ;
int pushCbcOsiDbl(CoinParam *param) ;
int pushCbcOsiKwd(CoinParam *param) ;
int pushCbcOsiHint(CoinParam *param) ;
}


#endif

