/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/
/*
  This file is part of cbc-generic.
*/

#ifndef CbcCbcParam_H
#define CbcCbcParam_H

/* \file CbcGenCbcParam.hpp
   \brief Declarations for parameters that act on a CbcModel object.
*/

/*
  $Id: CbcGenCbcParam.hpp 1173 2009-06-04 09:44:10Z forrest $
*/

/*! \class CbcCbcParam
    \brief Class for control parameters that act on a CbcModel object.

    Adds parameter type codes and push/pull functions to the generic parameter
    object.
*/

class CbcCbcParam : public CoinParam {

public:

    /*! \name Subtypes */
//@{

    /*! \enum CbcCbcParamCode
        \brief Enumeration for parameters that control a CbcModel object

      These are parameters that control the operation of a CbcModel object.
      CBCCBC_FIRSTPARAM and CBCCBC_LASTPARAM are markers to allow convenient
      separation of parameter groups.
    */
    typedef enum { CBCCBC_FIRSTPARAM = CbcGenParam::CBCGEN_LASTPARAM + 1,

                   ALLOWABLEGAP, COSTSTRATEGY,
                   CUTDEPTH, CUTOFF, CUTPASS, DIRECTION,
                   GAPRATIO,
                   INCREMENT, INFEASIBILITYWEIGHT, INTEGERTOLERANCE,
                   LOGLEVEL, MAXIMIZE, MAXNODES, MINIMIZE,
                   MIPOPTIONS, MOREMIPOPTIONS, NUMBERANALYZE,
                   NUMBERBEFORE, NUMBERMINI,
                   STRONGBRANCHING, TIMELIMIT_BAB,

                   CBCCBC_LASTPARAM

                 } CbcCbcParamCode ;

//@}

    /*! \name Constructors and Destructors

      Be careful how you specify parameters for the constructors! There's great
      potential for confusion.
    */
//@{
    /*! \brief Default constructor */

    CbcCbcParam() ;

    /*! \brief Constructor for a parameter with a double value

      The default value is 0.0. Be careful to clearly indicate that \p lower and
      \p upper are real (double) values to distinguish this constructor from the
      constructor for an integer parameter.
    */
    CbcCbcParam(CbcCbcParamCode code, std::string name, std::string help,
                double lower, double upper, double dflt = 0.0,
                bool display = true) ;

    /*! \brief Constructor for a parameter with an integer value

      The default value is 0.
    */
    CbcCbcParam(CbcCbcParamCode code, std::string name, std::string help,
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
    CbcCbcParam(CbcCbcParamCode code, std::string name, std::string help,
                std::string firstValue, int dflt, bool display = true) ;

    /*! \brief Constructor for a string parameter

      The default string value must be specified explicitly to distinguish
      a string constructor from an action parameter constructor.
    */

    CbcCbcParam(CbcCbcParamCode code, std::string name, std::string help,
                std::string dflt, bool display = true) ;

    /*! \brief Constructor for an action parameter */

    CbcCbcParam(CbcCbcParamCode code, std::string name, std::string help,
                bool display = true) ;

    /*! \brief Copy constructor */

    CbcCbcParam(const CbcCbcParam &orig) ;

    /*! \brief Clone */

    CbcCbcParam *clone() ;

    /*! \brief Assignment */

    CbcCbcParam &operator=(const CbcCbcParam &rhs) ;

    /*! \brief  Destructor */

    ~CbcCbcParam() ;

//@}

    /*! \name Methods to query and manipulate a parameter object */
//@{

    /*! \brief Get the parameter code  */

    inline CbcCbcParamCode paramCode() const {
        return (paramCode_) ;
    }

    /*! \brief Set the parameter code */

    inline void setParamCode(CbcCbcParamCode code) {
        paramCode_ = code ;
    }

    /*! \brief Get the underlying CbcModel object */

    inline CbcModel *obj() const {
        return (obj_) ;
    }

    /*! \brief Set the underlying CbcModel object */

    inline void setObj(CbcModel *obj) {
        obj_ = obj ;
    }

//@}

private:

    /*! \name Data */
//@{

    /// Parameter code
    CbcCbcParamCode paramCode_ ;

    /// CbcModel object
    CbcModel *obj_ ;

//@}

} ;


/*
  Declare the utility functions.
*/

namespace CbcCbcParamUtils {
void addCbcCbcParams(int &numParams, CoinParamVec &paramVec,
                     CbcModel *model) ;
void loadCbcParamObj(const CoinParamVec paramVec, int first, int last,
                     CbcModel *model) ;
void setCbcModelDefaults (CbcModel *model) ;

int pushCbcCbcDbl(CoinParam *param) ;
int pushCbcCbcInt(CoinParam *param) ;
}

#endif

