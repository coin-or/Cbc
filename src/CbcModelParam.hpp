/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

*/

#ifndef CbcModelParam_H
#define CbcModelParam_H

/* \file CbcModelParam.hpp
   \brief Declarations for parameters that act on a CbcModel object.
*/

/*
*/

/*! \class CbcModelParam
    \brief Class for control parameters that act on a CbcModel object.

    Adds parameter type codes and push/pull functions to the generic parameter
    object.
*/

class CbcModelParam : public CoinParam {

public:
  /*! \name Subtypes */
  //@{

  /*! \enum CbcModelParamCode
        \brief Enumeration for parameters that control a CbcModel object

      These are parameters that control the operation of a CbcModel object.
      CBCCBC_FIRSTPARAM and CBCCBC_LASTPARAM are markers to allow convenient
      separation of parameter groups.
    */
  enum CbcModelParamCode {
     CBCMODEL_FIRSTPARAM = CbcSolverParam::CBCSOLVER_LASTPARAM + 1,

     ALLOWABLEGAP,
     COSTSTRATEGY,
     CUTDEPTH,
     CUTOFF,
     CUTPASS,
     DIRECTION,
     GAPRATIO,
     INCREMENT,
     INFEASIBILITYWEIGHT,
     INTEGERTOLERANCE,
     LOGLEVEL,
     MAXIMIZE,
     MAXNODES,
     MAXNODESNOTIMPROVING, //Added
     MAXSECONDSNOTIMPROVING, //Added
     MAXSOLS, //Added
     MINIMIZE,
     MIPOPTIONS,
     MOREMIPOPTIONS,
     NUMBERANALYZE,
     NUMBERBEFORE,
     STRONGBRANCHING,
     TIMELIMIT,
     
     CBCMODEL_LASTPARAM
                 
  };

  enum OptimizationDirection {
      OptDirMaximize = 0,
      OptDirMinimize,
      OptDirZero,
      OptDireEndMarker
  };
   
  //@}

  /*! \name Constructors and Destructors

      Be careful how you specify parameters for the constructors! There's great
      potential for confusion.
    */
  //@{
  /*! \brief Default constructor */

  CbcModelParam();

  /*! \brief Constructor for a parameter with a double value

      The default value is 0.0. Be careful to clearly indicate that \p lower and
      \p upper are real (double) values to distinguish this constructor from the
      constructor for an integer parameter.
    */
  CbcModelParam(CbcModelParamCode code, std::string name, std::string help,
                double lower = -COIN_DBL_MAX, double upper = COIN_DBL_MAX,
                double defaultValue = 0.0, std::string longHelp = "",
                CoinDisplayPriority displayPriority = displayPriorityHigh);

  /*! \brief Constructor for a parameter with an integer value

      The default value is 0.
  */
  CbcModelParam(CbcModelParamCode code, std::string name, std::string help,
                int lower = -COIN_INT_MAX, int upper = COIN_INT_MAX,
                int defaultValue = 0, std::string longHelp = "",
                CoinDisplayPriority displayPriority = displayPriorityHigh);

  /*! \brief Constructor for a parameter with keyword values

      The string supplied as \p firstValue becomes the first keyword.
      Additional keywords can be added using appendKwd(). Keywords are numbered
      from zero. It's necessary to specify both the first keyword (\p
      firstValue) and the default keyword index (\p dflt) in order to
      distinguish this constructor from the string and action parameter
      constructors.
    */
  CbcModelParam(CbcModelParamCode code, std::string name, std::string help,
                std::string defaultKwd, int defaultMode, std::string longHelp = "",
                CoinDisplayPriority displayPriority = displayPriorityHigh);

  /*! \brief Constructor for a string parameter

      The default string value must be specified explicitly to distinguish
      a string constructor from an action parameter constructor.
    */

   CbcModelParam(CbcModelParamCode code, std::string name, std::string help,
                std::string defaultValue, std::string longHelp = "",
                 CoinDisplayPriority displayPriority = displayPriorityHigh);

  /*! \brief Constructor for an action parameter */

  // No defaults to resolve ambiguity
  CbcModelParam(CbcModelParamCode code, std::string name, std::string help, std::string longHelp,
                CoinDisplayPriority displayPriority);

  /*! \brief Copy constructor */

  CbcModelParam(const CbcModelParam &orig);

  /*! \brief Clone */

  CbcModelParam *clone();

  /*! \brief Assignment */

  CbcModelParam &operator=(const CbcModelParam &rhs);

  /*! \brief  Destructor */

  ~CbcModelParam();

  //@}

  /*! \name Methods to query and manipulate a parameter object */
  //@{

  /*! \brief Get the parameter code  */

  inline CbcModelParamCode paramCode() const
  {
    return (paramCode_);
  }

  /*! \brief Set the parameter code */

  inline void setParamCode(CbcModelParamCode code)
  {
    paramCode_ = code;
  }

  /*! \brief Get the underlying CbcModel object */

  inline CbcModel *obj() const
  {
    return (obj_);
  }

  /*! \brief Set the underlying CbcModel object */

  inline void setObj(CbcModel *obj)
  {
    obj_ = obj;
  }

  //@}

private:
  /*! \name Data */
  //@{

  /// Parameter code
  CbcModelParamCode paramCode_;

  /// CbcModel object
  CbcModel *obj_;

  //@}
};

/*
  Declare the utility functions.
*/

namespace CbcModelParamUtils {
   void addCbcModelParams(int &numParams, CoinParamVec &paramVec,
                          CbcModel *model);
   void loadCbcParamObj(const CoinParamVec paramVec, int first, int last,
                        CbcModel *model);
   void setCbcModelDefaults(CbcModel *model);
   
   int pushCbcModelDblParam(CoinParam &param);
   int pushCbcModelIntParam(CoinParam &param);
}

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
