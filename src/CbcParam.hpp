/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

*/

#ifndef CbcParam_H
#define CbcParam_H

#include "CoinParam.hpp"

/* \file CbcParam.hpp
   \brief Declarations for parameters that control the cbc-generic main
          program.
*/

/*
 */

class CbcModel;
class CbcSettings;

/*! \class CbcParam
    \brief Class for cbc-generic control parameters

    Adds parameter type codes and push/pull functions to the generic parameter
    object.
*/

class CbcParam : public CoinParam {

public:

  /*! \name Constructors and Destructors

      Be careful how you specify parameters for the constructors! There's great
      potential for confusion.
    */
  //@{
  /*! \brief Default constructor */

  CbcParam();

  /*! \brief Constructor for a parameter with a double value

      The default value is 0.0. Be careful to clearly indicate that \p lower and
      \p upper are real (double) values to distinguish this constructor from the
      constructor for an integer parameter.
    */
  CbcParam(int code, std::string name, std::string help,
           double lower = -COIN_DBL_MAX, double upper = COIN_DBL_MAX,
           double defaultValue = 0.0, std::string longHelp = "",
           CoinDisplayPriority display = displayPriorityHigh);

  /*! \brief Constructor for a parameter with an integer value

      The default value is 0.
    */
  CbcParam(int code, std::string name, std::string help,
           int lower = -COIN_INT_MAX, int upper = COIN_INT_MAX,
           int defaultValue = 0, std::string longHelp = "",
           CoinDisplayPriority display = displayPriorityHigh);

  /*! \brief Constructor for a parameter with keyword values

      The string supplied as \p firstValue becomes the first keyword.
      Additional keywords can be added using appendKwd(). Keywords are numbered
      from zero. It's necessary to specify both the first keyword (\p
      firstValue) and the default keyword index (\p dflt) in order to
      distinguish this constructor from the string and action parameter
      constructors.
    */
  CbcParam(int code, std::string name, std::string help,
           std::string defaultKwd, int defaultMode,
           std::string longHelp = "",
           CoinDisplayPriority display = displayPriorityHigh);

  /*! \brief Constructor for a string parameter

      The default string value must be specified explicitly to distinguish
      a string constructor from an action parameter constructor.
    */
  CbcParam(int code, std::string name, std::string help,
           std::string defaultValue, std::string longHelp = "",
           CoinDisplayPriority display = displayPriorityHigh);

  /*! \brief Constructor for an action parameter */
  // No defaults to resolve ambiguity
  CbcParam(int code, std::string name, std::string help,
           std::string longHelp, CoinDisplayPriority display);

  /*! \brief Copy constructor */
  CbcParam(const CbcParam &orig);

  /*! \brief Clone */
  CbcParam *clone();

  /*! \brief Assignment */
  CbcParam &operator=(const CbcParam &rhs);

  /*! \brief  Destructor */
  ~CbcParam();

  //@}

  /*! \name Methods to query and manipulate a parameter object */
  //@{

  /*! \brief Get the parameter code  */
  inline int paramCode() const { return (paramCode_); }

  /*! \brief Set the parameter code */
  inline void setParamCode(int code) { paramCode_ = code; }

  /*! \brief Get the associated Cbc settings object */
  inline CbcSettings *settings() const { return (settings_); }

  /*! \brief Set the associated Cbc settings object */
  inline void setSettings(CbcSettings *s) { settings_ = s; }

  /*! \brief Get the associated model */
  inline CbcModel *model() const { return (model_); }

  /*! \brief Set the associated Cbc settings object */
  inline void setModel(CbcModel *m) { model_ = m; }

  /*! \brief A hacky function to print some information about string parameters */
  std::string printString() const;

//@}

private:
  /*! \name Data */
  //@{

  /// Parameter code
  int paramCode_;

  /// Pointer to settings object
  CbcSettings *settings_;

  /// Pointer to model object
  CbcModel * model_;

  //@}
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
