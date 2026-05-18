/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

*/
#ifndef ClpParameters_H
#define ClpParameters_H

/* \file ClpParameters.hpp
   \brief Declarations for parameters of Clp.
*/

#include "ClpConfig.h"

#include "ClpParam.hpp"
#include "ClpSimplex.hpp"

// For now, we need this by default. It will be removed later.
#define CBC_CLUMSY_CODING

/* \brief Clp algorithm control class

  This class defines and stores the parameters used to control the operation
  of Clp.
*/

class CLPLIB_EXPORT ClpParameters {

public:
  /*! \name Constructors and destructors */
  //@{

  /*! \brief Constructors */
  ClpParameters(bool cbcMode = false);

  ClpParameters(int strategy, bool cbcMode = false);

  void init(int strategy);

  /*! \brief Destructor */
  ~ClpParameters();

  /*! \brief Copy constructor */
  ClpParameters(const ClpParameters &rhs);

  /*! \brief Assignment operator (deep-copies parameter vector) */
  ClpParameters &operator=(const ClpParameters &rhs);

  /*! \name Enumeration types used for Clp keyword parameters */
  //@{

   /*! \brief Codes to specify overall strategies */

   enum ClpStrategy { DefaultStrategy = 0 };

  /*! \brief Codes to specify one or off for binary parameters

     - ParamOff: Capability is switched off
     - ParamOn: Capability is switched on
   */

  enum OnOffMode { ParamOff = 0, ParamOn, ParamEndMarker };

  /*! \brief What parameters to print

    - displayAll:
    - displayLowHigh:
    - displayHigh:
    - displayEndMarker

   */

  enum CommandDisplayMode {
     displayAll = 0,
     displayLowHigh,
     displayHigh,
     displayEndMarker
  };
  //@}

  /*! \name Operators
      \brief Functions that define operators to allow access by [].
  */

  //@{

  ClpParam *operator[](std::size_t idx){
     return getParam(idx);
  }

  /*! \name Functions for Setting Up Parameters
      \brief Functions that populate the parameters objects.
  */

  //@{

  /*! set up the solver parameter vector */
  void addClpParams();
  void addClpStrParams();
  void addClpDirParams();
  void addClpFileParams();
  void addClpHelpParams();
  void addClpActionParams();
  void addClpKwdParams();
  void addClpDblParams();
  void addClpIntParams();
  void addClpBoolParams();

  /*! set up the default */
  void setDefaults(int strategy);

  //@{

  /*! \name Access functions
      \brief Functions that get and set data.
  */

  //@{
  /* \brief Get Clp solver parameter vector */
  inline CoinParamVec &paramVec() { return parameters_; }

  /* \brief Get specific Clp solver parameter object */
  inline ClpParam *getParam(int code) {
     return static_cast<ClpParam *>(parameters_[code]);
  }

  /* \brief Get value of parameter */
  void getParamVal(int code, std::string &value) {
     value = parameters_[code]->getVal(value);
  }
  void getParamVal(int code, double &value) {
     value = parameters_[code]->getVal(value);
  }
  void getParamVal(int code, int &value) {
     value = parameters_[code]->getVal(value);
  }

  /* \brief Set value of parameter */
   void setParamVal(int code, std::string value,
                    std::string *message = NULL,
                    CoinParam::ParamPushMode pMode = CoinParam::pushDefault) {
      parameters_[code]->setVal(value, message, pMode);
  }
   void setParamVal(int code, double value,
                    std::string *message = NULL,
                    CoinParam::ParamPushMode pMode = CoinParam::pushDefault) {
      parameters_[code]->setVal(value, message, pMode);
  }
   void setParamVal(int code, int value,
                    std::string *message = NULL,
                    CoinParam::ParamPushMode pMode = CoinParam::pushDefault) {
      parameters_[code]->setVal(value, message, pMode);
  }

  /* \brief Get version */
  inline std::string getVersion() { return CLP_VERSION; }

  /* \brief Get default directory */
  inline std::string getDefaultDirectory() { return dfltDirectory_; }

  /* \brief Set default directory */
   inline void setDefaultDirectory(std::string dir) { dfltDirectory_ = dir; }

  /*! \brief Set Clp model */
  inline void setModel(ClpSimplex *model) { model_ = model; }

  /*! \brief Get Clp model */
  inline ClpSimplex *getModel() const { return (model_); }

  /*! \brief Say whether in CbcMode */
  inline void setCbcMode(bool yesNo) { cbcMode_ = yesNo; }
#ifdef CBC_CLUMSY_CODING
  /*! \brief Synchronize Clp model - Int and Dbl */
  void synchronizeModel();

#endif
  //@{

   /*! \brief Returns index of first parameter that matches and number of
     matches overall. Returns CLP_INVALID if no match */
   int matches(std::string field, int &numberMatches);

private:

  /*! \brief Default directory prefix */
  std::string dfltDirectory_;

  /*! \brief The Cbc parameter vector (parameters stored by their index) */
  CoinParamVec parameters_;

  /*! \brief A pointer to the current ClpSimplex object */
  ClpSimplex *model_;

  /*! \brief To say if in Cbc mode */
  bool cbcMode_;

};

#endif
