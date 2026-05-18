/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

*/

#ifndef ClpParam_H
#define ClpParam_H

#include "ClpConfig.h"

#include "CoinParam.hpp"

class ClpParameters;

/* \file ClpParam.hpp
   \brief Declarations for parameters that control Clp.
*/

/*
 */

/*! \class ClpParam
    \brief Class for Clp control parameters

    Adds parameter type codes and push/pull functions to the generic parameter
    object.
*/

class CLPLIB_EXPORT ClpParam : public CoinParam {

public:

   /*! \name Enumeration types used to index parameters */
   //@{

   /*! \brief Codes to specify paramters */

   enum ClpParamCode {
      FIRSTPARAM = 0,

      // Help and Information Parameters (none of these are used in Cbc mode)
      FIRSTHELPPARAM,
      GENERALQUERY,
      FULLGENERALQUERY,
      HELP,
      LASTHELPPARAM,

      // ActionParameters
      FIRSTACTIONPARAM,
      ALLSLACK,
      BARRIER,
      BASISIN,
      BASISOUT,
      CLEARCUTS,
      DUALSIMPLEX,
      EITHERSIMPLEX,
      GUESS,
      NETLIB_BARRIER,
      NETLIB_DUAL,
      NETLIB_EITHER,
      NETLIB_PRIMAL,
      NETLIB_TUNE,
      NETWORK,
      PARAMETRICS,
      PLUSMINUS,
      PRIMALSIMPLEX,
      REALLY_SCALE,
      REVERSE,
      TIGHTEN,
      USERCLP,
      LASTCLPACTIONPARAM,
      // The remaining actions are duplicates of those in Cbc and not used
      // in Cbc mode.
      DUMMY,
      END,
      EXIT,
      STOP,
      QUIT,
      ENVIRONMENT,
      EXPORT,
      IMPORT,
      MAXIMIZE,
      MINIMIZE,
      OUTDUPROWS,
      PRINTSOL,
      READMODEL,
      READMODEL_OLD, // For backward compatibility
      READSOL,
      SOLVE,
      STATISTICS,
      UNITTEST,
      PRINTVERSION,
      WRITEGMPLSOL,
      WRITEGMPLSOL_OLD, // For backward compatibility
      WRITEMODEL,
      WRITEMODEL_OLD, // For backward compatibility
      WRITESOL,
      WRITESOL_OLD, // For backward compatibility
      WRITESOLBINARY,
      WRITESOLBINARY_OLD, // For backward compatibility
      LASTACTIONPARAM,

      // File Parameters
      FIRSTFILEPARAM,
      BASISFILE,
      PARAMETRICSFILE,
      LASTCLPFILEPARAM,
      // The remaining parameters are duplicates of those in Cbc and not used
      // in Cbc mode.
      EXPORTFILE,
      GMPLSOLFILE,
      IMPORTFILE,
      MODELFILE,
      SOLUTIONFILE,
      SOLUTIONBINARYFILE,
      LASTFILEPARAM,

      // Directory Paramaters (none of these are used in Cbc mode)
      FIRSTDIRECTORYPARAM,
      DIRECTORY,
      DIRSAMPLE,
      DIRNETLIB,
      DIRMIPLIB,
      LASTDIRECTORYPARAM,

      //String Parameters (none of these are used in Cbc mode)
      FIRSTSTRINGPARAM,
      PRINTMASK,
      LASTSTRINGPARAM,

      // On/Off Parameters
      FIRSTBOOLPARAM,
      AUTOSCALE,
      BUFFER_MODE,
      ERRORSALLOWED,
      KEEPNAMES,
      KKT,
      MESSAGES,
      PERTURBATION,
      PFI,
      SPARSEFACTOR,
      LASTBOOLPARAM,

      // Keyword Parameters
      FIRSTKWDPARAM,
      ABCWANTED,
      BARRIERSCALE,
      BIASLU,
      CHOLESKY,
      COMMANDPRINTLEVEL,
      CRASH,
      CROSSOVER,
      DIRECTION,
      DUALPIVOT,
      FACTORIZATION,
      GAMMA,
      INTPRINT,
      PRESOLVE,
      PRIMALPIVOT,
      SCALING,
      VECTOR,
      LASTKWDPARAM,

      // Integer Parameters
      FIRSTINTPARAM,
      CPP,
      DECOMPOSE_BLOCKS,
      DENSE,
      DUALIZE,
      IDIOT,
      MAXFACTOR,
      MAXITERATION,
      MORESPECIALOPTIONS,
      PERTVALUE,
      PRESOLVEPASS,
      PRESOLVEOPTIONS,
      RANDOMSEED,
      SPRINT,
      SLPVALUE,
      SMALLFACT,
      SPECIALOPTIONS,
      SUBSTITUTION,
      THREADS,
      VECTOR_MODE,
      LASTCLPINTPARAM,
      // The remaining parameters are duplicates of those in Cbc and not used
      // in Cbc mode.
      LOGLEVEL,
      OUTPUTFORMAT,
      PRINTOPTIONS,
      PROGRESSITER,
      VERBOSE,
      LASTINTPARAM,

      // Double Paramters
      FIRSTDBLPARAM,
      DUALBOUND,
      DUALTOLERANCE,
      FAKEBOUND,
      OBJSCALE,
      OBJSCALE2,
      PRESOLVETOLERANCE,
      PRIMALTOLERANCE,
      PRIMALWEIGHT,
      PROGRESS,
      PSI,
      RHSSCALE,
      TIMELIMIT,
      ZEROTOLERANCE,
      LASTDBLPARAM,

      INVALID,

      LASTPARAM
   };

   //@}

  /*! \name Constructors and Destructors

      Be careful how you specify parameters for the constructors! There's great
      potential for confusion.
    */
  //@{
  /*! \brief Default constructor */

  ClpParam();

  /*! \brief Constructor for a parameter with a double value

      The default value is 0.0. Be careful to clearly indicate that \p lower and
      \p upper are real (double) values to distinguish this constructor from the
      constructor for an integer parameter.
    */
  ClpParam(int code, std::string name, std::string help,
           double lower = -COIN_DBL_MAX, double upper = COIN_DBL_MAX,
           std::string longHelp = "",
           CoinDisplayPriority displayPriority = CoinParam::displayPriorityHigh);

  /*! \brief Constructor for a parameter with an integer value

    The default value is 0.
  */
  ClpParam(int code, std::string name, std::string help,
           int lower = -COIN_INT_MAX, int upper = COIN_INT_MAX,
           std::string longHelp = "",
           CoinDisplayPriority displayPriority = CoinParam::displayPriorityHigh);

  /*! \brief Constructor for a parameter with string values (or no value)
    Type is not optional to resolve ambiguity.

    The default value is "" for all such parameter types
  */
  ClpParam(int code, std::string name, CoinParam::CoinParamType type,
           std::string help, std::string longHelp = "",
           CoinDisplayPriority displayPriority = CoinParam::displayPriorityHigh);

  /*! \brief Copy constructor */
  ClpParam(const ClpParam &orig);

  /*! \brief Clone */
  ClpParam *clone();

  /*! \brief Assignment */
  ClpParam &operator=(const ClpParam &rhs);

  /*! \brief  Destructor */
  ~ClpParam();

  //@}

  /*! \name Methods to query and manipulate a parameter object */
  //@{

  /*! \brief Get the parameter code  */
  inline int paramCode() const { return (paramCode_); }

  /*! \brief Set the parameter code */
  inline void setParamCode(int code) { paramCode_ = code; }

  /*! \brief Get the enclosing ClpParameters object */
  inline ClpParameters *parameters() const { return (parameters_); }

  /*! \brief Set the enclosing ClpParameters object */
  inline void setParameters(ClpParameters *p) { parameters_ = p; }

  /*! \brief A hacky function to print some information about string parameters */
   std::string printString() const;

//@}

private:
  /*! \name Data */
  //@{

  /// Parameter code
  int paramCode_;

  /// Settings object
  ClpParameters *parameters_;

  //@}
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
