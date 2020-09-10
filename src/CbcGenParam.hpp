/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

*/
/*
  This file is part of cbc-generic.
*/

#ifndef CbcSolverParam_H
#define CbcSolverParam_H

/* \file CbcSolverParam.hpp
   \brief Declarations for parameters that control the cbc-generic main
	  program.
*/

/*
*/

class CbcSolverSettings;

/*! \class CbcSolverParam
    \brief Class for cbc-generic control parameters

    Adds parameter type codes and push/pull functions to the generic parameter
    object.
*/

class CbcSolverParam : public CoinParam {

public:
  /*! \name Subtypes */
  //@{

  /*! \enum CbcSolverParamCode
        \brief Enumeration for cbc-generic parameters

      These are parameters that control the operation of the cbc-generic main
      program by operating on a CbcSolverSettings object. CBCGEN_FIRSTPARAM and
      CBCGEN_LASTPARM are markers to allow convenient separation of parameter
      groups.
    */
  typedef enum {
     CBCGEN_FIRSTPARAM = 0,
     
     GENERALQUERY,
     FULLGENERALQUERY,
     ALLCOMMANDS, //Added
     HELP,
     ARTIFICIALCOST, //Added        
     BAB,
     BKPIVOTINGSTRATEGY, //Added        
     BKMAXCALLS, //Added        
     BKCLQEXTMETHOD, //Added
     BRANCHSTRATEGY, //Not used
     CLEARCUTS,
     CLIQUECUTS,
     CLQSTRENGTHENING, //Added
     COMBINE,
     COSTSTRATEGY,
     CPP,
     CPX, //Added
     CROSSOVER, //Added
     CSVSTATISTICS, //Added
     CUTDEPTH,
     CUTLENGTH, //Added
     CUTOFFCONSTRAINT, //Added
     CUTPASSINTREE, //Added
     CUTSTRATEGY,
     DEBUG,
     DEPTHMINIBAB, //Added
     DEXTRA3, //Added
     DEXTRA4, //Added
     DEXTRA5, //Added
     DINS, //Added
     DIRECTORY,
     DIRSAMPLE, //Added
     DIRNETLIB, //Added
     DIRMIPLIB, //Added
     DIVEOPT, //Added
     DIVEOPTSOLVES, //Added
     DIVINGC, //Added
     DIVINGF, //Added
     DIVINGG, //Added
     DIVINGL, //Added
     DIVINGP, //Added
     DIVINGS, //Added
     DIVINGV, //Added
     DJFIX,
     DOHEURISTIC, //Added
     DUMMY,
     DW, //Added
     ENVIRONMENT, //Added
     ERRORSALLOWED,
     EXIT,
     EXPERIMENT, //Added
     EXPORT,
     EXTRA1, //Added 
     EXTRA2, //Added
     EXTRA3, //Added
     EXTRA4, //Added
     EXTRA_VARIABLES, //Added
     FAKECUTOFF, //Added
     FAKEINCREMENT, //Added
     FLOWCUTS,
     FPUMP,
     FPUMPITS,
     FPUMPTUNE, //Added
     FPUMPTUNE2, //Added
     GMICUTS, //Added
     GMPL_SOLUTION, //Added
     GOMORYCUTS,
     GREEDY,
     HEURISTICSTRATEGY,
     HOPTIONS, //Added
     IMPORT,
     INTPRINT,
     KNAPSACKCUTS,
     LANDPCUTS, //Added
     LAGOMORYCUTS, //Added
     LATWOMIRCUTS, //Added
     LOCALTREE,
     LOGLEVEL,
     LPLOGLEVEL, //Added
     MAXSAVEDSOLS, //Added
     MAX_SLOW_CUTS, //Added
     MESSAGES,
     MIPLIB,
     MIPSTART, //Added
     MIRDCUTS,
     MOREMOREMIPOPTIONS, //Added
     MULTIPLEROOTS, //Added
     NAIVE, //Added
     NEXTBESTSOLUTION, //Added
     NODESTRATEGY, //Added
     ODDHOLECUTS, //Not used
     ODDWHEELCUTS, //Added
     ODDWEXTMETHOD, //Added
     ORBITAL,
     OUTDUPROWS,
     OUTPUTFORMAT,
     PIVOTANDFIX, //Added
     PIVOTANDCOMPLEMENT, //Added
     PREPROCESS,
     PREPROCNAMES, //Added
     PRINTMASK,
     PRINTOPTIONS,
     PRINTVERSION,
     PRIORITYIN,
     PROBINGCUTS,
     PROCESSTUNE, //Added
     PROXIMITY, //Added
     RANDOMSEED, //Added
     RANDROUND, //Added
     REDSPLITCUTS,
     REDSPLIT2CUTS, //Added
     RENS, //Added
     RESIDCAPCUTS, //Added
     RINS, //Added
     ROUNDING,
     SAVESOL, //Added
     SHOWUNIMP,
     SMALLBAB, //Added
     SOLUTION,
     SOLVECONTINUOUS,
     SOLVER, //Not used
     SOS,
     SOSPRIORITIZE, //Added
     STATISTICS, //Added
     STDIN,
     STRATEGY, //Added
     STRENGTHEN,
     STRONGSTRATEGY, //Added
     TESTOSI, //Added
     THREADS, //Added
     TIGHTENFACTOR,
     TIMEMODE, //Added
     TWOMIRCUTS,
     UNITTEST,
     USECGRAPH, //Added
     USERCBC,
     USESOLUTION,
     VERBOSE,
     VND, //Added
     VUBTRY, //Added
     ZEROHALFCUTS, //Added
     
     CBCGEN_LASTPARAM
     
  } CbcSolverParamCode;

  //@}

  /*! \name Constructors and Destructors

      Be careful how you specify parameters for the constructors! There's great
      potential for confusion.
    */
  //@{
  /*! \brief Default constructor */

  CbcSolverParam();

  /*! \brief Constructor for a parameter with a double value

      The default value is 0.0. Be careful to clearly indicate that \p lower and
      \p upper are real (double) values to distinguish this constructor from the
      constructor for an integer parameter.
    */
  CbcSolverParam(CbcSolverParamCode code, std::string name, std::string help,
              double lower, double upper, double dflt = 0.0,
              CoinDisplayPriority displayPriority = CoinDisplayPriority::high);

  /*! \brief Constructor for a parameter with an integer value

      The default value is 0.
    */
  CbcSolverParam(CbcSolverParamCode code, std::string name, std::string help,
              int lower, int upper, int dflt = 0,
              CoinDisplayPriority displayPriority = CoinDisplayPriority::high);

  /*! \brief Constructor for a parameter with keyword values

      The string supplied as \p firstValue becomes the first keyword.
      Additional keywords can be added using appendKwd(). Keywords are numbered
      from zero. It's necessary to specify both the first keyword (\p
      firstValue) and the default keyword index (\p dflt) in order to
      distinguish this constructor from the string and action parameter
      constructors.
    */
  CbcSolverParam(CbcSolverParamCode code, std::string name, std::string help,
              std::string firstValue, int dflt,
              CoinDisplayPriority displayPriority = CoinDisplayPriority::high);

  /*! \brief Constructor for a string parameter

      The default string value must be specified explicitly to distinguish
      a string constructor from an action parameter constructor.
    */

  CbcSolverParam(CbcSolverParamCode code, std::string name, std::string help,
              std::string dflt,
              CoinDisplayPriority displayPriority = CoinDisplayPriority::high);

  /*! \brief Constructor for an action parameter */

  CbcSolverParam(CbcSolverParamCode code, std::string name, std::string help,
              CoinDisplayPriority displayPriority = CoinDisplayPriority::high);

  /*! \brief Copy constructor */

  CbcSolverParam(const CbcSolverParam &orig);

  /*! \brief Clone */

  CbcSolverParam *clone();

  /*! \brief Assignment */

  CbcSolverParam &operator=(const CbcSolverParam &rhs);

  /*! \brief  Destructor */

  ~CbcSolverParam();

  //@}

  /*! \name Methods to query and manipulate a parameter object */
  //@{

  /*! \brief Get the parameter code  */

  inline CbcSolverParamCode paramCode() const
  {
    return (paramCode_);
  }

  /*! \brief Set the parameter code */

  inline void setParamCode(CbcSolverParamCode code)
  {
    paramCode_ = code;
  }

  /*! \brief Get the underlying cbc-generic control object */

  inline CbcSolverSettings *obj() const
  {
    return (obj_);
  }

  /*! \brief Set the underlying cbc-generic control object */

  inline void setObj(CbcSolverSettings *obj)
  {
    obj_ = obj;
  }

  //@}

private:
  /*! \name Data */
  //@{

  /// Parameter code
  CbcSolverParamCode paramCode_;

  /// cbc-generic control object
  CbcSolverSettings *obj_;

  //@}
};

/*
  Declare the utility functions.
*/

namespace CbcSolverParamUtils {
void addCbcSolverParams(int &numParams, CoinParamVec &paramVec,
  CbcSettings *cbcSettings);
void loadGenParamObj(const CoinParamVec paramVec, int first, int last,
  CbcSettings *cbcSetting);

void saveSolution(const OsiSolverInterface *osi, std::string fileName);
bool readSolution(std::string fileName,
  int &numRows, int &numCols, double &objVal,
  double **rowActivity, double **dualVars,
  double **primalVars, double **reducedCosts);

int doBaCParam(CoinParam *param);
int doDebugParam(CoinParam *param);
int doExitParam(CoinParam *param);
int doHelpParam(CoinParam *param);
int doImportParam(CoinParam *param);
int doPrintMaskParam(CoinParam *param);
int doNothingParam(CoinParam *param);
int doSolutionParam(CoinParam *param);
int doUnimplementedParam(CoinParam *param);
int doVersionParam(CoinParam *param);

int pushCbcSolverDblParam(CoinParam *param);
int pushCbcSolverIntParam(CoinParam *param);
int pushCbcSolverKwdParam(CoinParam *param);
int pushCbcSolverStrParam(CoinParam *param);

int pushCbcSolverCutParam(CoinParam *param);
}

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
