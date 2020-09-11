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
     CBCSOLVER_FIRSTPARAM = 0,
     
     GENERALQUERY,
     FULLGENERALQUERY,
     ALLCOMMANDS,
     HELP,
     ARTIFICIALCOST,        
     BAB,
     BKPIVOTINGSTRATEGY,        
     BKMAXCALLS,        
     BKCLQEXTMETHOD,
     BRANCHPRIORITY,
     BRANCHSTRATEGY, //Not used
     CLEARCUTS,
     CLIQUECUTS,
     CLQSTRENGTHENING,
     COMBINE,
     CPP,
     CPX,
     CROSSOVER,
     CSVSTATISTICS,
     CUTDEPTH,
     CUTLENGTH,
     CUTOFFCONSTRAINT,
     CUTPASSINTREE,
     CUTSTRATEGY,
     DEBUG,
     DEPTHMINIBAB,
     DEXTRA3,
     DEXTRA4,
     DEXTRA5,
     DINS,
     DIRECTORY,
     DIRSAMPLE,
     DIRNETLIB,
     DIRMIPLIB,
     DIVEOPT,
     DIVEOPTSOLVES,
     DIVINGC,
     DIVINGF,
     DIVINGG,
     DIVINGL,
     DIVINGP,
     DIVINGS,
     DIVINGV,
     DJFIX,
     DOHEURISTIC,
     DUMMY,
     DW,
     ENVIRONMENT,
     ERRORSALLOWED,
     EXIT,
     EXPERIMENT,
     EXPORT,
     EXTRA1, 
     EXTRA2,
     EXTRA3,
     EXTRA4,
     EXTRA_VARIABLES,
     FAKECUTOFF,
     FAKEINCREMENT,
     FLOWCUTS,
     FPUMP,
     FPUMPITS,
     FPUMPTUNE,
     FPUMPTUNE2,
     GMICUTS,
     GMPL_SOLUTION,
     GOMORYCUTS,
     GREEDY,
     HEURISTICSTRATEGY,
     HEUROPTIONS,
     IMPORT,
     INTPRINT,
     KNAPSACKCUTS,
     LANDPCUTS,
     LAGOMORYCUTS,
     LATWOMIRCUTS,
     LOCALTREE,
     LOGLEVEL,
     LPLOGLEVEL,
     MAXSAVEDSOLS,
     MAXSLOWCUTS,
     MESSAGES,
     MIPLIB,
     MIPSTART,
     MIRCUTS,
     MOREMOREMIPOPTIONS,
     MULTIPLEROOTS,
     NAIVE,
     NEXTBESTSOLUTION,
     NODESTRATEGY,
     ODDHOLECUTS, //Not used
     ODDWHEELCUTS,
     ODDWEXTMETHOD,
     ORBITAL,
     OUTDUPROWS,
     OUTPUTFORMAT,
     PIVOTANDFIX,
     PIVOTANDCOMPLEMENT,
     PREPROCESS,
     PREPROCNAMES,
     PRINTMASK,
     PRINTOPTIONS,
     PRINTVERSION,
     PRIORITYIN,
     PROBINGCUTS,
     PROCESSTUNE,
     PROXIMITY,
     RANDOMSEED,
     RANDROUND,
     REDSPLITCUTS,
     REDSPLIT2CUTS,
     RENS,
     RESIDCAPCUTS,
     RINS,
     ROUNDING,
     SAVESOL,
     SHOWUNIMP,
     SMALLBAB,
     SOLUTION,
     SOLVECONTINUOUS,
     SOLVER, //Not used
     SOS,
     SOSPRIORITIZE,
     STATISTICS,
     STDIN,
     STRATEGY,
     STRENGTHEN,
     STRONGSTRATEGY,
     TESTOSI,
     THREADS,
     TIGHTENFACTOR,
     TIMEMODE,
     TWOMIRCUTS,
     UNITTEST,
     USECGRAPH,
     USERCBC,
     USESOLUTION,
     VERBOSE,
     VND,
     VUBTRY,
     ZEROHALFCUTS,
     
     CBCSOLVER_LASTPARAM
     
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
                 CoinDisplayPriority display = displayPriorityHigh);

  /*! \brief Constructor for a parameter with an integer value

      The default value is 0.
    */
  CbcSolverParam(CbcSolverParamCode code, std::string name, std::string help,
                 int lower, int upper, int dflt = 0,
                 CoinDisplayPriority display = displayPriorityHigh);

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
                 CoinDisplayPriority display = displayPriorityHigh);

  /*! \brief Constructor for a string parameter

      The default string value must be specified explicitly to distinguish
      a string constructor from an action parameter constructor.
    */

  CbcSolverParam(CbcSolverParamCode code, std::string name, std::string help,
                 std::string dflt,
                 CoinDisplayPriority display = displayPriorityHigh);

  /*! \brief Constructor for an action parameter */

  CbcSolverParam(CbcSolverParamCode code, std::string name, std::string help,
                 CoinDisplayPriority display = displayPriorityHigh);

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
                        CbcSolverSettings *cbcSettings);
void loadGenParamObj(const CoinParamVec paramVec, int first, int last,
                     CbcSolverSettings *cbcSetting);

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
int pushCbcSolverBoolParam(CoinParam *param);

int pushCbcSolverHeurParam(CoinParam *param);
int pushCbcSolverCutParam(CoinParam *param);
}

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
