/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

*/

#ifndef CbcParam_H
#define CbcParam_H

#include "CoinParam.hpp"
#include "CbcConfig.h"

/* \file CbcParam.hpp
   \brief Declarations for parameters that control the cbc-generic main
          program.
*/

/*
 */

class CbcModel;
class CbcParameters;

/*! \class CbcParam
    \brief Class for cbc-generic control parameters

    Adds parameter type codes and push/pull functions to the generic parameter
    object.
*/

class CBCLIB_EXPORT CbcParam : public CoinParam {

public:

   /*! \name Enumeration types used to index parameters */
   //@{
   
   /*! \brief Codes to specify psramters */

   enum CbcParamCode {
      FIRSTPARAM = 0,
      
      // Unused paramters that we may delete
      FIRSTUNUSEDPARAM,
      BRANCHSTRATEGY,
      CLEARCUTS,
      SOLVER,
      LASTUNUSEDPARAM,
      
      // Help and Information Parameters
      FIRSTHELPPARAM,
      GENERALQUERY,
      FULLGENERALQUERY,
      HELP,
      LASTHELPPARAM,
      
      // Action Parameters
      FIRSTACTIONPARAM,
      BAB,
      DEBUG,
      ENVIRONMENT,
      END,
      EXIT,
      QUIT,
      STOP,
      EXPORT,
      IMPORT,
      MIPLIB,
      OUTDUPROWS,
      PRINTSOL,
      PRINTVERSION,
      READMIPSTART,
      READPRIORITIES,
      READSOL,
      READMODEL,
      READMODEL_OLD, // For backward compatibility
      SHOWUNIMP,
      SOLVECONTINUOUS,
      STATISTICS,
      STDIN,
      STRENGTHEN,
      UNITTEST,
      WRITEGMPLSOL,
      WRITEGMPLSOL_OLD, // For backward compatibility
      WRITEMODEL,
      WRITEMODEL_OLD, // For backward compatibility
      WRITENEXTSOL,
      WRITESOL,
      WRITESOL_OLD, // For backward compatibility
      WRITESOLBINARY,
      WRITESOLBINARY_OLD, // For backward compatibility
      WRITESTATS,
      LASTACTIONPARAM,
      
      // File Parameters
      FIRSTFILEPARAM,
      CSVSTATSFILE,
      DEBUGFILE,
      EXPORTFILE,
      IMPORTFILE,
      GMPLSOLFILE,
      MIPSTARTFILE,
      MODELFILE,
      NEXTSOLFILE,
      PRIORITYFILE,
      SOLUTIONFILE,
      SOLUTIONBINARYFILE,
      LASTFILEPARAM,

      // Directory Paramaters
      FIRSTDIRECTORYPARAM,
      DIRECTORY,
      DIRSAMPLE,
      DIRNETLIB,
      DIRMIPLIB,
      LASTDIRECTORYPARAM,

      // String Parameters
      FIRSTSTRINGPARAM,
      PRINTMASK,
      LASTSTRINGPARAM,
      
      // Cut Parameters
      FIRSTCUTPARAM,
      CLIQUECUTS,
      CUTSTRATEGY,
      FLOWCUTS,
      GMICUTS,
      GOMORYCUTS,
      KNAPSACKCUTS,
      LAGOMORYCUTS,
      LANDPCUTS,
      LATWOMIRCUTS,
      MIRCUTS,
      //ODDHOLECUTS, // Not used
      ODDWHEELCUTS,
      PROBINGCUTS,
      REDSPLITCUTS,
      REDSPLIT2CUTS,
      RESIDCAPCUTS,
      TWOMIRCUTS,
      ZEROHALFCUTS,
      LASTCUTPARAM,
      
      // Heuristic Parameters
      FIRSTHEURPARAM,
      COMBINE,
      CROSSOVER,
      DINS,
      DIVINGC,
      DIVINGF,
      DIVINGG,
      DIVINGL,
      DIVINGP,
      DIVINGS,
      DIVINGV,
      DW,
      FPUMP,
      GREEDY,
      HEURISTICSTRATEGY,
      LOCALTREE,
      NAIVE,
      PIVOTANDFIX,
      PIVOTANDCOMPLEMENT,
      PROXIMITY,
      RANDROUND,
      RENS,
      RINS,
      ROUNDING,
      VND,
      LASTHEURPARAM,
      
      // On/Off Parameters
      FIRSTBOOLPARAM,
      CPX,
      DOHEURISTIC,
      ERRORSALLOWED,
      MESSAGES,
      PREPROCNAMES,
      SOS,
      USESOLUTION,
      LASTBOOLPARAM,
      
      // Keyword Parameters
      FIRSTKWDPARAM,
      COMMANDPRINTLEVEL,
      CLQSTRENGTHENING,
      BRANCHPRIORITY,
      CUTOFFCONSTRAINT,
      INTPRINT,
      NODESTRATEGY,
      ORBITAL,
      PREPROCESS,
      SOSPRIORITIZE,
      STRATEGY,
      TIMEMODE,
      USECGRAPH,
      LASTKWDPARAM,
      
      // Integer Parameters
      FIRSTINTPARAM,
      BKPIVOTINGSTRATEGY,
      BKMAXCALLS,
      BKCLQEXTMETHOD,
      CPP,
      CUTDEPTH,
      CUTLENGTH,
      CUTPASSINTREE,
      DEPTHMINIBAB,
      DIVEOPT,
      DIVEOPTSOLVES,
      DUMMY,
      EXPERIMENT,
      EXTRA1,
      EXTRA2,
      EXTRA3,
      EXTRA4,
      EXTRAVARIABLES,
      FPUMPITS,
      FPUMPTUNE,
      FPUMPTUNE2,
      HEUROPTIONS,
      LOGLEVEL,
      LPLOGLEVEL,
      MAXHOTITS,
      MAXSAVEDSOLS,
      MAXSLOWCUTS,
      MOREMOREMIPOPTIONS,
      MULTIPLEROOTS,
      ODDWEXTMETHOD,
      OUTPUTFORMAT,
      PRINTOPTIONS,
      PROCESSTUNE,
      RANDOMSEED,
      STRONGSTRATEGY,
      TESTOSI,
      THREADS,
      USERCBC,
      VERBOSE,
      VUBTRY,
      LASTINTPARAM,
      
      // Double Parameters
      FIRSTDBLPARAM,
      ARTIFICIALCOST,
      DEXTRA3,
      DEXTRA4,
      DEXTRA5,
      DJFIX,
      FAKECUTOFF,
      FAKEINCREMENT,
      SMALLBAB,
      TIGHTENFACTOR,
      LASTDBLPARAM,
      
      FIRSTMODELPARAM,
      ALLOWABLEGAP,
      COSTSTRATEGY,
      //CUTDEPTH,
      CUTOFF,
      CUTPASS,
      DIRECTION,
      GAPRATIO,
      INCREMENT,
      INFEASIBILITYWEIGHT,
      INTEGERTOLERANCE,
      //LOGLEVEL,
      MAXIMIZE,
      MAXNODES,
      MAXNODESNOTIMPROVING,   // Added
      MAXSECONDSNOTIMPROVING, // Added
      MAXSOLS,                // Added
      MINIMIZE,
      MIPOPTIONS,
      MOREMIPOPTIONS,
      NUMBERANALYZE,
      NUMBERBEFORE,
      REVERSE,
      STRONGBRANCHING,
      TIMELIMIT,
      LASTMODELPARAM,
      
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

  CbcParam();

  /*! \brief Constructor for a parameter with a double value

      The default value is 0.0. Be careful to clearly indicate that \p lower and
      \p upper are real (double) values to distinguish this constructor from the
      constructor for an integer parameter.
    */
  CbcParam(int code, std::string name, std::string help,
           double lower = -COIN_DBL_MAX, double upper = COIN_DBL_MAX,
           std::string longHelp = "",
           CoinDisplayPriority display = displayPriorityHigh);

  /*! \brief Constructor for a parameter with an integer value

      The default value is 0.
    */
  CbcParam(int code, std::string name, std::string help,
           int lower = -COIN_INT_MAX, int upper = COIN_INT_MAX,
           std::string longHelp = "",
           CoinDisplayPriority display = displayPriorityHigh);

  /*! \brief Constructor for a parameter with string values (or no value)
    Type is not optional to resolve ambiguity.

    The default value is "" for all such parameter types
  */
  CbcParam(int code, std::string name, CoinParam::CoinParamType type,
           std::string help, std::string longHelp = "",
           CoinDisplayPriority displayPriority = CoinParam::displayPriorityHigh);
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

  /*! \brief Get the associated Cbc parameters object */
  inline CbcParameters *parameters() const { return (parameters_); }

  /*! \brief Set the associated Cbc parameters object */
  inline void setParameters(CbcParameters *s) { parameters_ = s; }

  /*! \brief Get the associated model */
  inline CbcModel *model() const { return (model_); }

  /*! \brief Set the associated Cbc model object */
  inline void setModel(CbcModel *m) { model_ = m; }

//@}

private:
  /*! \name Data */
  //@{

  /// Parameter code
  int paramCode_;

  /// Pointer to enclosing parameters object
  CbcParameters *parameters_;

  /// Pointer to model object
  CbcModel * model_;

  //@}
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
