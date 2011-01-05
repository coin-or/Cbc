/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/
/*
  This file is part of cbc-generic.
*/

#ifndef CbcGenParam_H
#define CbcGenParam_H

/* \file CbcGenParam.hpp
   \brief Declarations for parameters that control the cbc-generic main
	  program.
*/

/*
  $Id: CbcGenParam.hpp 1173 2009-06-04 09:44:10Z forrest $
*/

class CbcGenCtlBlk ;

/*! \class CbcGenParam
    \brief Class for cbc-generic control parameters

    Adds parameter type codes and push/pull functions to the generic parameter
    object.
*/

class CbcGenParam : public CoinParam {

public:

    /*! \name Subtypes */
//@{

    /*! \enum CbcGenParamCode
        \brief Enumeration for cbc-generic parameters

      These are parameters that control the operation of the cbc-generic main
      program by operating on a CbcGenCtlBlk object. CBCGEN_FIRSTPARAM and
      CBCGEN_LASTPARM are markers to allow convenient separation of parameter
      groups.
    */
    typedef enum { CBCGEN_FIRSTPARAM = 0,

                   GENERALQUERY, FULLGENERALQUERY, HELP,
                   BAB, CLEARCUTS, CLIQUECUTS, COMBINE, COSTSTRATEGY,
                   CPP, CUTDEPTH, CUTSTRATEGY,
                   DEBUG, DIRECTORY, DJFIX, DUMMY,
                   ERRORSALLOWED, EXIT, EXPORT,
                   FLOWCUTS, FPUMP, FPUMPITS, GOMORYCUTS, GREEDY, HEURISTICSTRATEGY,
                   IMPORT, INTPRINT, KNAPSACKCUTS, LOCALTREE, LOGLEVEL,
                   MESSAGES, MIPLIB, MIXEDCUTS, ODDHOLECUTS, OUTDUPROWS, OUTPUTFORMAT,
                   PREPROCESS, PRINTMASK, PRINTOPTIONS, PRINTVERSION, PRIORITYIN, PROBINGCUTS,
                   REDSPLITCUTS, ROUNDING,
                   SOLUTION, SOLVECONTINUOUS, SOLVER, SOS, STDIN, STRENGTHEN,
                   TIGHTENFACTOR, TWOMIRCUTS,
                   UNITTEST, USERCBC, USESOLUTION, VERBOSE, SHOWUNIMP,

                   CBCGEN_LASTPARAM

                 } CbcGenParamCode ;

//@}

    /*! \name Constructors and Destructors

      Be careful how you specify parameters for the constructors! There's great
      potential for confusion.
    */
//@{
    /*! \brief Default constructor */

    CbcGenParam() ;

    /*! \brief Constructor for a parameter with a double value

      The default value is 0.0. Be careful to clearly indicate that \p lower and
      \p upper are real (double) values to distinguish this constructor from the
      constructor for an integer parameter.
    */
    CbcGenParam(CbcGenParamCode code, std::string name, std::string help,
                double lower, double upper, double dflt = 0.0,
                bool display = true) ;

    /*! \brief Constructor for a parameter with an integer value

      The default value is 0.
    */
    CbcGenParam(CbcGenParamCode code, std::string name, std::string help,
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
    CbcGenParam(CbcGenParamCode code, std::string name, std::string help,
                std::string firstValue, int dflt, bool display = true) ;

    /*! \brief Constructor for a string parameter

      The default string value must be specified explicitly to distinguish
      a string constructor from an action parameter constructor.
    */

    CbcGenParam(CbcGenParamCode code, std::string name, std::string help,
                std::string dflt, bool display = true) ;

    /*! \brief Constructor for an action parameter */

    CbcGenParam(CbcGenParamCode code, std::string name, std::string help,
                bool display = true) ;

    /*! \brief Copy constructor */

    CbcGenParam(const CbcGenParam &orig) ;

    /*! \brief Clone */

    CbcGenParam *clone() ;

    /*! \brief Assignment */

    CbcGenParam &operator=(const CbcGenParam &rhs) ;

    /*! \brief  Destructor */

    ~CbcGenParam() ;

//@}

    /*! \name Methods to query and manipulate a parameter object */
//@{

    /*! \brief Get the parameter code  */

    inline CbcGenParamCode paramCode() const {
        return (paramCode_) ;
    }

    /*! \brief Set the parameter code */

    inline void setParamCode(CbcGenParamCode code) {
        paramCode_ = code ;
    }

    /*! \brief Get the underlying cbc-generic control object */

    inline CbcGenCtlBlk *obj() const {
        return (obj_) ;
    }

    /*! \brief Set the underlying cbc-generic control object */

    inline void setObj(CbcGenCtlBlk *obj) {
        obj_ = obj ;
    }

//@}


private:

    /*! \name Data */
//@{

    /// Parameter code
    CbcGenParamCode paramCode_ ;

    /// cbc-generic control object
    CbcGenCtlBlk *obj_ ;

//@}

} ;

/*
  Declare the utility functions.
*/

namespace CbcGenParamUtils {
void addCbcGenParams(int &numParams, CoinParamVec &paramVec,
                     CbcGenCtlBlk *ctlBlk) ;
void loadGenParamObj(const CoinParamVec paramVec, int first, int last,
                     CbcGenCtlBlk *ctlBlk) ;

void saveSolution(const OsiSolverInterface *osi, std::string fileName) ;
bool readSolution(std::string fileName,
                  int &numRows, int &numCols, double &objVal,
                  double **rowActivity, double **dualVars,
                  double **primalVars, double **reducedCosts) ;

int doBaCParam(CoinParam *param) ;
int doDebugParam(CoinParam *param) ;
int doExitParam(CoinParam *param) ;
int doHelpParam(CoinParam *param) ;
int doImportParam(CoinParam *param) ;
int doPrintMaskParam(CoinParam *param) ;
int doNothingParam(CoinParam *param) ;
int doSolutionParam(CoinParam *param) ;
int doUnimplementedParam(CoinParam *param) ;
int doVersionParam(CoinParam *param) ;

int pushCbcGenDblParam(CoinParam *param) ;
int pushCbcGenIntParam(CoinParam *param) ;
int pushCbcGenKwdParam(CoinParam *param) ;
int pushCbcGenStrParam(CoinParam *param) ;

int pushCbcGenCutParam(CoinParam *param) ;
}


#endif
