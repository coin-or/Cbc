// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef CbcOrClpParam_H
#define CbcOrClpParam_H
/**
   This has parameter handling stuff which can be shared between Cbc and Clp (and Dylp etc).

   This (and .cpp) should be copied so that it is the same in Cbc/Test and Clp/Test.
   I know this is not elegant but it seems simplest.
   
   It uses COIN_USE_CBC for parameters wanted by CBC
   It uses COIN_USE_CLP for parameters wanted by CLP (or CBC using CLP)
   It could use COIN_USE_DYLP for parameters wanted by DYLP
   It could use COIN_USE_DYLP_OR_CLP for parameters wanted by DYLP or CLP etc etc

 */
class OsiSolverInterface;
class CbcModel;
class ClpSimplex;
/*! \brief Parameter codes

  Parameter type ranges are allocated as follows
  <ul>
    <li>   1 -- 100	double parameters
    <li> 101 -- 200	integer parameters
    <li> 201 -- 250	string parameters
    <li> 251 -- 300	cuts etc(string but broken out for clarity)
    <li> 301 -- 400	`actions'
  </ul>

  `Actions' do not necessarily invoke an immediate action; it's just that they
  don't fit neatly into the parameters array.
  
  This coding scheme is in flux.
  NODESTRATEGY, BRANCHSTRATEGY, ADDCUTSSTRATEGY,
  CLEARCUTS, OSLSTUFF, CBCSTUFF are not used at present (03.10.24).
*/

enum CbcOrClpParameterType

  { GENERALQUERY = -100,FULLGENERALQUERY,
    
    PRIMALTOLERANCE = 1, DUALTOLERANCE,CUTOFF,TIMELIMIT,
    DUALBOUND, PRIMALWEIGHT,  OBJSCALE, RHSSCALE,
    
    INFEASIBILITYWEIGHT = 51, INTEGERTOLERANCE,INCREMENT,ALLOWABLEGAP,
    
    DJFIX = 81, GAPRATIO,TIGHTENFACTOR,

    SOLVERLOGLEVEL=101,
#ifndef COIN_USE_CBC 
    LOGLEVEL = 101, 
#endif
    MAXFACTOR,PERTVALUE,MAXITERATION,PRESOLVEPASS,IDIOT,SPRINT,
    OUTPUTFORMAT,SLPVALUE,PRESOLVEOPTIONS,PRINTOPTIONS,SPECIALOPTIONS,

    STRONGBRANCHING=151,CUTDEPTH, MAXNODES,
#ifdef COIN_USE_CBC 
    LOGLEVEL , 
#endif
  
    DIRECTION=201,DUALPIVOT,SCALING,ERRORSALLOWED,KEEPNAMES,SPARSEFACTOR,
    PRIMALPIVOT,PRESOLVE,CRASH,BIASLU,PERTURBATION,MESSAGES,AUTOSCALE,
    CHOLESKY,KKT,BARRIERSCALE,GAMMA,CROSSOVER,PFI,
    
    NODESTRATEGY = 251,BRANCHSTRATEGY,ADDCUTSSTRATEGY,
    GOMORYCUTS,PROBINGCUTS,KNAPSACKCUTS,ODDHOLECUTS,
    ROUNDING,SOLVER,CLIQUECUTS,COSTSTRATEGY,FLOWCUTS,MIXEDCUTS,
    TWOMIRCUTS,
    
    DIRECTORY=301,IMPORT,EXPORT,RESTORE,SAVE,DUALSIMPLEX,PRIMALSIMPLEX,
    MAXIMIZE,MINIMIZE,EXIT,STDIN,UNITTEST,NETLIB_DUAL,NETLIB_PRIMAL,SOLUTION,
    TIGHTEN,FAKEBOUND,HELP,PLUSMINUS,NETWORK,ALLSLACK,REVERSE,BARRIER,NETLIB_BARRIER,
    REALLY_SCALE,BASISIN,BASISOUT,SOLVECONTINUOUS,CLEARCUTS,VERSION,

    BAB=351,MIPLIB,

    OSLSTUFF = 401,CBCSTUFF,

    INVALID = 1000 } ;

#include <vector>
#include <string>

/// Very simple class for setting parameters

class CbcOrClpParam
{
public:

  /**@name Constructor and destructor */
  //@{
  /// Constructors
  CbcOrClpParam (  );
  CbcOrClpParam (std::string name, std::string help,
	   double lower, double upper, CbcOrClpParameterType type,bool display=true);
  CbcOrClpParam (std::string name, std::string help,
	   int lower, int upper, CbcOrClpParameterType type,bool display=true);
  // Other strings will be added by insert
  CbcOrClpParam (std::string name, std::string help, std::string firstValue,
	   CbcOrClpParameterType type,int defaultIndex=0,bool display=true);
  // Action
  CbcOrClpParam (std::string name, std::string help,
	   CbcOrClpParameterType type,int indexNumber=-1,bool display=true);
  /// Copy constructor. 
  CbcOrClpParam(const CbcOrClpParam &);
  /// Assignment operator. This copies the data
    CbcOrClpParam & operator=(const CbcOrClpParam & rhs);
  /// Destructor
  ~CbcOrClpParam (  );
  //@}

  /**@name stuff */
  //@{
  /// Insert string (only valid for keywords)
  void append(std::string keyWord);
  /// Adds one help line
  void addHelp(std::string keyWord);
  /// Returns name
  inline std::string  name(  ) const {
    return name_;
  };
  /// Returns short help
  inline std::string  shortHelp(  ) const {
    return shortHelp_;
  };
  /// Sets a double parameter (nonzero code if error)
  int setDoubleParameter(CbcModel & model, double value) ;
  /// Gets a double parameter
  double doubleParameter(CbcModel & model) const;
  /// Sets a int parameter (nonzero code if error)
  int setIntParameter(CbcModel & model, int value) ;
  /// Gets a int parameter
  int intParameter(CbcModel & model) const;
  /// Sets a double parameter (nonzero code if error)
  int setDoubleParameter(ClpSimplex * model, double value) ;
  /// Gets a double parameter
  double doubleParameter(ClpSimplex * model) const;
  /// Sets a int parameter (nonzero code if error)
  int setIntParameter(ClpSimplex * model, int value) ;
  /// Gets a int parameter
  int intParameter(ClpSimplex * model) const;
  /// Sets a double parameter (nonzero code if error)
  int setDoubleParameter(OsiSolverInterface * model, double value) ;
  /// Gets a double parameter
  double doubleParameter(OsiSolverInterface * model) const;
  /// Sets a int parameter (nonzero code if error)
  int setIntParameter(OsiSolverInterface * model, int value) ;
  /// Gets a int parameter
  int intParameter(OsiSolverInterface * model) const;
  /// Checks a double parameter (nonzero code if error)
  int checkDoubleParameter(double value) const;
  /// Returns name which could match
  std::string matchName (  ) const;
  /// Returns parameter option which matches (-1 if none)
  int parameterOption ( std::string check ) const;
  /// Prints parameter options
  void printOptions (  ) const;
  /// Returns current parameter option
  inline std::string currentOption (  ) const
  { return definedKeyWords_[currentKeyWord_]; };
  /// Sets current parameter option
  inline void setCurrentOption ( int value )
  { currentKeyWord_=value; };
  /// Sets int value
  void setIntValue ( int value );
    inline int intValue () const
  { return intValue_; };
  /// Sets double value
  void setDoubleValue ( double value );
  inline double doubleValue () const
  { return doubleValue_; };
  /// Sets string value
  void setStringValue ( std::string value );
  inline std::string stringValue () const
  { return stringValue_; };
  /// Returns 1 if matches minimum, 2 if matches less, 0 if not matched
  int matches (std::string input) const;
  /// type
  inline CbcOrClpParameterType type() const
  { return type_;};
  /// whether to display
  inline bool displayThis() const
  { return display_;};
  /// Set Long help
  inline void setLonghelp(const std::string help) 
  {longHelp_=help;};
  /// Print Long help
  void printLongHelp() const;
  /// Print action and string
  void printString() const;
  /// type for classification
  inline int indexNumber() const
  { return indexNumber_;};
private:
  /// gutsOfConstructor
  void gutsOfConstructor();
  //@}
////////////////// data //////////////////
private:

  /**@name data 
   We might as well throw all type data in - could derive?
  */
  //@{
  // Type see CbcOrClpParameterType
  CbcOrClpParameterType type_;
  /// If double == okay
  double lowerDoubleValue_;
  double upperDoubleValue_;
  /// If int == okay
  int lowerIntValue_;
  int upperIntValue_;
  // Length of name
  unsigned int lengthName_;
  // Minimum match
  unsigned int lengthMatch_;
  /// set of valid strings
  std::vector<std::string> definedKeyWords_;
  /// Name
  std::string name_;
  /// Short help
  std::string shortHelp_;
  /// Long help
  std::string longHelp_;
  /// Action
  CbcOrClpParameterType action_;
  /// Current keyWord (if a keyword parameter)
  int currentKeyWord_;
  /// Display on ?
  bool display_;
  /// Integer parameter - current value
  int intValue_;
  /// Double parameter - current value
  double doubleValue_;
  /// String parameter - current value
  std::string stringValue_;
  /// index number to use for display purposes
  int indexNumber_;
  //@}
};
/// Simple read stuff
std::string CoinReadNextField();

std::string CoinReadGetCommand(int argc, const char *argv[]);
std::string CoinReadGetString(int argc, const char *argv[]);
// valid 0 - okay, 1 bad, 2 not there
int CoinReadGetIntField(int argc, const char *argv[],int * valid);
double CoinReadGetDoubleField(int argc, const char *argv[],int * valid);
void CoinReadPrintit(const char * input);

#define CBCMAXPARAMETERS 200
/*
  Subroutine to establish the cbc parameter array. See the description of
  class CbcOrClpParam for details. Pulled from C..Main() for clarity. 
*/
void establishParams (int &numberParameters, CbcOrClpParam *const parameters);
// Given a parameter type - returns its number in list
int whichParam (CbcOrClpParameterType name, 
		int numberParameters, CbcOrClpParam *const parameters);
#endif	/* CbcOrClpParam_H */
