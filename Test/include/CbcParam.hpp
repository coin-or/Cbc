// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef CbcParam_H
#define CbcParam_H

#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"

/*! \brief Parameter codes

  Parameter type ranges are allocated as follows
  <ul>
    <li>   1 -- 100	double parameters
    <li> 101 -- 200	integer parameters
    <li> 201 -- 300	string parameters
    <li> 301 -- 400	`actions'
  </ul>

  `Actions' do not necessarily invoke an immediate action; it's just that they
  don't fit neatly into the parameters array.
  
  This coding scheme is in flux.
  NODESTRATEGY, BRANCHSTRATEGY, ADDCUTSSTRATEGY,
  CLEARCUTS, OSLSTUFF, CBCSTUFF are not used at present (03.10.24).
*/

enum CbcParameterType

{ GENERALQUERY = -100,
  
  PRIMALTOLERANCE = 1, DUALTOLERANCE,CUTOFF,TIMELIMIT,

  INFEASIBILITYWEIGHT = 51, INTEGERTOLERANCE,INCREMENT,ALLOWABLEGAP,

  DJFIX = 81, GAPRATIO,TIGHTENFACTOR,

  LOGLEVEL = 151, SOLVERLOGLEVEL, MAXNODES,STRONGBRANCHING,
  
  DIRECTION = 201,ERRORSALLOWED,KEEPNAMES,SCALING,

  NODESTRATEGY = 251,BRANCHSTRATEGY,ADDCUTSSTRATEGY,
  GOMORYCUTS,PROBINGCUTS,KNAPSACKCUTS,ODDHOLECUTS,PRESOLVE,
  ROUNDING,SOLVER,CLIQUECUTS,COSTSTRATEGY,FLOWCUTS,MIXEDCUTS,
  TWOMIRCUTS,
  
  DIRECTORY = 301,IMPORT,EXPORT,RESTORE,SAVE,SOLVECONTINUOUS,BAB,
    MAXIMIZE,MINIMIZE,EXIT,STDIN,UNITTEST,MIPLIB,SOLUTION,CLEARCUTS,
    VERSION,

  OSLSTUFF = 401,CBCSTUFF,

  INVALID = 1000 } ;


/// Very simple class for setting parameters

class CbcParam

{ public:

  /**@name Constructor and destructor */
  //@{
  /// Constructors
  CbcParam (  );
  CbcParam (std::string name, std::string help,
	   double lower, double upper, CbcParameterType type);
  CbcParam (std::string name, std::string help,
	   int lower, int upper, CbcParameterType type);
  // Other strings will be added by insert
  CbcParam (std::string name, std::string help, std::string defaultValue,
	   CbcParameterType type);
  // Action
  CbcParam (std::string name, std::string help,
	   CbcParameterType type);
  /// Copy constructor. 
  CbcParam(const CbcParam &);
  /// Assignment operator. This copies the data
    CbcParam & operator = (const CbcParam & rhs);
  /// Destructor
  ~CbcParam (  );
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
  int setDoubleParameter(OsiSolverInterface * model, double value) const;
  /// Gets a double parameter
  double doubleParameter(OsiSolverInterface * model) const;
  /// Sets a int parameter (nonzero code if error)
  int setIntParameter(OsiSolverInterface * model, int value) const;
  /// Gets a int parameter
  int intParameter(OsiSolverInterface * model) const;
  /// Sets a double parameter (nonzero code if error)
  int setDoubleParameter(CbcModel &model, double value) const;
  /// Gets a double parameter
  double doubleParameter(CbcModel &model) const;
  /// Checks a double parameter (nonzero code if error)
  int checkDoubleParameter(double value) const;
  /// Sets a int parameter (nonzero code if error)
  int setIntParameter(CbcModel &model, int value) const;
  /// Gets a int parameter
  int intParameter(CbcModel &model) const;
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
  { currentKeyWord_ = value; };
  /// Returns 1 if matches minimum, 2 if matches less, 0 if not matched
  int matches (std::string input) const;
  /// type
  inline CbcParameterType type() const
  { return type_;};
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
  // Type see CbcParameterType
  CbcParameterType type_;
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
  std::vector<std::string> longHelp_;
  /// Action
  CbcParameterType action_;
  /// Current keyWord (if a keyword parameter)
  int currentKeyWord_;
  //@}
};

#endif	/* CbcParam_H */
