// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef OsiCbcSolverInterface_H
#define OsiCbcSolverInterface_H

#include <string>
#include <cfloat>
#include <map>
#include "CbcModel.hpp"
#include "CoinPackedMatrix.hpp"
#include "OsiSolverInterface.hpp"
#include "CbcStrategy.hpp"
#include "CoinWarmStartBasis.hpp"

class OsiRowCut;
class OsiClpSolverInterface;
#ifndef COIN_DBL_MAX
static const double OsiCbcInfinity = DBL_MAX;
#else
static const double OsiCbcInfinity = COIN_DBL_MAX;
#endif

//#############################################################################

/** Cbc Solver Interface
    
Instantiation of OsiCbcSolverInterface for the Model Algorithm.

*/

class OsiCbcSolverInterface :
  virtual public OsiSolverInterface {
  friend int OsiCbcSolverInterfaceUnitTest(const std::string & mpsDir, const std::string & netlibDir);
  
public:
  //---------------------------------------------------------------------------
  /**@name Solve methods */
  //@{
  /// Solve initial LP relaxation 
  virtual void initialSolve();
  
  /// Resolve an LP relaxation after problem modification
  virtual void resolve();
  
  /// Invoke solver's built-in enumeration algorithm
  virtual void branchAndBound();
  //@}
  
  //---------------------------------------------------------------------------
  /**@name Parameter set/get methods
     
  The set methods return true if the parameter was set to the given value,
  false otherwise. There can be various reasons for failure: the given
  parameter is not applicable for the solver (e.g., refactorization
  frequency for the cbc algorithm), the parameter is not yet implemented
  for the solver or simply the value of the parameter is out of the range
  the solver accepts. If a parameter setting call returns false check the
  details of your solver.
  
  The get methods return true if the given parameter is applicable for the
  solver and is implemented. In this case the value of the parameter is
  returned in the second argument. Otherwise they return false.
  */
  //@{
  // Set an integer parameter
  bool setIntParam(OsiIntParam key, int value);
  // Set an double parameter
  bool setDblParam(OsiDblParam key, double value);
  // Set a string parameter
  bool setStrParam(OsiStrParam key, const std::string & value);
  // Get an integer parameter
  bool getIntParam(OsiIntParam key, int& value) const;
  // Get an double parameter
  bool getDblParam(OsiDblParam key, double& value) const;
  // Get a string parameter
  bool getStrParam(OsiStrParam key, std::string& value) const;
  // Set a hint parameter - overrides OsiSolverInterface
  virtual bool setHintParam(OsiHintParam key, bool yesNo=true,
                            OsiHintStrength strength=OsiHintTry,
                            void * otherInformation=NULL);
  /// Get a hint parameter
    virtual bool getHintParam(OsiHintParam key, bool& yesNo,
			      OsiHintStrength& strength,
			      void *& otherInformation) const;

  using OsiSolverInterface::getHintParam ;
  /// Get a hint parameter
    virtual bool getHintParam(OsiHintParam key, bool& yesNo,
			      OsiHintStrength& strength) const;
  //@}
  
  //---------------------------------------------------------------------------
  ///@name Methods returning info on how the solution process terminated
  //@{
  /// Are there a numerical difficulties?
  virtual bool isAbandoned() const;
  /// Is optimality proven?
  virtual bool isProvenOptimal() const;
  /// Is primal infeasiblity proven?
  virtual bool isProvenPrimalInfeasible() const;
  /// Is dual infeasiblity proven?
  virtual bool isProvenDualInfeasible() const;
  /// Is the given primal objective limit reached?
  virtual bool isPrimalObjectiveLimitReached() const;
  /// Is the given dual objective limit reached?
  virtual bool isDualObjectiveLimitReached() const;
  /// Iteration limit reached?
  virtual bool isIterationLimitReached() const;
  //@}
  
  //---------------------------------------------------------------------------
  /**@name WarmStart related methods */
  //@{
  
  /*! \brief Get an empty warm start object
    
  This routine returns an empty CoinWarmStartBasis object. Its purpose is
  to provide a way to give a client a warm start basis object of the
  appropriate type, which can resized and modified as desired.
  */
  
  virtual CoinWarmStart *getEmptyWarmStart () const;
  
  /// Get warmstarting information
  virtual CoinWarmStart* getWarmStart() const;
  /** Set warmstarting information. Return true/false depending on whether
      the warmstart information was accepted or not. */
  virtual bool setWarmStart(const CoinWarmStart* warmstart);
  //@}
  
  //---------------------------------------------------------------------------
  /**@name Hotstart related methods (primarily used in strong branching). <br>
     The user can create a hotstart (a snapshot) of the optimization process
     then reoptimize over and over again always starting from there.<br>
     <strong>NOTE</strong>: between hotstarted optimizations only
     bound changes are allowed. */
  //@{
  /// Create a hotstart point of the optimization process
  virtual void markHotStart();
  /// Optimize starting from the hotstart
  virtual void solveFromHotStart();
  /// Delete the snapshot
  virtual void unmarkHotStart();
  //@}
  
  //---------------------------------------------------------------------------
  /**@name Problem information methods
     
  These methods call the solver's query routines to return
  information about the problem referred to by the current object.
  Querying a problem that has no data associated with it result in
  zeros for the number of rows and columns, and NULL pointers from
  the methods that return vectors.
  
  Const pointers returned from any data-query method are valid as
  long as the data is unchanged and the solver is not called.
  */
  //@{
  /**@name Methods related to querying the input data */
  //@{
  /// Get number of columns
  virtual int getNumCols() const;
  
  /// Get number of rows
  virtual int getNumRows() const;
  
  /// Get number of nonzero elements
  virtual int getNumElements() const ;
  
  /// Get pointer to array[getNumCols()] of column lower bounds
  virtual const double * getColLower() const;
  
  /// Get pointer to array[getNumCols()] of column upper bounds
  virtual const double * getColUpper() const;
  
  /** Get pointer to array[getNumRows()] of row constraint senses.
      <ul>
      <li>'L' <= constraint
      <li>'E' =  constraint
      <li>'G' >= constraint
      <li>'R' ranged constraint
      <li>'N' free constraint
      </ul>
  */
  virtual const char * getRowSense() const;
  
  /** Get pointer to array[getNumRows()] of rows right-hand sides
      <ul>
      <li> if rowsense()[i] == 'L' then rhs()[i] == rowupper()[i]
      <li> if rowsense()[i] == 'G' then rhs()[i] == rowlower()[i]
      <li> if rowsense()[i] == 'R' then rhs()[i] == rowupper()[i]
      <li> if rowsense()[i] == 'N' then rhs()[i] == 0.0
      </ul>
  */
  virtual const double * getRightHandSide() const ;
  
  /** Get pointer to array[getNumRows()] of row ranges.
      <ul>
      <li> if rowsense()[i] == 'R' then
      rowrange()[i] == rowupper()[i] - rowlower()[i]
      <li> if rowsense()[i] != 'R' then
      rowrange()[i] is undefined
      </ul>
  */
  virtual const double * getRowRange() const ;
  
  /// Get pointer to array[getNumRows()] of row lower bounds
  virtual const double * getRowLower() const ;
  
  /// Get pointer to array[getNumRows()] of row upper bounds
  virtual const double * getRowUpper() const ;
  
  /// Get pointer to array[getNumCols()] of objective function coefficients
  virtual const double * getObjCoefficients() const; 
  
  /// Get objective function sense (1 for min (default), -1 for max)
  virtual double getObjSense() const ;
  
  /// Return true if column is continuous
  virtual bool isContinuous(int colNumber) const;
  
  
  /// Get pointer to row-wise copy of matrix
  virtual const CoinPackedMatrix * getMatrixByRow() const;
  
  /// Get pointer to column-wise copy of matrix
  virtual const CoinPackedMatrix * getMatrixByCol() const;
  
  /// Get solver's value for infinity
  virtual double getInfinity() const;
  //@}
  
  /**@name Methods related to querying the solution */
  //@{
  /// Get pointer to array[getNumCols()] of primal solution vector
  virtual const double * getColSolution() const; 
  
  /// Get pointer to array[getNumRows()] of dual prices
  virtual const double * getRowPrice() const;
  
  /// Get a pointer to array[getNumCols()] of reduced costs
  virtual const double * getReducedCost() const; 
  
  /** Get pointer to array[getNumRows()] of row activity levels (constraint
      matrix times the solution vector */
  virtual const double * getRowActivity() const; 
  
  /// Get objective function value
  virtual double getObjValue() const;
  
  /** Get how many iterations it took to solve the problem (whatever
      "iteration" mean to the solver. */
  virtual int getIterationCount() const ;
  
  /** Get as many dual rays as the solver can provide. (In case of proven
      primal infeasibility there should be at least one.)
      
      <strong>NOTE for implementers of solver interfaces:</strong> <br>
      The double pointers in the vector should point to arrays of length
      getNumRows() and they should be allocated via new[]. <br>
      
      <strong>NOTE for users of solver interfaces:</strong> <br>
      It is the user's responsibility to free the double pointers in the
      vector using delete[].
  */
  virtual std::vector<double*> getDualRays(int maxNumRays) const;
  /** Get as many primal rays as the solver can provide. (In case of proven
      dual infeasibility there should be at least one.)
      
      <strong>NOTE for implementers of solver interfaces:</strong> <br>
      The double pointers in the vector should point to arrays of length
      getNumCols() and they should be allocated via new[]. <br>
      
      <strong>NOTE for users of solver interfaces:</strong> <br>
      It is the user's responsibility to free the double pointers in the
      vector using delete[].
  */
  virtual std::vector<double*> getPrimalRays(int maxNumRays) const;
  
  //@}

  /*! \name Methods for row and column names.

    Because OsiCbc is a pass-through class, it's necessary to override any
    virtual method in order to be sure we catch an override by the underlying
    solver. See the OsiSolverInterface class documentation for detailed
    descriptions.
  */
  //@{

    /*! \brief Generate a standard name of the form Rnnnnnnn or Cnnnnnnn */

    virtual std::string dfltRowColName(char rc,
				 int ndx, unsigned digits = 7) const ;

    /*! \brief Return the name of the objective function */

    virtual std::string getObjName (unsigned maxLen = std::string::npos) const ;

    /*! \brief Set the name of the objective function */

    virtual void setObjName (std::string name) ;

    /*! \brief Return the name of the row.  */

    virtual std::string getRowName(int rowIndex,
				   unsigned maxLen = std::string::npos) const ;

    /*! \brief Return a pointer to a vector of row names */

    virtual const OsiNameVec &getRowNames() ;

    /*! \brief Set a row name */

    virtual void setRowName(int ndx, std::string name) ;

    /*! \brief Set multiple row names */

    virtual void setRowNames(OsiNameVec &srcNames,
		     int srcStart, int len, int tgtStart) ;

    /*! \brief Delete len row names starting at index tgtStart */

    virtual void deleteRowNames(int tgtStart, int len) ;
  
    /*! \brief Return the name of the column */

    virtual std::string getColName(int colIndex,
				   unsigned maxLen = std::string::npos) const ;

    /*! \brief Return a pointer to a vector of column names */

    virtual const OsiNameVec &getColNames() ;

    /*! \brief Set a column name */

    virtual void setColName(int ndx, std::string name) ;

    /*! \brief Set multiple column names */

    virtual void setColNames(OsiNameVec &srcNames,
		     int srcStart, int len, int tgtStart) ;

    /*! \brief Delete len column names starting at index tgtStart */
    virtual void deleteColNames(int tgtStart, int len) ;

  //@}

  //@}
  
  //---------------------------------------------------------------------------
  
  /**@name Problem modifying methods */
  //@{
  //-------------------------------------------------------------------------
  /**@name Changing bounds on variables and constraints */
  //@{
  /** Set an objective function coefficient */
  virtual void setObjCoeff( int elementIndex, double elementValue );

  using OsiSolverInterface::setColLower ;
  /** Set a single column lower bound<br>
      Use -DBL_MAX for -infinity. */
  virtual void setColLower( int elementIndex, double elementValue );
  
  using OsiSolverInterface::setColUpper ;
  /** Set a single column upper bound<br>
      Use DBL_MAX for infinity. */
  virtual void setColUpper( int elementIndex, double elementValue );
  
  /** Set a single column lower and upper bound */
  virtual void setColBounds( int elementIndex,
                             double lower, double upper );
  
  /** Set the bounds on a number of columns simultaneously<br>
      The default implementation just invokes setColLower() and
      setColUpper() over and over again.
      @param indexFirst,indexLast pointers to the beginning and after the
      end of the array of the indices of the variables whose
      <em>either</em> bound changes
      @param boundList the new lower/upper bound pairs for the variables
  */
  virtual void setColSetBounds(const int* indexFirst,
                               const int* indexLast,
                               const double* boundList);
  
  /** Set a single row lower bound<br>
      Use -DBL_MAX for -infinity. */
  virtual void setRowLower( int elementIndex, double elementValue );
  
  /** Set a single row upper bound<br>
      Use DBL_MAX for infinity. */
  virtual void setRowUpper( int elementIndex, double elementValue ) ;
  
  /** Set a single row lower and upper bound */
  virtual void setRowBounds( int elementIndex,
                             double lower, double upper ) ;
  
  /** Set the type of a single row<br> */
  virtual void setRowType(int index, char sense, double rightHandSide,
                          double range);
  
  /** Set the bounds on a number of rows simultaneously<br>
      The default implementation just invokes setRowLower() and
      setRowUpper() over and over again.
      @param indexFirst,indexLast pointers to the beginning and after the
      end of the array of the indices of the constraints whose
      <em>either</em> bound changes
      @param boundList the new lower/upper bound pairs for the constraints
  */
  virtual void setRowSetBounds(const int* indexFirst,
                               const int* indexLast,
                               const double* boundList);
  
  /** Set the type of a number of rows simultaneously<br>
      The default implementation just invokes setRowType()
      over and over again.
      @param indexFirst,indexLast pointers to the beginning and after the
      end of the array of the indices of the constraints whose
      <em>any</em> characteristics changes
      @param senseList the new senses
      @param rhsList   the new right hand sides
      @param rangeList the new ranges
  */
  virtual void setRowSetTypes(const int* indexFirst,
                              const int* indexLast,
                              const char* senseList,
                              const double* rhsList,
                              const double* rangeList);
  //@}
  
  //-------------------------------------------------------------------------
  /**@name Integrality related changing methods */
  //@{
  /** Set the index-th variable to be a continuous variable */
  virtual void setContinuous(int index);
  /** Set the index-th variable to be an integer variable */
  virtual void setInteger(int index);
  /** Set the variables listed in indices (which is of length len) to be
      continuous variables */
  virtual void setContinuous(const int* indices, int len);
  /** Set the variables listed in indices (which is of length len) to be
      integer variables */
  virtual void setInteger(const int* indices, int len);
  //@}
  
  //-------------------------------------------------------------------------
  /// Set objective function sense (1 for min (default), -1 for max,)
  virtual void setObjSense(double s ); 
  
  /** Set the primal solution column values
      
  colsol[numcols()] is an array of values of the problem column
  variables. These values are copied to memory owned by the
  solver object or the solver.  They will be returned as the
  result of colsol() until changed by another call to
  setColsol() or by a call to any solver routine.  Whether the
  solver makes use of the solution in any way is
  solver-dependent. 
  */
  virtual void setColSolution(const double * colsol);
  
  /** Set dual solution vector
      
  rowprice[numrows()] is an array of values of the problem row
  dual variables. These values are copied to memory owned by the
  solver object or the solver.  They will be returned as the
  result of rowprice() until changed by another call to
  setRowprice() or by a call to any solver routine.  Whether the
  solver makes use of the solution in any way is
  solver-dependent. 
  */
  virtual void setRowPrice(const double * rowprice);
  
  //-------------------------------------------------------------------------
  /**@name Methods to expand a problem.<br>
     Note that if a column is added then by default it will correspond to a
     continuous variable. */
  //@{
  using OsiSolverInterface::addCol ;
  /** */
  virtual void addCol(const CoinPackedVectorBase& vec,
                      const double collb, const double colub,   
                      const double obj);
  /** Add a column (primal variable) to the problem. */
  virtual void addCol(int numberElements, const int * rows, const double * elements,
                      const double collb, const double colub,   
                      const double obj) ;

  using OsiSolverInterface::addCols ;
  /** */
  virtual void addCols(const int numcols,
                       const CoinPackedVectorBase * const * cols,
                       const double* collb, const double* colub,   
                       const double* obj);
  /** */
  virtual void deleteCols(const int num, const int * colIndices);
  
  using OsiSolverInterface::addRow ;
  /** */
  virtual void addRow(const CoinPackedVectorBase& vec,
                      const double rowlb, const double rowub);
  /** */
  virtual void addRow(const CoinPackedVectorBase& vec,
                      const char rowsen, const double rowrhs,   
                      const double rowrng);

  using OsiSolverInterface::addRows ;
  /** */
  virtual void addRows(const int numrows,
                       const CoinPackedVectorBase * const * rows,
                       const double* rowlb, const double* rowub);
  /** */
  virtual void addRows(const int numrows,
                       const CoinPackedVectorBase * const * rows,
                       const char* rowsen, const double* rowrhs,   
                       const double* rowrng);
  /** */
  virtual void deleteRows(const int num, const int * rowIndices);
  
  //-----------------------------------------------------------------------
  /** Apply a collection of row cuts which are all effective.
      applyCuts seems to do one at a time which seems inefficient.
  */
  virtual void applyRowCuts(int numberCuts, const OsiRowCut * cuts);
  /** Apply a collection of row cuts which are all effective.
      applyCuts seems to do one at a time which seems inefficient.
      This uses array of pointers
  */
  virtual void applyRowCuts(int numberCuts, const OsiRowCut ** cuts);
  //@}
  //@}
  
  //---------------------------------------------------------------------------
  
public:
  
  /**@name Methods to input a problem */
  //@{
  /** Load in an problem by copying the arguments (the constraints on the
      rows are given by lower and upper bounds). If a pointer is 0 then the
      following values are the default:
      <ul>
      <li> <code>colub</code>: all columns have upper bound infinity
      <li> <code>collb</code>: all columns have lower bound 0 
      <li> <code>rowub</code>: all rows have upper bound infinity
      <li> <code>rowlb</code>: all rows have lower bound -infinity
      <li> <code>obj</code>: all variables have 0 objective coefficient
      </ul>
  */
  virtual void loadProblem(const CoinPackedMatrix& matrix,
                           const double* collb, const double* colub,   
                           const double* obj,
                           const double* rowlb, const double* rowub);
  
  /** Load in an problem by assuming ownership of the arguments (the
      constraints on the rows are given by lower and upper bounds). For
      default values see the previous method.  <br>
      <strong>WARNING</strong>: The arguments passed to this method will be
      freed using the C++ <code>delete</code> and <code>delete[]</code>
      functions. 
  */
  virtual void assignProblem(CoinPackedMatrix*& matrix,
    			     double*& collb, double*& colub, double*& obj,
    			     double*& rowlb, double*& rowub);
  
  /** Load in an problem by copying the arguments (the constraints on the
      rows are given by sense/rhs/range triplets). If a pointer is 0 then the
      following values are the default:
      <ul>
      <li> <code>colub</code>: all columns have upper bound infinity
      <li> <code>collb</code>: all columns have lower bound 0 
      <li> <code>obj</code>: all variables have 0 objective coefficient
      <li> <code>rowsen</code>: all rows are >=
      <li> <code>rowrhs</code>: all right hand sides are 0
      <li> <code>rowrng</code>: 0 for the ranged rows
      </ul>
  */
  virtual void loadProblem(const CoinPackedMatrix& matrix,
    			   const double* collb, const double* colub,
    			   const double* obj,
    			   const char* rowsen, const double* rowrhs,   
    			   const double* rowrng);
  
  /** Load in an problem by assuming ownership of the arguments (the
      constraints on the rows are given by sense/rhs/range triplets). For
      default values see the previous method. <br>
      <strong>WARNING</strong>: The arguments passed to this method will be
      freed using the C++ <code>delete</code> and <code>delete[]</code>
      functions. 
  */
  virtual void assignProblem(CoinPackedMatrix*& matrix,
    			     double*& collb, double*& colub, double*& obj,
    			     char*& rowsen, double*& rowrhs,
    			     double*& rowrng);
  
  /** Just like the other loadProblem() methods except that the matrix is
      given in a standard column major ordered format (without gaps). */
  virtual void loadProblem(const int numcols, const int numrows,
                           const CoinBigIndex * start, const int* index,
                           const double* value,
                           const double* collb, const double* colub,   
                           const double* obj,
                           const double* rowlb, const double* rowub);
  
  /** Just like the other loadProblem() methods except that the matrix is
      given in a standard column major ordered format (without gaps). */
  virtual void loadProblem(const int numcols, const int numrows,
                           const CoinBigIndex * start, const int* index,
                           const double* value,
                           const double* collb, const double* colub,   
                           const double* obj,
                           const char* rowsen, const double* rowrhs,   
                           const double* rowrng);

  using OsiSolverInterface::readMps ;
  /** Read an mps file from the given filename (defaults to Osi reader) - returns
      number of errors (see OsiMpsReader class) */
  virtual int readMps(const char *filename,
                      const char *extension = "mps") ;
  
  /** Write the problem into an mps file of the given filename.
      If objSense is non zero then -1.0 forces the code to write a
      maximization objective and +1.0 to write a minimization one.
      If 0.0 then solver can do what it wants */
  virtual void writeMps(const char *filename,
                        const char *extension = "mps",
                        double objSense=0.0) const;
  /** Write the problem into an mps file of the given filename,
      names may be null.  formatType is
      0 - normal
      1 - extra accuracy 
      2 - IEEE hex (later)
      
      Returns non-zero on I/O error
  */
  virtual int writeMpsNative(const char *filename, 
                             const char ** rowNames, const char ** columnNames,
                             int formatType=0,int numberAcross=2,
                             double objSense=0.0) const ;
  //@}
  
  /**@name Message handling (extra for Cbc messages).
     Normally I presume you would want the same language.
     If not then you could use underlying model pointer */
  //@{
  /// Set language
  void newLanguage(CoinMessages::Language language);
  void setLanguage(CoinMessages::Language language)
  {newLanguage(language);}
  //@}
  //---------------------------------------------------------------------------
  
  /**@name Cbc specific public interfaces */
  //@{
  /// Get pointer to Cbc model
  inline CbcModel * getModelPtr() const 
  { return modelPtr_;}
  /// Get pointer to underlying solver
  inline OsiSolverInterface * getRealSolverPtr() const 
  { return modelPtr_->solver();}
  /// Set cutoff bound on the objective function.
  inline void setCutoff(double value) 
  { modelPtr_->setCutoff(value);}
  /// Get the cutoff bound on the objective function - always as minimize
  inline double getCutoff() const
  { return modelPtr_->getCutoff();}
  /// Set the CbcModel::CbcMaxNumNode maximum node limit 
  inline void setMaximumNodes( int value)
  { modelPtr_->setMaximumNodes(value);}
  /// Get the CbcModel::CbcMaxNumNode maximum node limit
  inline int getMaximumNodes() const
  { return modelPtr_->getMaximumNodes();}
  /// Set the CbcModel::CbcMaxNumSol maximum number of solutions
  inline void setMaximumSolutions( int value) 
  { modelPtr_->setMaximumSolutions(value);}
  /// Get the CbcModel::CbcMaxNumSol maximum number of solutions 
  inline int getMaximumSolutions() const 
  { return modelPtr_->getMaximumSolutions();}
  /// Set the CbcModel::CbcMaximumSeconds maximum number of seconds 
  inline void setMaximumSeconds( double value) 
  { modelPtr_->setMaximumSeconds(value);}
  /// Get the CbcModel::CbcMaximumSeconds maximum number of seconds 
  inline double getMaximumSeconds() const 
  { return modelPtr_->getMaximumSeconds();}
  /// Node limit reached?
  inline bool isNodeLimitReached() const
  { return modelPtr_->isNodeLimitReached();}
  /// Solution limit reached?
  inline bool isSolutionLimitReached() const
  { return modelPtr_->isSolutionLimitReached();}
  /// Get how many Nodes it took to solve the problem.
  inline int getNodeCount() const
  { return modelPtr_->getNodeCount();}
    /// Final status of problem - 0 finished, 1 stopped, 2 difficulties
    inline int status() const
  { return modelPtr_->status();}
  /** Pass in a message handler
  
    It is the client's responsibility to destroy a message handler installed
    by this routine; it will not be destroyed when the solver interface is
    destroyed. 
  */
  virtual void passInMessageHandler(CoinMessageHandler * handler);
  //@}
  
  //---------------------------------------------------------------------------
  
  /**@name Constructors and destructors */
  //@{
  /// Default Constructor
  OsiCbcSolverInterface (OsiSolverInterface * solver=NULL,
                         CbcStrategy * strategy=NULL);
  
  /// Clone
  virtual OsiSolverInterface * clone(bool copyData = true) const;
  
  /// Copy constructor 
  OsiCbcSolverInterface (const OsiCbcSolverInterface &);
#if 0    
  /// Borrow constructor - only delete one copy
  OsiCbcSolverInterface (CbcModel * rhs, bool reallyOwn=false);
  
  /// Releases so won't error
  void releaseCbc();
#endif    
  /// Assignment operator 
  OsiCbcSolverInterface & operator=(const OsiCbcSolverInterface& rhs);
  
  /// Destructor 
  virtual ~OsiCbcSolverInterface ();
  
  //@}
  //---------------------------------------------------------------------------
  
protected:
  ///@name Protected methods
  //@{
  /** Apply a row cut (append to constraint matrix). */
  virtual void applyRowCut(const OsiRowCut& rc);
  
  /** Apply a column cut (adjust one or more bounds). */
  virtual void applyColCut(const OsiColCut& cc);
  //@}
  /**@name Protected member data */
  //@{
  /// Cbc model represented by this class instance
  mutable CbcModel * modelPtr_;
  //@}
};
// So unit test can find out if NDEBUG set
bool OsiCbcHasNDEBUG();

//#############################################################################
/** A function that tests the methods in the OsiCbcSolverInterface class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. Also, if this method is compiled with
    optimization, the compilation takes 10-15 minutes and the machine pages
    (has 256M core memory!)... */
int
OsiCbcSolverInterfaceUnitTest(const std::string & mpsDir, const std::string & netlibDir);

#endif
