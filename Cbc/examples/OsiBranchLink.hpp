// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef OsiBranchLink_H
#define OsiBranchLink_H

#include "OsiBranchingObject.hpp"

/** Define Special Linked Ordered Sets.

*/
class CoinWarmStartBasis;

class OsiOldLink : public OsiSOS {

public:

  // Default Constructor 
  OsiOldLink ();

  /** Useful constructor - A valid solution is if all variables are zero
      apart from k*numberLink to (k+1)*numberLink-1 where k is 0 through
      numberInSet-1.  The length of weights array is numberInSet.
      For this constructor the variables in matrix are the numberInSet*numberLink
      starting at first. If weights null then 0,1,2..
  */
  OsiOldLink (const OsiSolverInterface * solver, int numberMembers,
           int numberLinks, int first,
           const double * weights, int setNumber);
  /** Useful constructor - A valid solution is if all variables are zero
      apart from k*numberLink to (k+1)*numberLink-1 where k is 0 through
      numberInSet-1.  The length of weights array is numberInSet.
      For this constructor the variables are given by list - grouped.
      If weights null then 0,1,2..
  */
  OsiOldLink (const OsiSolverInterface * solver, int numberMembers,
           int numberLinks, int typeSOS, const int * which,
           const double * weights, int setNumber);
  
  // Copy constructor 
  OsiOldLink ( const OsiOldLink &);
   
  /// Clone
  virtual OsiObject * clone() const;

  // Assignment operator 
  OsiOldLink & operator=( const OsiOldLink& rhs);

  // Destructor 
  virtual ~OsiOldLink ();
  
  /// Infeasibility - large is 0.5
  virtual double infeasibility(const OsiBranchingInformation * info,int & whichWay) const;

  /** Set bounds to fix the variable at the current (integer) value.

    Given an integer value, set the lower and upper bounds to fix the
    variable. Returns amount it had to move variable.
  */
  virtual double feasibleRegion(OsiSolverInterface * solver, const OsiBranchingInformation * info) const;

  /** Creates a branching object

    The preferred direction is set by \p way, 0 for down, 1 for up.
  */
  virtual OsiBranchingObject * createBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) const;

  /// Redoes data when sequence numbers change
  virtual void resetSequenceEtc(int numberColumns, const int * originalColumns);

  /// Number of links for each member
  inline int numberLinks() const
  {return numberLinks_;};

  /** \brief Return true if object can take part in normal heuristics
  */
  virtual bool canDoHeuristics() const 
  {return false;};
  /** \brief Return true if branch should only bound variables
  */
  virtual bool boundBranch() const 
  {return false;};

private:
  /// data

  /// Number of links
  int numberLinks_;
};
/** Branching object for Linked ordered sets

 */
class OsiOldLinkBranchingObject : public OsiSOSBranchingObject {

public:

  // Default Constructor 
  OsiOldLinkBranchingObject ();

  // Useful constructor
  OsiOldLinkBranchingObject (OsiSolverInterface * solver,  const OsiOldLink * originalObject,
                          int way,
                          double separator);
  
  // Copy constructor 
  OsiOldLinkBranchingObject ( const OsiOldLinkBranchingObject &);
   
  // Assignment operator 
  OsiOldLinkBranchingObject & operator=( const OsiOldLinkBranchingObject& rhs);

  /// Clone
  virtual OsiBranchingObject * clone() const;

  // Destructor 
  virtual ~OsiOldLinkBranchingObject ();
  
  /// Does next branch and updates state
  virtual double branch(OsiSolverInterface * solver);

  /** \brief Print something about branch - only if log level high
  */
  virtual void print(const OsiSolverInterface * solver=NULL);
private:
  /// data
};
/** Define data for one link
    
*/


class OsiOneLink {

public:

  // Default Constructor 
  OsiOneLink ();

  /** Useful constructor - 
      
  */
  OsiOneLink (const OsiSolverInterface * solver, int xRow, int xColumn, int xyRow,
	      const char * functionString);
  
  // Copy constructor 
  OsiOneLink ( const OsiOneLink &);
   
  // Assignment operator 
  OsiOneLink & operator=( const OsiOneLink& rhs);

  // Destructor 
  virtual ~OsiOneLink ();
  
  /// data

  /// Row which defines x (if -1 then no x)
  int xRow_;
  /// Column which defines x
  int xColumn_;
  /// Output row
  int xyRow;
  /// Function
  std::string function_;
};
/** Define Special Linked Ordered Sets. New style

    members and weights may be stored in SOS object

    This is for y and x*f(y) and z*g(y) etc

*/


class OsiLink : public OsiSOS {

public:

  // Default Constructor 
  OsiLink ();

  /** Useful constructor -
      
  */
  OsiLink (const OsiSolverInterface * solver, int yRow,
	   int yColumn, double meshSize);
  
  // Copy constructor 
  OsiLink ( const OsiLink &);
   
  /// Clone
  virtual OsiObject * clone() const;

  // Assignment operator 
  OsiLink & operator=( const OsiLink& rhs);

  // Destructor 
  virtual ~OsiLink ();
  
  /// Infeasibility - large is 0.5
  virtual double infeasibility(const OsiBranchingInformation * info,int & whichWay) const;

  /** Set bounds to fix the variable at the current (integer) value.

    Given an integer value, set the lower and upper bounds to fix the
    variable. Returns amount it had to move variable.
  */
  virtual double feasibleRegion(OsiSolverInterface * solver, const OsiBranchingInformation * info) const;

  /** Creates a branching object

    The preferred direction is set by \p way, 0 for down, 1 for up.
  */
  virtual OsiBranchingObject * createBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) const;

  /// Redoes data when sequence numbers change
  virtual void resetSequenceEtc(int numberColumns, const int * originalColumns);

  /// Number of links for each member
  inline int numberLinks() const
  {return numberLinks_;};

  /** \brief Return true if object can take part in normal heuristics
  */
  virtual bool canDoHeuristics() const 
  {return false;};
  /** \brief Return true if branch should only bound variables
  */
  virtual bool boundBranch() const 
  {return false;};

private:
  /// data
  /// Current increment for y points
  double meshSize_;
  /// Links
  OsiOneLink * data_;
  /// Number of links
  int numberLinks_;
  /// Row which defines y
  int yRow_;
  /// Column which defines y
  int yColumn_;
};
/** Branching object for Linked ordered sets

 */
class OsiLinkBranchingObject : public OsiTwoWayBranchingObject {

public:

  // Default Constructor 
  OsiLinkBranchingObject ();

  // Useful constructor
  OsiLinkBranchingObject (OsiSolverInterface * solver,  const OsiLink * originalObject,
                          int way,
                          double separator);
  
  // Copy constructor 
  OsiLinkBranchingObject ( const OsiLinkBranchingObject &);
   
  // Assignment operator 
  OsiLinkBranchingObject & operator=( const OsiLinkBranchingObject& rhs);

  /// Clone
  virtual OsiBranchingObject * clone() const;

  // Destructor 
  virtual ~OsiLinkBranchingObject ();
  
  /// Does next branch and updates state
  virtual double branch(OsiSolverInterface * solver);

  /** \brief Print something about branch - only if log level high
  */
  virtual void print(const OsiSolverInterface * solver=NULL);
private:
  /// data
};
/** Define BiLinear objects

    This models x*y where one or both are integer

*/


class OsiBiLinear : public OsiObject2 {

public:

  // Default Constructor 
  OsiBiLinear ();

  /** Useful constructor - 
      This Adds in rows and variables to construct valid Linked Ordered Set
      Adds extra constraints to match other x/y
      So note not const solver
  */
  OsiBiLinear (OsiSolverInterface * solver, int xColumn,
	       int yColumn, int xyRow, double coefficient,
	       double xMesh, double yMesh,
	       int numberExistingObjects=0,const OsiObject ** objects=NULL );
  
  // Copy constructor 
  OsiBiLinear ( const OsiBiLinear &);
   
  /// Clone
  virtual OsiObject * clone() const;

  // Assignment operator 
  OsiBiLinear & operator=( const OsiBiLinear& rhs);

  // Destructor 
  virtual ~OsiBiLinear ();
  
  /// Infeasibility - large is 0.5
  virtual double infeasibility(const OsiBranchingInformation * info,int & whichWay) const;

  /** Set bounds to fix the variable at the current (integer) value.

    Given an integer value, set the lower and upper bounds to fix the
    variable. Returns amount it had to move variable.
  */
  virtual double feasibleRegion(OsiSolverInterface * solver, const OsiBranchingInformation * info) const;

  /** Creates a branching object

    The preferred direction is set by \p way, 0 for down, 1 for up.
  */
  virtual OsiBranchingObject * createBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) const;

  /// Redoes data when sequence numbers change
  virtual void resetSequenceEtc(int numberColumns, const int * originalColumns);

  // This does NOT set mutable stuff
  virtual double checkInfeasibility(const OsiBranchingInformation * info) const;
  
  /** \brief Return true if object can take part in normal heuristics
  */
  virtual bool canDoHeuristics() const 
  {return false;};
  /** \brief Return true if branch should only bound variables
  */
  virtual bool boundBranch() const 
  { return false;};
  /// X column
  inline int xColumn() const
  { return xColumn_;};
  /// Y column
  inline int yColumn() const
  { return yColumn_;};
  /// X satisfied if less than this away from mesh
  inline double xSatisfied() const
  { return xSatisfied_;};
  inline void setXSatisfied(double value)
  { xSatisfied_=value;};
  /// Y satisfied if less than this away from mesh
  inline double ySatisfied() const
  { return ySatisfied_;};
  inline void setYSatisfied(double value)
  { ySatisfied_=value;};
  /// XY satisfied if two version differ by less than this
  inline double xySatisfied() const
  { return xySatisfied_;};
  inline void setXYSatisfied(double value)
  { xySatisfied_=value;};
  /// 0 branch on either, 1 branch on x, 2 branch on y
  inline int branchingStrategy() const
  { return branchingStrategy_;};
  inline void setBranchingStrategy(int value)
  { branchingStrategy_=value;};
  /// Does work of branching
  void newBounds(OsiSolverInterface * solver, int way, short xOrY, double separator) const;
  /// Updates coefficients
  void updateCoefficients(const double * lower, const double * upper,
			  CoinPackedMatrix * matrix, CoinWarmStartBasis * basis) const;


private:
  /// data
  
  /// Coefficient
  double coefficient_;
  /// x mesh
  double xMeshSize_;
  /// y mesh
  double yMeshSize_;
  /// x satisfied if less than this away from mesh
  double xSatisfied_;
  /// y satisfied if less than this away from mesh
  double ySatisfied_;
  /// xy satisfied if less than this away from true
  double xySatisfied_;
  /// value of x or y to branch about
  mutable double xyBranchValue_;
  /// x column
  int xColumn_;
  /// y column 
  int yColumn_;
  /// First lambda (of 4)
  int firstLambda_;
  /// 0 branch on either, 1 branch on x, 2 branch on y
  int branchingStrategy_;
  /// x row
  int xRow_;
  /// y row (-1 if x*x)
  int yRow_;
  /// Output row
  int xyRow_;
  /// Convexity row
  int convexity_;
  /// Which chosen -1 none, 0 x, 1 y
  mutable short chosen_;
};
/** Branching object for BiLinear objects

 */
class OsiBiLinearBranchingObject : public OsiTwoWayBranchingObject {

public:

  // Default Constructor 
  OsiBiLinearBranchingObject ();

  // Useful constructor
  OsiBiLinearBranchingObject (OsiSolverInterface * solver,  const OsiBiLinear * originalObject,
                          int way,
                          double separator, int chosen);
  
  // Copy constructor 
  OsiBiLinearBranchingObject ( const OsiBiLinearBranchingObject &);
   
  // Assignment operator 
  OsiBiLinearBranchingObject & operator=( const OsiBiLinearBranchingObject& rhs);

  /// Clone
  virtual OsiBranchingObject * clone() const;

  // Destructor 
  virtual ~OsiBiLinearBranchingObject ();
  
  /// Does next branch and updates state
  virtual double branch(OsiSolverInterface * solver);

  /** \brief Print something about branch - only if log level high
  */
  virtual void print(const OsiSolverInterface * solver=NULL);
  /** \brief Return true if branch should only bound variables
  */
  virtual bool boundBranch() const 
  { return false;};
private:
  /// data
  /// 1 means branch on x, 2 branch on y
  short chosen_;
};
/// Define a single integer class - but one where you kep branching until fixed even if satsified


class OsiSimpleFixedInteger : public OsiSimpleInteger {

public:

  /// Default Constructor 
  OsiSimpleFixedInteger ();

  /// Useful constructor - passed solver index
  OsiSimpleFixedInteger (const OsiSolverInterface * solver, int iColumn);
  
  /// Useful constructor - passed solver index and original bounds
  OsiSimpleFixedInteger (int iColumn, double lower, double upper);
  
  /// Useful constructor - passed simple integer
  OsiSimpleFixedInteger (const OsiSimpleInteger &);
  
  /// Copy constructor 
  OsiSimpleFixedInteger ( const OsiSimpleFixedInteger &);
   
  /// Clone
  virtual OsiObject * clone() const;

  /// Assignment operator 
  OsiSimpleFixedInteger & operator=( const OsiSimpleFixedInteger& rhs);

  /// Destructor 
  virtual ~OsiSimpleFixedInteger ();
  
  /// Infeasibility - large is 0.5
  virtual double infeasibility(const OsiBranchingInformation * info, int & whichWay) const;
  /** Creates a branching object

    The preferred direction is set by \p way, 0 for down, 1 for up.
  */
  virtual OsiBranchingObject * createBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) const;
protected:
  /// data
  
};
#endif
