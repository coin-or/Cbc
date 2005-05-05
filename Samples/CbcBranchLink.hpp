// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcBranchLink_H
#define CbcBranchLink_H

#include "CbcBranchBase.hpp"

/** Define Special Linked Ordered Sets.

*/


class CbcLink : public CbcObject {

public:

  // Default Constructor 
  CbcLink ();

  /** Useful constructor - A valid solution is if all variables are zero
      apart from k*numberLink to (k+1)*numberLink-1 where k is 0 through
      numberInSet-1.  The length of weights array is numberInSet.
      For this simple version the variables in matrix are the numberInSet*numberLink
      starting at first. If weights null then 0,1,2..
  */
  CbcLink (CbcModel * model, int numberMembers,
           int numberLinks, int first,
           const double * weights, int setNumber);
  
  // Copy constructor 
  CbcLink ( const CbcLink &);
   
  /// Clone
  virtual CbcObject * clone() const;

  // Assignment operator 
  CbcLink & operator=( const CbcLink& rhs);

  // Destructor 
  ~CbcLink ();
  
  /// Infeasibility - large is 0.5
  virtual double infeasibility(int & preferredWay) const;

  /// This looks at solution and sets bounds to contain solution
  virtual void feasibleRegion();
  /// Creates a branching object
  virtual CbcBranchingObject * createBranch(int way) const;

  /// Number of members
  inline int numberMembers() const
  {return numberMembers_;};

  /// Number of links for each member
  inline int numberLinks() const
  {return numberLinks_;};

  /// First variable in matrix
  inline int first() const
  {return first_;};

  /** Array of weights */
  inline const double * weights() const
  { return weights_;};

private:
  /// data

  /// Weights
  double * weights_;

  /// Number of members
  int numberMembers_;
  /// Number of links
   int numberLinks_;
  /// First member
  int first_;
};
/** Branching object for Special ordered sets

    Variable_ is the set id number (redundant, as the object also holds a
    pointer to the set.
 */
class CbcLinkBranchingObject : public CbcBranchingObject {

public:

  // Default Constructor 
  CbcLinkBranchingObject ();

  // Useful constructor
  CbcLinkBranchingObject (CbcModel * model,  const CbcLink * set,
                          int way,
                          double separator);
  
  // Copy constructor 
  CbcLinkBranchingObject ( const CbcLinkBranchingObject &);
   
  // Assignment operator 
  CbcLinkBranchingObject & operator=( const CbcLinkBranchingObject& rhs);

  /// Clone
  virtual CbcBranchingObject * clone() const;

  // Destructor 
  virtual ~CbcLinkBranchingObject ();
  
  /// Does next branch and updates state
  virtual double branch(bool normalBranch=false);

  /** \brief Print something about branch - only if log level high
  */
  virtual void print(bool normalBranch);
private:
  /// data
  const CbcLink * set_;
  /// separator
  double separator_;
};
#endif
