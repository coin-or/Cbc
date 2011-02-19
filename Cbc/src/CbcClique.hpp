// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/9/2009-- carved out of CbcBranchActual

#ifndef CbcClique_H
#define CbcClique_H

/** \brief Branching object for cliques

  A clique is defined to be a set of binary variables where fixing any one
  variable to its `strong' value fixes all other variables. An example is the
  most common SOS1 construction: a set of binary variables x_j s.t.  SUM{j}
  x_j = 1.  Setting any one variable to 1 forces all other variables to 0.
  (See comments for CbcSOS below.)

  Other configurations are possible, however: Consider x1-x2+x3 <= 0.
  Setting x1 (x3) to 1 forces x2 to 1 and x3 (x1) to 0. Setting x2 to 0
  forces x1 and x3 to 0.

  The proper point of view to take when interpreting CbcClique is
  `generalisation of SOS1 on binary variables.' To get into the proper frame
  of mind, here's an example.
  
  Consider the following sequence, where x_j = (1-y_j):
  \verbatim
     x1 + x2 + x3 <=  1		all strong at 1
     x1 - y2 + x3 <=  0		y2 strong at 0; x1, x3 strong at 1
    -y1 - y2 + x3 <= -1		y1, y2 strong at 0, x3 strong at 1
    -y1 - y2 - y3 <= -2		all strong at 0
  \endverbatim
  The first line is a standard SOS1 on binary variables.
  
  Variables with +1 coefficients are `SOS-style' and variables with -1
  coefficients are `non-SOS-style'. So #numberNonSOSMembers_ simply tells you
  how many variables have -1 coefficients. The implicit rhs for a clique is
  1-numberNonSOSMembers_.
*/
class CbcClique : public CbcObject {

public:

    /// Default Constructor
    CbcClique ();

    /** Useful constructor (which are integer indices) slack can denote a slack
	in set.  If type == NULL then as if 1
    */
    CbcClique (CbcModel * model, int cliqueType, int numberMembers,
               const int * which, const char * type,
               int identifier, int slack = -1);

    /// Copy constructor
    CbcClique ( const CbcClique &);

    /// Clone
    virtual CbcObject * clone() const;

    /// Assignment operator
    CbcClique & operator=( const CbcClique& rhs);

    /// Destructor
    virtual ~CbcClique ();

    /// Infeasibility - large is 0.5
    virtual double infeasibility(const OsiBranchingInformation * info,
                                 int &preferredWay) const;

    using CbcObject::feasibleRegion ;
    /// This looks at solution and sets bounds to contain solution
    virtual void feasibleRegion();

    /// Creates a branching object
    virtual CbcBranchingObject * createCbcBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) ;
    /// Number of members
    inline int numberMembers() const {
        return numberMembers_;
    }
    /** \brief Number of variables with -1 coefficient
      
      Number of non-SOS members, i.e., fixing to zero is strong.
      See comments at head of class, and comments for #type_.
    */
    inline int numberNonSOSMembers() const {
        return numberNonSOSMembers_;
    }

    /// Members (indices in range 0 ... numberIntegers_-1)
    inline const int * members() const {
        return members_;
    }

    /*! \brief Type of each member, i.e., which way is strong.
  
      This also specifies whether a variable has a +1 or -1 coefficient.
	- 0 => -1 coefficient, 0 is strong value
	- 1 => +1 coefficient, 1 is strong value
      If unspecified, all coefficients are assumed to be positive.
      
      Indexed as 0 .. numberMembers_-1
    */
    inline char type(int index) const {
        if (type_) return type_[index];
        else return 1;
    }

    /// Clique type: 0 is <=, 1 is ==
    inline int cliqueType() const {
        return cliqueType_;
    }
    /// Redoes data when sequence numbers change
    virtual void redoSequenceEtc(CbcModel * model, int numberColumns, const int * originalColumns);

protected:
    /// data
    /// Number of members
    int numberMembers_;

    /// Number of Non SOS members i.e. fixing to zero is strong
    int numberNonSOSMembers_;

    /// Members (indices in range 0 ... numberIntegers_-1)
    int * members_;

    /** \brief Strong value for each member.

      This also specifies whether a variable has a +1 or -1 coefficient.
        - 0 => -1 coefficient, 0 is strong value
        - 1 => +1 coefficient, 1 is strong value
      If unspecified, all coefficients are assumed to be positive.
    
      Indexed as 0 .. numberMembers_-1
    */
    char * type_;

    /** \brief Clique type

      0 defines a <= relation, 1 an equality. The assumed value of the rhs is
      numberNonSOSMembers_+1. (See comments for the class.)
    */
    int cliqueType_;

    /** \brief Slack variable for the clique
  
      Identifies the slack variable for the clique (typically added to convert
      a <= relation to an equality). Value is sequence number within clique
      menbers.
    */
    int slack_;
};

/** Branching object for unordered cliques

    Intended for cliques which are long enough to make it worthwhile
    but <= 64 members.  There will also be ones for long cliques.

    Variable_ is the clique id number (redundant, as the object also holds a
    pointer to the clique.
 */
class CbcCliqueBranchingObject : public CbcBranchingObject {

public:

    // Default Constructor
    CbcCliqueBranchingObject ();

    // Useful constructor
    CbcCliqueBranchingObject (CbcModel * model,  const CbcClique * clique,
                              int way,
                              int numberOnDownSide, const int * down,
                              int numberOnUpSide, const int * up);

    // Copy constructor
    CbcCliqueBranchingObject ( const CbcCliqueBranchingObject &);

    // Assignment operator
    CbcCliqueBranchingObject & operator=( const CbcCliqueBranchingObject& rhs);

    /// Clone
    virtual CbcBranchingObject * clone() const;

    // Destructor
    virtual ~CbcCliqueBranchingObject ();

    using CbcBranchingObject::branch ;
    /// Does next branch and updates state
    virtual double branch();

    using CbcBranchingObject::print ;
    /** \brief Print something about branch - only if log level high
    */
    virtual void print();

    /** Return the type (an integer identifier) of \c this */
    virtual CbcBranchObjType type() const {
        return CliqueBranchObj;
    }

    /** Compare the original object of \c this with the original object of \c
        brObj. Assumes that there is an ordering of the original objects.
        This method should be invoked only if \c this and brObj are of the same
        type.
        Return negative/0/positive depending on whether \c this is
        smaller/same/larger than the argument.
    */
    virtual int compareOriginalObject(const CbcBranchingObject* brObj) const;

    /** Compare the \c this with \c brObj. \c this and \c brObj must be of the
        same type and must have the same original object, but they may have
        different feasible regions.
        Return the appropriate CbcRangeCompare value (first argument being the
        sub/superset if that's the case). In case of overlap (and if \c
        replaceIfOverlap is true) replace the current branching object with one
        whose feasible region is the overlap.
     */
    virtual CbcRangeCompare compareBranchingObject
    (const CbcBranchingObject* brObj, const bool replaceIfOverlap = false);

private:
    /// data
    const CbcClique * clique_;
    /// downMask - bit set to fix to weak bounds, not set to leave unfixed
    unsigned int downMask_[2];
    /// upMask - bit set to fix to weak bounds, not set to leave unfixed
    unsigned int upMask_[2];
};

/** Unordered Clique Branching Object class.
    These are for cliques which are > 64 members
    Variable is number of clique.
 */
class CbcLongCliqueBranchingObject : public CbcBranchingObject {

public:

    // Default Constructor
    CbcLongCliqueBranchingObject ();

    // Useful constructor
    CbcLongCliqueBranchingObject (CbcModel * model,  const CbcClique * clique,
                                  int way,
                                  int numberOnDownSide, const int * down,
                                  int numberOnUpSide, const int * up);

    // Copy constructor
    CbcLongCliqueBranchingObject ( const CbcLongCliqueBranchingObject &);

    // Assignment operator
    CbcLongCliqueBranchingObject & operator=( const CbcLongCliqueBranchingObject& rhs);

    /// Clone
    virtual CbcBranchingObject * clone() const;

    // Destructor
    virtual ~CbcLongCliqueBranchingObject ();

    using CbcBranchingObject::branch ;
    /// Does next branch and updates state
    virtual double branch();

    using CbcBranchingObject::print ;
    /** \brief Print something about branch - only if log level high
    */
    virtual void print();

    /** Return the type (an integer identifier) of \c this */
    virtual CbcBranchObjType type() const {
        return LongCliqueBranchObj;
    }

    /** Compare the original object of \c this with the original object of \c
        brObj. Assumes that there is an ordering of the original objects.
        This method should be invoked only if \c this and brObj are of the same
        type.
        Return negative/0/positive depending on whether \c this is
        smaller/same/larger than the argument.
    */
    virtual int compareOriginalObject(const CbcBranchingObject* brObj) const;

    /** Compare the \c this with \c brObj. \c this and \c brObj must be os the
        same type and must have the same original object, but they may have
        different feasible regions.
        Return the appropriate CbcRangeCompare value (first argument being the
        sub/superset if that's the case). In case of overlap (and if \c
        replaceIfOverlap is true) replace the current branching object with one
        whose feasible region is the overlap.
     */
    virtual CbcRangeCompare compareBranchingObject
    (const CbcBranchingObject* brObj, const bool replaceIfOverlap = false);

private:
    /// data
    const CbcClique * clique_;
    /// downMask - bit set to fix to weak bounds, not set to leave unfixed
    unsigned int * downMask_;
    /// upMask - bit set to fix to weak bounds, not set to leave unfixed
    unsigned int * upMask_;
};

#endif

