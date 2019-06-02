// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/9/2009-- carved out of CbcBranchActual

#ifndef CbcSOS_H
#define CbcSOS_H

/** \brief Branching object for Special Ordered Sets of type 1 and 2.

  SOS1 are an ordered set of variables where at most one variable can be
  non-zero. SOS1 are commonly defined with binary variables (interpreted as
  selection between alternatives) but this is not necessary.  An SOS1 with
  all binary variables is a special case of a clique (setting any one
  variable to 1 forces all others to 0).

  In theory, the implementation makes no assumptions about integrality in
  Type 1 sets. In practice, there are places where the code seems to have been
  written with a binary SOS mindset. Current development of SOS branching
  objects is proceeding in OsiSOS.

  SOS2 are an ordered set of variables in which at most two consecutive
  variables can be non-zero and must sum to 1 (interpreted as interpolation
  between two discrete values). By definition the variables are non-integer.
*/

class CbcSOS : public CbcObject {

public:
  // Default Constructor
  CbcSOS();

  /** \brief Constructor with SOS type and member information

    Type specifies SOS 1 or 2. Identifier is an arbitrary value.

    Which should be an array of variable indices with numberMembers entries.
    Weights can be used to assign arbitrary weights to variables, in the order
    they are specified in which. If no weights are provided, a default array of
    0, 1, 2, ... is generated.
	*/

  CbcSOS(CbcModel *model, int numberMembers,
    const int *which, const double *weights, int identifier,
    int type = 1);

  // Copy constructor
  CbcSOS(const CbcSOS &);

  /// Clone
  virtual CbcObject *clone() const;

  // Assignment operator
  CbcSOS &operator=(const CbcSOS &rhs);

  // Destructor
  virtual ~CbcSOS();

  /// Infeasibility - large is 0.5
  virtual double infeasibility(const OsiBranchingInformation *info,
    int &preferredWay) const;

  using CbcObject::feasibleRegion;
  /// This looks at solution and sets bounds to contain solution
  virtual void feasibleRegion();

  /// Creates a branching object
  virtual CbcBranchingObject *createCbcBranch(OsiSolverInterface *solver, const OsiBranchingInformation *info, int way);

  /** Pass in information on branch just done and create CbcObjectUpdateData instance.
        If object does not need data then backward pointer will be NULL.
        Assumes can get information from solver */
  virtual CbcObjectUpdateData createUpdateInformation(const OsiSolverInterface *solver,
    const CbcNode *node,
    const CbcBranchingObject *branchingObject);
  /// Update object by CbcObjectUpdateData
  virtual void updateInformation(const CbcObjectUpdateData &data);
  using CbcObject::solverBranch;
  /** Create an OsiSolverBranch object

        This returns NULL if branch not represented by bound changes
    */
  virtual OsiSolverBranch *solverBranch() const;
  /// Redoes data when sequence numbers change
  virtual void redoSequenceEtc(CbcModel *model, int numberColumns, const int *originalColumns);

  /// Construct an OsiSOS object
  OsiSOS *osiObject(const OsiSolverInterface *solver) const;
  /// Number of members
  inline int numberMembers() const
  {
    return numberMembers_;
  }

  /// Members (indices in range 0 ... numberColumns-1)
  inline const int *members() const
  {
    return members_;
  }

  /// SOS type
  inline int sosType() const
  {
    return sosType_;
  }
  /// Down number times
  inline int numberTimesDown() const
  {
    return numberTimesDown_;
  }
  /// Up number times
  inline int numberTimesUp() const
  {
    return numberTimesUp_;
  }

  /** Array of weights */
  inline const double *weights() const
  {
    return weights_;
  }

  /// Set number of members
  inline void setNumberMembers(int n)
  {
    numberMembers_ = n;
  }

  /// Members (indices in range 0 ... numberColumns-1)
  inline int *mutableMembers() const
  {
    return members_;
  }

  /** Array of weights */
  inline double *mutableWeights() const
  {
    return weights_;
  }

  /** \brief Return true if object can take part in normal heuristics
    */
  virtual bool canDoHeuristics() const
  {
    return (sosType_ == 1 && integerValued_);
  }
  /// Set whether set is integer valued or not
  inline void setIntegerValued(bool yesNo)
  {
    integerValued_ = yesNo;
  }

protected:
  /// data

  /// Members (indices in range 0 ... numberColumns-1)
  int *members_;
  /** \brief Weights for individual members

    Arbitrary weights for members. Can be used to attach meaning to variable
    values independent of objective coefficients. For example, if the SOS set
    comprises binary variables used to choose a facility of a given size, the
    weight could be the corresponding facilty size. Fractional values of the
    SOS variables can then be used to estimate ideal facility size.

    Weights cannot be completely arbitrary. From the code, they must be
    differ by at least 1.0e-7
  */

  double *weights_;
  /// Current pseudo-shadow price estimate down
  mutable double shadowEstimateDown_;
  /// Current pseudo-shadow price estimate up
  mutable double shadowEstimateUp_;
  /// Down pseudo ratio
  double downDynamicPseudoRatio_;
  /// Up pseudo ratio
  double upDynamicPseudoRatio_;
  /// Number of times we have gone down
  int numberTimesDown_;
  /// Number of times we have gone up
  int numberTimesUp_;
  /// Number of members
  int numberMembers_;
  /// SOS type
  int sosType_;
  /// Whether integer valued
  bool integerValued_;
  /// Whether odd values e.g. negative
  bool oddValues_;
};

/** Branching object for Special ordered sets

    Variable_ is the set id number (redundant, as the object also holds a
    pointer to the set.
 */
class CbcSOSBranchingObject : public CbcBranchingObject {

public:
  // Default Constructor
  CbcSOSBranchingObject();

  // Useful constructor
  CbcSOSBranchingObject(CbcModel *model, const CbcSOS *clique,
    int way,
    double separator);

  // Copy constructor
  CbcSOSBranchingObject(const CbcSOSBranchingObject &);

  // Assignment operator
  CbcSOSBranchingObject &operator=(const CbcSOSBranchingObject &rhs);

  /// Clone
  virtual CbcBranchingObject *clone() const;

  // Destructor
  virtual ~CbcSOSBranchingObject();

  using CbcBranchingObject::branch;
  /// Does next branch and updates state
  virtual double branch();
  /** Update bounds in solver as in 'branch' and update given bounds.
        branchState is -1 for 'down' +1 for 'up' */
  virtual void fix(OsiSolverInterface *solver,
    double *lower, double *upper,
    int branchState) const;

  /** Reset every information so that the branching object appears to point to
        the previous child. This method does not need to modify anything in any
        solver. */
  virtual void previousBranch()
  {
    CbcBranchingObject::previousBranch();
    computeNonzeroRange();
  }

  using CbcBranchingObject::print;
  /** \brief Print something about branch - only if log level high
    */
  virtual void print();

  /** Return the type (an integer identifier) of \c this */
  virtual CbcBranchObjType type() const
  {
    return SoSBranchObj;
  }

  /** Compare the original object of \c this with the original object of \c
        brObj. Assumes that there is an ordering of the original objects.
        This method should be invoked only if \c this and brObj are of the same
        type.
        Return negative/0/positive depending on whether \c this is
        smaller/same/larger than the argument.
    */
  virtual int compareOriginalObject(const CbcBranchingObject *brObj) const;

  /** Compare the \c this with \c brObj. \c this and \c brObj must be os the
        same type and must have the same original object, but they may have
        different feasible regions.
        Return the appropriate CbcRangeCompare value (first argument being the
        sub/superset if that's the case). In case of overlap (and if \c
        replaceIfOverlap is true) replace the current branching object with one
        whose feasible region is the overlap.
     */
  virtual CbcRangeCompare compareBranchingObject(const CbcBranchingObject *brObj, const bool replaceIfOverlap = false);

  /** Fill out the \c firstNonzero_ and \c lastNonzero_ data members */
  void computeNonzeroRange();

protected:
  /// data
  const CbcSOS *set_;
  /// separator
  double separator_;
  /** The following two members describe the range in the members_ of the
        original object that whose upper bound is not fixed to 0. This is not
        necessary for Cbc to function correctly, this is there for heuristics so
        that separate branching decisions on the same object can be pooled into
        one branching object. */
  int firstNonzero_;
  int lastNonzero_;
};
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
