// Edwin 11/9/2009-- carved out of CbcBranchActual
#ifndef CbcSOS_H
#define CbcSOS_H

/** Define Special Ordered Sets of type 1 and 2.  These do not have to be
    integer - so do not appear in lists of integers.

    which_ points directly to columns of matrix
*/


class CbcSOS : public CbcObject {

public:

    // Default Constructor
    CbcSOS ();

    /** Useful constructor - which are indices
        and  weights are also given.  If null then 0,1,2..
        type is SOS type
    */
    CbcSOS (CbcModel * model, int numberMembers,
            const int * which, const double * weights, int identifier,
            int type = 1);

    // Copy constructor
    CbcSOS ( const CbcSOS &);

    /// Clone
    virtual CbcObject * clone() const;

    // Assignment operator
    CbcSOS & operator=( const CbcSOS& rhs);

    // Destructor
    virtual ~CbcSOS ();

    /// Infeasibility - large is 0.5
    virtual double infeasibility(const OsiBranchingInformation * info,
                                 int &preferredWay) const;

    using CbcObject::feasibleRegion ;
    /// This looks at solution and sets bounds to contain solution
    virtual void feasibleRegion();

    /// Creates a branching object
    virtual CbcBranchingObject * createCbcBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) ;



    /** Pass in information on branch just done and create CbcObjectUpdateData instance.
        If object does not need data then backward pointer will be NULL.
        Assumes can get information from solver */
    virtual CbcObjectUpdateData createUpdateInformation(const OsiSolverInterface * solver,
            const CbcNode * node,
            const CbcBranchingObject * branchingObject);
    /// Update object by CbcObjectUpdateData
    virtual void updateInformation(const CbcObjectUpdateData & data) ;
    using CbcObject::solverBranch ;
    /** Create an OsiSolverBranch object

        This returns NULL if branch not represented by bound changes
    */
    virtual OsiSolverBranch * solverBranch() const;
    /// Redoes data when sequence numbers change
    virtual void redoSequenceEtc(CbcModel * model, int numberColumns, const int * originalColumns);

    /// Construct an OsiSOS object
    OsiSOS * osiObject(const OsiSolverInterface * solver) const;
    /// Number of members
    inline int numberMembers() const {
        return numberMembers_;
    }

    /// Members (indices in range 0 ... numberColumns-1)
    inline const int * members() const {
        return members_;
    }

    /// SOS type
    inline int sosType() const {
        return sosType_;
    }
    /// Down number times
    inline int numberTimesDown() const {
        return numberTimesDown_;
    }
    /// Up number times
    inline int numberTimesUp() const {
        return numberTimesUp_;
    }

    /** Array of weights */
    inline const double * weights() const {
        return weights_;
    }

    /// Set number of members
    inline void setNumberMembers(int n) {
        numberMembers_ = n;
    }

    /// Members (indices in range 0 ... numberColumns-1)
    inline int * mutableMembers() const {
        return members_;
    }

    /** Array of weights */
    inline double * mutableWeights() const {
        return weights_;
    }

    /** \brief Return true if object can take part in normal heuristics
    */
    virtual bool canDoHeuristics() const {
        return (sosType_ == 1 && integerValued_);
    }
    /// Set whether set is integer valued or not
    inline void setIntegerValued(bool yesNo) {
        integerValued_ = yesNo;
    }
private:
    /// data

    /// Members (indices in range 0 ... numberColumns-1)
    int * members_;
    /// Weights
    double * weights_;
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
};

/** Branching object for Special ordered sets

    Variable_ is the set id number (redundant, as the object also holds a
    pointer to the set.
 */
class CbcSOSBranchingObject : public CbcBranchingObject {

public:

    // Default Constructor
    CbcSOSBranchingObject ();

    // Useful constructor
    CbcSOSBranchingObject (CbcModel * model,  const CbcSOS * clique,
                           int way,
                           double separator);

    // Copy constructor
    CbcSOSBranchingObject ( const CbcSOSBranchingObject &);

    // Assignment operator
    CbcSOSBranchingObject & operator=( const CbcSOSBranchingObject& rhs);

    /// Clone
    virtual CbcBranchingObject * clone() const;

    // Destructor
    virtual ~CbcSOSBranchingObject ();

    using CbcBranchingObject::branch ;
    /// Does next branch and updates state
    virtual double branch();
    /** Update bounds in solver as in 'branch' and update given bounds.
        branchState is -1 for 'down' +1 for 'up' */
    virtual void fix(OsiSolverInterface * solver,
                     double * lower, double * upper,
                     int branchState) const ;

    /** Reset every information so that the branching object appears to point to
        the previous child. This method does not need to modify anything in any
        solver. */
    virtual void previousBranch() {
        CbcBranchingObject::previousBranch();
        computeNonzeroRange();
    }

    using CbcBranchingObject::print ;
    /** \brief Print something about branch - only if log level high
    */
    virtual void print();

    /** Return the type (an integer identifier) of \c this */
    virtual int type() const {
        return 104;
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

    /** Fill out the \c firstNonzero_ and \c lastNonzero_ data members */
    void computeNonzeroRange();

private:
    /// data
    const CbcSOS * set_;
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