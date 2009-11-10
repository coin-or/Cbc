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
#endif