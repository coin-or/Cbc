// Edwin 11/9/2009-- carved out of CbcBranchActual
/** Define an n-way class for variables.
    Only valid value is one at UB others at LB
    Normally 0-1
*/
#ifndef CbcNWay_H
#define CbcNWay_H

class CbcNWay : public CbcObject {

public:

    // Default Constructor
    CbcNWay ();

    /** Useful constructor (which are matrix indices)
    */
    CbcNWay (CbcModel * model, int numberMembers,
             const int * which, int identifier);

    // Copy constructor
    CbcNWay ( const CbcNWay &);

    /// Clone
    virtual CbcObject * clone() const;

    /// Assignment operator
    CbcNWay & operator=( const CbcNWay& rhs);

    /// Destructor
    virtual ~CbcNWay ();

    /// Set up a consequence for a single member
    void setConsequence(int iColumn, const CbcConsequence & consequence);

    /// Applies a consequence for a single member
    void applyConsequence(int iSequence, int state) const;

    /// Infeasibility - large is 0.5 (and 0.5 will give this)
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

    /// Members (indices in range 0 ... numberColumns-1)
    inline const int * members() const {
        return members_;
    }
    /// Redoes data when sequence numbers change
    virtual void redoSequenceEtc(CbcModel * model, int numberColumns, const int * originalColumns);

protected:
    /// data
    /// Number of members
    int numberMembers_;

    /// Members (indices in range 0 ... numberColumns-1)
    int * members_;
    /// Consequences (normally NULL)
    CbcConsequence ** consequence_;
};
#endif