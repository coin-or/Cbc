// Edwin 11/9/2009-- carved out of CbcBranchActual
#ifndef CbcClique_H
#define CbcClique_H

/// Define a clique class

class CbcClique : public CbcObject {

public:

    // Default Constructor
    CbcClique ();

    /** Useful constructor (which are integer indices)
        slack can denote a slack in set.
        If type == NULL then as if 1
    */
    CbcClique (CbcModel * model, int cliqueType, int numberMembers,
               const int * which, const char * type,
               int identifier, int slack = -1);

    // Copy constructor
    CbcClique ( const CbcClique &);

    /// Clone
    virtual CbcObject * clone() const;

    // Assignment operator
    CbcClique & operator=( const CbcClique& rhs);

    // Destructor
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

    /// Number of Non SOS members i.e. fixing to zero is strong
    inline int numberNonSOSMembers() const {
        return numberNonSOSMembers_;
    }

    /// Members (indices in range 0 ... numberIntegers_-1)
    inline const int * members() const {
        return members_;
    }

    /** Type of each member i.e. which way is strong 0=non SOS, 1 =SOS,
        index is 0 ... numberMembers_-1 */
    inline char type(int index) const {
        if (type_) return type_[index];
        else return 1;
    }

    /// Clique type - 0 <=, 1 ==
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

    /// Type of each member 0=SOS, 1 =clique
    char * type_;

    /// Clique type - 0 <=, 1 ==
    int cliqueType_;

    /// Which one is slack (if any) sequence within this set
    int slack_;
};
#endif