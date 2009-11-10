// Edwin 11/10/2009-- carved out of CbcBranchActual
#ifndef CbcGeneralDepth_H
#define CbcGeneralDepth_H

#include "CbcGeneral.hpp"

#ifdef COIN_HAS_CLP

/** Define a catch all class.
    This will create a list of subproblems using partial evaluation
*/
#include "ClpSimplex.hpp"
#include "ClpNode.hpp"


class CbcGeneralDepth : public CbcGeneral {

public:

    // Default Constructor
    CbcGeneralDepth ();

    /** Useful constructor
        Just needs to point to model.
        Initial version does evaluation to depth N
        This is stored in CbcModel but may be
        better here
    */
    CbcGeneralDepth (CbcModel * model, int maximumDepth);

    // Copy constructor
    CbcGeneralDepth ( const CbcGeneralDepth &);

    /// Clone
    virtual CbcObject * clone() const;

    // Assignment operator
    CbcGeneralDepth & operator=( const CbcGeneralDepth& rhs);

    // Destructor
    ~CbcGeneralDepth ();

    /// Infeasibility - large is 0.5
    virtual double infeasibility(const OsiBranchingInformation * info,
                                 int &preferredWay) const;

    using CbcObject::feasibleRegion ;
    /// This looks at solution and sets bounds to contain solution
    virtual void feasibleRegion();

    /// Creates a branching object
    virtual CbcBranchingObject * createCbcBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) ;
    /// Return maximum number of nodes
    inline int maximumNodes() const {
        return maximumNodes_;
    }
    /// Get maximum depth
    inline int maximumDepth() const {
        return maximumDepth_;
    }
    /// Set maximum depth
    inline void setMaximumDepth(int value) {
        maximumDepth_ = value;
    }
    /// Get which solution
    inline int whichSolution() const {
        return whichSolution_;
    }
    /// Get ClpNode info
    inline ClpNode * nodeInfo(int which) {
        return nodeInfo_->nodeInfo_[which];
    }

    /// Redoes data when sequence numbers change
    virtual void redoSequenceEtc(CbcModel * model, int numberColumns, const int * originalColumns);

protected:
    /// data
    /// Maximum depth
    int maximumDepth_;
    /// Maximum nodes
    int maximumNodes_;
    /// Which node has solution (or -1)
    mutable int whichSolution_;
    /// Number of valid nodes (including whichSolution_)
    mutable int numberNodes_;
    /// For solving nodes
    mutable ClpNodeStuff * nodeInfo_;
};
#endif //COIN_HAS_CLP
#endif