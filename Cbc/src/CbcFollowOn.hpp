// Edwin 11/10/2009-- carved out of CbcBranchActual
#ifndef CbcFollowOn_H
#define CbcFollowOn_H

#include "CbcBranchBase.hpp"
#include "CoinPackedMatrix.hpp"

/** Define a follow on class.
    The idea of this is that in air-crew scheduling problems crew may fly in on flight A
    and out on flight B or on some other flight.  A useful branch is one which on one side
    fixes all which go out on flight B to 0, while the other branch fixes all those that do NOT
    go out on flight B to 0.

    This branching rule should be in addition to normal rules and have a high priority.
*/

class CbcFollowOn : public CbcObject {

public:

    // Default Constructor
    CbcFollowOn ();

    /** Useful constructor
    */
    CbcFollowOn (CbcModel * model);

    // Copy constructor
    CbcFollowOn ( const CbcFollowOn &);

    /// Clone
    virtual CbcObject * clone() const;

    // Assignment operator
    CbcFollowOn & operator=( const CbcFollowOn& rhs);

    // Destructor
    ~CbcFollowOn ();

    /// Infeasibility - large is 0.5
    virtual double infeasibility(const OsiBranchingInformation * info,
                                 int &preferredWay) const;

    using CbcObject::feasibleRegion ;
    /// This looks at solution and sets bounds to contain solution
    virtual void feasibleRegion();

    /// Creates a branching object
    virtual CbcBranchingObject * createCbcBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) ;
    /// As some computation is needed in more than one place - returns row
    virtual int gutsOfFollowOn(int & otherRow, int & preferredWay) const;

protected:
    /// data
    /// Matrix
    CoinPackedMatrix matrix_;
    /// Matrix by row
    CoinPackedMatrix matrixByRow_;
    /// Possible rhs (if 0 then not possible)
    int * rhs_;
};

#endif