// Edwin 11/9/2009-- carved out of CbcBranchActual

#ifndef CbcSimpleInteger_H
#define CbcSimpleInteger_H

#include "CbcIntegerBranchingObject.hpp"
/// Define a single integer class


class CbcSimpleInteger : public CbcObject {

public:

    // Default Constructor
    CbcSimpleInteger ();

    // Useful constructor - passed model and index
    CbcSimpleInteger (CbcModel * model,  int iColumn, double breakEven = 0.5);

    // Useful constructor - passed model and Osi object
    CbcSimpleInteger (CbcModel * model,  const OsiSimpleInteger * object);

    // Copy constructor
    CbcSimpleInteger ( const CbcSimpleInteger &);

    /// Clone
    virtual CbcObject * clone() const;

    // Assignment operator
    CbcSimpleInteger & operator=( const CbcSimpleInteger& rhs);

    // Destructor
    virtual ~CbcSimpleInteger ();
    /// Construct an OsiSimpleInteger object
    OsiSimpleInteger * osiObject() const;
    /// Infeasibility - large is 0.5
    virtual double infeasibility(const OsiBranchingInformation * info,
                                 int &preferredWay) const;

    using CbcObject::feasibleRegion ;
    /** Set bounds to fix the variable at the current (integer) value.

      Given an integer value, set the lower and upper bounds to fix the
      variable. Returns amount it had to move variable.
    */
    virtual double feasibleRegion(OsiSolverInterface * solver, const OsiBranchingInformation * info) const;

    /** Create a branching object and indicate which way to branch first.

        The branching object has to know how to create branches (fix
        variables, etc.)
    */
    virtual CbcBranchingObject * createCbcBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) ;
    /// Fills in a created branching object
    void fillCreateBranch(CbcIntegerBranchingObject * branching, const OsiBranchingInformation * info, int way) ;

    using CbcObject::solverBranch ;
    /** Create an OsiSolverBranch object

        This returns NULL if branch not represented by bound changes
    */
    virtual OsiSolverBranch * solverBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info) const;

    /** Set bounds to fix the variable at the current (integer) value.

      Given an integer value, set the lower and upper bounds to fix the
      variable. The algorithm takes a bit of care in order to compensate for
      minor numerical inaccuracy.
    */
    virtual void feasibleRegion();

    /** Column number if single column object -1 otherwise,
        so returns >= 0
        Used by heuristics
    */
    virtual int columnNumber() const;
    /// Set column number
    inline void setColumnNumber(int value) {
        columnNumber_ = value;
    }

    /** Reset variable bounds to their original values.

      Bounds may be tightened, so it may be good to be able to set this info in object.
     */
    virtual void resetBounds(const OsiSolverInterface * solver) ;
    /**  Change column numbers after preprocessing
     */
    virtual void resetSequenceEtc(int numberColumns, const int * originalColumns) ;
    /// Original bounds
    inline double originalLowerBound() const {
        return originalLower_;
    }
    inline void setOriginalLowerBound(double value) {
        originalLower_ = value;
    }
    inline double originalUpperBound() const {
        return originalUpper_;
    }
    inline void setOriginalUpperBound(double value) {
        originalUpper_ = value;
    }
    /// Breakeven e.g 0.7 -> >= 0.7 go up first
    inline double breakEven() const {
        return breakEven_;
    }
    /// Set breakeven e.g 0.7 -> >= 0.7 go up first
    inline void setBreakEven(double value) {
        breakEven_ = value;
    }


protected:
    /// data

    /// Original lower bound
    double originalLower_;
    /// Original upper bound
    double originalUpper_;
    /// Breakeven i.e. >= this preferred is up
    double breakEven_;
    /// Column number in model
    int columnNumber_;
    /// If -1 down always chosen first, +1 up always, 0 normal
    int preferredWay_;
};
#endif