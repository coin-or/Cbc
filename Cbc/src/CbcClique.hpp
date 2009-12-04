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

#if 0
    // No need to override. Default works fine.
    /** Reset every information so that the branching object appears to point to
        the previous child. This method does not need to modify anything in any
        solver. */
    virtual void previousBranch();
#endif

    using CbcBranchingObject::print ;
    /** \brief Print something about branch - only if log level high
    */
    virtual void print();

    /** Return the type (an integer identifier) of \c this */
    virtual int type() const {
        return 102;
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

#if 0
    // No need to override. Default works fine.
    /** Reset every information so that the branching object appears to point to
        the previous child. This method does not need to modify anything in any
        solver. */
    virtual void previousBranch();
#endif

    using CbcBranchingObject::print ;
    /** \brief Print something about branch - only if log level high
    */
    virtual void print();

    /** Return the type (an integer identifier) of \c this */
    virtual int type() const {
        return 103;
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