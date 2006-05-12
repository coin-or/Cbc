// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcCountRowCut_H
#define CbcCountRowCut_H


class OsiCuts;
class OsiRowCut;
class CbcNodeInfo;

//#############################################################################
/** \brief OsiRowCut augmented with bookkeeping

  CbcCountRowCut is an OsiRowCut object augmented with bookkeeping
  information: a reference count and information that specifies the
  the generator that created the cut and the node to which it's associated.

  The general principles for handling the reference count are as follows:
  <ul>
    <li> Once it's determined how the node will branch, increment the
	 reference count under the assumption that all children will use
         all cuts currently tight at the node and will survive to be placed
	 in the search tree.
    <li> As this assumption is proven incorrect (a cut becomes loose, or a
	 child is fathomed), decrement the reference count accordingly.
  </ul>
  When all possible uses of a cut have been demonstrated to be unnecessary,
  the reference count (#numberPointingToThis_) will fall to zero. The
  CbcCountRowCut object (and its included OsiRowCut object) are then deleted.
*/

class CbcCountRowCut : public OsiRowCut {

public:
  
  /** @name Constructors & destructors */
  //@{

    /// Default Constructor 
    CbcCountRowCut ();

    /// `Copy' constructor using an OsiRowCut
    CbcCountRowCut ( const OsiRowCut &);  

    /// `Copy' constructor using an OsiRowCut and an CbcNodeInfo
    CbcCountRowCut(const OsiRowCut &, CbcNodeInfo *, int whichOne,
		   int whichGenerator=-1);  

    /** Destructor

      \note The destructor will reach out (via #owner_) and NULL the
      reference to the cut in the owner's
      \link CbcNodeInfo::cuts_ cuts_ \endlink list.
    */
    virtual ~CbcCountRowCut ();
  //@}

  /// Increment the number of references
  void increment(int change=1);

  /// Decrement the number of references and return the number left.
  int decrement(int change=1);

  /** \brief Set the information associating this cut with a node

    An CbcNodeInfo object and an index in the cut set of the node.
    For locally valid cuts, the node will be the  search tree node where the
    cut was generated. For globally valid cuts, it's the node where the cut
    was activated.
  */
  void setInfo(CbcNodeInfo *, int whichOne);

  /// Number of other CbcNodeInfo objects pointing to this row cut
  inline int numberPointingToThis()
  { return numberPointingToThis_;};

  /// Which generator for cuts - as user order
   int whichCutGenerator() const
  { return whichCutGenerator_;};

#ifdef CHECK_CUT_COUNTS
  // Just for printing sanity checks
  int tempNumber_;
#endif

private:

  /// Standard copy is illegal (reference counts would be incorrect)
  CbcCountRowCut(const CbcCountRowCut &);
  
  /// Standard assignment is illegal (reference counts would be incorrect)
  CbcCountRowCut & operator=(const CbcCountRowCut& rhs);

  /// Backward pointer to owning CbcNodeInfo
  CbcNodeInfo * owner_;
  
  /// Index of cut in owner's cut set
  /// (\link CbcNodeInfo::cuts_ cuts_ \endlink).
  int ownerCut_;
  
  /// Number of other CbcNodeInfo objects pointing to this cut
  int numberPointingToThis_;

  /// Which generator created this cut
  int whichCutGenerator_;

};

#endif
