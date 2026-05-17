// Copyright (C) 2010, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).
#ifndef CglZeroHalf_H
#define CglZeroHalf_H

#include <string>

#include "CglCutGenerator.hpp"
#include "CoinPackedMatrix.hpp"
#include "Cgl012cut.hpp" 

/** Zero Half Cut Generator Class

 This class generates zero half cuts via the following method:

 See - 

G. Andreello, A. Caprara, M. Fischetti,
 “Embedding Cuts in a Branch and Cut Framework: a Computational Study 
  with {0,1/2}-Cuts”, INFORMS Journal on Computing 19(2), 229-238, 2007.
 
*/

class CGLLIB_EXPORT CglZeroHalf : public CglCutGenerator {
   friend CGLLIB_EXPORT void CglZeroHalfUnitTest(const OsiSolverInterface * siP,
					 const std::string mpdDir );
 
public:

  /**@name Generate Cuts */
  //@{
  /** Generate zero half cuts for the model accessed through the solver interface. 
  Insert generated cuts into the cut set cs.
  */
  virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
			     const CglTreeInfo info = CglTreeInfo());
  //@}

  /**@name Sets and Gets */
  //@{
  /// Get flags
  inline int getFlags() const
  { return flags_;}
  /// Set flags
  inline void setFlags(int value)
  { flags_ = value;}
  /// Get active-node threshold for sparse separation graph use (0 means always sparse)
  inline int getSepGraphSparseThreshold() const
  { return cutInfo_.getSepGraphSparseThreshold(); }
  /// Set active-node threshold for sparse separation graph use (0 means always sparse)
  inline void setSepGraphSparseThreshold(int value)
  { cutInfo_.setSepGraphSparseThreshold(value); }
  /// Get maximum row pair count allowed before skipping the row (-1 disables)
  inline int getRowMaxPairCount() const
  { return cutInfo_.getRowMaxPairCount(); }
  /// Set maximum row pair count allowed before skipping the row (-1 disables)
  inline void setRowMaxPairCount(int value)
  { cutInfo_.setRowMaxPairCount(value); }
  /// Get maximum row fractional count allowed before skipping the row (-1 disables)
  inline int getRowMaxFractionalCount() const
  { return cutInfo_.getRowMaxFractionalCount(); }
  /// Set maximum row fractional count allowed before skipping the row (-1 disables)
  inline void setRowMaxFractionalCount(int value)
  { cutInfo_.setRowMaxFractionalCount(value); }
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglZeroHalf ();
 
  /// Copy constructor 
  CglZeroHalf (
    const CglZeroHalf &);

  /// Clone
  virtual CglCutGenerator * clone() const;

  /// Assignment operator 
  CglZeroHalf &
    operator=(
    const CglZeroHalf& rhs);
  
  /// Destructor 
  virtual
    ~CglZeroHalf ();
  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);
  /// This can be used to refresh any information
  virtual void refreshSolver(OsiSolverInterface * solver);
  //@}

private:
  
  // Private member methods
   
  /**@name Private methods */
  //@{
  //@}
  
  
  /**@name Private member data */
  //@{
  /// number of rows in the ILP matrix 
  int mr_;
  /// number of columns in the ILP matrix 
  int mc_;
  /// number of nonzero's in the ILP matrix 
  int mnz_;
  /// starting position of each row in arrays mtind and mtval 
  int *mtbeg_;
  /// number of entries of each row in arrays mtind and mtval 
  int *mtcnt_;
  /// column indices of the nonzero entries of the ILP matrix 
  int *mtind_;
  /// values of the nonzero entries of the ILP matrix 
  int *mtval_;
  /// lower bounds on the variables 
  int *vlb_;
  /// upper bounds on the variables 
  int *vub_;
  /// right hand sides of the constraints 
  int *mrhs_;
  /// senses of the constraints: 'L', 'G' or 'E' 
  char *msense_;
  /// Cgl012Cut object to make thread safe
  Cgl012Cut cutInfo_;
  /** Flags
      1 bit - global cuts 
  */
  int flags_;
  //@}
};
/// A simple Dijkstra shortest path - make better later
CGLLIB_EXPORT
#ifndef CGL_NEW_SHORT
void cglShortestPath(cgl_graph * graph, int source, int maximumLength);
#else
void cglShortestPath(auxiliary_graph * graph, int source, int maximumLength);
#endif
//#############################################################################
/** A function that tests the methods in the CglZeroHalf class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
CGLLIB_EXPORT
void CglZeroHalfUnitTest(const OsiSolverInterface * siP,
			       const std::string mpdDir );
  
#endif
