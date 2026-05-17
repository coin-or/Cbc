// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglKnapsackCover_H
#define CglKnapsackCover_H

#include <string>

#include "CglCutGenerator.hpp"
#include "CglTreeInfo.hpp"

/** Knapsack Cover Cut Generator Class */
class CGLLIB_EXPORT CglKnapsackCover : public CglCutGenerator {
   friend CGLLIB_EXPORT void CglKnapsackCoverUnitTest(const OsiSolverInterface * siP,
					const std::string mpdDir );

public:
   /** A method to set which rows should be tested for knapsack covers */
   void setTestedRowIndices(int num, const int* ind);

   /**@name Generate Cuts */
  //@{
  /** Generate knapsack cover cuts for the model of the solver interface, si. 
      Insert the generated cuts into OsiCut, cs.
  */
  virtual void generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			    const CglTreeInfo info = CglTreeInfo());
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglKnapsackCover ();
 
  /// Copy constructor 
  CglKnapsackCover (
    const CglKnapsackCover &);

  /// Clone
  virtual CglCutGenerator * clone() const;

  /// Assignment operator 
  CglKnapsackCover &
    operator=(
    const CglKnapsackCover& rhs);
  
  /// Destructor 
  virtual
    ~CglKnapsackCover ();
  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);
  /// This can be used to refresh any information
  virtual void refreshSolver(OsiSolverInterface * solver);
  //@}


  /**@name Sets and gets */
  //@{
  /// Set limit on number in knapsack
  inline void setMaxInKnapsack(int value)
           { if (value>0) maxInKnapsack_ = value;}
  /// get limit on number in knapsack
  inline int getMaxInKnapsack() const
           {return maxInKnapsack_;}
  /// Switch off expensive cuts
  inline void switchOffExpensive()
  { expensiveCuts_=false;}
  /// Switch on expensive cuts
  inline void switchOnExpensive()
  { expensiveCuts_=true;}
private:
  
 // Private member methods


  /**@name Private methods */
  //@{

  /** deriveAKnapsack 
                 returns 1 if it is able to derive
                 a (canonical) knapsack inequality
                in binary variables of the form ax<=b 
                 from the rowIndex-th  row in the model, 
                returns 0 otherwise.
  */
  int deriveAKnapsack(
    const OsiSolverInterface & si, 
    OsiCuts & cs,
    CoinPackedVector & krow,
    bool treatAsLRow,
    double & b,
    int *  complement,
    double *  xstar,
    int rowIndex,
    int numberElements,
    const int * index,
    const double * element);

  int deriveAKnapsack(
    const OsiSolverInterface & si, 
    OsiCuts & cs,
    CoinPackedVector & krow,
    double & b,
    int *  complement,
    double *  xstar,
    int rowIndex,
    const CoinPackedVectorBase & matrixRow);

  /** Find a violated minimal cover from 
 a canonical form knapsack inequality by
 solving the -most- violated cover problem
 and postprocess to ensure minimality
  */
  int findExactMostViolatedMinCover(
      int nCols, 
      int row,
      CoinPackedVector & krow,
      double b, 
      double *  xstar, 
      CoinPackedVector & cover,
      CoinPackedVector & remainder);

  /** Find the most violate minimum cover by solving the lp-relaxation of the
      most-violate-min-cover problem 
  */
  int findLPMostViolatedMinCover(
      int nCols,
      int row,
      CoinPackedVector & krow,
      double & b,
      double * xstar, 
      CoinPackedVector & cover,
      CoinPackedVector & remainder);
  
/// find a minimum cover by a simple greedy approach
  int findGreedyCover(
      int row,
      CoinPackedVector & krow,
      double & b,
      double * xstar,
      CoinPackedVector & cover,
      CoinPackedVector & remainder
      );

  /// lift the cover inequality
  int liftCoverCut(
     double & b,
     int nRowElem,
     CoinPackedVector & cover,
     CoinPackedVector & remainder,
     CoinPackedVector & cut );
 
  /// sequence-independent lift and uncomplement and add the resulting cut to the cut set
  int liftAndUncomplementAndAdd(
     double rowub,
     CoinPackedVector & krow,
     double & b,
     int * complement,
     int row,
     CoinPackedVector & cover,
     CoinPackedVector & remainder,
     OsiCuts & cs );

  /// sequence-dependent lift, uncomplement and add the resulting cut to the cut set
void seqLiftAndUncomplementAndAdd(
      int nCols,
      double * xstar, 
      int * complement,
      int row,
      int nRowElem,
      double & b,
      CoinPackedVector & cover,      // need not be violated
      CoinPackedVector & remainder,
      OsiCuts & cs );

  /// sequence-dependent lift binary variables either up or down, uncomplement and add to the cut set
void liftUpDownAndUncomplementAndAdd(
         int nCols,
         double * xstar, 
         int * complement,
         int row,
         int nRowElem,
         double & b,

         // the following 3 packed vectors partition the krow:
         CoinPackedVector & fracCover, // vars have frac soln values in lp relaxation
                                       // and form cover with the vars atOne
         CoinPackedVector & atOne,     // vars have soln value of 1 in lp relaxation
                                       // and together with fracCover form minimal (?) cover. 
         CoinPackedVector & remainder,
         OsiCuts & cs );

  /// find a cover using a variation of the logic found in OSL (w/o SOS)
  int findPseudoJohnAndEllisCover (
     int row,
     CoinPackedVector & krow,
     double & b,
     double * xstar,                     
     CoinPackedVector & cover,  
     CoinPackedVector & remainder);

  /// find a cover using the basic logic found in OSL (w/o SOS)
  int findJohnAndEllisCover (
     int row,
     CoinPackedVector & krow,
     double & b,
     double * xstar,                     
     CoinPackedVector & fracCover,  
     CoinPackedVector & atOnes,  
     CoinPackedVector & remainder);


  /** A C-style implementation of the Horowitz-Sahni exact solution 
   procedure for solving knapsack problem. 
   
   (ToDo: implement the more efficient dynamic programming approach)

   (Reference: Martello and Toth, Knapsack Problems, Wiley, 1990, p30.)
  */
  int exactSolveKnapsack(
      int n, 
      double c, 
      double const *pp, 
      double const *ww,
      double & z, 
      int * x);
  /// For testing gub stuff
  int gubifyCut(CoinPackedVector & cut);
public:
  /** Creates cliques for use by probing.
      Only cliques >= minimumSize and < maximumSize created
      Can also try and extend cliques as a result of probing (root node).
      Returns number of cliques found.
  */
  int createCliques( OsiSolverInterface & si, 
		    int minimumSize=2, int maximumSize=100, bool extendCliques=false);
private:
  /// Delete all clique information
  void deleteCliques();
  //@}

  // Private member data

  /**@name Private member data */
  //@{
  /// epsilon
  double epsilon_;  
  /// Tolerance to use for violation - bigger than epsilon_
  double epsilon2_;
  /// 1-epsilon
  double onetol_;  
  /// Maximum in knapsack
  int maxInKnapsack_;
  /** which rows to look at. If specified, only these rows will be considered
      for generating knapsack covers. Otherwise all rows will be tried */
  int numRowsToCheck_;
  int* rowsToCheck_;
  /// exactKnapsack can be expensive - this switches off some
  bool expensiveCuts_;
  /// Cliques
  /// **** TEMP so can reference from listing
  const OsiSolverInterface * solver_;
  int whichRow_;
  int * complement_;
  double * elements_;
  /// Number of cliques
  int numberCliques_;
  /// Clique type
  typedef struct {
    unsigned int equality:1; //  nonzero if clique is ==
  } CliqueType;
  CliqueType * cliqueType_;
  /// Start of each clique
  int * cliqueStart_;
  /// Entries for clique
  CliqueEntry * cliqueEntry_;
  /** Start of oneFixes cliques for a column in matrix or -1 if not
      in any clique */
  int * oneFixStart_;
  /** Start of zeroFixes cliques for a column in matrix or -1 if not
      in any clique */
  int * zeroFixStart_;
  /// End of fixes for a column
  int * endFixStart_;
  /// Clique numbers for one or zero fixes
  int * whichClique_;
  /// Number of columns
  int numberColumns_;
  /** For each column with nonzero in row copy this gives a clique "number".
      So first clique mentioned in row is always 0.  If no entries for row
      then no cliques.  If sequence > numberColumns then not in clique.
  */
  //CliqueEntry * cliqueRow_;
  /// cliqueRow_ starts for each row
  //int * cliqueRowStart_;
  //@}
};

//#############################################################################
/** A function that tests the methods in the CglKnapsackCover class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
CGLLIB_EXPORT
void CglKnapsackCoverUnitTest(const OsiSolverInterface * siP,
			      const std::string mpdDir );
  
#endif
