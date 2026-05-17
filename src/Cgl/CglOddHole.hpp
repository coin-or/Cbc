// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglOddHole_H
#define CglOddHole_H

#include <string>

#include "CglCutGenerator.hpp"

/** Odd Hole Cut Generator Class */
class CGLLIB_EXPORT CglOddHole : public CglCutGenerator {
   friend CGLLIB_EXPORT void CglOddHoleUnitTest(const OsiSolverInterface * siP,
				  const std::string mpdDir );
 
public:
    
  
  /**@name Generate Cuts */
  //@{
  /** Generate odd hole cuts for the model of the solver interface, si.
      This looks at all rows of type sum x(i) <= 1 (or == 1) (x 0-1)
      and sees if there is an odd cycle cut.  See Grotschel, Lovasz
      and Schrijver (1988) for method.
      This is then lifted by using the corresponding Chvatal cut i.e.
      Take all rows in cycle and add them together. RHS will be odd so  
      weaken all odd coefficients so 1.0 goes to 0.0 etc - then
      constraint is  sum even(j)*x(j) <= odd which can be replaced by
      sum (even(j)/2)*x(j) <= (odd-1.0)/2.
      A similar cut can be generated for sum x(i) >= 1.

      Insert the generated cuts into OsiCut, cs.

      This is only done for rows with unsatisfied 0-1 variables.  If there
      are many of these it will be slow.  Improvements would do a 
      randomized subset and also speed up shortest path algorithm used.

  */
  virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
			     const CglTreeInfo info = CglTreeInfo());
  //@}

  /**@name Create Row List */
  //@{
  /// Create a list of rows which might yield cuts
  /// this is to speed up process
  /// The possible parameter is a list to cut down search
  void createRowList( const OsiSolverInterface & si,
		      const int * possible=NULL);
  /// This version passes in a list - 1 marks possible
  void createRowList(int numberRows, const int * whichRow);
  //@}

  /**@name Create Clique List */
  //@{
  /// Create a list of extra row cliques which may not be in matrix
  /// At present these are classical cliques
  void createCliqueList(int numberCliques, const int * cliqueStart,
		     const int * cliqueMember);
  //@}

  /**@name Number Possibilities */
  //@{
  /// Returns how many rows might give odd hole cuts
  int numberPossible();
  //@}
  /**@name Gets and Sets */
  //@{
  /// Minimum violation
  double getMinimumViolation() const;
  void setMinimumViolation(double value);
  /// Minimum violation per entry
  double getMinimumViolationPer() const;
  void setMinimumViolationPer(double value);
  /// Maximum number of entries in a cut
  int getMaximumEntries() const;
  void setMaximumEntries(int value);
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglOddHole ();
 
  /// Copy constructor 
  CglOddHole (
    const CglOddHole &);

  /// Clone
  virtual CglCutGenerator * clone() const;

  /// Assignment operator 
  CglOddHole &
    operator=(
    const CglOddHole& rhs);
  
  /// Destructor 
  virtual
    ~CglOddHole ();

  /// This can be used to refresh any inforamtion
  virtual void refreshSolver(OsiSolverInterface * solver);
  //@}
      
private:
  
 // Private member methods


  /**@name Private methods */
  //@{
  /// Generate cuts from matrix copy and solution
  /// If packed true then <=1 rows, otherwise >=1 rows.
  void generateCuts(const OsiRowCutDebugger * debugger, 
		    const CoinPackedMatrix & rowCopy,
		    const double * solution, const double * dj,
		    OsiCuts & cs, const int * suitableRow,
		    const int * fixedColumn,const CglTreeInfo info,
		    bool packed);
  //@}

  // Private member data

  /**@name Private member data */
  //@{
  /// list of suitableRows
  int * suitableRows_;
  /// start of each clique
  int * startClique_;
  /// clique members
  int * member_;
  /// epsilon
  double epsilon_;  
  /// 1-epsilon
  double onetol_;
  /// Minimum violation
  double minimumViolation_;
  /// Minimum violation per entry
  double minimumViolationPer_;
  /// Maximum number of entries in a cut
  int maximumEntries_;
  /// number of rows when suitability tested
  int numberRows_;
  /// number of cliques
  int numberCliques_;
  //@}
};

//#############################################################################
/** A function that tests the methods in the CglOddHole class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
CGLLIB_EXPORT
void CglOddHoleUnitTest(const OsiSolverInterface * siP,
			const std::string mpdDir );
  
#endif
