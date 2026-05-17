// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglLiftAndProject_H
#define CglLiftAndProject_H

#include <string>

#include "CglCutGenerator.hpp"

/** Lift And Project Cut Generator Class */
class CGLLIB_EXPORT CglLiftAndProject : public CglCutGenerator {
   friend CGLLIB_EXPORT void CglLiftAndProjectUnitTest(const OsiSolverInterface * siP,
					const std::string mpdDir );

public:
  /**@name Generate Cuts */
  //@{
  /** Generate lift-and-project cuts for the 
      model of the solver interface, si. 
      Insert the generated cuts into OsiCut, cs.
  */
  virtual void generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			    const CglTreeInfo info = CglTreeInfo());

  /** Get the normalization : Either beta=+1 or beta=-1.
  */

  double getBeta() const {
    return beta_;
  }

  /** Set the normalization : Either beta=+1 or beta=-1.
      Default value is 1.
  */
  void setBeta(int oneOrMinusOne){
    if (oneOrMinusOne==1 || oneOrMinusOne==-1){
      beta_= static_cast<double>(oneOrMinusOne);
    }
    else {
      throw CoinError("Unallowable value. Beta must be 1 or -1",
		      "cutGeneration","CglLiftAndProject");
    }
  }

  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglLiftAndProject ();
 
  /// Copy constructor 
  CglLiftAndProject (
    const CglLiftAndProject &);

  /// Clone
  virtual CglCutGenerator * clone() const;

  /// Assignment operator 
  CglLiftAndProject &
    operator=(
    const CglLiftAndProject& rhs);
  
  /// Destructor 
  virtual
    ~CglLiftAndProject ();
  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);
  //@}

private:
  
 // Private member methods

  /**@name Private methods */
  //@{

  //@}

  // Private member data

  /**@name Private member data */
  //@{
  /// The normalization is beta_=1 or beta_=-1
  double beta_;  
  /// epsilon
  double epsilon_;  
  /// 1-epsilon
  double onetol_;  
  //@}
};

//#############################################################################
/** A function that tests the methods in the CglLiftAndProject class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
CGLLIB_EXPORT
void CglLiftAndProjectUnitTest(const OsiSolverInterface * siP,
			      const std::string mpdDir );
  
#endif
