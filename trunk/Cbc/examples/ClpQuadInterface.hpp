// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef ClpQuadInterface_H
#define ClpQuadInterface_H

#include "OsiClpSolverInterface.hpp"

//#############################################################################

/**

    This is to allow the user to replace initialSolve and resolve
*/

class ClpQuadInterface : public OsiClpSolverInterface {

public:
  //---------------------------------------------------------------------------
  /**@name Solve methods */
  //@{
    /// Solve initial LP relaxation 
    virtual void initialSolve();

    /// Resolve an LP relaxation after problem modification
    virtual void resolve();

  //@}


  /**@name Constructors and destructors */
  //@{
    /// Default Constructor
    ClpQuadInterface ();
    
    /// Clone
    virtual OsiSolverInterface * clone(bool CopyData=true) const;
    
    /// Copy constructor 
    ClpQuadInterface (const ClpQuadInterface &);
    
    /// Assignment operator 
    ClpQuadInterface & operator=(const ClpQuadInterface& rhs);
    
    /// Destructor 
    virtual ~ClpQuadInterface ();

  //@}


  /**@name Sets and Getss */
  //@{
  /** Setup fake objective.  It could also read an ampl .nl file
      or somehow get nonlinear  info */
  void initialize();
  /// Get objective function value (can't use default)
  virtual double getObjValue() const;

  //@}

  //---------------------------------------------------------------------------

private:
  
  /**@name Private member data */
  //@{
  /// True quadratic objective
  ClpObjective * quadraticObjective_;
  //@}
};

#endif
