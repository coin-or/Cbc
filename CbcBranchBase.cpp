// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cmath>
#include <cfloat>

#include "OsiSolverInterface.hpp"
#include "OsiSolverBranch.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcBranchBase.hpp"


// Default Constructor
CbcObject::CbcObject() 
  :model_(NULL),
   id_(-1),
   priority_(1000)
{
}

// Constructor from model
CbcObject::CbcObject(CbcModel * model)
:
  model_(model),
  id_(-1),
  priority_(1000)
{
}


// Destructor 
CbcObject::~CbcObject ()
{
}

// Copy constructor 
CbcObject::CbcObject ( const CbcObject & rhs)
{
  model_ = rhs.model_;
  id_ = rhs.id_;
  priority_ = rhs.priority_;
}

// Assignment operator 
CbcObject & 
CbcObject::operator=( const CbcObject& rhs)
{
  if (this!=&rhs) {
    model_ = rhs.model_;
    id_ = rhs.id_;
    priority_ = rhs.priority_;
  }
  return *this;
}

/* Returns floor and ceiling i.e. closest valid points
 */
void 
CbcObject::floorCeiling(double & floorValue, double & ceilingValue, double value,
			double tolerance) const
{
  if (fabs(floor(value+0.5)-value)>tolerance) {
    floorValue = floor(value);
  } else {
    floorValue = floor(value+0.5);
  }
  ceilingValue = floorValue+1.0;
}
// Return "up" estimate (default 1.0e-5)
double 
CbcObject::upEstimate() const
{
  return 1.0e-5;
}
// Return "down" estimate (default 1.0e-5)
double 
CbcObject::downEstimate() const
{
  return 1.0e-5;
}
// Column number if single column object -1 otherwise
int 
CbcObject::columnNumber() const
{
  // Say not as at present only used by heuristics
  return -1;
}
/* Create an OsiSolverBranch object
   
This returns NULL if branch not represented by bound changes
*/
OsiSolverBranch * 
CbcObject::solverBranch() const
{
  return NULL;
}
  
// Default Constructor 
CbcBranchingObject::CbcBranchingObject()
{
  model_=NULL;
  variable_=-1;
  way_=0;
  value_=0.0;
  numberBranchesLeft_=2;
}

// Useful constructor
CbcBranchingObject::CbcBranchingObject (CbcModel * model, int variable, int way , double value)
{
  model_= model;
  variable_=variable;
  way_=way;
  value_=value;
  numberBranchesLeft_=2;
}

// Copy constructor 
CbcBranchingObject::CbcBranchingObject ( const CbcBranchingObject & rhs)
{
  model_=rhs.model_;
  variable_=rhs.variable_;
  way_=rhs.way_;
  value_=rhs.value_;
  numberBranchesLeft_=rhs.numberBranchesLeft_;
}

// Assignment operator 
CbcBranchingObject & 
CbcBranchingObject::operator=( const CbcBranchingObject& rhs)
{
  if (this != &rhs) {
    model_=rhs.model_;
    variable_=rhs.variable_;
    way_=rhs.way_;
    value_=rhs.value_;
    numberBranchesLeft_=rhs.numberBranchesLeft_;
  }
  return *this;
}

// Destructor 
CbcBranchingObject::~CbcBranchingObject ()
{
}
// Default Constructor 
CbcBranchDecision::CbcBranchDecision ()
  : object_(NULL)
{
}

CbcBranchDecision::~CbcBranchDecision()
{
  delete object_;
}
/* Compare N branching objects. Return index of best
   and sets way of branching in chosen object.
   
   This routine is used only after strong branching.
   This is reccommended version as it can be more sophisticated
*/

int
CbcBranchDecision::bestBranch (CbcBranchingObject ** objects, int numberObjects,
			       int numberUnsatisfied,
			       double * changeUp, int * numberInfeasibilitiesUp,
			       double * changeDown, int * numberInfeasibilitiesDown,
			       double objectiveValue) 
{
  int bestWay=0;
  int whichObject = -1;
  if (numberObjects) {
    initialize(objects[0]->model()); 
    CbcBranchingObject * bestObject = NULL;
    for (int i = 0 ; i < numberObjects ; i++) {
      int betterWay = betterBranch(objects[i],
				   bestObject,
				   changeUp[i],
				   numberInfeasibilitiesUp [i],
				   changeDown[i],
				   numberInfeasibilitiesDown[i] );
      if (betterWay) {
	bestObject = objects[i];
	bestWay = betterWay;
	whichObject=i;
      }
    }
    // set way in best
    if (whichObject>=0)
      objects[whichObject]->way(bestWay);
    else
      printf("debug\n");
  }
  return whichObject;
}

// Default constructor
CbcConsequence::CbcConsequence()
{
}


// Destructor 
CbcConsequence::~CbcConsequence ()
{
}

// Copy constructor 
CbcConsequence::CbcConsequence ( const CbcConsequence & rhs)
{
}

// Assignment operator 
CbcConsequence & 
CbcConsequence::operator=( const CbcConsequence& rhs)
{
  if (this!=&rhs) {
  }
  return *this;
}

  
