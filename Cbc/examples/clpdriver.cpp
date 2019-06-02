// $Id: clpdriver.cpp 2101 2014-12-03 17:43:20Z forrest $
// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>
#include <iomanip>

#include "CoinPragma.hpp"
#include "CbcModel.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CbcSolver.hpp"

#include "CoinTime.hpp"

//#############################################################################

/************************************************************************

This main program shows how to take advantage of the standalone cbc in your program,
while still making major modifications.  
Clp has not got a similar function so there can be times when it is
useful to use Cbc for Lps.  Normally if a postsolved problem has small
infeasibilities then simplex is called to clean up.  For some problems
you may wish to kip this step (e.g. to save memory).
First it reads in a model from an mps file
Then it initializes the model with cbc defaults
Then it calls CbcMain1 passing all parameters apart from first but with 
an event handler to modify stuff.

************************************************************************/
#include "ClpEventHandler.hpp"
/** This is so user can trap events and do useful stuff.  

    ClpSimplex model_ is available as well as anything else you care 
    to pass in
*/

class MyEventHandler3 : public ClpEventHandler {

public:
  /**@name Overrides */
  //@{
  virtual int event(Event whichEvent);
  //@}

  /**@name Constructors, destructor etc*/
  //@{
  /** Default constructor. */
  MyEventHandler3();
  /// Constructor with pointer to model (redundant as setEventHandler does)
  MyEventHandler3(ClpSimplex *model);
  /** Destructor */
  virtual ~MyEventHandler3();
  /** The copy constructor. */
  MyEventHandler3(const MyEventHandler3 &rhs);
  /// Assignment
  MyEventHandler3 &operator=(const MyEventHandler3 &rhs);
  /// Clone
  virtual ClpEventHandler *clone() const;
  //@}

protected:
  // data goes here
};
//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
MyEventHandler3::MyEventHandler3()
  : ClpEventHandler()
{
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
MyEventHandler3::MyEventHandler3(const MyEventHandler3 &rhs)
  : ClpEventHandler(rhs)
{
}

// Constructor with pointer to model
MyEventHandler3::MyEventHandler3(ClpSimplex *model)
  : ClpEventHandler(model)
{
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
MyEventHandler3::~MyEventHandler3()
{
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
MyEventHandler3 &
MyEventHandler3::operator=(const MyEventHandler3 &rhs)
{
  if (this != &rhs) {
    ClpEventHandler::operator=(rhs);
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpEventHandler *MyEventHandler3::clone() const
{
  return new MyEventHandler3(*this);
}

int MyEventHandler3::event(Event whichEvent)
{
  // See if just after postsolve
  if (whichEvent == presolveAfterFirstSolve || whichEvent == presolveAfterSolve) {
    return -2; // skip clean up
  } else {
    return -1; // carry on
  }
}

static int dummyCallBack(CbcModel * /*model*/, int /*whereFrom*/)
{
  return 0;
}

int main(int argc, const char *argv[])
{

  OsiClpSolverInterface solver1;
  // Say we are keeping names (a bit slower this way)
  solver1.setIntParam(OsiNameDiscipline, 1);
  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string mpsFileName;
#if defined(SAMPLEDIR)
  mpsFileName = SAMPLEDIR "/p0033.mps";
#else
  if (argc < 2) {
    fprintf(stderr, "Do not know where to find sample MPS files.\n");
    exit(1);
  }
#endif
  if (argc >= 2)
    mpsFileName = argv[1];
  int numMpsReadErrors;
  if (strstr(mpsFileName.c_str(), ".mps"))
    numMpsReadErrors = solver1.readMps(mpsFileName.c_str(), "");
  else
    numMpsReadErrors = solver1.readLp(mpsFileName.c_str());
  if (numMpsReadErrors != 0) {
    printf("%d errors reading MPS file\n", numMpsReadErrors);
    return numMpsReadErrors;
  }
  ClpSimplex *model = solver1.getModelPtr();
  // Event handler
  MyEventHandler3 eventHandler(model);
  model->passInEventHandler(&eventHandler);
  // Messy code below copied from CbcSolver.cpp
  // Pass to Cbc initialize defaults
  CbcModel modelA(solver1);
  OsiClpSolverInterface *solver = dynamic_cast< OsiClpSolverInterface * >(modelA.solver());
  model = solver->getModelPtr();
  CbcSolverUsefulData cbcData;
#ifndef CBC_NO_INTERRUPT
  cbcData.useSignalHandler_ = true;
#endif
  cbcData.noPrinting_ = false;
  // initialize
  CbcMain0(modelA, cbcData);
  /* Now go into code for standalone solver
     Could copy arguments and add -quit at end to be safe
     but this will do
  */
  if (argc > 2) {
    CbcMain1(argc - 1, argv + 1, modelA, dummyCallBack , cbcData);
  } else {
    const char *argv2[] = { "clpdriver", "-solve", "-quit" };
    CbcMain1(3, argv2, modelA, dummyCallBack, cbcData);
  }

  if (!model->problemStatus()) {
    std::cout << "Objective value " << model->objectiveValue() << std::endl;
  } else {
    std::cout << "Infeasible!" << std::endl;
  }
  return 0;
}
