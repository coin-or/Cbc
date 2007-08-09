// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
   
#include "CbcConfig.h"
#include "CoinPragma.hpp"
#include "CbcModel.hpp"
#include "OsiClpSolverInterface.hpp"

#include <cassert>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <iostream>
  // define TEST_MESSAGE_HANDLER to check works on all messages
//#define TEST_MESSAGE_HANDLER
#ifdef TEST_MESSAGE_HANDLER
// This driver shows how to trap messages - this is just as in unitTest.cpp
// ****** THis code is similar to MyMessageHandler.hpp and MyMessagehandler.cpp
#include "CoinMessageHandler.hpp"

/** This just adds a model to CoinMessage and a void pointer so
    user can trap messages and do useful stuff.  
    This is used in Clp/Test/unitTest.cpp

    The file pointer is just there as an example of user stuff.

*/
class CbcModel;

class MyMessageHandler2 : public CoinMessageHandler {
  
public:
  /**@name Overrides */
  //@{
  virtual int print();
  //@}
  /**@name set and get */
  //@{
  /// Model
  const CbcModel * model() const;
  void setModel(CbcModel * model);
  //@}

  /**@name Constructors, destructor */
  //@{
  /** Default constructor. */
  MyMessageHandler2();
  /// Constructor with pointer to model
  MyMessageHandler2(CbcModel * model,
			   FILE * userPointer=NULL);
  /** Destructor */
  virtual ~MyMessageHandler2();
  //@}

  /**@name Copy method */
  //@{
  /** The copy constructor. */
  MyMessageHandler2(const MyMessageHandler2&);
  /** The copy constructor from an CoinSimplexMessageHandler. */
  MyMessageHandler2(const CoinMessageHandler&);
  
  MyMessageHandler2& operator=(const MyMessageHandler2&);
  /// Clone
  virtual CoinMessageHandler * clone() const ;
  //@}
   
    
protected:
  /**@name Data members
     The data members are protected to allow access for derived classes. */
  //@{
  /// Pointer back to model
  CbcModel * model_;
  //@}
};


//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
MyMessageHandler2::MyMessageHandler2 () 
  : CoinMessageHandler(),
    model_(NULL)
{
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
MyMessageHandler2::MyMessageHandler2 (const MyMessageHandler2 & rhs) 
: CoinMessageHandler(rhs),
    model_(rhs.model_)
{  
}

MyMessageHandler2::MyMessageHandler2 (const CoinMessageHandler & rhs) 
  : CoinMessageHandler(),
    model_(NULL)
{  
}

// Constructor with pointer to model
MyMessageHandler2::MyMessageHandler2(CbcModel * model,
               FILE * userPointer)
  : CoinMessageHandler(),
    model_(model)
{
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
MyMessageHandler2::~MyMessageHandler2 ()
{
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
MyMessageHandler2 &
MyMessageHandler2::operator=(const MyMessageHandler2& rhs)
{
  if (this != &rhs) {
    CoinMessageHandler::operator=(rhs);
    model_ = rhs.model_;
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CoinMessageHandler * MyMessageHandler2::clone() const
{
  return new MyMessageHandler2(*this);
}
int 
MyMessageHandler2::print()
{
  // Just add ==
  fprintf(fp_," == ");
  fprintf(fp_,"%s\n",messageBuffer_);
  return 0;
}
const CbcModel *
MyMessageHandler2::model() const
{
  return model_;
}
void 
MyMessageHandler2::setModel(CbcModel * model)
{
  model_ = model;
}
#endif
//#############################################################################
// To use USERCBC or USERCLP change 0 to 1 in defines and add in your fake main program(s) and any other code
//#define USER_HAS_FAKE_CBC 
//#define USER_HAS_FAKE_CLP 
#ifdef USER_HAS_FAKE_CBC
#endif
void fakeMain (ClpSimplex & model,OsiSolverInterface & osiSolver, CbcModel & babSolver)
{
#ifdef USER_HAS_FAKE_CBC
#else
  printf("Dummy user cbc code - model has %d rows and %d columns\n",
	 model.getNumRows(),model.getNumCols());
  // Reduce printout
  babSolver.solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
  // Do complete search
  babSolver.branchAndBound();
#endif
}
// Clp stuff
#ifdef USER_HAS_FAKE_CLP
#endif
void fakeMain2 (ClpSimplex & model,OsiClpSolverInterface & osiSolver,int options) {
#ifdef USER_HAS_FAKE_CLP
#else
  ClpSimplex * lpSolver = osiSolver.getModelPtr();
  printf("Dummy user clp code - model has %d rows and %d columns\n",
	 lpSolver->numberRows(),lpSolver->numberColumns());
  osiSolver.initialSolve();
#endif
}
//  End any fake main program
//#############################################################################
void CbcClpUnitTest (const CbcModel & saveModel);
int main (int argc, const char *argv[])
{
  OsiClpSolverInterface solver1;
  CbcModel model(solver1);
  // define TEST_MESSAGE_HANDLER at top of file to check works on all messages
#ifdef TEST_MESSAGE_HANDLER
  MyMessageHandler2 messageHandler(&model);
  std::cout<<"Testing derived message handler"<<std::endl;
  model.passInMessageHandler(&messageHandler);
  OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (model.solver());
  // Could use different handlers (if different log levels)
  clpSolver->passInMessageHandler(&messageHandler);
  //clpSolver->getModelPtr()->passInMessageHandler(&messageHandler);
#endif
  // initialize
  CbcMain0(model);
#ifdef TEST_MESSAGE_HANDLER
  // Set log levels same so can use one message handler
  clpSolver->messageHandler()->setLogLevel(1) ;
  model.messageHandler()->setLogLevel(1);
  // switch off some printing
  void setCbcOrClpPrinting(bool yesNo);
  setCbcOrClpPrinting(false);
#endif
  int returnCode = CbcMain1 (argc, argv,model);
  if (returnCode!=777) {
    //const CoinWarmStartBasis* debugws =
    //dynamic_cast <const CoinWarmStartBasis*>(model.solver()->getWarmStart()) ;
    //debugws->print() ;
    //delete debugws ; 
    return returnCode;
  } else {
    // do miplib
    CbcClpUnitTest(model);
    return 0;
  }
}
/*
  Version 1.00.00 November 16 2005.
  This is to stop me (JJF) messing about too much.
  Tuning changes should be noted here.
  The testing next version may be activated by CBC_NEXT_VERSION
  This applies to OsiClp, Clp etc
  Version 1.00.01 November 24 2005
  Added several classes for advanced users.  This can't affect code (if you don't use it)
  Made some tiny changes (for N way branching) which should not change anything.
  CbcNWay object class - for N way branching this also allows use of CbcConsequence class.
  CbcBranchAllDifferent object class - for branching on general integer variables
  to stop them having same value so branches are x >= y+1 and x <= y-1.
  Added two new Cgl classes - CglAllDifferent which does column fixing (too slowly)
  and CglStored which just has a list of cuts which can be activated.
  Modified preprocess option to SOS
  Version 1.00.02 December 9 2005
  Added use of CbcStrategy to do clean preprocessing
  Added use of referenceSolver for cleaner repetition of Cbc
  Version 1.01.00 February 2 2006
  Added first try at Ampl interface
  Made dummy program so real main can be called from other programs
*/
