/* $Id$ */
// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.

/*! \file CbcSolver.cpp
    \brief Main routine for the cbc stand-alone solver.
*/

#include "CbcConfig.h"
#include "CoinPragma.hpp"
#include "CbcModel.hpp"
#include "CbcOrClpParam.hpp"
#include "OsiClpSolverInterface.hpp"

/*
  We have the following compile-time symbols.

  NEW_STYLE_SOLVER	CoinSolve.cpp, CbcSolver.cpp

    Unclear what this does just yet. A value of 0 seems to be `old style
    solver'.


  CBC_OTHER_SOLVER	CoinSolve.cpp, CbcSolver.[cpp,hpp], CbcModel.cpp

    Usage in CbcSolver.hpp says `other solver' is Cplex (only).

    Here in CoinSolver, CBC_OTHER_SOLVER dominates NEW_STYLE_SOLVER.

    Usage in CbcModel is a fake; a small bit of code that's now `#if 0'.


  CPX_KEEP_RESULTS	CoinSolve.cpp, CbcSolver.cpp

    Unclear what this does just yet. The name seems clear, but how / what is
    affected is not. Defining this symbol forces CBC_OTHER_SOLVER.


  CLP_DEBUG_MALLOC

    This ties in with the functions clp_malloc, clp_free, and clp_memory,
    which are defined in CoinOslFactorization.cpp. (Right where you'd expect
    to find them, eh?).  Looks to be a relatively nice debugging wrapper for
    standard C malloc.  Calls standard C malloc/free directly if
    CLP_DEBUG_MALLOC is not defined.  Worth consideration for breaking out as
    a separate utility.  The hooks for new and delete defined here should be
    incorporated.

    Absolutely not thread safe --- lots of static variables.

    Hmmm ... is it still the case that standard C malloc and C++ new/delete
    do not play well together? 'Cause the hooks here for new and delete will
    not escape from this file.

*/


/*
  Allow (force?) use of cplex for something.
*/

#ifdef CPX_KEEP_RESULTS
#define CBC_OTHER_SOLVER 1
#endif
#if CBC_OTHER_SOLVER==1
#include "OsiCpxSolverInterface.hpp"
#endif

/*
  Hooks for a debugging wrapper for malloc/free. This bit of definition hooks
  C++ new / delete and diverts them into the debugging wrapper.
*/
//#define CLP_DEBUG_MALLOC
#ifdef CLP_DEBUG_MALLOC
/*extern "C" */void clp_memory(int type);
/*extern "C" */
void * clp_malloc(int length);
/*extern "C" */
void clp_free(void * array);
#include <malloc.h>
#include <exception>
#include <new>
void * operator new (size_t size) throw (std::bad_alloc)
{
    void * p = clp_malloc(size);
    return p;
}
void operator delete (void *p) throw()
{
    clp_free(p);
}
#endif		// CLP_DEBUG_MALLOC

#include <cassert>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <iostream>

/*
  The subtleties are unclear, but the gross action is that CBC_OTHER_SOLVER
  will overrule NEW_STYLE_SOLVER. NEW_STYLE_SOLVER undefined or 0 seems to
  mean `old style solver'.
*/
#ifndef NEW_STYLE_SOLVER
# define NEW_STYLE_SOLVER 0
#endif
#ifdef CBC_OTHER_SOLVER
# undef NEW_STYLE_SOLVER
# define NEW_STYLE_SOLVER 0
#endif



#if NEW_STYLE_SOLVER == 0
// define TEST_MESSAGE_HANDLER to check works on all messages
// #define TEST_MESSAGE_HANDLER
#ifdef TEST_MESSAGE_HANDLER
// This driver shows how to trap messages - this is just as in unitTest.cpp
// ****** THis code is similar to MyMessageHandler.hpp and MyMessagehandler.cpp
#include "CoinMessageHandler.hpp"

/** This just adds a model to CoinMessage and a void pointer so
    user can trap messages and do useful stuff.
    This is used in Clp/Test/unitTest.cpp

    The file pointer is just there as an example of user stuff.

  -- lh 071026 -- An accurate summary. Nothing is actually happening here
  except that messages will be prefixed with "==", which serves the purpose
  of demonstrating that this message handler is active. The extra parameters
  (CbcModel, FILE) are unused.

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
                      FILE * userPointer = NULL);
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
MyMessageHandler2::operator=(const MyMessageHandler2 & rhs)
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
    fprintf(fp_, " == ");
    fprintf(fp_, "%s\n", messageBuffer_);
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
#endif /* TEST_MESSAGE_HANDLER */

//#############################################################################

// To use USERCBC or USERCLP change 0 to 1 in defines and add in your fake main program(s) and any other code
//#define USER_HAS_FAKE_CBC
//#define USER_HAS_FAKE_CLP

#ifdef USER_HAS_FAKE_CBC
#endif
void fakeMain (ClpSimplex & model, OsiSolverInterface & /*osiSolver*/, CbcModel & babSolver)
{
#ifdef USER_HAS_FAKE_CBC
#else
    printf("Dummy user cbc code - model has %d rows and %d columns\n",
           model.getNumRows(), model.getNumCols());
    // Reduce printout
    babSolver.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
    // Do complete search
    babSolver.branchAndBound();
#endif
}

// Clp stuff
#ifdef USER_HAS_FAKE_CLP
#endif
void fakeMain2 (ClpSimplex & /*model*/,
                OsiClpSolverInterface & osiSolver,
                int /*options*/)
{
#ifdef USER_HAS_FAKE_CLP
#else
    ClpSimplex * lpSolver = osiSolver.getModelPtr();
    printf("Dummy user clp code - model has %d rows and %d columns\n",
           lpSolver->numberRows(), lpSolver->numberColumns());
    osiSolver.initialSolve();
#endif
}
//  End any fake main program

//#############################################################################

// void CbcClpUnitTest (const CbcModel & saveModel);

#ifdef CBC_STATISTICS
int osi_crunch = 0;
static int cbc_resolve = 0;
int osi_primal = 0;
int osi_dual = 0;
int osi_hot = 0;
void cbc_resolve_check(const OsiSolverInterface * solver)
{
    cbc_resolve++;
    printf("R %d stats %d %d %d\n",
           cbc_resolve, solver->getNumRows(), solver->getNumCols(),
           solver->getMatrixByCol()->getNumElements());
    if ((cbc_resolve % 1000) == 0)
        printf("RR %d resolve crunch %d primal %d dual %d hot %d\n",
               cbc_resolve, osi_crunch, osi_primal, osi_dual, osi_hot);
}
#endif

int main (int argc, const char *argv[])
{
    int returnCode = 0;
#ifdef CLP_DEBUG_MALLOC
    clp_memory(0);
#endif
    {
#ifndef CBC_OTHER_SOLVER
        OsiClpSolverInterface solver1;
#elif CBC_OTHER_SOLVER==1
        OsiCpxSolverInterface solver1;
#endif
        CbcModel model(solver1);

        // define TEST_MESSAGE_HANDLER at top of file to check works on all messages
#ifdef TEST_MESSAGE_HANDLER
        MyMessageHandler2 messageHandler(&model);
        std::cout << "Testing derived message handler" << std::endl;
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

        returnCode = CbcMain1 (argc, argv, model);
    }

#ifdef CLP_DEBUG_MALLOC
    clp_memory(1);
#endif

#ifdef CBC_STATISTICS
#endif

    if (returnCode != 777) {
        return returnCode;
    } else {
        return 0;
    }
}



#else	/* NEW_STYLE_SOLVER */

/*
  As best I can see, this is not yet fully functional. NEW_STYLE_SOLVER is
  normally undefined or forced to 0. See CbcSolver.hpp for jjf's thoughts on
  where this is going (a good direction, in my opinion [lh]).
*/

#include "CbcSolver.hpp"

/*
  ClpAmplStuff.cpp
*/
void addAmplToCbc(CbcSolver *);

int main (int argc, const char *argv[])
{
    int returnCode;

    /*
      // Only active if malloc switched on in CbcSolver.cpp

      Magic value 0 initialises the package. Anything else dumps statistics.
    */
# ifdef CLP_DEBUG_MALLOC
    clp_memory(0);
# endif

    // Open a block to ease memory leak checks.
    {
        OsiClpSolverInterface solver1;
        CbcSolver control(solver1);
        // initialize
        control.fillValuesInSolver();

#   ifdef COIN_HAS_ASL
        addAmplToCbc(&control);
#   endif

        returnCode = control.solve (argc, argv, 1);
    }
#ifdef CLP_DEBUG_MALLOC
    clp_memory(1);
#endif
    return returnCode;
}

#endif	/* NEW_STYLE_SOLVER */



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
