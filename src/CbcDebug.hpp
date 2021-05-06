#ifndef CbcDebug_H
#define CbcDebug_H

#ifdef CLP_DEBUG_MALLOC

/*extern "C" */ void clp_memory(int type);
/*extern "C" */
void *clp_malloc(int length);
/*extern "C" */
void clp_free(void *array);
#include <malloc.h>
#include <exception>
#include <new>

#endif

// define TEST_MESSAGE_HANDLER to check works on all messages
// #define TEST_MESSAGE_HANDLER
#ifdef TEST_MESSAGE_HANDLER
// This driver shows how to trap messages - this is just as in unitTest.cpp
// ****** THis code is similar to MyMessageHandler.hpp and MyMessagehandler.cpp
#include "CoinMessageHandler.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"

/** This just adds a model to CoinMessage and a void pointer so
    user can trap messages and do useful stuff.
    This is used in Clp/Test/unitTest.cpp

    The file pointer is just there as an example of user stuff.

  -- lh 071026 -- An accurate summary. Nothing is actually happening here
  except that messages will be prefixed with "==", which serves the purpose
  of demonstrating that this message handler is active. The extra parameters
  (CbcModel, FILE) are unused.

*/

class MyMessageHandler2 : public CoinMessageHandler {

public:
  /**@name Overrides */
  //@{
  virtual int print();
  //@}
  /**@name set and get */
  //@{
  /// Model
  const CbcModel *model() const;
  void setModel(CbcModel *model);
  //@}

  /**@name Constructors, destructor */
  //@{
  /** Default constructor. */
  MyMessageHandler2();
  /// Constructor with pointer to model
  MyMessageHandler2(CbcModel *model,
    FILE *userPointer = NULL);
  /** Destructor */
  virtual ~MyMessageHandler2();
  //@}

  /**@name Copy method */
  //@{
  /** The copy constructor. */
  MyMessageHandler2(const MyMessageHandler2 &);
  /** The copy constructor from an CoinSimplexMessageHandler. */
  MyMessageHandler2(const CoinMessageHandler &);

  MyMessageHandler2 &operator=(const MyMessageHandler2 &);
  /// Clone
  virtual CoinMessageHandler *clone() const;
  //@}

protected:
  /**@name Data members
       The data members are protected to allow access for derived classes. */
  //@{
  /// Pointer back to model
  CbcModel *model_;
  //@}
};

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
MyMessageHandler2::MyMessageHandler2()
  : CoinMessageHandler()
  , model_(NULL)
{
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
MyMessageHandler2::MyMessageHandler2(const MyMessageHandler2 &rhs)
  : CoinMessageHandler(rhs)
  , model_(rhs.model_)
{
}

MyMessageHandler2::MyMessageHandler2(const CoinMessageHandler &rhs)
  : CoinMessageHandler()
  , model_(NULL)
{
}

// Constructor with pointer to model
MyMessageHandler2::MyMessageHandler2(CbcModel *model,
  FILE *userPointer)
  : CoinMessageHandler()
  , model_(model)
{
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
MyMessageHandler2::~MyMessageHandler2()
{
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
MyMessageHandler2 &
MyMessageHandler2::operator=(const MyMessageHandler2 &rhs)
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
CoinMessageHandler *MyMessageHandler2::clone() const
{
  return new MyMessageHandler2(*this);
}
int MyMessageHandler2::print()
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
void MyMessageHandler2::setModel(CbcModel *model)
{
  model_ = model;
}

void test_message_handler(CbcModel &model);

#endif /* TEST_MESSAGE_HANDLER */

#endif
