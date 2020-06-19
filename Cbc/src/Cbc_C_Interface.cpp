// $Id$
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cmath>
#include <cfloat>
#include <cctype>

#include "CoinPragma.hpp"
//#include "CoinHelperFunctions.hpp"
//#include "CoinPackedMatrix.hpp"
#include "CoinTime.hpp"

#include "CbcModel.hpp"
#include "CbcSolver.hpp"
#include "CbcBranchActual.hpp"

#include "CoinMessageHandler.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CglCutGenerator.hpp"

//  bobe including extras.h to get strdup()
#if defined(__MWERKS__)
// #include <extras.h>  // bobe 06-02-14
#endif

// Get C stuff but with extern C
#define CBC_EXTERN_C
#include "Coin_C_defines.h"

#define CbcGetProperty(T, prop) \
  COINLIBAPI T COINLINKAGE      \
    Cbc_##prop(Cbc_Model *m)    \
  {                             \
    return m->model_->prop();   \
  }

#define CbcSetSolverProperty(T, prop)          \
  COINLIBAPI void COINLINKAGE                  \
    Cbc_##prop(Cbc_Model *m, int index, T val) \
  {                                            \
    m->model_->solver()->prop(index, val);     \
  }

const int VERBOSE = 0;

// cut generator to accept callbacks in CBC
//
class CglCallback : public CglCutGenerator
{
    public:
        CglCallback();
        
        cbc_cut_callback cut_callback_;
        void *appdata;
        //CbcModel *model;

        /// Copy constructor
        CglCallback(const CglCallback& rhs);

        /// Clone
        virtual CglCutGenerator * clone() const;

        virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
                const CglTreeInfo info = CglTreeInfo() );

        virtual ~CglCallback();
    private:
};


CglCallback::CglCallback()
    : cut_callback_(NULL),
    appdata(NULL)
{
}

CglCallback::CglCallback(const CglCallback& rhs)
{
    this->cut_callback_ = rhs.cut_callback_;
    this->appdata = rhs.appdata;
}

CglCutGenerator* CglCallback::clone() const
{
    CglCallback *cglcb = new CglCallback();
    cglcb->cut_callback_ = this->cut_callback_;
    cglcb->appdata = this->appdata;

    return static_cast<CglCutGenerator*>(cglcb);
}

void CglCallback::generateCuts( const OsiSolverInterface &si, OsiCuts &cs, const CglTreeInfo info )
{
  this->cut_callback_( (OsiSolverInterface *) &si, &cs, this->appdata );
}

CglCallback::~CglCallback()
{

}

// To allow call backs
class Cbc_MessageHandler
  : public CoinMessageHandler {

public:
  /**@name Overrides */
  //@{
  virtual int print();
  //@}
  /**@name set and get */
  //@{
  /// Model
  const Cbc_Model *model() const;
  void setModel(Cbc_Model *model);
  /// Call back
  void setCallBack(cbc_callback callback);
  //@}

  /**@name Constructors, destructor */
  //@{
  /** Default constructor. */
  Cbc_MessageHandler();
  /// Constructor with pointer to model
  Cbc_MessageHandler(Cbc_Model *model,
    FILE *userPointer = NULL);
  /** Destructor */
  virtual ~Cbc_MessageHandler();
  //@}

  /**@name Copy method */
  //@{
  /** The copy constructor. */
  Cbc_MessageHandler(const Cbc_MessageHandler &);
  /** The copy constructor from an CoinSimplexMessageHandler. */
  Cbc_MessageHandler(const CoinMessageHandler &);

  Cbc_MessageHandler &operator=(const Cbc_MessageHandler &);
  /// Clone
  virtual CoinMessageHandler *clone() const;
  //@}

protected:
  /**@name Data members
       The data members are protected to allow access for derived classes. */
  //@{
  /// Pointer back to model
  Cbc_Model *model_;
  /// call back
  cbc_callback callback_;
  //@}
};

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
Cbc_MessageHandler::Cbc_MessageHandler()
  : CoinMessageHandler()
  , model_(NULL)
  , callback_(NULL)
{
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
Cbc_MessageHandler::Cbc_MessageHandler(const Cbc_MessageHandler &rhs)
  : CoinMessageHandler(rhs)
  , model_(rhs.model_)
  , callback_(rhs.callback_)
{
}

Cbc_MessageHandler::Cbc_MessageHandler(const CoinMessageHandler &rhs)
  : CoinMessageHandler(rhs)
  , model_(NULL)
  , callback_(NULL)
{
}

// Constructor with pointer to model
Cbc_MessageHandler::Cbc_MessageHandler(Cbc_Model *model,
  FILE * /*userPointer*/)
  : CoinMessageHandler()
  , model_(model)
  , callback_(NULL)
{
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
Cbc_MessageHandler::~Cbc_MessageHandler()
{
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
Cbc_MessageHandler &
Cbc_MessageHandler::operator=(const Cbc_MessageHandler &rhs)
{
  if (this != &rhs) {
    CoinMessageHandler::operator=(rhs);
    model_ = rhs.model_;
    callback_ = rhs.callback_;
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CoinMessageHandler *Cbc_MessageHandler::clone() const
{
  return new Cbc_MessageHandler(*this);
}
int Cbc_MessageHandler::print()
{
  if (callback_) {
    int messageNumber = currentMessage().externalNumber();
    if (currentSource() != "Cbc")
      messageNumber += 1000000;
    int i;
    int nDouble = numberDoubleFields();
    assert(nDouble <= 200);
    double vDouble[200];
    for (i = 0; i < nDouble; i++)
      vDouble[i] = doubleValue(i);
    int nInt = numberIntFields();
    assert(nInt <= 200);
    int vInt[200];
    for (i = 0; i < nInt; i++)
      vInt[i] = intValue(i);
    int nString = numberStringFields();
    assert(nString <= 200);
    char *vString[200];
    for (i = 0; i < nString; i++) {
      std::string value = stringValue(i);
      vString[i] = CoinStrdup(value.c_str());
    }
    callback_(model_, messageNumber,
      nDouble, vDouble,
      nInt, vInt,
      nString, vString);
    for (i = 0; i < nString; i++)
      free(vString[i]);
  }
  return CoinMessageHandler::print();
  return 0;
}
const Cbc_Model *
Cbc_MessageHandler::model() const
{
  return model_;
}
void Cbc_MessageHandler::setModel(Cbc_Model *model)
{
  model_ = model;
}
// Call back
void Cbc_MessageHandler::setCallBack(cbc_callback callback)
{
  callback_ = callback;
}
/**
  *
  *  C Interface Routines
  *
  */
#include "Cbc_C_Interface.h"
#include <string>
#include <stdio.h>
#include <iostream>

#if defined(__MWERKS__)
#pragma export on
#endif

/* Version */
COINLIBAPI const char *COINLINKAGE Cbc_getVersion()
{
  return CBC_VERSION;
}

// flushes buffers of new variables
static void Cbc_flush( Cbc_Model *model )
{
  if (model->nCols)
  {
    OsiSolverInterface *solver = model->model_->solver();

    int *starts = new int[model->nCols+1];
    for ( int i=0 ; (i<model->nCols+1) ; ++i )
      starts[i] = 0;

    int idx = 0; double coef = 0.0;

    int colsBefore = solver->getNumCols();
    
    solver->addCols( model->nCols, starts, &idx, &coef, model->cLB, model->cUB, model->cObj );

    for ( int i=0 ; i<model->nCols; ++i )
      if (model->cInt[i])
        solver->setInteger( colsBefore+i );

    for ( int i=0 ; i<model->nCols; ++i )
      solver->setColName( colsBefore+i, std::string(model->cNames+model->cNameStart[i]) );

    model->nCols = 0;

    delete[] starts;
  }
}

static void Cbc_checkSpaceColBuffer( Cbc_Model *model, int additionlNameSpace )
{
  // initialize buffer
  if ( model->colSpace == 0 )
  {
    // initial buffer allocation
    model->colSpace = 8192;
    int c = model->colSpace;
    model->nCols = 0;
    model->cNameSpace = 16384;

    model->cNameStart = (int *) malloc( sizeof(int)*c );
    assert( model->cNameStart );
    model->cNameStart[0] = 0;

    model->cInt = (char *) malloc( sizeof(char)*c );
    assert( model->cInt );

    model->cNames = (char *) malloc( sizeof(char)*model->cNameSpace );
    assert( model->cNames );

    model->cLB = (double *) malloc( sizeof(double)*c );
    assert( model->cLB );

    model->cUB = (double *)malloc( sizeof(double)*c );
    assert( model->cUB );

    model->cObj = (double *)malloc( sizeof(double)*c );
    assert( model->cObj );
  }
  else
  {
    // check buffer space
    if (model->nCols+2 >= model->colSpace)
    {
      model->colSpace *= 2;
      int c = model->colSpace;

      model->cNameStart = (int *) realloc( model->cNameStart, sizeof(int)*c );
      assert( model->cNameStart );

      model->cInt = (char *) realloc( model->cInt, sizeof(char)*c );
      assert( model->cInt );

      model->cLB = (double *) realloc( model->cLB, sizeof(double)*c );
      assert( model->cLB );

      model->cUB = (double *) realloc( model->cUB, sizeof(double)*c );
      assert( model->cUB );

      model->cObj = (double *) realloc( model->cObj, sizeof(double)*c );
      assert( model->cObj );
    }
    // check string buffer space
    int slen = additionlNameSpace + 1;
    int reqsize = slen + model->cNameStart[model->nCols]+1;
    if (reqsize>model->cNameSpace)
    {
      model->cNameSpace *= 2;
      model->cNames = (char *) realloc( model->cNames, sizeof(char)*model->cNameSpace );
    }
  }
}

static void Cbc_addColBuffer( Cbc_Model *model, 
    const char *name, double lb, double ub, double obj, 
    char isInteger )
{
  Cbc_checkSpaceColBuffer( model, 512 );
  int p = model->nCols;
  model->cInt[p] = isInteger;
  model->cLB[p] = lb;
  model->cUB[p] = ub;
  model->cObj[p] = obj;
  
  int ps = model->cNameStart[p];
  strcpy( model->cNames+ps, name );
  int len = strlen(name);

  model->nCols++;
  model->cNameStart[model->nCols] = ps + len + 1;
}

static void Cbc_deleteColBuffer( Cbc_Model *model )
{
  if ( model->colSpace > 0 )
  {
    free(model->cNameStart);
    free(model->cInt);
    free(model->cNames);
    free(model->cLB);
    free(model->cUB);
    free(model->cObj);
  }
}


/* Default Cbc_Model constructor */
COINLIBAPI Cbc_Model *COINLINKAGE
Cbc_newModel()
{
  const char prefix[] = "Cbc_C_Interface::Cbc_newModel(): ";
  //  const int  VERBOSE = 1;
  if (VERBOSE > 0)
    printf("%s begin\n", prefix);

  Cbc_Model *model = new Cbc_Model();
  OsiClpSolverInterface solver1; // will be release at the end of the scope, CbcModel clones it
  model->model_ = new CbcModel(solver1);
  model->solver_ = dynamic_cast< OsiClpSolverInterface * >(model->model_->solver());
  model->cbcData = new CbcSolverUsefulData();
  CbcMain0(*model->model_, *model->cbcData);
  model->handler_ = NULL;
  model->cbcData->noPrinting_ = false;
  model->relax_ = 0;

  // initialize columns buffer
  model->colSpace = 0;
  model->nCols = 0;
  model->cNameSpace = 0;
  model->cNameStart = NULL;
  model->cInt = NULL;
  model->cNames= NULL;
  model->cLB = NULL;
  model->cUB = NULL;
  model->cObj = NULL;

  if (VERBOSE > 0)
    printf("%s return\n", prefix);
  return model;
}
/* Cbc_Model Destructor */
COINLIBAPI void COINLINKAGE
Cbc_deleteModel(Cbc_Model *model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_deleteModel(): ";
  //  const int  VERBOSE = 1;
  if (VERBOSE > 0)
    printf("%s begin\n", prefix);
  fflush(stdout);

  Cbc_deleteColBuffer(model);

  if (VERBOSE > 1)
    printf("%s delete model->model_\n", prefix);
  fflush(stdout);
  delete model->model_;

  if (VERBOSE > 1)
    printf("%s delete model->handler_\n", prefix);
  fflush(stdout);
  delete model->handler_;

  delete model->cbcData;

  if (VERBOSE > 1)
    printf("%s delete model\n", prefix);
  fflush(stdout);
  delete model;

  if (VERBOSE > 0)
    printf("%s return\n", prefix);
  fflush(stdout);
}

/* Loads a problem (the constraints on the
    rows are given by lower and upper bounds). If a pointer is NULL then the
    following values are the default:
    <ul>
    <li> <code>colub</code>: all columns have upper bound infinity
    <li> <code>collb</code>: all columns have lower bound 0
    <li> <code>rowub</code>: all rows have upper bound infinity
    <li> <code>rowlb</code>: all rows have lower bound -infinity
    <li> <code>obj</code>: all variables have 0 objective coefficient
    </ul>

   Just like the other loadProblem() method except that the matrix is
   given in a standard column major ordered format (without gaps).
*/
COINLIBAPI void COINLINKAGE
Cbc_loadProblem(Cbc_Model *model, const int numcols, const int numrows,
  const CoinBigIndex *start, const int *index,
  const double *value,
  const double *collb, const double *colub,
  const double *obj,
  const double *rowlb, const double *rowub)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_loadProblem(): ";
  //  const int  VERBOSE = 2;
  if (VERBOSE > 0)
    printf("%s begin\n", prefix);

  OsiSolverInterface *solver = model->model_->solver();

  if (VERBOSE > 1) {
    printf("%s numcols = %i, numrows = %i\n",
      prefix, numcols, numrows);
    printf("%s model = %p, start = %p, index = %p, value = %p\n",
      prefix, static_cast< void * >(model), static_cast< const void * >(start),
      static_cast< const void * >(index), static_cast< const void * >(value));
    printf("%s collb = %p, colub = %p, obj = %p, rowlb = %p, rowub = %p\n",
      prefix, static_cast< const void * >(collb),
      static_cast< const void * >(colub), static_cast< const void * >(obj),
      static_cast< const void * >(rowlb), static_cast< const void * >(rowub));
  }

  if (VERBOSE > 1)
    printf("%s Calling solver->loadProblem()\n", prefix);
  fflush(stdout);

  solver->loadProblem(numcols, numrows, start, index, value,
    collb, colub, obj, rowlb, rowub);
  if (VERBOSE > 1)
    printf("%s Finished solver->loadProblem()\n", prefix);
  fflush(stdout);

  if (VERBOSE > 0)
    printf("%s return\n", prefix);
} //  Cbc_loadProblem()

/* Read an mps file from the given filename */
COINLIBAPI int COINLINKAGE
Cbc_readMps(Cbc_Model *model, const char *filename)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_readMps(): ";
  //  const int  VERBOSE = 2;
  if (VERBOSE > 0)
    printf("%s begin\n", prefix);
  if (VERBOSE > 1)
    printf("%s filename = '%s'\n", prefix, filename);

  int result = 1;
  result = model->model_->solver()->readMps(filename);
  assert(result == 0);

  if (VERBOSE > 0)
    printf("%s return %i\n", prefix, result);
  return result;
}
/* Write an mps file from the given filename */
COINLIBAPI void COINLINKAGE
Cbc_writeMps(Cbc_Model *model, const char *filename)
{
  Cbc_flush(model);

  const char prefix[] = "Cbc_C_Interface::Cbc_writeMps(): ";
  //  const int  VERBOSE = 2;
  if (VERBOSE > 0)
    printf("%s begin\n", prefix);
  if (VERBOSE > 1)
    printf("%s filename = '%s'\n", prefix, filename);

  model->model_->solver()->writeMps(filename, "mps", Cbc_getObjSense(model));

  if (VERBOSE > 0)
    printf("%s return\n", prefix);
  return;
}

/* Write an mps file from the given filename */
COINLIBAPI void COINLINKAGE
Cbc_writeLp(Cbc_Model *model, const char *filename)
{
  Cbc_flush(model);

  const char prefix[] = "Cbc_C_Interface::Cbc_writeLp(): ";
  //  const int  VERBOSE = 2;
  if (VERBOSE > 0)
    printf("%s begin\n", prefix);
  if (VERBOSE > 1)
    printf("%s filename = '%s'\n", prefix, filename);

  char outFile[512];
  strncpy(outFile, filename, 511);

  char *s = NULL;
  if ((s = strstr(outFile, ".lp"))) {
    if (s != outFile) // not at the start
      *s = '\0';
  }

  model->model_->solver()->writeLp(outFile);

  if (VERBOSE > 0)
    printf("%s return\n", prefix);
  return;
}

/* Read an lp file from the given filename */
COINLIBAPI int COINLINKAGE
Cbc_readLp(Cbc_Model *model, const char *filename)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_readLp(): ";
  //  const int  VERBOSE = 2;
  if (VERBOSE > 0)
    printf("%s begin\n", prefix);
  if (VERBOSE > 1)
    printf("%s filename = '%s'\n", prefix, filename);
  int result = 1;
  result = model->model_->solver()->readLp(filename);
  assert(result == 0);

  if (VERBOSE > 0)
    printf("%s return %i\n", prefix, result);
  return result;
}

COINLIBAPI void COINLINKAGE
Cbc_setInitialSolution(Cbc_Model *model, const double *sol)
{
  Cbc_flush(model); 
  int n = Cbc_getNumCols(model);
  // We need to manually compute the objective here for some reason
  const double *objvec = Cbc_getObjCoefficients(model);
  double objval = 0;
  for (int i = 0; i < n; i++) {
    objval += objvec[i] * sol[i];
  }
  model->model_->setBestSolution(sol, n, objval, true);
}

COINLIBAPI void COINLINKAGE
Cbc_setParameter(Cbc_Model *model, const char *name, const char *value)
{
  // checking if parameter is not included with another value
  // if this is the case just replacing this value
  std::string argname=std::string("-")+name;
  for ( int i=0 ; (i<((int)model->cmdargs_.size())-1) ; ++i )
  {
    if (argname==model->cmdargs_[i])
    {
      model->cmdargs_[i+1] = std::string(value);
      return;
    }
  }
  
  model->cmdargs_.push_back(argname);
  model->cmdargs_.push_back(value);
}

/* Fills in array with problem name  */
COINLIBAPI void COINLINKAGE
Cbc_problemName(Cbc_Model *model, int maxNumberCharacters, char *array)
{
  std::string name;
  model->model_->solver()->getStrParam(OsiProbName, name);
  strncpy(array, name.c_str(), maxNumberCharacters);
}
/* Sets problem name.  Must have \0 at end.  */
COINLIBAPI int COINLINKAGE
Cbc_setProblemName(Cbc_Model *model, const char *array)
{
  bool result = false;
  result = model->model_->solver()->setStrParam(OsiProbName, array);

  return (result) ? 1 : 0;
}

CbcGetProperty(int, status)

  CbcGetProperty(int, secondaryStatus)

  /* Number of elements in matrix */
  COINLIBAPI int COINLINKAGE
  Cbc_getNumElements(Cbc_Model *model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getNumElements(): ";
  //  const int  VERBOSE = 1;
  if (VERBOSE > 0)
    printf("%s begin\n", prefix);

  int result = 0;
  result = model->model_->getNumElements();

  if (VERBOSE > 0)
    printf("%s return %i\n", prefix, result);
  return result;
}

COINLIBAPI int COINLINKAGE
Cbc_getNumIntegers(Cbc_Model *model)
{
  Cbc_flush(model); 
  return model->model_->solver()->getNumIntegers();
}

// Column starts in matrix
COINLIBAPI const CoinBigIndex *COINLINKAGE
Cbc_getVectorStarts(Cbc_Model *model)
{
  const CoinPackedMatrix *matrix = NULL;
  matrix = model->model_->solver()->getMatrixByCol();
  return (matrix == NULL) ? NULL : matrix->getVectorStarts();
}
// Row indices in matrix
COINLIBAPI const int *COINLINKAGE
Cbc_getIndices(Cbc_Model *model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getIndices(): ";
  //  const int  VERBOSE = 1;
  if (VERBOSE > 0)
    printf("%s begin\n", prefix);

  const int *result = NULL;
  const CoinPackedMatrix *matrix = NULL;
  matrix = model->model_->solver()->getMatrixByCol();
  result = (matrix == NULL) ? NULL : matrix->getIndices();

  if (VERBOSE > 0)
    printf("%s return %p\n", prefix, static_cast< const void * >(result));
  return result;
}

// Element values in matrix
COINLIBAPI const double *COINLINKAGE
Cbc_getElements(Cbc_Model *model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_getElements(): ";
  //  const int  VERBOSE = 1;
  if (VERBOSE > 0)
    printf("%s begin\n", prefix);

  const double *result = NULL;
  const CoinPackedMatrix *matrix = NULL;
  matrix = model->model_->solver()->getMatrixByCol();
  result = (matrix == NULL) ? NULL : matrix->getElements();

  if (VERBOSE > 0)
    printf("%s return %p\n", prefix, static_cast< const void * >(result));
  return result;
}
// ======================================================================

/* Pass in Callback function */
COINLIBAPI void COINLINKAGE
Cbc_registerCallBack(Cbc_Model *model,
  cbc_callback userCallBack)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_registerCallBack(): ";
  //  const int  VERBOSE = 1;
  if (VERBOSE > 0)
    printf("%s begin\n", prefix);

  // Will be copy of users one
  delete model->handler_;
  model->handler_ = new Cbc_MessageHandler(*(model->model_->messageHandler()));
  model->handler_->setCallBack(userCallBack);
  model->handler_->setModel(model);
  model->model_->passInMessageHandler(model->handler_);

  if (VERBOSE > 0)
    printf("%s return\n", prefix);
}
/* Unset Callback function */
COINLIBAPI void COINLINKAGE
Cbc_clearCallBack(Cbc_Model *model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_clearCallBack(): ";
  //  const int  VERBOSE = 1;
  if (VERBOSE > 0)
    printf("%s begin\n", prefix);

  delete model->handler_;
  model->handler_ = NULL;

  if (VERBOSE > 0)
    printf("%s return\n", prefix);
}
/* length of names (0 means no names0 */
COINLIBAPI size_t COINLINKAGE
Cbc_maxNameLength(Cbc_Model *model)
{
  size_t result = 0;
  OsiSolverInterface::OsiNameVec const &rownames = model->model_->solver()->getRowNames();
  for (size_t i = 0; i < rownames.size(); i++) {
    if (rownames[i].length() > result)
      result = rownames[i].length();
  }
  OsiSolverInterface::OsiNameVec const &colnames = model->model_->solver()->getColNames();
  for (size_t i = 0; i < colnames.size(); i++) {
    if (colnames[i].length() > result)
      result = colnames[i].length();
  }
  return result;
}
COINLIBAPI void COINLINKAGE
Cbc_getRowName(Cbc_Model *model, int iRow, char *name, size_t maxLength)
{
  std::string rowname = model->model_->solver()->getRowName(iRow);
  strncpy(name, rowname.c_str(), maxLength);
  name[maxLength - 1] = '\0';
}

COINLIBAPI void COINLINKAGE
Cbc_getColName(Cbc_Model *model, int iColumn, char *name, size_t maxLength)
{
  assert( iColumn >= 0 );
  assert( iColumn < Cbc_getNumCols(model) );

  Cbc_flush(model);

  std::string colname = model->model_->solver()->getColName(iColumn);
  strncpy(name, colname.c_str(), maxLength);
  name[maxLength - 1] = '\0';
}

COINLIBAPI void COINLINKAGE
Cbc_setColName(Cbc_Model *model, int iColumn, const char *name)
{
  Cbc_flush(model);
  model->model_->solver()->setColName(iColumn, name);
}

COINLIBAPI void COINLINKAGE
Cbc_setRowName(Cbc_Model *model, int iRow, const char *name)
{
  model->model_->solver()->setRowName(iRow, name);
}

COINLIBAPI void COINLINKAGE
Cbc_setSolveRelax(Cbc_Model *model, char solveOnlyRelax)
{
  model->relax_ = solveOnlyRelax;
}

COINLIBAPI int COINLINKAGE
Cbc_solve(Cbc_Model *model)
{
  Cbc_flush( model );

  OsiSolverInterface *solver = model->solver_;
  if (solver->getNumIntegers() == 0 || model->relax_ == 1) {
    if (solver->basisIsAvailable()) {
      solver->resolve();
    } else {
      solver->initialSolve();
    }

    if (solver->isProvenOptimal())
      return 0;

    return 1;
  }

  const char prefix[] = "Cbc_C_Interface::Cbc_solve(): ";
  int result = 0;
  std::vector< const char * > argv;
  argv.push_back("Cbc_C_Interface");
  for (size_t i = 0; i < model->cmdargs_.size(); i++) {
    argv.push_back(model->cmdargs_[i].c_str());
  }
  argv.push_back("-solve");
  argv.push_back("-quit");
  try {

    CbcMain1((int)argv.size(), &argv[0], *model->model_, NULL, *model->cbcData);
  } catch (CoinError e) {
    printf("%s ERROR: %s::%s, %s\n", prefix,
      e.className().c_str(), e.methodName().c_str(), e.message().c_str());
  }
  result = model->model_->status();

  return result;
}

COINLIBAPI void COINLINKAGE Cbc_addCutCallback( 
    Cbc_Model *model, cbc_cut_callback cutcb, 
    const char *name, void *appData )
{
  assert( model != NULL );
  assert( model->model_ != NULL );

  CbcModel *cbcModel = model->model_;
  cbcModel->setKeepNamesPreproc(true);

  CglCallback cglCb;
  cglCb.appdata = appData;
  cglCb.cut_callback_ = cutcb;

  cbcModel->addCutGenerator( &cglCb, 1, name );
}

/* Sum of primal infeasibilities */
COINLIBAPI double COINLINKAGE
Cbc_sumPrimalInfeasibilities(Cbc_Model * /*model*/)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_sumPrimalInfeasibilities(): ";
  //  const int  VERBOSE = 1;
  if (VERBOSE > 0)
    printf("%s begin\n", prefix);

  double result = 0;
  // cannot find names in Cbc, Osi, or OsiClp
  //tbd result = model->model_->sumPrimalInfeasibilities();
  if (VERBOSE > 0)
    printf("%s WARNING:  NOT IMPLEMENTED\n", prefix);

  if (VERBOSE > 0)
    printf("%s return %g\n", prefix, result);
  return result;
}
/* Number of primal infeasibilities */
COINLIBAPI int COINLINKAGE
Cbc_numberPrimalInfeasibilities(Cbc_Model * /*model*/)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_numberPrimalInfeasibilities(): ";
  //  const int  VERBOSE = 1;
  if (VERBOSE > 0)
    printf("%s begin\n", prefix);

  int result = 0;
  //tbd  result = model->model_->getContinuousInfeasibilities();
  if (VERBOSE > 0)
    printf("%s WARNING:  NOT IMPLEMENTED\n", prefix);

  if (VERBOSE > 0)
    printf("%s return %i\n", prefix, result);
  return result;
}

/** Call this to really test if a valid solution can be feasible
    Solution is number columns in size.
    If fixVariables true then bounds of continuous solver updated.
    Returns objective value (worse than cutoff if not feasible)
*/
COINLIBAPI void COINLINKAGE
Cbc_checkSolution(Cbc_Model * /*model*/)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_checkSolution(): ";
  //  const int  VERBOSE = 1;
  if (VERBOSE > 0)
    printf("%s begin\n", prefix);

  // see CbcModel::checkSolution(double cutoff, const double * solution,
  //	       bool fixVariables);
  //  model->model_->checkSolution();

  if (VERBOSE > 0)
    printf("%s return\n", prefix);
  return;
}

COINLIBAPI int COINLINKAGE
Cbc_getNumCols(Cbc_Model *model)
{
  return model->model_->solver()->getNumCols() + model->nCols;
}

CbcGetProperty(int, getNumRows)
CbcGetProperty(int, getIterationCount)

  /** Number of non-zero entries in a row */
  COINLIBAPI int COINLINKAGE
  Cbc_getRowNz(Cbc_Model *model, int row)
{
  const CoinPackedMatrix *cpmRow = model->model_->solver()->getMatrixByRow();
  return cpmRow->getVectorLengths()[row];
}

/** Indices of variables that appear on this row */
COINLIBAPI const int *COINLINKAGE
Cbc_getRowIndices(Cbc_Model *model, int row)
{
  const CoinPackedMatrix *cpmRow = model->model_->solver()->getMatrixByRow();
  const CoinBigIndex *starts = cpmRow->getVectorStarts();
  const int *ridx = cpmRow->getIndices() + starts[row];
  return ridx;
}

/** Coefficients of variables that appear on this row */
COINLIBAPI const double *COINLINKAGE
Cbc_getRowCoeffs(Cbc_Model *model, int row)
{
  const CoinPackedMatrix *cpmRow = model->model_->solver()->getMatrixByRow();
  const CoinBigIndex *starts = cpmRow->getVectorStarts();
  const double *rcoef = cpmRow->getElements() + starts[row];
  return rcoef;
}

/** Number of non-zero entries in a column */
COINLIBAPI int COINLINKAGE
Cbc_getColNz(Cbc_Model *model, int col)
{
  const CoinPackedMatrix *cpmCol = model->model_->solver()->getMatrixByCol();
  return cpmCol->getVectorLengths()[col];
}

/** Indices of rows that a column appears */
COINLIBAPI const int *COINLINKAGE
Cbc_getColIndices(Cbc_Model *model, int col)
{
  const CoinPackedMatrix *cpmCol = model->model_->solver()->getMatrixByCol();
  const CoinBigIndex *starts = cpmCol->getVectorStarts();
  const int *cidx = cpmCol->getIndices() + starts[col];
  return cidx;
}

/** Coefficients that a column appear in rows */
COINLIBAPI const double *COINLINKAGE
Cbc_getColCoeffs(Cbc_Model *model, int col)
{
  const CoinPackedMatrix *cpmCol = model->model_->solver()->getMatrixByCol();
  const CoinBigIndex *starts = cpmCol->getVectorStarts();
  const double *rcoef = cpmCol->getElements() + starts[col];
  return rcoef;
}

/** Right hand side of a row */
COINLIBAPI double COINLINKAGE
Cbc_getRowRHS(Cbc_Model *model, int row)
{
  return model->model_->solver()->getRightHandSide()[row];
}

/** Sense a row */
COINLIBAPI char COINLINKAGE
Cbc_getRowSense(Cbc_Model *model, int row)
{
  return model->model_->solver()->getRowSense()[row];
}

/** Are there a numerical difficulties? */
COINLIBAPI int COINLINKAGE
Cbc_isAbandoned(Cbc_Model *model)
{
  if (Cbc_getNumIntegers(model) == 0 || model->relax_ == 1)
    return model->solver_->isAbandoned();
  else
    return model->model_->isAbandoned();
}

/** Is optimality proven? */
COINLIBAPI int COINLINKAGE
Cbc_isProvenOptimal(Cbc_Model *model)
{
  if (Cbc_getNumIntegers(model) == 0 || model->relax_ == 1)
    return model->solver_->isProvenOptimal();
  else
    return model->model_->isProvenOptimal();
}

COINLIBAPI int COINLINKAGE
Cbc_isProvenInfeasible(Cbc_Model *model)
{
  if (Cbc_getNumIntegers(model) == 0 || model->relax_ == 1)
    return (model->solver_->isProvenDualInfeasible() || model->solver_->isProvenPrimalInfeasible());
  else
    return model->model_->isProvenInfeasible();
}

COINLIBAPI double COINLINKAGE
Cbc_getObjValue(Cbc_Model *model)
{
  if (Cbc_getNumIntegers(model) == 0 || model->relax_ == 1)
    return (model->solver_->getObjValue());

  return model->model_->getObjValue();
}

COINLIBAPI const double *COINLINKAGE
Cbc_getReducedCost(Cbc_Model *model)
{
  if (Cbc_getNumIntegers(model) == 0 || model->relax_ == 1)
    return (model->solver_->getReducedCost());

  return model->model_->getReducedCost();
}

COINLIBAPI int COINLINKAGE
Cbc_numberSavedSolutions(Cbc_Model *model)
{
  return model->model_->numberSavedSolutions();
}

COINLIBAPI const double *COINLINKAGE
Cbc_savedSolution(Cbc_Model *model, int whichSol)
{
  return model->model_->savedSolution(whichSol);
}

COINLIBAPI double COINLINKAGE
Cbc_savedSolutionObj(Cbc_Model *model, int whichSol)
{
  return model->model_->savedSolutionObjective(whichSol);
}

COINLIBAPI const double *COINLINKAGE
Cbc_getColSolution(Cbc_Model *model)
{
  if (Cbc_getNumIntegers(model) == 0 || model->relax_ == 1)
    return (model->solver_->getColSolution());

  return model->model_->getColSolution();
}

CbcGetProperty(int, isContinuousUnbounded)
  CbcGetProperty(int, isNodeLimitReached)
    CbcGetProperty(int, isSecondsLimitReached)
      CbcGetProperty(int, isSolutionLimitReached)
        CbcGetProperty(int, isInitialSolveAbandoned)
          CbcGetProperty(int, isInitialSolveProvenOptimal)
            CbcGetProperty(int, isInitialSolveProvenPrimalInfeasible)

              CbcGetProperty(double, getObjSense)

                COINLIBAPI void COINLINKAGE
  Cbc_setObjSense(Cbc_Model *model, double sense)
{
  Cbc_flush(model);
  model->model_->setObjSense(sense);
}

CbcGetProperty(const double *, getRowActivity)

CbcGetProperty(const double *, getRowLower)
CbcSetSolverProperty(double, setRowLower)
CbcGetProperty(const double *, getRowUpper)
CbcSetSolverProperty(double, setRowUpper)

COINLIBAPI const double *COINLINKAGE
Cbc_getColLower(Cbc_Model *model)
{
  Cbc_flush(model); 
  return model->model_->solver()->getColLower();
}

COINLIBAPI const double *COINLINKAGE
Cbc_getColUpper(Cbc_Model *model)
{
  Cbc_flush(model); 
  return model->model_->solver()->getColUpper();
}

CbcGetProperty(double, getBestPossibleObjValue)

COINLIBAPI const double *COINLINKAGE
Cbc_getObjCoefficients(Cbc_Model *model)
{
  Cbc_flush(model); 
  return model->model_->solver()->getObjCoefficients();
}

COINLIBAPI void COINLINKAGE
Cbc_setObjCoeff(Cbc_Model *model, int index, double value)
{
  Cbc_flush( model );
  model->model_->solver()->setObjCoeff( index, value );
}

COINLIBAPI void COINLINKAGE
Cbc_setColLower(Cbc_Model *model, int index, double value)
{
  Cbc_flush(model);
  model->model_->solver()->setColLower( index, value );
}

COINLIBAPI void COINLINKAGE
Cbc_setColUpper(Cbc_Model *model, int index, double value)
{
  Cbc_flush(model);
  model->model_->solver()->setColUpper( index, value );
}


COINLIBAPI double *COINLINKAGE
Cbc_bestSolution(Cbc_Model *model)
{
  return model->model_->bestSolution();
}

/* Print model */
COINLIBAPI void COINLINKAGE
Cbc_printModel(Cbc_Model *model, const char *argPrefix)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_printModel(): ";
  const int VERBOSE = 4;
  if (VERBOSE > 0)
    printf("%s begin\n", prefix);

  CbcModel *cbc_model = model->model_;
  int numrows = cbc_model->getNumRows();
  int numcols = cbc_model->getNumCols();
  int numelem = cbc_model->getNumElements();
  const CoinPackedMatrix *matrix = cbc_model->solver()->getMatrixByCol();
  const CoinBigIndex *start = matrix->getVectorStarts();
  const int *index = matrix->getIndices();
  const double *value = matrix->getElements();
  const double *collb = cbc_model->getColLower();
  const double *colub = cbc_model->getColUpper();
  const double *obj = cbc_model->getObjCoefficients();
  const double *rowlb = cbc_model->getRowLower();
  const double *rowub = cbc_model->getRowUpper();

  printf("%s numcols = %i, numrows = %i, numelem = %i\n",
    argPrefix, numcols, numrows, numelem);
  printf("%s model = %p, start = %p, index = %p, value = %p\n",
    argPrefix, static_cast< void * >(model), static_cast< const void * >(start),
    static_cast< const void * >(index), static_cast< const void * >(value));
  matrix->dumpMatrix(NULL);
  {
    int i;
    for (i = 0; i <= numcols; i++)
      printf("%s start[%i] = %i\n", argPrefix, i, start[i]);
    for (i = 0; i < numelem; i++)
      printf("%s index[%i] = %i, value[%i] = %g\n",
        argPrefix, i, index[i], i, value[i]);
  }

  printf("%s collb = %p, colub = %p, obj = %p, rowlb = %p, rowub = %p\n",
    argPrefix, static_cast< const void * >(collb),
    static_cast< const void * >(colub), static_cast< const void * >(obj),
    static_cast< const void * >(rowlb), static_cast< const void * >(rowub));
  printf("%s optimization direction = %g\n", argPrefix, Cbc_getObjSense(model));
  printf("  (1 - minimize, -1 - maximize, 0 - ignore)\n");
  {
    int i;
    for (i = 0; i < numcols; i++)
      printf("%s collb[%i] = %g, colub[%i] = %g, obj[%i] = %g\n",
        argPrefix, i, collb[i], i, colub[i], i, obj[i]);
    for (i = 0; i < numrows; i++)
      printf("%s rowlb[%i] = %g, rowub[%i] = %g\n",
        argPrefix, i, rowlb[i], i, rowub[i]);
  }

  if (VERBOSE > 0)
    printf("%s return\n", prefix);
} // Cbc_printModel()

COINLIBAPI int COINLINKAGE
Cbc_isInteger(Cbc_Model *model, int i)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_isInteger(): ";
  //  const int  VERBOSE = 1;
  if (VERBOSE > 0)
    printf("%s begin\n", prefix);

  Cbc_flush(model); 

  bool result = false;
  result = model->model_->isInteger(i);

  if (VERBOSE > 0)
    printf("%s return %i\n", prefix, result);
  return (result) ? 1 : 0;
}

CbcGetProperty(int, getNodeCount)

/** Return a copy of this model */
COINLIBAPI Cbc_Model *COINLINKAGE
Cbc_clone(Cbc_Model *model)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_clone(): ";
  //  const int  VERBOSE = 1;
  if (VERBOSE > 0)
    printf("%s begin\n", prefix);

  Cbc_flush(model); 
  Cbc_Model *result = new Cbc_Model();
  result->model_ = new CbcModel(*(model->model_));
  result->solver_ = dynamic_cast< OsiClpSolverInterface * >(result->model_->solver());
  result->cbcData = new CbcSolverUsefulData();
  result->handler_ = NULL;
  result->cmdargs_ = model->cmdargs_;
  result->relax_ = model->relax_;
  result->cbcData->noPrinting_ = model->cbcData->noPrinting_;

  result->colSpace = 0;
  result->nCols = 0;
  result->cNameSpace = 0;
  result->cNameStart = NULL;
  result->cInt = NULL;
  result->cNames= NULL;
  result->cLB = NULL;
  result->cUB = NULL;
  result->cObj = NULL;

  if (VERBOSE > 0)
    printf("%s return\n", prefix);
  return result;
}
/** Set this the variable to be continuous */
COINLIBAPI void COINLINKAGE
Cbc_setContinuous(Cbc_Model *model, int iColumn)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setContinuous(): ";
  //  const int  VERBOSE = 1;
  if (VERBOSE > 0)
    printf("%s begin\n", prefix);

  Cbc_flush(model);

  model->model_->solver()->setContinuous(iColumn);

  if (VERBOSE > 0)
    printf("%s return\n", prefix);
}
/** Set this the variable to be integer */
COINLIBAPI void COINLINKAGE
Cbc_setInteger(Cbc_Model *model, int iColumn)
{
  const char prefix[] = "Cbc_C_Interface::Cbc_setContinuous(): ";
  //  const int  VERBOSE = 1;
  if (VERBOSE > 0)
    printf("%s begin\n", prefix);

  Cbc_flush(model);

  model->model_->solver()->setInteger(iColumn);

  if (VERBOSE > 0)
    printf("%s return\n", prefix);
}

/** Adds a new column */
COINLIBAPI void COINLINKAGE
Cbc_addCol(Cbc_Model *model, const char *name, double lb,
  double ub, double obj, char isInteger,
  int nz, int *rows, double *coefs)
{
  OsiSolverInterface *solver = model->model_->solver();

  if ( nz==0 )
  {
    Cbc_addColBuffer( model, name, lb, ub, obj, isInteger ); 
  }
  else
  {
    solver->addCol(nz, rows, coefs, lb, ub, obj, std::string(name));
    if (isInteger)
      solver->setInteger(solver->getNumCols() - 1);
  }
}

/** Adds a new row */
COINLIBAPI void COINLINKAGE
Cbc_addRow(Cbc_Model *model, const char *name, int nz,
  const int *cols, const double *coefs, char sense, double rhs)
{
  Cbc_flush(model);
  OsiSolverInterface *solver = model->model_->solver();
  double rowLB = -DBL_MAX, rowUB = DBL_MAX;
  switch (toupper(sense)) {
  case '=':
    rowLB = rowUB = rhs;
    break;
  case 'E':
    rowLB = rowUB = rhs;
    break;
  case '<':
    rowUB = rhs;
    break;
  case 'L':
    rowUB = rhs;
    break;
  case '>':
    rowLB = rhs;
    break;
  case 'G':
    rowLB = rhs;
    break;
  default:
    fprintf(stderr, "unknown row sense %c.", toupper(sense));
    abort();
  }
  solver->addRow(nz, cols, coefs, rowLB, rowUB);
  solver->setRowName(solver->getNumRows() - 1, std::string(name));
}

/** Add SOS constraints to the model using row-order matrix */

COINLIBAPI void COINLINKAGE
Cbc_addSOS(Cbc_Model *model, int numRows, const int *rowStarts,
  const int *colIndices, const double *weights, const int type)
{
  Cbc_flush(model);
  const char prefix[] = "Cbc_C_Interface::Cbc_addSOS(): ";
  //const int  VERBOSE = 4;
  if (VERBOSE > 0)
    printf("%sbegin\n", prefix);

  if (VERBOSE > 0)
    printf("%s numRows = %i\n", prefix, numRows);

  int row, i;
  const int *colIndex;
  const double *colWeight;

  // loop on rows and count number of objects according to numWeights>0
  int numObjects = 0;
  for (row = 0; row < numRows; row++) {
    if (VERBOSE > 2) {
      printf("%s row = %i\n", prefix, row);
      printf("%s rowStarts[%i] = %i\n", prefix, row, rowStarts[row]);
      printf("%s rowStarts[%i+1] = %i\n", prefix, row, rowStarts[row + 1]);
      fflush(stdout);
    }
    const int numWeights = rowStarts[row + 1] - rowStarts[row];
    if (VERBOSE > 2)
      printf("%s  numWeights = %i\n", prefix, numWeights);
    if (numWeights > 0)
      numObjects++;
  }

  // make objects
  CbcObject **objects = new CbcObject *[numObjects];
  //  if (VERBOSE>1) printf("%s numObjects = %i, objects = %X\n",prefix,numObjects,objects);

  // loop on rows and make an object when numWeights>0
  int objNum = 0;
  for (row = 0; row < numRows; row++) {
    if (VERBOSE > 2) {
      printf("%s row = %i\n", prefix, row);
      printf("%s rowStarts[%i] = %i\n", prefix, row, rowStarts[row]);
      printf("%s rowStarts[%i+1] = %i\n", prefix, row, rowStarts[row + 1]);
    }
    const int numWeights = rowStarts[row + 1] - rowStarts[row];
    if (VERBOSE > 2)
      printf("%s  numWeights = %i\n", prefix, numWeights);
    colIndex = colIndices + rowStarts[row];
    colWeight = weights + rowStarts[row];
    if (numWeights > 0) {
      // Make a CbcSOS and assign it to objects
      if (VERBOSE > 3) {
        for (i = 0; i < numWeights; i++) {
          printf("%s  colIndex [%i] = %i\n", prefix, i, colIndex[i]);
          printf("%s  colWeight[%i] = %f\n", prefix, i, colWeight[i]);
        }
        fflush(stdout);
      }
      objects[objNum] = new CbcSOS(model->model_, (int)(numWeights),
        (const int *)colIndex, (const double *)colWeight, (int)objNum, (int)type);
      //      if (VERBOSE>2) printf("%s objects[%i] = %X\n",prefix,objNum,objects[objNum]);
      if (objects[objNum] == NULL) {
        printf("%s ERROR: objects[%i] == NULL\n", prefix, objNum);
        fflush(stdout);
        assert(objects[objNum] != NULL);
      }
      objNum++;
    }
  }
  if (VERBOSE > 2) {
    printf("%s calling addObjects()\n", prefix);

    //    printf("%s numObjects = %i, objects = %X\n",prefix,numObjects,objects);
    //    for (row=0; row<numObjects; row++)
    //      printf("%s  objects[%i] = %X\n",prefix,row,objects[row]);
  }
  fflush(stdout);
  model->model_->addObjects(numObjects, objects);
  if (VERBOSE > 1)
    printf("%s finished addObjects()\n", prefix);

  for (objNum = 0; objNum < numObjects; objNum++)
    delete objects[objNum];
  delete[] objects;

  if (VERBOSE > 0)
    printf("%sreturn\n", prefix);
  return;
}

COINLIBAPI void COINLINKAGE
Cbc_setMIPStart(Cbc_Model *model, int count, const char **colNames, const double colValues[])
{
  Cbc_flush(model); 
  model->model_->setMIPStart(count, colNames, colValues);
}

COINLIBAPI void COINLINKAGE
Cbc_setMIPStartI(Cbc_Model *model, int count, const int colIdxs[], const double colValues[])
{
  Cbc_flush(model); 
  CbcModel *cbcModel = model->model_;
  OsiSolverInterface *solver = cbcModel->solver();

  int charSpace = count;
  for (int i = 0; (i < count); ++i)
    charSpace += solver->getColName(colIdxs[i]).size();

  char *allChars = new char[charSpace];
  char *s = allChars;
  char **names = new char *[count];
  for (int i = 0; (i < count); ++i) {
    names[i] = s;
    strcpy(s, solver->getColName(colIdxs[i]).c_str());
    s += solver->getColName(colIdxs[i]).size() + 1;
  }

  cbcModel->setMIPStart(count, (const char **)names, colValues);

  delete[] names;
  delete[] allChars;
}

/** Print the solution */
COINLIBAPI void COINLINKAGE
Cbc_printSolution(Cbc_Model *model)
{
  {
    //
    //  Now to print out row solution.  The methods used return const
    //  pointers - which is of course much more virtuous.
    //
    //  This version just does non-zero columns
    //

    // * Rows

    int numberRows = Cbc_getNumRows(model);
    int iRow;

    const double *rowPrimal = Cbc_getRowActivity(model);
    const double *rowLower = Cbc_getRowLower(model);
    const double *rowUpper = Cbc_getRowUpper(model);
    printf("--------------------------------------\n");

    // * If we have not kept names (parameter to readMps) this will be 0
    //    assert(Cbc_lengthNames(model));

    printf("                       Primal          Lower         Upper\n");
    for (iRow = 0; iRow < numberRows; iRow++) {
      double value;
      value = rowPrimal[iRow];
      if (value > 1.0e-8 || value < -1.0e-8) {
        char name[20];
        //      	Cbc_columnName(model,iColumn,name);
        sprintf(name, "ROW%5i", iRow);
        printf("%6d %8s", iRow, name);
        printf(" %13g", rowPrimal[iRow]);
        printf(" %13g", rowLower[iRow]);
        printf(" %13g", rowUpper[iRow]);
        printf("\n");
      }
    }
    printf("--------------------------------------\n");
  }
  {
    //
    //  Now to print out column solution.  The methods used return const
    //  pointers - which is of course much more virtuous.
    //
    //  This version just does non-zero columns
    //
    //

    // * Columns

    int numberColumns = Cbc_getNumCols(model);
    int iColumn;

    const double *columnPrimal = Cbc_getColSolution(model);
    const double *columnLower = Cbc_getColLower(model);
    const double *columnUpper = Cbc_getColUpper(model);
    const double *columnObjective = Cbc_getObjCoefficients(model);

    printf("--------------------------------------\n");

    // * If we have not kept names (parameter to readMps) this will be 0
    //    assert(Cbc_lengthNames(model));

    printf("                       Primal          Lower         Upper          Cost     isInteger\n");
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value;
      value = columnPrimal[iColumn];
      if (value > 1.0e-8 || value < -1.0e-8) {
        char name[20];
        //      	Cbc_columnName(model,iColumn,name);
        sprintf(name, "COL%5i", iColumn);
        printf("%6d %8s", iColumn, name);
        printf(" %13g", columnPrimal[iColumn]);
        printf(" %13g", columnLower[iColumn]);
        printf(" %13g", columnUpper[iColumn]);
        printf(" %13g", columnObjective[iColumn]);
        printf(" %13i", Cbc_isInteger(model, iColumn));
        printf("\n");
      }
    }
    printf("--------------------------------------\n");
  }
  if (0)
    Cbc_printModel(model, "cbc::main(): ");
  return;
}

COINLIBAPI int COINLINKAGE
Osi_getNumCols( void *osi )
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;
  return osiSolver->getNumCols();
}

/** @brief Returns column name in OsiSolverInterface object */
COINLIBAPI void COINLINKAGE
Osi_getColName( void *osi, int i, char *name, int maxLen )
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;
  strncpy( name, osiSolver->getColName(i).c_str(), maxLen );
}

/** @brief Returns column lower bounds in OsiSolverInterface object */
COINLIBAPI const double * COINLINKAGE
Osi_getColLower( void *osi )
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;
  return osiSolver->getColLower();
}

/** @brief Returns column upper bounds in OsiSolverInterface object */
COINLIBAPI const double * COINLINKAGE
Osi_getColUpper( void *osi )
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;
  return osiSolver->getColUpper();
}

/** @brief Returns integrality information for columns in OsiSolverInterface object */
COINLIBAPI int COINLINKAGE
Osi_isInteger( void *osi, int col )
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;
  return osiSolver->isInteger(col);
}

/** @brief Returns number of rows in OsiSolverInterface object */
COINLIBAPI int COINLINKAGE
Osi_getNumRows( void *osi )
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;
  return osiSolver->getNumRows();
}

COINLIBAPI int COINLINKAGE
Osi_getRowNz(void *osi, int row)
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;

  const CoinPackedMatrix *cpmRow = osiSolver->getMatrixByRow();
  return cpmRow->getVectorLengths()[row];
}

/** @brief Indices of variables that appear on a row */
COINLIBAPI const int *COINLINKAGE
Osi_getRowIndices(void *osi, int row)
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;

  const CoinPackedMatrix *cpmRow = osiSolver->getMatrixByRow();
  const CoinBigIndex *starts = cpmRow->getVectorStarts();
  const int *ridx = cpmRow->getIndices() + starts[row];
  return ridx;
}

/** @brief Coefficients of variables that appear on this row  */
COINLIBAPI const double *COINLINKAGE
Osi_getRowCoeffs(void *osi, int row)
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;

  const CoinPackedMatrix *cpmRow = osiSolver->getMatrixByRow();
  const CoinBigIndex *starts = cpmRow->getVectorStarts();
  const double *rcoef = cpmRow->getElements() + starts[row];
  return rcoef;
}

/** @brief Right hand side of a row  */
COINLIBAPI double COINLINKAGE
Osi_getRowRHS(void *osi, int row)
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;
  return osiSolver->getRightHandSide()[row];
}

/** @brief Sense a row 
     * @param model problem object 
     * @param row row index
     * @return row sense: E for =, L for <=, G for >= and R for ranged row
     **/
COINLIBAPI char COINLINKAGE
Osi_getRowSense(void *osi, int row)
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;

  return osiSolver->getRowSense()[row];
}

COINLIBAPI void COINLINKAGE
OsiCuts_addRowCut( void *osiCuts, int nz, const int *idx, const double *coef, char sense, double rhs )
{
  sense = toupper(sense);
  OsiCuts *oc = (OsiCuts *) osiCuts;

  OsiRowCut orc;
  orc.setRow( nz, idx, coef );


  orc.setLb(-DBL_MAX);
  orc.setUb(DBL_MAX);

  switch (toupper(sense)) {
  case '=':
    orc.setLb(rhs);
    orc.setUb(rhs);
    break;
  case 'E':
    orc.setLb(rhs);
    orc.setUb(rhs);
    break;
  case '<':
    orc.setUb(rhs);
    break;
  case 'L':
    orc.setUb(rhs);
    break;
  case '>':
    orc.setLb(rhs);
    break;
  case 'G':
    orc.setLb(rhs);
    break;
  default:
    fprintf(stderr, "unknown row sense %c.", toupper(sense));
    abort();
  }

  oc->insert(orc);
}



/** @brief Returns solution vector in OsiSolverInterface object */
COINLIBAPI const double * COINLINKAGE
Osi_getColSolution( void *osi )
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;

  return osiSolver->getColSolution();
}

COINLIBAPI double COINLINKAGE
Cbc_getAllowableGap(Cbc_Model* model)
{
  return model->model_->getAllowableGap();
}

COINLIBAPI void COINLINKAGE
Cbc_setAllowableGap(Cbc_Model* model, double allowedGap)
{
  model->model_->setAllowableGap(allowedGap);
}

COINLIBAPI double COINLINKAGE
Cbc_getAllowableFractionGap(Cbc_Model* model)
{
  return model->model_->getAllowableFractionGap();
}

COINLIBAPI void COINLINKAGE
Cbc_setAllowableFractionGap(Cbc_Model* model, double allowedFracionGap)
{
  model->model_->setAllowableFractionGap(allowedFracionGap);
}

/** returns the maximum number of nodes that can be explored in the search tree
 */
COINLIBAPI int COINLINKAGE
Cbc_getMaximumNodes(Cbc_Model *model)
{
  return model->model_->getMaximumNodes();
}

/** sets the maximum number of nodes that can be explored in the search tree
 */
COINLIBAPI void COINLINKAGE
Cbc_setMaximumNodes(Cbc_Model *model, int maxNodes)
{
  model->model_->setMaximumNodes(maxNodes);
}

/** returns solution limit for the search process
 */
COINLIBAPI int COINLINKAGE
Cbc_getMaximumSolutions(Cbc_Model *model)
{
  return model->model_->getMaximumSolutions();
}

/** sets a solution limit as a stopping criterion 
 */
COINLIBAPI void COINLINKAGE
Cbc_setMaximumSolutions(Cbc_Model *model, int maxSolutions)
{
  model->model_->setMaximumSolutions(maxSolutions);
}

/** returns the current log leven
 */
COINLIBAPI int COINLINKAGE
Cbc_getLogLevel(Cbc_Model *model)
{
  return model->model_->logLevel();
}

/** sets the log level
 */
COINLIBAPI void COINLINKAGE
Cbc_setLogLevel(Cbc_Model *model, int logLevel)
{
  model->model_->setLogLevel(logLevel);
}



COINLIBAPI double COINLINKAGE
Cbc_getCutoff(Cbc_Model* model)
{
  return model->model_->getCutoff();
}

COINLIBAPI void COINLINKAGE
Cbc_setCutoff(Cbc_Model* model, double cutoff)
{
  model->model_->setCutoff(cutoff);
}

COINLIBAPI double COINLINKAGE
Cbc_getAllowablePercentageGap(Cbc_Model* model)
{
  return model->model_->getAllowablePercentageGap();
}

COINLIBAPI void COINLINKAGE
Cbc_setAllowablePercentageGap(Cbc_Model* model,
    double allowedPercentageGap)
{
  model->model_->setAllowablePercentageGap(allowedPercentageGap);
}

COINLIBAPI double COINLINKAGE
Cbc_getMaximumSeconds(Cbc_Model *model)
{
  return model->model_->getMaximumSeconds();
}

COINLIBAPI void COINLINKAGE
Cbc_setMaximumSeconds(Cbc_Model *model, double maxSeconds)
{
  model->model_->setMaximumSeconds(maxSeconds);
}

#if defined(__MWERKS__)
#pragma export off
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
