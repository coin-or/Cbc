// $Id$
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cmath>
#include <cfloat>
#include <cctype>
#include <map>
#include <string>
#ifdef CBC_THREAD
#include <pthread.h>
#endif

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
#include "CbcCutGenerator.hpp"

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

typedef std::map< std::string, int > NameIndex;

// cut generator to accept callbacks in CBC
//
class CglCallback : public CglCutGenerator
{
    public:
        CglCallback();

        cbc_cut_callback cut_callback_;
        void *appdata;
#ifdef CBC_THREAD
        pthread_mutex_t *cbcMutex;
#endif

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
#ifdef CBC_THREAD
    ,cbcMutex(NULL)
#endif

{
}

CglCallback::CglCallback(const CglCallback& rhs)
{
    this->cut_callback_ = rhs.cut_callback_;
    this->appdata = rhs.appdata;
#ifdef CBC_THREAD
    this->cbcMutex = rhs.cbcMutex;
#endif
}

CglCutGenerator* CglCallback::clone() const
{
    CglCallback *cglcb = new CglCallback();
    cglcb->cut_callback_ = this->cut_callback_;
    cglcb->appdata = this->appdata;
#ifdef CBC_THREAD
    cglcb->cbcMutex = this->cbcMutex;
#endif

    return static_cast<CglCutGenerator*>(cglcb);
}

void CglCallback::generateCuts( const OsiSolverInterface &si, OsiCuts &cs, const CglTreeInfo info )
{
#ifdef CBC_THREAD
    assert(this->cbcMutex);
    pthread_mutex_lock((this->cbcMutex));
#endif

  this->cut_callback_( (OsiSolverInterface *) &si, &cs, this->appdata );


#ifdef CBC_THREAD
    pthread_mutex_unlock((this->cbcMutex));
#endif
}

CglCallback::~CglCallback()
{

}

class Cbc_EventHandler : public CbcEventHandler {
public:
  virtual CbcAction event(CbcEvent whichEvent);
  Cbc_EventHandler();
  Cbc_EventHandler(CbcModel *model);
  virtual ~Cbc_EventHandler();
  Cbc_EventHandler(const Cbc_EventHandler &rhs);
  Cbc_EventHandler &operator=(const Cbc_EventHandler &rhs);
  virtual CbcEventHandler *clone() const;

  double bestCost;
  double bestBound;
  cbc_incumbent_callback inc_callback;
  cbc_progress_callback progr_callback;
  void *appData;
  void *pgAppData;
#ifdef CBC_THREAD
  pthread_mutex_t *cbcMutex;
#endif
};

Cbc_EventHandler::Cbc_EventHandler()
  : CbcEventHandler()
  , bestCost(COIN_DBL_MAX)
  , bestBound(COIN_DBL_MIN)
  , inc_callback(NULL)
  , progr_callback(NULL)
  , appData(NULL)
  , pgAppData(NULL)
#ifdef CBC_THREAD
  , cbcMutex(NULL)
#endif
{
}

Cbc_EventHandler::Cbc_EventHandler(const Cbc_EventHandler &rhs)
  : CbcEventHandler(rhs)
  , bestCost(rhs.bestCost)
  , bestBound(rhs.bestBound)
  , inc_callback(rhs.inc_callback)
  , progr_callback(rhs.progr_callback)
  , appData(rhs.appData)
  , pgAppData(rhs.pgAppData)
#ifdef CBC_THREAD
  , cbcMutex(rhs.cbcMutex)
#endif

{
}

Cbc_EventHandler::Cbc_EventHandler(CbcModel *model)
  : CbcEventHandler(model)
  , bestCost(COIN_DBL_MAX)
  , bestBound(COIN_DBL_MIN)
  , inc_callback(NULL)
  , progr_callback(NULL)
  , appData(NULL)
  , pgAppData(NULL)
#ifdef CBC_THREAD
  , cbcMutex(NULL)
#endif

{
}

Cbc_EventHandler::~Cbc_EventHandler()
{
}

Cbc_EventHandler &Cbc_EventHandler::operator=(const Cbc_EventHandler &rhs)
{
  if (this != &rhs) {
    CbcEventHandler::operator=(rhs);
    this->bestCost = rhs.bestCost;
    this->bestBound = rhs.bestBound;
    this->inc_callback = rhs.inc_callback;
    this->progr_callback = rhs.progr_callback;
    this->appData = rhs.appData;
    this->pgAppData = rhs.pgAppData;
#ifdef CBC_THREAD
    this->cbcMutex = rhs.cbcMutex;
#endif
  }
  return *this;
}

CbcEventHandler *Cbc_EventHandler::clone() const
{
  return new Cbc_EventHandler(*this);
}

CbcEventHandler::CbcAction Cbc_EventHandler::event(CbcEvent whichEvent)
{
  // if in sub tree carry on
  if ((model_->specialOptions() & 2048) == 0) {
    if ((whichEvent == solution || whichEvent == heuristicSolution)) {

      OsiSolverInterface *origSolver = model_->solver();
      const OsiSolverInterface *pps = model_->postProcessedSolver(1);

      const OsiSolverInterface *solver = pps ? pps : origSolver;

#ifdef CBC_THREAD
    assert(this->cbcMutex);
    pthread_mutex_lock((this->cbcMutex));
#endif
      if (bestCost >= solver->getObjValue() + 1e-10) {
          bestCost = solver->getObjValue();

        if (this->inc_callback != NULL) {
          int charSize = 0, nNZ = 0;
          const double *x = solver->getColSolution();
          for (int i = 0; (i < solver->getNumCols()); ++i) {
            if (fabs(x[i]) <= 1e-7)
                continue;
            charSize += solver->getColName(i).size()+1;
            ++nNZ;
          } // checking non zero cols
          char **cnames = new char*[nNZ];
          double *xv = new double[nNZ];
          cnames[0] = new char[charSize];

          int cnz = 0;
          for (int i = 0; (i < solver->getNumCols()); ++i)
          {
            if (fabs(x[i]) <= 1e-7)
              continue;
            if (cnz>=1)
              cnames[cnz] = cnames[cnz-1] + solver->getColName(i).size()+1;
            cnz++;
          }

          cnz = 0;
          for (int i = 0; (i < solver->getNumCols()); ++i)
          {
            if (fabs(x[i]) <= 1e-7)
                continue;
            strcpy(cnames[cnz], solver->getColName(i).c_str());
            xv[cnz] = x[i];

            cnz++;
          }

          this->inc_callback(model_, bestCost, nNZ, cnames, xv, this->appData);

          delete[] xv;
          delete[] cnames[0];
          delete[] cnames;
        } // incumbent callback
        if (this->progr_callback) {
          const CbcModel *model = this->getModel();
          double secs = model->getCurrentSeconds();
          double lb = model->getBestPossibleObjValue();
          this->progr_callback((void *) model, 0, model->getNodeCount(), "sol", secs, lb, bestCost, 0, NULL, this->pgAppData);
        } // incumbent callback
      } // improved cost
#ifdef CBC_THREAD
    assert(this->cbcMutex);
    pthread_mutex_unlock((this->cbcMutex));
#endif
    } // solution found

    if (this->progr_callback)
    {
#ifdef CBC_THREAD
    assert(this->cbcMutex);
    pthread_mutex_lock((this->cbcMutex));
#endif
      if (this->getModel()->getBestPossibleObjValue() >= this->bestBound+1e-5)
      {
        this->bestBound = this->getModel()->getBestPossibleObjValue();
        const CbcModel *model = this->getModel();
        double secs = model->getCurrentSeconds();
        this->progr_callback((void *) model, 0, model->getNodeCount(), "bnd", secs, bestBound, bestCost, 0, NULL, this->pgAppData);
      } // improved lower bound
#ifdef CBC_THREAD
    assert(this->cbcMutex);
    pthread_mutex_unlock((this->cbcMutex));
#endif
    } // has progress callback
  } // not in subtree

  return noAction;
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

enum FlushContents
{
  FCColumns,
  FCRows,
  FCBoth
};

// flushes buffers of new variables
static void Cbc_flush( Cbc_Model *model, enum FlushContents fc = FCBoth )
{
  OsiSolverInterface *solver = model->model_->solver();

  if (fc == FCBoth || fc == FCColumns)
  {
    if (model->nCols)
    {
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
  if (fc == FCRows || fc == FCBoth)
  {
    if (model->nRows)
    {
      int rowsBefore = solver->getNumRows();
      solver->addRows(model->nRows, model->rStart, model->rIdx, model->rCoef, model->rLB, model->rUB);

      for ( int i=0 ; i<model->nRows; ++i )
      {
        const int rIdx = rowsBefore+i;
        const std::string rName = std::string(model->rNames+model->rNameStart[i]);
        solver->setRowName(rIdx, rName);
        if (model->rowNameIndex)
        {
          NameIndex &rowNameIndex = *((NameIndex  *)model->rowNameIndex);
          rowNameIndex[rName] = rIdx;
        }
      }

      model->nRows = 0;
    }
  }
} // flush cols, rows or both

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

static void Cbc_checkSpaceRowBuffer(Cbc_Model *model, int nzRow, int rowNameLen)
{
  if (model->rowSpace == 0)
  {
    // allocating buffer
    model->rowSpace = 8192;
    model->rStart = (int *)malloc(sizeof(int)*model->rowSpace);
    model->rStart[0] = 0;
    model->rLB = (double *)malloc(sizeof(double)*model->rowSpace);
    model->rUB = (double *)malloc(sizeof(double)*model->rowSpace);
    model->rNameStart = (int *)malloc(sizeof(int)*model->rowSpace);
    model->rNameStart[0] = 0;

    model->rElementsSpace = std::max(131072, nzRow * 2);
    model->rIdx = (int *)malloc(sizeof(int)*model->rElementsSpace);
    model->rCoef = (double *)malloc(sizeof(double)*model->rElementsSpace);

    model->rNameSpace = 131072;
    model->rNames = (char *)malloc(sizeof(char)*model->rNameSpace);
  }
  else
  {
    // checking if some resize is needed
    if (model->nRows+2 >= model->rowSpace)
    {
      // checking row space
      if (model->rowSpace < 1048576)
      {
        model->rowSpace *= 2;
        model->rStart = (int *)realloc(model->rStart, sizeof(int)*model->rowSpace);
        model->rLB = (double *)realloc(model->rLB, sizeof(double)*model->rowSpace);
        model->rUB = (double *)realloc(model->rUB, sizeof(double)*model->rowSpace);
        model->rNameStart = (int *)realloc(model->rNameStart, sizeof(int)*model->rowSpace);
      }
      else
      {
        Cbc_flush(model, FCRows);
      }
    } // rowSpace

    if (model->rStart[model->nRows]+nzRow+1 >= model->rElementsSpace)
    {
      if (model->rElementsSpace<4194304 || nzRow>=4194304)
      {
        model->rElementsSpace *= 2;
        model->rElementsSpace = std::max(model->rElementsSpace, nzRow*2);
        model->rIdx = (int *)realloc(model->rIdx, sizeof(int)*model->rElementsSpace);
        model->rCoef = (double *)realloc(model->rCoef, sizeof(double)*model->rElementsSpace);
      }
      else
      {
        Cbc_flush(model, FCRows);
      }
    } // elements space

    if (model->rNameStart[model->nRows]+rowNameLen+2 >= model->rNameSpace)
    {
      if (model->rNameSpace < 8388608)
      {
        model->rNameSpace *= 2;
        model->rNames = (char *)realloc(model->rNames, sizeof(char)*model->rNameSpace);
      }
      else
      {
        Cbc_flush(model, FCRows);
      }
    } // row names resize
  } // resizing
} // Cbc_checkSpaceRowBuffer

static void Cbc_addRowBuffer(
    Cbc_Model *model,
    int nz,
    const int *rIdx,
    const double *rCoef,
    double rLB,
    double rUB,
    const char *rName)
{
  int nameLen = strlen(rName);
  Cbc_checkSpaceRowBuffer(model, nz, nameLen);
  const int st = model->rStart[model->nRows];

  model->rLB[model->nRows] = rLB;
  model->rUB[model->nRows] = rUB;
  memcpy(model->rIdx + st, rIdx, sizeof(int)*nz);
  memcpy(model->rCoef + st, rCoef, sizeof(double)*nz);

  char *spcName = model->rNames + model->rNameStart[model->nRows];
  strcpy(spcName, rName);

  model->nRows++;
  model->rStart[model->nRows] = st + nz;
  model->rNameStart[model->nRows] = model->rNameStart[model->nRows-1] + nameLen + 1;
}

static void Cbc_deleteRowBuffer(Cbc_Model *model)
{
  if (model->rowSpace)
  {
    free(model->rStart);
    free(model->rLB);
    free(model->rUB);
    free(model->rNameStart);
    free(model->rIdx);
    free(model->rCoef);
    free(model->rNames);
  }
}

static void Cbc_iniBuffer(Cbc_Model *model) 
{
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

  // initialize rows buffer
  model->rowSpace = 0;
  model->nRows = 0;
  model->rNameSpace = 0;
  model->rNameStart = 0;
  model->rNames = NULL;
  model->rLB = NULL;
  model->rUB = NULL;
  model->rElementsSpace = 0;
  model->rStart = NULL;
  model->rIdx = NULL;
  model->rCoef = NULL;
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
  model->inc_callback = NULL;
  model->progr_callback = NULL;
  model->icAppData = NULL;
  model->pgrAppData = NULL;
  model->colNameIndex = NULL;
  model->rowNameIndex = NULL;

  Cbc_iniBuffer(model);

#ifdef CBC_THREAD
  pthread_mutex_init(&(model->cbcMutex), NULL);
#endif

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
  Cbc_deleteRowBuffer(model);

  #ifdef CBC_THREAD
  pthread_mutex_destroy(&(model->cbcMutex));
#endif


  if (model->colNameIndex)
  {
    NameIndex *m = (NameIndex *)model->colNameIndex;
    delete m;
    m = (NameIndex *)model->rowNameIndex;
    assert( m != NULL );
    delete m;
  }

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

/* should be called after reading a new problem */
static void fillAllNameIndexes(Cbc_Model *model)
{
  if (!model->colNameIndex)
    return;

  OsiSolverInterface *solver = model->model_->solver();
  NameIndex &colNameIndex = *((NameIndex  *)model->colNameIndex);
  colNameIndex.clear();
  NameIndex &rowNameIndex = *((NameIndex  *)model->rowNameIndex);
  rowNameIndex.clear();

  for ( int i=0 ; (i<solver->getNumCols()) ; ++i )
    colNameIndex[solver->getColName(i)] = i;

  for ( int i=0 ; (i<solver->getNumRows()) ; ++i )
    rowNameIndex[solver->getRowName(i)] = i;
}

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
  OsiSolverInterface *solver = model->model_->solver();
  result = solver->readMps(filename);
  assert(result == 0);

  Cbc_deleteColBuffer(model);
  Cbc_deleteRowBuffer(model);
  Cbc_iniBuffer(model);

  fillAllNameIndexes(model);

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
  OsiSolverInterface *solver = model->model_->solver();
  result = solver->readLp(filename);
  assert(result == 0);

  Cbc_deleteColBuffer(model);
  Cbc_deleteRowBuffer(model);
  Cbc_iniBuffer(model);

  fillAllNameIndexes(model);

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
  Cbc_flush(model);

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
  Cbc_flush(model, FCColumns);
  return model->model_->solver()->getNumIntegers();
}

// Column starts in matrix
COINLIBAPI const CoinBigIndex *COINLINKAGE
Cbc_getVectorStarts(Cbc_Model *model)
{
  Cbc_flush(model);
  const CoinPackedMatrix *matrix = NULL;
  matrix = model->model_->solver()->getMatrixByCol();
  return (matrix == NULL) ? NULL : matrix->getVectorStarts();
}
// Row indices in matrix
COINLIBAPI const int *COINLINKAGE
Cbc_getIndices(Cbc_Model *model)
{
  Cbc_flush(model, FCRows);
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
  Cbc_flush(model, FCRows);
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
  Cbc_flush(model, FCRows);
  std::string rowname = model->model_->solver()->getRowName(iRow);
  strncpy(name, rowname.c_str(), maxLength);
  name[maxLength - 1] = '\0';
}

COINLIBAPI void COINLINKAGE
Cbc_getColName(Cbc_Model *model, int iColumn, char *name, size_t maxLength)
{
  assert( iColumn >= 0 );
  assert( iColumn < Cbc_getNumCols(model) );

  Cbc_flush(model, FCColumns);

  std::string colname = model->model_->solver()->getColName(iColumn);
  strncpy(name, colname.c_str(), maxLength);
  name[maxLength - 1] = '\0';
}

COINLIBAPI void COINLINKAGE
Cbc_setColName(Cbc_Model *model, int iColumn, const char *name)
{
  Cbc_flush(model);
  OsiSolverInterface *solver = model->model_->solver();
  std::string previousName = solver->getColName(iColumn);
  solver->setColName(iColumn, name);

  if (!model->colNameIndex)
    return;
  NameIndex &colNameIndex = *((NameIndex  *)model->colNameIndex);
  NameIndex::iterator it = colNameIndex.find(previousName);
  if (it!=colNameIndex.end())
    colNameIndex.erase(it);
  colNameIndex[name] = iColumn;
}

COINLIBAPI void COINLINKAGE
Cbc_setRowName(Cbc_Model *model, int iRow, const char *name)
{
  Cbc_flush(model, FCRows);
  OsiSolverInterface *solver = model->model_->solver();
  std::string previousName = solver->getRowName(iRow);
  solver->setRowName(iRow, name);

  if (!model->rowNameIndex)
    return;
  NameIndex &rowNameIndex = *((NameIndex  *)model->rowNameIndex);
  NameIndex::iterator it = rowNameIndex.find(std::string(previousName));
  if (it!=rowNameIndex.end())
    rowNameIndex.erase(it);
  rowNameIndex[name] = iRow;
}

COINLIBAPI void COINLINKAGE
Cbc_setSolveRelax(Cbc_Model *model, char solveOnlyRelax)
{
  model->relax_ = solveOnlyRelax;
}

static int cbc_callb(CbcModel *cbcModel, int whereFrom) {
  Cbc_EventHandler *cbc_eh = dynamic_cast<Cbc_EventHandler *>(cbcModel->getEventHandler());
  if (cbc_eh == NULL)
    return 0;

  if (cbc_eh->progr_callback != NULL)
  {
    switch (whereFrom)
    {
      case 1:
        cbc_eh->progr_callback(cbcModel, 0, 0, "root lp relax", cbcModel->getCurrentSeconds(),
                               cbcModel->getBestPossibleObjValue(), DBL_MAX, 0, NULL, cbc_eh->pgAppData);
        break;
      case 2:
        cbc_eh->progr_callback(cbcModel, 0, 0, "pre-proc", cbcModel->getCurrentSeconds(),
                               cbcModel->getBestPossibleObjValue(), cbcModel->getObjValue(), 0, NULL, cbc_eh->pgAppData);
        break;
    } // where from
  } // has progress callback

  return 0;
}

COINLIBAPI int COINLINKAGE
Cbc_solve(Cbc_Model *model)
{
  Cbc_flush( model );

  OsiSolverInterface *solver = model->solver_;
  if (solver->getNumIntegers() == 0 || model->relax_ == 1) {

    solver->messageHandler()->setLogLevel( model->model_->logLevel() );

#ifdef COIN_HAS_CLP
    double maxTime = Cbc_getMaximumSeconds(model);
    if (maxTime != DBL_MAX) {
      OsiClpSolverInterface *clpSolver
        = dynamic_cast< OsiClpSolverInterface * >(solver);
      if (clpSolver) {
        bool useElapsedTime = model->model_->useElapsedTime();
        if (useElapsedTime)
          clpSolver->getModelPtr()->setMaximumWallSeconds(maxTime);
        else
          clpSolver->getModelPtr()->setMaximumSeconds(maxTime);
      } // clp solver
    } // set time limit
#endif

    if (solver->basisIsAvailable()) {
      solver->resolve();
    } else {
      solver->initialSolve();
    } // initial solve

    if (solver->isProvenOptimal())
      return 0;

    return 1;
  } // solve only lp relaxation

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

    Cbc_EventHandler *cbc_eh = NULL;
    if (model->inc_callback!=NULL || model->progr_callback!=NULL)
    {
      cbc_eh = new Cbc_EventHandler(model->model_);
#ifdef CBC_THREAD
      cbc_eh->cbcMutex = &(model->cbcMutex);
#endif

      if (model->inc_callback) {
        cbc_eh->inc_callback = model->inc_callback;
        cbc_eh->appData = model->icAppData;
      }
      if (model->progr_callback) {
        cbc_eh->progr_callback = model->progr_callback;
        cbc_eh->pgAppData = model->pgrAppData;
      }

      model->model_->passInEventHandler(cbc_eh);
    }

    CbcMain1((int)argv.size(), &argv[0], *model->model_, cbc_callb, *model->cbcData);

    if (cbc_eh)
    {
      delete cbc_eh;
      cbc_eh = NULL;
    }
  } catch (CoinError e) {
    printf("%s ERROR: %s::%s, %s\n", prefix,
      e.className().c_str(), e.methodName().c_str(), e.message().c_str());
  }
  result = model->model_->status();


  return result;
}

COINLIBAPI void COINLINKAGE Cbc_addIncCallback(
    Cbc_Model *model, cbc_incumbent_callback inccb,
    void *appData )
{
  model->inc_callback = inccb;
  model->icAppData = appData;
}

COINLIBAPI void COINLINKAGE Cbc_addProgrCallback(
  Cbc_Model *model, cbc_progress_callback prgcbc,
  void *appData)
{
  model->progr_callback = prgcbc;
  model->pgrAppData = appData;
}

COINLIBAPI void COINLINKAGE Cbc_addCutCallback(
    Cbc_Model *model, cbc_cut_callback cutcb,
    const char *name, void *appData )
{
  bool deleteCb = false;
  assert( model != NULL );
  assert( model->model_ != NULL );

  CbcModel *cbcModel = model->model_;
  cbcModel->setKeepNamesPreproc(true);

  CglCallback *cglCb = NULL;
  for ( int i=0 ; (i<cbcModel->numberCutGenerators()) ; ++i )
  {
    CbcCutGenerator *ccb = cbcModel->cutGenerators()[i];
    CglCallback *t = dynamic_cast<CglCallback *>(ccb->generator());
    if (t) {
      cglCb = t;
      break;
    }
  }
  if (cglCb==NULL)
  {
    cglCb = new CglCallback();
    deleteCb = true;
  }

  cglCb->appdata = appData;
  cglCb->cut_callback_ = cutcb;
#ifdef CBC_THREAD
  cglCb->cbcMutex = &(model->cbcMutex);
#endif

  cbcModel->addCutGenerator( cglCb, 1, name );

  if (deleteCb)
    delete cglCb;
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

COINLIBAPI int COINLINKAGE
Cbc_getNumRows(Cbc_Model *model)
{
  return model->model_->solver()->getNumRows() + model->nRows;
}


CbcGetProperty(int, getIterationCount)

/** Number of non-zero entries in a row */
COINLIBAPI int COINLINKAGE
Cbc_getRowNz(Cbc_Model *model, int row)
{
  Cbc_flush(model, FCRows);
  const CoinPackedMatrix *cpmRow = model->model_->solver()->getMatrixByRow();
  return cpmRow->getVectorLengths()[row];
}

/** Indices of variables that appear on this row */
COINLIBAPI const int *COINLINKAGE
Cbc_getRowIndices(Cbc_Model *model, int row)
{
  Cbc_flush(model, FCRows);
  const CoinPackedMatrix *cpmRow = model->model_->solver()->getMatrixByRow();
  const CoinBigIndex *starts = cpmRow->getVectorStarts();
  const int *ridx = cpmRow->getIndices() + starts[row];
  return ridx;
}

/** Coefficients of variables that appear on this row */
COINLIBAPI const double *COINLINKAGE
Cbc_getRowCoeffs(Cbc_Model *model, int row)
{
  Cbc_flush(model, FCRows);
  const CoinPackedMatrix *cpmRow = model->model_->solver()->getMatrixByRow();
  const CoinBigIndex *starts = cpmRow->getVectorStarts();
  const double *rcoef = cpmRow->getElements() + starts[row];
  return rcoef;
}

/** Number of non-zero entries in a column */
COINLIBAPI int COINLINKAGE
Cbc_getColNz(Cbc_Model *model, int col)
{
  Cbc_flush(model);
  const CoinPackedMatrix *cpmCol = model->model_->solver()->getMatrixByCol();
  return cpmCol->getVectorLengths()[col];
}

/** Indices of rows that a column appears */
COINLIBAPI const int *COINLINKAGE
Cbc_getColIndices(Cbc_Model *model, int col)
{
  Cbc_flush(model);
  const CoinPackedMatrix *cpmCol = model->model_->solver()->getMatrixByCol();
  const CoinBigIndex *starts = cpmCol->getVectorStarts();
  const int *cidx = cpmCol->getIndices() + starts[col];
  return cidx;
}

/** Coefficients that a column appear in rows */
COINLIBAPI const double *COINLINKAGE
Cbc_getColCoeffs(Cbc_Model *model, int col)
{
  Cbc_flush(model);
  const CoinPackedMatrix *cpmCol = model->model_->solver()->getMatrixByCol();
  const CoinBigIndex *starts = cpmCol->getVectorStarts();
  const double *rcoef = cpmCol->getElements() + starts[col];
  return rcoef;
}

/** Right hand side of a row */
COINLIBAPI double COINLINKAGE
Cbc_getRowRHS(Cbc_Model *model, int row)
{
  Cbc_flush(model, FCRows);
  return model->model_->solver()->getRightHandSide()[row];
}

/** Sense a row */
COINLIBAPI char COINLINKAGE
Cbc_getRowSense(Cbc_Model *model, int row)
{
  Cbc_flush(model, FCRows);
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

COINLIBAPI const double *COINLINKAGE
Cbc_getRowPrice(Cbc_Model *model)
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
  Cbc_flush(model, FCColumns);
  return model->model_->solver()->getColLower();
}

COINLIBAPI const double *COINLINKAGE
Cbc_getColUpper(Cbc_Model *model)
{
  Cbc_flush(model, FCColumns);
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
  Cbc_flush(model, FCColumns);
  model->model_->solver()->setColLower( index, value );
}

COINLIBAPI void COINLINKAGE
Cbc_setColUpper(Cbc_Model *model, int index, double value)
{
  Cbc_flush(model, FCColumns);
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
  Cbc_flush(model);
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

  Cbc_flush(model, FCColumns);

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
  result->inc_callback = model->inc_callback;
  result->progr_callback = model->progr_callback;
  result->icAppData = model->icAppData;
  result->pgrAppData = model->pgrAppData;
  if (model->colNameIndex)
    Cbc_storeNameIndexes(result, 1);

  result->colSpace = 0;
  result->nCols = 0;
  result->cNameSpace = 0;
  result->cNameStart = NULL;
  result->cInt = NULL;
  result->cNames= NULL;
  result->cLB = NULL;
  result->cUB = NULL;
  result->cObj = NULL;

  model->rowSpace = 0;
  model->nRows = 0;
  model->rNameSpace = 0;
  model->rNameStart = 0;
  model->rNames = NULL;
  model->rLB = NULL;
  model->rUB = NULL;
  model->rElementsSpace = 0;
  model->rStart = NULL;
  model->rIdx = NULL;
  model->rCoef = NULL;

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

  Cbc_flush(model, FCColumns);

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

  Cbc_flush(model, FCColumns);

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

  if (model->colNameIndex)
  {
    NameIndex &colNameIndex = *((NameIndex  *)model->colNameIndex);
    colNameIndex[std::string(name)] = Cbc_getNumCols(model)-1;
  }
}

/** Adds a new row */
COINLIBAPI void COINLINKAGE
Cbc_addRow(Cbc_Model *model, const char *name, int nz,
  const int *cols, const double *coefs, char sense, double rhs)
{
  Cbc_flush(model, FCColumns);
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
    fprintf(stderr, "unknow row sense %c.", toupper(sense));
    abort();
  }

  Cbc_addRowBuffer(model, nz, cols, coefs, rowLB, rowUB, name);
}

COINLIBAPI void COINLINKAGE
Osi_addRow(void *osi, const char *name, int nz,
  const int *cols, const double *coefs, char sense, double rhs)
{
  OsiSolverInterface *solver = (OsiSolverInterface *) osi;
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
    fprintf(stderr, "unknow row sense %c.", toupper(sense));
    abort();
  }

  solver->addRow(nz, cols, coefs, rowLB, rowUB);
  solver->setRowName(solver->getNumRows()-1, std::string(name));
}

COINLIBAPI void COINLINKAGE
Osi_setObjSense(void *osi, double sense)
{
  OsiSolverInterface *solver = (OsiSolverInterface *) osi;
  solver->setObjSense(sense);
}

COINLIBAPI void COINLINKAGE
Cbc_deleteRows(Cbc_Model *model, int numRows, const int rows[])
{
  Cbc_flush(model, FCRows);
  OsiSolverInterface *solver = model->model_->solver();

  if (model->rowNameIndex)
  {
    NameIndex &rowNameIndex = *((NameIndex  *)model->rowNameIndex);
    for ( int i=0 ; i<numRows; ++i )
      rowNameIndex.erase(solver->getRowName(rows[i]));
  }

  solver->deleteRows(numRows, rows);
}

COINLIBAPI void COINLINKAGE
Cbc_deleteCols(Cbc_Model *model, int numCols, const int cols[])
{
  Cbc_flush(model, FCColumns);
  OsiSolverInterface *solver = model->model_->solver();

  if (model->colNameIndex)
  {
    NameIndex &colNameIndex = *((NameIndex  *)model->colNameIndex);
    for ( int i=0 ; i<numCols; ++i )
      colNameIndex.erase(solver->getColName(cols[i]));
  }

  solver->deleteCols(numCols, cols);
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

/** @brief Creates a new OsiClpSolverInterface and returns a pointer to an OsiSolverInterface object */
COINLIBAPI void * COINLINKAGE
Osi_newSolver()
{
  OsiClpSolverInterface *clp = new OsiClpSolverInterface();

  return dynamic_cast<OsiSolverInterface *>(clp);
}

COINLIBAPI void COINLINKAGE
Osi_setObjCoef(void *osi, int index, double obj)
{
  OsiSolverInterface *osis = (OsiSolverInterface *)osi;
  osis->setObjCoeff( index, obj );
}

/** @brief Solves initial LP relaxation */
COINLIBAPI void COINLINKAGE
Osi_initialSolve(void *osi)
{
  OsiSolverInterface *osis = ( OsiSolverInterface *)osi;
  osis->initialSolve();
}

/** @brief Reoptimizes linear program  */
COINLIBAPI void COINLINKAGE
Osi_resolve(void *osi)
{
  OsiSolverInterface *osis = ( OsiSolverInterface *)osi;
  osis->resolve();
}

/** @brief Performs branch and bound */
COINLIBAPI void COINLINKAGE
Osi_branchAndBound(void *osi)
{
  OsiSolverInterface *osis = ( OsiSolverInterface *)osi;
  osis->branchAndBound();
}

// solution query methods

/** @brief Checks if optimization was abandoned */
COINLIBAPI char COINLINKAGE
Osi_isAbandoned(void *osi)
{
  OsiSolverInterface *osis = ( OsiSolverInterface *)osi;
  return (char)osis->isAbandoned();
}

/** @brief Checks if optimal solution was found */
COINLIBAPI char COINLINKAGE
Osi_isProvenOptimal(void *osi)
{
  OsiSolverInterface *osis = ( OsiSolverInterface *)osi;
  return (char)osis->isProvenOptimal();
}

/** @brief Checks if problem is primal infeasible */
COINLIBAPI char COINLINKAGE
Osi_isProvenPrimalInfeasible(void *osi)
{
  OsiSolverInterface *osis = ( OsiSolverInterface *)osi;
  return (char)osis->isProvenPrimalInfeasible();
}

/** @brief Checks if problem is dual infeasible */
COINLIBAPI char COINLINKAGE
Osi_isProvenDualInfeasible(void *osi)
{
  OsiSolverInterface *osis = ( OsiSolverInterface *)osi;
  return (char)osis->isProvenDualInfeasible();
}

/** @brief Checks if primal objective limit was reached */
COINLIBAPI char COINLINKAGE
Osi_isPrimalObjectiveLimitReached(void *osi)
{
  OsiSolverInterface *osis = ( OsiSolverInterface *)osi;
  return (char)osis->isPrimalObjectiveLimitReached();
}

/** @brief Checks if dual objective limit was reached */
COINLIBAPI char COINLINKAGE
Osi_isDualObjectiveLimitReached(void *osi)
{
  OsiSolverInterface *osis = ( OsiSolverInterface *)osi;
  return (char)osis->isDualObjectiveLimitReached();
}

/** @brief Checks if iteration limit was reached */
COINLIBAPI char COINLINKAGE
Osi_isIterationLimitReached(void *osi)
{
  OsiSolverInterface *osis = ( OsiSolverInterface *)osi;
  return (char)osis->isIterationLimitReached();
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

/** @brief Returns number non-zeros in the constraint matrix */
COINLIBAPI int COINLINKAGE
Osi_getNumNz( void *osi )
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;
  return osiSolver->getNumElements();
}

/** @brief Returns number integer/binary variables */
COINLIBAPI int COINLINKAGE
Osi_getNumIntegers( void *osi )
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;
  return osiSolver->getNumIntegers();
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
    fprintf(stderr, "unknow row sense %c.", toupper(sense));
    abort();
  }

  oc->insert(orc);
}

/** @brief Sets a variable to integer */
COINLIBAPI void COINLINKAGE
Osi_setInteger(void *osi, int index)
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;

  osiSolver->setInteger(index);
}

/** @brief Sets a variable to continuous */
COINLIBAPI void COINLINKAGE
Osi_setContinuous(void *osi, int index)
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;

  osiSolver->setContinuous(index);
}

COINLIBAPI int COINLINKAGE
Osi_getColNz(void *osi, int col)
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;
  const CoinPackedMatrix *cpmCol = osiSolver->getMatrixByCol();
  return cpmCol->getVectorLengths()[col];
}

/** @brief Indices of rows that a column appears 
     *
     * @param model problem object 
     * @param col column index
     * @return indices of rows that this column appears
     **/
COINLIBAPI const int *COINLINKAGE
Osi_getColIndices(void *osi, int col)
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;
  const CoinPackedMatrix *cpmCol = osiSolver->getMatrixByCol();
  const CoinBigIndex *starts = cpmCol->getVectorStarts();
  const int *cidx = cpmCol->getIndices() + starts[col];
  return cidx;
}

/** @brief Coefficients that a column appear in rows 
     *
     * @param model problem object 
     * @param col column index
     * @return coefficients of this column in rows
     **/
COINLIBAPI const double *COINLINKAGE
Osi_getColCoeffs(void *osi, int col)
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;

  const CoinPackedMatrix *cpmCol = osiSolver->getMatrixByCol();
  const CoinBigIndex *starts = cpmCol->getVectorStarts();
  const double *rcoef = cpmCol->getElements() + starts[col];
  return rcoef;
}

/** @brief Returns solution vector in OsiSolverInterface object */
COINLIBAPI const double * COINLINKAGE
Osi_getColSolution( void *osi )
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;

  return osiSolver->getColSolution();
}

/** Adds a new column */
COINLIBAPI void COINLINKAGE
Osi_addCol(void *osi, const char *name, double lb,
  double ub, double obj, char isInteger,
  int nz, int *rows, double *coefs)
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;

    osiSolver->addCol(nz, rows, coefs, lb, ub, obj, std::string(name));
    if (isInteger)
      osiSolver->setInteger(osiSolver->getNumCols() - 1);
}

/** @brief Returns vector of reduced costs */
COINLIBAPI const double * COINLINKAGE
Osi_getReducedCost( void *osi )
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;

  return osiSolver->getColSolution();
}

/** @brief Returns vector dual variables */
COINLIBAPI const double * COINLINKAGE
Osi_getRowPrice( void *osi )
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;

  return osiSolver->getRowPrice();
}

COINLIBAPI double COINLINKAGE
Osi_getObjValue( void *osi )
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;

  return osiSolver->getObjValue();
}

/** @brief Sets column upper bound */
COINLIBAPI void COINLINKAGE
Osi_setColUpper (void *osi, int elementIndex, double ub)
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;
  osiSolver->setColUpper(elementIndex, ub);
}

/** @brief Sets column upper bound */
COINLIBAPI void COINLINKAGE
Osi_setColLower(void *osi, int elementIndex, double lb)
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;
  osiSolver->setColLower(elementIndex, lb);
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

COINLIBAPI void COINLINKAGE
Cbc_storeNameIndexes(Cbc_Model *model, char _store)
{
  if (_store)
  {
    if (model->colNameIndex==NULL)
    {
      assert(model->rowNameIndex==NULL);
      model->colNameIndex = new NameIndex();
      model->rowNameIndex = new NameIndex();
    }
  }
  else
  {
    if (model->colNameIndex!=NULL)
    {
      NameIndex *m = (NameIndex *)model->colNameIndex;
      delete m;
      m = (NameIndex *)model->rowNameIndex;
      assert( m != NULL );
      delete m;

      model->colNameIndex = model->rowNameIndex = NULL;
    }
  }
}

COINLIBAPI int COINLINKAGE
Cbc_getColNameIndex(Cbc_Model *model, const char *name)
{
  if (!model->colNameIndex)
  {
    fprintf(stderr, "Call Cbc_storeNameIndex to enable name index search.");
    abort();
  }

  OsiSolverInterface *solver = model->model_->solver();
  NameIndex &colNameIndex = *((NameIndex  *)model->colNameIndex);
  NameIndex::iterator it = colNameIndex.find(std::string(name));
  if (it == colNameIndex.end())
    return -1;

  return it->second;
}

COINLIBAPI int COINLINKAGE
Cbc_getRowNameIndex(Cbc_Model *model, const char *name)
{
  if (!model->rowNameIndex)
  {
    fprintf(stderr, "Call Cbc_storeNameIndex to enable name index search.");
    abort();
  }

  OsiSolverInterface *solver = model->model_->solver();
  NameIndex &rowNameIndex = *((NameIndex  *)model->rowNameIndex);
  NameIndex::iterator it = rowNameIndex.find(std::string(name));
  if (it == rowNameIndex.end())
    return -1;

  return it->second;
}


#if defined(__MWERKS__)
#pragma export off
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
