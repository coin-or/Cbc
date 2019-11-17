// $Id$
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cmath>
#include <sstream>
#include <cfloat>
#include <climits>
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
#include "CglStored.hpp"
#include "CbcCutGenerator.hpp"
#include "ClpDualRowSteepest.hpp"
#include "ClpDualRowDantzig.hpp"
#include "ClpPEDualRowSteepest.hpp"
#include "ClpPEDualRowDantzig.hpp"
#include <OsiAuxInfo.hpp>

using namespace std;

static char **to_char_vec( const vector< string > names );
static void *xmalloc( const size_t size );
static void *xrealloc( void *ptr, const size_t newSize );

#define VALIDATE_ROW_INDEX(iRow, model)  \
      if (iRow<0 || iRow >= Cbc_getNumRows(model)) { \
        fprintf( stderr, "Invalid row index (%d), valid range is [0,%d). At %s:%d\n", \
            iRow, model->solver_->getNumRows()+model->nRows, __FILE__, __LINE__ ); \
        fflush(stderr); \
        abort(); \
      } 

#define VALIDATE_COL_INDEX(iColumn , model)  \
  if ( iColumn<0 || iColumn >= Cbc_getNumCols(model) ) { \
    fprintf( stderr, "Invalid column index (%d), valid range is [0,%d). At %s:%d\n", \
        iColumn, model->solver_->getNumCols()+model->nCols, __FILE__, __LINE__ ); \
    fflush(stderr); \
    abort(); \
  }

/* to check if solution (and of which type)
 * is available */
enum OptimizationTask {
  ModelNotOptimized      = 0,
  ContinuousOptimization = 1,
  IntegerOptimization    = 2
};

/**
  *
  *  C Interface Routines
  *
  */
#include "Cbc_C_Interface.h"


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

// command line arguments related to the lp relaxation solution
//static char cmdLPoptions[][80] = 
//  {"pertv", "idiot", "log", "seconds", "primalt", "dualt", "zerot", 
//   "pretol", "psi", "maxit", "crash", "scal" };
//static int nCmcLpOptions = 12;

struct Cbc_Model {
  /**
   * Problem is stored here: before optimizing 
   * this should be cloned because CBC modifies
   * directly this object dution the solution process
   **/
  OsiClpSolverInterface *solver_;

  vector< string > vcbcOptions; // to keep the order of the options
  map< string, string > cbcOptions; // to quickly check current value of option

  char relax_;

  // a new CbcModel needs to be created at every 
  // integer optimization, can be ignored for 
  // continuous optimization
  CbcModel *cbcModel_;

  // buffer for columns
  int colSpace;
  int nCols;
  int cNameSpace;
  int *cNameStart;
  char *cInt;
  char *cNames;
  double *cLB;
  double *cUB;
  double *cObj;

  vector< double > *iniSol;
  double iniObj;

  // buffer for rows
  int rowSpace;
  int nRows;
  int rNameSpace;
  int *rNameStart;
  char *rNames;
  double *rLB;
  double *rUB;
  int rElementsSpace;
  int *rStart;
  int *rIdx;
  double *rCoef;

  // for fast search of columns
  void *colNameIndex;
  void *rowNameIndex;

  cbc_incumbent_callback inc_callback;
  cbc_progress_callback progr_callback;
  cbc_callback userCallBack;

  cbc_cut_callback cut_callback;
  string cutCBName;
  void *cutCBData;
  int cutCBhowOften;
  char cutCBAtSol;

  enum OptimizationTask lastOptimization;


  /**
   * Incumbent callback callback data
   **/
  void *icAppData;

  /**
   * Progress callback callback data
   **/
  void *pgrAppData;

#ifdef CBC_THREAD
  pthread_mutex_t cbcMutex;
#endif


  int nSos;
  int sosCap;
  int sosSize;
  int sosElCap;
  int sosElSize;
  int *sosRowStart;
  int *sosType;
  int *sosEl;
  double *sosElWeight;

  int nColsMS;
  char **colNamesMS;
  double *colValuesMS;
  int charSpaceMS;

  // parameters
  enum LPMethod lp_method;
  enum DualPivot dualp;

  // lazy constraints
  CglStored *lazyConstrs;

  int int_param[N_INT_PARAMS];
  double dbl_param[N_DBL_PARAMS];
};

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
            charSize += (int)solver->getColName(i).size()+1;
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
  OsiSolverInterface *solver = model->solver_;

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

    model->cNameStart = (int *) xmalloc( sizeof(int)*c );
    assert( model->cNameStart );
    model->cNameStart[0] = 0;

    model->cInt = (char *) xmalloc( sizeof(char)*c );
    assert( model->cInt );

    model->cNames = (char *) xmalloc( sizeof(char)*model->cNameSpace );
    assert( model->cNames );

    model->cLB = (double *) xmalloc( sizeof(double)*c );
    assert( model->cLB );

    model->cUB = (double *)xmalloc( sizeof(double)*c );
    assert( model->cUB );

    model->cObj = (double *)xmalloc( sizeof(double)*c );
    assert( model->cObj );
  }
  else
  {
    // check buffer space
    if (model->nCols+2 >= model->colSpace)
    {
      model->colSpace *= 2;
      int c = model->colSpace;

      model->cNameStart = (int *) xrealloc( model->cNameStart, sizeof(int)*c );
      assert( model->cNameStart );

      model->cInt = (char *) xrealloc( model->cInt, sizeof(char)*c );
      assert( model->cInt );

      model->cLB = (double *) xrealloc( model->cLB, sizeof(double)*c );
      assert( model->cLB );

      model->cUB = (double *) xrealloc( model->cUB, sizeof(double)*c );
      assert( model->cUB );

      model->cObj = (double *) xrealloc( model->cObj, sizeof(double)*c );
      assert( model->cObj );
    }
    // check string buffer space
    int slen = additionlNameSpace + 1;
    int reqsize = slen + model->cNameStart[model->nCols]+1;
    if (reqsize>model->cNameSpace)
    {
      model->cNameSpace *= 2;
      model->cNames = (char *) xrealloc( model->cNames, sizeof(char)*model->cNameSpace );
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
  int len = (int)strlen(name);

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
    model->rStart = (int *)xmalloc(sizeof(int)*model->rowSpace);
    model->rStart[0] = 0;
    model->rLB = (double *)xmalloc(sizeof(double)*model->rowSpace);
    model->rUB = (double *)xmalloc(sizeof(double)*model->rowSpace);
    model->rNameStart = (int *)xmalloc(sizeof(int)*model->rowSpace);
    model->rNameStart[0] = 0;

    model->rElementsSpace = std::max(131072, nzRow * 2);
    model->rIdx = (int *)xmalloc(sizeof(int)*model->rElementsSpace);
    model->rCoef = (double *)xmalloc(sizeof(double)*model->rElementsSpace);

    model->rNameSpace = 131072;
    model->rNames = (char *)xmalloc(sizeof(char)*model->rNameSpace);
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
        model->rStart = (int *)xrealloc(model->rStart, sizeof(int)*model->rowSpace);
        model->rLB = (double *)xrealloc(model->rLB, sizeof(double)*model->rowSpace);
        model->rUB = (double *)xrealloc(model->rUB, sizeof(double)*model->rowSpace);
        model->rNameStart = (int *)xrealloc(model->rNameStart, sizeof(int)*model->rowSpace);
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
        model->rIdx = (int *)xrealloc(model->rIdx, sizeof(int)*model->rElementsSpace);
        model->rCoef = (double *)xrealloc(model->rCoef, sizeof(double)*model->rElementsSpace);
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
        model->rNames = (char *)xrealloc(model->rNames, sizeof(char)*model->rNameSpace);
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
  int nameLen = (int)strlen(rName);
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

  // SOS buffer
  model->nSos = 0;
  model->sosCap = 0;
  model->sosSize = 0;
  model->sosElCap = 0;
  model->sosElSize = 0;
  model->sosRowStart = NULL;
  model->sosEl = NULL;
  model->sosElWeight = NULL;
  model->sosType = NULL;

  model->nColsMS = 0;
  model->colNamesMS = NULL;
  model->colValuesMS = NULL;
}

static void Cbc_iniParams( Cbc_Model *model );

/* Default Cbc_Model constructor */
COINLIBAPI Cbc_Model *COINLINKAGE
Cbc_newModel()
{
  Cbc_Model *model = new Cbc_Model();

  Cbc_iniParams(model);
  
  model->solver_ = new OsiClpSolverInterface();
  model->relax_ = 0;

  model->cbcModel_ = NULL;

  model->inc_callback = NULL;
  model->progr_callback = NULL;
  model->userCallBack = NULL;

  model->cut_callback = NULL;
  model->cutCBData = NULL;
  model->cutCBhowOften = -1;
  model->cutCBAtSol = 0;
  
  model->lastOptimization = ModelNotOptimized;

  model->icAppData = NULL;
  model->pgrAppData = NULL;
  model->colNameIndex = NULL;
  model->rowNameIndex = NULL;
  model->iniObj = DBL_MAX;

  model->iniSol = NULL;


  Cbc_iniBuffer(model);

#ifdef CBC_THREAD
  pthread_mutex_init(&(model->cbcMutex), NULL);
#endif

  model->lazyConstrs = NULL;

  return model;
}

/* Cbc_Model Destructor */
COINLIBAPI void COINLINKAGE
Cbc_deleteModel(Cbc_Model *model)
{
  Cbc_deleteColBuffer(model);
  Cbc_deleteRowBuffer(model);

  if (model->sosRowStart) {
    free(model->sosRowStart);
    free(model->sosEl);
    free(model->sosElWeight);
    free(model->sosType);
  }

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

  if (model->cbcModel_)
    delete model->cbcModel_;

  if (model->solver_)
    delete model->solver_;

  if (model->iniSol)
    delete model->iniSol;

  if (model->nColsMS) {
    if (model->colNamesMS) {
      free(model->colNamesMS[0]);
      free(model->colNamesMS);
    }
    free(model->colValuesMS);
  }

  if (model->lazyConstrs)
    delete model->lazyConstrs;

  delete model;
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
  OsiSolverInterface *solver = model->solver_;

  solver->loadProblem(numcols, numrows, start, index, value,
    collb, colub, obj, rowlb, rowub);
} //  Cbc_loadProblem()

/* should be called after reading a new problem */
static void fillAllNameIndexes(Cbc_Model *model)
{
  if (!model->colNameIndex)
    return;

  OsiSolverInterface *solver = model->solver_;
  NameIndex &colNameIndex = *((NameIndex  *)model->colNameIndex);
  colNameIndex.clear();
  NameIndex &rowNameIndex = *((NameIndex  *)model->rowNameIndex);
  rowNameIndex.clear();

  for ( int i=0 ; (i<solver->getNumCols()) ; ++i )
    colNameIndex[solver->getColName(i)] = i;

  for ( int i=0 ; (i<solver->getNumRows()) ; ++i )
    rowNameIndex[solver->getRowName(i)] = i;
}

/** Reads an MPS file
 *
 * @param model problem object
 * @param fileName file name
 **/
COINLIBAPI int COINLINKAGE
Cbc_readMps(Cbc_Model *model, const char *filename)
{
  int result = 1;
  OsiSolverInterface *solver = model->solver_;
  result = solver->readMps(filename);
  assert(result == 0);

  Cbc_deleteColBuffer(model);
  Cbc_deleteRowBuffer(model);
  Cbc_iniBuffer(model);

  fillAllNameIndexes(model);

  return result;
}

/** Writes an MPS file
 *
 * @param model problem object
 * @param fileName file name
 **/
COINLIBAPI void COINLINKAGE
Cbc_writeMps(Cbc_Model *model, const char *filename)
{
  Cbc_flush(model);

  model->solver_->writeMps(filename, "mps", Cbc_getObjSense(model));
}

/** Writes an LP file
 *
 * @param model problem object
 * @param fileName file name */
COINLIBAPI void COINLINKAGE
Cbc_writeLp(Cbc_Model *model, const char *filename)
{
  Cbc_flush(model);

  char outFile[512];
  strncpy(outFile, filename, 511);

  char *s = NULL;
  if ((s = strstr(outFile, ".lp"))) {
    if (s != outFile) // not at the start
      *s = '\0';
  }

  model->solver_->writeLp(outFile);
}

/** Reads an LP file
 *
 * @param model problem object
 * @param fileName file name */
COINLIBAPI int COINLINKAGE
Cbc_readLp(Cbc_Model *model, const char *filename)
{
  int result = 1;
  OsiSolverInterface *solver = model->solver_;
  result = solver->readLp(filename);
  assert(result == 0);

  Cbc_deleteColBuffer(model);
  Cbc_deleteRowBuffer(model);
  Cbc_iniBuffer(model);

  fillAllNameIndexes(model);

  return result;
}

/* Sets an initial feasible solution
 *
 * @param model problem object
 **/
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

  if (model->iniSol) {
    model->iniSol->resize( Cbc_getNumCols(model) );
    double *iniSol = &((*model->iniSol)[0]);
    memcpy( iniSol, sol, sizeof(double)*Cbc_getNumCols(model) );
  } else {
    model->iniSol = new vector<double>(sol, sol+n);
  }

  model->iniObj = objval;
}


COINLIBAPI void COINLINKAGE
Cbc_setIntParam(Cbc_Model *model, enum IntParam which, const int val) {
  model->int_param[which] = val;
}

COINLIBAPI void COINLINKAGE
Cbc_setDblParam(Cbc_Model *model, enum DblParam which, const double val) {
  model->dbl_param[which] = val;
}

COINLIBAPI void COINLINKAGE
Cbc_setParameter(Cbc_Model *model, const char *name, const char *value)
{
  if (model->cbcOptions.find(string(name))==model->cbcOptions.end())
    model->vcbcOptions.push_back(string(name));
  model->cbcOptions[name] = string(value);
}

/* Fills in array with problem name  */
COINLIBAPI void COINLINKAGE
Cbc_problemName(Cbc_Model *model, int maxNumberCharacters, char *array)
{
  std::string name;
  model->solver_->getStrParam(OsiProbName, name);
  strncpy(array, name.c_str(), maxNumberCharacters);
}
/* Sets problem name.  Must have \0 at end.  */
COINLIBAPI int COINLINKAGE
Cbc_setProblemName(Cbc_Model *model, const char *array)
{
  bool result = false;
  result = model->solver_->setStrParam(OsiProbName, array);

  return (result) ? 1 : 0;
}

COINLIBAPI int COINLINKAGE Cbc_status(Cbc_Model *model) {
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Status not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      fprintf( stderr, "Cbc_status only available for MIP models.\n");
      abort();
      break;
    case IntegerOptimization:
      return model->cbcModel_->status();
  }

  return INT_MAX;
}


COINLIBAPI int COINLINKAGE
Cbc_secondaryStatus(Cbc_Model *model) {
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Status not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      fprintf( stderr, "Cbc_status only available for MIP models.\n");
      abort();
      break;
    case IntegerOptimization:
      return model->cbcModel_->secondaryStatus();
  }

  return INT_MAX;
}

/* Number of elements in matrix */
COINLIBAPI int COINLINKAGE
Cbc_getNumElements(Cbc_Model *model)
{
  Cbc_flush(model);

  int result = 0;
  result = model->solver_->getNumElements();

  return result;
}

COINLIBAPI int COINLINKAGE
Cbc_getNumIntegers(Cbc_Model *model)
{
  Cbc_flush(model, FCColumns);

  return model->solver_->getNumIntegers();
}

// Column starts in matrix
COINLIBAPI const CoinBigIndex *COINLINKAGE
Cbc_getVectorStarts(Cbc_Model *model)
{
  Cbc_flush(model);

  const CoinPackedMatrix *matrix = NULL;
  matrix = model->solver_->getMatrixByCol();
  return (matrix == NULL) ? NULL : matrix->getVectorStarts();
}

// Row indices in matrix
COINLIBAPI const int *COINLINKAGE
Cbc_getIndices(Cbc_Model *model)
{
  Cbc_flush(model, FCRows);

  const int *result = NULL;
  const CoinPackedMatrix *matrix = NULL;
  matrix = model->solver_->getMatrixByCol();
  result = (matrix == NULL) ? NULL : matrix->getIndices();

  return result;
}

// Element values in matrix
COINLIBAPI const double *COINLINKAGE
Cbc_getElements(Cbc_Model *model)
{
  Cbc_flush(model, FCRows);

  const double *result = NULL;
  const CoinPackedMatrix *matrix = NULL;
  matrix = model->solver_->getMatrixByCol();
  result = (matrix == NULL) ? NULL : matrix->getElements();

  return result;
}
// ======================================================================

/* Pass in Callback function */
COINLIBAPI void COINLINKAGE
Cbc_registerCallBack(Cbc_Model *model,
  cbc_callback userCallBack)
{
  model->userCallBack = userCallBack;
}
/* Unset Callback function */
COINLIBAPI void COINLINKAGE
Cbc_clearCallBack(Cbc_Model *model)
{
  model->userCallBack = NULL;
}

/* length of names (0 means no names0 */
COINLIBAPI size_t COINLINKAGE
Cbc_maxNameLength(Cbc_Model *model)
{
  size_t result = 0;
  OsiSolverInterface::OsiNameVec const &rownames = model->solver_->getRowNames();
  for (size_t i = 0; i < rownames.size(); i++) {
    if (rownames[i].length() > result)
      result = rownames[i].length();
  }
  OsiSolverInterface::OsiNameVec const &colnames = model->solver_->getColNames();
  for (size_t i = 0; i < colnames.size(); i++) {
    if (colnames[i].length() > result)
      result = colnames[i].length();
  }

  // go trough buffered names also
  for ( size_t i=0 ; (i<(size_t)model->nCols) ; ++i )
    result = max( result, strlen(model->cNames+model->cNameStart[i]) );

  for ( size_t i=0 ; (i<(size_t)model->nRows) ; ++i )
    result = max( result, strlen(model->rNames+model->rNameStart[i]) );

  return result;
}

COINLIBAPI void COINLINKAGE
Cbc_getRowName(Cbc_Model *model, int iRow, char *name, size_t maxLength)
{
  VALIDATE_ROW_INDEX( iRow, model );

  if (iRow < model->solver_->getNumRows()) {
    std::string rowname = model->solver_->getRowName(iRow);
    strncpy(name, rowname.c_str(), maxLength);
    name[maxLength - 1] = '\0';
  } else {
    int idxRowBuffer = iRow - model->solver_->getNumRows();
    strncpy( name, model->rNames+model->rNameStart[idxRowBuffer], maxLength );
  }
}

COINLIBAPI void COINLINKAGE
Cbc_getColName(Cbc_Model *model, int iColumn, char *name, size_t maxLength)
{
  VALIDATE_COL_INDEX( iColumn, model );

  if (iColumn < model->solver_->getNumCols()) {
    std::string colname = model->solver_->getColName(iColumn);
    strncpy(name, colname.c_str(), maxLength);
    name[maxLength - 1] = '\0';
  } else {
    int idxColBuffer = iColumn - model->solver_->getNumCols();
    strncpy( name, model->cNames+model->cNameStart[idxColBuffer], maxLength );
  }

}

COINLIBAPI void COINLINKAGE
Cbc_setColName(Cbc_Model *model, int iColumn, const char *name)
{
  VALIDATE_COL_INDEX( iColumn, model );

  Cbc_flush(model);
  OsiSolverInterface *solver = model->solver_;
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
  VALIDATE_ROW_INDEX( iRow, model );

  Cbc_flush(model, FCRows);
  OsiSolverInterface *solver = model->solver_;
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

// adds all sos objects to the current cbcModel_ object
static void Cbc_addAllSOS( Cbc_Model *model );

// adds mipstart if available
static void Cbc_addMS( Cbc_Model *model );


COINLIBAPI int COINLINKAGE
Cbc_solve(Cbc_Model *model)
{
  Cbc_flush( model );

  OsiSolverInterface *solver = model->solver_;

  solver->setDblParam( OsiPrimalTolerance, model->dbl_param[DBL_PARAM_PRIMAL_TOL]);
  solver->setDblParam( OsiDualTolerance, model->dbl_param[DBL_PARAM_DUAL_TOL]);

  OsiClpSolverInterface *clpSolver = dynamic_cast< OsiClpSolverInterface * >(solver);
  assert(clpSolver);
  ClpSimplex *clps = clpSolver->getModelPtr();
  assert(clps);
  clps->setPerturbation(model->int_param[INT_PARAM_PERT_VALUE]);
  double maxTime = Cbc_getMaximumSeconds(model);
  if (maxTime != COIN_DBL_MAX)
    clps->setMaximumWallSeconds(maxTime);
  solver->messageHandler()->setLogLevel( model->int_param[INT_PARAM_LOG_LEVEL] );

  if (solver->getNumIntegers() == 0 || model->relax_ == 1) {
    model->lastOptimization = ContinuousOptimization;

    if (solver->basisIsAvailable()) {
      solver->resolve();
      return 0;
    } // resolve
  } // solve only lp relaxation

  /* for integer or linear optimization starting with LP relaxation */
  ClpSolve clpOptions;
  switch (model->lp_method) {
    case LPM_Auto:
      break;
    case LPM_Primal:
      clpOptions.setSolveType( ClpSolve::usePrimal );
      break;
    case LPM_Dual:
      clpOptions.setSolveType( ClpSolve::useDual );
      break;
    case LPM_Barrier:
      clpOptions.setSolveType( ClpSolve::useBarrier );
      clpOptions.setSpecialOption(4, 4);
      break;
  }
  if (model->int_param[INT_PARAM_IDIOT] != -1) {
    clpOptions.setSpecialOption(1, 2, model->int_param[INT_PARAM_IDIOT]);
  }
  clpSolver->setSolveOptions(clpOptions);

  switch (model->dualp) {
    case DP_Auto:
      {
        ClpDualRowSteepest asteep(3);
        clps->setDualRowPivotAlgorithm(asteep);
        break;
      }
    case DP_Dantzig:
      {
        ClpDualRowDantzig dantzig;
        clps->setDualRowPivotAlgorithm(dantzig);
        break;
      }
    case DP_Partial:
      {
        ClpDualRowSteepest bsteep(2);
        clps->setDualRowPivotAlgorithm(bsteep);
        break;
      }
    case DP_Steepest:
      {
        ClpDualRowSteepest csteep;
        clps->setDualRowPivotAlgorithm(csteep);
        break;
      }
    case DP_PESteepest:
      {
        ClpPEDualRowSteepest p(model->dbl_param[DBL_PARAM_PSI]);
        clps->setDualRowPivotAlgorithm(p);
        break;
      }
  }

  solver->initialSolve();

  if (solver->isProvenPrimalInfeasible() || solver->isProvenDualInfeasible() ||
      solver->isAbandoned() || solver->isIterationLimitReached()) {
    model->lastOptimization = ContinuousOptimization;
    return 0;
  }

  if (solver->getNumIntegers() == 0 || model->relax_ == 1) {
    return 0;
  }

  // MIP Optimization
  if (model->cbcModel_) {
    delete model->cbcModel_;
  }

  OsiClpSolverInterface *linearProgram = dynamic_cast<OsiClpSolverInterface *>( model->solver_->clone() );
  model->lastOptimization = IntegerOptimization;
  CbcModel *cbcModel = model->cbcModel_ = new CbcModel( *linearProgram );

  try {
    if (model->lazyConstrs) {
      cbcModel->addCutGenerator(model->lazyConstrs, 1, "Stored LazyConstraints", true, 1);
      model->cutCBAtSol = 1;
    }

    Cbc_EventHandler *cbc_eh = NULL;
    if (model->inc_callback!=NULL || model->progr_callback!=NULL)
    {
      cbc_eh = new Cbc_EventHandler(model->cbcModel_);
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

      cbcModel->passInEventHandler(cbc_eh);
    }

    if (model->iniSol)
      cbcModel->setBestSolution(&((*model->iniSol)[0]), Cbc_getNumCols(model), model->iniObj, true);

    // add cut generator if necessary
    if (model->cut_callback) {
      cbcModel->setKeepNamesPreproc(true);

      CglCallback cglCb;
      cglCb.appdata = model->cutCBData;
      cglCb.cut_callback_ = model->cut_callback;
#ifdef CBC_THREAD
      cglCb.cbcMutex = &(model->cbcMutex);
#endif
      cbcModel->addCutGenerator( &cglCb, model->cutCBhowOften, model->cutCBName.c_str(), true, model->cutCBAtSol );
    }

    Cbc_MessageHandler *cbcmh  = NULL;

    if (model->userCallBack) {
      cbcmh = new Cbc_MessageHandler(*cbcModel->messageHandler());
      cbcmh->setCallBack(model->userCallBack);
      cbcmh->setModel(model);
      cbcModel->passInMessageHandler(cbcmh);
    }

    CbcSolverUsefulData cbcData;
    CbcMain0(*cbcModel, cbcData);

    cbcModel->solver()->setDblParam( OsiPrimalTolerance, model->dbl_param[DBL_PARAM_PRIMAL_TOL] );
    cbcModel->solver()->setDblParam( OsiDualTolerance, model->dbl_param[DBL_PARAM_DUAL_TOL] );
    // adds SOSs if any
    Cbc_addAllSOS(model);

    // adds MIPStart if any
    Cbc_addMS( model );

    // parameters
    if (model->dbl_param[DBL_PARAM_TIME_LIMIT] != COIN_DBL_MAX)
      cbcModel->setMaximumSeconds( model->dbl_param[DBL_PARAM_TIME_LIMIT] );
    if ( model->int_param[INT_PARAM_MAX_SOLS] != INT_MAX && model->int_param[INT_PARAM_MAX_SOLS] != -1 )
      cbcModel->setMaximumSolutions( model->int_param[INT_PARAM_MAX_SOLS] );
    cbcModel->setAllowableGap( model->dbl_param[DBL_PARAM_ALLOWABLE_GAP] );
    cbcModel->setAllowableFractionGap( model->dbl_param[DBL_PARAM_GAP_RATIO] );
    if ( model->int_param[INT_PARAM_MAX_NODES] != INT_MAX )
      cbcModel->setMaximumNodes( model->int_param[INT_PARAM_MAX_NODES] );
    cbcModel->setLogLevel( model->int_param[INT_PARAM_LOG_LEVEL] );
    if ( model->dbl_param[DBL_PARAM_CUTOFF] != COIN_DBL_MAX )
      cbcModel->setCutoff( model->dbl_param[DBL_PARAM_CUTOFF] );

    std::vector< string > argv;
    argv.push_back("Cbc_C_Interface");

    /*
    char custom = 0;

    if (model->int_param[INT_PARAM_PERT_VALUE] != 50) {
      custom = 1;
      char str[256];
      sprintf(str, "-pertv=%d", model->int_param[INT_PARAM_PERT_VALUE]);
      argv.push_back(str);
    }
    if (model->int_param[INT_PARAM_IDIOT] != -1) {
      custom = 1;
      char str[256];
      sprintf(str, "-idiot=%d", model->int_param[INT_PARAM_IDIOT]);
      argv.push_back(str);
    }
    if (model->dbl_param[DBL_PARAM_TIME_LIMIT] != COIN_DBL_MAX) {
      char str[256];
      sprintf(str, "-seco=%g", model->dbl_param[DBL_PARAM_TIME_LIMIT]);
      argv.push_back(str);
    }  */


    // LP options should be included here
    /*
    switch (model->dualp) {
      case DP_Auto:
        break;
      case DP_Partial:
        argv.push_back("-dualp=partial");
        custom = 1;
        break;
      case DP_PESteepest:
        custom = 1;
        char str[256];
        sprintf(str, "-psi=%g", model->dbl_param[DBL_PARAM_PSI]);
        argv.push_back(str);
        argv.push_back("-dualp=pesteep");
        break;
      case DP_Dantzig:
        custom = 1;
        argv.push_back("-dualp=dantzig");
        break;
      case DP_Steepest:
        custom = 1;
        argv.push_back("-dualp=steep");
        break;
    }

    switch (model->lp_method) {
      case LPM_Auto:
        if (custom == 0)
          argv.push_back("-guess");
          argv.push_back("-initialSolve");
        break;
      case LPM_Primal:
        argv.push_back("-primals");
        break;
      case LPM_Dual:
        argv.push_back("-duals");
        break;
      case LPM_Barrier:
        argv.push_back("-chol=univ");
        argv.push_back("-barrier");
        break;

    } */

    for ( size_t i=0 ; i<model->vcbcOptions.size() ; ++i ) {
      string param = model->vcbcOptions[i];
      string val = model->cbcOptions[param];
      if (val.size()) {
        stringstream ss;
        ss << "-" << param << "=" << val;
        argv.push_back(ss.str().c_str());
      } else {
        stringstream ss;
        ss << "-" << param;
        argv.push_back(ss.str());
      }
    }

    argv.push_back("-solve");
    argv.push_back("-quit");

    cbcData.noPrinting_= false;

    char **charCbcOpts = to_char_vec(argv);
    const int nargs = (int) argv.size();
    const char **args = (const char **)charCbcOpts;

    OsiBabSolver defaultC;
    if (model->cutCBAtSol) {
      defaultC.setSolverType(4);
      model->cbcModel_->solver()->setAuxiliaryInfo(&defaultC);
      model->cbcModel_->passInSolverCharacteristics(&defaultC);
    }

    CbcMain1( nargs, args, *model->cbcModel_, cbc_callb, cbcData );

    free(charCbcOpts);
    delete linearProgram;

    if (cbc_eh)
      delete cbc_eh;

    if (cbcmh)
      delete cbcmh;
  } catch (CoinError &e) {
    fprintf( stderr, "%s ERROR: %s::%s, %s\n", "Cbc_solve",
      e.className().c_str(), e.methodName().c_str(), e.message().c_str());
    abort();
  }
  
  return cbcModel->status();
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
    Cbc_Model *model, 
    cbc_cut_callback cutcb, 
    const char *name, 
    void *appData, 
    int howOften,
    char atSolution )
{
  model->cut_callback = cutcb;
  model->cutCBName = string(name);
  model->cutCBData = appData;
  model->cutCBhowOften = howOften;
  model->cutCBAtSol = atSolution;
}

COINLIBAPI int COINLINKAGE
Cbc_getNumCols(Cbc_Model *model)
{
  return model->solver_->getNumCols() + model->nCols;
}

COINLIBAPI int COINLINKAGE
Cbc_getNumRows(Cbc_Model *model)
{
  return model->solver_->getNumRows() + model->nRows;
}


COINLIBAPI int COINLINKAGE
Cbc_getIterationCount(Cbc_Model *model) {
  return model->cbcModel_->getIterationCount();
}

/** Number of non-zero entries in a row */
COINLIBAPI int COINLINKAGE
Cbc_getRowNz(Cbc_Model *model, int row)
{
  VALIDATE_ROW_INDEX( row, model);

  if (row<model->solver_->getNumRows()) {
    const CoinPackedMatrix *cpmRow = model->solver_->getMatrixByRow();
    return cpmRow->getVectorLengths()[row];
  } else {
    int idxRowBuffer = row - model->solver_->getNumRows();
    return model->rStart[idxRowBuffer+1]-model->rStart[idxRowBuffer];
  }
}

/** Indices of variables that appear on this row */
COINLIBAPI const int *COINLINKAGE
Cbc_getRowIndices(Cbc_Model *model, int row)
{
  VALIDATE_ROW_INDEX( row, model);

  if (row<model->solver_->getNumRows()) {
    const CoinPackedMatrix *cpmRow = model->solver_->getMatrixByRow();
    const CoinBigIndex *starts = cpmRow->getVectorStarts();
    const int *ridx = cpmRow->getIndices() + starts[row];
    return ridx;
  } else {
    int idxRowBuffer = row - model->solver_->getNumRows();
    return model->rIdx + model->rStart[idxRowBuffer];
  }
}

/** Coefficients of variables that appear on this row */
COINLIBAPI const double *COINLINKAGE
Cbc_getRowCoeffs(Cbc_Model *model, int row)
{
  VALIDATE_ROW_INDEX( row, model);

  if (row<model->solver_->getNumRows()) {
    const CoinPackedMatrix *cpmRow = model->solver_->getMatrixByRow();
    const CoinBigIndex *starts = cpmRow->getVectorStarts();
    const double *rcoef = cpmRow->getElements() + starts[row];
    return rcoef;
  } else {
    int idxRowBuffer = row - model->solver_->getNumRows();
    return model->rCoef + model->rStart[idxRowBuffer];
  }
}

/** Number of non-zero entries in a column */
COINLIBAPI int COINLINKAGE
Cbc_getColNz(Cbc_Model *model, int col)
{
  VALIDATE_COL_INDEX( col, model );

  Cbc_flush(model);
  const CoinPackedMatrix *cpmCol = model->solver_->getMatrixByCol();
  return cpmCol->getVectorLengths()[col];
}

/** Indices of rows that a column appears */
COINLIBAPI const int *COINLINKAGE
Cbc_getColIndices(Cbc_Model *model, int col)
{
  VALIDATE_COL_INDEX( col, model );

  Cbc_flush(model);
  const CoinPackedMatrix *cpmCol = model->solver_->getMatrixByCol();
  const CoinBigIndex *starts = cpmCol->getVectorStarts();
  const int *cidx = cpmCol->getIndices() + starts[col];
  return cidx;
}

/** Coefficients that a column appear in rows */
COINLIBAPI const double *COINLINKAGE
Cbc_getColCoeffs(Cbc_Model *model, int col)
{
  VALIDATE_COL_INDEX( col, model );

  Cbc_flush(model);

  const CoinPackedMatrix *cpmCol = model->solver_->getMatrixByCol();
  const CoinBigIndex *starts = cpmCol->getVectorStarts();
  const double *rcoef = cpmCol->getElements() + starts[col];
  return rcoef;
}

/** Right hand side of a row */
COINLIBAPI double COINLINKAGE
Cbc_getRowRHS(Cbc_Model *model, int row)
{
  VALIDATE_ROW_INDEX( row, model );

  if (row<model->solver_->getNumRows()) {
    return model->solver_->getRightHandSide()[row];
  } else {
    int idxRowBuffer = row - model->solver_->getNumRows();
    if (model->rUB[idxRowBuffer] < COIN_DBL_MAX) 
      return model->rUB[idxRowBuffer];
    else 
      return model->rLB[idxRowBuffer];
  }
}

/** Sense a row */
COINLIBAPI char COINLINKAGE
Cbc_getRowSense(Cbc_Model *model, int row)
{
  VALIDATE_ROW_INDEX( row, model );

  if (row<model->solver_->getNumRows()) {
    return model->solver_->getRowSense()[row];
  } else {
    int idxRowBuffer = row - model->solver_->getNumRows();
    if (fabs(model->rLB[idxRowBuffer]-model->rUB[idxRowBuffer]) <= 1e-15)
      return 'E';
    if (model->rUB[idxRowBuffer] == COIN_DBL_MAX) 
      return 'G';
    if (model->rLB[idxRowBuffer] == -COIN_DBL_MAX) 
      return 'L';

    return 'R';
  }
}

/** Are there a numerical difficulties? */
COINLIBAPI int COINLINKAGE
Cbc_isAbandoned(Cbc_Model *model)
{
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Information not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      return model->solver_->isAbandoned();
      break;
    case IntegerOptimization:
      return model->cbcModel_->isAbandoned();
  }

  return false;
}

/** Is optimality proven? */
COINLIBAPI int COINLINKAGE
Cbc_isProvenOptimal(Cbc_Model *model)
{
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Information not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      return model->solver_->isProvenOptimal();
    case IntegerOptimization:
      return model->cbcModel_->isProvenOptimal();
  }

  return false;
}

COINLIBAPI int COINLINKAGE
Cbc_isProvenInfeasible(Cbc_Model *model)
{
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Information not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      return model->solver_->isProvenPrimalInfeasible();
    case IntegerOptimization:
      if (model->cbcModel_->status() == -1)
        return model->solver_->isProvenPrimalInfeasible();
      return model->cbcModel_->isProvenInfeasible();
  }

  return false;
}

COINLIBAPI double COINLINKAGE
Cbc_getObjValue(Cbc_Model *model)
{
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Information not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      return model->solver_->getObjValue();
    case IntegerOptimization:
      return model->cbcModel_->getObjValue();
  }

  return COIN_DBL_MAX;
}

COINLIBAPI const double *COINLINKAGE
Cbc_getReducedCost(Cbc_Model *model)
{
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Information not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      return model->solver_->getReducedCost();
    case IntegerOptimization:
      fprintf( stderr, "Information only available when optimizing continuous models.\n");
      abort();
      break;
  }

  return NULL;
}

COINLIBAPI const double *COINLINKAGE
Cbc_getRowPrice(Cbc_Model *model)
{
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Information not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      return model->solver_->getRowPrice();
    case IntegerOptimization:
      fprintf( stderr, "Information only available when optimizing continuous models.\n");
      abort();
      break;
  }

  return NULL;
}

COINLIBAPI int COINLINKAGE
Cbc_numberSavedSolutions(Cbc_Model *model)
{
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Information not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      if (model->solver_->isProvenOptimal())
        return 1;
      else
        return 0;
    case IntegerOptimization:
      return model->cbcModel_->numberSavedSolutions();
  }

  return 0;
}

COINLIBAPI const double *COINLINKAGE
Cbc_savedSolution(Cbc_Model *model, int whichSol)
{
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Information not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      fprintf( stderr, "Information only available when optimizing integer models.\n");
      abort();
      break;
    case IntegerOptimization:
      return model->cbcModel_->savedSolution(whichSol);
  }

  return NULL;
}

COINLIBAPI double COINLINKAGE
Cbc_savedSolutionObj(Cbc_Model *model, int whichSol)
{
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Information not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      fprintf( stderr, "Information only available when optimizing integer models.\n");
      abort();
      break;
    case IntegerOptimization:
      return model->cbcModel_->savedSolutionObjective(whichSol);
  }

  return COIN_DBL_MAX;
}

COINLIBAPI const double *COINLINKAGE
Cbc_getColSolution(Cbc_Model *model)
{
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Information not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      return model->solver_->getColSolution();
    case IntegerOptimization:
      return model->cbcModel_->getColSolution();
  }

  return NULL;
}


COINLIBAPI int COINLINKAGE
Cbc_isContinuousUnbounded(Cbc_Model *model) {
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Information not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      return model->solver_->isProvenDualInfeasible();
    case IntegerOptimization:
      if (model->cbcModel_->status() == -1)
        return model->solver_->isProvenDualInfeasible();

      return model->cbcModel_->isContinuousUnbounded();
  }

  return false;
}


COINLIBAPI int COINLINKAGE
Cbc_isNodeLimitReached(Cbc_Model *model) {
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Information not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      fprintf( stderr, "Information not available when optimizing continuous models.\n");
      abort();
      break;
    case IntegerOptimization:
      return model->cbcModel_->isNodeLimitReached();
  }

  return false;
}


COINLIBAPI int COINLINKAGE
Cbc_isSecondsLimitReached(Cbc_Model *model) {
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Information not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      fprintf( stderr, "Information only available when optimizing integer models.\n");
      abort();
      break;
    case IntegerOptimization:
      return model->cbcModel_->isSecondsLimitReached();
  }

  return false;

}

COINLIBAPI int COINLINKAGE
Cbc_isSolutionLimitReached(Cbc_Model *model) {
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Information not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      fprintf( stderr, "Information only available when optimizing integer models.\n");
      abort();
      break;
    case IntegerOptimization:
      return model->cbcModel_->isSolutionLimitReached();
  }

  return false;
}




COINLIBAPI int COINLINKAGE
Cbc_isInitialSolveAbandoned(Cbc_Model *model) {
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Information not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      fprintf( stderr, "Information only available when optimizing integer models.\n");
      abort();
      break;
    case IntegerOptimization:
      return model->cbcModel_->isInitialSolveAbandoned();
  }

  return false;
}

COINLIBAPI int COINLINKAGE
Cbc_isInitialSolveProvenOptimal(Cbc_Model *model) {
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Information not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      fprintf( stderr, "Information only available when optimizing integer models.\n");
      abort();
      break;
    case IntegerOptimization:
      return model->cbcModel_->isInitialSolveProvenOptimal();
  }

  return false;
}


COINLIBAPI int COINLINKAGE
Cbc_isInitialSolveProvenPrimalInfeasible(Cbc_Model *model) {
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Information not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      fprintf( stderr, "Information only available when optimizing integer models.\n");
      abort();
      break;
    case IntegerOptimization:
      return model->cbcModel_->isInitialSolveProvenPrimalInfeasible();
  }

  return false;

}

COINLIBAPI double COINLINKAGE
Cbc_getObjSense(Cbc_Model *model) {
  return model->solver_->getObjSense();
}

COINLIBAPI void COINLINKAGE
Cbc_setObjSense(Cbc_Model *model, double sense)
{
  Cbc_flush(model, FCColumns);
  model->solver_->setObjSense(sense);
}

COINLIBAPI void COINLINKAGE
Cbc_setRowLower(Cbc_Model *model, int index, double value)
{
  Cbc_flush(model, FCRows);
  OsiSolverInterface *solver = model->solver_;
  solver->setRowLower(index, value);
}

COINLIBAPI void COINLINKAGE
Cbc_setRowUpper(Cbc_Model *model, int index, double value)
{
  Cbc_flush(model, FCRows);
  OsiSolverInterface *solver = model->solver_;
  solver->setRowUpper(index, value);
}

/** @brief Constraint lower bounds 
  *
  * @param model problem object 
  * @return vector with lower bounds of constraints
  **/
COINLIBAPI const double *COINLINKAGE
Cbc_getRowLower(Cbc_Model *model)
{
  Cbc_flush(model, FCRows);
  OsiSolverInterface *solver = model->solver_;
  return solver->getRowLower();
}

/** @brief Constraint upper bounds 
  *
  * @param model problem object 
  * @return constraint upper bounds
  **/
COINLIBAPI const double *COINLINKAGE
Cbc_getRowUpper(Cbc_Model *model)
{
  Cbc_flush(model, FCRows);
  OsiSolverInterface *solver = model->solver_;
  return solver->getRowUpper();
}


COINLIBAPI const double *COINLINKAGE
Cbc_getRowActivity(Cbc_Model *model) {
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Information not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      return model->solver_->getRowActivity();
    case IntegerOptimization:
      return model->cbcModel_->getRowActivity();
  }

  return NULL;

}

COINLIBAPI const double *COINLINKAGE
Cbc_getColLower(Cbc_Model *model)
{
  Cbc_flush(model, FCColumns);
  return model->solver_->getColLower();
}

COINLIBAPI const double *COINLINKAGE
Cbc_getColUpper(Cbc_Model *model)
{
  Cbc_flush(model, FCColumns);
  return model->solver_->getColUpper();
}


COINLIBAPI double COINLINKAGE
Cbc_getBestPossibleObjValue(Cbc_Model *model) {
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Information not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      return model->solver_->getObjValue();
    case IntegerOptimization:
      return model->cbcModel_->getBestPossibleObjValue();
  }

  return COIN_DBL_MIN;
}

COINLIBAPI const double *COINLINKAGE
Cbc_getObjCoefficients(Cbc_Model *model)
{
  Cbc_flush(model, FCColumns);
  return model->solver_->getObjCoefficients();
}

COINLIBAPI void COINLINKAGE
Cbc_setObjCoeff(Cbc_Model *model, int index, double value)
{
  Cbc_flush( model, FCColumns );
  model->solver_->setObjCoeff( index, value );
}

COINLIBAPI void COINLINKAGE
Cbc_setColLower(Cbc_Model *model, int index, double value)
{
  Cbc_flush(model, FCColumns);
  model->solver_->setColLower( index, value );
}

COINLIBAPI void COINLINKAGE
Cbc_setColUpper(Cbc_Model *model, int index, double value)
{
  Cbc_flush(model, FCColumns);
  model->solver_->setColUpper( index, value );
}


COINLIBAPI const double *COINLINKAGE
Cbc_bestSolution(Cbc_Model *model)
{
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Information not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      return model->solver_->getColSolution();
    case IntegerOptimization:
      return model->cbcModel_->bestSolution();
  }

  return NULL;
}


COINLIBAPI int COINLINKAGE
Cbc_isInteger(Cbc_Model *model, int i)
{
  Cbc_flush(model);
  VALIDATE_COL_INDEX( i, model );
  return model->solver_->isInteger(i);
}


COINLIBAPI int COINLINKAGE
Cbc_getNodeCount(Cbc_Model *model) {
  switch (model->lastOptimization) {
    case ModelNotOptimized:
      fprintf( stderr, "Information not available, model was not optimized yet.\n");
      abort();
      break;
    case ContinuousOptimization:
      fprintf( stderr, "Information only available when optimizing integer models.\n");
      abort();
      break;
    case IntegerOptimization:
      return model->cbcModel_->getNodeCount();
  }

  return false;

}


/** Return a copy of this model */
COINLIBAPI Cbc_Model *COINLINKAGE
Cbc_clone(Cbc_Model *model)
{
  Cbc_flush(model);
  Cbc_Model *result = new Cbc_Model();

  result->solver_ = dynamic_cast<OsiClpSolverInterface *>(model->solver_->clone());

  if (model->cbcModel_)
    result->cbcModel_ = model->cbcModel_->clone(true);
  else
    result->cbcModel_ = NULL;

  result->relax_ = model->relax_;

  result->inc_callback = model->inc_callback;
  result->progr_callback = model->progr_callback;
  result->userCallBack = model->userCallBack;

  result->cut_callback = model->cut_callback;
  result->cutCBData = model->cutCBData;
  result->cutCBhowOften = model->cutCBhowOften;
  result->cutCBAtSol = model->cutCBAtSol;


  result->lastOptimization = model->lastOptimization;

  result->icAppData = model->icAppData;
  result->pgrAppData = model->pgrAppData;
  result->colNameIndex = NULL;
  result->rowNameIndex = NULL;
  if (model->colNameIndex)
    Cbc_storeNameIndexes(result, 1);

  Cbc_iniBuffer(result);

  if (model->iniSol) {
    result->iniSol = new vector<double>( model->iniSol->begin(), model->iniSol->end() );
    result->iniObj = model->iniObj;
  }
  else
  {
    result->iniSol = NULL;
    result->iniObj = COIN_DBL_MAX;
  }

  if (model->nColsMS) {
    result->nColsMS = model->nColsMS;
    result->colNamesMS = (char **) xmalloc( sizeof(char *)*model->nColsMS );
    result->charSpaceMS = model->charSpaceMS;
    result->colNamesMS[0] = (char *) xmalloc( result->charSpaceMS );
    memcpy( result->colNamesMS[0], model->colNamesMS[0], model->charSpaceMS );
    for ( int i=1 ; (i<model->nColsMS) ; ++i )
      result->colNamesMS[i] = result->colNamesMS[i-1] + strlen(result->colNamesMS[i-1]);
  }
  else
  {
    result->nColsMS = 0;
    result->colNamesMS = NULL;
    result->colValuesMS = NULL;
    result->charSpaceMS = 0;
  }

  if ( model->nSos ) {
    result->nSos = model->nSos;
    result->sosSize = model->sosSize;
    result->sosElSize = model->sosElSize;

    result->sosCap = result->nSos;
    result->sosElCap = model->sosElCap;

    result->sosRowStart = (int*) xmalloc( sizeof(int)*(result->sosCap+1) );
    result->sosEl = (int*) xmalloc( sizeof(int)*(result->sosElCap) );
    result->sosElWeight = (double*) xmalloc( sizeof(double)*(result->sosElCap) );

    memcpy( result->sosRowStart, model->sosRowStart, sizeof(int)*(result->nSos+1) );
    memcpy( result->sosEl, model->sosEl, sizeof(int)*(result->sosElSize) );
    memcpy( result->sosElWeight, model->sosElWeight, sizeof(double)*(result->sosElSize) );
  }

#ifdef CBC_THREAD
  pthread_mutex_init(&(result->cbcMutex), NULL);
#endif

  // copying parameters
  result->lp_method = model->lp_method;
  result->dualp = model->dualp;

  memcpy(result->int_param, model->int_param, sizeof(result->int_param));
  memcpy(result->dbl_param, model->dbl_param, sizeof(result->dbl_param));

  return result;
}

/** Set this the variable to be continuous */
COINLIBAPI void COINLINKAGE
Cbc_setContinuous(Cbc_Model *model, int iColumn)
{
  Cbc_flush(model, FCColumns);

  model->solver_->setContinuous(iColumn);
}

/** Set this the variable to be integer */
COINLIBAPI void COINLINKAGE
Cbc_setInteger(Cbc_Model *model, int iColumn)
{
  Cbc_flush(model, FCColumns);

  model->solver_->setInteger(iColumn);
}

/** Adds a new column */
COINLIBAPI void COINLINKAGE
Cbc_addCol(Cbc_Model *model, const char *name, double lb,
  double ub, double obj, char isInteger,
  int nz, int *rows, double *coefs)
{
  OsiSolverInterface *solver = model->solver_;

  if ( nz==0 )
  {
    Cbc_addColBuffer( model, name, lb, ub, obj, isInteger );
  }
  else
  {
    Cbc_flush(model, FCRows);
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
Cbc_addLazyConstraint(Cbc_Model *model, int nz,
  int *cols, double *coefs, char sense, double rhs)
{
  if (model->lazyConstrs==NULL)
    model->lazyConstrs = new CglStored();

  OsiRowCut orc;
  orc.setRow( nz, cols, coefs );

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

  model->lazyConstrs->addCut(orc);
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
  OsiSolverInterface *solver = model->solver_;

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
  OsiSolverInterface *solver = model->solver_;

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
  int newEl = rowStarts[numRows] - rowStarts[0];

  if (numRows + model->nSos > model->sosCap) {
    int prevCap = model->sosCap;
    if (prevCap) {
      model->sosCap = std::max( 2*model->sosCap, numRows + model->nSos);
      model->sosRowStart = (int *) xrealloc(model->sosRowStart, sizeof(int)*(model->sosCap+1) );
      model->sosType = (int *) xrealloc(model->sosRowStart, sizeof(int)*(model->sosCap) );
    } else {
      model->sosCap = max(1024, numRows);
      model->sosRowStart = (int *) xmalloc(sizeof(int)*(model->sosCap+1) );
      model->sosType = (int *) xmalloc(sizeof(int)*(model->sosCap) );
      model->sosRowStart[0] = 0;
      model->sosElCap = std::max(8192, newEl);
      model->sosEl = (int *) xmalloc( sizeof(int)*model->sosElCap );
      model->sosElWeight = (double *) xmalloc( sizeof(double)*model->sosElCap );
      model->sosElSize = 0;
      model->nSos = 0;
    }
  }

  for ( int i=0 ; i<numRows ; ++i )
    model->sosType[model->nSos+i] = type;

  if ( model->sosElSize + newEl > model->sosElCap ) {
    model->sosElCap = max( 2*model->sosElCap, newEl );
    model->sosEl = (int *) xrealloc( model->sosEl, sizeof(int)*model->sosElCap );
    model->sosElWeight  = (double *) xrealloc( model->sosElWeight, sizeof(double)*model->sosElCap );
  }

  memcpy( model->sosEl + model->sosElSize, colIndices, sizeof(int)*newEl );
  memcpy( model->sosElWeight + model->sosElSize, weights, sizeof(double)*newEl );

  for ( int i=0 ; (i<numRows) ; ++i ) {
    int size = rowStarts[i+1] - rowStarts[i];
    model->sosRowStart[model->nSos+1] = model->sosRowStart[model->nSos] + size;
    model->nSos++;
  }

  model->sosElSize += newEl;
}

COINLIBAPI void COINLINKAGE
Cbc_setMIPStart(Cbc_Model *model, int count, const char **colNames, const double colValues[])
{
  if (model->nColsMS) {
    if (model->colNamesMS) {
      free( model->colNamesMS[0]);
      free( model->colNamesMS);
    }
    free( model->colValuesMS );
  }

  int nameSpace = 0;
  for ( int i=0 ; (i<count) ; ++i )
    nameSpace += strlen(colNames[i]);
  nameSpace += count;

  model->colValuesMS = (double *) xmalloc( sizeof(double)*count );
  model->colNamesMS = (char **) xmalloc( sizeof(char*)*count );
  model->charSpaceMS = sizeof(char)*nameSpace;
  model->colNamesMS[0] = (char *) xmalloc( model->charSpaceMS );


  for ( int i=1 ; (i<count) ; ++i )
    model->colNamesMS[i] = model->colNamesMS[i-1] + 1 + strlen(colNames[i-1]);

  for ( int i=0 ; (i<count) ; ++i )
    strcpy( model->colNamesMS[i], colNames[i] );

  memcpy(model->colValuesMS, colValues, sizeof(double)*count );
}

COINLIBAPI void COINLINKAGE
Cbc_setMIPStartI(Cbc_Model *model, int count, const int colIdxs[], const double colValues[])
{
  OsiSolverInterface *solver = model->solver_;

  if (model->nColsMS) {
    if (model->colNamesMS) {
      free( model->colNamesMS[0]);
      free( model->colNamesMS);
    }
    free( model->colValuesMS );
  }

  int nameSpace = 0;
  for ( int i=0 ; (i<count) ; ++i )
    nameSpace += solver->getColName(colIdxs[i]).size();
  nameSpace += count;

  model->colValuesMS = (double *) xmalloc( sizeof(double)*count );
  model->colNamesMS = (char **) xmalloc( sizeof(char*)*count );
  model->charSpaceMS = sizeof(char)*nameSpace;
  model->colNamesMS[0] = (char *) xmalloc( model->charSpaceMS );

  for ( int i=1 ; (i<count) ; ++i )
    model->colNamesMS[i] = model->colNamesMS[i-1] + 1 + solver->getColName(colIdxs[i-1]).size();

  for ( int i=0 ; (i<count) ; ++i ) {
    strcpy( model->colNamesMS[i], solver->getColName(colIdxs[i]).c_str() );
    model->colValuesMS[i] = colValues[i];
  }

  memcpy(model->colValuesMS, colValues, sizeof(double)*count );
  model->nColsMS = count;
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

/** @brief integer tolerance */
COINLIBAPI double COINLINKAGE
Osi_getIntegerTolerance(void *osi)
{
  OsiSolverInterface *osiSolver = (OsiSolverInterface *) osi;
  return osiSolver->getIntegerTolerance();
}

COINLIBAPI double COINLINKAGE
Cbc_getAllowableGap(Cbc_Model* model)
{
  return model->dbl_param[DBL_PARAM_ALLOWABLE_GAP];
}

COINLIBAPI void COINLINKAGE
Cbc_setAllowableGap(Cbc_Model* model, double allowedGap)
{
  model->dbl_param[DBL_PARAM_ALLOWABLE_GAP] = allowedGap;
}

COINLIBAPI double COINLINKAGE
Cbc_getAllowableFractionGap(Cbc_Model* model)
{
  return model->dbl_param[DBL_PARAM_GAP_RATIO];
}

COINLIBAPI void COINLINKAGE
Cbc_setAllowableFractionGap(Cbc_Model* model, double allowedFracionGap)
{
  model->dbl_param[DBL_PARAM_GAP_RATIO] = allowedFracionGap;
}

/** returns the maximum number of nodes that can be explored in the search tree
 */
COINLIBAPI int COINLINKAGE
Cbc_getMaximumNodes(Cbc_Model *model)
{
  return model->int_param[INT_PARAM_MAX_NODES];
}

/** sets the maximum number of nodes that can be explored in the search tree
 */
COINLIBAPI void COINLINKAGE
Cbc_setMaximumNodes(Cbc_Model *model, int maxNodes)
{
  model->int_param[INT_PARAM_MAX_NODES] = maxNodes;
}

/** returns solution limit for the search process
 */
COINLIBAPI int COINLINKAGE
Cbc_getMaximumSolutions(Cbc_Model *model)
{
  return model->int_param[INT_PARAM_MAX_SOLS];
}

/** sets a solution limit as a stopping criterion
 */
COINLIBAPI void COINLINKAGE
Cbc_setMaximumSolutions(Cbc_Model *model, int maxSolutions)
{
  model->int_param[INT_PARAM_MAX_SOLS] = maxSolutions;
}

/** returns the current log leven
 */
COINLIBAPI int COINLINKAGE
Cbc_getLogLevel(Cbc_Model *model)
{
  return model->int_param[INT_PARAM_LOG_LEVEL];
}

/** sets the log level
 */
COINLIBAPI void COINLINKAGE
Cbc_setLogLevel(Cbc_Model *model, int logLevel)
{
  model->int_param[INT_PARAM_LOG_LEVEL] = logLevel;
}

/** gets the tolerance for infeasibility in the LP solver
 */
COINLIBAPI double COINLINKAGE
Cbc_getPrimalTolerance(Cbc_Model *model)
{
  return model->dbl_param[DBL_PARAM_PRIMAL_TOL];
}

/** sets the tolerance for infeasibility in the LP solver
 */
COINLIBAPI void COINLINKAGE
Cbc_setPrimalTolerance(Cbc_Model *model, double tol)
{
  model->dbl_param[DBL_PARAM_PRIMAL_TOL] = tol;
}

/** gets the tolerance for optimality in the LP solver
 */
COINLIBAPI double COINLINKAGE
Cbc_getDualTolerance(Cbc_Model *model)
{
  return model->dbl_param[DBL_PARAM_DUAL_TOL];
}

/** sets the tolerance for optimality in the LP solver
 */
COINLIBAPI void COINLINKAGE
Cbc_setDualTolerance(Cbc_Model *model, double tol)
{
  model->dbl_param[DBL_PARAM_DUAL_TOL] = tol;
}

COINLIBAPI double COINLINKAGE
Cbc_getCutoff(Cbc_Model* model)
{
  return model->dbl_param[DBL_PARAM_CUTOFF];
}

COINLIBAPI void COINLINKAGE
Cbc_setCutoff(Cbc_Model* model, double cutoff)
{
  model->dbl_param[DBL_PARAM_CUTOFF] = cutoff;
}

COINLIBAPI void COINLINKAGE
Cbc_setLPmethod(Cbc_Model *model, enum LPMethod lpm ) {
  model->lp_method = lpm;
}

COINLIBAPI void COINLINKAGE
Cbc_setDualPivot(Cbc_Model *model, enum DualPivot dp ) {
  model->dualp = dp;
}

COINLIBAPI double COINLINKAGE
Cbc_getMaximumSeconds(Cbc_Model *model)
{
  return model->dbl_param[DBL_PARAM_TIME_LIMIT];
}

COINLIBAPI void COINLINKAGE
Cbc_setMaximumSeconds(Cbc_Model *model, double maxSeconds)
{
  model->dbl_param[DBL_PARAM_TIME_LIMIT] = maxSeconds;
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

  NameIndex &rowNameIndex = *((NameIndex  *)model->rowNameIndex);
  NameIndex::iterator it = rowNameIndex.find(std::string(name));
  if (it == rowNameIndex.end())
    return -1;

  return it->second;
}

static char **to_char_vec( const vector< string > names )
{
    size_t spaceVec = (sizeof(char*)*names.size());
    size_t totLen = names.size(); // all \0
    for ( int i=0 ; (i<(int)names.size()) ; ++i )
        totLen += names[i].size();
    totLen *= sizeof(char);

    char **r = (char **)xmalloc(spaceVec+totLen);
    assert( r );
    r[0] = (char *)(r + names.size());
    for ( size_t i=1 ; (i<names.size()) ; ++i )
        r[i] = r[i-1] + names[i-1].size() + 1;

    for ( size_t i=0 ; (i<names.size()) ; ++i )
        strcpy(r[i], names[i].c_str());

    return r;
}

static void *xmalloc( const size_t size )
{
   void *result = malloc( size );
   if (!result)
   {
      fprintf(stderr, "No more memory available. Trying to allocate %zu bytes.", size);
      abort();
   }

   return result;
}


static void *xrealloc( void *ptr, const size_t newSize ) {
  void *res = realloc( ptr, newSize );

  if (!res) {
      fprintf(stderr, "No more memory available. Trying to allocate %zu bytes.", newSize);
      abort();
   }

  return res;
}

void Cbc_addAllSOS( Cbc_Model *model ) {
  if (model->nSos == 0)
    return;

  CbcModel *cbcModel = model->cbcModel_;

  vector< CbcObject *> objects;
  objects.reserve( model->nSos );
  for ( int i=0 ; i<model->nSos ; ++i ) {
    objects.push_back(
        new CbcSOS( 
            cbcModel, 
            model->sosRowStart[i+1]-model->sosRowStart[i],
            model->sosEl + model->sosRowStart[i],
            model->sosElWeight + model->sosRowStart[i],
            (int)objects.size(),
            model->sosType[i]
          ) 
        ); // add in objects
  }

  cbcModel->addObjects( (int) objects.size(), &objects[0] );

  for ( int i=0 ; i<model->nSos ; ++i ) 
    delete objects[i];
}

static void Cbc_addMS( Cbc_Model *model ) {
  if ( model->nColsMS == 0 )
    return;

  CbcModel *cbc = model->cbcModel_;

  cbc->setMIPStart( model->nColsMS, (const char **)model->colNamesMS, model->colValuesMS );
}

void Cbc_iniParams( Cbc_Model *model ) {
  model->lp_method = LPM_Auto;
  model->dualp = DP_Auto;
  memset(model->int_param, 0, sizeof(model->int_param) );
  for ( int i=0 ; (i<N_DBL_PARAMS) ; ++i )
    model->dbl_param[i] = 0.0;

  model->int_param[INT_PARAM_PERT_VALUE]       =       50;
  model->int_param[INT_PARAM_IDIOT]            =       -1;
  model->int_param[INT_PARAM_STRONG_BRANCHING] =        5;
  model->int_param[INT_PARAM_CUT_DEPTH]        =       -1;
  model->int_param[INT_PARAM_MAX_NODES]        =  INT_MAX;
  model->int_param[INT_PARAM_NUMBER_BEFORE]    =        5;
  model->int_param[INT_PARAM_FPUMP_ITS]        =       30;
  model->int_param[INT_PARAM_MAX_SOLS]         =       -1;
  model->int_param[INT_PARAM_CUT_PASS_IN_TREE] =        1;
  model->int_param[INT_PARAM_LOG_LEVEL]        =        1;
  model->int_param[INT_PARAM_MAX_SAVED_SOLS]   =       -1;
  model->int_param[INT_PARAM_MULTIPLE_ROOTS]   =        0;

  model->dbl_param[DBL_PARAM_PRIMAL_TOL]       =          1e-6;
  model->dbl_param[DBL_PARAM_DUAL_TOL]         =          1e-6;
  model->dbl_param[DBL_PARAM_ZERO_TOL]         =         1e-20;
  model->dbl_param[DBL_PARAM_INT_TOL]          =          1e-6;
  model->dbl_param[DBL_PARAM_PRESOLVE_TOL]     =          1e-8;
  model->dbl_param[DBL_PARAM_TIME_LIMIT]       =  COIN_DBL_MAX;
  model->dbl_param[DBL_PARAM_PSI]              =          -1.0;
  model->dbl_param[DBL_PARAM_CUTOFF]           =  COIN_DBL_MAX;
  model->dbl_param[DBL_PARAM_ALLOWABLE_GAP]    =         1e-10;
  model->dbl_param[DBL_PARAM_GAP_RATIO]        =        0.0001;
}

#if defined(__MWERKS__)
#pragma export off
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
