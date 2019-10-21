// * -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */

/**
  The original version of this file is taken from the MiniZinc source code
  www.minizinc.org
  */

/*
 *  Main authors:
 *     Gleb Belov <gleb.belov@monash.edu>
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifdef _MSC_VER 
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <chrono>
#include <stdexcept>

using namespace std;

#include "MIP_osicbc_wrap.hh"
#include "utils.hh"

#include "CbcSolver.hpp"
#include "CbcConfig.h"
#include "ClpConfig.h"
#include "CbcEventHandler.hpp"
#include "CglPreProcess.hpp"
#include "CglCutGenerator.hpp"
#include "CoinSignal.hpp"

#define WANT_SOLUTION

MIP_wrapper* createCBCWrapper(const MIP_wrapper::Options* opt) {
  return new MIP_osicbc_wrapper(opt);
}

string MIP_osicbc_wrapper::getDescription() {
  string v = "MIP wrapper for COIN-BC ";
  v += CBC_VERSION;                     // E.g., 2.9 stable or 2.9.7 latest release
  v += ",  using CLP ";
  v += CLP_VERSION;
  v += "  Compiled  " __DATE__ "  " __TIME__;
  return v;
}

string MIP_osicbc_wrapper::getVersion() {
  return string(CBC_VERSION)+"/"+string(CLP_VERSION);
}

string MIP_osicbc_wrapper::getId() {
  return "coin-bc";
}

string MIP_osicbc_wrapper::getName() {
  return "COIN-BC";
}

vector<string> MIP_osicbc_wrapper::getTags() {
  return {"mip","float","api","osicbc","coinbc","cbc"};
}

vector<string> MIP_osicbc_wrapper::getStdFlags() {
  return {"-a", "-p", "-s"};
}

void MIP_wrapper::Options::printHelp(ostream& os) {
  os
  << "\nMIP wrapper options:" << std::endl
  // -s                  print statistics
  //            << "  --readParam <file>  read OSICBC parameters from file
  //               << "--writeParam <file> write OSICBC parameters to file
  //               << "--tuneParam         instruct OSICBC to tune parameters instead of solving
  << "  --solver-flags \"args\"\n"
     "    command-line args passed to callCbc, e.g., \"-cuts off -preprocess off -passc 1\"." << std::endl
     //  \"-preprocess off\" recommended in 2.9.6
  << "  --writeModel <file>" << endl
  << "    write model to <file> (.mps)" << std::endl
  << "  -a, --all\n    print intermediate solutions for optimization problems\n"
     "    (not from FeasPump. Can be slow.)" << std::endl
   << "  -p <N>\n    use N threads, default: 1. CBC should be configured with --enable-cbc-parallel" << std::endl
//   << "--nomippresolve     disable MIP presolving   NOT IMPL" << std::endl
  << "  --solver-time-limit <N>\n    stop search after N milliseconds" << std::endl
//   << "--workmem <N>       maximal amount of RAM used, MB" << std::endl
//   << "--readParam <file>  read OSICBC parameters from file" << std::endl
//   << "--writeParam <file> write OSICBC parameters to file" << std::endl
//   << "--tuneParam         instruct OSICBC to tune parameters instead of solving   NOT IMPL"

  << "  --absGap <n>\n    absolute gap |primal-dual| to stop" << std::endl
  << "  --relGap <n>\n    relative gap |primal-dual|/<solver-dep> to stop. Default 1e-8, set <0 to use backend's default" << std::endl
  << "  --intTol <n>\n    integrality tolerance for a variable. Default 1e-8" << std::endl
//   << "--objDiff <n>       objective function discretization. Default 1.0" << std::endl

  << std::endl;
}

bool MIP_wrapper::Options::processOption(int& i, std::vector<std::string>& argv) {
  MiniZinc::CLOParser cop( i, argv );
  string sTmp;
  if ( string(argv[i])=="-a"
      || string(argv[i])=="--all"
      || string(argv[i])=="--all-solutions" ) {
    flag_all_solutions = true;
  } else if (string(argv[i])=="-f") {
//     std::cerr << "  Flag -f: ignoring fixed strategy anyway." << std::endl;
  } else if ( cop.get( "--writeModel", &sExportModel ) ) {
  } else if ( cop.get( "-p", &nThreads ) ) {
  } else if ( cop.get( "--solver-time-limit", &nTimeout ) ) {
  } else if ( cop.get( "--workmem", &nWorkMemLimit ) ) {
  } else if ( cop.get( "--readParam", &sReadParams ) ) {
  } else if ( cop.get( "--writeParam", &sWriteParams ) ) {
  } else if ( cop.get( "--solver-flags",
    &sTmp ) ) {
    cmdOptions += ' ';
    cmdOptions += sTmp;
  } else if ( cop.get( "--absGap", &absGap ) ) {
  } else if ( cop.get( "--relGap", &relGap ) ) {
  } else if ( cop.get( "--intTol", &intTol ) ) {
//   } else if ( cop.get( "--objDiff", &objDiff ) ) {
  } else
    return false;
  return true;
}

void MIP_osicbc_wrapper::wrap_assert(bool cond, string msg, bool fTerm)
{
   if ( !cond ) {
//       strcpy(osicbc_buffer, "[NO ERROR STRING GIVEN]");
//       CBCgeterrorstring (env, status, osicbc_buffer);
      string msgAll = ("  MIP_osicbc_wrapper runtime error:  " + msg + "  " + osicbc_buffer);
      cerr << msgAll << endl;
      if (fTerm) {
        cerr << "TERMINATING." << endl;
        throw runtime_error(msgAll);
      }
   }
}

void MIP_osicbc_wrapper::doAddVars
  (size_t n, double* obj, double* lb, double* ub, MIP_wrapper::VarType* vt, string *names)
{
  /// Convert var types:
//   vector<char> ctype(n);
//   vector<char*> pcNames(n);
  CoinPackedVector cpv;
  vector<CoinPackedVectorBase*> pCpv(n, &cpv);
  osi.addCols(n, pCpv.data(), lb, ub, obj);   // setting integer & names later
//   status = CBCnewcols (env, lp, n, obj, lb, ub, &ctype[0], &pcNames[0]);
//   wrap_assert( !status,  "Failed to declare variables." );
}

void MIP_osicbc_wrapper::addRow
  (int nnz, const int* rmatind, const double* rmatval, MIP_wrapper::LinConType sense,
   double rhs, int mask, const string& rowName)
{
  /// Convert var types:
  double rlb=rhs, rub=rhs;
    switch (sense) {
      case LQ:
        rlb = -osi.getInfinity();
        break;
      case EQ:
        break;
      case GQ:
        rub = osi.getInfinity();
        break;
      default:
        throw runtime_error("  MIP_wrapper: unknown constraint type");
    }
  // ignoring mask for now.  TODO
  // 1-by-1 too slow:
//   try {
//     CoinPackedVector cpv(nnz, rmatind, rmatval);
//     osi.addRow(cpv, rlb, rub);
//   } catch (const CoinError& err) {
//     cerr << "  COIN-OR Error: " << err.message() << endl;
//     throw runtime_error(err.message());
//   }
  /// Segfault:
//   rowStarts.push_back(columns.size());
//   columns.insert(columns.end(), rmatind, rmatind + nnz);
//   element.insert(element.end(), rmatval, rmatval + nnz);
  rows.push_back(CoinPackedVector(nnz, rmatind, rmatval));
  rowlb.push_back(rlb);
  rowub.push_back(rub);
}


bool MIP_osicbc_wrapper::addWarmStart( const std::vector<VarId>& vars, const std::vector<double> vals ) {
  assert( vars.size()==vals.size() );
  static_assert( sizeof(VarId)==sizeof(int), "VarId should be (u)int currently" );
  for (int i=0; i<vars.size(); ++i)
    warmstart[vars[i]] = vals[i];
  return true;
}

/// SolutionCallback ------------------------------------------------------------------------
/// OSICBC ensures thread-safety?? TODO
/// Event handling copied from examples/interrupt.cpp, Cbc 2.9.8 rev 2272

/************************************************************************

This main program shows how to take advantage of the standalone cbc in your program,
while still making major modifications.
This is very like driver4 but allows interrupts in clp for faster stopping
It would be up to user to clean up output as stopping in Clp seems to
give correct results but can think it is stopping in an odd way.
To make cleaner would need more events defined (in Cbc AND Clp)
First it reads in an integer model from an mps file
Then it initializes the integer model with cbc defaults
Then it calls CbcMain1 passing all parameters apart from first but with callBack to modify stuff
Finally it prints solution

************************************************************************/
/* Meaning of whereFrom:
   1 after initial solve by dualsimplex etc
   2 after preprocessing
   3 just before branchAndBound (so user can override)
   4 just after branchAndBound (before postprocessing)
   5 after postprocessing
*/
/* Meaning of model status is as normal
   status
      -1 before branchAndBound
      0 finished - check isProvenOptimal or isProvenInfeasible to see if solution found
      (or check value of best solution)
      1 stopped - on maxnodes, maxsols, maxtime
      2 difficulties so run was abandoned
      (5 event user programmed event occurred) 

      cbc secondary status of problem
        -1 unset (status_ will also be -1)
  0 search completed with solution
  1 linear relaxation not feasible (or worse than cutoff)
  2 stopped on gap
  3 stopped on nodes
  4 stopped on time
  5 stopped on user event
  6 stopped on solutions
  7 linear relaxation unbounded

   but initially check if status is 0 and secondary status is 1 -> infeasible
   or you can check solver status.
*/
/* Return non-zero to return quickly */   
//static int callBack(CbcModel * model, int whereFrom)
//{
//  int returnCode=0;
//  switch (whereFrom) {
//  case 1:
//  case 2:
//    if (!model->status()&&model->secondaryStatus())
//      returnCode=1;
//    break;
//  case 3:
//    {
//      //CbcCompareUser compare;
//      //model->setNodeComparison(compare);
//    }
//    break;
//  case 4:
//    // If not good enough could skip postprocessing
//    break;
//  case 5:
//    break;
//  default:
//    abort();
//  }
//  return returnCode;
//}
static int cancelAsap=0;
/*
  0 - not yet in Cbc
  1 - in Cbc with new signal handler
  2 - ending Cbc
*/
static int statusOfCbc=0;
static CoinSighandler_t saveSignal = static_cast<CoinSighandler_t> (0);

extern "C" {
     static void
#if defined(_MSC_VER)
     __cdecl
#endif // _MSC_VER
     signal_handler(int /*whichSignal*/)
     {
       cancelAsap=3;
       return;
     }
}
/** This is so user can trap events and do useful stuff.  

    CbcModel model_ is available as well as anything else you care 
    to pass in
*/

struct EventUserInfo {
  MIP_wrapper::CBUserInfo* pCbui=0;
  CglPreProcess* pPP=0;
};

extern CglPreProcess * cbcPreProcessPointer;
class MyEventHandler3 : public CbcEventHandler {
  
public:
  /**@name Overrides */
  //@{
  virtual CbcAction event(CbcEvent whichEvent);
  //@}

  /**@name Constructors, destructor etc*/
  //@{
  /** Default constructor. */
  MyEventHandler3(EventUserInfo& u_);
  /// Constructor with pointer to model (redundant as setEventHandler does)
  MyEventHandler3(CbcModel * model, EventUserInfo& u_);
  /** Destructor */
  virtual ~MyEventHandler3();
  /** The copy constructor. */
  MyEventHandler3(const MyEventHandler3 & rhs);
  /// Assignment
  MyEventHandler3& operator=(const MyEventHandler3 & rhs);
  /// Clone
  virtual CbcEventHandler * clone() const ;
  //@}
   
    
protected:
  // data goes here
  EventUserInfo ui;
  double bestSolutionValue_ = DBL_MAX;  // always min
};
//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
MyEventHandler3::MyEventHandler3 (EventUserInfo& u_) 
  : CbcEventHandler(), ui(u_)
{
  assert(0);
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
MyEventHandler3::MyEventHandler3 (const MyEventHandler3 & rhs) 
: CbcEventHandler(rhs), ui(rhs.ui)
{  
  
}

// Constructor with pointer to model
MyEventHandler3::MyEventHandler3(CbcModel * model, EventUserInfo& u_)
  : CbcEventHandler(model), ui(u_)
{
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
MyEventHandler3::~MyEventHandler3 ()
{
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
MyEventHandler3 &
MyEventHandler3::operator=(const MyEventHandler3& rhs)
{
  if (this != &rhs) {
    CbcEventHandler::operator=(rhs);
  }
  ui = rhs.ui;
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CbcEventHandler * MyEventHandler3::clone() const
{
  return new MyEventHandler3(*this);
}

CbcEventHandler::CbcAction 
MyEventHandler3::event(CbcEvent whichEvent)
{
  if(!statusOfCbc) {
    // override signal handler
    // register signal handler
    saveSignal = signal(SIGINT, signal_handler);
    statusOfCbc=1;
  }
  if ( (cancelAsap&2)!=0 ) {
    //     printf("Cbc got cancel\n");
    // switch off Clp cancel
    cancelAsap &= 2;
    return stop;
  }
  // If in sub tree carry on
  if (!model_->parentModel()) {
    if (whichEvent==endSearch&&statusOfCbc==1) {
      // switch off cancel
      cancelAsap=0;
      // restore signal handler
      signal(SIGINT, saveSignal);
      statusOfCbc=2;
    }
    if (whichEvent==solution||whichEvent==heuristicSolution) {
      // John Forrest  27.2.16:
      // check not duplicate
      if (model_->getObjValue()<bestSolutionValue_) {
        bestSolutionValue_ = model_->getObjValue();
        // If preprocessing was done solution will be to processed model
        //       int numberColumns = model_->getNumCols();
        const double * bestSolution = model_->bestSolution();
        assert (bestSolution);
        //       printf("value of solution is %g\n",model_->getObjValue());
        
        // Trying to obtain solution for the original model:
        assert( model_ && model_->solver() );
        // double objOffset=0;
        // model_->solver()->getDblParam(OsiObjOffset, objOffset);
        double objVal = (model_->getObjValue() );  //- objOffset);   John Forrest suggested to remove, 17.11.17
        double bestBnd = (model_->getBestPossibleObjValue() );  //- objOffset);
        if ( 0!=cbcPreProcessPointer ) {
          if ( OsiSolverInterface* cbcPreOrig = cbcPreProcessPointer->originalModel() ) {
            objVal *= cbcPreOrig->getObjSense();
            bestBnd *= cbcPreOrig->getObjSense();
          }
        } else {
          objVal *= model_->getObjSense();
          bestBnd *= model_->getObjSense();
        }
        OsiSolverInterface* origModel=0;
        if ( 0!=cbcPreProcessPointer && 0!=model_->continuousSolver() ) {
          OsiSolverInterface * solver = (model_->continuousSolver()->clone());
          //       ? model_->continuousSolver()->clone()
          //       : model_->continuousSolver()->clone();
          int numberColumns = solver->getNumCols();
          for (int i=0;i<numberColumns;i++) {
            if (solver->isInteger(i)) {
              solver->setColLower(i,bestSolution[i]);
              solver->setColUpper(i,bestSolution[i]);
            }
          }
          solver->resolve();
          cbcPreProcessPointer->postProcess( *solver, false );
          delete solver;
          origModel = cbcPreProcessPointer->originalModel();
          ui.pCbui->pOutput->x = origModel->getColSolution();
        } else {
          origModel = model_->solver();
          ui.pCbui->pOutput->x = bestSolution;
        }
        if ( ui.pCbui->fVerb )
          cerr
          << " % OBJ VAL RAW: " << model_->getObjValue()
          << "  OBJ VAL ORIG(?): " << objVal
          << " % BND RAW: " << model_->getBestPossibleObjValue()
          << "  BND ORIG(?): " << bestBnd
          //         << "  &prepro: " << cbcPreProcessPointer
          //         << "  &model_._solver(): " << model_->solver()
          << "  orig NCols: " << ui.pCbui->pOutput->nCols
          << "  prepro NCols:  " << model_->getNumCols()
          ;
        assert( origModel->getNumCols() == ui.pCbui->pOutput->nCols );
        if ( ui.pCbui->fVerb ) {
          if ( ui.pCbui->pOutput->nObjVarIndex>=0 )
            cerr
            << "  objVAR: " << ui.pCbui->pOutput->x[ui.pCbui->pOutput->nObjVarIndex];
          cerr << endl;
        }
        ui.pCbui->pOutput->objVal = objVal;
        //         origModel->getObjValue();
        ui.pCbui->pOutput->status = MIP_wrapper::SAT;
        ui.pCbui->pOutput->statusName = "feasible from a callback";
        ui.pCbui->pOutput->bestBound = bestBnd;
        ui.pCbui->pOutput->dWallTime = std::chrono::duration<double>(
                                                                     std::chrono::steady_clock::now() - ui.pCbui->pOutput->dWallTime0).count();
        ui.pCbui->pOutput->dCPUTime = model_->getCurrentSeconds();
        ui.pCbui->pOutput->nNodes = model_->getNodeCount();
        ui.pCbui->pOutput->nOpenNodes = -1; // model_->getNodeCount2();
        
        /// Call the user function:
        if (ui.pCbui->solcbfn) {
          (*(ui.pCbui->solcbfn))(*(ui.pCbui->pOutput), ui.pCbui->psi);
          ui.pCbui->printed = true;
        }
        return noAction; // carry on
      } else {
        return noAction; // carry on
      }
    } else {
      return noAction;
    }
  } else {
    return noAction; // carry on
  }
}

/** This is so user can trap events and do useful stuff.  

    ClpSimplex model_ is available as well as anything else you care 
    to pass in
*/
class MyEventHandler4 : public ClpEventHandler {
  
public:
  /**@name Overrides */
  //@{
  virtual int event(Event whichEvent);
  //@}

  /**@name Constructors, destructor etc*/
  //@{
  /** Default constructor. */
  MyEventHandler4();
  /// Constructor with pointer to model (redundant as setEventHandler does)
  MyEventHandler4(ClpSimplex * model);
  /** Destructor */
  virtual ~MyEventHandler4();
  /** The copy constructor. */
  MyEventHandler4(const MyEventHandler4 & rhs);
  /// Assignment
  MyEventHandler4& operator=(const MyEventHandler4 & rhs);
  /// Clone
  virtual ClpEventHandler * clone() const ;
  //@}
   
    
protected:
  // data goes here
};
//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
MyEventHandler4::MyEventHandler4 () 
  : ClpEventHandler()
{
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
MyEventHandler4::MyEventHandler4 (const MyEventHandler4 & rhs) 
: ClpEventHandler(rhs)
{  
}

// Constructor with pointer to model
MyEventHandler4::MyEventHandler4(ClpSimplex * model)
  : ClpEventHandler(model)
{
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
MyEventHandler4::~MyEventHandler4 ()
{
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
MyEventHandler4 &
MyEventHandler4::operator=(const MyEventHandler4& rhs)
{
  if (this != &rhs) {
    ClpEventHandler::operator=(rhs);
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpEventHandler * MyEventHandler4::clone() const
{
  return new MyEventHandler4(*this);
}

int
MyEventHandler4::event(Event whichEvent)
{
  if ( (cancelAsap&1)!=0 ) {
//     printf("Clp got cancel\n");
    return 5;
  } else {
    return -1;
  }
}
// end SolutionCallback ---------------------------------------------------------------------

MIP_osicbc_wrapper::Status MIP_osicbc_wrapper::convertStatus(CbcModel *pModel)
{
  Status s = Status::UNKNOWN;
   /* Converting the status. */
   if (pModel->isProvenOptimal()) {
       s = Status::OPT;
       output.statusName = "Optimal";
//        wrap_assert(osi., "Optimality reported but pool empty?", false);
   } else if (pModel->isProvenInfeasible()) {
       s = Status::UNSAT;
       output.statusName = "Infeasible";
   } else if (pModel->isProvenDualInfeasible()) {
       s = Status::UNBND;
       output.statusName = "Dual infeasible";
//        s = Status::UNSATorUNBND;
   } else if   // wrong: (pModel->getColSolution())
     (fabs(pModel->getObjValue()) < 1e50)
   {
       s = Status::SAT;
       output.statusName = "Feasible";
   } else if (pModel->isAbandoned()) {     // AFTER feas-ty
       s = Status::__ERROR;
       output.statusName = "Abandoned";
   } else {
     s = Status::UNKNOWN;
     output.statusName = "Unknown";
   }
   return s;
}

MIP_osicbc_wrapper::Status MIP_osicbc_wrapper::convertStatus()
{
  Status s = Status::UNKNOWN;
   /* Converting the status. */
   if (osi.isProvenOptimal()) {
       s = Status::OPT;
       output.statusName = "Optimal";
//        wrap_assert(osi., "Optimality reported but pool empty?", false);
   } else if (osi.isProvenPrimalInfeasible()) {
       s = Status::UNSAT;
       output.statusName = "Infeasible";
   } else if (osi.isProvenDualInfeasible()) {
       s = Status::UNBND;
       output.statusName = "Dual infeasible";
//        s = Status::UNSATorUNBND;
   } else if (osi.isAbandoned()) {
       s = Status::__ERROR;
       output.statusName = "Abandoned";
   } else if   // wrong: (pModel->getColSolution())
     (fabs(osi.getObjValue()) < osi.getInfinity())
   {
       s = Status::SAT;
       output.statusName = "Feasible";
       cout << " getSolverObjValue(as minim) == " << osi.getObjValue() << endl;
   }
   else {
     s = Status::UNKNOWN;
     output.statusName = "Unknown";
   }
   return s;
}


void MIP_osicbc_wrapper::solve() {  // Move into ancestor?
  if ( options->flag_all_solutions && 0==nProbType )
    cerr << "WARNING. --all-solutions for SAT problems not implemented." << endl;
  try {
    /// Not using CoinPackedMatrix any more, so need to add all constraints at once:
    /// But this gives segf:
//     osi.addRows(rowStarts.size(), rowStarts.data(),
//                 columns.data(), element.data(), rowlb.data(), rowub.data());
    /// So:
    MIP_wrapper::addPhase1Vars();         // only now
    if (fVerbose) {
      cerr << getDescription() << std::endl;
      cerr << "  MIP_osicbc_wrapper: adding constraints physically..." << flush;
    }
    vector<CoinPackedVectorBase*> pRows(rowlb.size());
    for (int i=0; i<rowlb.size(); ++i)
      pRows[i] = &rows[i];
    osi.addRows(rowlb.size(), pRows.data(), rowlb.data(), rowub.data());
//     rowStarts.clear();
//     columns.clear();
//     element.clear();
    pRows.clear();
    rows.clear();
    rowlb.clear();
    rowub.clear();
    if (fVerbose)
      cerr << " done." << endl;
  /////////////// Last-minute solver options //////////////////
//       osi->loadProblem(*matrix, 
    {
      std::vector<VarId> integer_vars;
      for(unsigned int i=0; i<colObj.size(); i++) {
        if(REAL != colTypes[i]
//           && is_used[i]
        ) {
          integer_vars.push_back(i);
        }
      }
      osi.setInteger(integer_vars.data(), integer_vars.size());
    }
    if(options->sExportModel.size()) {
      // Not implemented for OsiClp:
//       osi.setColNames(colNames, 0, colObj.size(), 0);
      vector<const char*> colN(colObj.size());
      for (int j=0; j<colNames.size(); ++j)
        colN[j] = colNames[j].c_str();
      osi.writeMpsNative(options->sExportModel.c_str(), 0, colN.data());
    }
    
    // Tell solver to return fast if presolve or initial solve infeasible
    osi.getModelPtr()->setMoreSpecialOptions(3);
    // allow Clp to handle interrupts
    MyEventHandler4 clpEventHandler;
    osi.getModelPtr()->passInEventHandler(&clpEventHandler);

  /* switch on/off output to the screen */
    class NullCoinMessageHandler : public CoinMessageHandler {
      int print() {
        return 0;
      }
      void checkSeverity() {
      }
    } nullHandler;

//     CbcSolver control(osi);
//     // initialize   ???????
//     control.fillValuesInSolver();  
//     CbcModel * pModel = control.model();
    if ( fVerbose )
      cerr << " Model creation..." << endl;
    
// #define __USE_CbcSolver__  -- not linked rev2274
    /// FOR WARMSTART
    for (const auto& vv: warmstart) {
      osi.setColName(vv.first, colNames[vv.first]);
    }
#ifdef __USE_CbcSolver__
    CbcSolver control(osi);
    // initialize
    control.fillValuesInSolver();
    CbcModel& model = *control.model();
#else
//    CbcModel model(osi);

    /// Changing so we can just assign solver
    CbcModel model;
    {
      OsiSolverInterface *pOsi = &osi;
      model.assignSolver(pOsi, false);
      model.setModelOwnsSolver(false);
    }

#endif
//     CbcSolver control(osi);
//     control.solve();
    if ( options->absGap>=0.0 )
      model.setAllowableGap( options->absGap );
    if ( options->relGap>=0.0 )
      model.setAllowableFractionGap( options->relGap );
    if ( options->intTol>=0.0 )
      model.setIntegerTolerance( options->intTol );
//     model.setCutoffIncrement( objDiff );

    /// WARMSTART
    {
      std::vector< std::pair< std::string, double > > mipstart;
      for (const auto& vv: warmstart) {
        mipstart.push_back( std::make_pair(colNames[vv.first], vv.second) );
      }
      warmstart.clear();
      model.setMIPStart(mipstart);
    }
    
    CoinMessageHandler msgStderr(stderr);

    class StderrCoinMessageHandler : public CoinMessageHandler {
      int print() {
        cerr << messageBuffer_ << endl;
		    return 0;
      }
      void checkSeverity() {
      }
    } stderrHandler;
    
    if(fVerbose) {
//        osi.messageHandler()->setLogLevel(1);
//        osi.getModelPtr()->setLogLevel(1);
//        osi.getRealSolverPtr()->messageHandler()->setLogLevel(0);
      // DOES NOT WORK:                                                     TODO
//        model.passInMessageHandler( &stderrHandler );
       msgStderr.setLogLevel( 0, 1 );
       model.passInMessageHandler( &msgStderr );
//        model.setLogLevel(1);
//        model.solver()->messageHandler()->setLogLevel(0);
    } else {
      model.passInMessageHandler(&nullHandler);
      model.messageHandler()->setLogLevel(0);
      model.setLogLevel(0);
      model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
//       osi.passInMessageHandler(&nullHandler);
//       osi.messageHandler()->setLogLevel(0);
//       osi.setHintParam(OsiDoReducePrint, true, OsiHintTry);
    }

    if(options->nTimeout != 0) {
//       osi.setMaximumSeconds(nTimeout);
      model.setMaximumSeconds(static_cast<double>(options->nTimeout)/1000.0);
    }

   /// TODO
//     if(all_solutions && obj.getImpl()) {
//       IloNum lastObjVal = (obj.getSense() == IloObjective::Minimize ) ?
//       _iloosicbc->use(SolutionCallback(_iloenv, lastObjVal, *this));
      // Turn off OSICBC logging

    
   /// Solution callback
   output.nCols = colObj.size();
//    x.resize(output.nCols);
//    output.x = &x[0];

   if (options->flag_all_solutions && cbui.solcbfn) {
     // Event handler. Should be after CbcMain0()?
     EventUserInfo ui;
     ui.pCbui = &cbui;
//      ui.pPP = 0;
     MyEventHandler3 eventHandler(&model, ui);
     model.passInEventHandler(&eventHandler);
   }

   /// Cuts needed
   if ( cbui.cutcbfn ) {
     /// This class is passed to CBC to organize cut callbacks
     /// We need original solutions here (combinatorial cuts)
     class CutCallback : public CglCutGenerator {
       MIP_wrapper::CBUserInfo& cbui;
     public:
       CutCallback(MIP_wrapper::CBUserInfo& ui) : cbui(ui) { }
       CglCutGenerator* clone() const override { return new CutCallback(cbui); }
       /// Make sure this overrides but we might need to compile this with old CBC as well
       bool needsOriginalModel() const /*override*/ { return true; }
       void generateCuts(const OsiSolverInterface &si, OsiCuts &cs,
                         const CglTreeInfo info = CglTreeInfo()) override {
         cbui.pOutput->nCols = si.getNumCols();
         MZN_ASSERT_HARD_MSG(cbui.pOutput->nCols == ((MIP_wrapper*)(cbui.wrapper))->colNames.size(),
                    "CBC cut callback: current model is different? Ncols=" << cbui.pOutput->nCols
                             << ", originally " << ((MIP_wrapper*)(cbui.wrapper))->colNames.size()
                             << ". If you have an old version of CBC, to use combinatorial cuts"
                                " run with --cbcArgs '-preprocess off'" );
         cbui.pOutput->x = si.getColSolution(); // change the pointer?
         MIP_wrapper::CutInput cuts;
         cbui.cutcbfn( *cbui.pOutput, cuts, cbui.psi, info.options&128 );
                                     // options&128 is nonzero when we have an integer candidate
         for (const auto& cut: cuts) {                  // Convert cut sense
           OsiRowCut rc;
           switch (cut.sense) {
           case LQ:
             rc.setUb(cut.rhs);
             break;
           case GQ:
             rc.setLb(cut.rhs);
             break;
           default:
             assert(EQ==cut.sense);
             rc.setLb(cut.rhs);
             rc.setUb(cut.rhs);
           }
           rc.setRow(cut.rmatind.size(), cut.rmatind.data(), cut.rmatval.data());
           cs.insertIfNotDuplicate(rc);
         }
       }
     };
     CutCallback ccb(cbui);
     model.addCutGenerator(&ccb, 1,                        // how often
                           "MZN_cuts",
                           true,                            // "normal" (what's this?)
                           true);                           // true: also at solution
                                                // TODO @solution only for lazy cuts (see cut.mask)
   }

   cbc_cmdOptions = options->cmdOptions;
   if ( 1<options->nThreads ) {
    ostringstream oss;
    oss << " -threads " << options->nThreads;
    cbc_cmdOptions += oss.str();
   }
   cbc_cmdOptions += " -solve";
   cbc_cmdOptions += " -quit";

   cbui.pOutput->dWallTime0 = output.dWallTime0 =
     std::chrono::steady_clock::now();
   output.dCPUTime = clock();

   /* OLD: Optimize the problem and obtain solution. */
//       model.branchAndBound();
//       osi.branchAndBound();

      /// TAKEN FORM DRIVER3.CPP, seems to use most features:
//      CbcMain0(model);
//      CbcCbcParamUtils::setCbcModelDefaults(model) ;
//       const char * argv2[]={"mzn-cbc","-solve","-quit"};
//        CbcMain1(3,argv2,model);
#ifdef __USE_CbcSolver__
  if (fVerbose)
    cerr << "  Calling control.solve() with options '" << options->cbc_cmdOptions << "'..." << endl;
  control.solve (options->cbc_cmdOptions.c_str(), 1);
#else
#define __USE_callCbc1__
#ifdef __USE_callCbc1__
    if (fVerbose)
      cerr << "  Calling callCbc with options '" << cbc_cmdOptions << "'..." << endl;
    callCbc(cbc_cmdOptions, model);
//     callCbc1(cbc_cmdOptions, model, callBack);
    // What is callBack() for?    TODO
#else
  CbcMain0(model);
  // should be here?
//   // Event handler
//    EventUserInfo ui;
//    MyEventHandler3 eventHandler( &model, ui );
//    model.passInEventHandler(&eventHandler);
  /* Now go into code for standalone solver
     Could copy arguments and add -quit at end to be safe
     but this will do
  */
  vector<string> argvS;
  MiniZinc::split(cbc_cmdOptions, argvS);
  vector<const char*> argv;
  MiniZinc::vecString2vecPChar(argvS, argv);
  if (fVerbose)
    cerr << "  Calling CbcMain1 with options '" << cbc_cmdOptions << "'..." << endl;
  CbcMain1(argv.size(),argv.data(),model,callBack);
#endif
#endif
  
    output.dWallTime = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - output.dWallTime0).count();
    output.dCPUTime = (clock() - output.dCPUTime) / CLOCKS_PER_SEC;

    output.status = convertStatus(&model);
  //    output.status = convertStatus();

    /// Continuing to fill the output object:
    if (Status::OPT == output.status || Status::SAT ==output.status) {
      output.objVal = model.getObjValue();
  //      output.objVal = osi.getObjValue();

        /* The size of the problem should be obtained by asking OSICBC what
            the actual size is, rather than using what was passed to CBCcopylp.
            cur_numrows and cur_numcols store the current number of rows and
            columns, respectively.  */   // ?????????????? TODO

        int cur_numcols = model.getNumCols ();
  //       int cur_numcols = osi.getNumCols ();
        assert(cur_numcols == colObj.size());
        
        wrap_assert(model.getColSolution(), "Failed to get variable values.");
        x.assign( model.getColSolution(), model.getColSolution() + cur_numcols ); // ColSolution();
        output.x = x.data();
  //       output.x = osi.getColSolution();
      if (cbui.solcbfn && (!options->flag_all_solutions || !cbui.printed)) {
        cbui.solcbfn(output, cbui.psi);
      }
    }
    output.bestBound = model.getBestPossibleObjValue();
  //    output.bestBound = -1;
    output.nNodes = model.getNodeCount();
  //    output.nNodes = osi.getNodeCount();
    output.nOpenNodes = -1;
    
  } catch (CoinError& err) {
    err.print(true);
  }
}

void MIP_osicbc_wrapper::setObjSense(int s)
{
  osi.setObjSense(-s);
}

/*

try the following for example:

CbcMain0(model); 
const char * argv2[]={"driver4","-cuts","off" ,"-preprocess","off","-passc","1","-solve","-quit"};
CbcMain1(9,argv2,model);

you can add any feature you want to argv2 ...

if you want to add cuts yourself, or heuristics, do the following:

  OsiSolverInterface *solver2 = osi;
  CglPreProcess *process = new CglPreProcess;
  solver2 = process->preProcess(*solver,false,2);

    CbcModel model1(*solver2);
  
    model1.initialSolve();
    
  //==============================================
 
  CglProbing generator1;
  generator1.setUsingObjective(true);
  generator1.setMaxPass(1);
  generator1.setMaxPassRoot(5);
  generator1.setMaxProbe(10);
  generator1.setMaxProbeRoot(1000);
  generator1.setMaxLook(50);
  generator1.setMaxLookRoot(500);
  generator1.setMaxElements(200);
  generator1.setRowCuts(3);

  CglGomory generator2;
  generator2.setLimit(300);

  CglKnapsackCover generator3;

  CglRedSplit generator4;
  generator4.setLimit(200);

  CglClique generator5;
  generator5.setStarCliqueReport(false);
  generator5.setRowCliqueReport(false);

  CglMixedIntegerRounding2 mixedGen;
  CglFlowCover flowGen;
  
  CglGMI cut1;
  CglMixedIntegerRounding2 cut2;
  CglOddHole cut3;
  CglSimpleRounding cut4;
  CglResidualCapacity cut5;
  CglTwomir cut6;
  CglZeroHalf cut7;
  
  model1.addCutGenerator(&generator1,-1,"Probing");
  model1.addCutGenerator(&generator2,-1,"Gomory");
  model1.addCutGenerator(&generator3,-1,"Knapsack");
  model1.addCutGenerator(&generator4,-1,"RedSplit");
  model1.addCutGenerator(&generator5,-1,"Clique");
  model1.addCutGenerator(&flowGen,-1,"FlowCover");
  model1.addCutGenerator(&mixedGen,-1,"MixedIntegerRounding");
  model1.addCutGenerator(&cut1,-1,"GMI");
  model1.addCutGenerator(&cut2,-1,"MixedIntegerRounding2");
  model1.addCutGenerator(&cut3,-1,"OddHole");
  model1.addCutGenerator(&cut4,-1,"SimpleRounding");
  model1.addCutGenerator(&cut5,-1,"ResidualCapacity");
  model1.addCutGenerator(&cut6,-1,"Twomir");
  model1.addCutGenerator(&cut7,-1,"ZeroHalf");
  
 

  CbcRounding heuristic1(model1);
  CbcHeuristicLocal heuristic2(model1);

  
 model1.addHeuristic(&heuristic1);
 model1.addHeuristic(&heuristic2);

 
  
  
    model1.setMaximumCutPassesAtRoot(50); 
    model1.setMaximumCutPasses(1000);
  
  
  
  model1.branchAndBound();
  
  
  OsiSolverInterface * solver3;
  
  process->postProcess(*model1.solver());
  
  solver3 = solver;
  
 or, use the default strategy:

CbcStrategyDefault strategy(5);
model1.setStrategy(strategy);
 

  
  
  

On Sun, Oct 11, 2015 at 8:38 PM, Gleb Belov <gleb.belov@monash.edu> wrote:

    Hi,

    I am trying to call Cbc 2.9.6 from my program. When using the tutorial-style approach

    OsiClpSolverInterface osi;
    osi.add .......
    CbcModel model(osi);
    model.branchAndBound();

    there seem to be no cuts and other stuff applied. When using the method from the examples,

    CbcMain0(model);
    const char * argv2[]={"driver4","-solve","-quit"};
    CbcMain1(3,argv2,model);

    there are cuts applied, but obviously different (less aggressive) to the standalone Cbc executable. I also tried CbcSolver class but its method solve() is not found by the linker. So what is the 'standard' way of using the 'default' add-ons?

    Moreover. The attached example crashes both in the standalone Cbc and in the CbcCmain0/1 variant after a few minutes.

    Thanks

    _______________________________________________
    Cbc mailing list
    Cbc@list.coin-or.org
    http://list.coin-or.org/mailman/listinfo/cbc


    

Hi, what is currently good way to have a solution callback in Cbc? the 
interrupt example shows 2 ways, don't know which is right.

Moreover, it says that the solution would be given for the preprocessed 
model. Is it possible to produce one for the original? Is it possible to 
call other functions from inside, such as number of nodes, dual bound?

Thanks

From john.forrest at fastercoin.com  Thu Oct  8 10:34:15 2015
From: john.forrest at fastercoin.com (John Forrest)
Date: Thu, 8 Oct 2015 15:34:15 +0100
Subject: [Cbc] Solution callbacks
In-Reply-To: <5615F778.9020601@monash.edu>
References: <5615F778.9020601@monash.edu>
Message-ID: <56167EE7.6000607@fastercoin.com>

Gleb,

On 08/10/15 05:56, Gleb Belov wrote:
> Hi, what is currently good way to have a solution callback in Cbc? the 
> interrupt example shows 2 ways, don't know which is right.
>

It is the event handling code you would be using.
> Moreover, it says that the solution would be given for the 
> preprocessed model. Is it possible to produce one for the original? 

At present no.  In principle not difficult.  First the callback function 
would have to be modified to get passed the CglPreProcess object - 
easy.  Then in event handler you could make a copy of object and 
postsolve (you need a copy as postsolve deletes data).
> Is it possible to call other functions from inside, such as number of 
> nodes, dual bound?

Yes - you have CbcModel * model_ so things like that are available (or 
could easily be made available)

>
> Thanks
>

John Forrest


 */
