/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */

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

#ifndef __MIP_WRAPPER__
#define __MIP_WRAPPER__

#include <vector>
//#include <map>
#include <unordered_map>
#include <string>
#include <iostream>
#include <sstream>
#include <cassert>
#include <chrono>


/// Facilitate lhs computation of a cut
inline
double computeSparse( int n, const int* ind, const double* coef, const double* dense, int nVarsDense ) {
  assert( ind && coef && dense );
  double val=0.0;
  for ( int i=0; i<n; ++i ) {
    assert( ind[i]>=0 );
    assert( ind[i]<nVarsDense );
    val += coef[i] * dense[ ind [ i ] ];
  }
  return val;
}

/// An abstract MIP wrapper.
/// Does not include MZN stuff so can be used independently
/// although it's limited to the MZN solving needs.
class MIP_wrapper {

public:
  class Options  {
  public:
    bool fVerbose=false;
    int nThreads=1;
    std::string sExportModel;
    int nTimeout=0;
    long int nSolLimit = -1;
    double nWorkMemLimit=-1;
    std::string sReadParams;
    std::string sWriteParams;
    bool flag_all_solutions = false;

    double absGap=-1;
    double relGap=1e-8;
    double intTol=1e-8;
//    double objDiff=1.0;

    std::string cmdOptions;              // optiomns directly for the solver

    bool processOption(int& i, std::vector<std::string>& argv);
    static void printHelp(std::ostream& );
  };
  const Options* options=NULL;

public:
    typedef int VarId;    // CPLEX uses int
    enum VarType {
      REAL,
      INT,
      BINARY
    };
    enum LinConType {
      LQ = -1,
      EQ = 0,
      GQ = 1
    };
    
    // CPLEX 12.6.2 advises anti-symmetry constraints to be user+lazy
    static const int MaskConsType_Normal = 1;
    /// User cut. Only cuts off fractional points, no integer feasible points
    static const int MaskConsType_Usercut = 2;
    /// Lazy cut. Can cut off otherwise feasible integer solutions.
    /// Callback should be able to produce previously generated cuts again if needed [Gurobi]
    static const int MaskConsType_Lazy = 4;
    enum Status { OPT, SAT, UNSAT, UNBND, UNSATorUNBND, UNKNOWN, __ERROR };
  public:
    /// Columns for SCIP upfront and with obj coefs:
    std::vector<double> colObj, colLB, colUB;
    std::vector<VarType> colTypes;
    std::vector<std::string> colNames;
//     , rowLB, rowUB, elements;
//     veci whichInt
//     , starts, column;
//     double objUB;
//     double qpu;
    
  public:
    /// Parameter
    bool fVerbose = true;
    
    int nProbType = -2;  // +-1: max/min; 0: sat

  public:
    struct Output {
      Status status;
      std::string statusName="Untouched";
      double objVal = 1e308;
      double bestBound = 1e308;
      int nCols = 0;
      int nObjVarIndex=-1;
      const double *x = 0;
      int nNodes=0;
      int nOpenNodes=0;
      double dWallTime = 0.0;
      std::chrono::time_point<std::chrono::steady_clock> dWallTime0;
      double dCPUTime = 0;
      std::clock_t cCPUTime0 = 0;
    };      
    Output output;

    /// General cut definition, could be used for addRow() too
    class CutDef {
      CutDef() { }
    public:
      CutDef( LinConType s, int m ) : sense( s ), mask( m ) { }
      std::vector<int> rmatind;
      std::vector<double> rmatval;
      LinConType sense=LQ;
      double rhs=0.0;
      int mask = 0; // need to know what type of cuts are registered before solve()  TODO
      std::string rowName = "";
      void addVar( int i, double c ) {
        rmatind.push_back( i );
        rmatval.push_back( c );
      }
      double computeViol( const double* x, int nCols ) {
        double lhs = computeSparse( static_cast<int>(rmatind.size()), rmatind.data(), rmatval.data(), x, nCols );
        if ( LQ==sense ) {
          return lhs-rhs;
        } else if ( GQ==sense ) {
          return rhs-lhs;
        } else
          assert( 0 );
        return 0.0;
      }
    };
    /// Result of cut generation: cut callback fills one
    typedef std::vector<CutDef> CutInput;
    
  public:

    /// Generic cut generator
    /// Callback should be able to produce previously generated cuts again if needed [Gurobi]
    class CutGen {
    public:
      virtual ~CutGen() { }
      /// Say what type of cuts
      virtual int getMask() = 0;
      /// Adds new cuts to the 2nd parameter
      virtual void generate(const MIP_wrapper::Output&, MIP_wrapper::CutInput&) = 0;
      virtual void print( std::ostream& ) { }
    };

  private:

    /// The cut generators
    std::vector< CutGen* > cutGenerators;
  public:

    /// This should be used to provide a cut generator
    void registerCutGenerator( CutGen* pCG ) {
      cbui.cutMask |= pCG->getMask();
      cutGenerators.push_back( pCG );
      provideCutCallback(HandleCutCallback, this);
    }

  protected:
    /// solution callback handler, the wrapper might not have these callbacks implemented
    typedef void (*SolCallbackFn)(const Output& , void* );
    /// cut callback handler, the wrapper might not have these callbacks implemented
    typedef void (*CutCallbackFn)(const Output& , CutInput& ,
                                  void* ,          // the extra info provided
                  bool fMIPSol  // if with a MIP feas sol candidate - lazy cuts only
                                 );
    static void HandleCutCallback(const MIP_wrapper::Output& out, MIP_wrapper::CutInput& in,
                           void* pp, bool fMIPSol) {
      // multi-threading? Assumed to be cared by the solver. TODO
      MIP_wrapper* pw = static_cast<MIP_wrapper*>( pp );
      assert(pw);
      assert(&out);
      assert(&in);
      pw->genCuts( out, in, fMIPSol );
    }
    /// This function calls all cut generators appropriate
    void genCuts(const Output& slvOut, CutInput& cutsIn, bool fMIPSol) {
      for ( auto& pCG : cutGenerators ) {
        if ( !fMIPSol || pCG->getMask() & MIP_wrapper::MaskConsType_Lazy )
          pCG->generate( slvOut, cutsIn );
      }
      /// Select some most violated? TODO
    }
  public:
    struct CBUserInfo {
      MIP_wrapper* wrapper = 0;
      MIP_wrapper::Output* pOutput=0;
      MIP_wrapper::Output* pCutOutput=0;
      void *psi=0;  // external info. Intended to keep MIP_solverinstance
      SolCallbackFn solcbfn=0;
      CutCallbackFn cutcbfn=0;
      /// Union of all flags used for the registered callback cuts
      /// See MaskConstrType_..
      /// Solvers need to know this
      /// In MIP_solverinstance, class CutGen defines getMask() which should return that
      int cutMask = 0; // can be any combination of User/Lazy
      bool fVerb = false;              // used in Gurobi
      bool printed = false;            // whether any solution was output
      double nTimeoutFeas = -1.0;      // >=0 => stop that long after 1st feas
      double nTime1Feas = -1e100;      // time of the 1st feas
    };
    CBUserInfo cbui;

  public:
    MIP_wrapper(const Options* opt)
      : options(opt), fVerbose(opt->fVerbose) {
      cbui.wrapper = this;
      cbui.pOutput = &output;
    }
    virtual ~MIP_wrapper() {
      for (int i=cutGenerators.size(); i--; ) {
        delete cutGenerators[i];
      }
    }

    /// re-create solver object. Called from the base class constructor
//     virtual void resetModel() { };
    
//     virtual void printVersion(ostream& os) { os << "Abstract MIP wrapper"; }
//     virtual void printHelp(ostream& ) { }

  public:
    bool fPhase1Over = false;
  private:
    /// adding a variable just internally (in Phase 1 only that). Not to be used directly.
    virtual VarId addVarLocal(double obj, double lb, double ub, 
                             VarType vt, const std::string& name="") {
//       cerr << "  addVarLocal: colObj.size() == " << colObj.size()
//         << " obj == " <<obj
//         << " lb == " << lb
//         << " ub == " << ub
//         << " vt == " << vt
//         << " nm == " << name
//         << endl;
      colObj.push_back(obj);
      colLB.push_back(lb);
      colUB.push_back(ub);
      colTypes.push_back(vt);
      colNames.push_back(name);
      return static_cast<VarId>(colObj.size()-1);
    }
    /// add the given var to the solver. Asserts all previous are added. Phase >=2. No direct use
    virtual void addVar(int j) {
      assert(j == getNCols());
      assert(fPhase1Over);
      doAddVars(1, &colObj[j], &colLB[j], &colUB[j], &colTypes[j], &colNames[j]);
    }
    /// actual adding new variables to the solver. "Updates" the model (e.g., Gurobi). No direct use
    virtual void doAddVars(size_t n, double *obj, double *lb, double *ub,
      VarType *vt, std::string *names) = 0;

  public:
    /// debugging stuff
//     set<double> sLitValues;
    std::unordered_map<double, VarId> sLitValues;
    
    void setProbType( int t ) { nProbType=t; }
    
    /// adding a variable, at once to the solver
    virtual VarId addVar(double obj, double lb, double ub, 
                             VarType vt, const std::string& name="") {
//       cerr << "  AddVar: " << lb << ":   ";
      VarId res = addVarLocal(obj, lb, ub, vt, name);
      if (fPhase1Over)
        addVar(res);
      return res;
    }
    int nLitVars=0;
    /// adding a literal as a variable. Should not happen in feasible models
    virtual VarId addLitVar(double v) {
      // Cannot do this: at least CBC does not support duplicated indexes    TODO??
//       ++nLitVars;
//       auto itFound = sLitValues.find(v);
//       if (sLitValues.end() != itFound)
//         return itFound->second;
      std::ostringstream oss;
      oss << "lit_" << v << "__" << (nLitVars++);
      std::string name = oss.str();
      size_t pos = name.find('.');
      if (std::string::npos != pos)
        name.replace(pos, 1, "p");
      VarId res = addVarLocal(0.0, v, v, REAL, name);
      if (fPhase1Over)
        addVar(res);
//       cerr << "  AddLitVar " << v << "   (PROBABLY WRONG)" << endl;
      sLitValues[v] = res;
      return res;
    }
    /// adding all local variables upfront. Makes sure it's called only once
    virtual void addPhase1Vars() {
      assert(0 == getNColsModel());
      assert(! fPhase1Over);
      if (fVerbose)
        std::cerr << "  MIP_wrapper: adding the " << colObj.size() << " Phase-1 variables..." << std::flush;
      doAddVars(colObj.size(), &colObj[0], &colLB[0], &colUB[0], &colTypes[0], &colNames[0]);
      if (fVerbose)
        std::cerr << " done." << std::endl;
      fPhase1Over = true;    // SCIP needs after adding
    }

    /// var bounds
    void setVarBounds( int iVar, double lb, double ub )  {
      if (fPhase1Over)
        doSetVarBounds(iVar, lb, ub);
      colLB.at(iVar) = lb;
      colUB.at(iVar) = ub;
    }
    void setVarLB( int iVar, double lb )  { setVarBounds(iVar, lb, colUB.at(iVar)); }
    void setVarUB( int iVar, double ub )  { setVarBounds(iVar, colLB.at(iVar), ub); }
  protected:
    /// This goes into the actual solver after phase 1
    virtual void doSetVarBounds( int iVar, double lb, double ub )  {
      throw std::runtime_error("This MIP wrapper has not implemented doSetVarBounds()");
    }
  public:
    /// adding a linear constraint
    virtual void addRow(int nnz, const int *rmatind, const double* rmatval,
                        LinConType sense, double rhs,
                        int mask = MaskConsType_Normal,
                        const std::string& rowName = "") = 0;
    /// Indicator constraint: x[iBVar]==bVal -> lin constr
    virtual void addIndicatorConstraint(int iBVar, int bVal, int nnz, int *rmatind, double* rmatval,
                        LinConType sense, double rhs,
                        std::string rowName = "") { throw std::runtime_error("Indicator constraints not supported. "); }
                
    /// Bounds disj for SCIP
    virtual void addBoundsDisj(int n, double *fUB, double *bnd, int* vars,
                               int nF, double *fUBF, double *bndF, int* varsF,
                        std::string rowName = "") { throw std::runtime_error("Bounds disjunctions not supported. "); }

    /// Cumulative, currently SCIP only
    virtual void addCumulative(int nnz, int *rmatind, double* d, double* r, double b, std::string rowName="")
    { throw std::runtime_error("Cumulative constraints not supported. "); }
    /// 0: model-defined level, 1: free, 2: uniform search
    virtual int getFreeSearch() { return 1; }
    /// Return 0 if ignoring searches
    virtual bool addSearch( const std::vector<VarId>& vars, const std::vector<int> pri ) {
      return false;
    }
    /// Return 0 if ignoring warm starts    
    virtual bool addWarmStart( const std::vector<VarId>& vars, const std::vector<double> vals ) {
      return false;
    }
                        
    int nAddedRows = 0;   // for name counting
    int nIndicatorConstr = 0;
    /// adding an implication
//     virtual void addImpl() = 0;
    virtual void setObjSense(int s) = 0;   // +/-1 for max/min
    
    virtual double getInfBound() = 0;
                        
    virtual int getNCols() = 0;
    virtual int getNColsModel() { return getNCols(); }  // from the solver
    virtual int getNRows() = 0;
                        
//     void setObjUB(double ub) { objUB = ub; }
//     void addQPUniform(double c) { qpu = c; } // also sets problem type to MIQP unless c=0

  protected:
    /// Set solution callback. Thread-safety??
    /// solution callback handler, the wrapper might not have these callbacks implemented
    virtual void provideSolutionCallback(SolCallbackFn cbfn, void* info) {
      assert(cbfn);
      cbui.pOutput = &output;
      cbui.psi = info;
      cbui.solcbfn = cbfn;
    }
    /// solution callback handler, the wrapper might not have these callbacks implemented
    virtual void provideCutCallback(CutCallbackFn cbfn, void* info) {
      assert(cbfn);
      cbui.pCutOutput = 0;  // &outpCuts;   thread-safety is to be provided by the solver
      cbui.psi = info;
      cbui.cutcbfn = cbfn;
    }

  public:
    virtual void solve() = 0; 
    
    /// OUTPUT, should also work in a callback
    virtual const double* getValues() = 0;
    virtual double getObjValue() = 0;
    virtual double getBestBound() = 0;
    virtual double getWallTimeElapsed() { return output.dWallTime; }
    virtual double getCPUTime() = 0;
    
    virtual Status getStatus() = 0;
    virtual std::string getStatusName() = 0;

     virtual int getNNodes() = 0;
     virtual int getNOpen() = 0;

    /// Default MZN library for MIP
    static std::string getMznLib();

  };

/// Wrapper factory
/// Pointer not owned
MIP_wrapper* createCBCWrapper(const MIP_wrapper::Options* opt);

#endif  // __MIP_WRAPPER__
