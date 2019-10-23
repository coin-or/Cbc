// Copyright (C) 2019, Gleb.Belov@monash.edu
// This code is licensed under the terms of the Eclipse Public License (EPL).
// Much of the Cbc interface code has been taken from MiniZinc (www.minizinc.org)

/************************************************************************
  This TSP solver can use MTZ formulation and/or Subtour Elimination cuts
  Run with -h to see options

  The example code can be reused in several ways.
    - runTSPApp() handles the complete application
      including command-line parameters.
    - solveTSP() is interfaced by structs describing TSP instance.
    - addMIPCircuit() just adds the circuit constraint on given variables
      (MTZ formulation and/or SEC cut handler), so it can be used
      to add the TSP subproblem in a larger model.

  Those declarations can be placed into .h files as the interface.
*/

#include <sstream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <string>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <chrono>
#include <exception>
#define OLD_BOOST_JSON
#ifdef OLD_BOOST_JSON
  #include <boost/property_tree/ptree.hpp>
  #include <boost/property_tree/json_parser.hpp>
  namespace pt = boost::property_tree;
#endif
#define COMPILE_BOOST_MINCUT
#ifdef COMPILE_BOOST_MINCUT
  #include <boost/graph/adjacency_list.hpp>
  #include <boost/graph/graph_traits.hpp>
  #include <boost/graph/one_bit_color_map.hpp>
  #include <boost/graph/stoer_wagner_min_cut.hpp>
  #include <boost/property_map/property_map.hpp>
  #include <boost/typeof/typeof.hpp>
#endif

#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CoinError.hpp"
#include "CbcSolver.hpp"
#include "CbcConfig.h"
#include "ClpConfig.h"
#include "CbcEventHandler.hpp"
#include "CglPreProcess.hpp"
#include "CglCutGenerator.hpp"
#include "CoinSignal.hpp"

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// BEGIN INTERFACE ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
/// This function interfaces the whole application for possible use as a library
int runTSPApp(int argc, const char* const *argv);

/// An abstract MIP wrapper.
/// Does not include MZN stuff so can be used independently
/// although it's limited to the MZN solving needs.
class MIP_wrapper {

public:
  class Options  {
  public:
    bool fVerbose=true;
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
    /// Facilitate lhs computation of a cut
    static
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
/// Input/output structs for solveTSP()
struct TSPIO {
  struct Input {
    struct InstanceData {
      std::vector<double> dist;           // TSP distance matrix
      int n;
      bool check() const {
        return n*n == dist.size();
      }
    };
    InstanceData dat;
    struct Param {
      int fMTZ=1;                      // whether to use MTZ formulation
      int fSEC=1;                      // whether to use SEC cuts

      bool fVerbose=true;
    };
    Param prm;
    /// Strictly speaking, this is already method's options
    MIP_wrapper::Options opt_mip;
  };
  struct Output {
    std::vector<int> succ;              // the successor list
    double obj;
  };
};
/// The function for reading instance data
void readTSP(std::istream& is, TSPIO::Input::InstanceData& dat, bool fJSON);
/// The function solving TSP
int solveTSP(const TSPIO::Input& inp, TSPIO::Output& outp);

/// The "global constraint macro" adding MTZ reformulation and/or SEC handler
/// to a given model
struct MIPCircuitDef {
  MIP_wrapper* wrp;
  int n;                           // number of nodes
  const int* flowVars;             // binary flow variable indices x11, x12, ...
  bool fMTZ, fSEC;
};
void addMIPCircuit(const MIPCircuitDef& cdef);

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// END INTERFACE ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////// int main() ////////////////////////////////////////
int main(int argc, const char **argv)
{
  try {
    return runTSPApp(argc, argv);
  } catch (const CoinError& e) {
    e.print();
    if (e.lineNumber() >= 0)
      std::cerr << "This was from a CoinAssert" << std::endl;
    exit(-1);
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    exit(-1);
  }
  return 0;
}


#ifdef MZN_HAS_LLROUND
#include <cmath>
namespace MiniZinc {
  inline
  long long int round_to_longlong(double v) {
    return ::llround(v);
  }
}
#else
namespace MiniZinc {
  inline
  long long int round_to_longlong(double v) {
    return static_cast<long long int>(v < 0 ? v-0.5 : v+0.5);
  }
}
#endif

namespace MiniZinc {

// #define __MZN_PRINTATONCE__
#ifdef __MZN_PRINTATONCE__
  #define __MZN_PRINT_SRCLOC(e1, e2) \
    std::cerr << '\n' << __FILE__ << ": " << __LINE__ << " (" << __func__ \
     << "): not " << e1 << ":  " << std::flush; \
     std::cerr << e2 << std::endl
#else
  #define __MZN_PRINT_SRCLOC(e1, e2)
#endif
#define MZN_ASSERT_HARD( c ) \
   do { if ( !(c) ) { __MZN_PRINT_SRCLOC( #c, "" ); throw std::runtime_error( #c ); } } while (0)
#define MZN_ASSERT_HARD_MSG( c, e ) \
   do { if ( !(c) ) { __MZN_PRINT_SRCLOC( #c, e ); \
     std::ostringstream oss; oss << "not " << #c << ":  " << e; \
     throw std::runtime_error( oss.str() ); } } while (0)

  inline bool beginswith(const std::string& s, const std::string& t) {
    return s.compare(0, t.length(), t)==0;
  }

  /// https://stackoverflow.com/questions/874134/find-out-if-string-ends-with-another-string-in-c
  inline bool endswith(std::string const &fullString, std::string const &ending) {
      if (fullString.length() >= ending.length()) {
          return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
      } else {
          return false;
      }
  }

  inline void checkIOStatus( bool fOk, std::string msg, bool fHard=1 )
  {
    if ( !fOk ) {
#ifdef _MSC_VER
      char errBuf[1024];
      strerror_s(errBuf, sizeof(errBuf), errno);
#else
      char* errBuf = strerror(errno);
#endif
      std::cerr << "\n  " << msg
        << ":   " << errBuf << "." << std::endl;
      MZN_ASSERT_HARD_MSG ( !fHard, msg << ": " << errBuf );
    }
  }

  template <class T> inline bool assignStr(T*, const std::string ) { return false; }
  template<> inline bool assignStr(std::string* pS, const std::string s ) {
    *pS = s;
    return true;
  }

  /// A simple per-cmdline option parser
  class CLOParser {
    int& i;              // current item
    std::vector<std::string>& argv;

  public:
    CLOParser( int& ii, std::vector<std::string>& av )
      : i(ii), argv(av) { }
    template <class Value=int>
    inline bool get(  const char* names, // space-separated option list
                      Value* pResult=nullptr, // pointer to value storage
                      bool fValueOptional=false // if pResult, for non-string values
                ) {
      return getOption( names, pResult, fValueOptional );
    }
    template <class Value=int>
    inline bool getOption(  const char* names, // space-separated option list
                            Value* pResult=nullptr, // pointer to value storage
                            bool fValueOptional=false // if pResult, for non-string values
                ) {
      assert(0 == strchr(names, ','));
      assert(0 == strchr(names, ';'));
      if( i>=argv.size() )
        return false;
      std::string arg( argv[i] );
      /// Separate keywords
      std::string keyword;
      std::istringstream iss( names );
      while ( iss >> keyword ) {
        if ( ((2<keyword.size() || 0==pResult) && arg!=keyword) ||  // exact cmp
          (0!=arg.compare( 0, keyword.size(), keyword )) )           // truncated cmp
          continue;
        /// Process it
        bool combinedArg = false; // whether arg and value are combined in one string (like -Ggecode)
        if ( keyword.size() < arg.size() ) {
          if ( 0==pResult )
            continue;
          combinedArg = true;
          arg.erase( 0, keyword.size() );
        } else {
          if ( 0==pResult )
            return true;
          i++;
          if( i>=argv.size() ) {
            --i;
            return fValueOptional;
          }
          arg = argv[i];
        }
        assert( pResult );
        if ( assignStr( pResult, arg ) )
          return true;
        std::istringstream iss( arg );
        Value tmp;
        if ( !( iss >> tmp ) ) {
          if (!combinedArg)
            --i;
          if ( fValueOptional ) {
            return true;
          }
          // Not print because another agent can handle this option
//           cerr << "\nBad value for " << keyword << ": " << arg << endl;
          return false;
        }
        *pResult = tmp;
        return true;
      }
      return false;
    }
  };  // class CLOParser

  /// This class prints a value if non-0 and adds comma if not 1st time
  class HadOne {
    bool fHadOne=false;
  public:
    template <class N>
    std::string operator()(const N& val, const char* descr=0) {
      std::ostringstream oss;
      if ( val ) {
        if ( fHadOne )
          oss << ", ";
        fHadOne=true;
        oss << val;
        if ( descr )
          oss << descr;
      }
      return oss.str();
    }
    void reset() { fHadOne=false; }
    operator bool() const { return fHadOne; }
    bool operator!() const { return !fHadOne; }
  };

  /// Split a string into words
  /// Add the words into the given vector
  inline void split(const std::string& str, std::vector<std::string>& words) {
    std::istringstream iss(str);
    std::string buf;
    while (iss) {
      iss >> buf;
      words.push_back(buf);
    }
  }

  /// Puts the strings' c_str()s into the 2nd argument.
  /// The latter is only valid as long as the former isn't changed.
  inline void vecString2vecPChar(const std::vector<std::string>& vS, std::vector<const char*>& vPC) {
    vPC.resize(vS.size());
    for ( size_t i=0; i<vS.size(); ++i ) {
      vPC[i] = vS[i].c_str();
    }
  }

}

/// SEC cut generator for circuit
class MIP_SECCutGen : public MIP_wrapper::CutGen {
  MIP_SECCutGen() { }
  MIP_wrapper* pMIP=0;
public:
  MIP_SECCutGen( MIP_wrapper* pw ) : pMIP(pw) { }
  /// Say what type of cuts
  virtual int getMask() { return MIP_wrapper::MaskConsType_Lazy | MIP_wrapper::MaskConsType_Usercut; }
  /// The flow variables
  std::vector<MIP_wrapper::VarId> varXij;
  int nN=0;                        // N nodes
  /// returns error message if fails
  std::string validate() const;
  void generate(const MIP_wrapper::Output&, MIP_wrapper::CutInput&);
  void print( std::ostream& );
};


////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// runTSPApp() ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
int runTSPApp(int argc, const char* const *argv) {
  class TSPApp {
    const int d_argc;
    std::vector<std::string> d_argv;
    int d_result=0;                         // from the solver
  public:
    int run() {
      if (readParams()) {
        readInput();
        solve();
        writeResult();
      }
      return d_result;
    }
  private:
    TSPIO::Input d_inp;
    TSPIO::Output d_outp;
    std::string d_inpFile;

  protected:
    bool readParams() {
      if (d_argc < 2)
        throw std::runtime_error("Need at least an input filename");
      for (int i=1; i<d_argc; ++i) {
        MiniZinc::CLOParser cop( i, d_argv );
        if (cop.get( "-h --help" )) {
          printHelp();
          return false;
        }
        if (cop.get( "--version -v" )) {
          d_inp.prm.fVerbose = d_inp.opt_mip.fVerbose = true;         // version printed from the CBC wrapper
        } else if (cop.get("--MTZ --mtz -MTZ --fMTZ", &d_inp.prm.fMTZ)) {
        } else if (cop.get("--SEC --sec -SEC --fSEC", &d_inp.prm.fSEC)) {
        } else
          if (d_inp.opt_mip.processOption(i, d_argv)) {
          } else
            d_inpFile = d_argv[i];
      }
      return true;
    }
    void printHelp() {
      std::cout
          << "CBC TSP solver." << std::endl
          << "Usage: exe <options> input_file." << std::endl
        << "\nGeneral options:" << std::endl
        << "  --help, -h\n    Print this help message." << std::endl
        << "  --version\n    Print version information." << std::endl
        << "  -v\n    Print progress/log statements." << std::endl
        << "  --MTZ <n>\n    MTZ reformulation on(1)*/off(0)." << std::endl
        << "  --SEC <n>\n    Subtour Elimination Cuts on(1)*/off(0)." << std::endl
           ;

      MIP_wrapper::Options::printHelp(std::cout);
    }
    void readInput() {
      if (getParams().fVerbose)
        std::cout << "Reading instance from '" << d_inpFile << "'..." << std::endl;
      std::ifstream ifs(d_inpFile);
      readTSP(ifs, d_inp.dat,
              MiniZinc::endswith(d_inpFile, ".json") || MiniZinc::endswith(d_inpFile, ".JSON"));
      if (!ifs.good())
        throw std::runtime_error("I/O error.");
      if (getParams().fVerbose)
        std::cout << "Got instance with " << getData().n
                  << " nodes." << std::endl;
      if (!ifs.good())
        throw std::runtime_error("I/O error.");
    }
    void solve() {
      d_result = solveTSP(d_inp, d_outp);
    }
    void writeResult() {

    }

    const TSPIO::Input::Param& getParams() const { return d_inp.prm; }
    const TSPIO::Input::InstanceData& getData() const { return d_inp.dat; }
  public:
    TSPApp(int argc, const char* const *argv)
      : d_argc(argc), d_argv(argv, argv+argc) { }
  };
  TSPApp tspApp(argc, argv);
  return tspApp.run();
}

void readTSP_dist(std::istream& is, TSPIO::Input::InstanceData& dat) {
  is >> dat.n;
  dat.dist.resize(dat.n * dat.n);
  for (int ij=0; (ij<dat.n*dat.n); ++ij) {
    is >> dat.dist[ij];
  }
}
void readTSP_json(std::istream& is, TSPIO::Input::InstanceData& dat) {
  pt::ptree pTree;
  pt::read_json(is, pTree);
  /// Could understand several input styles
  boost::optional<const pt::ptree&> optEucl =
      ((const pt::ptree&)(pTree)).get_child_optional("TSP__EUCLIDIAN");
  if (optEucl) {
    /// Point coordinates
    const pt::ptree& pX = optEucl.get().get_child("X");
    const pt::ptree& pY = optEucl.get().get_child("Y");
    std::vector<double> X, Y;
    for (pt::ptree::const_iterator it=pX.begin(); it!=pX.end(); ++it)
      X.push_back(it->second.get_value<double>());
    for (pt::ptree::const_iterator it=pY.begin(); it!=pY.end(); ++it)
      Y.push_back(it->second.get_value<double>());
    dat.n = X.size();
    if (dat.n != Y.size())
      throw std::runtime_error("TSP__EUCLIDIAN: sizes of X and Y differ");
    dat.dist.resize(dat.n * dat.n);
    /// Converting to dist + cost
    for (size_t i=0; i<dat.n; ++i)
      for (size_t j=0; j<dat.n; ++j)
        dat.dist[i*dat.n + j] =
            std::sqrt(std::pow(X[i]-X[j], 2.0) + std::pow(Y[j]-Y[i], 2.0));
  } else
    throw std::runtime_error("None of input styles recognized.");
}
void readTSP(std::istream& is, TSPIO::Input::InstanceData& dat, bool fJSON) {
  if (fJSON)
    readTSP_json(is, dat);
  else
    readTSP_dist(is, dat);
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// solveTSP() ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

int solveTSP(const TSPIO::Input& inp, TSPIO::Output& outp) {
  class CBC_TSP_Solver {
    const TSPIO::Input& d_inp;
    TSPIO::Output& d_outp;
    int d_result=0;
  public:
    ////////////////////////// The high-level TSP solver logic //////////////////////////////////
    int run() {
      initModel();
      addFlowVars();                      // includes objective coefficients
      addCircuitGlobalConstraint();
      solve();
      getResult();
      return d_result;
    }
  private:
    /// MIP interface
    MIP_wrapper* d_mip=NULL;

    std::vector<int> xFlow;               // indices of the flow variables
  protected:
    void initModel() {
      getMIP()->setObjSense(-1);
    }
    /// Flow vars define the objective
    void addFlowVars() {
      for (int i = 0; i < getData().n; ++i) {
        for (int j = 0; j < getData().n; ++j) {
          std::ostringstream oss;
          oss << "TSP" << getData().n << "x_" << i << '_' << j;
          xFlow.push_back( getMIP()->addVar(getData().dist[i*getData().n + j],
                           0, i!=j ? 1 : 0,            // upper bound 0 for i==j
              MIP_wrapper::BINARY, oss.str()) );
        }
      }
    }
    void addCircuitGlobalConstraint() {
      addMIPCircuit( { getMIP(), getData().n, xFlow.data(),
                       getParams().fMTZ, getParams().fSEC } );
    }
    void solve() {
      getMIP()->solve();
    }
    void getResult() {

    }

    //////////////////////////// Internal stuff //////////////////////////////////
    const TSPIO::Input::Param& getParams() const { return d_inp.prm; }
    const TSPIO::Input::InstanceData& getData() const { return d_inp.dat; }
    TSPIO::Output& getOutp() { return d_outp; }
//    OsiSolverInterface& getOsi() { return d_osi; }
//    CbcModel& getModel() { return d_model; }
    MIP_wrapper* getMIP() { return d_mip; }
  public:
    CBC_TSP_Solver(const TSPIO::Input& inp, TSPIO::Output& outp)
      : d_inp(inp), d_outp(outp), d_mip(createCBCWrapper(&inp.opt_mip)) { }
    ~CBC_TSP_Solver() { delete d_mip; }
  };
  CBC_TSP_Solver cbcTSP(inp, outp);
  return cbcTSP.run();
}


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// addCircuit() ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

/// The encoding of MTZ follows fzn_circuit.mzn in the MiniZinc linear library
void addMIPCircuit(const MIPCircuitDef& cdef) {
  class MIPCircuit {
    const MIPCircuitDef& cdef;
  public:
    void run() {
      if (!cdef.fMTZ && !cdef.fSEC)
        throw std::logic_error("addMIPCircuit(): neither MTZ nor SEC requested");
      addMTZ();
      if (cdef.fSEC)
        addSECHandler();
    }
  private:
    std::vector<int> xOrder;
  protected:
    void addMTZ() {
      setVarBounds();
      addFlowConstraints();
      if (cdef.fMTZ) {
        addMTZOrderVars();
        addMTZOrderConstraints();
      }
      addRedundantAntiLoops();
    }
    void setVarBounds() {
      for (int i = 0; i < getN(); ++i) {
        getMIP()->setVarUB(getX(i, i), 0);
      }
    }
    /// add alldiff()
    void addFlowConstraints() {
      const std::vector<double> all1(getN(), 1.0);               // n-vec of 1's
      for (int i = 0; i < getN(); ++i) {
        std::vector<int> colInd(getN());                   // indices of column i
        for (int j = 0; j < getN(); ++j) {
          colInd[j] = getX(j, i);
        }
        /// Outflow
        getMIP()->addRow(getN(), cdef.flowVars+i*getN(), all1.data(),
                         MIP_wrapper::EQ, 1.0);
        /// Inflow
        getMIP()->addRow(getN(), colInd.data(), all1.data(),
                         MIP_wrapper::EQ, 1.0);
      }
    }
    void addMTZOrderVars() {
      xOrder.resize(getN()-1);
      for (int i = 1; i < getN(); ++i) {
        std::ostringstream oss;
        oss << "TSP" << getN() << "_MTZOrder_" << i;
        xOrder[i-1] = getMIP()->addVar(0.0, 2.0, getN(),
                                     MIP_wrapper::INT,          // integer seems better!
                                     oss.str());
      }
    }
    void addMTZOrderConstraints() {
      const std::vector<double> coefs {{1.0, -1.0, getN()-1, getN()-3}};
      for (int i = 1; i < getN(); ++i) {
        for (int j = 1; j < getN(); ++j) {
          if (i!=j) {                          // TODO only where ub(xij)>0
            std::vector<int> vars {{xOrder.at(i-1), xOrder.at(j-1),
                                   getX(i, j),
                    getX(j, i)                // the Desrochers & Laporte '91 term
                              // seems only efficient with all 2-node SECs added apriori
                                   }};
            getMIP()->addRow(coefs.size(), vars.data(), coefs.data(),
                             MIP_wrapper::LQ, getN()-2);
          }
        }
      }
    }
    /// These are size-2 SECs
    void addRedundantAntiLoops() {
      const std::vector<double> coefs {{1.0, 1.0}};
      for (int i = 0; i < getN(); ++i) {
        for (int j = i+1; j < getN(); ++j) {
          std::vector<int> vars {{getX(i, j), getX(j, i)}};
          getMIP()->addRow(coefs.size(), vars.data(), coefs.data(),
                           MIP_wrapper::LQ, 1.0);
        }
      }
    }
    void addSECHandler() {
      MIP_SECCutGen* pCG = new MIP_SECCutGen(getMIP());
      pCG->nN = getN();
      pCG->varXij.assign(cdef.flowVars, cdef.flowVars+getN()*getN());
      const std::string sVld = pCG->validate();
      if (sVld.size()) {
        throw std::runtime_error( sVld );
      }
      getMIP()->registerCutGenerator(pCG);
    }

    /// Helper funcs
    int getN() const { return cdef.n; }
    /// Index of the X[i, j] flow variable
    int getX(int i, int j) const { return cdef.flowVars[i*cdef.n + j]; }
    MIP_wrapper* getMIP() { return cdef.wrp; }
  public:
    MIPCircuit(const MIPCircuitDef& cd) : cdef(cd) { }
  };
  MIPCircuit mipc(cdef);
  mipc.run();
}



///////////////////////////////////////////////////////////////////////////
/// Interfaces && technical

class MIP_osicbc_wrapper : public MIP_wrapper {
//     OsiCbcSolverInterface osi;   // deprecated in Cbc 2.9.6
    OsiClpSolverInterface osi;
//     CoinPackedMatrix* matrix = 0;
    int             error;
    std::string          osicbc_buffer;   // [CBC_MESSAGEBUFSIZE];
//     string          osicbc_status_buffer; // [CBC_MESSAGEBUFSIZE];

    std::vector<double> x;

    // To add constraints:
//     vector<int> rowStarts, columns;
    std::vector<CoinPackedVector> rows;
    std::vector<double> //element,
      rowlb, rowub;

    std::unordered_map<VarId, double> warmstart;           // this accumulates warmstart infos

    std::string cbc_cmdOptions;

  public:

    MIP_osicbc_wrapper(const Options* opt) : MIP_wrapper(opt) { openOSICBC(); }
    virtual ~MIP_osicbc_wrapper() { closeOSICBC(); }

    static std::string getDescription();
    static std::string getVersion();
    static std::string getId(void);
    static std::string getName(void);
    static std::vector<std::string> getTags(void);
    static std::vector<std::string> getStdFlags(void);
    static std::string needDllFlag(void) { return ""; }

    void printVersion(std::ostream& );
    void printHelp(std::ostream& );
//       Statistics& getStatistics() { return _statistics; }

//      IloConstraintArray *userCuts, *lazyConstraints;

    /// derived should overload and call the ancestor
//     virtual void cleanup();
    void openOSICBC() { }
    void closeOSICBC() { }

    /// actual adding new variables to the solver
    virtual void doAddVars(size_t n, double *obj, double *lb, double *ub,
      VarType *vt, std::string *names);

    void addPhase1Vars() {
      if (fVerbose)
        std::cerr << "  MIP_osicbc_wrapper: delaying physical addition of variables..." << std::endl;
    }

    /// adding a linear constraint
    virtual void addRow(int nnz, const int *rmatind, const double* rmatval,
                        LinConType sense, double rhs,
                        int mask = MaskConsType_Normal,
                        const std::string& rowName = "");
    /// adding an implication
//     virtual void addImpl() = 0;

    virtual bool addWarmStart( const std::vector<VarId>& vars, const std::vector<double> vals );

    virtual void setObjSense(int s);   // +/-1 for max/min

    virtual double getInfBound() { return osi.getInfinity(); }

    virtual int getNCols() {
      int nc = osi.getNumCols();
      return nc ? nc : colLB.size();
    }
    virtual int getNColsModel() {
      return osi.getNumCols();
    }
    virtual int getNRows() {
      if (rowlb.size())
        return rowlb.size();
      return osi.getNumRows();
    }

//     void setObjUB(double ub) { objUB = ub; }
//     void addQPUniform(double c) { qpu = c; } // also sets problem type to MIQP unless c=0

    virtual void solve();

    /// OUTPUT:
    virtual const double* getValues() { return output.x; }
    virtual double getObjValue() { return output.objVal; }
    virtual double getBestBound() { return output.bestBound; }
    virtual double getCPUTime() { return output.dCPUTime; }

    virtual Status getStatus()  { return output.status; }
    virtual std::string getStatusName() { return output.statusName; }

     virtual int getNNodes() { return output.nNodes; }
     virtual int getNOpen() { return output.nOpenNodes; }

//     virtual int getNNodes() = 0;
//     virtual double getTime() = 0;

  protected:
//     OsiSolverInterface& getOsiSolver(void) { return osi; }

    void wrap_assert(bool , std::string , bool fTerm=true);

    /// Need to consider the 100 status codes in OSICBC and change with every version? TODO
    Status convertStatus(CbcModel *pModel);
    Status convertStatus();
};


using namespace std;


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
           rc.setGloballyValid(true);
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


namespace Algorithms {

  /// An interface to some min-cut algorithm for undirected graphs
  class MinCut {
  public:
    /// INPUT
    int nNodes=0;
    std::vector< std::pair< int, int > > edges;
    std::vector< double > weights;
    /// OUTPUT
    std::vector< bool > parities;
    double wMinCut=1e100;
    /// Invocation
    void solve();
  };

}  // namespace Algorithms

std::string MIP_SECCutGen::validate() const {
  std::ostringstream oss;
  /// Check that diagonal flows are 0
  for (int i=0; i<nN; ++i)
    if (pMIP->colUB[varXij[i*nN+i]] > 0.0)
      oss << "SECutGen with " << nN
          << " cities: diagonal flow " << (i+1)
          << " has UB=" << pMIP->colUB[varXij[i*nN+i]]
          << "\n";
  return oss.str();
}

void MIP_SECCutGen::generate(const MIP_wrapper::Output& slvOut, MIP_wrapper::CutInput& cutsIn) {
  assert( pMIP );
  /// Extract graph, converting to undirected
  typedef std::map< std::pair< int, int >, double > TMapFlow;
  TMapFlow mapFlow;
  for ( int i=0; i<nN; ++i ) {
    for ( int j=0; j<nN; ++j ) {
      const double xij = slvOut.x[ varXij[ nN*i + j ] ];
      if ( i==j )
        MZN_ASSERT_HARD_MSG( 1e-4 > std::fabs(xij), "circuit: X[" << (i+1) << ", " << (j+1) << "]==" << xij );
      MZN_ASSERT_HARD_MSG( -1e-4 < xij && 1.0+1e-4 > xij,           // adjusted from 1e-6 to 1e-4 for CBC. 7.8.19
                           "circuit: X[" << (i+1) << ", " << (j+1) << "]==" << xij );
      if ( 1e-4 <= xij ) {
        mapFlow[ std::make_pair( std::min(i,j), std::max(i,j) ) ] += xij;
      }
    }
  }
  /// Invoking Min Cut
//   cerr << "  MIN CUT... " << flush;
  Algorithms::MinCut mc;
  mc.nNodes = nN;
  mc.edges.reserve( mapFlow.size() );
  mc.weights.reserve( mapFlow.size() );
  for ( const auto& mf: mapFlow ) {
    mc.edges.push_back( mf.first );
    mc.weights.push_back( mf.second );
  }
  mc.solve();
  /// Check if violation
  if ( mc.wMinCut <= 1.999 ) {
    MIP_wrapper::CutDef cut( MIP_wrapper::GQ, MIP_wrapper::MaskConsType_Lazy | MIP_wrapper::MaskConsType_Usercut );
    cut.rhs = 1.0;
    int nCutSize=0;
    constexpr int nElemPrint = 20;
    // cerr << "  CUT: [ ";
    for ( int i=0; i<nN; ++i )
    if ( mc.parities[i] ) {
      ++nCutSize;
      //if ( nCutSize<=nElemPrint )
      //  cerr << (i+1) << ", ";
      //else if ( nCutSize==nElemPrint+1 )
      //  cerr << "...";
      for ( int j=0; j<nN; ++j )
      if ( !mc.parities[j] ) {
        cut.addVar( varXij[ nN*i + j ], 1.0 );
      }
    }
    // cerr << "]. " << flush;
    double dViol = cut.computeViol( slvOut.x, slvOut.nCols );
    if ( dViol > 0.0001 ) {   // ?? PARAM?  TODO. See also min cut value required
      cutsIn.push_back( cut );
      std::cerr << "  SEC cut: violation=" << dViol
        << "  N NODES: " << nN
        << "  |X|: : " << nCutSize
        << std::endl;
    } else {
      MZN_ASSERT_HARD_MSG( 0, "  SEC cut: N nodes = " << nN << ": violation = " << dViol
        << ": too small compared to the min-cut value " << (2.0-mc.wMinCut)
      );
    }
  }
}

void MIP_SECCutGen::print(std::ostream&) {
}


///////////////////////////// MIN CUT ///////////////////////////////
void Algorithms::MinCut::solve() {
#ifndef COMPILE_BOOST_MINCUT
  throw std::runtime_error( "MIP/circuit: Subtour Elimination Constraints: MinCut::solve not compiled. "
  "Compile with COMPILE_BOOST_MINCUT (needs boost), or use -D nSECcuts=0 for flattening." );
#else
  typedef std::pair< unsigned long, unsigned long > edge_t;

  // A graphic of the min-cut is available at <http://www.boost.org/doc/libs/release/libs/graph/doc/stoer_wagner_imgs/stoer_wagner.cpp.gif>
  using namespace std;

  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
    boost::no_property, boost::property<boost::edge_weight_t, double> > undirected_graph;
  typedef boost::property_map<undirected_graph, boost::edge_weight_t>::type weight_map_type;
  typedef boost::property_traits<weight_map_type>::value_type weight_type;

  // define the 16 edges of the graph. {3, 4} means an undirected edge between vertices 3 and 4.
//   edge_t edges[] = {{3, 4}, {3, 6}, {3, 5}, {0, 4}, {0, 1}, {0, 6}, {0, 7},
//     {0, 5}, {0, 2}, {4, 1}, {1, 6}, {1, 5}, {6, 7}, {7, 5}, {5, 2}, {3, 4}};
  std::vector< edge_t > edges_array;  // = edges;
//   std::vector< weight_type > ws_array;

  edges_array.reserve( edges.size() );
  for ( const auto& edg: edges )
    edges_array.push_back( std::pair<unsigned long, unsigned long>( edg.first, edg.second ) );

  // for each of the 16 edges, define the associated edge weight. ws[i] is the weight for the edge
  // that is described by edges[i].
//   weight_type ws[] = {0.2, 3.1, 1.3, 3.7, 1.5, 2.6, 6.44, 1.26, 8.77, 1.29, 1.95, 80.74, 2.23, 1.94, 1.23, 4.957};

  // construct the graph object. 8 is the number of vertices, which are numbered from 0
  // through 7, and 16 is the number of edges.
  undirected_graph g(edges_array.data(), edges_array.data() + edges_array.size(),
                     weights.data(), nNodes, edges_array.size());

  // define a property map, `parities`, that will store a boolean value for each vertex.
  // Vertices that have the same parity after `stoer_wagner_min_cut` runs are on the same side of the min-cut.
  BOOST_AUTO(parities00, boost::make_one_bit_color_map(num_vertices(g), get(boost::vertex_index, g)));

  // run the Stoer-Wagner algorithm to obtain the min-cut weight. `parities` is also filled in.
  wMinCut = boost::stoer_wagner_min_cut(g, get(boost::edge_weight, g), boost::parity_map(parities00));

  assert( nNodes==num_vertices(g) );

//   parities = parities00;
  parities.clear();
  parities.reserve( nNodes );
  for (size_t i = 0; i < num_vertices(g); ++i) {
    parities.push_back (get(parities00, i));
  }
#endif
}
