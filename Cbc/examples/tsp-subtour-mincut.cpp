// Copyright (C) 2019, Gleb.Belov@monash.edu
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>
#include <iomanip>
#include <sstream>
#include <vector>
#include <exception>
#define OLD_BOOST_JSON
#ifdef OLD_BOOST_JSON
  #include <boost/property_tree/ptree.hpp>
  #include <boost/property_tree/json_parser.hpp>
  namespace pt = boost::property_tree;
#endif

#include "CoinError.hpp"

#include "MIP_wrap.hh"
#include "mip_seccutgen.h"
#include "utils.hh"

/************************************************************************

  This TSP solver can use MTZ formulation and/or Subtour Elimination cuts

  It can be used in 3 ways.
    - runTSPApp() handles the complete application
      including command-line parameters.
    - solveTSP() is interfaced by structs describing TSP instance.
    - addCircuit() just adds the circuit constraint on given variables
      (MTZ formulation and/or SEC cut handler), so it can be used
      to add the TSP subproblem in a larger model.

  These functions can be placed into a .h as the interface.
*/

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// BEGIN INTERFACE ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
/// This function interfaces the whole application for possible use as a library
int runTSPApp(int argc, const char* const *argv);

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
void readTSP(std::istream& is, TSPIO::Input::InstanceData& dat);
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
      readTSP(ifs, d_inp.dat);
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

void readTSP(std::istream& is, TSPIO::Input::InstanceData& dat) {
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


