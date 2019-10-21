 
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

#ifndef __MIP_OSICBC_WRAPPER_H__
#define __MIP_OSICBC_WRAPPER_H__

#include "MIP_wrap.hh"
                    // CMakeLists.txt needs OSICBC_HOME defined
// #include <coin/CoinPackedVector.hpp>
// #include <coin/CoinPackedMatrix.hpp>
// #include <coin/CoinShallowPackedVector.hpp>
// #include <coin/CoinTime.hpp>
// #include <coin/OsiSolverInterface.hpp>
//  #include <coin/OsiCbcSolverInterface.hpp>
#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"
// #include <coin/CbcSolver.hpp>


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

#endif  // __MIP_OSICBC_WRAPPER_H__
