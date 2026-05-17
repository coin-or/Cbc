//-----------------------------------------------------------------------------
// name:     Cgl Lifed Simple Generalized Flow Cover Cut Generator
// author:   Yan Xu                email: Yan.Xu@sas.com
//           Jeff Linderoth        email: jtl3@lehigh.edu
//           Martin Savelsberg     email: martin.savelsbergh@isye.gatech.edu
// date:     05/01/2003
// comments: please scan this file for '???' and read the comments
//-----------------------------------------------------------------------------
// Copyright (C) 2003, Yan Xu, Jeff Linderoth, Martin Savelsberg and others. 
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>
#include <iostream>

//#include "CoinPackedMatrix.hpp"
#include "CglFlowCover.hpp"

//--------------------------------------------------------------------------

void
CglFlowCoverUnitTest(
  const OsiSolverInterface * baseSiP,
  const std::string mpsDir )
{
  // Test default constructor
  {
    CglFlowCover aGenerator;
    assert (aGenerator.getMaxNumCuts() >= 2000);
  }
  
  // Test copy & assignment
  {
    CglFlowCover rhs;
    {
      CglFlowCover bGenerator;
      bGenerator.setMaxNumCuts(100);
      CglFlowCover cGenerator(bGenerator);
      rhs=bGenerator;
    }
  }

  {
    OsiCuts osicuts1;
    CglFlowCover test;

    OsiSolverInterface  * siP = baseSiP->clone();
    
    int nRowCuts;
    std::string fn(mpsDir+"egout");
    std::string fn2 = mpsDir+"egout.mps";
    FILE *in_f = fopen(fn2.c_str(), "r");
    if(in_f == NULL) {
      std::cout<<"Can not open file "<<fn2<<std::endl<<"Skip test of CglFlowCover::generateCuts()"<<std::endl;
    }
    else {
      fclose(in_f);
      siP->readMps(fn.c_str(), "mps");
 
      // Check the preprocess
      test.flowPreprocess(*siP);

#ifdef CGL_DEBUG
      test.printVubs(std::cout);
#endif

      // Test generating cuts
      siP->initialSolve();
      double lpRelax = siP->getObjValue();
      test.generateCuts(*siP, osicuts1);
      nRowCuts = osicuts1.sizeRowCuts();

#ifdef CGL_DEBUG      
      std::cout<<"There are " << nRowCuts << " flow cuts" << std::endl;
      int i;
      for (i = 0; i < nRowCuts; i++){
	OsiRowCut rcut;
	CoinPackedVector rpv;
	const double* colsol = siP->getColSolution();
	rcut = osicuts1.rowCut(i);
	rpv = rcut.row();
	const int n = rpv.getNumElements();
	const int* indices = rpv.getIndices();
	double* elements = rpv.getElements();
	double sum2 = 0.0;
	int k = 0;
	double lb = rcut.lb();
	double ub = rcut.ub();
	for (k = 0; k < n; ++k){
	  int column = indices[k];
	  sum2 += colsol[column] * elements[k];
	}
	if (sum2 > ub + 1.0e-7 ||sum2 < lb - 1.0e-7) {
	  std::cout << "Cut " << i <<" lb " << lb << " solution " << sum2 
		    <<" ub " << ub << std::endl;
	  for (k = 0; k < n; ++k){
	    int column = indices[k];
	    std::cout << "(col=" << column << ", el = " << elements[k] 
		      << ", sol = " << colsol[column] << ") ";
	  }
	  std::cout <<std::endl;
	}
      }
#endif

    // Test generating cuts again
      OsiCuts osicuts2;
      test.generateCuts(*siP, osicuts2);
      OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(osicuts2);
      siP->resolve();
      nRowCuts = osicuts2.sizeRowCuts();
      std::cout<<"There are " << nRowCuts << " flow cuts" << std::endl;

#ifdef CGL_DEBUG
      for (i = 0; i < nRowCuts; i++){
	OsiRowCut rcut;
	CoinPackedVector rpv;
	const double* colsol = siP->getColSolution();
	rcut = osicuts2.rowCut(i);
	rpv = rcut.row();
	const int n = rpv.getNumElements();
	const int* indices = rpv.getIndices();
	double* elements = rpv.getElements();
	double sum2 = 0.0;
	int k = 0;
	double lb = rcut.lb();
	double ub = rcut.ub();
	for (k = 0; k < n; ++k){
	  int column = indices[k];
	  sum2 += colsol[column] * elements[k];
	}
	if (sum2 > ub + 1.0e-7 ||sum2 < lb - 1.0e-7) {
	  std::cout << "Cut " << i <<" lb " << lb << " solution " << sum2 
		    <<" ub " << ub << std::endl;
	  for (k = 0; k < n; ++k){
	    int column = indices[k];
	    std::cout << "(col=" << column << ", el = " << elements[k] 
		      << ", sol = " << colsol[column] << ") ";
	  }
	  std::cout <<std::endl;
	}
      }
#endif
      assert(osicuts2.sizeRowCuts() > 0);
      rc = siP->applyCuts(osicuts2);
      
      siP->resolve();
      
      double lpRelaxAfter= siP->getObjValue(); 
      std::cout<<"Initial LP value: "<<lpRelax<<std::endl;
      std::cout<<"LP value with cuts: "<<lpRelaxAfter<<std::endl;
      assert( lpRelax < lpRelaxAfter );
      assert(lpRelaxAfter < 569);
    }
    delete siP;
  }
}
