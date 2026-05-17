// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdio>

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

#include "CoinPragma.hpp"
#include "CoinPackedMatrix.hpp"
#include "CglProbing.hpp"


//--------------------------------------------------------------------------
// test EKKsolution methods.
void
CglProbingUnitTest(
  const OsiSolverInterface * baseSiP,
  const std::string mpsDir )
{
# ifdef CGL_DEBUG
  int i ;	// define just once
# endif
  CoinRelFltEq eq(0.000001);

  // Test default constructor
  {
    CglProbing aGenerator;
  }
  
  // Test copy & assignment
  {
    CglProbing rhs;
    {
      CglProbing bGenerator;
      CglProbing cGenerator(bGenerator);
      rhs=bGenerator;
    }
  }

  {
    OsiCuts osicuts;
    CglProbing test1;
    OsiSolverInterface  * siP = baseSiP->clone();
    int nColCuts;
    int nRowCuts;
    
    std::string fn = mpsDir+"p0033";
    siP->readMps(fn.c_str(),"mps");
    siP->initialSolve();
    // just unsatisfied variables
    test1.generateCuts(*siP,osicuts);
    nColCuts = osicuts.sizeColCuts();
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" probing cuts"<<std::endl;
    {
      std::cout<<"there are "<<nColCuts<<" probing column cuts"<<std::endl;

#ifdef CGL_DEBUG
      const double * lo = siP->getColLower();
      const double * up = siP->getColUpper();
      for (i=0; i<nColCuts; i++){
	OsiColCut ccut;
	CoinPackedVector cpv;
	ccut = osicuts.colCut(i);
	cpv = ccut.lbs();
	int n = cpv.getNumElements();
        int j;
	const int * indices = cpv.getIndices();
	double* elements = cpv.getElements();
	for (j=0;j<n;j++) {
	  int icol=indices[j];
	  if (elements[j]>lo[icol])
	    std::cout<<"Can increase lb on "<<icol<<" from "<<lo[icol]<<
	      " to "<<elements[j]<<std::endl;
	}
	cpv = ccut.ubs();
	n = cpv.getNumElements();
	indices = cpv.getIndices();
	elements = cpv.getElements();

	for (j=0;j<n;j++) {
	  int icol=indices[j];
	  if (elements[j]<up[icol])
	    std::cout<<"Can decrease ub on "<<icol<<" from "<<up[icol]<<
	      " to "<<elements[j]<<std::endl;
	}
      }

#endif

    }

#ifdef CGL_DEBUG
    for (i=0; i<nRowCuts; i++){
      OsiRowCut rcut;
      CoinPackedVector rpv;
      const double * colsol = siP->getColSolution();
      rcut = osicuts.rowCut(i);
      rpv = rcut.row();
      const int n = rpv.getNumElements();
      const int * indices = rpv.getIndices();
      double* elements = rpv.getElements();
      double sum2=0.0;
      int k=0;
      double lb=rcut.lb();
      double ub=rcut.ub();
      for (k=0; k<n; k++){
	int column=indices[k];
	sum2 += colsol[column]*elements[k];
      }
      if (sum2 >ub + 1.0e-7 ||sum2 < lb - 1.0e-7) {
	std::cout<<"Cut "<<i<<" lb "<<lb<<" solution "<<sum2<<" ub "<<ub<<std::endl;
	for (k=0; k<n; k++){
	  int column=indices[k];
	  std::cout<<"(col="<<column<<",el="<<elements[k]<<",sol="<<
	    colsol[column]<<") ";
	}
	std::cout <<std::endl;
      }
    }
#endif

    if (nRowCuts==1) {
      CoinPackedVector check;
      int index[] = {6,32};
      double el[] = {1,1};
      check.setVector(2,index,el);
      // sort Elements in increasing order
      CoinPackedVector rpv=osicuts.rowCut(0).row();
      assert (rpv.getNumElements()==2);
      rpv.sortIncrIndex();
      assert (check==rpv);
      assert (osicuts.rowCut(0).lb()==1.0);
    }
    // now all variables
    osicuts=OsiCuts();
    test1.setMode(2);
    test1.setRowCuts(3);
    test1.generateCuts(*siP,osicuts);
    nColCuts = osicuts.sizeColCuts();
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" probing cuts"<<std::endl;
    {
      std::cout<<"there are "<<nColCuts<<" probing column cuts"<<std::endl;

#ifdef CGL_DEBUG
      const double * lo = siP->getColLower();
      const double * up = siP->getColUpper();
      for (i=0; i<nColCuts; i++){
	OsiColCut ccut;
	CoinPackedVector cpv;
	ccut = osicuts.colCut(i);
	cpv = ccut.lbs();
	int n = cpv.getNumElements();
        int j;
	const int * indices = cpv.getIndices();
	double* elements = cpv.getElements();
	for (j=0;j<n;j++) {
	  int icol=indices[j];
	  if (elements[j]>lo[icol])
	    std::cout<<"Can increase lb on "<<icol<<" from "<<lo[icol]<<
	      " to "<<elements[j]<<std::endl;
	}
	cpv = ccut.ubs();
	n = cpv.getNumElements();
	indices = cpv.getIndices();
	elements = cpv.getElements();
	for (j=0;j<n;j++) {
	  int icol=indices[j];
	  if (elements[j]<up[icol])
	    std::cout<<"Can decrease ub on "<<icol<<" from "<<up[icol]<<
	      " to "<<elements[j]<<std::endl;
	}
      }
#endif

    }

#ifdef CGL_DEBUG
    for (i=0; i<nRowCuts; i++){
      OsiRowCut rcut;
      CoinPackedVector rpv;
      const double * colsol = siP->getColSolution();
      rcut = osicuts.rowCut(i);
      rpv = rcut.row();
      const int n = rpv.getNumElements();
      const int * indices = rpv.getIndices();
      double* elements = rpv.getElements();
      double sum2=0.0;
      int k=0;
      double lb=rcut.lb();
      double ub=rcut.ub();
      for (k=0; k<n; k++){
	int column=indices[k];
	sum2 += colsol[column]*elements[k];
      }
      if (sum2 >ub + 1.0e-7 ||sum2 < lb - 1.0e-7) {
	std::cout<<"Cut "<<i<<" lb "<<lb<<" solution "<<sum2<<" ub "<<ub<<std::endl;
	for (k=0; k<n; k++){
	  int column=indices[k];
	  std::cout<<"(col="<<column<<",el="<<elements[k]<<",sol="<<
	    colsol[column]<<") ";
	}
	std::cout <<std::endl;
      }
    }
#endif

    assert (osicuts.sizeRowCuts()>=4);
    delete siP;
  }

}

