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
#include "OsiCuts.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CglGomory.hpp"

void convertColumnCuts(OsiCuts & osicuts)
{
  int nColCuts = osicuts.sizeColCuts();
  if (nColCuts) {
    std::cout<<"There are "<<nColCuts<<" gomory column cuts! - converting"
	     <<std::endl;
    // convert to row cut
    for (int i=0;i<nColCuts;i++) {
      const OsiColCut * cut = osicuts.colCutPtr(i);
      const CoinPackedVector lbs = cut->lbs();
      const CoinPackedVector ubs = cut->ubs();
      int ncols;
      const int * cols;
      const double * values;
      double one = -1.0;
      ncols = lbs.getNumElements();
      cols = lbs.getIndices();
      values = lbs.getElements();
      for (int j=0;j<ncols;j++) {
	int jColumn = cols[j];
	double value = values[j];
	// convert to <=
	OsiRowCut rc;
	rc.setRow(1,&jColumn,&one,false);
	rc.setLb(-COIN_DBL_MAX);
	rc.setUb(-value);
	osicuts.insertIfNotDuplicate(rc);
      }
      one = 1.0;
      ncols = ubs.getNumElements();
      cols = ubs.getIndices();
      values = ubs.getElements();
      for (int j=0;j<ncols;j++) {
	int jColumn = cols[j];
	double value = values[j];
	OsiRowCut rc;
	rc.setRow(1,&jColumn,&one,false);
	rc.setLb(-COIN_DBL_MAX);
	rc.setUb(value);
	osicuts.insertIfNotDuplicate(rc);
      }
    }
  } 
}
//--------------------------------------------------------------------------
// ** At present this does not use any solver
void
CglGomoryUnitTest(
  const OsiSolverInterface * baseSiP,
  const std::string mpsDir )
{
  CoinRelFltEq eq(0.000001);

  // Test default constructor
  {
    CglGomory aGenerator;
    assert (aGenerator.getLimit()==50);
    assert (aGenerator.getAway()==0.05);
  }
  
  // Test copy & assignment etc
  {
    CglGomory rhs;
    {
      CglGomory bGenerator;
      bGenerator.setLimit(99);
      bGenerator.setAway(0.2);
      CglGomory cGenerator(bGenerator);
      rhs=bGenerator;
      assert (rhs.getLimit()==99);
      assert (rhs.getAway()==0.2);
    }
  }

  // Test explicit form - all integer (pg 125 Wolsey)
  if (1) {
    OsiCuts osicuts;
    CglGomory test1;
    int i;
    int nOldCuts=0,nRowCuts;
 
    // matrix data
    //deliberate hiccup of 2 between 0 and 1
    CoinBigIndex start[5]={0,4,7,8,9};
    int length[5]={2,3,1,1,1};
    int rows[11]={0,2,-1,-1,0,1,2,0,1,2};
    double elements[11]={7.0,2.0,1.0e10,1.0e10,-2.0,1.0,-2.0,1,1,1};
    CoinPackedMatrix matrix(true,3,5,8,elements,rows,start,length);
    
    // rim data (objective not used just yet)
    double rowLower[5]={14.0,3.0,3.0,1.0e10,1.0e10};
    double rowUpper[5]={14.0,3.0,3.0,-1.0e10,-1.0e10};
    double colLower[7]={0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    double colUpper[7]={100.0,100.0,100.0,100.0,100.0,100.0,100.0};
  
    // integer
    char intVar[7]={2,2,2,2,2,2,2};

    // basis 1
    int rowBasis1[3]={-1,-1,-1};
    int colBasis1[5]={1,1,-1,-1,1};
    CoinWarmStartBasis warm;
    warm.setSize(5,3);
    for (i=0;i<3;i++) {
      if (rowBasis1[i]<0) {
	warm.setArtifStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setArtifStatus(i,CoinWarmStartBasis::basic);
      }
    }
    for (i=0;i<5;i++) {
      if (colBasis1[i]<0) {
	warm.setStructStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setStructStatus(i,CoinWarmStartBasis::basic);
      }
    }

    // solution 1
    double colsol1[5]={20.0/7.0,3.0,0.0,0.0,23.0/7.0};
    test1.generateCuts(NULL, osicuts, matrix,
		       /*objective,*/ colsol1,
		 colLower, colUpper,
		 rowLower, rowUpper, intVar, &warm);
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" gomory cuts"<<std::endl;
    assert (nRowCuts==2);
    // cuts always <=
    int testCut=0; // test first cut as stronger
    double rhs=-6.0;
    double testCut1[5]={0.0,0.0,-1.0,-2.0,0.0};
    double * cut = testCut1;
    double * colsol = colsol1;
    for (i=nOldCuts; i<nRowCuts; i++){
      OsiRowCut rcut;
      CoinPackedVector rpv;
      rcut = osicuts.rowCut(i);
      rpv = rcut.row();
      const int n = rpv.getNumElements();
      const int * indices = rpv.getIndices();
      double* elements = rpv.getElements();
      double sum2=0.0;
      int k=0;
      for (k=0; k<n; k++){
	int column=indices[k];
	sum2 += colsol[column]*elements[k];
      }

      double ub=rcut.ub();

#ifdef CGL_DEBUG
      double lb=rcut.lb();
      if (sum2 >ub + 1.0e-7 ||sum2 < lb - 1.0e-7) {
	std::cout<<"Cut "<<i<<" lb "<<lb<<" solution "<<sum2<<" ub "<<ub<<std::endl;
	for (k=0; k<n; k++){
	  int column=indices[k];
	  std::cout<<"(col="<<column<<",el="<<elements[k]<<",sol="<<
	    colsol[column]<<") ";
	}
	std::cout <<std::endl;
      }
#endif

      if (i-nOldCuts==testCut) {
	assert( eq(rhs,ub));
	assert(n==2);
	for (k=0; k<n; k++){
	  int column=indices[k];
	  assert (eq(cut[column],elements[k]));
	}
	// add cut
	// explicit slack
	matrix.setDimensions(-1,6);
	rpv.insert(5,1.0*7.0); // to get cut in book
	rowLower[3]=ub;
	rowUpper[3]=ub;
	matrix.appendRow(rpv);
      }
    }
    nOldCuts=nRowCuts;
    // basis 2
    int rowBasis2[4]={-1,-1,-1,-1};
    int colBasis2[6]={1,1,1,1,-1,-1};
    warm.setSize(6,4);
    for (i=0;i<4;i++) {
      if (rowBasis2[i]<0) {
	warm.setArtifStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setArtifStatus(i,CoinWarmStartBasis::basic);
      }
    }
    for (i=0;i<6;i++) {
      if (colBasis2[i]<0) {
	warm.setStructStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setStructStatus(i,CoinWarmStartBasis::basic);
      }
    }

    // solution 2
    double colsol2[6]={2.0,0.5,1.0,2.5,0.0,0.0};
    test1.generateCuts(NULL, osicuts, matrix,
		       /*objective,*/ colsol2,
		 colLower, colUpper,
		 rowLower, rowUpper, intVar, &warm);
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" gomory cuts"<<std::endl;
    //assert (nRowCuts-nOldCuts==2);
    // cuts always <=
    testCut=0; // test first cut as stronger
    rhs=-1.0;
    double testCut2[6]={0.0,0.0,0.0,0.0,-1.0,0.0};
    cut = testCut2;
    colsol = colsol2;
    for (i=nOldCuts; i<nRowCuts; i++){
      OsiRowCut rcut;
      CoinPackedVector rpv;
      rcut = osicuts.rowCut(i);
      rpv = rcut.row();
      const int n = rpv.getNumElements();
      const int * indices = rpv.getIndices();
      double* elements = rpv.getElements();
      double sum2=0.0;
      int k=0;
      for (k=0; k<n; k++){
	int column=indices[k];
	sum2 += colsol[column]*elements[k];
      }

      double ub=rcut.ub();

#ifdef CGL_DEBUG
      double lb=rcut.lb();
      if (sum2 >ub + 1.0e-7 ||sum2 < lb - 1.0e-7) {
	std::cout<<"Cut "<<i<<" lb "<<lb<<" solution "<<sum2<<" ub "<<ub<<std::endl;
	for (k=0; k<n; k++){
	  int column=indices[k];
	  std::cout<<"(col="<<column<<",el="<<elements[k]<<",sol="<<
	    colsol[column]<<") ";
	}
	std::cout <<std::endl;
      }
#endif

      if (i-nOldCuts==testCut) {
	assert( eq(rhs,ub));
	assert(n==1);
	for (k=0; k<n; k++){
	  int column=indices[k];
	  assert (eq(cut[column],elements[k]));
	}
	// add cut
	// explicit slack
	matrix.setDimensions(-1,7);
	rpv.insert(6,1.0);
	rowLower[4]=ub;
	rowUpper[4]=ub;
	matrix.appendRow(rpv);
      }
    }
    nOldCuts=nRowCuts;
    // basis 3
    int rowBasis3[5]={-1,-1,-1,-1,-1};
    int colBasis3[7]={1,1,1,1,1,-1,-1};
    warm.setSize(7,5);
    for (i=0;i<5;i++) {
      if (rowBasis3[i]<0) {
	warm.setArtifStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setArtifStatus(i,CoinWarmStartBasis::basic);
      }
    }
    for (i=0;i<7;i++) {
      if (colBasis3[i]<0) {
	warm.setStructStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setStructStatus(i,CoinWarmStartBasis::basic);
      }
    }

    // solution 3
    double colsol3[7]={2.0,1.0,2.0,2.0,1.0,0.0,0.0};
    test1.generateCuts(NULL, osicuts, matrix,
		       /*objective,*/ colsol3,
		 colLower, colUpper,
		 rowLower, rowUpper, intVar, &warm);
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" gomory cuts"<<std::endl;
    assert (nRowCuts==nOldCuts);
    
  }
  // Test explicit form - this time with x4 flipped
  if (1) {
    OsiCuts osicuts;
    CglGomory test1;
    int i;
    int nOldCuts=0,nRowCuts;
 
    // matrix data
    //deliberate hiccup of 2 between 0 and 1
    CoinBigIndex start[5]={0,4,7,8,9};
    int length[5]={2,3,1,1,1};
    int rows[11]={0,2,-1,-1,0,1,2,0,1,2};
    double elements[11]={7.0,2.0,1.0e10,1.0e10,-2.0,1.0,-2.0,1,-1,1};
    CoinPackedMatrix matrix(true,3,5,8,elements,rows,start,length);
    
    // rim data (objective not used just yet)
    double rowLower[5]={14.0,-5.0,3.0,1.0e10,1.0e10};
    double rowUpper[5]={14.0,-5.0,3.0,-1.0e10,-1.0e10};
    double colLower[7]={0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    double colUpper[7]={100.0,100.0,100.0,8.0,100.0,100.0,100.0};
  
    // integer
    char intVar[7]={2,2,2,2,2,2,2};

    // basis 1
    int rowBasis1[3]={-1,-1,-1};
    int colBasis1[5]={1,1,-1,-1,1};
    CoinWarmStartBasis warm;
    warm.setSize(5,3);
    for (i=0;i<3;i++) {
      if (rowBasis1[i]<0) {
	warm.setArtifStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setArtifStatus(i,CoinWarmStartBasis::basic);
      }
    }
    for (i=0;i<5;i++) {
      if (colBasis1[i]<0) {
	warm.setStructStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setStructStatus(i,CoinWarmStartBasis::basic);
      }
    }

    // solution 1
    double colsol1[5]={20.0/7.0,3.0,0.0,8.0,23.0/7.0};
    test1.generateCuts(NULL, osicuts, matrix,
		       /*objective,*/ colsol1,
		 colLower, colUpper,
		 rowLower, rowUpper, intVar, &warm);
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" gomory cuts"<<std::endl;
    assert (nRowCuts==2);
    // cuts always <=
    int testCut=0; // test first cut as stronger
    double rhs=10.0;
    double testCut1[5]={0.0,0.0,-1.0,2.0,0.0};
    double * cut = testCut1;
    double * colsol = colsol1;
    for (i=nOldCuts; i<nRowCuts; i++){
      OsiRowCut rcut;
      CoinPackedVector rpv;
      rcut = osicuts.rowCut(i);
      rpv = rcut.row();
      const int n = rpv.getNumElements();
      const int * indices = rpv.getIndices();
      double* elements = rpv.getElements();
      double sum2=0.0;
      int k=0;
      for (k=0; k<n; k++){
	int column=indices[k];
	sum2 += colsol[column]*elements[k];
      }

      double ub=rcut.ub();

#ifdef CGL_DEBUG
      double lb=rcut.lb();
      if (sum2 >ub + 1.0e-7 ||sum2 < lb - 1.0e-7) {
	std::cout<<"Cut "<<i<<" lb "<<lb<<" solution "<<sum2<<" ub "<<ub<<std::endl;
	for (k=0; k<n; k++){
	  int column=indices[k];
	  std::cout<<"(col="<<column<<",el="<<elements[k]<<",sol="<<
	    colsol[column]<<") ";
	}
	std::cout <<std::endl;
      }
#endif

      if (i-nOldCuts==testCut) {
	assert( eq(rhs,ub));
	assert(n==2);
	for (k=0; k<n; k++){
	  int column=indices[k];
	  assert (eq(cut[column],elements[k]));
	}
	// add cut
	// explicit slack
	matrix.setDimensions(-1,6);
	rpv.insert(5,1.0*7.0); // to get cut in book
	rowLower[3]=ub;
	rowUpper[3]=ub;
	matrix.appendRow(rpv);
      }
    }
    nOldCuts=nRowCuts;
    // basis 2
    int rowBasis2[4]={-1,-1,-1,-1};
    int colBasis2[6]={1,1,1,1,-1,-1};
    warm.setSize(6,4);
    for (i=0;i<4;i++) {
      if (rowBasis2[i]<0) {
	warm.setArtifStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setArtifStatus(i,CoinWarmStartBasis::basic);
      }
    }
    for (i=0;i<6;i++) {
      if (colBasis2[i]<0) {
	warm.setStructStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setStructStatus(i,CoinWarmStartBasis::basic);
      }
    }

    // solution 2
    double colsol2[6]={2.0,0.5,1.0,5.5,0.0,0.0};
    test1.generateCuts(NULL, osicuts, matrix,
		       /*objective,*/ colsol2,
		 colLower, colUpper,
		 rowLower, rowUpper, intVar, &warm);
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" gomory cuts"<<std::endl;
    //assert (nRowCuts-nOldCuts==2);
    // cuts always <=
    testCut=0; // test first cut as stronger
    rhs=-1.0;
    double testCut2[6]={0.0,0.0,0.0,0.0,-1.0,0.0};
    cut = testCut2;
    colsol = colsol2;
    for (i=nOldCuts; i<nRowCuts; i++){
      OsiRowCut rcut;
      CoinPackedVector rpv;
      rcut = osicuts.rowCut(i);
      rpv = rcut.row();
      const int n = rpv.getNumElements();
      const int * indices = rpv.getIndices();
      double* elements = rpv.getElements();
      double sum2=0.0;
      int k=0;
      for (k=0; k<n; k++){
	int column=indices[k];
	sum2 += colsol[column]*elements[k];
      }

      double ub=rcut.ub();

#ifdef CGL_DEBUG
      double lb=rcut.lb();
      if (sum2 >ub + 1.0e-7 ||sum2 < lb - 1.0e-7) {
	std::cout<<"Cut "<<i<<" lb "<<lb<<" solution "<<sum2<<" ub "<<ub<<std::endl;
	for (k=0; k<n; k++){
	  int column=indices[k];
	  std::cout<<"(col="<<column<<",el="<<elements[k]<<",sol="<<
	    colsol[column]<<") ";
	}
	std::cout <<std::endl;
      }
#endif

      if (i-nOldCuts==testCut) {
	assert( eq(rhs,ub));
	assert(n==1);
	for (k=0; k<n; k++){
	  int column=indices[k];
	  assert (eq(cut[column],elements[k]));
	}
	// add cut
	// explicit slack
	matrix.setDimensions(-1,7);
	rpv.insert(6,1.0);
	rowLower[4]=ub;
	rowUpper[4]=ub;
	matrix.appendRow(rpv);
      }
    }
    nOldCuts=nRowCuts;
    // basis 3
    int rowBasis3[5]={-1,-1,-1,-1,-1};
    int colBasis3[7]={1,1,1,1,1,-1,-1};
    warm.setSize(7,5);
    for (i=0;i<5;i++) {
      if (rowBasis3[i]<0) {
	warm.setArtifStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setArtifStatus(i,CoinWarmStartBasis::basic);
      }
    }
    for (i=0;i<7;i++) {
      if (colBasis3[i]<0) {
	warm.setStructStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setStructStatus(i,CoinWarmStartBasis::basic);
      }
    }

    // solution 3
    double colsol3[7]={2.0,1.0,2.0,6.0,1.0,0.0,0.0};
    test1.generateCuts(NULL, osicuts, matrix,
		       /*objective,*/ colsol3,
		 colLower, colUpper,
		 rowLower, rowUpper, intVar, &warm);
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" gomory cuts"<<std::endl;
    assert (nRowCuts==nOldCuts);
    
  }
  // Test with slacks 
  if (1) {
    OsiCuts osicuts;
    CglGomory test1;
    int i;
    int nOldCuts=0,nRowCuts;
 
    // matrix data
    //deliberate hiccup of 2 between 0 and 1
    CoinBigIndex start[5]={0,4};
    int length[5]={2,3};
    int rows[11]={0,2,-1,-1,0,1,2};
    double elements[11]={7.0,2.0,1.0e10,1.0e10,-2.0,1.0,-2.0};
    CoinPackedMatrix matrix(true,3,2,5,elements,rows,start,length);
    
    // rim data (objective not used just yet)
    double rowLower[5]={-1.0e10,-1.0e10,-1.0e10,1.0e10,1.0e10};
    double rowUpper[5]={14.0,3.0,3.0,-1.0e10,-1.0e10};
    double colLower[2]={0.0,0.0};
    double colUpper[2]={100.0,100.0};
  
    // integer
    char intVar[2]={2,2};

    // basis 1
    int rowBasis1[3]={-1,-1,1};
    int colBasis1[2]={1,1};
    CoinWarmStartBasis warm;
    warm.setSize(2,3);
    for (i=0;i<3;i++) {
      if (rowBasis1[i]<0) {
	warm.setArtifStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setArtifStatus(i,CoinWarmStartBasis::basic);
      }
    }
    for (i=0;i<2;i++) {
      if (colBasis1[i]<0) {
	warm.setStructStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setStructStatus(i,CoinWarmStartBasis::basic);
      }
    }

    // solution 1
    double colsol1[2]={20.0/7.0,3.0};
    test1.generateCuts(NULL, osicuts, matrix,
		       /* objective,*/ colsol1,
		 colLower, colUpper,
		 rowLower, rowUpper, intVar, &warm);
    // new version may create a column cut if just one element
    convertColumnCuts(osicuts);
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" gomory cuts"<<std::endl;
    assert (nRowCuts==1);
    // cuts always <=
    int testCut=0; // test first cut as stronger
    double rhs=2.0;
    double testCut1[2]={1.0,0.0};
    double * cut = testCut1;
    double * colsol = colsol1;
    for (i=nOldCuts; i<nRowCuts; i++){
      OsiRowCut rcut;
      CoinPackedVector rpv;
      rcut = osicuts.rowCut(i);
      rpv = rcut.row();
      const int n = rpv.getNumElements();
      const int * indices = rpv.getIndices();
      double* elements = rpv.getElements();
      double sum2=0.0;
      int k=0;
      for (k=0; k<n; k++){
	int column=indices[k];
	sum2 += colsol[column]*elements[k];
      }

      double ub=rcut.ub();

#ifdef CGL_DEBUG
      double lb=rcut.lb();
      if (sum2 >ub + 1.0e-7 ||sum2 < lb - 1.0e-7) {
	std::cout<<"Cut "<<i<<" lb "<<lb<<" solution "<<sum2<<" ub "<<ub<<std::endl;
	for (k=0; k<n; k++){
	  int column=indices[k];
	  std::cout<<"(col="<<column<<",el="<<elements[k]<<",sol="<<
	    colsol[column]<<") ";
	}
	std::cout <<std::endl;
      }
#endif

      if (i-nOldCuts==testCut) {
	assert( eq(rhs,ub));
	assert(n==1);
	for (k=0; k<n; k++){
	  int column=indices[k];
	  assert (eq(cut[column],elements[k]));
	}
	// add cut
	rowLower[3]=-1.0e100;
	rowUpper[3]=ub;
	matrix.appendRow(rpv);
      }
    }
    nOldCuts=nRowCuts;
    // basis 2
    int rowBasis2[4]={1,1,-1,-1};
    int colBasis2[2]={1,1};
    warm.setSize(2,4);
    for (i=0;i<4;i++) {
      if (rowBasis2[i]<0) {
	warm.setArtifStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setArtifStatus(i,CoinWarmStartBasis::basic);
      }
    }
    for (i=0;i<2;i++) {
      if (colBasis2[i]<0) {
	warm.setStructStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setStructStatus(i,CoinWarmStartBasis::basic);
      }
    }

    // solution 2
    double colsol2[2]={2.0,0.5};
    test1.generateCuts(NULL, osicuts, matrix,
		       /*objective,*/ colsol2,
		 colLower, colUpper,
		 rowLower, rowUpper, intVar, &warm);
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" gomory cuts"<<std::endl;
    assert (nRowCuts-nOldCuts==1);
    // cuts always <=
    testCut=0; // test first cut as stronger
    rhs=1.0;
    double testCut2[2]={1.0,-1.0};
    cut = testCut2;
    colsol = colsol2;
    for (i=nOldCuts; i<nRowCuts; i++){
      OsiRowCut rcut;
      CoinPackedVector rpv;
      rcut = osicuts.rowCut(i);
      rpv = rcut.row();
      const int n = rpv.getNumElements();
      const int * indices = rpv.getIndices();
      double* elements = rpv.getElements();
      double sum2=0.0;
      int k=0;
      for (k=0; k<n; k++){
	int column=indices[k];
	sum2 += colsol[column]*elements[k];
      }

      double ub=rcut.ub();

#ifdef CGL_DEBUG
      double lb=rcut.lb();
      if (sum2 >ub + 1.0e-7 ||sum2 < lb - 1.0e-7) {
	std::cout<<"Cut "<<i<<" lb "<<lb<<" solution "<<sum2<<" ub "<<ub<<std::endl;
	for (k=0; k<n; k++){
	  int column=indices[k];
	  std::cout<<"(col="<<column<<",el="<<elements[k]<<",sol="<<
	    colsol[column]<<") ";
	}
	std::cout <<std::endl;
      }
#endif

      if (i-nOldCuts==testCut) {
	assert( eq(rhs,ub));
	assert(n==2);
	for (k=0; k<n; k++){
	  int column=indices[k];
	  assert (eq(cut[column],elements[k]));
	}
	// add cut
	rowLower[4]=-1.0e100;
	rowUpper[4]=ub;
	matrix.appendRow(rpv);
      }
    }
    nOldCuts=nRowCuts;
    // basis 3
    int rowBasis3[5]={1,1,1,-1,-1};
    int colBasis3[2]={1,1};
    warm.setSize(2,5);
    for (i=0;i<5;i++) {
      if (rowBasis3[i]<0) {
	warm.setArtifStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setArtifStatus(i,CoinWarmStartBasis::basic);
      }
    }
    for (i=0;i<2;i++) {
      if (colBasis3[i]<0) {
	warm.setStructStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setStructStatus(i,CoinWarmStartBasis::basic);
      }
    }

    // solution 3
    double colsol3[2]={2.0,1.0};
    test1.generateCuts(NULL, osicuts, matrix,
		       /*objective,*/ colsol3,
		 colLower, colUpper,
		 rowLower, rowUpper, intVar, &warm);
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" gomory cuts"<<std::endl;
    assert (nRowCuts==nOldCuts);
    
  }
  // swap some rows to G
  if (1) {
    OsiCuts osicuts;
    CglGomory test1;
    int i;
    int nOldCuts=0,nRowCuts;
 
    // matrix data
    //deliberate hiccup of 2 between 0 and 1
    CoinBigIndex start[5]={0,4};
    int length[5]={2,3};
    int rows[11]={0,2,-1,-1,0,1,2};
    double elements[11]={-7.0,-2.0,1.0e10,1.0e10,+2.0,1.0,+2.0};
    CoinPackedMatrix matrix(true,3,2,5,elements,rows,start,length);
    
    // rim data (objective not used just yet)
    double rowUpper[5]={1.0e10,3.0,1.0e10,-1.0e10,-1.0e10};
    double rowLower[5]={-14.0,-1.0e10,-3.0,1.0e10,1.0e10};
    double colLower[2]={0.0,0.0};
    double colUpper[2]={100.0,100.0};
  
    // integer
    char intVar[2]={2,2};

    // basis 1
    int rowBasis1[3]={-1,-1,1};
    int colBasis1[2]={1,1};
    CoinWarmStartBasis warm;
    warm.setSize(2,3);
    for (i=0;i<3;i++) {
      if (rowBasis1[i]<0) {
	warm.setArtifStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setArtifStatus(i,CoinWarmStartBasis::basic);
      }
    }
    for (i=0;i<2;i++) {
      if (colBasis1[i]<0) {
	warm.setStructStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setStructStatus(i,CoinWarmStartBasis::basic);
      }
    }

    // solution 1
    double colsol1[2]={20.0/7.0,3.0};
    test1.generateCuts(NULL, osicuts, matrix,
		 /*objective,*/ colsol1,
		 colLower, colUpper,
		 rowLower, rowUpper, intVar, &warm);
    convertColumnCuts(osicuts);
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" gomory cuts"<<std::endl;
    assert (nRowCuts==1);
    // cuts always <=
    int testCut=0; // test first cut as stronger
    double rhs=2.0;
    double testCut1[2]={1.0,0.0};
    double * cut = testCut1;
    double * colsol = colsol1;
    for (i=nOldCuts; i<nRowCuts; i++){
      OsiRowCut rcut;
      CoinPackedVector rpv;
      rcut = osicuts.rowCut(i);
      rpv = rcut.row();
      const int n = rpv.getNumElements();
      const int * indices = rpv.getIndices();
      double* elements = rpv.getElements();
      double sum2=0.0;
      int k=0;
      for (k=0; k<n; k++){
	int column=indices[k];
	sum2 += colsol[column]*elements[k];
      }

      double ub=rcut.ub();

#ifdef CGL_DEBUG
      double lb=rcut.lb();
      if (sum2 >ub + 1.0e-7 ||sum2 < lb - 1.0e-7) {
	std::cout<<"Cut "<<i<<" lb "<<lb<<" solution "<<sum2<<" ub "<<ub<<std::endl;
	for (k=0; k<n; k++){
	  int column=indices[k];
	  std::cout<<"(col="<<column<<",el="<<elements[k]<<",sol="<<
	    colsol[column]<<") ";
	}
	std::cout <<std::endl;
      }
#endif

      if (i-nOldCuts==testCut) {
	assert( eq(rhs,ub));
	assert(n==1);
	for (k=0; k<n; k++){
	  int column=indices[k];
	  assert (eq(cut[column],elements[k]));
	}
	// add cut
	rowLower[3]=-1.0e100;
	rowUpper[3]=ub;
	matrix.appendRow(rpv);
      }
    }
    nOldCuts=nRowCuts;
    // basis 2
    int rowBasis2[4]={1,1,-1,-1};
    int colBasis2[2]={1,1};
    warm.setSize(2,4);
    for (i=0;i<4;i++) {
      if (rowBasis2[i]<0) {
	warm.setArtifStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setArtifStatus(i,CoinWarmStartBasis::basic);
      }
    }
    for (i=0;i<2;i++) {
      if (colBasis2[i]<0) {
	warm.setStructStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setStructStatus(i,CoinWarmStartBasis::basic);
      }
    }

    // solution 2
    double colsol2[2]={2.0,0.5};
    test1.generateCuts(NULL, osicuts, matrix,
		 /*objective,*/ colsol2,
		 colLower, colUpper,
		 rowLower, rowUpper, intVar, &warm);
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" gomory cuts"<<std::endl;
    assert (nRowCuts-nOldCuts==1);
    // cuts always <=
    testCut=0; // test first cut as stronger
    rhs=1.0;
    double testCut2[2]={1.0,-1.0};
    cut = testCut2;
    colsol = colsol2;
    for (i=nOldCuts; i<nRowCuts; i++){
      OsiRowCut rcut;
      CoinPackedVector rpv;
      rcut = osicuts.rowCut(i);
      rpv = rcut.row();
      const int n = rpv.getNumElements();
      const int * indices = rpv.getIndices();
      double* elements = rpv.getElements();
      double sum2=0.0;
      int k=0;
      for (k=0; k<n; k++){
	int column=indices[k];
	sum2 += colsol[column]*elements[k];
      }

      double ub=rcut.ub();

#ifdef CGL_DEBUG
      double lb=rcut.lb();
      if (sum2 >ub + 1.0e-7 ||sum2 < lb - 1.0e-7) {
	std::cout<<"Cut "<<i<<" lb "<<lb<<" solution "<<sum2<<" ub "<<ub<<std::endl;
	for (k=0; k<n; k++){
	  int column=indices[k];
	  std::cout<<"(col="<<column<<",el="<<elements[k]<<",sol="<<
	    colsol[column]<<") ";
	}
	std::cout <<std::endl;
      }
#endif

      if (i-nOldCuts==testCut) {
	assert( eq(rhs,ub));
	assert(n==2);
	for (k=0; k<n; k++){
	  int column=indices[k];
	  assert (eq(cut[column],elements[k]));
	}
	// add cut
	rowLower[4]=-1.0e100;
	rowUpper[4]=ub;
	matrix.appendRow(rpv);
      }
    }
    nOldCuts=nRowCuts;
    // basis 3
    int rowBasis3[5]={1,1,1,-1,-1};
    int colBasis3[2]={1,1};
    warm.setSize(2,5);
    for (i=0;i<5;i++) {
      if (rowBasis3[i]<0) {
	warm.setArtifStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setArtifStatus(i,CoinWarmStartBasis::basic);
      }
    }
    for (i=0;i<2;i++) {
      if (colBasis3[i]<0) {
	warm.setStructStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setStructStatus(i,CoinWarmStartBasis::basic);
      }
    }

    // solution 3
    double colsol3[2]={2.0,1.0};
    test1.generateCuts(NULL, osicuts, matrix,
		 /*objective,*/ colsol3,
		 colLower, colUpper,
		 rowLower, rowUpper, intVar, &warm);
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" gomory cuts"<<std::endl;
    assert (nRowCuts==nOldCuts);
    
  }


  // NOW mixed integer gomory cuts

  // Test explicit form - (pg 130 Wolsey)
  // Some arrays left same size as previously although not used in full
  if (1) {
    OsiCuts osicuts;
    CglGomory test1;
    int i;
    int nOldCuts=0,nRowCuts;
 
    // matrix data
    //deliberate hiccup of 2 between 0 and 1
    CoinBigIndex start[5]={0,4,7,8,9};
    int length[5]={2,3,1,1,1};
    int rows[11]={0,2,-1,-1,0,1,2,0,1,2};
    double elements[11]={7.0,2.0,1.0e10,1.0e10,-2.0,1.0,-2.0,1,1,1};
    CoinPackedMatrix matrix(true,3,5,8,elements,rows,start,length);
    
    // rim data (objective not used just yet)
    double rowLower[5]={14.0,3.0,3.0,1.0e10,1.0e10};
    double rowUpper[5]={14.0,3.0,3.0,-1.0e10,-1.0e10};
    double colLower[7]={0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    double colUpper[7]={100.0,100.0,100.0,100.0,100.0,100.0,100.0};
  
    // integer
    char intVar[7]={2,0,0,0,0,0,0};

    // basis 1
    int rowBasis1[3]={-1,-1,-1};
    int colBasis1[5]={1,1,-1,-1,1};
    CoinWarmStartBasis warm;
    warm.setSize(5,3);
    for (i=0;i<3;i++) {
      if (rowBasis1[i]<0) {
	warm.setArtifStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setArtifStatus(i,CoinWarmStartBasis::basic);
      }
    }
    for (i=0;i<5;i++) {
      if (colBasis1[i]<0) {
	warm.setStructStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setStructStatus(i,CoinWarmStartBasis::basic);
      }
    }

    // solution 1
    double colsol1[5]={20.0/7.0,3.0,0.0,0.0,23.0/7.0};
    test1.generateCuts(NULL, osicuts, matrix,
		 /*objective,*/ colsol1,
		 colLower, colUpper,
		 rowLower, rowUpper, intVar, &warm);
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" gomory cuts"<<std::endl;
    assert (nRowCuts==1);
    // cuts always <=
    int testCut=0; // test first cut as stronger
    double rhs=-6.0/7.0;
    double testCut1[5]={0.0,0.0,-1.0/7.0,-2.0/7.0,0.0};
    double * cut = testCut1;
    double * colsol = colsol1;
    for (i=nOldCuts; i<nRowCuts; i++){
      OsiRowCut rcut;
      CoinPackedVector rpv;
      rcut = osicuts.rowCut(i);
      rpv = rcut.row();
      const int n = rpv.getNumElements();
      const int * indices = rpv.getIndices();
      double* elements = rpv.getElements();
      double sum2=0.0;
      int k=0;
      for (k=0; k<n; k++){
	int column=indices[k];
	sum2 += colsol[column]*elements[k];
      }

      double ub=rcut.ub();

#ifdef CGL_DEBUG
      double lb=rcut.lb();
      if (sum2 >ub + 1.0e-7 ||sum2 < lb - 1.0e-7) {
	std::cout<<"Cut "<<i<<" lb "<<lb<<" solution "<<sum2<<" ub "<<ub<<std::endl;
	for (k=0; k<n; k++){
	  int column=indices[k];
	  std::cout<<"(col="<<column<<",el="<<elements[k]<<",sol="<<
	    colsol[column]<<") ";
	}
	std::cout <<std::endl;
      }
#endif

      if (i-nOldCuts==testCut) {
	assert( rhs*ub>=0.0);
	assert(n==2);
	double ratio = cut[indices[0]]/elements[0];
	assert (eq(ub*ratio,rhs));
	for (k=0; k<n; k++){
	  int column=indices[k];
	  assert (eq(cut[column],ratio*elements[k]));
	}
	// add cut
	// explicit slack
	matrix.setDimensions(-1,6);
	rpv.insert(5,1.0); // to get cut in book
	rowLower[3]=ub;
	rowUpper[3]=ub;
	matrix.appendRow(rpv);
      }
    }
    nOldCuts=nRowCuts;
    // basis 2
    int rowBasis2[4]={-1,-1,-1,-1};
    int colBasis2[6]={1,1,1,1,-1,-1};
    warm.setSize(6,4);
    for (i=0;i<4;i++) {
      if (rowBasis2[i]<0) {
	warm.setArtifStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setArtifStatus(i,CoinWarmStartBasis::basic);
      }
    }
    for (i=0;i<6;i++) {
      if (colBasis2[i]<0) {
	warm.setStructStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setStructStatus(i,CoinWarmStartBasis::basic);
      }
    }

    // solution 2
    double colsol2[6]={2.0,0.5,1.0,2.5,0.0,0.0};
    test1.generateCuts(NULL, osicuts, matrix,
		 /*objective,*/ colsol2,
		 colLower, colUpper,
		 rowLower, rowUpper, intVar, &warm);
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" gomory cuts"<<std::endl;
    assert (nRowCuts==nOldCuts);
    
  }
  // Test explicit form - this time with x4 flipped 
  if (1) {
    OsiCuts osicuts;
    CglGomory test1;
    int i;
    int nOldCuts=0,nRowCuts;
 
    // matrix data
    //deliberate hiccup of 2 between 0 and 1
    CoinBigIndex start[5]={0,4,7,8,9};
    int length[5]={2,3,1,1,1};
    int rows[11]={0,2,-1,-1,0,1,2,0,1,2};
    double elements[11]={7.0,2.0,1.0e10,1.0e10,-2.0,1.0,-2.0,1,-1,1};
    CoinPackedMatrix matrix(true,3,5,8,elements,rows,start,length);
    
    // rim data (objective not used just yet)
    double rowLower[5]={14.0,-5.0,3.0,1.0e10,1.0e10};
    double rowUpper[5]={14.0,-5.0,3.0,-1.0e10,-1.0e10};
    double colLower[7]={0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    double colUpper[7]={100.0,100.0,100.0,8.0,100.0,100.0,100.0};
  
    // integer
    char intVar[7]={2,0,0,0,0,0,0};

    // basis 1
    int rowBasis1[3]={-1,-1,-1};
    int colBasis1[5]={1,1,-1,-1,1};
    CoinWarmStartBasis warm;
    warm.setSize(5,3);
    for (i=0;i<3;i++) {
      if (rowBasis1[i]<0) {
	warm.setArtifStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setArtifStatus(i,CoinWarmStartBasis::basic);
      }
    }
    for (i=0;i<5;i++) {
      if (colBasis1[i]<0) {
	warm.setStructStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setStructStatus(i,CoinWarmStartBasis::basic);
      }
    }

    // solution 1
    double colsol1[5]={20.0/7.0,3.0,0.0,8.0,23.0/7.0};
    test1.generateCuts(NULL, osicuts, matrix,
		 /*objective,*/ colsol1,
		 colLower, colUpper,
		 rowLower, rowUpper, intVar, &warm);
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" gomory cuts"<<std::endl;
    assert (nRowCuts==1);
    // cuts always <=
    int testCut=0; 
    double rhs=10.0/7.0;
    double testCut1[5]={0.0,0.0,-1.0/7.0,2.0/7.0,0.0};
    double * cut = testCut1;
    double * colsol = colsol1;
    for (i=nOldCuts; i<nRowCuts; i++){
      OsiRowCut rcut;
      CoinPackedVector rpv;
      rcut = osicuts.rowCut(i);
      rpv = rcut.row();
      const int n = rpv.getNumElements();
      const int * indices = rpv.getIndices();
      double* elements = rpv.getElements();
      double sum2=0.0;
      int k=0;
      for (k=0; k<n; k++){
	int column=indices[k];
	sum2 += colsol[column]*elements[k];
      }

      double ub=rcut.ub();

#ifdef CGL_DEBUG
      double lb=rcut.lb();
      if (sum2 >ub + 1.0e-7 ||sum2 < lb - 1.0e-7) {
	std::cout<<"Cut "<<i<<" lb "<<lb<<" solution "<<sum2<<" ub "<<ub<<std::endl;
	for (k=0; k<n; k++){
	  int column=indices[k];
	  std::cout<<"(col="<<column<<",el="<<elements[k]<<",sol="<<
	    colsol[column]<<") ";
	}
	std::cout <<std::endl;
      }
#endif

      if (i-nOldCuts==testCut) {
	assert( rhs*ub>=0.0);
	assert(n==2);
	double ratio = cut[indices[0]]/elements[0];
	assert (eq(ub*ratio,rhs));
	for (k=0; k<n; k++){
	  int column=indices[k];
	  assert (eq(cut[column],ratio*elements[k]));
	}
	// add cut
	// explicit slack
	matrix.setDimensions(-1,6);
	rpv.insert(5,1.0); // to get cut in book
	rowLower[3]=ub;
	rowUpper[3]=ub;
	matrix.appendRow(rpv);
      }
    }
    nOldCuts=nRowCuts;
    // basis 2
    int rowBasis2[4]={-1,-1,-1,-1};
    int colBasis2[6]={1,1,1,1,-1,-1};
    warm.setSize(6,4);
    for (i=0;i<4;i++) {
      if (rowBasis2[i]<0) {
	warm.setArtifStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setArtifStatus(i,CoinWarmStartBasis::basic);
      }
    }
    for (i=0;i<6;i++) {
      if (colBasis2[i]<0) {
	warm.setStructStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setStructStatus(i,CoinWarmStartBasis::basic);
      }
    }

    // solution 2
    double colsol2[6]={2.0,0.5,1.0,5.5,0.0,0.0};
    test1.generateCuts(NULL, osicuts, matrix,
		 /*objective,*/ colsol2,
		 colLower, colUpper,
		 rowLower, rowUpper, intVar, &warm);
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" gomory cuts"<<std::endl;
    assert (nRowCuts==nOldCuts);
    
  }
  // Test with slacks 
  if (1) {
    OsiCuts osicuts;
    CglGomory test1;
    int i;
    int nOldCuts=0,nRowCuts;
 
    // matrix data
    //deliberate hiccup of 2 between 0 and 1
    CoinBigIndex start[5]={0,4};
    int length[5]={2,3};
    int rows[11]={0,2,-1,-1,0,1,2};
    double elements[11]={7.0,2.0,1.0e10,1.0e10,-2.0,1.0,-2.0};
    CoinPackedMatrix matrix(true,3,2,5,elements,rows,start,length);
    
    // rim data (objective not used just yet)
    double rowLower[5]={-1.0e10,-1.0e10,-1.0e10,1.0e10,1.0e10};
    double rowUpper[5]={14.0,3.0,3.0,-1.0e10,-1.0e10};
    double colLower[2]={0.0,0.0};
    double colUpper[2]={100.0,100.0};
  
    // integer
    char intVar[2]={2,0};

    // basis 1
    int rowBasis1[3]={-1,-1,1};
    int colBasis1[2]={1,1};
    CoinWarmStartBasis warm;
    warm.setSize(2,3);
    for (i=0;i<3;i++) {
      if (rowBasis1[i]<0) {
	warm.setArtifStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setArtifStatus(i,CoinWarmStartBasis::basic);
      }
    }
    for (i=0;i<2;i++) {
      if (colBasis1[i]<0) {
	warm.setStructStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setStructStatus(i,CoinWarmStartBasis::basic);
      }
    }

    // solution 1
    double colsol1[2]={20.0/7.0,3.0};
    test1.generateCuts(NULL, osicuts, matrix,
		 /*objective,*/ colsol1,
		 colLower, colUpper,
		 rowLower, rowUpper, intVar, &warm);
    convertColumnCuts(osicuts);
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" gomory cuts"<<std::endl;
    assert (nRowCuts==1);
    // cuts always <=
    int testCut=0; // test first cut as stronger
    double rhs=2.0;
    double testCut1[2]={1.0,0.0};
    double * cut = testCut1;
    double * colsol = colsol1;
    for (i=nOldCuts; i<nRowCuts; i++){
      OsiRowCut rcut;
      CoinPackedVector rpv;
      rcut = osicuts.rowCut(i);
      rpv = rcut.row();
      const int n = rpv.getNumElements();
      const int * indices = rpv.getIndices();
      double* elements = rpv.getElements();
      double sum2=0.0;
      int k=0;
      for (k=0; k<n; k++){
	int column=indices[k];
	sum2 += colsol[column]*elements[k];
      }

      double ub=rcut.ub();

#ifdef CGL_DEBUG
      double lb=rcut.lb();
      if (sum2 >ub + 1.0e-7 ||sum2 < lb - 1.0e-7) {
	std::cout<<"Cut "<<i<<" lb "<<lb<<" solution "<<sum2<<" ub "<<ub<<std::endl;
	for (k=0; k<n; k++){
	  int column=indices[k];
	  std::cout<<"(col="<<column<<",el="<<elements[k]<<",sol="<<
	    colsol[column]<<") ";
	}
	std::cout <<std::endl;
      }
#endif

      if (i-nOldCuts==testCut) {
	assert( eq(rhs,ub));
	assert(n==1);
	for (k=0; k<n; k++){
	  int column=indices[k];
	  assert (eq(cut[column],elements[k]));
	}
	// add cut
	rowLower[3]=-1.0e100;
	rowUpper[3]=ub;
	matrix.appendRow(rpv);
      }
    }
    nOldCuts=nRowCuts;
    // basis 2
    int rowBasis2[4]={1,1,-1,-1};
    int colBasis2[2]={1,1};
    warm.setSize(2,4);
    for (i=0;i<4;i++) {
      if (rowBasis2[i]<0) {
	warm.setArtifStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setArtifStatus(i,CoinWarmStartBasis::basic);
      }
    }
    for (i=0;i<2;i++) {
      if (colBasis2[i]<0) {
	warm.setStructStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setStructStatus(i,CoinWarmStartBasis::basic);
      }
    }

    // solution 2
    double colsol2[2]={2.0,0.5};
    test1.generateCuts(NULL, osicuts, matrix,
		       /*objective,*/ colsol2,
		 colLower, colUpper,
		 rowLower, rowUpper, intVar, &warm);
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" gomory cuts"<<std::endl;
    assert (nRowCuts==nOldCuts);
    
  }
  // swap some rows to G
  if (1) {
    OsiCuts osicuts;
    CglGomory test1;
    int i;
    int nOldCuts=0,nRowCuts;
 
    // matrix data
    //deliberate hiccup of 2 between 0 and 1
    CoinBigIndex start[5]={0,4};
    int length[5]={2,3};
    int rows[11]={0,2,-1,-1,0,1,2};
    double elements[11]={-7.0,-2.0,1.0e10,1.0e10,+2.0,1.0,+2.0};
    CoinPackedMatrix matrix(true,3,2,5,elements,rows,start,length);
    
    // rim data (objective not used just yet)
    double rowUpper[5]={1.0e10,3.0,1.0e10,-1.0e10,-1.0e10};
    double rowLower[5]={-14.0,-1.0e10,-3.0,1.0e10,1.0e10};
    double colLower[2]={0.0,0.0};
    double colUpper[2]={100.0,100.0};
  
    // integer
    char intVar[2]={2,0};

    // basis 1
    int rowBasis1[3]={-1,-1,1};
    int colBasis1[2]={1,1};
    CoinWarmStartBasis warm;
    warm.setSize(2,3);
    for (i=0;i<3;i++) {
      if (rowBasis1[i]<0) {
	warm.setArtifStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setArtifStatus(i,CoinWarmStartBasis::basic);
      }
    }
    for (i=0;i<2;i++) {
      if (colBasis1[i]<0) {
	warm.setStructStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setStructStatus(i,CoinWarmStartBasis::basic);
      }
    }

    // solution 1
    double colsol1[2]={20.0/7.0,3.0};
    test1.generateCuts(NULL, osicuts, matrix,
		       /*objective,*/ colsol1,
		 colLower, colUpper,
		 rowLower, rowUpper, intVar, &warm);
    convertColumnCuts(osicuts);
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" gomory cuts"<<std::endl;
    assert (nRowCuts==1);
    // cuts always <=
    int testCut=0; // test first cut as stronger
    double rhs=2.0;
    double testCut1[2]={1.0,0.0};
    double * cut = testCut1;
    double * colsol = colsol1;
    for (i=nOldCuts; i<nRowCuts; i++){
      OsiRowCut rcut;
      CoinPackedVector rpv;
      rcut = osicuts.rowCut(i);
      rpv = rcut.row();
      const int n = rpv.getNumElements();
      const int * indices = rpv.getIndices();
      double* elements = rpv.getElements();
      double sum2=0.0;
      int k=0;
      for (k=0; k<n; k++){
	int column=indices[k];
	sum2 += colsol[column]*elements[k];
      }

      double ub=rcut.ub();

#ifdef CGL_DEBUG
      double lb=rcut.lb();
      if (sum2 >ub + 1.0e-7 ||sum2 < lb - 1.0e-7) {
	std::cout<<"Cut "<<i<<" lb "<<lb<<" solution "<<sum2<<" ub "<<ub<<std::endl;
	for (k=0; k<n; k++){
	  int column=indices[k];
	  std::cout<<"(col="<<column<<",el="<<elements[k]<<",sol="<<
	    colsol[column]<<") ";
	}
	std::cout <<std::endl;
      }
#endif

      if (i-nOldCuts==testCut) {
	assert( eq(rhs,ub));
	assert(n==1);
	for (k=0; k<n; k++){
	  int column=indices[k];
	  assert (eq(cut[column],elements[k]));
	}
	// add cut
	rowLower[3]=-1.0e100;
	rowUpper[3]=ub;
	matrix.appendRow(rpv);
      }
    }
    nOldCuts=nRowCuts;
    // basis 2
    int rowBasis2[4]={1,1,-1,-1};
    int colBasis2[2]={1,1};
    warm.setSize(2,4);
    for (i=0;i<4;i++) {
      if (rowBasis2[i]<0) {
	warm.setArtifStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setArtifStatus(i,CoinWarmStartBasis::basic);
      }
    }
    for (i=0;i<2;i++) {
      if (colBasis2[i]<0) {
	warm.setStructStatus(i,CoinWarmStartBasis::atLowerBound);
      } else {
	warm.setStructStatus(i,CoinWarmStartBasis::basic);
      }
    }

    // solution 2
    double colsol2[2]={2.0,0.5};
    test1.generateCuts(NULL, osicuts, matrix,
		 /*objective,*/ colsol2,
		 colLower, colUpper,
		 rowLower, rowUpper, intVar, &warm);
    nRowCuts = osicuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" gomory cuts"<<std::endl;
    assert (nRowCuts==nOldCuts);
    
  }

  // Miplib3 problem p0033
  if (1) {
    // Setup
    OsiSolverInterface  * siP = baseSiP->clone();
    std::string fn(mpsDir+"p0033");
    siP->readMps(fn.c_str(),"mps");
    siP->activateRowCutDebugger("p0033");
    CglGomory test;

    // Solve the LP relaxation of the model and
    // print out ofv for sake of comparison 
    siP->initialSolve();
    double lpRelaxBefore=siP->getObjValue();
    std::cout<<"Initial LP value: "<<lpRelaxBefore<<std::endl;
    assert( eq(lpRelaxBefore, 2520.5717391304347) );

    // Fails with OsiCpx, OsiXpr:
    /**********
    double mycs[] = {0, 1, 0, 0, -2.0837010502455788e-19, 1, 0, 0, 1,
		       0.021739130434782594, 0.35652173913043478, 
		       -6.7220534694101275e-18, 5.3125906451789717e-18, 
		       1, 0, 1.9298798670241979e-17, 0, 0, 0,
		       7.8875708048320448e-18, 0.5, 0, 
		       0.85999999999999999, 1, 1, 0.57999999999999996,
		       1, 0, 1, 0, 0.25, 0, 0.67500000000000004};
    siP->setColSolution(mycs);
    ****/

    OsiCuts cuts;    
    
    // Test generateCuts method
    test.generateCuts(*siP,cuts);
    int nRowCuts = cuts.sizeRowCuts();
    std::cout<<"There are "<<nRowCuts<<" Gomory cuts"<<std::endl;
    assert(cuts.sizeRowCuts() > 0);
    OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cuts);
    
    siP->resolve();
    double lpRelaxAfter=siP->getObjValue(); 
    std::cout<<"LP value with cuts: "<<lpRelaxAfter<<std::endl;
    //assert( eq(lpRelaxAfter, 2592.1908295194507) );
    assert( lpRelaxAfter> 2545.0 );
    assert( lpRelaxBefore < lpRelaxAfter );
    assert(lpRelaxAfter < 3089.1);
    
    delete siP;
  } 
}

