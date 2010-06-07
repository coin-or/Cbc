// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#ifdef NDEBUG
#undef NDEBUG
#endif

#include "OsiConfig.h"

#include <cassert>
//#include <cstdlib>
//#include <cstdio>
//#include <iostream>

#include "OsiCbcSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "CoinMessage.hpp"
#include "CoinModel.hpp"
#ifdef COIN_HAS_OSL
#include "OsiOslSolverInterface.hpp"
#endif

//#############################################################################

namespace {
CoinPackedMatrix &BuildExmip1Mtx ()
/*
  Simple function to build a packed matrix for the exmip1 example used in
  tests. The function exists solely to hide the intermediate variables.
  Probably could be written as an initialised declaration.
  See COIN/Mps/Sample/exmip1.mps for a human-readable presentation.

  Ordered triples seem easiest. They're listed in row-major order.
*/

{ int rowndxs[] = { 0, 0, 0, 0, 0,
		    1, 1,
		    2, 2,
		    3, 3,
		    4, 4, 4 } ;
  int colndxs[] = { 0, 1, 3, 4, 7,
		    1, 2,
		    2, 5,
		    3, 6,
		    0, 4, 7 } ;
  double coeffs[] = { 3.0, 1.0, -2.0, -1.0, -1.0,
		      2.0, 1.1,
		      1.0, 1.0,
		      2.8, -1.2,
		      5.6, 1.0, 1.9 } ;

  static CoinPackedMatrix exmip1mtx =
    CoinPackedMatrix(true,&rowndxs[0],&colndxs[0],&coeffs[0],14) ;

  return (exmip1mtx) ; }
}

//--------------------------------------------------------------------------
// test solution methods.
int
OsiCbcSolverInterfaceUnitTest(const std::string & mpsDir, const std::string & netlibDir)
{
  
  
  {    
    CoinRelFltEq eq;
    OsiCbcSolverInterface m;
    std::string fn = mpsDir+"exmip1";
    m.readMps(fn.c_str(),"mps");

    {
      OsiCbcSolverInterface im;    
      
      assert( im.getNumCols() == 0 ); 
      
      assert( im.getModelPtr()!=NULL );
    }
    
    // Test copy constructor and assignment operator
    {
      OsiCbcSolverInterface lhs;
      {      
        OsiCbcSolverInterface im(m);        
        
        OsiCbcSolverInterface imC1(im);
        assert( imC1.getModelPtr()!=im.getModelPtr() );
        assert( imC1.getNumCols() == im.getNumCols() );
        assert( imC1.getNumRows() == im.getNumRows() );   
        
        OsiCbcSolverInterface imC2(im);
        assert( imC2.getModelPtr()!=im.getModelPtr() );
        assert( imC2.getNumCols() == im.getNumCols() );
        assert( imC2.getNumRows() == im.getNumRows() );  
        
        assert( imC2.getModelPtr()!=imC1.getModelPtr() );
        
        lhs=imC2;
      }
      // Test that lhs has correct values even though rhs has gone out of scope
      
      assert( lhs.getModelPtr() != m.getModelPtr() );
      assert( lhs.getNumCols() == m.getNumCols() );
      assert( lhs.getNumRows() == m.getNumRows() );      
    }
    // Test clone
    {
      OsiCbcSolverInterface oslSi(m);
      OsiSolverInterface * siPtr = &oslSi;
      OsiSolverInterface * siClone = siPtr->clone();
      OsiCbcSolverInterface * oslClone = dynamic_cast<OsiCbcSolverInterface*>(siClone);
      assert( oslClone != NULL );
      assert( oslClone->getModelPtr() != oslSi.getModelPtr() );
      assert( oslClone->getNumRows() == oslSi.getNumRows() );
      assert( oslClone->getNumCols() == m.getNumCols() );
      
      delete siClone;
    }
  
    // test infinity
    {
      OsiCbcSolverInterface si;
      assert( eq(si.getInfinity(),OsiCbcInfinity));
    }     
    
    // Test setting solution
    {
      OsiCbcSolverInterface m1(m);
      int i;

      double * cs = new double[m1.getNumCols()];
      for ( i = 0;  i < m1.getNumCols();  i++ ) 
        cs[i] = i + .5;
      m1.setColSolution(cs);
      for ( i = 0;  i < m1.getNumCols();  i++ ) 
        assert(m1.getColSolution()[i] == i + .5);
      
      double * rs = new double[m1.getNumRows()];
      for ( i = 0;  i < m1.getNumRows();  i++ ) 
        rs[i] = i - .5;
      m1.setRowPrice(rs);
      for ( i = 0;  i < m1.getNumRows();  i++ ) 
        assert(m1.getRowPrice()[i] == i - .5);

      delete [] cs;
      delete [] rs;
    }
    
    
    // Test fraction Indices
    {
      OsiCbcSolverInterface fim;
      std::string fn = mpsDir+"exmip1";
      fim.readMps(fn.c_str(),"mps");
      // exmip1.mps has 2 integer variables with index 2 & 3
      assert(  fim.isContinuous(0) );
      assert(  fim.isContinuous(1) );
      assert( !fim.isContinuous(2) );
      assert( !fim.isContinuous(3) );
      assert(  fim.isContinuous(4) );
      
      assert( !fim.isInteger(0) );
      assert( !fim.isInteger(1) );
      assert(  fim.isInteger(2) );
      assert(  fim.isInteger(3) );
      assert( !fim.isInteger(4) );
      
      assert( !fim.isBinary(0) );
      assert( !fim.isBinary(1) );
      assert(  fim.isBinary(2) );
      assert(  fim.isBinary(3) );
      assert( !fim.isBinary(4) );
      
      assert( !fim.isIntegerNonBinary(0) );
      assert( !fim.isIntegerNonBinary(1) );
      assert( !fim.isIntegerNonBinary(2) );
      assert( !fim.isIntegerNonBinary(3) );
      assert( !fim.isIntegerNonBinary(4) );

    }
    // Test some catches
      if (!OsiCbcHasNDEBUG())
    { bool thrown ;

      thrown = false ;
      OsiCbcSolverInterface solver;
      try {
        solver.setObjCoeff(0,0.0);
      }
      catch (CoinError e) {
        std::cout<<"Correct throw"<<std::endl;
	thrown = true ;
      }
      assert( thrown == true ) ;

      std::string fn = mpsDir+"exmip1";
      solver.readMps(fn.c_str(),"mps");
      try {
        solver.setObjCoeff(0,0.0);
      }
      catch (CoinError e) {
        std::cout<<"** Incorrect throw"<<std::endl;
        abort();
      }

      thrown = false ;
      try {
        int index[]={0,20};
        double value[]={0.0,0.0,0.0,0.0};
        solver.setColSetBounds(index,index+2,value);
      }
      catch (CoinError e) {
        std::cout<<"Correct throw"<<std::endl;
	thrown = true ;
      }
      assert( thrown == true ) ;
    }
    // Test apply cuts method
    {      
      OsiCbcSolverInterface im(m);
      OsiCuts cuts;
      
      // Generate some cuts 
      {
        // Get number of rows and columns in model
        int nr=im.getNumRows();
        int nc=im.getNumCols();
        assert( nr == 5 );
        assert( nc == 8 );
        
        // Generate a valid row cut from thin air
        int c;
        {
          int *inx = new int[nc];
          for (c=0;c<nc;c++) inx[c]=c;
          double *el = new double[nc];
          for (c=0;c<nc;c++) el[c]=((double)c)*((double)c);
          
          OsiRowCut rc;
          rc.setRow(nc,inx,el);
          rc.setLb(-100.);
          rc.setUb(100.);
          rc.setEffectiveness(22);
          
          cuts.insert(rc);
          delete[]el;
          delete[]inx;
        }
        
        // Generate valid col cut from thin air
        {
          const double * oslColLB = im.getColLower();
          const double * oslColUB = im.getColUpper();
          int *inx = new int[nc];
          for (c=0;c<nc;c++) inx[c]=c;
          double *lb = new double[nc];
          double *ub = new double[nc];
          for (c=0;c<nc;c++) lb[c]=oslColLB[c]+0.001;
          for (c=0;c<nc;c++) ub[c]=oslColUB[c]-0.001;
          
          OsiColCut cc;
          cc.setLbs(nc,inx,lb);
          cc.setUbs(nc,inx,ub);
          
          cuts.insert(cc);
          delete [] ub;
          delete [] lb;
          delete [] inx;
        }
        
        {
          // Generate a row and column cut which are ineffective
          OsiRowCut * rcP= new OsiRowCut;
          rcP->setEffectiveness(-1.);
          cuts.insert(rcP);
          assert(rcP==NULL);
          
          OsiColCut * ccP= new OsiColCut;
          ccP->setEffectiveness(-12.);
          cuts.insert(ccP);
          assert(ccP==NULL);
        }
        {
          //Generate inconsistent Row cut
          OsiRowCut rc;
          const int ne=1;
          int inx[ne]={-10};
          double el[ne]={2.5};
          rc.setRow(ne,inx,el);
          rc.setLb(3.);
          rc.setUb(4.);
          assert(!rc.consistent());
          cuts.insert(rc);
        }
        {
          //Generate inconsistent col cut
          OsiColCut cc;
          const int ne=1;
          int inx[ne]={-10};
          double el[ne]={2.5};
          cc.setUbs(ne,inx,el);
          assert(!cc.consistent());
          cuts.insert(cc);
        }
        {
          // Generate row cut which is inconsistent for model m
          OsiRowCut rc;
          const int ne=1;
          int inx[ne]={10};
          double el[ne]={2.5};
          rc.setRow(ne,inx,el);
          assert(rc.consistent());
          assert(!rc.consistent(im));
          cuts.insert(rc);
        }
        {
          // Generate col cut which is inconsistent for model m
          OsiColCut cc;
          const int ne=1;
          int inx[ne]={30};
          double el[ne]={2.0};
          cc.setLbs(ne,inx,el);
          assert(cc.consistent());
          assert(!cc.consistent(im));
          cuts.insert(cc);
        }
        {
          // Generate col cut which is infeasible
          OsiColCut cc;
          const int ne=1;
          int inx[ne]={0};
          double el[ne]={2.0};
          cc.setUbs(ne,inx,el);
          cc.setEffectiveness(1000.);
          assert(cc.consistent());
          assert(cc.consistent(im));
          assert(cc.infeasible(im));
          cuts.insert(cc);
        }
      }
      assert(cuts.sizeRowCuts()==4);
      assert(cuts.sizeColCuts()==5);
      
      OsiSolverInterface::ApplyCutsReturnCode rc = im.applyCuts(cuts);
      assert( rc.getNumIneffective() == 2 );
      assert( rc.getNumApplied() == 2 );
      assert( rc.getNumInfeasible() == 1 );
      assert( rc.getNumInconsistentWrtIntegerModel() == 2 );
      assert( rc.getNumInconsistent() == 2 );
      assert( cuts.sizeCuts() == rc.getNumIneffective() +
        rc.getNumApplied() +
        rc.getNumInfeasible() +
        rc.getNumInconsistentWrtIntegerModel() +
        rc.getNumInconsistent() );
    }
    
    {    
      OsiCbcSolverInterface oslSi(m);
      int nc = oslSi.getNumCols();
      int nr = oslSi.getNumRows();
      const double * cl = oslSi.getColLower();
      const double * cu = oslSi.getColUpper();
      const double * rl = oslSi.getRowLower();
      const double * ru = oslSi.getRowUpper();
      assert( nc == 8 );
      assert( nr == 5 );
      assert( eq(cl[0],2.5) );
      assert( eq(cl[1],0.0) );
      assert( eq(cu[1],4.1) );
      assert( eq(cu[2],1.0) );
      assert( eq(rl[0],2.5) );
      assert( eq(rl[4],3.0) );
      assert( eq(ru[1],2.1) );
      assert( eq(ru[4],15.0) );
      
      const double * cs = oslSi.getColSolution();
      assert( eq(cs[0],2.5) );
      assert( eq(cs[7],0.0) );
      
      assert( !eq(cl[3],1.2345) );
      oslSi.setColLower( 3, 1.2345 );
      assert( eq(oslSi.getColLower()[3],1.2345) );
      
      assert( !eq(cu[4],10.2345) );
      oslSi.setColUpper( 4, 10.2345 );
      assert( eq(oslSi.getColUpper()[4],10.2345) );
      // LH: Objective will depend on how underlying solver constructs and
      // LH: maintains initial solution.
      double objValue = oslSi.getObjValue();
      assert( eq(objValue,3.5) || eq(objValue,10.5) );

      assert( eq( oslSi.getObjCoefficients()[0],  1.0) );
      assert( eq( oslSi.getObjCoefficients()[1],  0.0) );
      assert( eq( oslSi.getObjCoefficients()[2],  0.0) );
      assert( eq( oslSi.getObjCoefficients()[3],  0.0) );
      assert( eq( oslSi.getObjCoefficients()[4],  2.0) );
      assert( eq( oslSi.getObjCoefficients()[5],  0.0) );
      assert( eq( oslSi.getObjCoefficients()[6],  0.0) );
      assert( eq( oslSi.getObjCoefficients()[7], -1.0) );
    }
    
    // Test matrixByRow method
    { 
      const OsiCbcSolverInterface si(m);
      const CoinPackedMatrix * smP = si.getMatrixByRow();
      // LL:      const OsiCbcPackedMatrix * osmP = dynamic_cast<const OsiCbcPackedMatrix*>(smP);
      // LL: assert( osmP!=NULL );

#ifdef OSICBC_TEST_MTX_STRUCTURE

      CoinRelFltEq eq;
      const double * ev = smP->getElements();
      assert( eq(ev[0],   3.0) );
      assert( eq(ev[1],   1.0) );
      assert( eq(ev[2],  -2.0) );
      assert( eq(ev[3],  -1.0) );
      assert( eq(ev[4],  -1.0) );
      assert( eq(ev[5],   2.0) );
      assert( eq(ev[6],   1.1) );
      assert( eq(ev[7],   1.0) );
      assert( eq(ev[8],   1.0) );
      assert( eq(ev[9],   2.8) );
      assert( eq(ev[10], -1.2) );
      assert( eq(ev[11],  5.6) );
      assert( eq(ev[12],  1.0) );
      assert( eq(ev[13],  1.9) );
      
      const CoinBigIndex * mi = smP->getVectorStarts();
      assert( mi[0]==0 );
      assert( mi[1]==5 );
      assert( mi[2]==7 );
      assert( mi[3]==9 );
      assert( mi[4]==11 );
      assert( mi[5]==14 );
      
      const int * ei = smP->getIndices();
      assert( ei[0]  ==  0 );
      assert( ei[1]  ==  1 );
      assert( ei[2]  ==  3 );
      assert( ei[3]  ==  4 );
      assert( ei[4]  ==  7 );
      assert( ei[5]  ==  1 );
      assert( ei[6]  ==  2 );
      assert( ei[7]  ==  2 );
      assert( ei[8]  ==  5 );
      assert( ei[9]  ==  3 );
      assert( ei[10] ==  6 );
      assert( ei[11] ==  0 );
      assert( ei[12] ==  4 );
      assert( ei[13] ==  7 );    

#else	// OSICBC_TEST_MTX_STRUCTURE

      CoinPackedMatrix exmip1Mtx ;
      exmip1Mtx.reverseOrderedCopyOf(BuildExmip1Mtx()) ;
      assert( exmip1Mtx.isEquivalent(*smP) ) ;

#endif	// OSICBC_TEST_MTX_STRUCTURE

      assert( smP->getMajorDim() == 5 );
      assert( smP->getMinorDim() == 8 );
      assert( smP->getNumElements() == 14 );
      assert( smP->getSizeVectorStarts() == 6 );
    }

    // Test adding several cuts, and handling of a coefficient of infinity
    // in the constraint matrix.
    {
      OsiCbcSolverInterface fim;
      std::string fn = mpsDir+"exmip1";
      fim.readMps(fn.c_str(),"mps");
      // exmip1.mps has 2 integer variables with index 2 & 3
      fim.initialSolve();
      OsiRowCut cuts[3];
      
      
      // Generate one ineffective cut plus two trivial cuts
      int c;
      int nc = fim.getNumCols();
      int *inx = new int[nc];
      for (c=0;c<nc;c++) inx[c]=c;
      double *el = new double[nc];
      for (c=0;c<nc;c++) el[c]=1.0e-50+((double)c)*((double)c);
      
      cuts[0].setRow(nc,inx,el);
      cuts[0].setLb(-100.);
      cuts[0].setUb(500.);
      cuts[0].setEffectiveness(22);
      el[4]=0.0; // to get inf later
      
      for (c=2;c<4;c++) {
        el[0]=1.0;
        inx[0]=c;
        cuts[c-1].setRow(1,inx,el);
        cuts[c-1].setLb(1.);
        cuts[c-1].setUb(100.);
        cuts[c-1].setEffectiveness(c);
      }
      fim.writeMps("x1.mps");
      fim.applyRowCuts(3,cuts);
      fim.writeMps("x2.mps");
      // resolve - should get message about zero elements
      fim.resolve();
      fim.writeMps("x3.mps");
      // check integer solution
      const double * cs = fim.getColSolution();
      CoinRelFltEq eq;
      assert( eq(cs[2],   1.0) );
      assert( eq(cs[3],   1.0) );
      // check will find invalid matrix
      el[0]=1.0/el[4];
      inx[0]=0;
      cuts[0].setRow(nc,inx,el);
      cuts[0].setLb(-100.);
      cuts[0].setUb(500.);
      cuts[0].setEffectiveness(22);
      fim.applyRowCut(cuts[0]);
      // resolve - should get message about zero elements
      fim.resolve();
      assert (fim.isAbandoned());
      delete[]el;
      delete[]inx;
    }

        // Test matrixByCol method
    {
  
      const OsiCbcSolverInterface si(m);
      const CoinPackedMatrix * smP = si.getMatrixByCol();
      // LL:      const OsiCbcPackedMatrix * osmP = dynamic_cast<const OsiCbcPackedMatrix*>(smP);
      // LL: assert( osmP!=NULL );
      
#ifdef OSICBC_TEST_MTX_STRUCTURE

      CoinRelFltEq eq;
      const double * ev = smP->getElements();
      assert( eq(ev[0],   3.0) );
      assert( eq(ev[1],   5.6) );
      assert( eq(ev[2],   1.0) );
      assert( eq(ev[3],   2.0) );
      assert( eq(ev[4],   1.1) );
      assert( eq(ev[5],   1.0) );
      assert( eq(ev[6],  -2.0) );
      assert( eq(ev[7],   2.8) );
      assert( eq(ev[8],  -1.0) );
      assert( eq(ev[9],   1.0) );
      assert( eq(ev[10],  1.0) );
      assert( eq(ev[11], -1.2) );
      assert( eq(ev[12], -1.0) );
      assert( eq(ev[13],  1.9) );
      
      const CoinBigIndex * mi = smP->getVectorStarts();
      assert( mi[0]==0 );
      assert( mi[1]==2 );
      assert( mi[2]==4 );
      assert( mi[3]==6 );
      assert( mi[4]==8 );
      assert( mi[5]==10 );
      assert( mi[6]==11 );
      assert( mi[7]==12 );
      assert( mi[8]==14 );
      
      const int * ei = smP->getIndices();
      assert( ei[0]  ==  0 );
      assert( ei[1]  ==  4 );
      assert( ei[2]  ==  0 );
      assert( ei[3]  ==  1 );
      assert( ei[4]  ==  1 );
      assert( ei[5]  ==  2 );
      assert( ei[6]  ==  0 );
      assert( ei[7]  ==  3 );
      assert( ei[8]  ==  0 );
      assert( ei[9]  ==  4 );
      assert( ei[10] ==  2 );
      assert( ei[11] ==  3 );
      assert( ei[12] ==  0 );
      assert( ei[13] ==  4 );    
      
#else // OSICBC_TEST_MTX_STRUCTURE

      CoinPackedMatrix &exmip1Mtx = BuildExmip1Mtx() ;
      assert( exmip1Mtx.isEquivalent(*smP) ) ;

#endif	// OSICBC_TEST_MTX_STRUCTURE     

      assert( smP->getMajorDim() == 8 ); 
      assert( smP->getMinorDim() == 5 );
      assert( smP->getNumElements() == 14 );
      assert( smP->getSizeVectorStarts() == 9 );
    }

    //--------------
    // Test rowsense, rhs, rowrange, matrixByRow, solver assignment
    {
      OsiCbcSolverInterface lhs;
      {      

        OsiCbcSolverInterface siC1(m);     
        const char   * siC1rs  = siC1.getRowSense();
        assert( siC1rs[0]=='G' );
        assert( siC1rs[1]=='L' );
        assert( siC1rs[2]=='E' );
        assert( siC1rs[3]=='R' );
        assert( siC1rs[4]=='R' );
        
        const double * siC1rhs = siC1.getRightHandSide();
        assert( eq(siC1rhs[0],2.5) );
        assert( eq(siC1rhs[1],2.1) );
        assert( eq(siC1rhs[2],4.0) );
        assert( eq(siC1rhs[3],5.0) );
        assert( eq(siC1rhs[4],15.) ); 
        
        const double * siC1rr  = siC1.getRowRange();
        assert( eq(siC1rr[0],0.0) );
        assert( eq(siC1rr[1],0.0) );
        assert( eq(siC1rr[2],0.0) );
        assert( eq(siC1rr[3],5.0-1.8) );
        assert( eq(siC1rr[4],15.0-3.0) );
        
        const CoinPackedMatrix * siC1mbr = siC1.getMatrixByRow();
        assert( siC1mbr != NULL );

#ifdef OSICBC_TEST_MTX_STRUCTURE

        const double * ev = siC1mbr->getElements();
        assert( eq(ev[0],   3.0) );
        assert( eq(ev[1],   1.0) );
        assert( eq(ev[2],  -2.0) );
        assert( eq(ev[3],  -1.0) );
        assert( eq(ev[4],  -1.0) );
        assert( eq(ev[5],   2.0) );
        assert( eq(ev[6],   1.1) );
        assert( eq(ev[7],   1.0) );
        assert( eq(ev[8],   1.0) );
        assert( eq(ev[9],   2.8) );
        assert( eq(ev[10], -1.2) );
        assert( eq(ev[11],  5.6) );
        assert( eq(ev[12],  1.0) );
        assert( eq(ev[13],  1.9) );
        
        const CoinBigIndex * mi = siC1mbr->getVectorStarts();
        assert( mi[0]==0 );
        assert( mi[1]==5 );
        assert( mi[2]==7 );
        assert( mi[3]==9 );
        assert( mi[4]==11 );
        assert( mi[5]==14 );
        
        const int * ei = siC1mbr->getIndices();
        assert( ei[0]  ==  0 );
        assert( ei[1]  ==  1 );
        assert( ei[2]  ==  3 );
        assert( ei[3]  ==  4 );
        assert( ei[4]  ==  7 );
        assert( ei[5]  ==  1 );
        assert( ei[6]  ==  2 );
        assert( ei[7]  ==  2 );
        assert( ei[8]  ==  5 );
        assert( ei[9]  ==  3 );
        assert( ei[10] ==  6 );
        assert( ei[11] ==  0 );
        assert( ei[12] ==  4 );
        assert( ei[13] ==  7 );    

#else	// OSICBC_TEST_MTX_STRUCTURE

	CoinPackedMatrix exmip1Mtx ;
	exmip1Mtx.reverseOrderedCopyOf(BuildExmip1Mtx()) ;
	assert( exmip1Mtx.isEquivalent(*siC1mbr) ) ;

#endif	// OSICBC_TEST_MTX_STRUCTURE

        assert( siC1mbr->getMajorDim() == 5 ); 
	assert( siC1mbr->getMinorDim() == 8 );
        assert( siC1mbr->getNumElements() == 14 );
	assert( siC1mbr->getSizeVectorStarts()==6 );

        assert( siC1rs  == siC1.getRowSense() );
        assert( siC1rhs == siC1.getRightHandSide() );
        assert( siC1rr  == siC1.getRowRange() );

        // Change OSL Model by adding free row
        OsiRowCut rc;
        rc.setLb(-DBL_MAX);
        rc.setUb( DBL_MAX);
        OsiCuts cuts;
        cuts.insert(rc);
        siC1.applyCuts(cuts);
             
        
        siC1rs  = siC1.getRowSense();
        assert( siC1rs[0]=='G' );
        assert( siC1rs[1]=='L' );
        assert( siC1rs[2]=='E' );
        assert( siC1rs[3]=='R' );
        assert( siC1rs[4]=='R' );
        assert( siC1rs[5]=='N' );

        siC1rhs = siC1.getRightHandSide();
        assert( eq(siC1rhs[0],2.5) );
        assert( eq(siC1rhs[1],2.1) );
        assert( eq(siC1rhs[2],4.0) );
        assert( eq(siC1rhs[3],5.0) );
        assert( eq(siC1rhs[4],15.) ); 
        assert( eq(siC1rhs[5],0.0 ) ); 

        siC1rr  = siC1.getRowRange();
        assert( eq(siC1rr[0],0.0) );
        assert( eq(siC1rr[1],0.0) );
        assert( eq(siC1rr[2],0.0) );
        assert( eq(siC1rr[3],5.0-1.8) );
        assert( eq(siC1rr[4],15.0-3.0) );
        assert( eq(siC1rr[5],0.0) );
    
        lhs=siC1;
      }
      // Test that lhs has correct values even though siC1 has gone out of scope    
      
      const char * lhsrs  = lhs.getRowSense();
      assert( lhsrs[0]=='G' );
      assert( lhsrs[1]=='L' );
      assert( lhsrs[2]=='E' );
      assert( lhsrs[3]=='R' );
      assert( lhsrs[4]=='R' );
      assert( lhsrs[5]=='N' );
      
      const double * lhsrhs = lhs.getRightHandSide();
      assert( eq(lhsrhs[0],2.5) );
      assert( eq(lhsrhs[1],2.1) );
      assert( eq(lhsrhs[2],4.0) );
      assert( eq(lhsrhs[3],5.0) );
      assert( eq(lhsrhs[4],15.) ); 
      assert( eq(lhsrhs[5],0.0) ); 
      
      const double *lhsrr  = lhs.getRowRange();
      assert( eq(lhsrr[0],0.0) );
      assert( eq(lhsrr[1],0.0) );
      assert( eq(lhsrr[2],0.0) );
      assert( eq(lhsrr[3],5.0-1.8) );
      assert( eq(lhsrr[4],15.0-3.0) );
      assert( eq(lhsrr[5],0.0) );      
      
      const CoinPackedMatrix * lhsmbr = lhs.getMatrixByRow();
      assert( lhsmbr != NULL );       

#ifdef OSICBC_TEST_MTX_STRUCTURE

      const double * ev = lhsmbr->getElements();
      assert( eq(ev[0],   3.0) );
      assert( eq(ev[1],   1.0) );
      assert( eq(ev[2],  -2.0) );
      assert( eq(ev[3],  -1.0) );
      assert( eq(ev[4],  -1.0) );
      assert( eq(ev[5],   2.0) );
      assert( eq(ev[6],   1.1) );
      assert( eq(ev[7],   1.0) );
      assert( eq(ev[8],   1.0) );
      assert( eq(ev[9],   2.8) );
      assert( eq(ev[10], -1.2) );
      assert( eq(ev[11],  5.6) );
      assert( eq(ev[12],  1.0) );
      assert( eq(ev[13],  1.9) );
      
      const CoinBigIndex * mi = lhsmbr->getVectorStarts();
      assert( mi[0]==0 );
      assert( mi[1]==5 );
      assert( mi[2]==7 );
      assert( mi[3]==9 );
      assert( mi[4]==11 );
      assert( mi[5]==14 );
      
      const int * ei = lhsmbr->getIndices();
      assert( ei[0]  ==  0 );
      assert( ei[1]  ==  1 );
      assert( ei[2]  ==  3 );
      assert( ei[3]  ==  4 );
      assert( ei[4]  ==  7 );
      assert( ei[5]  ==  1 );
      assert( ei[6]  ==  2 );
      assert( ei[7]  ==  2 );
      assert( ei[8]  ==  5 );
      assert( ei[9]  ==  3 );
      assert( ei[10] ==  6 );
      assert( ei[11] ==  0 );
      assert( ei[12] ==  4 );
      assert( ei[13] ==  7 );    

#else	// OSICBC_TEST_MTX_STRUCTURE

/*
  This admittedly looks bogus, but it's the equivalent operation on the matrix
  for inserting a cut of the form -Inf <= +Inf (i.e., a cut with no
  coefficients).
*/
      CoinPackedMatrix exmip1Mtx ;
      exmip1Mtx.reverseOrderedCopyOf(BuildExmip1Mtx()) ;
      CoinPackedVector freeRow ;
      exmip1Mtx.appendRow(freeRow) ;
      assert( exmip1Mtx.isEquivalent(*lhsmbr) ) ;

#endif	// OSICBC_TEST_MTX_STRUCTURE
      
      int md = lhsmbr->getMajorDim();
      assert(  md == 6 ); 
      assert( lhsmbr->getNumElements() == 14 );
    }
    
  }

  // Test add/delete columns
  {    
    OsiCbcSolverInterface m;
    std::string fn = mpsDir+"p0033";
    m.readMps(fn.c_str(),"mps");
    double inf = m.getInfinity();

    CoinPackedVector c0;
    c0.insert(0, 4);
    c0.insert(1, 1);
    m.addCol(c0, 0, inf, 3);
    m.initialSolve();
    double objValue = m.getObjValue();
    CoinRelFltEq eq(1.0e-2);
    assert( eq(objValue,2520.57) );
    // Try deleting first column that's nonbasic at lower bound (0).
    int * d = new int[1];
    CoinWarmStartBasis *cwsb =
	dynamic_cast<CoinWarmStartBasis *>(m.getWarmStart()) ;
    assert(cwsb) ;
    CoinWarmStartBasis::Status stati ;
    int iCol ;
    for (iCol = 0 ;  iCol < cwsb->getNumStructural() ; iCol++)
    { stati = cwsb->getStructStatus(iCol) ;
      if (stati == CoinWarmStartBasis::atLowerBound) break ; }
    d[0]=iCol;
    m.deleteCols(1,d);
    delete [] d;
    delete cwsb;
    d=NULL;
    m.resolve();
    objValue = m.getObjValue();
    assert( eq(objValue,2520.57) );
    // Try deleting column we added. If basic, go to initialSolve as deleting
    // basic variable trashes basis required for warm start.
    iCol = m.getNumCols()-1;
    cwsb = dynamic_cast<CoinWarmStartBasis *>(m.getWarmStart()) ;
    stati =  cwsb->getStructStatus(iCol) ;
    delete cwsb;
    m.deleteCols(1,&iCol);
    if (stati == CoinWarmStartBasis::basic)
    { m.initialSolve() ; }
    else
    { m.resolve(); }
    objValue = m.getObjValue();
    assert( eq(objValue,2520.57) );

  }
  // Test matt
  if (fopen("../Cbc/matt.mps","r")) {    
    OsiCbcSolverInterface m;
    m.readMps("../Cbc/matt","mps");
    m.setHintParam(OsiDoPresolveInResolve, true, OsiHintDo);
    m.resolve();
    
    std::vector<double *> rays = m.getDualRays(1);
    std::cout << "Dual Ray: " << std::endl;
    for(int i = 0; i < m.getNumRows(); i++){
      if(fabs(rays[0][i]) > 0.00001)
        std::cout << i << " : " << rays[0][i] << std::endl;
    }
    
    std::cout << "isProvenOptimal = " << m.isProvenOptimal() << std::endl;
    std::cout << "isProvenPrimalInfeasible = " << m.isProvenPrimalInfeasible()
         << std::endl;
    
    delete [] rays[0];
    
  }

  // Build a model
  {    
    OsiCbcSolverInterface model;
    std::string fn = mpsDir+"p0033";
    model.readMps(fn.c_str(),"mps");
    // Point to data
    int numberRows = model.getNumRows();
    const double * rowLower = model.getRowLower();
    const double * rowUpper = model.getRowUpper();
    int numberColumns = model.getNumCols();
    const double * columnLower = model.getColLower();
    const double * columnUpper = model.getColUpper();
    const double * columnObjective = model.getObjCoefficients();
    // get row copy
    CoinPackedMatrix rowCopy = *model.getMatrixByRow();
    const int * column = rowCopy.getIndices();
    const int * rowLength = rowCopy.getVectorLengths();
    const CoinBigIndex * rowStart = rowCopy.getVectorStarts();
    const double * element = rowCopy.getElements();
    
    // solve
    model.initialSolve();
    // Now build new model
    CoinModel build;
    // Row bounds
    int iRow;
    for (iRow=0;iRow<numberRows;iRow++) {
      build.setRowBounds(iRow,rowLower[iRow],rowUpper[iRow]);
    }
    // Column bounds and objective
    int iColumn;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      build.setColumnLower(iColumn,columnLower[iColumn]);
      build.setColumnUpper(iColumn,columnUpper[iColumn]);
      build.setObjective(iColumn,columnObjective[iColumn]);
    }
    // Adds elements one by one by row (backwards by row)
    for (iRow=numberRows-1;iRow>=0;iRow--) {
      int start = rowStart[iRow];
      for (int j=start;j<start+rowLength[iRow];j++) 
        build(iRow,column[j],element[j]);
    }
    // Now create Model
    OsiCbcSolverInterface model2;
    model2.loadFromCoinModel(build);
    model2.initialSolve();
    // Save - should be continuous
    model2.writeMps("continuous");
    int * whichInteger = new int[numberColumns];
    for (iColumn=0;iColumn<numberColumns;iColumn++) 
      whichInteger[iColumn]=iColumn;
    // mark as integer
    model2.setInteger(whichInteger,numberColumns);
    delete [] whichInteger;
    // save - should be integer
    model2.writeMps("integer");
    
    // Now do with strings attached
    // Save build to show how to go over rows
    CoinModel saveBuild = build;
    build = CoinModel();
    // Column bounds
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      build.setColumnLower(iColumn,columnLower[iColumn]);
      build.setColumnUpper(iColumn,columnUpper[iColumn]);
    }
    // Objective - half the columns as is and half with multiplier of "1.0+multiplier"
    // Pick up from saveBuild (for no reason at all)
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      double value = saveBuild.objective(iColumn);
      if (iColumn*2<numberColumns) {
        build.setObjective(iColumn,columnObjective[iColumn]);
      } else {
        // create as string
        char temp[100];
        sprintf(temp,"%g + abs(%g*multiplier)",value,value);
        build.setObjective(iColumn,temp);
      }
    }
    // It then adds rows one by one but for half the rows sets their values
    //      with multiplier of "1.0+1.5*multiplier"
    for (iRow=0;iRow<numberRows;iRow++) {
      if (iRow*2<numberRows) {
        // add row in simple way
        int start = rowStart[iRow];
        build.addRow(rowLength[iRow],column+start,element+start,
                     rowLower[iRow],rowUpper[iRow]);
      } else {
        // As we have to add one by one let's get from saveBuild
        CoinModelLink triple=saveBuild.firstInRow(iRow);
        while (triple.column()>=0) {
          int iColumn = triple.column();
          if (iColumn*2<numberColumns) {
            // just value as normal
            build(iRow,triple.column(),triple.value());
          } else {
            // create as string
            char temp[100];
            sprintf(temp,"%g + (1.5*%g*multiplier)",triple.value(), triple.value());
            build(iRow,iColumn,temp);
          }
          triple=saveBuild.next(triple);
        }
        // but remember to do rhs
        build.setRowLower(iRow,rowLower[iRow]);
        build.setRowUpper(iRow,rowUpper[iRow]);
      }
    }
    // If small switch on error printing
    if (numberColumns<50)
      build.setLogLevel(1);
    int numberErrors=model2.loadFromCoinModel(build);
    // should fail as we never set multiplier
    assert (numberErrors);
    build.associateElement("multiplier",0.0);
    numberErrors=model2.loadFromCoinModel(build);
    assert (!numberErrors);
    model2.initialSolve();
    // It then loops with multiplier going from 0.0 to 2.0 in increments of 0.1
    for (double multiplier=0.0;multiplier<2.0;multiplier+= 0.1) {
      build.associateElement("multiplier",multiplier);
      numberErrors=model2.loadFromCoinModel(build,true);
      assert (!numberErrors);
      model2.resolve();
    }
  }
  // branch and bound
  {    
    OsiCbcSolverInterface m;
    std::string fn = mpsDir+"p0033";
    m.readMps(fn.c_str(),"mps");
    m.initialSolve();
    //m.messageHandler()->setLogLevel(0);
    m.getModelPtr()->messageHandler()->setLogLevel(0);
    m.branchAndBound();
  }
  // branch and bound using CbcModel!!!!!!!
  {    
    OsiCbcSolverInterface mm;
    OsiCbcSolverInterface m(&mm);
    std::string fn = mpsDir+"p0033";
    m.readMps(fn.c_str(),"mps");
    m.initialSolve();
    m.branchAndBound();
  }
#ifdef COIN_HAS_OSL
  // branch and bound using OSL
  {    
    OsiOslSolverInterface mmm;
    OsiCbcSolverInterface mm(&mmm);
    CbcStrategyNull strategy;
    OsiCbcSolverInterface m(&mm,&strategy);
    std::string fn = mpsDir+"p0033";
    m.readMps(fn.c_str(),"mps");
    //m.initialSolve();
    m.branchAndBound();
  }
#endif
  int errCnt = 0;
  // Do common solverInterface testing 
  {
    OsiCbcSolverInterface m;
    errCnt += OsiSolverInterfaceCommonUnitTest(&m, mpsDir,netlibDir);
  }
  {
    OsiCbcSolverInterface mm;
    OsiCbcSolverInterface m(&mm);
    errCnt += OsiSolverInterfaceCommonUnitTest(&m, mpsDir,netlibDir);
  }
#ifdef COIN_HAS_OSL
  {
    OsiOslSolverInterface mm;
    OsiCbcSolverInterface m(&mm);
    errCnt += OsiSolverInterfaceCommonUnitTest(&m, mpsDir,netlibDir);
  }
#endif
  return errCnt;
}
