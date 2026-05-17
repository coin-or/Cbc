// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdio>

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

#include "CoinPragma.hpp"
#include "CglKnapsackCover.hpp"
#include "CoinPackedMatrix.hpp"

//--------------------------------------------------------------------------
void
CglKnapsackCoverUnitTest(
  const OsiSolverInterface * baseSiP,
  const std::string mpsDir )
{
  int i;
  CoinRelFltEq eq(0.000001);

  // Test default constructor
  {
    CglKnapsackCover kccGenerator;
  }
  
  // Test copy & assignment
  {
    CglKnapsackCover rhs;
    {
      CglKnapsackCover kccGenerator;
      CglKnapsackCover cgC(kccGenerator);
      rhs=kccGenerator;
    }
  }


  // test exactSolveKnapsack
  {  
    CglKnapsackCover kccg;
    const int n=7;
    double c=50;
    double p[n] = {70,20,39,37,7,5,10};
    double w[n] = {31, 10, 20, 19, 4, 3, 6};
    double z;
    int x[n];
    int exactsol = kccg.exactSolveKnapsack(n, c, p, w, z, x);
    assert(exactsol==1);
    assert (z == 107);
    assert (x[0]==1);
    assert (x[1]==0);
    assert (x[2]==0);
    assert (x[3]==1);
    assert (x[4]==0);
    assert (x[5]==0);
    assert (x[6]==0);
  }

  /*
  // Testcase /u/rlh/osl2/mps/scOneInt.mps
  // Model has 3 continous, 2 binary, and 1 general
  // integer variable.
  {
    OsiSolverInterface  * siP = baseSiP->clone();
    int * complement=NULL;
    double * xstar=NULL;

    siP->readMps("../Mps/scOneInt","mps");
    CglKnapsackCover kccg;
    int nCols=siP->getNumCols();
    
    // Test the siP methods for detecting
    // variable type
    int numCont=0, numBinary=0, numIntNonBinary=0, numInt=0;
    for (int thisCol=0; thisCol<nCols; thisCol++) {
      if ( siP->isContinuous(thisCol) ) numCont++;
      if ( siP->isBinary(thisCol) ) numBinary++;
      if ( siP->isIntegerNonBinary(thisCol) ) numIntNonBinary++;
      if ( siP->isInteger(thisCol) ) numInt++;
    }
    assert(numCont==3);
    assert(numBinary==2);
    assert(numIntNonBinary==1);
    assert(numInt==3);
    
    
    // Test initializeCutGenerator
    siP->initialSolve();
    assert(xstar !=NULL);
    for (i=0; i<nCols; i++){
      assert(complement[i]==0);
    }
    int nRows=siP->getNumRows();
    for (i=0; i<nRows; i++){
    int vectorsize = siP->getMatrixByRow()->vectorSize(i);
    assert(vectorsize==2);
    }
    
    kccg.cleanUpCutGenerator(complement,xstar);
    delete siP;
  }
  */  
  
  // Testcase /u/rlh/osl2/mps/tp3.mps
  // Models has 3 cols, 3 rows
  // Row 0 yields a knapsack, others do not.
  {
    // setup
    OsiSolverInterface  * siP = baseSiP->clone();
    std::string fn(mpsDir+"tp3");
    siP->readMps(fn.c_str(),"mps");     
    // All integer variables should be binary.
    // Assert that this is true.
    for ( i = 0;  i < siP->getNumCols();  i++ )
      if ( siP->isInteger(i) ) 
        assert(siP->getColUpper()[i]==1.0 && siP->isBinary(i));  
    OsiCuts cs;
    CoinPackedVector krow;
    double b=0;
    int nCols=siP->getNumCols();
    int * complement=new int [nCols];
    double * xstar=new double [nCols];

    CglKnapsackCover kccg;

    // solve LP relaxation
    // a "must" before calling initialization
    siP->initialSolve();
    double lpRelaxBefore=siP->getObjValue();
    std::cout<<"Initial LP value: "<<lpRelaxBefore<<std::endl;
    assert( eq(siP->getObjValue(), 97.185) );
    double mycs[] = {.627, .667558333333, .038};
    siP->setColSolution(mycs);
    const double *colsol = siP->getColSolution(); 
    int k;
    for (k=0; k<nCols; k++){
      xstar[k]=colsol[k];
      complement[k]=0;
    }
    
    // test deriveAKnapsack
    int rind = ( siP->getRowSense()[0] == 'N' ) ? 1 : 0;
    const CoinShallowPackedVector reqdBySunCC = siP->getMatrixByRow()->getVector(rind) ;
    int deriveaknap = kccg.deriveAKnapsack(*siP, cs, krow,b,complement,xstar,rind,reqdBySunCC);
    assert(deriveaknap ==1);
    assert(complement[0]==0);
    assert(complement[1]==1);
    assert(complement[2]==1);
    int inx[3] = {0,1,2};
    double el[3] = {161, 120, 68};
    CoinPackedVector r;
    r.setVector(3,inx,el);
    assert (krow == r);
    //assert (b == 183.0); ????? but x1 and x2 at 1 is valid 
    
    // test findGreedyCover 
    CoinPackedVector cover,remainder;
#if 0
    int findgreedy =  kccg.findGreedyCover( 0, krow, b, xstar, cover, remainder );
    assert( findgreedy == 1 );
    int coveri = cover.getNumElements();
    assert( cover.getNumElements() == 2);
    coveri = cover.getIndices()[0];
    assert( cover.getIndices()[0] == 0);
    assert( cover.getIndices()[1] == 1);
    assert( cover.getElements()[0] == 161.0);
    assert( cover.getElements()[1] == 120.0);
    assert( remainder.getNumElements() == 1);
    assert( remainder.getIndices()[0] == 2);
    assert( remainder.getElements()[0] == 68.0);

    // test liftCoverCut
    CoinPackedVector cut;
    double * rowupper = ekk_rowupper(model);
    double cutRhs = cover.getNumElements() - 1.0;
    kccg.liftCoverCut(b, krow.getNumElements(),
      cover, remainder,
      cut);
    assert ( cut.getNumElements() == 3 );
    assert ( cut.getIndices()[0] == 0 );
    assert ( cut.getIndices()[1] == 1 );
    assert ( cut.getIndices()[2] == 2 );
    assert( cut.getElements()[0] == 1 );
    assert( cut.getElements()[1] == 1 );
    assert( eq(cut.getElements()[2], 0.087719) );
    
    // test liftAndUncomplementAndAdd
    OsiCuts cuts;    
    kccg.liftAndUncomplementAndAdd(*siP.getRowUpper()[0],krow,b,complement,0,
      cover,remainder,cuts);   
    int sizerowcuts = cuts.sizeRowCuts();
    assert ( sizerowcuts== 1 );
    OsiRowCut testRowCut = cuts.rowCut(0);
    CoinPackedVector testRowPV = testRowCut.row(); 
    OsiRowCut sampleRowCut;
    const int sampleSize = 3;
    int sampleCols[sampleSize]={0,1,2};
    double sampleElems[sampleSize]={1.0,-1.0,-0.087719};
    sampleRowCut.setRow(sampleSize,sampleCols,sampleElems);
    sampleRowCut.setLb(-DBL_MAX);
    sampleRowCut.setUb(-0.087719);
    bool equiv =  testRowPV.equivalent(sampleRowCut.row(),CoinRelFltEq(1.0e-05) );
    assert ( equiv );
#endif
    
    // test find PseudoJohnAndEllisCover
    cover.setVector(0,NULL, NULL);
    remainder.setVector(0,NULL,NULL);

    rind = ( siP->getRowSense()[0] == 'N' ) ? 1 : 0;
    int findPJE =  kccg.findPseudoJohnAndEllisCover( rind, krow, 
						     b, xstar, cover, remainder );
    assert( findPJE == 1 );
    assert ( cover.getIndices()[0] == 0 );
    assert ( cover.getIndices()[1] == 2 );
    assert ( cover.getElements()[0] == 161 );    
    assert ( cover.getElements()[1] == 68 );    
    assert ( remainder.getIndices()[0] == 1 );
    assert ( remainder.getElements()[0] == 120 );    
    OsiCuts cuts;    
    kccg.liftAndUncomplementAndAdd((*siP).getRowUpper()[rind],krow,b, complement, rind,
      cover,remainder,cuts);   
    assert (cuts.sizeRowCuts() == 1 );

    OsiRowCut testRowCut = cuts.rowCut(0);
    CoinPackedVector testRowPV = testRowCut.row();


    const int sampleSize = 3;
    int sampleCols[sampleSize]={0,1,2};
    double sampleElems[sampleSize]={1.0, -1.0, -1.0};
    OsiRowCut sampleRowCut;
    sampleRowCut.setRow(sampleSize,sampleCols,sampleElems);
    sampleRowCut.setLb(-COIN_DBL_MAX);
    sampleRowCut.setUb(-1.0);
    
    // test for 'close enough'
    assert( testRowPV.isEquivalent(sampleRowCut.row(),CoinRelFltEq(1.0e-05) ) );
    // Reset complement & test next row
    for (i=0; i<nCols; i++){
      complement[i]=0;
    }

    rind++;
    const CoinShallowPackedVector reqdBySunCC2 = siP->getMatrixByRow()->getVector(rind) ;
    deriveaknap = kccg.deriveAKnapsack(*siP,cuts,krow,b,complement,xstar,rind,reqdBySunCC2);
    assert(deriveaknap==0);
    
    // Reset complement & test next row
    for (i=0; i<nCols; i++){
      complement[i]=0;
    }
    const CoinShallowPackedVector reqdBySunCC3 = siP->getMatrixByRow()->getVector(2) ;
    deriveaknap = kccg.deriveAKnapsack(*siP,cuts,krow,b,complement,xstar,2,
				       reqdBySunCC3);
    assert(deriveaknap == 0);
    
    // Clean up
    delete [] complement;
    delete [] xstar;
    
    delete siP;
  }

#if 0
  // Testcase /u/rlh/osl2/mps/tp4.mps
  // Models has 6 cols, 1 knapsack row and 
  // 3 rows explicily bounding variables
  // Row 0 yields a knapsack cover cut 
  // using findGreedyCover which moves the 
  // LP objective function value.
  {
    // Setup
    EKKContext * env=ekk_initializeContext();
    EKKModel * model = ekk_newModel(env,"");
    OsiSolverInterface si(model);
    ekk_importModel(model, "tp4.mps");
    CglKnapsackCover kccg;
    kccg.ekk_validateIntType(si);     
    
    // Solve the LP relaxation of the model and
    // print out ofv for sake of comparison 
    ekk_allSlackBasis(model);
    ekk_crash(model,1); 
    ekk_primalSimplex(model,1);
    double lpRelaxBefore=ekk_getRobjvalue(model);
#ifdef CGL_DEBUG
    printf("\n\nOrig LP min=%f\n",lpRelaxBefore);
#endif
    
    // Determine if lp sol is ip optimal
    // Note: no ekk_function to do this
    int nCols=ekk_getInumcols(model);
    double * optLpSol = ekk_colsol(model);
    int ipOpt = 1;
    i=0;
    while (i++<nCols && ipOpt){
      if(optLpSol[i] < 1.0-1.0e-08 && optLpSol[i]> 1.0e-08) ipOpt = 0;
    }
    
    if (ipOpt){
#ifdef CGL_DEBUG
      printf("Lp solution is within ip optimality tolerance\n");
#endif
    }    
    else {
      OsiSolverInterface iModel(model);
      OsiCuts cuts;    
      
      // Test generateCuts method
      kccg.generateCuts(iModel,cuts);
      OsiSolverInterface::ApplyCutsReturnCode rc = iModel.applyCuts(cuts);
      
      ekk_mergeBlocks(model,1);         
      ekk_dualSimplex(model);
      double lpRelaxAfter=ekk_getRobjvalue(model); 
#ifdef CGL_DEBUG
      printf("\n\nFinal LP min=%f\n",lpRelaxAfter);
#endif
      assert( lpRelaxBefore < lpRelaxAfter );
      
      // This may need to be updated as other 
      // minimal cover finders are added
      assert( cuts.sizeRowCuts() == 1 );
      OsiRowCut testRowCut = cuts.rowCut(0);
      CoinPackedVector testRowPV = testRowCut.row();
      
      OsiRowCut sampleRowCut;
      const int sampleSize = 6;
      int sampleCols[sampleSize]={0,1,2,3,4,5};
      double sampleElems[sampleSize]={1.0,1.0,1.0,1.0,0.5, 2.0};
      sampleRowCut.setRow(sampleSize,sampleCols,sampleElems);
      sampleRowCut.setLb(-DBL_MAX);
      sampleRowCut.setUb(3.0);
      bool equiv = testRowPV.equivalent(sampleRowCut.row(),CoinRelFltEq(1.0e-05) );
      assert( testRowPV.equivalent(sampleRowCut.row(),CoinRelFltEq(1.0e-05) ) );
    }
    
    // Exit out of OSL
    ekk_deleteModel(model);
    ekk_endContext(env);
    
  }
#endif


  // Testcase /u/rlh/osl2/mps/tp5.mps
  // Models has 6 cols, 1 knapsack row and 
  // 3 rows explicily bounding variables
  // Row 0 yields a knapsack cover cut 
  // using findGreedyCover which moves the 
  // LP objective function value.
  {
    // Setup
    OsiSolverInterface  * siP = baseSiP->clone();
    std::string fn(mpsDir+"tp5");
    siP->readMps(fn.c_str(),"mps");
    // All integer variables should be binary.
    // Assert that this is true.
    for ( i = 0;  i < siP->getNumCols();  i++ )
      if ( siP->isInteger(i) ) 
        assert(siP->getColUpper()[i]==1.0 && siP->isBinary(i));  
    CglKnapsackCover kccg;
    
    // Solve the LP relaxation of the model and
    // print out ofv for sake of comparison 
    siP->initialSolve();
    double lpRelaxBefore=siP->getObjValue();
    assert( eq(lpRelaxBefore, -51.66666666667) );
    double mycs[] = {.8999999999, .899999999999, .89999999999, 1.110223e-16, .5166666666667, 0};
    siP->setColSolution(mycs);
#ifdef CGL_DEBUG
    printf("\n\nOrig LP min=%f\n",lpRelaxBefore);
#endif
    
    // Determine if lp sol is 0/1 optimal
    int nCols=siP->getNumCols();
    const double * optLpSol = siP->getColSolution();
    bool ipOpt = true;
    i=0;
    while (i++<nCols && ipOpt){
      if(optLpSol[i] > kccg.epsilon_ && optLpSol[i] < kccg.onetol_) ipOpt = false;
    }
    
    if (ipOpt){
#ifdef CGL_DEBUG
      printf("Lp solution is within ip optimality tolerance\n");
#endif
    }    
    else {
      // set up
      OsiCuts cuts;    
      CoinPackedVector krow;
      double b=0.0;
      int * complement=new int[nCols];
      double * xstar=new double[nCols];
      // initialize cut generator
      const double *colsol = siP->getColSolution(); 
      for (i=0; i<nCols; i++){
	xstar[i]=colsol[i];
	complement[i]=0;
      }
      int row = ( siP->getRowSense()[0] == 'N' ) ? 1 : 0;
      // transform row into canonical knapsack form
      const CoinShallowPackedVector reqdBySunCC = siP->getMatrixByRow()->getVector(row) ;
      if (kccg.deriveAKnapsack(*siP, cuts, krow, b, complement, xstar, row,reqdBySunCC)){
        CoinPackedVector cover, remainder;  
        // apply greedy logic to detect violated minimal cover inequalities
        if (kccg.findGreedyCover(row, krow, b, xstar, cover, remainder) == 1){
          // lift, uncomplements, and add cut to cut set
          kccg.liftAndUncomplementAndAdd((*siP).getRowUpper()[row],krow, b, complement, row, cover, remainder, cuts);   
        }  
        // reset optimal column solution (xstar) information in OSL     
        const double * rowupper = siP->getRowUpper();
	int k;
        if (fabs(b-rowupper[row]) > 1.0e-05) {
          for(k=0; k<krow.getNumElements(); k++) {
            if (complement[krow.getIndices()[k]]){
              xstar[krow.getIndices()[k]]= 1.0-xstar[krow.getIndices()[k]];
              complement[krow.getIndices()[k]]=0;
            }
          }
        }  
        // clean up
        delete [] complement;
	delete [] xstar;
      }
      // apply the cuts
      OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cuts);
      
      siP->resolve();
      double lpRelaxAfter=siP->getObjValue();
      assert( eq(lpRelaxAfter, -30.0) );
#ifdef CGL_DEBUG
      printf("\n\nFinal LP min=%f\n",lpRelaxAfter);
#endif
      // test that expected cut was detected
      assert( lpRelaxBefore < lpRelaxAfter );
      assert( cuts.sizeRowCuts() == 1 );
      OsiRowCut testRowCut = cuts.rowCut(0);
      CoinPackedVector testRowPV = testRowCut.row();
      OsiRowCut sampleRowCut;
      const int sampleSize = 6;
      int sampleCols[sampleSize]={0,1,2,3,4,5};
      double sampleElems[sampleSize]={1.0,1.0,1.0,0.25,1.0,2.0};
      sampleRowCut.setRow(sampleSize,sampleCols,sampleElems);
      sampleRowCut.setLb(-COIN_DBL_MAX);
      sampleRowCut.setUb(3.0);
      assert(testRowPV.isEquivalent(sampleRowCut.row(),CoinRelFltEq(1.0e-05)));
    }
    
    delete siP;
  }
 

  // Testcase /u/rlh/osl2/mps/p0033
  // Miplib3 problem p0033
  // Test that no cuts chop off the optimal solution
  {
    // Setup
    OsiSolverInterface  * siP = baseSiP->clone();
    std::string fn(mpsDir+"p0033");
    siP->readMps(fn.c_str(),"mps");
    // All integer variables should be binary.
    // Assert that this is true.
    for ( i = 0;  i < siP->getNumCols();  i++ )
      if ( siP->isInteger(i) ) 
        assert(siP->getColUpper()[i]==1.0 && siP->isBinary(i));  
    int nCols=siP->getNumCols();
    CglKnapsackCover kccg;

    // Solve the LP relaxation of the model and
    // print out ofv for sake of comparison 
    siP->initialSolve();
    double lpRelaxBefore=siP->getObjValue();
    assert( eq(lpRelaxBefore, 2520.5717391304347) );
    double mycs[] = {0, 1, 0, 0, -2.0837010502455788e-19, 1, 0, 0, 1,
		       0.021739130434782594, 0.35652173913043478, 
		       -6.7220534694101275e-18, 5.3125906451789717e-18, 
		       1, 0, 1.9298798670241979e-17, 0, 0, 0,
		       7.8875708048320448e-18, 0.5, 0, 
		       0.85999999999999999, 1, 1, 0.57999999999999996,
		       1, 0, 1, 0, 0.25, 0, 0.67500000000000004};
    siP->setColSolution(mycs);
#ifdef CGL_DEBUG
    printf("\n\nOrig LP min=%f\n",lpRelaxBefore);
#endif
    
    OsiCuts cuts;    
    
    // Test generateCuts method
    kccg.generateCuts(*siP,cuts);
    OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cuts);
    
    siP->resolve();
    double lpRelaxAfter=siP->getObjValue(); 
    assert( eq(lpRelaxAfter, 2829.0597826086955) );
#ifdef CGL_DEBUG
    printf("\n\nOrig LP min=%f\n",lpRelaxBefore);
    printf("\n\nFinal LP min=%f\n",lpRelaxAfter);
#endif
    assert( lpRelaxBefore < lpRelaxAfter );
    
    // the CoinPackedVector p0033 is the optimal
    // IP solution to the miplib problem p0033
    int objIndices[14] = { 
       0,  6,  7,  9, 13, 17, 18,
      22, 24, 25, 26, 27, 28, 29 };
    CoinPackedVector p0033(14,objIndices,1.0);

    // Sanity check
    const double *  objective=siP->getObjCoefficients();
    double ofv =0 ;
    int r;
    for (r=0; r<nCols; r++){
      ofv=ofv + p0033[r]*objective[r];
    }
    CoinRelFltEq eq;
    assert( eq(ofv,3089.0) );

    int nRowCuts = cuts.sizeRowCuts();
    OsiRowCut rcut;
    CoinPackedVector rpv;
    for (i=0; i<nRowCuts; i++){
      rcut = cuts.rowCut(i);
      rpv = rcut.row();
      double p0033Sum = (rpv*p0033).sum();
      assert (p0033Sum <= rcut.ub() );
    }
  
    delete siP;
  } 

  // if a debug file is there then look at it
  {
    FILE * fp = fopen("knapsack.debug","r");
    if (fp) {
      int ncol,nel;
      double up;
      int x = fscanf(fp,"%d %d %lg",&ncol,&nel,&up);
      if (x<=0)
	throw("bad fscanf");
      printf("%d columns, %d elements, upper %g\n",ncol,nel,up);
      double * sol1 = new double[nel];
      double * el1 = new double[nel];
      int * col1 = new int[nel];
      CoinBigIndex * start = new CoinBigIndex [ncol+1];
      memset(start,0,ncol*sizeof(CoinBigIndex ));
      int * row = new int[nel];
      int i;
      for (i=0;i<nel;i++) {
	x=fscanf(fp,"%d %lg %lg",col1+i,el1+i,sol1+i);
	if (x<=0)
	  throw("bad fscanf");
	printf("[%d, e=%g, v=%g] ",col1[i],el1[i],sol1[i]);
	start[col1[i]]=1;
	row[i]=0;
      }
      printf("\n");
      // Setup
      OsiSolverInterface  * siP = baseSiP->clone();
      
      double lo=-1.0e30;
      double * upper = new double[ncol];
      start[ncol]=nel;
      CoinBigIndex last=0;
      for (i=0;i<ncol;i++) {
	upper[i]=1.0;
	CoinBigIndex marked=start[i];
	start[i]=last;
	if (marked)
	  last++;
      }
      siP->loadProblem(ncol,1,start,row,el1,NULL,upper,NULL,&lo,&up);
      // use upper for solution
      memset(upper,0,ncol*sizeof(double));
      for (i=0;i<nel;i++) {
	int icol=col1[i];
	upper[icol]=sol1[i];
	siP->setInteger(icol);
      }
      siP->setColSolution(upper);
      delete [] sol1;
      delete [] el1;
      delete [] col1;
      delete [] start;
      delete [] row;
      delete [] upper;
      CglKnapsackCover kccg;
      
      OsiCuts cuts;    
      
      // Test generateCuts method
      kccg.generateCuts(*siP,cuts);
      // print out and compare to known cuts
      int numberCuts = cuts.sizeRowCuts();
      if (numberCuts) {
	for (i=0;i<numberCuts;i++) {
	  OsiRowCut * thisCut = cuts.rowCutPtr(i);
	  int n=thisCut->row().getNumElements();
	  printf("Cut %d has %d entries, rhs %g %g =>",i,n,thisCut->lb(),
		 thisCut->ub());
	  int j;
	  const int * index = thisCut->row().getIndices();
	  const double * element = thisCut->row().getElements();
	  for (j=0;j<n;j++) {
	    printf(" (%d,%g)",index[j],element[j]);
	  }
	  printf("\n");
	}
      }
      fclose(fp);
    }
  }

  // Testcase /u/rlh/osl2/mps/p0201
  // Miplib3 problem p0282
  // Test that no cuts chop off the optimal ip solution
  {
    // Setup
    OsiSolverInterface  * siP = baseSiP->clone();
    std::string fn(mpsDir+"p0201");
    siP->readMps(fn.c_str(),"mps");
    // All integer variables should be binary.
    // Assert that this is true.
    for ( i = 0;  i < siP->getNumCols();  i++ )
      if ( siP->isInteger(i) ) 
        assert(siP->getColUpper()[i]==1.0 && siP->isBinary(i));    

    const int nCols=siP->getNumCols();
    CglKnapsackCover kccg;
    
    // Solve the LP relaxation of the model and
    // print out ofv for sake of comparisn 
    siP->initialSolve();
    double lpRelaxBefore=siP->getObjValue();
    assert( eq(lpRelaxBefore, 6875.) );
    double mycs[] =
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
       0, 0.5, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.5, 
       0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.5, 0, 0, 
       0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.5, 0, 0, 0, 0.5, 
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.5, 0, 0, 0, 0.5, 0, 0, 
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.5, 0, 0, 0, 0.5, 0, 0, 0, 0, 
       0, 0, 0, 0, 0, 0, 0, 0, 1, 0.5, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 
       0, 0, 0, 0, 0, 0, 1, 0.5, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 
       0, 0, 0, 0, 1, 0.5, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
       0, 0, 1, 0.5, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
       1};
    siP->setColSolution(mycs);
#ifdef CGL_DEBUG
    printf("\n\nOrig LP min=%f\n",lpRelaxBefore);
#endif
    
    OsiCuts cuts;    
    
    // Test generateCuts method
    kccg.generateCuts(*siP,cuts);
    OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cuts);
    
    siP->resolve();
    double lpRelaxAfter=siP->getObjValue(); 
    assert( eq(lpRelaxAfter, 7125) );
#ifdef CGL_DEBUG
    printf("\n\nOrig LP min=%f\n",lpRelaxBefore);
    printf("\n\nFinal LP min=%f\n",lpRelaxAfter);
#endif
    assert( lpRelaxBefore < lpRelaxAfter );
 
    // Optimal IP solution to p0201    
    int objIndices[22] = { 8, 10,  21,  38,  39,  56,
      60,   74, 79,  92, 94, 110, 111, 128, 132, 146, 
      151,164, 166, 182,183, 200 };
    CoinPackedVector p0201(22,objIndices,1.0);
    
    // Sanity check
    const double *  objective=siP->getObjCoefficients();
    double ofv =0 ;
    int r;
    for (r=0; r<nCols; r++){
      ofv=ofv + p0201[r]*objective[r];
    }
    CoinRelFltEq eq;
    assert( eq(ofv,7615.0) );
    //printf("p0201 optimal ofv = %g\n",ofv); 

    int nRowCuts = cuts.sizeRowCuts();
    OsiRowCut rcut;
    CoinPackedVector rpv;
    for (i=0; i<nRowCuts; i++){
      rcut = cuts.rowCut(i);
      rpv = rcut.row();
      double p0201Sum = (rpv*p0201).sum();
      assert (p0201Sum <= rcut.ub() );
    }
  
    delete siP;
  } 

 
  // see if I get the same covers that N&W get
  {
    OsiSolverInterface * siP=baseSiP->clone();
    std::string fn(mpsDir+"nw460");
    siP->readMps(fn.c_str(),"mps");   
    // All integer variables should be binary.
    // Assert that this is true.
    for ( i = 0;  i < siP->getNumCols();  i++ )
      if ( siP->isInteger(i) ) 
        assert(siP->getColUpper()[i]==1.0 && siP->isBinary(i));  
    CglKnapsackCover kccg;
    
    // Solve the LP relaxation of the model and
    // print out ofv for sake of comparison 
    siP->initialSolve();
    double lpRelaxBefore=siP->getObjValue();
    assert( eq(lpRelaxBefore, -225.68951787852194) );
    double mycs[] = {0.7099213482046447, 0, 0.34185802225477174, 1, 1, 0, 1, 1, 0};
    siP->setColSolution(mycs);

    OsiCuts cuts;    
    
    // Test generateCuts method
    kccg.generateCuts(*siP,cuts);
    OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cuts);
    
    siP->resolve();
    double lpRelaxAfter=siP->getObjValue(); 
    assert( eq(lpRelaxAfter, -176) );
#ifdef CGL_DEBUG
    printf("\n\nOrig LP min=%f\n",lpRelaxBefore);
    printf("\n\nFinal LP min=%f\n",lpRelaxAfter);
#endif
#ifdef MJS
    assert( lpRelaxBefore < lpRelaxAfter );
#endif    
    
    int nRowCuts = cuts.sizeRowCuts();
    OsiRowCut rcut;
    CoinPackedVector rpv;
    for (i=0; i<nRowCuts; i++){
      rcut = cuts.rowCut(i);
      rpv = rcut.row();
      int j;
      printf("Row cut number %i has rhs = %g\n",i,rcut.ub());
      for (j=0; j<rpv.getNumElements(); j++){
        printf("index %i, element %g\n", rpv.getIndices()[j], rpv.getElements()[j]);
      }
      printf("\n");
    }
    delete siP; 
  }

  // Debugging: try "exmip1.mps"
  {
    // Setup
    OsiSolverInterface  * siP = baseSiP->clone();
    std::string fn(mpsDir+"exmip1");
    siP->readMps(fn.c_str(),"mps");   
    // All integer variables should be binary.
    // Assert that this is true.
    for ( i = 0;  i < siP->getNumCols();  i++ )
      if ( siP->isInteger(i) ) 
        assert(siP->getColUpper()[i]==1.0 && siP->isBinary(i));  
    CglKnapsackCover kccg;
    
    // Solve the LP relaxation of the model and
    // print out ofv for sake of comparison 
    siP->initialSolve();
    double lpRelaxBefore=siP->getObjValue();
    assert( eq(lpRelaxBefore, 3.2368421052631575) );
    double mycs[] = {2.5, 0, 0, 0.6428571428571429, 0.5, 4, 0, 0.26315789473684253};
    siP->setColSolution(mycs);
    // Test generateCuts method
    OsiCuts cuts;    
    kccg.generateCuts(*siP,cuts);
    OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cuts);
    
    siP->resolve();
    double lpRelaxAfter=siP->getObjValue();
    assert( eq(lpRelaxAfter, 3.2368421052631575) );
#ifdef CGL_DEBUG
    printf("\n\nOrig LP min=%f\n",lpRelaxBefore);
    printf("\n\nFinal LP min=%f\n",lpRelaxAfter);
#endif
    assert( lpRelaxBefore <= lpRelaxAfter );

    delete siP;
  } 

#ifdef CGL_DEBUG
  // See what findLPMostViolatedMinCover for knapsack with 2 elements does
  {
    int nCols = 2;
    int row = 1;
    CoinPackedVector krow;
    double e[2] = {5,10};
    int ii[2] = {0,1};
    krow.setVector(nCols,ii,e);
    double b=11;
    double xstar[2] = {.2,.9};
    CoinPackedVector cover;
    CoinPackedVector remainder;
    CglKnapsackCover kccg;
    kccg.findLPMostViolatedMinCover(nCols, row, krow, b, xstar, cover, remainder);
    printf("num in cover = %i\n",cover.getNumElements());
    int j;
    for (j=0; j<cover.getNumElements(); j++){
      printf(" index %i element % g\n", cover.getIndices()[j], cover.getElements()[j]);
    }
  }
#endif 

#ifdef CGL_DEBUG
  // see what findLPMostViolatedMinCover does
  {
    int nCols = 5;
    int row = 1;
    CoinPackedVector krow;
    double e[5] = {1,1,1,1,10};
    int ii[5] = {0,1,2,3,4};
    krow.setVector(nCols,ii,e);
    double b=11;
    double xstar[5] = {.9,.9,1,1,.1};
    CoinPackedVector cover;
    CoinPackedVector remainder;
    CglKnapsackCover kccg;
    kccg.findLPMostViolatedMinCover(nCols, row, krow, b, xstar, cover, remainder);
    printf("num in cover = %i\n",cover.getNumElements());
    int j;
    for (j=0; j<cover.getNumElements(); j++){
      printf(" index %i element % g\n", cover.getIndices()[j], cover.getElements()[j]);
    }
  }
#endif

}

