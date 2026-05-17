// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <cfloat>
#include <iostream>

#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"
#include "CglLiftAndProject.hpp"
#include "CoinPackedVector.hpp"
#include "CoinSort.hpp"
#include "CoinPackedMatrix.hpp"

//-----------------------------------------------------------------------------
// Generate Lift-and-Project cuts
//------------------------------------------------------------------- 
void CglLiftAndProject::generateCuts(const OsiSolverInterface& si, OsiCuts& cs,
				     const CglTreeInfo /*info*/)
{
  // Assumes the mixed 0-1 problem 
  //
  //   min {cx: <Atilde,x> >= btilde} 
  //
  // is in canonical form with all bounds,
  // including x_t>=0, -x_t>=-1 for x_t binary,
  // explicitly stated in the constraint matrix. 
  // See ~/COIN/Examples/Cgl2/cgl2.cpp 
  // for a general purpose "convert" function. 

  // Reference [BCC]: Balas, Ceria, and Corneujols,
  // "A lift-and-project cutting plane algorithm
  // for mixed 0-1 program", Math Prog 58, (1993) 
  // 295-324.

  // This implementation uses Normalization 1.

  // Given canonical problem and
  // the lp-relaxation solution, x,
  // the LAP cut generator attempts to construct
  // a cut for every x_j such that 0<x_j<1
  // [BCC:307]
 

  // x_j is the strictly fractional binary variable
  // the cut is generated from
  int j = 0; 

  // Get basic problem information
  // let Atilde be an m by n matrix
  const int m = si.getNumRows(); 
  const int n = si.getNumCols(); 
  const double * x = si.getColSolution();

  // Remember - Atildes may have gaps..
  const CoinPackedMatrix * Atilde = si.getMatrixByRow();
  const double * AtildeElements =  Atilde->getElements();
  const int * AtildeIndices =  Atilde->getIndices();
  const CoinBigIndex * AtildeStarts = Atilde->getVectorStarts();
  const int * AtildeLengths = Atilde->getVectorLengths();  
  const CoinBigIndex AtildeFullSize = AtildeStarts[m];
  const double * btilde = si.getRowLower();

  // Set up memory for system (10) [BCC:307]
  // (the problem over the norm intersected 
  //  with the polar cone)
  // 
  // min <<x^T,Atilde^T>,u> + x_ju_0
  // s.t.
  //     <B,w> = (0,...,0,beta_,beta)^T
  //        w  is nonneg for all but the
  //           last two entries, which are free.
  // where 
  // w = (u,v,v_0,u_0)in BCC notation 
  //      u and v are m-vectors; u,v >=0
  //      v_0 and u_0 are free-scalars, and
  //  
  // B = Atilde^T  -Atilde^T  -e_j e_j
  //     btilde^T   e_0^T      0   0
  //     e_0^T      btilde^T   1   0

  // ^T indicates Transpose
  // e_0 is a (AtildeNCols x 1) vector of all zeros 
  // e_j is e_0 with a 1 in the jth position

  // Storing B in column order. B is a (n+2 x 2m+2) matrix 
  // But need to allow for possible gaps in Atilde.
  // At each iteration, only need to change 2 cols and objfunc
  // Sane design of OsiSolverInterface does not permit mucking
  // with matrix.
  // Because we must delete and add cols to alter matrix,
  // and we can only add columns on the end of the matrix
  // put the v_0 and u_0 columns on the end.
  // rather than as described in [BCC]
 
  // Initially allocating B with space for v_0 and u_O cols
  // but not populating, for efficiency.

  // B without u_0 and v_0 is a (n+2 x 2m) size matrix.

  int twoM = 2*m;
  int BNumRows = n+2;
  int BNumCols = twoM+2;
  CoinBigIndex BFullSize = 2*AtildeFullSize+twoM+3;
  double * BElements = new double[BFullSize];
  int * BIndices = new int[BFullSize];
  CoinBigIndex * BStarts = new CoinBigIndex [BNumCols+1];
  int * BLengths = new int[BNumCols];


  int i, k=0;
  CoinBigIndex ij;
  int nPlus1=n+1;
  CoinBigIndex offset = AtildeStarts[m]+m;
  for (i=0; i<m; i++){
    for (ij=AtildeStarts[i];ij<AtildeStarts[i]+AtildeLengths[i];ij++){
      BElements[k]=AtildeElements[ij];
      BElements[k+offset]=-AtildeElements[ij];
      BIndices[k]= AtildeIndices[ij];
      BIndices[k+offset]= AtildeIndices[ij];

      k++;
    }
    BElements[k]=btilde[i];
    BElements[k+offset]=btilde[i];
    BIndices[k]=n;
    BIndices[k+offset]=nPlus1;
    BStarts[i]= AtildeStarts[i]+i;
    BStarts[i+m]=offset+BStarts[i];// = AtildeStarts[m]+m+AtildeStarts[i]+i
    BLengths[i]= AtildeLengths[i]+1;
    BLengths[i+m]= AtildeLengths[i]+1;
    k++;
  }

  BStarts[twoM]=BStarts[twoM-1]+BLengths[twoM-1];

  // Cols that will be deleted each iteration
  int BNumColsLessOne=BNumCols-1;
  int BNumColsLessTwo=BNumCols-2;
  const int delCols[2] = {BNumColsLessOne, BNumColsLessTwo};

  // Set lower bound on u and v
  // u_0, v_0 will be reset as free
  const double solverINFINITY = si.getInfinity();
  double * BColLowers = new double[BNumCols];
  double * BColUppers = new double[BNumCols];
  CoinFillN(BColLowers,BNumCols,0.0);  
  CoinFillN(BColUppers,BNumCols,solverINFINITY); 

  // Set row lowers and uppers.
  // The rhs is zero, for but the last two rows.
  // For these the rhs is beta_
  double * BRowLowers = new double[BNumRows];
  double * BRowUppers = new double[BNumRows];
  CoinFillN(BRowLowers,BNumRows,0.0);  
  CoinFillN(BRowUppers,BNumRows,0.0);
  BRowLowers[BNumRows-2]=beta_;
  BRowUppers[BNumRows-2]=beta_;
  BRowLowers[BNumRows-1]=beta_;
  BRowUppers[BNumRows-1]=beta_;


  // Calculate base objective <<x^T,Atilde^T>,u>
  // Note: at each iteration coefficient u_0
  //       changes to <x^T,e_j>
  //       w=(u,v,beta,v_0,u_0) size 2m+3
  //       So, BOjective[2m+2]=x[j]
  double * BObjective= new double[BNumCols];
  double * Atildex = new double[m];
  CoinFillN(BObjective,BNumCols,0.0);
  Atilde->times(x,Atildex); // Atildex is size m, x is size n
  CoinDisjointCopyN(Atildex,m,BObjective); 

  // Number of cols and size of Elements vector
  // in B without the v_0 and u_0 cols
  CoinBigIndex BFullSizeLessThree = BFullSize-3;

  // Load B matrix into a column orders CoinPackedMatrix
  CoinPackedMatrix * BMatrix = new CoinPackedMatrix(true, BNumRows,
						  BNumColsLessTwo, 
						  BFullSizeLessThree,
						  BElements,BIndices, 
						  BStarts,BLengths);
  // Assign problem into a solver interface 
  // Note: coneSi will cleanup the memory itself
  OsiSolverInterface * coneSi = si.clone(false);
  coneSi->assignProblem (BMatrix, BColLowers, BColUppers, 
		      BObjective,
		      BRowLowers, BRowUppers);

  // Problem sense should default to "min" by default, 
  // but just to be virtuous...
  coneSi->setObjSense(1.0);

  // The plot outline from here on down:
  // coneSi has been assigned B without the u_0 and v_0 columns
  // Calculate base objective <<x^T,Atilde^T>,u>
  // bool haveWarmStart = false;
  // For (j=0; j<n, j++)
  //   if (!isBinary(x_j) || x_j<=0 || x_j>=1) continue;
  //   // IMPROVEME: if(haveWarmStart) check if j attractive
  //   add {-e_j,0,-1} matrix column for v_0
  //   add {e_j,0,0} matrix column for u_0
  //   objective coefficient for u_0 is  x_j 
  //   if (haveWarmStart) 
  //      set warmstart info
  //   solve min{objw:Bw=0; w>=0,except v_0, u_0 free}
  //   if (bounded)
  //      get warmstart info
  //      haveWarmStart=true;
  //      ustar = optimal u solution
  //      ustar_0 = optimal u_0 solution
  //      alpha^T= <ustar^T,Atilde> -ustar_0e_j^T
  //      (double check <alpha^T,x> >= beta_ should be violated)
  //      add <alpha^T,x> >= beta_ to cutset 
  //   endif
  //   delete column for u_0 // this deletes all column info.
  //   delete column for v_0
  // endFor
  // clean up memory
  // return 0;

  int * nVectorIndices = new int[n];
  CoinIotaN(nVectorIndices, n, 0);

  bool haveWarmStart = false;
  bool equalObj1, equalObj2;
  CoinRelFltEq eq;

  double v_0Elements[2] = {-1,1};
  double u_0Elements[1] = {1};

  CoinWarmStart * warmStart = 0;

  double * ustar = new double[m];
  CoinFillN(ustar, m, 0.0);

  double* alpha = new double[n];
  CoinFillN(alpha, n, 0.0);
  const char * intVar = si.getColType();

  for (j=0;j<n;j++){
    if (intVar[i]!=1) continue; // Better to ask coneSi? No! 
                                   // coneSi has no binInfo.
    equalObj1=eq(x[j],0);
    equalObj2=eq(x[j],1);
    if (equalObj1 || equalObj2) continue;
    // IMPROVEME: if (haveWarmStart) check if j attractive;

    // AskLL:wanted to declare u_0 and v_0 packedVec outside loop
    // and setIndices, but didn't see a method to do that(?)
    // (Could "insert". Seems inefficient)
    int v_0Indices[2]={j,nPlus1};
    int u_0Indices[1]={j};
    // 
    CoinPackedVector  v_0(2,v_0Indices,v_0Elements,false);
    CoinPackedVector  u_0(1,u_0Indices,u_0Elements,false);

#if CGL_DEBUG
    const CoinPackedMatrix *see1 = coneSi->getMatrixByRow();
#endif

    coneSi->addCol(v_0,-solverINFINITY,solverINFINITY,0);
    coneSi->addCol(u_0,-solverINFINITY,solverINFINITY,x[j]);
    if(haveWarmStart) {
      coneSi->setWarmStart(warmStart);
      coneSi->resolve();
    }
    else {

#if CGL_DEBUG
      const CoinPackedMatrix *see2 = coneSi->getMatrixByRow();
#endif

      coneSi->initialSolve();
    }
    if(coneSi->isProvenOptimal()){
      warmStart = coneSi->getWarmStart();
      haveWarmStart=true;
      const double * wstar = coneSi->getColSolution();
      CoinDisjointCopyN(wstar, m, ustar);
      Atilde->transposeTimes(ustar,alpha);
      alpha[j]+=wstar[BNumCols-1]; 
      
#if debug
      int p;
      double sum;
      for(p=0;p<n;p++)sum+=alpha[p]*x[p];
      if (sum<=beta_){
	throw CoinError("Cut not violated",
			"cutGeneration",
			"CglLiftAndProject");
      }
#endif

      // add <alpha^T,x> >= beta_ to cutset
      OsiRowCut rc;
      rc.setRow(n,nVectorIndices,alpha);
      rc.setLb(beta_);
      rc.setUb(solverINFINITY);
      cs.insertIfNotDuplicate(rc);
    }
    // delete col for u_o and v_0
    coneSi->deleteCols(2,delCols);

    // clean up memory
  }
  // clean up
  delete [] alpha;
  delete [] ustar;
  delete [] nVectorIndices;
  // BMatrix, BColLowers,BColUppers, BObjective, BRowLowers, BRowUppers
  // are all freed by OsiSolverInterface destructor (?)
  delete [] BLengths;
  delete [] BStarts;
  delete [] BIndices;
  delete [] BElements;
}

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglLiftAndProject::CglLiftAndProject ()
:
CglCutGenerator(),
beta_(1),
epsilon_(1.0e-08),
onetol_(1-epsilon_)
{
  // nothing to do here
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglLiftAndProject::CglLiftAndProject (const CglLiftAndProject & source) :
   CglCutGenerator(source),
   beta_(source.beta_),
   epsilon_(source.epsilon_),
   onetol_(source.onetol_)
{
  // Nothing to do here
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglLiftAndProject::clone() const
{
  return new CglLiftAndProject(*this);
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CglLiftAndProject::~CglLiftAndProject ()
{
  // Nothing to do here
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CglLiftAndProject &
CglLiftAndProject::operator=(
                                         const CglLiftAndProject& rhs)
{
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    beta_=rhs.beta_;
    epsilon_=rhs.epsilon_;
    onetol_=rhs.onetol_;
  }
  return *this;
}
// Create C++ lines to get to current state
std::string
CglLiftAndProject::generateCpp( FILE * fp) 
{
  CglLiftAndProject other;
  fprintf(fp,"0#include \"CglLiftAndProject.hpp\"\n");
  fprintf(fp,"3  CglLiftAndProject liftAndProject;\n");
  if (beta_!=other.beta_)
    fprintf(fp,"3  liftAndProject.setBeta(%d);\n",static_cast<int> (beta_));
  else
    fprintf(fp,"4  liftAndProject.setBeta(%d);\n",static_cast<int> (beta_));
  fprintf(fp,"3  liftAndProject.setAggressiveness(%d);\n",getAggressiveness());
  if (getAggressiveness()!=other.getAggressiveness())
    fprintf(fp,"3  liftAndProject.setAggressiveness(%d);\n",getAggressiveness());
  else
    fprintf(fp,"4  liftAndProject.setAggressiveness(%d);\n",getAggressiveness());
  return "liftAndProject";
}
