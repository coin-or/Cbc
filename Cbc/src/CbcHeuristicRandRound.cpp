// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <vector.h>

#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcHeuristicRandRound.hpp"
#include "OsiClpSolverInterface.hpp"
#include  "CoinTime.hpp"

// Default Constructor
CbcHeuristicRandRound::CbcHeuristicRandRound() 
  :CbcHeuristic()
{
}

// Constructor with model - assumed before cuts

CbcHeuristicRandRound::CbcHeuristicRandRound(CbcModel & model)
  :CbcHeuristic(model)
{
}

// Destructor 
CbcHeuristicRandRound::~CbcHeuristicRandRound ()
{
}

// Clone
CbcHeuristic *
CbcHeuristicRandRound::clone() const
{
  return new CbcHeuristicRandRound(*this);
}
// Create C++ lines to get to current state
void 
CbcHeuristicRandRound::generateCpp( FILE * fp) 
{
  CbcHeuristicRandRound other;
  fprintf(fp,"0#include \"CbcHeuristicRandRound.hpp\"\n");
  fprintf(fp,"3  CbcHeuristicRandRound heuristicPFX(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp,"heuristicPFX");
  fprintf(fp,"3  cbcModel->addHeuristic(&heuristicPFX);\n");
}

// Copy constructor 
CbcHeuristicRandRound::CbcHeuristicRandRound(const CbcHeuristicRandRound & rhs)
:
  CbcHeuristic(rhs)
{
}

// Assignment operator 
CbcHeuristicRandRound & 
CbcHeuristicRandRound::operator=( const CbcHeuristicRandRound& rhs)
{
  if (this!=&rhs) {
    CbcHeuristic::operator=(rhs);
  }
  return *this;
}

// Resets stuff if model changes
void 
CbcHeuristicRandRound::resetModel(CbcModel * model)
{
  //CbcHeuristic::resetModel(model);
}
/*
  Randomized Rounding Heuristic
  Returns 1 if solution, 0 if not */
int
CbcHeuristicRandRound::solution(double & solutionValue,
			 double * betterSolution)
{
  std::cout << "Lucky you! You're in the Randomized Rounding Heuristic" << std::endl;
  // The struct should be moved to member data

  typedef struct {
    int numberSolutions;
    int maximumSolutions;
    int numberColumns;
    double ** solution;
    int * numberUnsatisfied;
  } clpSolution;

  double start = CoinCpuTime();
  numCouldRun_++; // Todo: Ask JJHF what this for. 

  OsiClpSolverInterface * clpSolver 
    = dynamic_cast<OsiClpSolverInterface *> (model_->solver());
  assert (clpSolver);
  ClpSimplex * simplex = clpSolver->getModelPtr();

  // Initialize the structure holding the solutions
  clpSolution solutions;
  // Set typeStruct field of ClpTrustedData struct to one.
  // This tells Clp it's "Mahdi!" 
  ClpTrustedData trustedSolutions;
  trustedSolutions.typeStruct = 1;
  trustedSolutions.data = &solutions;
  solutions.numberSolutions=0;
  solutions.maximumSolutions=0;
  solutions.numberColumns=simplex->numberColumns();
  solutions.solution=NULL;
  solutions.numberUnsatisfied=NULL;
  simplex->setTrustedUserPointer(&trustedSolutions); // rlh: old was
					      // "userPointer"
   // Solve from all slack to get some points
  simplex->allSlackBasis();
  simplex->primal();

  // rlh: Mahdi's code
  // -------------------------------------------------
  // Get the problem information
  
  // - Get the number of cols and rows
  int numCols = clpSolver->getNumCols();
  int numRows = clpSolver->getNumRows();

  // - Get the right hand side of the rows 
  const double * rhs = clpSolver->getRightHandSide();

  // Find the integer variables
  // rlh: consider using columnType 
  // One if not continuous (i.e., bin or gen int)
  bool * varClassInt = new bool[numCols];
  for(int i=0; i<numCols; i++)
    {
      if(clpSolver->isContinuous(i))
	varClassInt[i] = 0;
      else
	varClassInt[i] = 1;
    }

  // -Get the rows sense
  const char * rowSense;
  rowSense = clpSolver->getRowSense();

  // -Get the objective coefficients
  const double *objCoefficients = clpSolver->getObjCoefficients();
  double *originalObjCoeff = new double [numCols];
  for(int i=0; i<numCols;i++)
    originalObjCoeff[i] = objCoefficients[i];

  // -Get the matrix of the problem
  double ** matrix = new double * [numRows];
  for(int i = 0; i < numRows; i++)
    {
      matrix[i] = new double[numCols];
      for(int j = 0; j < numCols; j++)
	matrix[i][j] = 0;
    }
  
  const CoinPackedMatrix* matrixByRow = clpSolver->getMatrixByRow();
  const double * matrixElements = matrixByRow->getElements();
  const int * matrixIndices = matrixByRow->getIndices();
  const int * matrixStarts = matrixByRow->getVectorStarts();
  for(int j = 0; j < numRows; j++)
    {
      for(int i = matrixStarts[j]; i < matrixStarts[j+1]; i++)
	{
	  matrix[j][matrixIndices[i]] = matrixElements[i];
	}
    }
  
  double * newObj = new double [numCols];
  srand ( time(NULL) +1);
  int randNum;

  // Shuffle the rows:
  // Put the rows in a random order 
  // so that the optimal solution is a different corner point than the
  // starting point. 
  int * index = new int [numRows];
  for(int i=0; i<numRows; i++)
    index[i]=i;
  for(int i=0; i<numRows; i++)
    {
      int temp = index[i];
      int randNumTemp = i+(rand()%(numRows-i));
      index[i]=index[randNumTemp];
      index[randNumTemp]=temp;
    }
  
  // Start finding corner points by iteratively doing the following:
  // - find randomly tilted objective
  // - solve
  for(int i=0; i<numRows; i++)
    {
      // TODO: that 10,000 could be a param in the member data
      if(solutions.numberSolutions  > 10000)
	break;
      randNum = rand()%10 + 1;
      randNum = fmod(randNum, 2);
      for(int j=0; j<numCols; j++)
	{
	  // for row i and column j vary the coefficient "a bit"
	  if(randNum == 1)		
	    // if the element is zero, then round the coefficient down to 0.1
	    if(fabs(matrix[index[i]][j]) < 1e-6)
	      newObj[j] = 0.1;
	    else
	      // if the element is nonzero, then increase it "a bit"
	      newObj[j] = matrix[index[i]][j] * 1.1;
	  else
	    // if randnum is 2, then
	    // if the element is zero, then round the coefficient down
	    // to NEGATIVE 0.1
	    if(fabs(matrix[index[i]][j]) < 1e-6)
	      newObj[j] = -0.1;
	    else
	      // if the element is nonzero, then DEcrease it "a bit"
	      newObj[j] = matrix[index[i]][j] * 0.9;
	}
      // Use the new "tilted" objective
      clpSolver->setObjective(newObj);

      // Based on the row sense, we decide whether to max or min
      if(rowSense[i] == 'L')
	clpSolver->setObjSense(-1);
      else
	clpSolver->setObjSense(1);

      // Solve with primal simplex
      clpSolver->getModelPtr()->primal(1);
      printf("---------------------------------------------------------------- %d\n", i);
    }
  // Iteratively do this process until...
  // either you reach the max number of corner points (aka 10K)
  // or all the rows have been used as an objective.
  
  // Look at solutions
  int numberSolutions = solutions.numberSolutions;
  const char * integerInfo = simplex->integerInformation();
  const double * columnLower = simplex->columnLower();
  const double * columnUpper = simplex->columnUpper();
  printf("there are %d solutions\n",numberSolutions);

  // Up to here we have all the corner points
  // Now we need to do the random walks and roundings
  
  double ** cornerPoints = new double * [numberSolutions];
  for(int j=0;j<numberSolutions;j++)
    cornerPoints[j] = solutions.solution[j];
  
  bool feasibility = 1;
  double bestObj = 1e30;
  vector< vector <double> > feasibles;
  int numFeasibles = 0;
  
  // Check the feasibility of the corner points
  int numCornerPoints = numberSolutions;

  rhs = clpSolver->getRightHandSide();
  rowSense = clpSolver->getRowSense();
  
  for(int i=0; i<numCornerPoints; i++)
    {
      //get the objective value for this this point
      double objValue = 0;
      for(int k=0; k<numCols; k++)
	objValue += cornerPoints[i][k] * originalObjCoeff[k];
      
      if(objValue < bestObj)
	{
	  feasibility = 1;
	  for(int j=0; j<numCols; j++)
	    {
	      if(varClassInt[j])
		{
		  double closest=floor(cornerPoints[i][j]+0.5);
		  if(fabs(cornerPoints[i][j] - closest)>1e-6)
		    {
		      feasibility = 0;
		      break;
		    }
		}
	    }
	  if(feasibility)
	      { 
		for(int irow = 0; irow < numRows; irow++)
		  {
		    double lhs = 0;
		    for(int j = 0; j <numCols; j++)
		      {
			lhs += matrix[irow][j] * cornerPoints[i][j];
		      }
		    if(rowSense[irow] == 'L' && lhs > rhs[irow] + 1e-6)
		      {
			feasibility = 0;
			break;
		      }
		    if(rowSense[irow] == 'G' && lhs < rhs[irow] - 1e-6)
		      {
			feasibility = 0;
			break;
		      }
		    if(rowSense[irow] == 'E' && (lhs - rhs[irow] > 1e-6 || lhs - rhs[irow] < -1e-6)) 
		      {
			feasibility = 0;
			break;
		      }
		  }
	      }
	  
	  if(feasibility)
	    {		
	      numFeasibles++;
		feasibles.push_back(vector <double> (numCols));
		for(int k=0; k<numCols; k++)
		  feasibles[numFeasibles-1][k] = cornerPoints[i][k];
		printf("obj: %f\n", objValue);
		if(objValue < bestObj)
		  bestObj = objValue;
	    }
	}
    }	
  int numFeasibleCorners;
  numFeasibleCorners = numFeasibles;
  //find the center of gravity of the corner points as the first random point
  double * rp = new double[numCols];
  for(int i=0; i<numCols; i++)
    {
      rp[i] = 0;
      for(int j=0; j < numCornerPoints; j++)
	{
	    rp[i] += cornerPoints[j][i];
	}
      rp[i] = rp[i]/numCornerPoints;
    }
  
  //-------------------------------------------
  //main loop:
  // -generate the next random point
  // -round the random point
  // -check the feasibility of the random point
  //-------------------------------------------
  
  srand ( time(NULL) +1);
  int numRandomPoints = 0;
  while(numRandomPoints < 50000)
    {
      numRandomPoints++;
      //generate the next random point
      int randomIndex = rand()%numCornerPoints;
      double random = rand();
      random = random/RAND_MAX;
      for(int i=0; i<numCols; i++)
	{
	  rp[i] = (random * (cornerPoints[randomIndex][i] - rp[i])) + rp[i];
	}
      
      //CRISP ROUNDING
      //round the random point just generated
      double * roundRp = new double[numCols];
      for(int i=0; i<numCols; i++)
	{
	  roundRp[i] = rp[i];
	  if(varClassInt[i])
	    {
	      if(rp[i]>=0)
		{				
		  if(fmod(rp[i],1) > 0.5)
		    roundRp[i] = floor(rp[i]) + 1;
		  else
		    roundRp[i] = floor(rp[i]);
		}
	      else
		{
		  if(abs(fmod(rp[i],1)) > 0.5)
		    roundRp[i] = floor(rp[i]);
		  else
		    roundRp[i] = floor(rp[i])+1;
		  
		}
	    }
	}
      
      
      //SOFT ROUNDING
      // Look at original files for the "how to" on soft rounding 
      
      
      //check the feasibility of the rounded random point
      // -Check the feasibility
      // -Get the rows sense
      rowSense = clpSolver->getRowSense();
      rhs = clpSolver->getRightHandSide();
      
      //get the objective value for this feasible point
      double objValue = 0;
      for(int i=0; i<numCols; i++)
	objValue += roundRp[i] * originalObjCoeff[i];
      
      if(objValue<bestObj)
	{
	  feasibility = 1;
	  for(int i = 0; i < numRows; i++)
	    {
	      double lhs = 0;
	      for(int j = 0; j <numCols; j++)
		{
		  lhs += matrix[i][j] * roundRp[j];
		}
	      if(rowSense[i] == 'L' && lhs > rhs[i] + 1e-6)
		{
		  feasibility = 0;
		  break;
		}
	      if(rowSense[i] == 'G' && lhs < rhs[i] - 1e-6)
		{
		  feasibility = 0;
		  break;
		}
	      if(rowSense[i] == 'E' && (lhs - rhs[i] > 1e-6 || lhs - rhs[i] < -1e-6)) 
		{
		  feasibility = 0;
		  break;
		}
	    }
	  if(feasibility)
	    {
	      printf("Feasible Found!!\n");
	      printf("%.2f\n", CoinCpuTime()-start);
	      numFeasibles++;
	      feasibles.push_back(vector <double> (numCols));
	      for(int i=0; i<numCols; i++)
		feasibles[numFeasibles-1][i] = roundRp[i];
	      printf("obj: %f\n", objValue);
	      if(objValue < bestObj)
		bestObj = objValue;
	    }
	}
    }
  printf("Number of Feasible Corners: %d\n", numFeasibleCorners);
  printf("Number of Feasibles Found: %d\n", numFeasibles);
  if(numFeasibles > 0)
    printf("Best Objective: %f\n", bestObj);
  printf("time: %.2f\n",CoinCpuTime()-start);
  
  if (numFeasibles == 0) return 0; 
  
  // We found something better
  solutionValue = bestObj;
  for (int k = 0; k<numCols; k++){
    betterSolution[k] =  feasibles[numFeasibles-1][k];
  }  
  std::cout << "See you soon! You're leaving the Randomized Rounding Heuristic" << std::endl;
  return 1;

}
// update model
void CbcHeuristicRandRound::setModel(CbcModel * model)
{
  // probably same as resetModel
}

  
