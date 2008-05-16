// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CbcHeuristicDivePseudoCost.hpp"
#include "CbcStrategy.hpp"

// Default Constructor
CbcHeuristicDivePseudoCost::CbcHeuristicDivePseudoCost() 
  :CbcHeuristicDive()
{
}

// Constructor from model
CbcHeuristicDivePseudoCost::CbcHeuristicDivePseudoCost(CbcModel & model)
  :CbcHeuristicDive(model)
{
}

// Destructor 
CbcHeuristicDivePseudoCost::~CbcHeuristicDivePseudoCost ()
{
}

// Clone
CbcHeuristicDivePseudoCost *
CbcHeuristicDivePseudoCost::clone() const
{
  return new CbcHeuristicDivePseudoCost(*this);
}

// Create C++ lines to get to current state
void 
CbcHeuristicDivePseudoCost::generateCpp( FILE * fp) 
{
  CbcHeuristicDivePseudoCost other;
  fprintf(fp,"0#include \"CbcHeuristicDivePseudoCost.hpp\"\n");
  fprintf(fp,"3  CbcHeuristicDivePseudoCost heuristicDivePseudoCost(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp,"heuristicDivePseudoCost");
  fprintf(fp,"3  cbcModel->addHeuristic(&heuristicDivePseudoCost);\n");
}

// Copy constructor 
CbcHeuristicDivePseudoCost::CbcHeuristicDivePseudoCost(const CbcHeuristicDivePseudoCost & rhs)
:
  CbcHeuristicDive(rhs)
{
}

// Assignment operator 
CbcHeuristicDivePseudoCost & 
CbcHeuristicDivePseudoCost::operator=( const CbcHeuristicDivePseudoCost& rhs)
{
  if (this!=&rhs) {
    CbcHeuristicDive::operator=(rhs);
  }
  return *this;
}

bool
CbcHeuristicDivePseudoCost::selectVariableToBranch(OsiSolverInterface* solver,
						   const double* newSolution,
						   int& bestColumn,
						   int& bestRound)
{
  int numberIntegers = model_->numberIntegers();
  const int * integerVariable = model_->integerVariable();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);

  // get the LP relaxation solution at the root node
  double * rootNodeLPSol = model_->continuousSolution();

  // get pseudo costs
  double * pseudoCostDown = new double[numberIntegers];
  double * pseudoCostUp = new double[numberIntegers];
  model_->fillPseudoCosts(pseudoCostDown, pseudoCostUp);

  bestColumn = -1;
  bestRound = -1; // -1 rounds down, +1 rounds up
  double bestScore = -1.0;
  bool allTriviallyRoundableSoFar = true;
  for (int i=0; i<numberIntegers; i++) {
    int iColumn = integerVariable[i];
    double rootValue=rootNodeLPSol[iColumn];
    double value=newSolution[iColumn];
    double fraction=value-floor(value);
    int round = 0;
    if (fabs(floor(value+0.5)-value)>integerTolerance) {
      if (allTriviallyRoundableSoFar||(downLocks_[i]>0&&upLocks_[i]>0)) {

	if (allTriviallyRoundableSoFar&&downLocks_[i]>0&&upLocks_[i]>0) {
	  allTriviallyRoundableSoFar = false;
	  bestScore = -1.0;
	}

	double pCostDown = pseudoCostDown[i];
	double pCostUp = pseudoCostUp[i];
	assert(pCostDown >= 0.0 && pCostUp >= 0.0);

	if(allTriviallyRoundableSoFar&&downLocks_[i]==0&&upLocks_[i]>0)
	    round = 1;
	else if(allTriviallyRoundableSoFar&&downLocks_[i]>0&&upLocks_[i]==0)
	  round = -1;
	else if(value - rootValue < -0.4)
	  round = -1;
	else if(value - rootValue > 0.4)
	  round = 1;
	else if(fraction < 0.3)
	  round = -1;
	else if(fraction > 0.7)
	  round = 1;
	else if(pCostDown < pCostUp)
	  round = -1;
	else
	  round = 1;

	// calculate score
	double score;
	if(round == 1)
	  score = fraction * (pCostDown + 1.0) / (pCostUp + 1.0);
	else
	  score = (1.0 - fraction) * (pCostUp + 1.0) / (pCostDown + 1.0);

	// if variable is binary, increase its chance of being selected
	if(solver->isBinary(iColumn))
	  score *= 1000.0;
	
	if(score > bestScore) {
	  bestColumn = iColumn;
	  bestScore = score;
	  bestRound = round;
	}
      }
    }
  }

  delete [] pseudoCostDown;
  delete [] pseudoCostUp;

  return allTriviallyRoundableSoFar;
}
