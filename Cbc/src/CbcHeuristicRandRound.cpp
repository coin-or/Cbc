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

#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcHeuristicRandRound.hpp"

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
  Comments needed
  Returns 1 if solution, 0 if not */
int
CbcHeuristicRandRound::solution(double & solutionValue,
			 double * betterSolution)
{

  numCouldRun_++;
  int returnCode = 0;
  return returnCode;
}
// update model
void CbcHeuristicRandRound::setModel(CbcModel * model)
{
  // probably same as resetModel
}

  
