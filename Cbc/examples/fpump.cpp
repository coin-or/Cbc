/* $Id$ */
// Copyright (C) 2009, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CbcModel.hpp"
#include "CbcHeuristicFPump.hpp"
#include "CbcHeuristicRINS.hpp"
#include "CbcHeuristicLocal.hpp"
#include "CbcHeuristicDiveCoefficient.hpp"
#include "CbcHeuristicDiveFractional.hpp"

// Using as solver
#include "OsiClpSolverInterface.hpp"

int main (int argc, const char *argv[])
{
  OsiClpSolverInterface solver1;
  // Read in example model
  // and assert that it is a clean model
  std::string mpsFileName = "../../Data/Sample/p0033.mps";
  if (argc>=2) mpsFileName = argv[1];
  int numMpsReadErrors = solver1.readMps(mpsFileName.c_str(),"");
  assert(numMpsReadErrors==0);
  // make sure perturbation on (may want to make stronger??)
  solver1.getModelPtr()->setPerturbation(50);
  solver1.getModelPtr()->defaultFactorizationFrequency();

  solver1.initialSolve();
  double objValue = solver1.getObjValue();
  // Pass data and solver to CbcModel 
  CbcModel model(solver1);

  // uncomment to reduce printout
  //model.setLogLevel(1);
  //model.solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
  CbcHeuristicDiveCoefficient heuristicDC(model);
  heuristicDC.setHeuristicName("DiveCoefficient");
  // allow to exit if close enough to optimum
  heuristicDC.setSwitches(1);
  model.addHeuristic(&heuristicDC);
  CbcHeuristicDiveFractional heuristicDF(model);
  heuristicDF.setHeuristicName("DiveFractional");
  heuristicDF.setSwitches(1);
  model.addHeuristic(&heuristicDF);
  CbcHeuristicFPump pump(model);
  // allow to exit if close enough to optimum
  // plus some dubious options
  pump.setSwitches(1+4+8);
  pump.setMaximumTime(200);
  pump.setMaximumPasses(100);
  pump.setMaximumRetries(1); 
  pump.setFixOnReducedCosts(1);
  pump.setHeuristicName("Feasibility pump");
  pump.setFractionSmall(1.0);
  pump.setWhen(13);
  pump.setFakeCutoff(objValue+0.01*fabs(objValue));
  pump.setReducedCostMultiplier(0.1);
  pump.setFeasibilityPumpOptions(80);
  model.addHeuristic(&pump);
  pump.setFakeCutoff(objValue+0.05*fabs(objValue));
  pump.setFeasibilityPumpOptions(80);
  model.addHeuristic(&pump);
  pump.setFakeCutoff(objValue+0.01*fabs(objValue));
  pump.setReducedCostMultiplier(0.1);
  pump.setFeasibilityPumpOptions(80);
  model.addHeuristic(&pump);
  pump.setFakeCutoff(objValue+0.05*fabs(objValue));
  pump.setReducedCostMultiplier(1.0);
  pump.setFeasibilityPumpOptions(80);
  pump.setMaximumTime(200);
  model.addHeuristic(&pump);
  CbcHeuristicRINS rins(model);
  rins.setHeuristicName("RINS");
  rins.setFractionSmall(0.5);
  rins.setDecayFactor(5.0);
  model.addHeuristic(&rins) ;
  CbcHeuristicLocal local(model);
  local.setHeuristicName("LOCAL");
  local.setFractionSmall(0.5);
  local.setSearchType(1);
  model.addHeuristic(&local) ;
  model.setThreadMode(8);
  model.setNumberThreads(2);
  model.setMaximumNodes(-1);
  // Say what is good enough for solution
  model.setAllowableFractionGap(0.005);
  // Do root
  model.branchAndBound();
#if 0
  /* Print solution.  CbcModel clones solver so we
     need to get current copy */
  int numberColumns = model.solver()->getNumCols();
    
  const double * solution = model.solver()->getColSolution();
    
  for (int iColumn=0;iColumn<numberColumns;iColumn++) {
    double value=solution[iColumn];
    if (fabs(value)>1.0e-7&&model.solver()->isInteger(iColumn)) 
      printf("%d has value %g\n",iColumn,value);
  }
#endif
  return 0;
}    
