// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <iomanip>


// For Branch and bound
#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcStrategy.hpp"
#include "CbcCutGenerator.hpp"
#include "CoinTime.hpp"
#include "CbcHeuristicLocal.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "CbcHeuristicFPump.hpp"
#include "CbcHeuristicRINS.hpp"

//#############################################################################

/** User class
 */

class CbcStrategyUser : public CbcStrategyDefault {
public:

  // Default Constructor 
  CbcStrategyUser (bool cutsOnlyAtRoot=true,
                      int numberStrong=5,
                      int numberBeforeTrust=0,
                      int printLevel=0);

  // Copy constructor 
  CbcStrategyUser ( const CbcStrategyUser &);
   
  // Destructor 
  ~CbcStrategyUser ();
  
  /// Clone
  virtual CbcStrategy * clone() const;
  /// Setup heuristics
  virtual void setupHeuristics(CbcModel & model);

protected:
  // Data

private:
  /// Illegal Assignment operator 
  CbcStrategyUser & operator=(const CbcStrategyUser& rhs);
};

// Default Constructor
CbcStrategyUser::CbcStrategyUser(bool cutsOnlyAtRoot,
                                       int numberStrong,
                                       int numberBeforeTrust,
                                       int printLevel)
  :CbcStrategyDefault(cutsOnlyAtRoot,numberStrong,
		      numberBeforeTrust,printLevel)
{
}


// Destructor 
CbcStrategyUser::~CbcStrategyUser ()
{
}

// Clone
CbcStrategy *
CbcStrategyUser::clone() const
{
  return new CbcStrategyUser(*this);
}

// Copy constructor 
CbcStrategyUser::CbcStrategyUser(const CbcStrategyUser & rhs)
:
  CbcStrategyDefault(rhs)
{
}
// Setup heuristics
void 
CbcStrategyUser::setupHeuristics(CbcModel & model)
{
  // FPump done first as it only works if no solution
  int numberHeuristics = model.numberHeuristics();
  int iHeuristic;
  bool found;
  found=false;
  for (iHeuristic=0;iHeuristic<numberHeuristics;iHeuristic++) {
    CbcHeuristic * heuristic = model.heuristic(iHeuristic);
    CbcHeuristicFPump * cgl = dynamic_cast<CbcHeuristicFPump *>(heuristic);
    if (cgl) {
      found=true;
      break;
    }
  }
  if (!found) {
    CbcHeuristicFPump heuristic(model);
    heuristic.setMaximumPasses(20);
    double value = model.solver()->getObjSense()*model.solver()->getObjValue();
    // also set increment
    heuristic.setAbsoluteIncrement(0.005*(fabs(value)+1.0e-12));
    heuristic.setAccumulate(0);
    heuristic.setMaximumRetries(2);
    heuristic.setWhen(13);
    heuristic.setHeuristicName("feasibility pump");
    model.addHeuristic(&heuristic);
  }
  // Allow CbcStrategyDefault heuristics
  CbcStrategyDefault::setupHeuristics(model);

  // Add if not there 
  
  numberHeuristics = model.numberHeuristics();
  found=false;
  for (iHeuristic=0;iHeuristic<numberHeuristics;iHeuristic++) {
    CbcHeuristic * heuristic = model.heuristic(iHeuristic);
    CbcHeuristicLocal * cgl = dynamic_cast<CbcHeuristicLocal *>(heuristic);
    if (cgl) {
      found=true;
      break;
    }
  }
  if (!found) {
    CbcHeuristicLocal heuristic(model);
    heuristic.setSearchType(1);
    heuristic.setHeuristicName("join solutions");
    model.addHeuristic(&heuristic);
  }
  
  // Add if not there 
  
  numberHeuristics = model.numberHeuristics();
  found=false;
  for (iHeuristic=0;iHeuristic<numberHeuristics;iHeuristic++) {
    CbcHeuristic * heuristic = model.heuristic(iHeuristic);
    CbcHeuristicGreedyCover * cgl = dynamic_cast<CbcHeuristicGreedyCover *>(heuristic);
    if (cgl) {
      found=true;
      break;
    }
  }
  if (!found) {
    CbcHeuristicGreedyCover heuristic(model);
    heuristic.setHeuristicName("greedy cover");
    model.addHeuristic(&heuristic);
  }
  
  // Add if not there 
  
  numberHeuristics = model.numberHeuristics();
  found=false;
  for (iHeuristic=0;iHeuristic<numberHeuristics;iHeuristic++) {
    CbcHeuristic * heuristic = model.heuristic(iHeuristic);
    CbcHeuristicGreedyEquality * cgl = dynamic_cast<CbcHeuristicGreedyEquality *>(heuristic);
    if (cgl) {
      found=true;
      break;
    }
  }
  if (!found) {
    CbcHeuristicGreedyEquality heuristic(model);
    heuristic.setHeuristicName("greedy equality");
    model.addHeuristic(&heuristic);
  }
  
  // Add if not there 
  
  numberHeuristics = model.numberHeuristics();
  found=false;
  for (iHeuristic=0;iHeuristic<numberHeuristics;iHeuristic++) {
    CbcHeuristic * heuristic = model.heuristic(iHeuristic);
    CbcHeuristicRINS * cgl = dynamic_cast<CbcHeuristicRINS *>(heuristic);
    if (cgl) {
      found=true;
      break;
    }
  }
  if (!found) {
    CbcHeuristicRINS heuristic(model);
    heuristic.setHeuristicName("RINS");
    model.addHeuristic(&heuristic);
  }
}

/************************************************************************

This main program reads in an integer model from an mps file.

It then sets up a strategy and solves a problem

This inherits from CbcStrategyUser and adds in other heuristics
************************************************************************/

int main (int argc, const char *argv[])
{

  // Define your favorite OsiSolver
  
  OsiClpSolverInterface solver1;

  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string mpsFileName = "../../Data/Sample/p0033.mps";
  if (argc>=2) mpsFileName = argv[1];
  int numMpsReadErrors = solver1.readMps(mpsFileName.c_str(),"");
  assert(numMpsReadErrors==0);
  double time1 = CoinCpuTime();

  solver1.initialSolve();
  // Reduce printout
  solver1.setHintParam(OsiDoReducePrint,true,OsiHintTry);
  // See if we want preprocessing
  OsiSolverInterface * solver2=&solver1;
  CbcModel model(solver1);
  // Stop after 20 minutes
  int minutes=20;
  std::cout<<"Stopping after "<<minutes<<" minutes"<<std::endl;
  model.setDblParam(CbcModel::CbcMaximumSeconds,60.0*minutes);
  ////////////////////////////////////////////////////////////

  // Default strategy will leave cut generators as they exist already
  // so cutsOnlyAtRoot (1) true
  // numberStrong (2) is 5 (default)
  // numberBeforeTrust (3) is 0 (default)
  // printLevel (4) defaults (0)
  // So same as CbcStrategyUser strategy(true,5,0,0);
  // NO - add more cuts CbcStrategyUser strategy;
  CbcStrategyUser strategy(false);
  // Set up pre-processing if wanted
  strategy.setupPreProcessing(1);
  model.setStrategy(strategy);

  ////////////////////////////////////////////////////////////
  // Do complete search
  
  model.branchAndBound();

  std::cout<<mpsFileName<<" took "<<CoinCpuTime()-time1<<" seconds, "
	   <<model.getNodeCount()<<" nodes with objective "
	   <<model.getObjValue()
	   <<(!model.status() ? " Finished" : " Not finished")
	   <<std::endl;

  // Print more statistics
  std::cout<<"Cuts at root node changed objective from "<<model.getContinuousObjective()
	   <<" to "<<model.rootObjectiveAfterCuts()<<std::endl;

  for (int iGenerator=0;iGenerator<model.numberCutGenerators();iGenerator++) {
    CbcCutGenerator * generator = model.cutGenerator(iGenerator);
    std::cout<<generator->cutGeneratorName()<<" was tried "
	     <<generator->numberTimesEntered()<<" times and created "
	     <<generator->numberCutsInTotal()<<" cuts of which "
	     <<generator->numberCutsActive()<<" were active after adding rounds of cuts";
    if (generator->timing())
      std::cout<<" ( "<<generator->timeInCutGenerator()<<" seconds)"<<std::endl;
    else
      std::cout<<std::endl;
  }
  // Print solution if finished - we can't get names from Osi! - so get from OsiClp

  if (model.getMinimizationObjValue()<1.0e50) {
    OsiSolverInterface * solver = model.solver();
    int numberColumns = solver->getNumCols();
    
    const double * solution = solver->getColSolution();

    // Get names from solver1 (as OsiSolverInterface may lose)
    std::vector<std::string> columnNames = *solver1.getModelPtr()->columnNames();
    
    int iColumn;
    std::cout<<std::setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14);
    
    std::cout<<"--------------------------------------"<<std::endl;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      double value=solution[iColumn];
      if (fabs(value)>1.0e-7&&solver->isInteger(iColumn)) 
	std::cout<<std::setw(6)<<iColumn<<" "
                 <<columnNames[iColumn]<<" "
                 <<value<<std::endl;
    }
    std::cout<<"--------------------------------------"<<std::endl;
  
    std::cout<<std::resetiosflags(std::ios::fixed|std::ios::showpoint|std::ios::scientific);
  }
  return 0;
}    
