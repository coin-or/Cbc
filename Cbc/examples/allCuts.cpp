// $Id$
// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>
#include <iomanip>

#include "CoinPragma.hpp"
// For Branch and bound
#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcBranchUser.hpp"
#include "CbcCompareUser.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcStrategy.hpp"
#include "CbcHeuristicLocal.hpp"
#include "OsiClpSolverInterface.hpp"

// Cuts

#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglRedSplit.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding2.hpp"
// Preprocessing
#include "CglPreProcess.hpp"

// For saying about solution validity
#include "OsiAuxInfo.hpp"

// Heuristics (But any will have to be special)

#include "CbcHeuristic.hpp"

#include "CoinTime.hpp"
// Need stored cuts

#include "CglStored.hpp"

/** Stored Cut Generator Class */
class CglStoredUser : public CglStored {

public:
  /**@name Generate Cuts */
  //@{
  /** Generate Mixed Integer Stored cuts for the model of the 
      solver interface, si.

      Insert the generated cuts into OsiCut, cs.

      This generator just looks at previously stored cuts
      and inserts any that are violated by enough
  */
  virtual void generateCuts(const OsiSolverInterface &si, OsiCuts &cs,
    const CglTreeInfo info = CglTreeInfo()) const;
  //@}

  /**@name Cut stuff */
  //@{
  OsiRowCut *mutableRowCutPointer(int index)
  {
    return cuts_.rowCutPtr(index);
  }
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor
  CglStoredUser();

  /// Copy constructor
  CglStoredUser(const CglStoredUser &rhs);

  /// Clone
  virtual CglCutGenerator *clone() const;

  /// Assignment operator
  CglStoredUser &
  operator=(const CglStoredUser &rhs);

  /// Destructor
  virtual ~CglStoredUser();
  //@}

protected:
  // Protected member methods

  // Protected member data

  /**@name Protected member data */
  //@{
  /** Don't add any more cuts after this number passes (per node)
      unless looks integer feasible
  */
  int numberPasses_;
  //@}
};
//-------------------------------------------------------------------
// Generate Stored cuts
//-------------------------------------------------------------------
void CglStoredUser::generateCuts(const OsiSolverInterface &si, OsiCuts &cs,
  const CglTreeInfo info) const
{
  // Get basic problem information
  const double *solution = si.getColSolution();
  if (info.inTree && info.pass > numberPasses_) {
    // only continue if integer feasible
    int numberColumns = si.getNumCols();
    int i;
    const double *colUpper = si.getColUpper();
    const double *colLower = si.getColLower();
    int numberAway = 0;
    for (i = 0; i < numberColumns; i++) {
      double value = solution[i];
      // In case slightly away from bounds
      value = CoinMax(colLower[i], value);
      value = CoinMin(colUpper[i], value);
      if (si.isInteger(i) && fabs(value - fabs(value + 0.5)) > 1.0e-5)
        numberAway++;
    }
    if (numberAway)
      return; // let code branch
  }
  int numberRowCuts = cuts_.sizeRowCuts();
  for (int i = 0; i < numberRowCuts; i++) {
    const OsiRowCut *rowCutPointer = cuts_.rowCutPtr(i);
    double violation = rowCutPointer->violated(solution);
    if (violation >= requiredViolation_)
      cs.insert(*rowCutPointer);
  }
}

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
CglStoredUser::CglStoredUser()
  : CglStored()
  , numberPasses_(5)
{
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
CglStoredUser::CglStoredUser(const CglStoredUser &source)
  : CglStored(source)
  , numberPasses_(source.numberPasses_)
{
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglStoredUser::clone() const
{
  return new CglStoredUser(*this);
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
CglStoredUser::~CglStoredUser()
{
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
CglStoredUser &
CglStoredUser::operator=(const CglStoredUser &rhs)
{
  if (this != &rhs) {
    CglStored::operator=(rhs);
    numberPasses_ = rhs.numberPasses_;
  }
  return *this;
}
// Class to disallow strong branching solutions
#include "CbcFeasibilityBase.hpp"
class CbcFeasibilityNoStrong : public CbcFeasibilityBase {
public:
  // Default Constructor
  CbcFeasibilityNoStrong() {};

  virtual ~CbcFeasibilityNoStrong() {};
  // Copy constructor
  CbcFeasibilityNoStrong(const CbcFeasibilityNoStrong &rhs) {};

  // Assignment operator
  CbcFeasibilityNoStrong &operator=(const CbcFeasibilityNoStrong &rhs)
  {
    return *this;
  };

  /// Clone
  virtual CbcFeasibilityBase *clone() const
  {
    return new CbcFeasibilityNoStrong();
  };

  /**
     On input mode:
     0 - called after a solve but before any cuts
     -1 - called after strong branching
     Returns :
     0 - no opinion
     -1 pretend infeasible
     1 pretend integer solution
  */
  virtual int feasible(CbcModel *model, int mode)
  {
    return mode;
  };
};

//#############################################################################

/************************************************************************

This main program reads in an integer model from an mps file.

It makes L or G rows into cuts

************************************************************************/

int main(int argc, const char *argv[])
{

  // Define your favorite OsiSolver

  OsiClpSolverInterface solver1;

  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string mpsFileName;
#if defined(SAMPLEDIR)
  mpsFileName = SAMPLEDIR "/p0033.mps";
#else
  if (argc < 2) {
    fprintf(stderr, "Do not know where to find sample MPS files.\n");
    exit(1);
  }
#endif
  if (argc >= 2)
    mpsFileName = argv[1];
  int numMpsReadErrors = solver1.readMps(mpsFileName.c_str(), "");
  if (numMpsReadErrors != 0) {
    printf("%d errors reading MPS file\n", numMpsReadErrors);
    return numMpsReadErrors;
  }
  double time1 = CoinCpuTime();
  OsiClpSolverInterface solverSave = solver1;

  /* Options are:
     preprocess to do preprocessing
     time in minutes
     if 2 parameters and numeric taken as time
  */
  bool preProcess = false;
  double minutes = -1.0;
  int nGoodParam = 0;
  for (int iParam = 2; iParam < argc; iParam++) {
    if (!strcmp(argv[iParam], "preprocess")) {
      preProcess = true;
      nGoodParam++;
    } else if (!strcmp(argv[iParam], "time")) {
      if (iParam + 1 < argc && isdigit(argv[iParam + 1][0])) {
        minutes = atof(argv[iParam + 1]);
        if (minutes >= 0.0) {
          nGoodParam += 2;
          iParam++; // skip time
        }
      }
    }
  }
  if (nGoodParam == 0 && argc == 3 && isdigit(argv[2][0])) {
    // If time is given then stop after that number of minutes
    minutes = atof(argv[2]);
    if (minutes >= 0.0)
      nGoodParam = 1;
  }
  if (nGoodParam != argc - 2 && argc >= 2) {
    printf("Usage <file> [preprocess] [time <minutes>] or <file> <minutes>\n");
    exit(1);
  }
  // Reduce printout
  solver1.setHintParam(OsiDoReducePrint, true, OsiHintTry);
  // See if we want preprocessing
  OsiSolverInterface *solver2 = &solver1;
  CglPreProcess process;
  // Never do preprocessing until dual tests out as can fix incorrectly
  preProcess = false;
  if (preProcess) {
    /* Do not try and produce equality cliques and
       do up to 5 passes */
    solver2 = process.preProcess(solver1, false, 5);
    if (!solver2) {
      printf("Pre-processing says infeasible\n");
      exit(2);
    }
    solver2->resolve();
  }
  // Turn L rows into cuts
  CglStoredUser stored;
  {
    int numberRows = solver2->getNumRows();

    int *whichRow = new int[numberRows];
    // get row copy
    const CoinPackedMatrix *rowCopy = solver2->getMatrixByRow();
    const int *column = rowCopy->getIndices();
    const int *rowLength = rowCopy->getVectorLengths();
    const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
    const double *rowLower = solver2->getRowLower();
    const double *rowUpper = solver2->getRowUpper();
    const double *element = rowCopy->getElements();
    int iRow, nDelete = 0;
    for (iRow = 0; iRow < numberRows; iRow++) {
      if (rowLower[iRow] < -1.0e20 || rowUpper[iRow] > 1.0e20) {
        // take out
        whichRow[nDelete++] = iRow;
      }
    }
    // leave some rows to avoid empty problem (Gomory does not like)
    nDelete = CoinMax(CoinMin(nDelete, numberRows - 5), 0);
    for (int jRow = 0; jRow < nDelete; jRow++) {
      iRow = whichRow[jRow];
      int start = rowStart[iRow];
      stored.addCut(rowLower[iRow], rowUpper[iRow], rowLength[iRow],
        column + start, element + start);
    }
    /* The following is problem specific.
     Normally cuts are deleted if slack on cut basic.
     On some problems you may wish to leave cuts in as long
     as slack value zero
  */
    int numberCuts = stored.sizeRowCuts();
    for (int iCut = 0; iCut < numberCuts; iCut++) {
      //stored.mutableRowCutPointer(iCut)->setEffectiveness(1.0e50);
    }
    solver2->deleteRows(nDelete, whichRow);
    delete[] whichRow;
  }
  CbcModel model(*solver2);
  model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
  // Set up some cut generators and defaults
  // Probing first as gets tight bounds on continuous

  CglProbing generator1;
  generator1.setUsingObjective(true);
  generator1.setMaxPass(1);
  generator1.setMaxPassRoot(5);
  // Number of unsatisfied variables to look at
  generator1.setMaxProbe(10);
  generator1.setMaxProbeRoot(1000);
  // How far to follow the consequences
  generator1.setMaxLook(50);
  generator1.setMaxLookRoot(500);
  // Only look at rows with fewer than this number of elements
  generator1.setMaxElements(200);
  generator1.setRowCuts(3);

  CglGomory generator2;
  // try larger limit
  generator2.setLimit(300);

  CglKnapsackCover generator3;

  CglRedSplit generator4;
  // try larger limit
  generator4.setLimit(200);

  CglClique generator5;
  generator5.setStarCliqueReport(false);
  generator5.setRowCliqueReport(false);

  CglMixedIntegerRounding2 mixedGen;
  CglFlowCover flowGen;

  // Add in generators
  // Experiment with -1 and -99 etc
  // This is just for one particular model
  model.addCutGenerator(&generator1, -1, "Probing");
  //model.addCutGenerator(&generator2,-1,"Gomory");
  model.addCutGenerator(&generator2, 1, "Gomory");
  model.addCutGenerator(&generator3, -1, "Knapsack");
  // model.addCutGenerator(&generator4,-1,"RedSplit");
  //model.addCutGenerator(&generator5,-1,"Clique");
  model.addCutGenerator(&generator5, 1, "Clique");
  model.addCutGenerator(&flowGen, -1, "FlowCover");
  model.addCutGenerator(&mixedGen, -1, "MixedIntegerRounding");
  // Add stored cuts (making sure at all depths)
  model.addCutGenerator(&stored, 1, "Stored", true, false, false, -100, 1, -1);

  int numberGenerators = model.numberCutGenerators();
  int iGenerator;
  // Say we want timings
  for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
    CbcCutGenerator *generator = model.cutGenerator(iGenerator);
    generator->setTiming(true);
  }
  OsiClpSolverInterface *osiclp = dynamic_cast< OsiClpSolverInterface * >(model.solver());
  // go faster stripes
  if (osiclp) {
    if (osiclp->getNumRows() < 300 && osiclp->getNumCols() < 500) {
      //osiclp->setupForRepeatedUse(2,0);
      osiclp->setupForRepeatedUse(0, 0);
    }
    // Don't allow dual stuff
    osiclp->setSpecialOptions(osiclp->specialOptions() | 262144);
  }
  // Uncommenting this should switch off all CBC messages
  // model.messagesPointer()->setDetailMessages(10,10000,NULL);
  // No heuristics
  // Do initial solve to continuous
  model.initialSolve();
  /*  You need the next few lines -
      a) so that cut generator will always be called again if it generated cuts
      b) it is known that matrix is not enough to define problem so do cuts even
         if it looks integer feasible at continuous optimum.
      c) a solution found by strong branching will be ignored.
      d) don't recompute a solution once found
  */
  // Make sure cut generator called correctly (a)
  iGenerator = numberGenerators - 1;
  model.cutGenerator(iGenerator)->setMustCallAgain(true);
  // Say cuts needed at continuous (b)
  OsiBabSolver oddCuts;
  oddCuts.setSolverType(4);
  // owing to bug must set after initialSolve
  model.passInSolverCharacteristics(&oddCuts);
  // Say no to all solutions by strong branching (c)
  CbcFeasibilityNoStrong noStrong;
  model.setProblemFeasibility(noStrong);
  // Say don't recompute solution d)
  model.setSpecialOptions(4);

  // Could tune more
  double objValue = model.solver()->getObjSense() * model.solver()->getObjValue();
  double minimumDropA = CoinMin(1.0, fabs(objValue) * 1.0e-3 + 1.0e-4);
  double minimumDrop = fabs(objValue) * 1.0e-4 + 1.0e-4;
  printf("min drop %g (A %g)\n", minimumDrop, minimumDropA);
  model.setMinimumDrop(minimumDrop);

  if (model.getNumCols() < 500)
    model.setMaximumCutPassesAtRoot(-100); // always do 100 if possible
  else if (model.getNumCols() < 5000)
    model.setMaximumCutPassesAtRoot(100); // use minimum drop
  else
    model.setMaximumCutPassesAtRoot(20);
  model.setMaximumCutPasses(10);
  //model.setMaximumCutPasses(2);

  // Switch off strong branching if wanted
  // model.setNumberStrong(0);
  // Do more strong branching if small
  if (model.getNumCols() < 5000)
    model.setNumberStrong(10);
  model.setNumberStrong(20);
  //model.setNumberStrong(5);
  model.setNumberBeforeTrust(5);

  model.solver()->setIntParam(OsiMaxNumIterationHotStart, 100);

  // If time is given then stop after that number of minutes
  if (minutes >= 0.0) {
    std::cout << "Stopping after " << minutes << " minutes" << std::endl;
    model.setDblParam(CbcModel::CbcMaximumSeconds, 60.0 * minutes);
  }
  // Switch off most output
  if (model.getNumCols() < 30000) {
    model.messageHandler()->setLogLevel(1);
    //model.solver()->messageHandler()->setLogLevel(0);
  } else {
    model.messageHandler()->setLogLevel(2);
    model.solver()->messageHandler()->setLogLevel(1);
  }
  //model.messageHandler()->setLogLevel(2);
  //model.solver()->messageHandler()->setLogLevel(2);
  //model.setPrintFrequency(50);
  //#define DEBUG_CUTS
#ifdef DEBUG_CUTS
  // Set up debugger by name (only if no preprocesing)
  if (!preProcess) {
    std::string problemName;
    model.solver()->getStrParam(OsiProbName, problemName);
    model.solver()->activateRowCutDebugger(problemName.c_str());
  }
#endif
  // Do complete search

  model.branchAndBound();

  std::cout << mpsFileName << " took " << CoinCpuTime() - time1 << " seconds, "
            << model.getNodeCount() << " nodes with objective "
            << model.getObjValue()
            << (!model.status() ? " Finished" : " Not finished")
            << std::endl;

  // Print more statistics
  std::cout << "Cuts at root node changed objective from " << model.getContinuousObjective()
            << " to " << model.rootObjectiveAfterCuts() << std::endl;

  for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
    CbcCutGenerator *generator = model.cutGenerator(iGenerator);
    std::cout << generator->cutGeneratorName() << " was tried "
              << generator->numberTimesEntered() << " times and created "
              << generator->numberCutsInTotal() << " cuts of which "
              << generator->numberCutsActive() << " were active after adding rounds of cuts";
    if (generator->timing())
      std::cout << " ( " << generator->timeInCutGenerator() << " seconds)" << std::endl;
    else
      std::cout << std::endl;
  }
  // Print solution if finished - we can't get names from Osi! - so get from OsiClp

  if (model.getMinimizationObjValue() < 1.0e50) {
    // post process
    OsiSolverInterface *solver;
    if (preProcess) {
      process.postProcess(*model.solver());
      // Solution now back in solver1
      solver = &solver1;
    } else {
      solver = model.solver();
    }
    int numberColumns = solver->getNumCols();

    const double *solution = solver->getColSolution();

    // Get names from solver1 (as OsiSolverInterface may lose)
    std::vector< std::string > columnNames = *solver1.getModelPtr()->columnNames();

    int iColumn;
    std::cout << std::setiosflags(std::ios::fixed | std::ios::showpoint) << std::setw(14);

    std::cout << "--------------------------------------" << std::endl;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value = solution[iColumn];
      if (fabs(value) > 1.0e-7 && solver->isInteger(iColumn)) {
        std::cout << std::setw(6) << iColumn << " "
                  << columnNames[iColumn] << " "
                  << value << std::endl;
        solverSave.setColLower(iColumn, value);
        solverSave.setColUpper(iColumn, value);
      }
    }
    std::cout << "--------------------------------------" << std::endl;

    std::cout << std::resetiosflags(std::ios::fixed | std::ios::showpoint | std::ios::scientific);
    solverSave.initialSolve();
  }
  return 0;
}
