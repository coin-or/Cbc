// $Id$
// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>

#include "CoinPragma.hpp"

// For Branch and bound
#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcCutGenerator.hpp"
#include "CoinHelperFunctions.hpp"
#include "CbcStrategy.hpp"

// Need stored cuts

#include "CglStored.hpp"

// For saying about solution validity
#include "OsiAuxInfo.hpp"

// Time
#include "CoinTime.hpp"
// Class to disallow strong branching solutions
#include "CbcFeasibilityBase.hpp"
class CbcFeasibilityNoStrong : public CbcFeasibilityBase {
public:
  // Default Constructor
  CbcFeasibilityNoStrong() {}

  virtual ~CbcFeasibilityNoStrong() {}
  // Copy constructor
  CbcFeasibilityNoStrong(const CbcFeasibilityNoStrong &rhs) {}

  // Assignment operator
  CbcFeasibilityNoStrong &operator=(const CbcFeasibilityNoStrong &rhs)
  {
    return *this;
  }

  /// Clone
  virtual CbcFeasibilityBase *clone() const
  {
    return new CbcFeasibilityNoStrong();
  }

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
  }
};

/************************************************************************

This main program solves the following 0-1 problem:

min -x0 - 2x1 - 3x2 - 4x3

subject to 
 
x0 + x1 + x2 + x3 <= 2

and quadratic constraints with positive random numbers 

It does it creating extra yij variables and constraints xi + xj -1 <= yij
and putting quadratic elements on y

The extra constraints are treated as stored cuts.

This is to show how to keep branching even if we have a solution

************************************************************************/

int main(int argc, const char *argv[])
{

  // Define a Solver which inherits from OsiClpsolverInterface -> OsiSolverInterface

  OsiClpSolverInterface solver1;

  int nX = 4;

  int nY = (nX * (nX - 1) / 2);
  // All columns
  double *obj = new double[nX + nY];
  double *clo = new double[nX + nY];
  double *cup = new double[nX + nY];
  int i;
  for (i = 0; i < nX; i++) {
    obj[i] = -(i + 1);
    clo[i] = 0.0;
    cup[i] = 1.0;
  }
  for (i = nX; i < nX + nY; i++) {
    obj[i] = 0.0;
    clo[i] = 0.0;
    cup[i] = 1.0;
  }
  // Just ordinary rows
  int nRow = 1 + nX;
  double *rlo = new double[nRow];
  double *rup = new double[nRow];
  for (i = 0; i < nRow; i++) {
    rlo[i] = -COIN_DBL_MAX;
    rup[i] = 1.0;
  }
  // and first row
  rup[0] = nX / 2.0;
  // Matrix
  int nEl = nX + nX * nX;
  int *row = new int[nEl];
  int *col = new int[nEl];
  double *el = new double[nEl];
  // X
  nEl = 0;
  // May need scale to make plausible
  double scaleFactor = 1.0;
  for (i = 0; i < nX; i++) {
    row[nEl] = 0;
    col[nEl] = i;
    el[nEl++] = 1.0;
    // and diagonal
    row[nEl] = i + 1;
    col[nEl] = i;
    double value = CoinDrand48() * scaleFactor;
    // make reasonable (so multiples of 0.000001)
    value *= 1.0e6;
    int iValue = (int)(value + 1.0);
    value = iValue;
    value *= 1.0e-6;
    el[nEl++] = value;
  }
  // Y
  nY = nX;
  // And stored cuts
  CglStored stored;
  double cutEls[3] = { 1.0, 1.0, -1.0 };
  int cutIndices[3];
  for (i = 0; i < nX; i++) {
    cutIndices[0] = i;
    for (int j = i + 1; j < nX; j++) {
      cutIndices[1] = j;
      cutIndices[2] = nY;
      // add cut
      stored.addCut(-COIN_DBL_MAX, 1.0, 3, cutIndices, cutEls);
      row[nEl] = i + 1;
      col[nEl] = nY;
      double value = CoinDrand48() * scaleFactor;
      // multiply to make ones with most negative objective violated
      // make reasonable (so multiples of 0.000001)
      value *= 1.0e6 + 1.0e6 * j;
      int iValue = (int)(value + 1.0);
      value = iValue;
      value *= 1.0e-6;
      el[nEl++] = value;
      // and other
      if (i != j) {
        row[nEl] = j + 1;
        col[nEl] = nY;
        el[nEl++] = value;
      }
      nY++;
    }
  }
  // Create model
  CoinPackedMatrix matrix(true, row, col, el, nEl);
  solver1.loadProblem(matrix, clo, cup, obj, rlo, rup);
  delete[] obj;
  delete[] clo;
  delete[] cup;
  delete[] rlo;
  delete[] rup;
  delete[] row;
  delete[] col;
  delete[] el;
  // Integers
  for (i = 0; i < nX; i++)
    solver1.setInteger(i);
  // Reduce printout
  solver1.setHintParam(OsiDoReducePrint, true, OsiHintTry);
  // This clones solver
  CbcModel model(solver1);
  // Add stored cuts (making sure at all depths)
  model.addCutGenerator(&stored, 1, "Stored", true, false, false, -100, 1, -1);
  /*  You need the next few lines -
      a) so that cut generator will always be called again if it generated cuts
      b) it is known that matrix is not enough to define problem so do cuts even
         if it looks integer feasible at continuous optimum.
      c) a solution found by strong branching will be ignored.
      d) don't recompute a solution once found
  */
  // Make sure cut generator called correctly (a)
  model.cutGenerator(0)->setMustCallAgain(true);
  // Say cuts needed at continuous (b)
  OsiBabSolver oddCuts;
  oddCuts.setSolverType(4);
  model.passInSolverCharacteristics(&oddCuts);
  // Say no to all solutions by strong branching (c)
  CbcFeasibilityNoStrong noStrong;
  model.setProblemFeasibility(noStrong);
  // Say don't recompute solution d)
  model.setSpecialOptions(4);

  double time1 = CoinCpuTime();

  // Do complete search

  model.branchAndBound();

  std::cout << argv[1] << " took " << CoinCpuTime() - time1 << " seconds, "
            << model.getNodeCount() << " nodes with objective "
            << model.getObjValue()
            << (!model.status() ? " Finished" : " Not finished")
            << std::endl;

  // Print solution if finished - we can't get names from Osi!

  if (!model.status() && model.getMinimizationObjValue() < 1.0e50) {
    int numberColumns = model.solver()->getNumCols();

    //const double * solution = model.bestSolution();
    const double *solution = model.solver()->getColSolution();

    int iColumn;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value = solution[iColumn];
      if (fabs(value) > 1.0e-7 && model.solver()->isInteger(iColumn))
        printf("Column %d has value %g\n", iColumn, value);
    }
  }
  return 0;
}
