// $Id$
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>
#include <iomanip>

#include "CoinPragma.hpp"

// For Branch and bound
#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CoinModel.hpp"
// For Linked Ordered Sets
#include "CbcBranchLink.hpp"
#include "OsiClpSolverInterface.hpp"

#include "CoinTime.hpp"

/************************************************************************

This shows how we can define a new branching method to solve problems with
nonlinearities and discontinuities.

We are going to solve the problem

minimize  abs  ( 1.0/6.931 - x1*x4/x2*x3)

where the variables are integral between 12 and 60.
See E.Sangren, "Nonlinear Integer and Discrete Programming in 
Mechanical Design Optimization". Trans. ASME, J. Mech Design 112, 223-229, 1990

One could try to use logarithms to make the problem separable but that leads to a 
weak formulation.  Instaed we are going to use linked
special ordered sets.  The generalization with column generation can be even more powerful
but is not yet in CBC.

The idea is simple:

A linear variable is a convex combination of its lower bound and upper bound!
If x must lie between 12 and 60 then we can substitute for x  as x == 12.0*xl + 60.0*xu where
xl + xu == 1.0.  At first this looks cumbersome but if we have xl12, xl13, ... xl60 and corresponding
xu and yl and yu then we can write:

x == sum 12.0*xl[i] + 60.0* xu[i] where sum xl[i] + xu[i] == 1.0
and 
x*y == 12.0*12.0*xl12 + 12.0*60.0*xu12 + 13.0*12.0*xl13 + 13.0*60.0*x13 ....
                 + 12.0*60*.0xl60 + 60.0*60.0*xu60

And now x*y is correct if x is integer and xl[i], xu[i] are only nonzero for one i.
Note that this would have worked just as easily for y**2 or any clean function of y.

So this is just like a special ordered set of type 1 but on two sets simultaneously.
The idea is even more powerful if we want other functions on y as we can branch on all
sets simultaneously.
Also note that convexity requirements for any non-linear functions are not needed.

So we need a new branching method to do that - see CbcBranchLink.?pp

We are going to need a CbcBranchLink method to see whether we are satisfied etc and also to
create another branching object which knows how to fix variables.  We might be able to use an
existing method for the latter but let us create two methods CbcLink and 
CbcLinkBranchingObject.

For CbcLink we will need the following methods:
Constructot/Destructor
infeasibility - returns 0.0 if feasible otherwise some measure of infeasibility
feasibleRegion - sets bounds to contain current solution
createBranch - creates a CbcLinkBranchingObject

For CbcLinkBranchingObject we need:
Constructor/Destructor
branch - does actual fixing
print - optional for debug purposes.

The easiest way to do this is to cut and paste from CbcBranchActual to get current
SOS stuff and then modify that.

************************************************************************/

int main(int argc, const char *argv[])
{

  OsiClpSolverInterface solver1;

  /*
    We are going to treat x1 and x2 as integer and x3 and x4 as a set.
    We define two new variables y1 == x1*x4 and y2 == x2*x3.
    We define a variable z == x1*x4/x2*x3 so y2*z == y1
    (we will treat y2 as a set)
    Then we have objective - minimize w1 + w2 where
    w1 - w2 = 1.0/6.931 - z

    The model would be a lot smaller if we had column generation.
  */
  // Create model
  CoinModel build;
  // Keep values of all variables for reporting purposes even if not necessary
  /*
    z is first, then x then y1,y2 then w1,w2 
    then y1 stuff, y2 stuff and finally y2 -> z stuff.
    For rows same but 2 per y then rest of z stuff
  */
  int loInt = 12;
  int hiInt = 60;
  int ybaseA = 5, ybaseB = 9, ylen = hiInt - loInt + 1;
  int base = ybaseB + 2 * 2 * ylen;
  int yylen = hiInt * hiInt - loInt * loInt + 1;
  int zbase = 10;
  int i;
  // Do single variables
  double value[] = { 1.0, 1.0 };
  int row[2];
  /* z - obviously we can't choose bounds too tight but we need bounds
     so choose 20% off as obviously feasible.
     fastest way to solve would be too run for a few seconds to get
     tighter bounds then re-formulate and solve.  */
  double loose = 0.2;
  double loZ = (1 - loose) * (1.0 / 6.931), hiZ = (1 + loose) * (1.0 / 6.931);
  row[0] = 0; // for reporting
  row[1] = zbase + 1; // for real use
  build.addColumn(2, row, value, loZ, hiZ, 0.0);
  // x
  for (i = 0; i < 4; i++) {
    row[0] = i + 1;
    build.addColumn(1, row, value, loInt, hiInt, 0.0);
    // we don't need to say x2, x3 integer but won't hurt
    build.setInteger(i + 1);
  }
  // y
  for (i = 0; i < 2; i++) {
    // y from x*x, and convexity
    row[0] = ybaseA + 2 * i;
    if (i == 0)
      row[1] = zbase + 2; // yb*z == ya
    else
      row[1] = zbase - 1; // to feed into z
    build.addColumn(2, row, value, loInt * loInt, hiInt * hiInt, 0.0);
    // we don't need to say integer but won't hurt
    build.setInteger(ybaseA + i);
  }
  // skip z convexity put w in final equation
  row[0] = zbase + 1;
  build.addColumn(1, row, value, 0.0, 1.0, 1.0);
  value[0] = -1.0;
  build.addColumn(1, row, value, 0.0, 1.0, 1.0);
  // Do columns so we know where each is
  for (i = ybaseB; i < base + (2 * yylen); i++)
    build.setColumnBounds(i, 0.0, 1.0);
  // Now do rows
  // z definition
  build.setRowBounds(0, 0.0, 0.0);
  for (i = 0; i < yylen; i++) {
    // l
    build.setElement(0, base + 2 * i, -loZ);
    // u
    build.setElement(0, base + 2 * i + 1, -hiZ);
  }
  // x
  for (i = 0; i < 2; i++) {
    int iVarRow = 1 + i;
    int iSetRow = 4 - i; // as it is x1*x4 and x2*x3
    build.setRowBounds(iVarRow, 0.0, 0.0);
    build.setRowBounds(iSetRow, 0.0, 0.0);
    int j;
    int base2 = ybaseB + 2 * ylen * i;
    for (j = 0; j < ylen; j++) {
      // l
      build.setElement(iVarRow, base2 + 2 * j, -loInt);
      build.setElement(iSetRow, base2 + 2 * j, -loInt - j);
      // u
      build.setElement(iVarRow, base2 + 2 * j + 1, -hiInt);
      build.setElement(iSetRow, base2 + 2 * j + 1, -loInt - j);
    }
  }
  // y
  for (i = 0; i < 2; i++) {
    int iRow = 5 + 2 * i;
    int iConvex = iRow + 1;
    build.setRowBounds(iRow, 0.0, 0.0);
    build.setRowBounds(iConvex, 1.0, 1.0);
    int j;
    int base2 = ybaseB + 2 * ylen * i;
    for (j = 0; j < ylen; j++) {
      // l
      build.setElement(iRow, base2 + 2 * j, -loInt * (j + loInt));
      build.setElement(iConvex, base2 + 2 * j, 1.0);
      // u
      build.setElement(iRow, base2 + 2 * j + 1, -hiInt * (j + loInt));
      build.setElement(iConvex, base2 + 2 * j + 1, 1.0);
    }
  }
  // row that feeds into z and convexity
  build.setRowBounds(zbase - 1, 0.0, 0.0);
  build.setRowBounds(zbase, 1.0, 1.0);
  for (i = 0; i < yylen; i++) {
    // l
    build.setElement(zbase - 1, base + 2 * i, -(i + loInt * loInt));
    build.setElement(zbase, base + 2 * i, 1.0);
    // u
    build.setElement(zbase - 1, base + 2 * i + 1, -(i + loInt * loInt));
    build.setElement(zbase, base + 2 * i + 1, 1.0);
  }
  // and real equation rhs
  build.setRowBounds(zbase + 1, 1.0 / 6.931, 1.0 / 6.931);
  // z*y
  build.setRowBounds(zbase + 2, 0.0, 0.0);
  for (i = 0; i < yylen; i++) {
    // l
    build.setElement(zbase + 2, base + 2 * i, -(i + loInt * loInt) * loZ);
    // u
    build.setElement(zbase + 2, base + 2 * i + 1, -(i + loInt * loInt) * hiZ);
  }
  // And finally two more rows to break symmetry
  build.setRowBounds(zbase + 3, -COIN_DBL_MAX, 0.0);
  build.setElement(zbase + 3, 1, 1.0);
  build.setElement(zbase + 3, 4, -1.0);
  build.setRowBounds(zbase + 4, -COIN_DBL_MAX, 0.0);
  build.setElement(zbase + 4, 2, 1.0);
  build.setElement(zbase + 4, 3, -1.0);
  solver1.loadFromCoinModel(build);
  // To make CbcBranchLink simpler assume that all variables with same i are consecutive

  double time1 = CoinCpuTime();
  solver1.initialSolve();
  solver1.writeMps("bad");
  CbcModel model(solver1);
  model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
  model.solver()->setHintParam(OsiDoScale, false, OsiHintTry);

  CbcObject **objects = new CbcObject *[3];
  /* Format is number in sets, number in each link, first variable in matrix)
      and then a weight for each in set to say where to branch.  
      In this case use NULL to say 0,1,2 ...
      Finally a set number as ID.
  */
  objects[0] = new CbcLink(&model, ylen, 2, ybaseB, NULL, 0);
  objects[0]->setPriority(10);
  objects[1] = new CbcLink(&model, ylen, 2, ybaseB + 2 * ylen, NULL, 0);
  objects[1]->setPriority(20);
  objects[2] = new CbcLink(&model, yylen, 2, base, NULL, 0);
  objects[2]->setPriority(1);
  model.addObjects(3, objects);
  for (i = 0; i < 3; i++)
    delete objects[i];
  delete[] objects;
  model.messageHandler()->setLogLevel(1);
  // Do complete search

  model.setDblParam(CbcModel::CbcMaximumSeconds, 1200.0);
  model.setDblParam(CbcModel::CbcCutoffIncrement, 1.0e-8);
  model.branchAndBound();

  std::cout << "took " << CoinCpuTime() - time1 << " seconds, "
            << model.getNodeCount() << " nodes with objective "
            << model.getObjValue()
            << (!model.status() ? " Finished" : " Not finished")
            << std::endl;

  if (model.getMinimizationObjValue() < 1.0e50) {

    const double *solution = model.bestSolution();
    int numberColumns = model.solver()->getNumCols();
    double x1 = solution[1];
    double x2 = solution[2];
    double x3 = solution[3];
    double x4 = solution[4];
    printf("Optimal solution %g %g %g %g\n", x1, x2, x3, x4);
    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value = solution[iColumn];
      if (fabs(value) > 1.0e-7)
        std::cout << iColumn << " " << value << std::endl;
    }
  }
  return 0;
}
