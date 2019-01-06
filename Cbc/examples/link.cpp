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

minimize  10.0*x + y + z - w

where x, y are continuous between 1 and 10, z can take the values 0, 0.1, 0.2 up to 1.0
and w can take the values 0 or 1.  There is one constraint

w  <= x*z**2 + y*sqrt(z) 

One could try to use logarithms to make the problem separable but that is a very
weak formulation as we want to branch on z directly.  The answer is the concept of linked
special ordered sets.  The generalization with column generation can be even more powerful
but here we are limiting z to discrete values to avoid column generation.

The idea is simple:

A linear variable is a convex combination of its lower bound and upper bound!
If x must lie between 2 and 10 then we can substitute for x  as x == 2.0*xl + 10.0*xu where
xl + xu == 1.0.  At first this looks cumbersome but if we have xl0, xl1, ... xl10 and corresponding
xu and yl and yu then we can write:

x == sum 2.0*xl[i] + 10.0* xu[i] where sum xl[i] + xu[i] == 1.0
and 
x*z**2 == 0.02*xl1 + 0.1*xu1 + 0.08*xl2 + 0.4*xu2 .... + 2.0*xl10 + 10.0*xu10

with similar substitutions for y and y*sqrt(z)

And now the problem is satisfied if w is 0 or 1 and xl[i], xu[i], yl[i] and yu[i] are only
nonzero for one i.

So this is just like a special ordered set of type 1 but on four sets simultaneously.
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

  // Create model
  CoinModel build;
  // Keep x,y and z for reporting purposes in rows 0,1,2
  // Do these
  double value = 1.0;
  int row = -1;
  // x
  row = 0;
  build.addColumn(1, &row, &value, 1.0, 10.0, 10.0);
  // y
  row = 1;
  build.addColumn(1, &row, &value, 1.0, 10.0, 1.0);
  // z
  row = 2;
  build.addColumn(1, &row, &value, 0.0, 1.0, 1.0);
  // w
  row = 3;
  build.addColumn(1, &row, &value, 0.0, 1.0, -1.0);
  build.setInteger(3);
  // Do columns so we know where each is
  int i;
  for (i = 4; i < 4 + 44; i++)
    build.setColumnBounds(i, 0.0, 1.0);
  // Now do rows
  // x
  build.setRowBounds(0, 0.0, 0.0);
  for (i = 0; i < 11; i++) {
    // xl
    build.setElement(0, 4 + 4 * i, -1.0);
    // xu
    build.setElement(0, 4 + 4 * i + 1, -10.0);
  }
  // y
  build.setRowBounds(1, 0.0, 0.0);
  for (i = 0; i < 11; i++) {
    // yl
    build.setElement(1, 4 + 4 * i + 2, -1.0);
    // yu
    build.setElement(1, 4 + 4 * i + 3, -10.0);
  }
  // z - just use x part
  build.setRowBounds(2, 0.0, 0.0);
  for (i = 0; i < 11; i++) {
    // xl
    build.setElement(2, 4 + 4 * i, -0.1 * i);
    // xu
    build.setElement(2, 4 + 4 * i + 1, -0.1 * i);
  }
  // w  <= x*z**2 + y* sqrt(z)
  build.setRowBounds(3, -COIN_DBL_MAX, 0.0);
  for (i = 0; i < 11; i++) {
    double value = 0.1 * i;
    // xl * z**2
    build.setElement(3, 4 + 4 * i, -1.0 * value * value);
    // xu * z**2
    build.setElement(3, 4 + 4 * i + 1, -10.0 * value * value);
    // yl * sqrt(z)
    build.setElement(3, 4 + 4 * i + 2, -1.0 * sqrt(value));
    // yu * sqrt(z)
    build.setElement(3, 4 + 4 * i + 3, -10.0 * sqrt(value));
  }
  // and convexity for x and y
  // x
  build.setRowBounds(4, 1.0, 1.0);
  for (i = 0; i < 11; i++) {
    // xl
    build.setElement(4, 4 + 4 * i, 1.0);
    // xu
    build.setElement(4, 4 + 4 * i + 1, 1.0);
  }
  // y
  build.setRowBounds(5, 1.0, 1.0);
  for (i = 0; i < 11; i++) {
    // yl
    build.setElement(5, 4 + 4 * i + 2, 1.0);
    // yu
    build.setElement(5, 4 + 4 * i + 3, 1.0);
  }
  solver1.loadFromCoinModel(build);
  // To make CbcBranchLink simpler assume that all variables with same i are consecutive

  double time1 = CoinCpuTime();
  solver1.initialSolve();
  solver1.writeMps("bad");
  CbcModel model(solver1);
  model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
  // Although just one set - code as if more
  CbcObject **objects = new CbcObject *[1];
  /* Format is number in sets, number in each link, first variable in matrix)
      and then a weight for each in set to say where to branch.  
      Finally a set number as ID.
  */
  double where[] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
  objects[0] = new CbcLink(&model, 11, 4, 4, where, 0);
  model.addObjects(1, objects);
  delete objects[0];
  delete[] objects;
  // Do complete search

  model.branchAndBound();

  std::cout << "took " << CoinCpuTime() - time1 << " seconds, "
            << model.getNodeCount() << " nodes with objective "
            << model.getObjValue()
            << (!model.status() ? " Finished" : " Not finished")
            << std::endl;

  if (model.getMinimizationObjValue() < 1.0e50) {

    const double *solution = model.bestSolution();
    // check correct
    int which = -1;
    for (int i = 0; i < 11; i++) {
      for (int j = 4 + 4 * i; j < 4 + 4 * i + 4; j++) {
        double value = solution[j];
        if (fabs(value) > 1.0e-7) {
          if (which == -1)
            which = i;
          else
            assert(which == i);
        }
      }
    }
    double x = solution[0];
    double y = solution[1];
    double z = solution[2];
    double w = solution[3];
    // check z
    assert(fabs(z - 0.1 * ((double)which)) < 1.0e-7);
    printf("Optimal solution when x is %g, y %g, z %g and w %g\n",
      x, y, z, w);
    printf("solution should be %g\n", 10.0 * x + y + z - w);
  }
  return 0;
}
