# The CBC Model Class

## Overview

The main class in CBC is `CbcModel`. The `CbcModel` class is where most
of the parameter setting is done. The absolute minimum number of actions
taken with `CbcModel` is two,

  - `CbcModel(OsiSolverInterface & linearSolver)` as constructor, and

  - `branchAndBound()` for solving the problem.

## Simple Branch-and-Bound Example

The first sample program shows how to perform simple branch-and-bound
with CBC. This program is short enough to present in full. Most of the
remaining examples will take the form of small code fragments. The
complete code for all the examples in this Guide can be found in the CBC
Samples directory, `Cbc/examples`.

```
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CbcModel.hpp"

// Using CLP as the solver
#include "OsiClpSolverInterface.hpp"

int main (int argc, const char *argv[])
{
  OsiClpSolverInterface solver1;

  // Read in example model in MPS file format
  // and assert that it is a clean model
  int numMpsReadErrors = solver1.readMps("../../Data/Sample/p0033.mps","");
  assert(numMpsReadErrors==0);

  // Pass the solver with the problem to be solved to CbcModel
  CbcModel model(solver1);

  // Do complete search
  model.branchAndBound();

  // Print the solution.  CbcModel clones the solver so we
  //  need to get current copy from the CbcModel
  int numberColumns = model.solver()->getNumCols();

  const double * solution = model.solver()->getColSolution();

  for (int iColumn=0;iColumn<numberColumns;iColumn++) {
    double value=solution[iColumn];
    if (fabs(value)>1.0e-7&&model.solver()->isInteger(iColumn))
      printf("%d has value %g\n",iColumn,value);
   }
  return 0;
}
```

This program creates a `OsiClpSolverInterface` solver interface
(i.e., `solver1`), and reads an MPS file. If there are no errors, the
program passes the problem to `CbcModel` which solves the problem using
the branch-and-bound algorithm. The part of the program which solves the
problem is very small (one line!) but before that one line, the LP solver
(i.e., `solver1`) had to be created and populated with the problem.
After that one line, the results were printed out.

## The Relationship Between OSI and CBC

The above program illustrates the dependency
of CBC on the `OsiSolverInterface` class. The constructor of `CbcModel`
takes a pointer to an `OsiSolverInterface` (i.e., a solver). The
`CbcModel` clones the solver, and uses its own instance of the solver.
The `CbcModel`'s solver and the original solver (e.g., `solver1`) are
not in sync unless the user synchronizes them. The user can always
access the `CbcModel`'s solver through the `model()` method. To
synchronize the two solvers, explicitly refreshing the original, e.g.,

```
solver1 = model.solver();
```

`CbcModel`'s method `solve()` returns a pointer to CBC's cloned solver.

For convenience, many of the OSI methods to access problem data have
identical method names in `CbcModel`. (It's just more convenient to type
`model.getNumCols()` rather than `model.solver()->getNumCols()`). The
`CbcModel` refreshes its solver at certain logical points during the
algorithm. At these points, the information from the `CbcModel` `model`
will match the information from the `model.solver()`. Elsewhere, the
information may vary. For instance, the OSI method `getColSolution()`
will contain the best solution so far, while the `CbcModel` method may
not. In this case, it is safer to use `CbcModel::bestSolution()`.

While all the OSI methods have equivalent methods
in `CbcModel`, there are some OSI methods which do not. For example, if
the program produced a lot of undesired output, one might add the line

```
model.solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
```

to reduce the output. There is no `setHintParam()` method in `CbcModel`.

## Getting Solution Information

Optimality can be checked through a call to `model.isProvenOptimal()`.
Also available are `isProvenInfeasible()`, `isSolutionLimitReached()`,
`isNodeLimitReached()` or the feared `isAbandoned()`. There is also
`int status()` which returns 0 if finished (which includes the case
when the algorithm is finished because it has been proved infeasible), 1
if stopped by user, and 2 if difficulties arose.

In addition to these `CbcModel` methods, solution values can be accessed
via OSI methods. The OSI methods pick up the current solution in the
`CBCModel`. The current solution will match the best solution found so
far if called after `branchAndBound()` and a solution was
found.

| Purpose                    | Name                              | Notes                                                                                                                                                      |
| -------------------------- | --------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Primal column solution     | `const double * getColSolution()` | The OSI method will return the best solution found thus far, unless none has been found. It is safer to use `CbcModel` version, `CbcModel::bestSolution()` |
| Dual row solution          | `const double * getRowPrice()`    | Identical `CbcModel` version available, `CbcModel::getRowPrice()`.                                                                                         |
| Primal row solution        | `const double * getRowActivity()` | Identical `CbcModel` version available, `CbcModel::getRowActivity()`.                                                                                      |
| Dual column solution       | `const double * getReducedCost()` | Identical `CbcModel` version available, `CbcModel::gtReducedCost()`.                                                                                       |
| Number of rows in model    | `int getNumRows()`                | Identical `CbcModel` version available, `CbcModel::getNumRows()`. Note: the number of rows can change due to cuts.                                         |
| Number of columns in model | `int getNumCols()`                | Identical `CbcModel` version available, `CbcModel::getNumCols()`.                                                                                          |

## Useful Set and Get Methods in `CbcModel`

Most of the parameter setting in CBC is done through `CbcModel` methods.
The most commonly used set and get methods are the following:

| Method(s)                                                                                                                                                                                                                                             | Description                                                                                                                                                |
| ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `bool setMaximumNodes(int value)` `int getMaximumNodes() const` `bool setMaximumSeconds(double value)` `double getMaximumSeconds()` `bool setMaximumSolutions(double value)` `double getMaximumSolutions() const`                                     | These set methods tell CBC to stop after a given number of nodes, seconds, or solutions is reached. The get methods return the corresponding values.       |
| `bool setIntegerTolerance(double value) const` `double getIntegerTolerance() const`                                                                                                                                                                   | An integer variable is deemed to be at an integral value if it is no further than this `value` (tolerance) away.                                           |
| `bool setAllowableGap(double value)` `double getAllowableGap() const` `bool setAllowablePercentageGap(double value)` `double getAllowablePercentageGap() const` `bool setAllowableFractionGap(double value)` `double getAllowableFractionGap() const` | `CbcModel` returns if the gap between the best known solution and the best possible solution is less than this `value`, or as a percentage, or a fraction. |
| ` void setNumberStrong(double value)  ` ` int numberStrong() const  `                                                                                                                                                                                 | These methods set or get the maximum number of candidates at a node to be evaluated for strong branching.                                                  |
| ` void setPrintFrequency(int value)  ` `int printFrequency() const`                                                                                                                                                                                   | Controls the number of nodes evaluated between status prints. Print frequency has a very slight overhead, if `value` is small.                             |
| `int getNodeCount() const`                                                                                                                                                                                                                            | Returns number of nodes evaluated in the search.                                                                                                           |
| `int numberRowsAtContinuous() const`                                                                                                                                                                                                                  | Returns number of rows at continuous                                                                                                                       |
| `int  numberIntegers() const` `const int * integerVariable() const`                                                                                                                                                                                   | Returns number of integer variables and an array specifying them.                                                                                          |
| `bool isBinary(int colIndex) const` `bool isContinuous(int colIndex) const` `bool isInteger(int colIndex) const`                                                                                                                                      | Returns information on variable `colIndex`. OSI methods can be used to set these attributes (before handing the model to `CbcModel`).                      |
| `double getObjValue() const`                                                                                                                                                                                                                          | This method returns the best objective value so far.                                                                                                       |
| `double getCurrentObjValue() const`                                                                                                                                                                                                                   | This method returns the current objective value.                                                                                                           |
| `const double * getObjCoefficients() const`                                                                                                                                                                                                           | This method return the objective coefficients.                                                                                                             |
| `const double * getRowLower() const` `const double * getRowUpper() const` `const double * getColLower() const` `const double * getColUpper() const`                                                                                                   | These methods return the lower and upper bounds on row and column activities.                                                                              |
| `const CoinPackMatrix * getMatrixByRow() const`                                                                                                                                                                                                       | This method returns a pointer to a row copy of matrix stored as a `CoinPackedMatrix` which can be further examined.                                        |
| `const CoinPackMatrix * getMatrixByCol() const`                                                                                                                                                                                                       | This method returns a pointer to a column copy of matrix stored as a `CoinPackedMatrix` which can be further examined.                                     |
| `CoinBigIndex getNumElements() const` (`CoinBigIndex` is a `typedef` which in most cases is the same as `int`)                                                                                                                                        | Returns the number of nonzero elements in the problem matrix.                                                                                              |
| `void setObjSense(double value)` `double getObjSense() const`                                                                                                                                                                                         | These methods set and get the objective sense. The parameter `value` should be +1 to minimize and -1 to maximize.                                          |

## Impacting the Solution Process

`CbcModel` is extremely flexible and customizable. The class structure
of CBC is designed to make the most commonly desired customizations of
branch-and-cut possible. These include:

  - selecting the next node to consider in the search tree,
  - determining which variable to branch on,
  - using heuristics to generate MIP-feasible solutions quickly,
  - including cut generation when solving the LP-relaxations, and
  - invoking customized subproblem solvers.

To enable this flexibility, `CbcModel` uses other classes in CBC (some
of which are virtual and may have multiple instances). Not all classes
are created equal. The two tables below list in alphabetical order the
classes used by `CbcModel` that are of most interest and of least
interest.

| Class name           | Description                                                                                                                                             | Notes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| -------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `CbcCompareBase`     | Controls which node on the tree is selected.                                                                                                            | The default is `CbcCompareDefault`. Other comparison classes in `CbcCompareActual.hpp` include `CbcCompareDepth` and `CbcCompareObjective`. Experimenting with these classes and creating new compare classes is easy.                                                                                                                                                                                                                                                                              |
| `CbcCutGenerator`    | A wrapper for `CglCutGenerator` with additional data to control when the cut generator is invoked during the tree search.                               | Other than knowing how to add a cut generator to `CbcModel`, there is not much the average user needs to know about this class. However, sophisticated users can implement their own cut generators.                                                                                                                                                                                                                                                                                                |
| `CbcHeuristic`       | Heuristic that attempts to generate valid MIP-solutions leading to good upper bounds.                                                                   | Specialized heuristics can dramatically improve branch-and-cut performance. As many different heuristics as desired can be used in CBC. Advanced users should consider implementing custom heuristics when tackling difficult problems.                                                                                                                                                                                                                                                             |
| `CbcObject`          | Defines what it means for a variable to be satisfied. Used in branching.                                                                                | Virtual class. CBC's concept of branching is based on the idea of an "object". An object has (i) a feasible region, (ii) can be evaluated for infeasibility, (iii) can be branched on, e.g., a method of generating a branching object, which defines an up branch and a down branch, and (iv) allows comparsion of the effect of branching. Instances of objects include `CbcSimpleInteger`, `CbcSimpleIntegerPseudoCosts`, `CbcClique`, `CbcSOS` (type 1 and 2), `CbcFollowOn`, and `CbcLotsize`. |
| `OsiSolverInterface` | Defines the LP solver being used and the LP model. Normally a pointer to the desired `OsiSolverInteface` is passed to `CbcModel` before branch and cut. | Virtual class. The user instantiates the solver interface of their choice, e.g., `OsiClpSolverInterface`.                                                                                                                                                                                                                                                                                                                                                                                           |

There is not much about the following classes that the average user needs to know about.

| Class name           | Description                                                                                                        | Notes                                                                                                                       |
| -------------------- | ------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------- |
| `CbcBranchDecision`  | Used in choosing which variable to branch on, however, most of the work is done by the definitions in `CbcObject`. | Defaults to `CbcBranchDefaultDecision`.                                                                                     |
| `CbcCountRowCut`     | Interface to `OsiRowCut`. It counts the usage so cuts can gracefully vanish.                                       | See `OsiRowCut` for more details.                                                                                           |
| `CbcNode`            | Controls which variable/entity is selected to be branch on.                                                        | Controlled via `CbcModel` parameters. Information from `CbcNode` can be useful in creating customized node selection rules. |
| `CbcNodeInfo`        | Contains data on bounds, basis, etc. for one node of the search tree.                                              | Header is located in `CbcNode.hpp`.                                                                                         |
| `CbcTree`            | Defines how the search tree is stored.                                                                             | This class can be changed but it is not likely to be modified.                                                              |
| `CoinMessageHandler` | Deals with message handling                                                                                        | The user can inherit from `CoinMessageHandler` to specialize message handling.                                              |
| `CoinWarmStartBasis` | Basis representation to be used by solver                                                                          |                                                                                                                             |
