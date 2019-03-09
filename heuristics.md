## Getting Good Bounds in CBC

# CbcHeuristic - Heuristic Methods

In practice, it is very useful to get a good solution reasonably fast.
Any MIP-feasible solution produces an upper bound, and a good bound will
greatly reduce the run time. Good solutions can satisfy the user on very
large problems where a complete search is impossible. Obviously,
heuristics are problem dependent, although some do have more general
use. At present there is only one heuristic in CBC itself,
`CbcRounding`. Hopefully, the number will grow. Other heuristics are in
the `Cbc/examples` directory. A heuristic tries to obtain a solution
to the original problem so it only needs to consider the original rows
and does not have to use the current bounds. CBC provides an abstract
base class `CbcHeuristic` and a rounding heuristic in CBC.

This chapter describes how to build a greedy heuristic for a set
covering problem, e.g., the miplib problem fast0507. A more general (and
efficient) version of the heuristic is in `CbcHeuristicGreedy.hpp` and
`CbcHeuristicGreedy.cpp` located in the `Cbc/examples` directory.

The greedy heuristic will leave all variables taking value one at this
node of the tree at value one, and will initially set all other variable
to value zero. All variables are then sorted in order of their cost
divided by the number of entries in rows which are not yet covered. (We
may randomize that value a bit so that ties will be broken in different
ways on different runs of the heuristic.) The best one is choosen, and
set to one. The process is repeated. Because this is a set covering
problem (i.e., all constraints are GE), the heuristic is guaranteed to
find a solution (but not necessarily an improved solution). The speed of
the heuristic could be improved by just redoing those affected, but for
illustrative purposes we will keep it simple.(The speed could also be
improved if all elements are 1.0).

The key `CbcHeuristic` method is `int solution(double & solutionValue,
double * betterSolution)`. The `solution()` method returns 0 if no
solution found, and returns 1 if a solution is found, in which case it
fills in the objective value and primal solution. The code in
`CbcHeuristicGreedy.cpp` is a little more complicated than this
following example. For instance, the code here assumes all variables are
integer. The important bit of data is a copy of the matrix (stored by
column) before any cuts have been made. The data used are bounds,
objective and the matrix plus two work arrays.

```
OsiSolverInterface * solver = model_->solver(); // Get solver from CbcModel
const double * columnLower = solver->getColLower(); // Column Bounds
const double * columnUpper = solver->getColUpper();
const double * rowLower = solver->getRowLower(); // We know we only need lower bounds
const double * solution = solver->getColSolution();
const double * objective = solver->getObjCoefficients(); // In code we also use min/max
double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
double primalTolerance;
solver->getDblParam(OsiPrimalTolerance,primalTolerance);
int numberRows = originalNumberRows_; // This is number of rows when matrix was passed in
// Column copy of matrix (before cuts)
const double * element = matrix_.getElements();
const int * row = matrix_.getIndices();
const CoinBigIndex * columnStart = matrix_.getVectorStarts();
const int * columnLength = matrix_.getVectorLengths();

// Get solution array for heuristic solution
int numberColumns = solver->getNumCols();
double * newSolution = new double [numberColumns];
// And to sum row activities
double * rowActivity = new double[numberRows];
```

The `newSolution` is then initialized to the rounded down solution:
```
for (iColumn=0;iColumn<numberColumns;iColumn++) {
  CoinBigIndex j;
  double value = solution[iColumn];
  // Round down integer
  if (fabs(floor(value+0.5)-value)<integerTolerance)
    value=floor(CoinMax(value+1.0e-3,columnLower[iColumn]));
  // make sure clean
  value = CoinMin(value,columnUpper[iColumn]);
  value = CoinMax(value,columnLower[iColumn]);
  newSolution[iColumn]=value;
  if (value) {
    double cost = direction * objective[iColumn];
    newSolutionValue += value*cost;
    for (j=columnStart[iColumn];j<columnStart[iColumn]+columnLength[iColumn];j++) {
      int iRow=row[j];
      rowActivity[iRow] += value*element[j];
    }
  }
}
```

At this point some row activities may be below their lower bound. To
correct this infeasibility, the variable which is cheapest in reducing
the sum of infeasibilities is found and updated, and the process
repeats. This is a finite process. (Theimplementation could be faster,
but is kept simple for illustrative purposes.)

```
while (true) {
  // Get column with best ratio
  int bestColumn=-1;
  double bestRatio=COIN_DBL_MAX;
  for (int iColumn=0;iColumn<numberColumns;iColumn++) {
    CoinBigIndex j;
    double value = newSolution[iColumn];
    double cost = direction * objective[iColumn];
    // we could use original upper rather than current
    if (value+0.99<columnUpper[iColumn]) {
      double sum=0.0; // Compute how much we will reduce infeasibility by
      for (j=columnStart[iColumn];
           j<columnStart[iColumn]+columnLength[iColumn];j++) {
        int iRow=row[j];
        double gap = rowLower[iRow]-rowActivity[iRow];
        if (gap>1.0e-7)
          sum += CoinMin(element[j],gap);
        if (element[j]+rowActivity[iRow]<rowLower[iRow]+1.0e-7)
          sum += element[j];
      }
      if (sum>0.0) {
        double ratio = (cost/sum)*(1.0+0.1*CoinDrand48());
        if (ratio<bestRatio) {
          bestRatio=ratio;
          bestColumn=iColumn;
        }
      }
    }
  }
  if (bestColumn<0)
    break; // we have finished
  // Increase chosen column
  newSolution[bestColumn] += 1.0;
  double cost = direction * objective[bestColumn];
  newSolutionValue += cost;
  for (CoinBigIndex j=columnStart[bestColumn];
    j<columnStart[bestColumn]+columnLength[bestColumn];j++) {
    int iRow = row[j];
    rowActivity[iRow] += element[j];
  }
}
```

A solution value of `newSolution` is compared to the best solution
value. If `newSolution` is an improvement, its feasibility is validated.
```
returnCode=0; // 0 means no good solution
if (newSolutionValue<solutionValue) { // minimization
  // check feasible
  memset(rowActivity,0,numberRows*sizeof(double));
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    CoinBigIndex j;
    double value = newSolution[iColumn];
    if (value) {
      for (j=columnStart[iColumn];
           j<columnStart[iColumn]+columnLength[iColumn];j++) {
        int iRow=row[j];
        rowActivity[iRow] += value*element[j];
      }
    }
  }
  // check was approximately feasible
  bool feasible=true;
  for (iRow=0;iRow<numberRows;iRow++) {
    if(rowActivity[iRow]<rowLower[iRow])
      if (rowActivity[iRow]<rowLower[iRow]-10.0*primalTolerance)
        feasible = false;
  }
  if (feasible) {
    // new solution
    memcpy(betterSolution,newSolution,numberColumns*sizeof(double));
    solutionValue = newSolutionValue;
    // We have good solution
    returnCode=1;
  }
}
```
