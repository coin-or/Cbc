# Advanced Solver Uses

## Creating a Solver via Inheritance

CBC uses a generic `OsiSolverInterface` and its `resolve` capability.
This does not give much flexibility so advanced users can inherit from
their interface of choice. This section illustrates how to implement
such a solver for a long thin problem, e.g., fast0507 again. As with the
other examples in the Guide, the sample code is not guaranteed to be the
fastest way to solve the problem. The main purpose of the example is to
illustrate techniques. The full source is in `CbcSolver2.hpp` and
`CbcSolver2.cpp` located in the CBC Samples directory.

The method `initialSolve` is called a few times in CBC, and provides a
convenient starting point. The `modelPtr_` derives from
`OsiClpSolverInterface`.

```
// modelPtr_ is of type ClpSimplex *
modelPtr_->setLogLevel(1); // switch on a bit of printout
modelPtr_->scaling(0); // We don't want scaling for fast0507
setBasis(basis_,modelPtr_); // Put basis into ClpSimplex
// Do long thin by sprint
ClpSolve options;
options.setSolveType(ClpSolve::usePrimalorSprint);
options.setPresolveType(ClpSolve::presolveOff);
options.setSpecialOption(1,3,15); // Do 15 sprint iterations
modelPtr_->initialSolve(options); // solve problem
basis_ = getBasis(modelPtr_); // save basis
modelPtr_->setLogLevel(0); // switch off printout
```

The `resolve()` method is more complicated than `initialSolve()`. The
main pieces of data are a counter `count_` which is incremented each
solve and an integer array `node_` which stores the last time a variable
was active in a solution. For the first few times, the normal Dual
Simplex is called and `node_` array is updated.

```
if (count_<10) {
  OsiClpSolverInterface::resolve(); // Normal resolve
  if (modelPtr_->status()==0) {
    count_++; // feasible - save any nonzero or basic
    const double * solution = modelPtr_->primalColumnSolution();
    for (int i=0;i<numberColumns;i++) {
      if (solution[i]>1.0e-6||modelPtr_->getStatus(i)==ClpSimplex::basic) {
        node_[i]=CoinMax(count_,node_[i]);
        howMany_[i]++;
      }
    }
  } else {
    printf("infeasible early on\n");
  }
}
```

After the first few solves, only those variables which took part in a
solution in the last so many solves are used. As fast0507 is a set
covering problem, any rows which are already covered can be taken out.
```
int * whichRow = new int[numberRows]; // Array to say which rows used
int * whichColumn = new int [numberColumns]; // Array to say which columns used
int i;
const double * lower = modelPtr_->columnLower();
const double * upper = modelPtr_->columnUpper();
setBasis(basis_,modelPtr_); // Set basis
int nNewCol=0; // Number of columns in small model
// Column copy of matrix
const double * element = modelPtr_->matrix()->getElements();
const int * row = modelPtr_->matrix()->getIndices();
const CoinBigIndex * columnStart = modelPtr_->matrix()->getVectorStarts();
const int * columnLength = modelPtr_->matrix()->getVectorLengths();

int * rowActivity = new int[numberRows]; // Number of columns with entries in each row
memset(rowActivity,0,numberRows*sizeof(int));
int * rowActivity2 = new int[numberRows]; // Lower bound on row activity for each row
memset(rowActivity2,0,numberRows*sizeof(int));
char * mark = (char *) modelPtr_->dualColumnSolution(); // Get some space to mark columns
memset(mark,0,numberColumns);
for (i=0;i<numberColumns;i++) {
  bool choose = (node_[i]>count_-memory_&&node_[i]>0); // Choose if used recently
  // Take if used recently or active in some sense
  if ((choose&&upper[i])||(modelPtr_->getStatus(i)!=ClpSimplex::atLowerBound&&
       modelPtr_->getStatus(i)!=ClpSimplex::isFixed)||lower[i]>0.0) {
    mark[i]=1; // mark as used
    whichColumn[nNewCol++]=i; // add to list
    CoinBigIndex j;
    double value = upper[i];
    if (value) {
      for (j=columnStart[i];
           j<columnStart[i]+columnLength[i];j++) {
        int iRow=row[j];
        assert (element[j]==1.0);
        rowActivity[iRow] ++; // This variable can cover this row
      }
      if (lower[i]>0.0) {
        for (j=columnStart[i];
             j<columnStart[i]+columnLength[i];j++) {
          int iRow=row[j];
          rowActivity2[iRow] ++; // This row redundant
        }
      }
    }
  }
}
int nOK=0; // Use to count rows which can be covered
int nNewRow=0; // Use to make list of rows needed
for (i=0;i<numberRows;i++) {
  if (rowActivity[i])
    nOK++;
  if (!rowActivity2[i])
    whichRow[nNewRow++]=i; // not satisfied
  else
    modelPtr_->setRowStatus(i,ClpSimplex::basic); // make slack basic
}
if (nOK<numberRows) {
  // The variables we have do not cover rows - see if we can find any that do
  for (i=0;i<numberColumns;i++) {
    if (!mark[i]&&upper[i]) {
      CoinBigIndex j;
      int good=0;
      for (j=columnStart[i];
           j<columnStart[i]+columnLength[i];j++) {
        int iRow=row[j];
        if (!rowActivity[iRow]) {
          rowActivity[iRow] ++;
          good++;
        }
      }
      if (good) {
        nOK+=good; // This covers - put in list
        whichColumn[nNewCol++]=i;
      }
    }
  }
}
delete [] rowActivity;
delete [] rowActivity2;
if (nOK<numberRows) {
  // By inspection the problem is infeasible - no need to solve
  modelPtr_->setProblemStatus(1);
  delete [] whichRow;
  delete [] whichColumn;
  printf("infeasible by inspection\n");
  return;
}
// Now make up a small model with the right rows and columns
ClpSimplex *  temp = new ClpSimplex(modelPtr_,nNewRow,whichRow,nNewCol,whichColumn);
```

If the variables cover the rows, then the problem is feasible (no cuts
are being used). If the rows were equality constraints, then this might
not be the case. More work would be needed. After the solution, the
reduct costs are checked. If any reduced costs are negative, the code
goes back to the full problem and cleans up with Primal Simplex.

```
temp->setDualObjectiveLimit(1.0e50); // Switch off dual cutoff as problem is restricted
temp->dual(); // solve
double * solution = modelPtr_->primalColumnSolution(); // put back solution
const double * solution2 = temp->primalColumnSolution();
memset(solution,0,numberColumns*sizeof(double));
for (i=0;i<nNewCol;i++) {
  int iColumn = whichColumn[i];
  solution[iColumn]=solution2[i];
  modelPtr_->setStatus(iColumn,temp->getStatus(i));
}
double * rowSolution = modelPtr_->primalRowSolution();
const double * rowSolution2 = temp->primalRowSolution();
double * dual = modelPtr_->dualRowSolution();
const double * dual2 = temp->dualRowSolution();
memset(dual,0,numberRows*sizeof(double));
for (i=0;i<nNewRow;i++) {
  int iRow=whichRow[i];
  modelPtr_->setRowStatus(iRow,temp->getRowStatus(i));
  rowSolution[iRow]=rowSolution2[i];
  dual[iRow]=dual2[i];
}
// See if optimal
double * dj = modelPtr_->dualColumnSolution();
// get reduced cost for large problem
// this assumes minimization
memcpy(dj,modelPtr_->objective(),numberColumns*sizeof(double));
modelPtr_->transposeTimes(-1.0,dual,dj);
modelPtr_->setObjectiveValue(temp->objectiveValue());
modelPtr_->setProblemStatus(0);
int nBad=0;

for (i=0;i<numberColumns;i++) {
  if (modelPtr_->getStatus(i)==ClpSimplex::atLowerBound
      &&upper[i]>lower[i]&&dj[i]<-1.0e-5)
    nBad++;
}
// If necessary clean up with primal (and save some statistics)
if (nBad) {
  timesBad_++;
  modelPtr_->primal(1);
  iterationsBad_ += modelPtr_->numberIterations();
}
```

The array `node_` is updated, as for the first few solves. To give some
idea of the effect of this tactic, the problem fast0507 has 63,009
variables but the small problem never has more than 4,000 variables. In
only about ten percent of solves was it necessary to resolve, and then
the average number of iterations on full problem was less than 20.

## Quadratic MIP

To give another example - again only for illustrative purposes -- it is
possible to do quadratic MIP with CBC. In this case, we make `resolve`
the same as `initialSolve`. The full code is in `ClpQuadInterface.hpp`
and `ClpQuadInterface.cpp` located in the CBC Samples directory.

```
// save cutoff
double cutoff = modelPtr_->dualObjectiveLimit();
modelPtr_->setDualObjectiveLimit(1.0e50);
modelPtr_->scaling(0);
modelPtr_->setLogLevel(0);
// solve with no objective to get feasible solution
setBasis(basis_,modelPtr_);
modelPtr_->dual();
basis_ = getBasis(modelPtr_);
modelPtr_->setDualObjectiveLimit(cutoff);
if (modelPtr_->problemStatus())
  return; // problem was infeasible
// Now pass in quadratic objective
ClpObjective * saveObjective  = modelPtr_->objectiveAsObject();
modelPtr_->setObjectivePointer(quadraticObjective_);
modelPtr_->primal();
modelPtr_->setDualObjectiveLimit(cutoff);
if (modelPtr_->objectiveValue()>cutoff)
  modelPtr_->setProblemStatus(1);
modelPtr_->setObjectivePointer(saveObjective);
```
