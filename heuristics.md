# Getting Good Bounds in CBC

## CbcHeuristic - Heuristic Methods

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

## Notes on the Feasibility Pump Implementation

References are
- M. Fischetti, L. Bertacco, and A. Lodi. A feasibility pump heuristic for general mixed-integer problems. Technical report, Universita di Bologna - D.E.I.S. - Operations Research, 2005.
- M. Fischetti, F. Glover, and A. Lodi. The feasibility pump. Mathematical Programming, 2005.

The basic idea (much simplified) is that you start with the relaxed continuous solution and then keep changing the objective function to try and minimize the sum of integer infeasibilities.
If this goes to zero then you have a solution.

So, how long do you try?  If you get a solution, can you get a better one?
Even if you fail to get a solution, can you get any benefit?
The answer to the second question is yes - by adding a constraint that says the objective must be better than this solution.
The answer to the third question is also yes.
By continually changing the objective, you are making the solution go all over the place but there will be variables which despite all this remain fixed at a value.
Maybe if you fix these variables and preprocess the resulting problem you will get a much smaller problem on which it will be worth doing a few nodes of branch and cut.

By default the Feasibility Pump heuristic is run at the start of CBC's branch and cut, but there are several ways to fine tune the heuristic.
Before I answer any of those questions, I should mention one or two undocumented features.
There are several parameters you can set to terminate the search early.
One is `maxNodes`. If you do `maxN??` in Cbc it will say that the valid range is -1 to 2147483647.
If you set maxNodes to 0, cbc will apply all cuts at root node and then terminate, but if you set it to -1, then the cuts will not be computed.
This can save some time if all you want to do is find a heuristic solution.
Also if `allowableGap` or `ratioGap` is set and a heuristic finds a solution which satisfies the gap criterion, then the code will also terminate before doing the node 0 cut computations.
Finally I should mention the `doh` option.  The full name of the parameter is `doHeuristic`.
Normally Cbc does preprocessing and then enters branch and cut, where the first thing to be done is the heuristics phase.
If a valid value for the objective is known then sometimes the preprocessing can do a better job as it can fix variables on reduced costs.
One way of doing this is to run some heuristics, then do preprocessing and go into branch and cut.
To do the feasibility pump then do preprocessing and then do feasibility pump inside branch and cut one would do
```
cbc ....mps -feas both -doh -solve
```

Anyway back to tuning the feasibility pump.
I should mention here that I will give the full name of the option but you can abbreviate it e.g. `passFeasibilityPump` is accepted if you enter `passF`.
There are three options which affect the heuristic:

 * `passfeasibilityPump` - this says how many passes of feasibility pump to do.
   If less than 200 and `hoptions` flag not set then the heuristic will modify this number depending on how things are going.

 * `pumpTune` - this is a very complex parameter as more and more was loaded onto it.
    It is numeric of the form `fedddcba` where
    - `a` is 0,1,2,3.  This is only used if the code is going to fix variables and try a mini branch and cut.
      - 0 and 1 have the same meaning and just says - fix all integer variables which have remained at one bound all the time.
      - 2 says - also fix those general integer variables which have stayed at an interior integral value.
      - 3 says do all that and fix all continuous variables which have stayed at a bound.
      The default is 3.
    - `b` is used to fine tune how first solve is done. Even I have no clue what this does!
    - `c`: at present the only useful values are 0 and 1.
      1 says to add a constraint forcing the objective to be better than some initial value.
      If the parameter `dextra1` is not set (i.e. is still zero) then this value is 5% more than continuous solution value.
      If dextra1 is set then this is used.
    - `ddd` says how many solutions to try and get.
      If 0 then the heuristic will exit after the first solution.
      If n then after n+1 solutions.
      After the first time a constraint is added forcing a better solution.
    - `e` can be 0,1,4,5. 
       If 1 or 5 then the mini branch and cut phase is entered.
       If 4 or 5 then a different method is used to deal with general integer variables.
       The default method is to treat general integer variables as if they were continuous.
       When all 0-1 variables are satisfied then if any general are unsatisfied then all 0-1 variables are fixed and a few nodes of branch and bound are tried to satisfy the general integer variables.
       If 4/5 is set then the original method of Fischetti and Lodi is used which involves adding one constraint and variable for every general integer variable which is being pushed away from bound.
    - `f`: if set then a fraction of the original objective is added in.
      This fraction will decay at each pass.  Larger `f` gives a larger fraction.
    Whew!  So you think you have finished with options - no way.

 * `hOptions` (heuristic options).
   The first part of this applies to all heuristics, but only if one of the gap parameters are set to say end search early if gap is reached.
   If the 1 bit is set i.e. `hOptions` is odd then the heuristic phase is exited immediately the gap is achieved.
   The default behavior is to do all heuristics and then check.
   If the 2 bit is set then do exactly the number of passes specified.
   Normally if number of passes is 30 then the code will do up to 30 passes each major iteration.
   If the 4 bit is set and an initial cutoff is given then the code will relax this cutoff every 50 passes - so that an optimistic guess can be corrected.
   If the 8 bit is set then on the second and subsequent major iterations the cutoff will be changed if it does not look as is the code will find a solution.
