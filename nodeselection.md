# Selecting the Next Node in the Search Tree

## CbcCompare - Comparison Methods

The order in which the nodes of the search tree are explored can
strongly influence the performance of branch-and-cut algorithms. CBC
give users complete control over the search order. The search order is
controlled via the `CbcCompare...` class. CBC provides an abstract base
class, `CbcCompareBase`, and several commonly used instances:

| Class name            | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| --------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `CbcCompareDepth`     | This will always choose the node deepest in tree. It gives minimum tree size but may take a long time to find the best solution.                                                                                                                                                                                                                                                                                                                                                                                            |
| `CbcCompareObjective` | This will always choose the node with the best objective value. This may give a very large tree. It is likely that the first solution found will be the best and the search should finish soon after the first solution is found.                                                                                                                                                                                                                                                                                           |
| `CbcCompareDefault`   | This is designed to do a mostly depth-first search until a solution has been found. It then use estimates that are designed to give a slightly better solution. If a reasonable number of nodes have been explored (or a reasonable number of solutions found), then this class will adopt a breadth-first search (i.e., making a comparison based strictly on objective function values) unless the tree is very large when it will revert to depth-first search. A better description of `CbcCompareUser` is given below. |
| `CbcCompareEstimate`  | When pseudo costs are invoked, they can be used to guess a solution. This class uses the guessed solution.                                                                                                                                                                                                                                                                                                                                                                                                                  |

It is relatively simple for an experienced user to create new compare
class instances. The code in the example below describes how to
build a new comparison class and the reasoning behind it. The complete
source can be found in `CbcCompareUser.hpp` and `CbcCompareUser.cpp`,
located in the CBC Samples directory. The key
method in `CbcCompare` is `bool test(CbcNode* x, CbcNode* y))` which
returns `true` if node `y` is preferred over node `x`. In the `test()`
method, information from `CbcNode` can easily be used.
The following table lists some commonly used methods to access
information at a node.

| Method Name                            | Description                                                                                                                  |
| -------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| `double objectiveValue() const`        | Value of objective at the node.                                                                                              |
| `int numberUnsatisfied() const`        | Number of unsatisfied integers (assuming branching object is an integer - otherwise it might be number of unsatisfied sets). |
| `int depth() const`                    | Depth of the node in the search tree.                                                                                        |
| `double guessedObjectiveValue() const` | If user was setting this (e.g., if using pseudo costs).                                                                      |
| `int way() const`                      | The way which branching would next occur from this node (for more advanced use).                                             |
| `int variable() const`                 | The branching "variable" (associated with the `CbcBranchingObject` -- for more advanced use).                                |

The node desired in the tree is often a function of the how the search
is progressing. In the design of CBC, there is no information on the
state of the tree. CBC is designed so that the method
`newSolution()` is called whenever a solution is found and the method
`every1000Nodes()` is called every 1000 nodes. When these methods are
called, the user has the opportunity to modify the behavior of `test()`
by adjusting their common variables (e.g., `weight_`). Because `CbcNode`
has a pointer to the model, the user can also influence the search
through actions such as changing the maximum time CBC is allowed, once a
solution has been found (e.g., `CbcModel::setMaximumSeconds(double
value)`). In `CbcCompareUser.cpp` of the `Cbc/examples` directory,
four items of data are used.

  1. The number of solutions found so far
  2. The size of the tree (defined to be the number of active nodes)
  3. A weight, `weight_`, which is initialized to -1.0
  4. A saved value of weight, `saveWeight_` (for when weight is set back to -1.0 for special reason)

The full code for the `CbcCompareUser::test()` method is the following:

```
// Returns true if y better than x
bool
CbcCompareUser::test (CbcNode * x, CbcNode * y)
{
  if (weight_==-1.0) {
    // before solution
    if (x->numberUnsatisfied() > y->numberUnsatisfied())
      return true;
    else if (x->numberUnsatisfied() < y->numberUnsatisfied())
      return false;
    else
      return x->depth() < y->depth();
  } else {
    // after solution.
    // note: if weight_=0, comparison is based
    //       solely on objective value
    double weight = CoinMax(weight_,0.0);
    return x->objectiveValue()+ weight*x->numberUnsatisfied() >
      y->objectiveValue() + weight*y->numberUnsatisfied();
  }
}
```

Initially, `weight`\_ is -1.0 and the search is biased towards depth
first. In fact, `test()` prefers `y` if `y` has fewer unsatisfied
variables. In the case of a tie, `test()` prefers the node with the
greater depth in tree. Once a solution is found, `newSolution()` is
called. The method `newSolution()` interacts with `test()` by means of
the variable `weight_`. If the solution was achieved by branching, a
calculation is made to determine the cost per unsatisfied integer
variable to go from the continuous solution to an integer solution. The
variable `weight_` is then set to aim at a slightly better solution.
From then on, `test()` returns `true` if it seems that `y` will lead to
a better solution than `x`. This source for `newSolution()` is the following:

```
// This allows the test() method to change behavior by resetting weight_.
// It is called after each new solution is found.
void
CbcCompareUser::newSolution(CbcModel * model,
                   double objectiveAtContinuous,
                   int numberInfeasibilitiesAtContinuous)
{
  if (model->getSolutionCount()==model->getNumberHeuristicSolutions())
    return; // solution was found by rounding so ignore it.

  // set weight_ to get close to this solution
  double costPerInteger =
    (model->getObjValue()-objectiveAtContinuous)/
    ((double) numberInfeasibilitiesAtContinuous);
  weight_ = 0.98*costPerInteger;
  saveWeight_=weight_;
  numberSolutions_++;
  if (numberSolutions_>5)
    weight_ =0.0; // comparison in test() will be
                  // based strictly on objective value.
}
```

As the search progresses, the comparison can be modified. If many nodes
(or many solutions) have been generated, then `weight_` is set to 0.0
leading to a breadth-first search. Breadth-first search can lead to an
enormous tree. If the tree size is exceeds 10000, it may be desirable to
return to a search biased towards depth first. Changing the behavior in
this manner is done by the method `every1000Nodes` shown next:

```
// This allows the test() method to change behavior every so often
bool
CbcCompareUser::every1000Nodes(CbcModel * model, int numberNodes)
{
  if (numberNodes>10000)
    weight_ =0.0; // compare nodes based on objective value
  else if (numberNodes==1000&&weight_==-2.0)
    weight_=-1.0; // Go to depth first
  // get size of tree
  treeSize_ = model->tree()->size();
  if (treeSize_>10000) {
    // set weight to reduce size most of time
    if (treeSize_>20000)
      weight_=-1.0;
    else if ((numberNodes%4000)!=0)
      weight_=-1.0;
    else
      weight_=saveWeight_;
  }
  return numberNodes==11000; // resort if first time
}
```
