# Introduction

The COIN-OR Branch and Cut solver (CBC) is an open-source
mixed-integer program (MIP) solver written in C++. CBC is intended to be
used primarily as a callable library to create customized branch-and-cut
solvers. A basic, stand-alone executable version is also available. CBC
is an active open-source project led by John Forrest at www.coin-or.org.

## Prerequisites

The primary users of CBC are expected to be developers implementing
customized branch-and-cut algorithms in C++ using CBC as a library.
Consequently, this document assumes a working knowledge of
[C++](http://www.cplusplus.com/doc/tutorial/), including basic
object-oriented programming terminology, and familiarity with the
fundamental concepts of [linear programming](http://carbon.cudenver.edu/~hgreenbe/courseware/LPshort/intro.html)
and [mixed integer programming](http://carbon.cudenver.edu/~hgreenbe/courseware/MIP/intro.html).

CBC relies other parts of the COIN-OR repository. CBC needs an LP solver
and relies the COIN-OR Open Solver Inteface (OSI) to communicate with the
user's choice of solver. Any LP solver with an OSI interface can be used
with CBC. The LP solver expected to be used most commonly is COIN-OR's
native linear program solver, CLP. For cut generators, CBC relies on the
COIN-OR Cut Generation Library (CGL). Any cut generator written to CGL
standards can be used with CBC. Some of the cut generators in CGL rely
on other parts of COIN, e.g., CGL's Gomory cut generator rely on the
factorization functionality of `CoinFactorization`. This document
assumes basic familiarity with OSI and CGL.

Technically speaking, CBC assesses the solver (and sometime the model
and data it contains) through an `OSISolverInterface`. For the sake of
simplicity, we will refer to the `OsiSolverInterface` as "the solver" in
this document, rather than "the standard application programming
interface to the solver." We hope any confusion caused by blurring this
distinction will be mitigated by the shorter sentences.

In summary, readers should have the following prerequisites:

  - C++ knowledge,
  - LP and MIP fundamentals, and
  - OSI familiarity.

Unless otherwise stated, we will assume the problem being optimized is a
minimization problem. The terms "model" and "problem" are used
synonymously.

## Branch-and-Cut Overview

Before examining CBC in more detail, we tersely describe the basic
branch-and-cut algorithm by way of example, (which should really be
called branch-and-cut-and-bound) and show the major C++ class(es) in CBC
related to each step. The major CBC classes, labeled (A) through (F),
are described in the table below.

- Step 1. (Bound) Given a MIP model to minimize where some variables must
  take on integer values (e.g., 0, 1, or 2), relax the integrality
  requirements (e.g., consider each "integer" variable to be continuous
  with a lower bound of 0.0 and an upper bound of 2.0). Solve the
  resulting linear model with an LP solver to obtain a lower bound on the
  MIP's objective function value. If the optimal LP solution has integer
  values for the MIP's integer variables, we are finished. Any
  MIP-feasible solution provides an upper bound on the objective value.
  The upper bound equals the lower bound; the solution is optimal.

- Step 2. (Branch) Otherwise, there exists an "integer" variable with a
  non-integral value. Choose one non-integral variable (e.g., with value
  1.3) (A)(B) and branch. Create two nodes, one with the branching
  variable having an upper bound of 1.0, and the other with the branching
  variable having a lower bound of 2.0. Add the two nodes to the search
  tree.

  While (search tree is not empty)
  - Step 3. (Choose Node) Pick a node off the tree (C)(D)
  - Step 4. (Re-optimize LP) Create an LP relaxation and solve.
  - Step 5. (Bound) Interrogate the optimal LP solution, and try to prune
    the node by one of the following.
    - LP is infeasible, prune the node.
    - Else, the optimal LP solution value of the node exceeds the current
      upper bound, prune the node.
    - Else, the optimal LP solution of the node does not exceed the
      current upper bound and the solution is feasible to the MIP. Update
      the upper bound, and the best known MIP solution, and prune the node
      by optimality.
  - Step 6. (Branch) If we were unable to prune the node, then branch.
    Choose one non-integral variable to branch on (A)(B). Create two nodes
    and add them to the search tree.

This is the outline of a "branch-and-bound" algorithm. If in optimizing
the linear programs, we use cuts to tighten the LP relaxations (E)(F),
then we have a "branch-and-cut" algorithm. (Note, if cuts are only used
in Step 1, the method is called a "cut-and-branch"
algorithm.)

| Note | Class name         | Description                                                                                                                                                                                                                                                             |
| ---- | ------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| (A)  | `CbcBranch...`     | These classes define the nature of MIP's discontinuity. The simplest discontinuity is a variable which must take an integral value. Other types of discontinuities exist, e.g., lot-sizing variables.                                                                   |
| (B)  | `CbcNode`          | This class decides which variable/entity to branch on next. Even advanced users will probably only interact with this class by setting `CbcModel` parameters ( e.g., priorities).                                                                                       |
| (C)  | `CbcTree`          | All unsolved models can be thought of as being nodes on a tree where each node (model) can branch two or more times. The user should not need to be concerned with this class.                                                                                          |
| (D)  | `CbcCompare...`    | These classes are used in determine which of the unexplored nodes in the tree to consider next. These classes are very small simple classes that can be tailored to suit the problem.                                                                                   |
| (E)  | `CglCutGenerators` | Any cut generator from CGL can be used in CBC. The cut generators are passed to CBC with parameters which modify when each generator will be tried. All cut generators should be tried to determine which are effective. Few users will write their own cut generators. |
| (F)  | `CbcHeuristics`    | Heuristics are very important for obtaining valid solutions quickly. Some heuristics are available, but this is an area where it is useful and interesting to write specialized ones.                                                                                   |

There are a number of resources available to help new CBC users get
started. This document is designed to be used in conjunction with the
files in the Samples subdirectory of the main CBC directory
(`Cbc/examples`). The Samples illustrate how to use CBC and may also
serve as useful starting points for user projects. In the event that
either this document or the available [Doxygen content](http://www.coin-or.org/Doxygen/Cbc)
conflicts with the observed behavior of the source code, the comments in
the Cbc header files are the ultimate reference.
