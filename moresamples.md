# CBC's Samples Directory

The CBC distribution includes a number of `.cpp` sample files. Users are
encouraged to use them as starting points for their own CBC projects.
The files can be found in the `Cbc/examples` directory. Most of them can be built by

    make DRIVER=name

which produces an executable `testit`. Below is a list of some of the
most useful sample files with a short description for each file.

Basic Samples:

| Source file | Description                                                                                                                                                                                                                                                                                            |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| MINIMUMCPP  | This is a CBC "Hello, world" program. It reads a problem in MPS file format, and solves the problem using simple branch-and-bound.                                                                                                                                                                     |
| SAMPLE2CPP  | This is designed to be a file that a user could modify to get a useful driver program for his or her project. In particular, it demonstrates the use of CGL's preprocess functionality. It uses `CbcBranchUser.cpp`, `CbcCompareUser.cpp` and `CbcHeuristicUser.cpp` with corresponding `*.hpp` files. |

Advanced Samples:

| Source file | Description                                                                                                                                                                                                                                                                                                                      |
| ----------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| CREWCPP     | This sample shows the use of advanced branching and a use of priorities. It uses `CbcCompareUser.cpp` with corresponding `*.hpp` files.                                                                                                                                                                                          |
| LONGTHINCPP | This sample shows the advanced use of a solver. It also has coding for a greedy heuristic. The solver is given in `CbcSolver2.hpp` and `CbcSolver2.cpp`. The heuristic is given in `CbcHeuristicGreedy.hpp` and `CbcHeuristicGreedy.cpp`. It uses `CbcBranchUser.cpp` and `CbcCompareUser.cpp` with corresponding `*.hpp` files. |
| QMIPCPP     | This solves a quadratic MIP. It is to show advanced use of a solver. The solver is given in `ClpQuadInterface.hpp` and `ClpQuadInterface.cpp`. It uses `CbcBranchUser.cpp` and `CbcCompareUser.cpp` with corresponding `*.hpp` files.                                                                                            |
| SOSCPP      | This artificially creates a Special Ordered Set problem.                                                                                                                                                                                                                                                                         |
| LOTSIZECPP  | This artificially creates a Lot Sizing problem.                                                                                                                                                                                                                                                                                  |
