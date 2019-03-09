# Messages

Messages and codes passed by CBC are listed in the tables below. For a
complete list, see `Cbc/src/CbcMessages.cpp`. The notation used is the
same as for the `printf` in the C programming language.

  - `%s` is a string
  - `%d` is an integer
  - `%g` or `%f` is a floating point value

There are several log levels. Setting the log level to be `i` produces
the log messages for level `i` and all levels less than `i`.

  - Logging Level 0: Switches off all CBC messages, but one.
  - Logging Level 1: The default.
  - Logging Level 2: Substantial amount of information, e.g., message 15
    is generated once per node. Can be useful when the evaluation at
    each node is slow.
  - Logging Level 3: Tremendous amount of information, e.g., multiple
    messages per node.

CBC Messages Passed At Logging Level 0:

| Code |  | Text and notes                         |
| ---- |  | -------------------------------------- |
| 3007 |  | `No integer variables - nothing to do` |

CBC Messages Passed At or Above Logging Level 1:

| Code |  | Text and notes                                                                                           |
| ---- |  | ---------------------------------------------------------------------------------------------------------|
| 1    |  | `Search completed - best objective %g, took %d iterations and %d nodes`                                  |
| 3    |  | `Exiting on maximum nodes`                                                                               |
| 4    |  | `Integer solution of %g found after %d iterations and %d nodes`                                          |
| 5    |  | `Partial search - best objective %g (best possible %g), took %d iterations and %d nodes`                 |
| 6    |  | `The LP relaxation is infeasible or too expensive`                                                       |
| 9    |  | `Objective coefficients multiple of %g`                                                                  |
| 10   |  | `After %d nodes, %d on tree, %g best solution, best possible %g`                                         |
| 11   |  | `Exiting as integer gap of %g less than %g or %g%%`                                                      |
| 12   |  | `Integer solution of %g found by heuristic after %d iterations and %d nodes`                             |
| 13   |  | `At root node, %d cuts changed objective from %g to %g in %d passes`                                     |
| 14   |  | `Cut generator %d (%s) - %d row cuts (%d active), %d column cuts %? in %g seconds - new frequency is %d` |
| 16   |  | `Integer solution of %g found by strong branching after %d iterations and %d nodes`                      |
| 17   |  | `%d solved, %d variables fixed, %d tightened`                                                            |
| 18   |  | `After tightenVubs, %d variables fixed, %d tightened`                                                    |
| 19   |  | `Exiting on maximum solutions`                                                                           |
| 20   |  | `Exiting on maximum time`                                                                                |
| 23   |  | `Cutoff set to %g - equivalent to best solution of %g`                                                   |
| 24   |  | `Integer solution of %g found by subtree after %d iterations and %d nodes`                               |
| 26   |  | `Setting priorities for objects %d to %d inclusive (out of %d)`                                          |
| 3008 |  | `Strong branching is fixing too many variables, too expensively!`                                        |

CBC Messages Passed At or Above Logging Level 2:

| Code |  | Text and notes                                                   |
| ---- |  | -----------------------------------------------------------------|
| 15   |  | `Node %d Obj %g Unsat %d depth %d`                               |
| 21   |  | `On closer inspection node is infeasible`                        |
| 22   |  | `On closer inspection objective value of %g above cutoff of %g`  |
| 23   |  | `Allowing solution, even though largest row infeasibility is %g` |

CBC Messages Passed At or Above Logging Level 3:

| Code |  | Text and notes                                                  |
| ---- |  | ----------------------------------------------------------------|
| 7    |  | `Strong branching on %d (%d), down %g (%d) up %g (%d) value %g` |
| 25   |  | `%d cleanup iterations before strong branching`                 |
