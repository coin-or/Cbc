# Experiment

Let's execute experiments on the set of instances in ~/inst/super .

Our plan is to evaluate the time taken by cbc with the default LP solver (-initialSolve) versus the version which tries to "guess" (-guess) the best parameters to solve the initial LP relaxation (i.e. primal simplex/dual simplex and other parameters). Our goal is to evaluate the initial LP relaxation *time* and *correctness*. To check correctnes, you can compare the LP optimal objective value produced by CBC with the optimal LP objective values for LP relaxation stored in ~/inst/super (there is a tsv file there I think).

Check in the CBC code if it is prepared to run with initialSolve and with guess.

As most of the instances consume less that 1GB or memory I think we can run 100 instances at time in this machine which has 128 cores and 256 GB of RAM. Let's use gnu parallel to run all the experiments. 

We can use a 2 hour time limit for solving the LP relaxation. Let's kill the execution if it does not respects the time limit at 2 hours and 30 minutes. 

At the end, we want to generate a tsv file with the results:

instance,method,time,killed,obj,matches

where:

instance: instance name
method: initialSolve/guess
time: time execution took
killed: if time limit was not respected and we had to kill the execution
obj: optimal LP relaxation objective value reported by cbc
matches: if optimal lp matches the expected saved LP result

one sample command line to call cbc would be:

cbc instance.mps.gz -seco 7200 -initialSolve -solu solution_instance.sol -quit

please note that time limit must be informed before the command to solve (initialSolve/guess) as these are actions. -solu allows to save the solution for verification purposes. -quit says CBC to not continue the search (no branch and bound).


