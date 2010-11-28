#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2010            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

SOLVER=$1
LPS=$2
NAME=$3
TIMELIMIT=$4
SOLFILE=$5
MIPGAP=$6

$SOLVER -import $NAME -sec $TIMELIMIT -ratio $MIPGAP -solve -solution $SOLFILE

if test -f $SOLFILE
then
    # translate SCIP solution format into format for solution checker.
    #  The SOLFILE format is a very simple format where in each line 
    #  we have a <variable, value> pair, separated by spaces. 
    #  A variable name of =obj= is used to store the objective value 
    #  of the solution, as computed by the solver. A variable name of 
    #  =infeas= can be used to indicate that an instance is infeasible.
    awk -f parse_cbc_sol.awk $SOLFILE
fi


