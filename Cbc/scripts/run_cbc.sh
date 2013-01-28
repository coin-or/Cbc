#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2010            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id$

SOLVER=$1
BINNAME=$2
NAME=$3
TIMELIMIT=$4
SOLFILE=$5
THREADS=$6
MIPGAP=$7

TMPFILE=results/check.$SOLVER.tmp

if test $THREADS != 0
then
    $BINNAME -import $NAME -sec $TIMELIMIT -threads $THREADS -ratio $MIPGAP -solve -solution $SOLFILE
else
    $BINNAME -import $NAME -sec $TIMELIMIT -ratio $MIPGAP -solve -solution $SOLFILE
fi

if test -f $SOLFILE
then
    # translate CBC solution format into format for solution checker.
    #  The SOLFILE format is a very simple format where in each line 
    #  we have a <variable, value> pair, separated by spaces. 
    #  A variable name of =obj= is used to store the objective value 
    #  of the solution, as computed by the solver. A variable name of 
    #  =infeas= can be used to indicate that an instance is infeasible.
    awk '
    BEGIN{
	infeasible = 0;
        nointsol = 0;
    }
    ($3 != "objective" && $3 != "gap" && $3 != "time" && $2 != "infeasible"){
	if (!infeasible){
		printf ("%s %s \n", $2, $3);
	    }
    }
    ($3 == "objective" && $1 != "Infeasible" && $2 != "infeasible"){
	printf ("=obj= %s \n", $5);
    }
    ($3 == "gap"){
	printf ("=obj= %s \n", $8);
    }
    ($3 == "time"){
        if ($5 == "integer"){
           printf ("=nointsol= \n");
           nointsol = 1;
        }else{
           printf ("=obj= %s \n", $7);
        }
    }   
    ($1 == "Infeasible" || $2 == "infeasible"){
	printf ("=infeas= \n");
	infeasible = 1;
    }' $SOLFILE | tee $TMPFILE
    mv $TMPFILE $SOLFILE
fi
