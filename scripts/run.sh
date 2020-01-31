#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2010            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id$

SHELL=$1
BINNAME=$2
TSTNAME=$3
TIMELIMIT=$4
HARDMEMLIMIT=$5
THREADS=$6

# construct paths
MIPLIBPATH=`pwd`
BINPATH=$MIPLIBPATH/bin
CHECKERPATH=$MIPLIBPATH/checker
RESULTSPATH=$MIPLIBPATH/results
SCRIPTPATH=$MIPLIBPATH/scripts
TSTPATH=$MIPLIBPATH/testset

# check if the solver link (binary) exists
if test ! -e $BINPATH/$BINNAME
then
    echo "ERROR: solver link <$BINNAME> does not exist in <bin> folder; see bin/README"
    exit;
fi

# check if the test set file/link exists
if test ! -e $TSTPATH/$TSTNAME.test
then
    echo "ERROR: test set file/link <$TSTNAME.test> does not exist in <testset> folder"
    exit;
fi

# grep solver name 
SOLVER=`echo $BINNAME | sed 's/\([a-zA-Z0-9_-]*\).*/\1/g'`

# check if the result folder exist. if not create the result folder
if test ! -e $RESULTSPATH
then
    mkdir $RESULTSPATH
fi

# construct name of output, results, and temporary solution file  
BASENAME=$RESULTSPATH/$TSTNAME.$BINNAME
OUTFILE=$BASENAME.out
RESFILE=$BASENAME.res
SOLFILE=$BASENAME.sol

# absolut tolerance for checking linear constraints and objective value
LINTOL=1e-4 
# absolut tolerance for checking integrality constraints 
INTTOL=1e-4 

# Note that the MIP gap (gap between primal and dual solution) is not
# uniqly defined through all solvers. For example, there is a difference
# between SCIP and CPLEX. All solver, however, have the some behaviour in
# case of a MIP gap of 0.0. 
MIPGAP=0.0

# post system information and current time into the output file
uname -a > $OUTFILE
date >> $OUTFILE

# convert hard memory limit to kilo bytes and post it into the output file
HARDMEMLIMIT=`expr $HARDMEMLIMIT \* 1024`
echo "hard mem limit: $HARDMEMLIMIT k" >> $OUTFILE

# loop over all instance names which are listed in the test set file name
for i in `cat $TSTPATH/$TSTNAME.test` 
do 
    # check if the current instance exists 
    if test -f $i
    then
        echo @01 $i ===========     
        echo -----------------------------
        date
        echo -----------------------------
        TIMESTART=`date +"%s"`
	echo @03 $TIMESTART
	$SHELL -c " ulimit -v $HARDMEMLIMIT k; ulimit -f 2000000; $SCRIPTPATH/run_$SOLVER.sh $SOLVER $BINPATH/$BINNAME $i $TIMELIMIT $SOLFILE $THREADS $MIPGAP"
	echo 
        TIMEEND=`date +"%s"`
	echo @04 $TIMEEND
	echo @05 $TIMELIMIT
	# check if a solution file was written
	if test -e $SOLFILE
	then
	    # check if the link to the solution checker exists
	    if test -f "$CHECKERPATH/bin/solchecker" 
	    then
	    	echo 
	    	$SHELL -c " $CHECKERPATH/bin/solchecker $i $SOLFILE $LINTOL $INTTOL"  
	    	echo
	    else
		echo WARNING: solution cannot be checked because solution checker is missing 
	    fi 
	fi
        echo -----------------------------
        date
        echo -----------------------------
	echo
        echo =ready=
    else
        echo @02 FILE NOT FOUND: $i ===========
    fi
done 2>&1 | tee -a $OUTFILE

date >> $OUTFILE

if test -e $SOLFILE
then
    rm $SOLFILE
fi

awk -f $SCRIPTPATH/parse.awk -f  $SCRIPTPATH/parse_$SOLVER.awk $OUTFILE | tee $RESFILE
