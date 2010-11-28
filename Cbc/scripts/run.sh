#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2010            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: run.sh,v 1.1 2010/10/12 07:43:55 bzfwolte Exp $

SOLVER=$1
LPS=$2
TSTNAME=$3
TIMELIMIT=$4
HARDMEMLIMIT=$5

if test ! -e results
then
    mkdir results
fi

OUTFILE=results/check.$SOLVER.$LPS.$TSTNAME.out
RESFILE=results/check.$SOLVER.$LPS.$TSTNAME.res
SOLFILE=results/check.$SOLVER.$LPS.$TSTNAME.sol

CHECKTOL=-4 # short for 1e-04
MIPGAP=0.0

uname -a > $OUTFILE
date >> $OUTFILE

HARDMEMLIMIT=`expr $HARDMEMLIMIT \* 1024`
echo "hard mem limit: $HARDMEMLIMIT k" >> $OUTFILE

for i in `cat $TSTNAME.test` 
do 
    if test -f $i
    then
        echo @01 $i ===========     
        echo -----------------------------
        date
        echo -----------------------------
        TIMESTART=`date +"%s"`
	echo @03 $TIMESTART
	bash -c " ulimit -v $HARDMEMLIMIT k; ulimit -f 200000; ./run_$SOLVER.sh $SOLVER $LPS $i $TIMELIMIT $SOLFILE $MIPGAP"
        TIMEEND=`date +"%s"`
	echo @04 $TIMEEND
	echo @05 $TIMELIMIT
	if test -f $SOLFILE
	then
	    echo ""
	    # bash -c " ./solchecker $i $SOLFILE $CHECKTOL"  
	    echo ""
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

if test -f $SOLFILE
then
    rm $SOLFILE
fi

awk -f parse.awk -f parse_$SOLVER.awk $OUTFILE | tee $RESFILE
