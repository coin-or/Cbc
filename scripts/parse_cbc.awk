#!/usr/bin/awk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2010            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id$

# set all solver specific data:
#  solver ["?"]
#  solverversion ["?"]
#  lps ["none"]
#  lpsversion ["-"]
#  bbnodes [0]
#  db [-infty]
#  pb [+infty]
#  aborted [1]
#  timeout [0]

# The solver name
BEGIN {
   solver = "CBC";
   gap = 0;
}

# The solver version 
/^Version:/ { 
   version = $2; 
}

/^Revision Number:/ {
   revision = $3;
   solverversion = version "-" revision
}

# The results
/^Result/ {
   if ($3 == "Optimal"){
      if ($7 == "gap"){
	 gap = 1;
      }else{
	 gap = 0
      }
      aborted = 0;
      timeout = 0;
   }
   if ($5 == "infeasible"){
      pb = +infty;
      db = +infty;
      aborted = 0;
      timeout = 0;
   }else if ($5 == "unbounded"){
      pb = -infty;
      db = -infty;
      aborted = 0;
      timeout = 0;
   }else if ($3 == "Stopped"){
      if ($5 == "time"){
	 timeout = 1;
	 aborted = 0;
      }
   }else if ($3 == "Difficulties"){
      aborted = 1
   }
}

/^Objective value:/ {
   pb = $3;
   if (!gap){
      db = pb;
   }
}

/^Lower bound:/ {
   db = $3;
}

/^Enumerated nodes:/ {
   bbnodes = $3
}

/errors on input/ {
  read_error = 1;
}