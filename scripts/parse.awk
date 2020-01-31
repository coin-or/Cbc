#!/usr/bin/awk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2010            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id$

function abs(x)
{
   return x < 0 ? -x : x;
}
function min(x,y)
{
   return (x) < (y) ? (x) : (y);
}
function max(x,y)
{
   return (x) > (y) ? (x) : (y);
}
BEGIN {
   printf("----------------------------+----------------+----------------+------+-------+-------+--------+---------\n");
   printf("Name                        |   Dual Bound   |  Primal Bound  | Gap%% | Nodes |  Time | Status | Solution \n");
   printf("----------------------------+----------------+----------------+------+-------+-------+--------+---------\n");

   infty = +1e+20;
   eps = 1e-04;
   largegap = 1e+04;

   # initialize summary data
   nsolved = 0;
   nstopped = 0;
   nfailed = 0;

   # initialize data to be set in parse_<solver>.awk
   solver = "?";
   solverversion = "?";
   solverremark = "";
}
# instance name 
/^@01/ { 
   n  = split ($2, a, "/");
   m = split(a[n], b, ".");
   prob = b[1];
   if( b[m] == "gz" || b[m] == "z" || b[m] == "GZ" || b[m] == "Z" )
      m--;
   for( i = 2; i < m; ++i )
      prob = prob "." b[i];

   # initialize data to be set in parse.awk
   timelimit = 0;
   starttime = 0.0;
   endtime = 0.0;
   time = 0.0;

   # initialize data to be set parse_<solver>.awk
   bbnodes = 0;
   pb = +infty;
   db = -infty;
   aborted = 1;
   timeout = 0;
   solstatus = "none";
   read_error = 0;
}
# time
/@03/ { starttime = $2; }
/@04/ { endtime = $2; }
/@05/ { timelimit = $2; }
# solution status
/Read SOL:/ { 
   solstatus = "--";
}
/Check SOL:/ { 
   intcheck = $4;
   conscheck = $6;
   objcheck = $8;
   if( intcheck && conscheck && objcheck ) 
      solstatus = "ok";
   else
      solstatus = "fail";
}
/^=ready=/ {
   # measure wallclock time externaly rounded up to the next second 
   time = max(1, endtime - starttime);

   if( timelimit > 0 )
   {
      # report time limit as time for instances stopped by time limit
      if( timeout )
	 time = timelimit;

      # report time limit as time for aborted instances 
      if( aborted )
	 time = timelimit;

      # report time limit as time for instances that exceeded the time limit but did not stop
      time = min(time, timelimit);
   }

   # determine solving status
   status = "";
   if( aborted && !read_error) 
     status = "abort";
   else if (aborted && read_error)
     status = "noread";
   else if( timeout )
     status = "stopped";
   else
     status = "ok";

   # determine overall status from solving status and solution status:

   # instance solved correctly (including case that no solution was found) 
   if( status == "ok" && (solstatus == "ok" || solstatus == "--") )
      nsolved++;
   # incorrect solving process or infeasible solution (including errors with solution checker)
   else if( status == "abort" || (solstatus == "fail" || solstatus == "error") )
      nfailed++;
   # stopped due to imposed limits
   else if ( status == "stopped" )
      nstopped++;
   else
     nnoread++;
 
  # compute gap
   temp = pb;
   pb = 1.0*temp;
   temp = db;
   db = 1.0*temp;

   if( abs(pb - db) < eps && pb < +infty ) 
      gap = 0.0;
   else if( abs(db) < eps )
      gap = -1.0;
   else if( pb*db < 0.0 )
      gap = -1.0;
   else if( abs(db) >= +infty )
      gap = -1.0;
   else if( abs(pb) >= +infty )
      gap = -1.0;
   else
      gap = 100.0*abs((pb-db)/db);

   if( gap < 0.0 )
      gapstr = "    --";
   else if( gap < largegap )
      gapstr = sprintf("%6.1f", gap);
   else
      gapstr = " Large";

   printf("%-28s %16.9g %16.9g %6s %7d %7d %8s %9s\n", prob, db, pb, gapstr, bbnodes, time, status, solstatus);
}
END {
   printf("----------------------------+----------------+----------------+------+-------+-------+--------+---------\n");
   printf("\n");
   printf("solved/stopped/noread/failed: %d/%d/%d/%d\n", nsolved, nstopped, nnoread, nfailed);
   printf("\n");
   printf("@02 timelimit: %g\n", timelimit);
   printf("@01 %s(%s)%s\n", solver, solverversion, solverremark);
}