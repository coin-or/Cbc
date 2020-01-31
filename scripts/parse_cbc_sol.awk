#!/usr/bin/awk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2010            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
/^Stopped/ {
   if( NF > 7 )
      exit;
   printf ("=obj= %s \n", $7);
   next;
}
/^Optimal/ {
   printf ("=obj= %s \n", $5);
   next;
}
/^Infeasible/ {
   printf ("=infeas= \n");
   exit;
}
/^Integer/ {
   if( $2 == "infeasible")
      printf ("=infeas= \n");
   exit;
}
//{
   printf ("%s %s \n", $2, $3);
}


