#!/usr/bin/awk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2010            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
BEGIN{
   printf("\nSolution: \n\n");
   infeasible = 0;
}
($3 != "objective" && $3 != "gap" && $3 != "time"){
   if (!infeasible){
      printf ("%s %s \n", $2, $3);
   }
}
($3 == "objective" && $1 != "Infeasible"){
      printf ("=obj= %s \n", $5);
}
($3 == "gap"){
      printf ("=obj= %s \n", $8);
}
($3 == "time"){
      printf ("=obj= %s \n", $7);
}   
($1 == "Infeasible"){
   printf ("=infeas= \n");
   infeasible = 1;
}

