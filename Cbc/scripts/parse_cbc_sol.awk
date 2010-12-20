#!/usr/bin/awk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2010            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
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
}

