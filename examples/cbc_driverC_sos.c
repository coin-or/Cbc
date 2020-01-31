/*
  Copyright (C) 2004-2007 EPRI Corporation and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/

/* This example shows the use of the "C" interface for CBC. */

#include "Cbc_C_Interface.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* prototypes */
void printSolution(Cbc_Model *cbc_model);
Cbc_Model * getDefaultModel(int argc, const char *argv[]);


/* Call back function - just says whenever it gets Clp0005 or Coin0002 */
/* TODO: It seems that Cbc gives callbacks but not Coin */
static void callBack(Cbc_Model * model, int messageNumber,
                     int nDouble, const double * vDouble,
                     int nInt, const int * vInt,
                     int nString, char ** vString) 
{
  const char prefix[] = "cbc_driverC_sos.cpp::callBack(): ";
  const int  VERBOSE = 4;
  if (VERBOSE>0) printf("%s begin\n",prefix);
  if (VERBOSE>1) printf("%s messageNumber %i\n",prefix,messageNumber);

  if (messageNumber==1000002) {
    /* Coin0002 */
    assert (nString==1&&nInt==3);
    printf("Name of problem is %s\n",vString[0]);
    printf("row %d col %d el %d\n",vInt[0],vInt[1],vInt[2]);
  } else if (messageNumber==5) {
    /* Clp0005 */
    int i;
    assert (nInt==4&&nDouble==3); /* they may not all print */
    for (i=0;i<3;i++)
      printf("%d %g\n",vInt[i],vDouble[i]);
  }

  if (VERBOSE>0) printf("%s return\n",prefix);
}

/**
* Get default model 
*/
Cbc_Model * getDefaultModel(int argc, const char *argv[])
{
  const char prefix[] = "cbc_driverC_sos.cpp::getDefaultModel(): ";
  const int  VERBOSE = 4;
  Cbc_Model  *model;
  int status; 

  if (VERBOSE>0) printf("%s begin\n",prefix);
  model = Cbc_newModel(); 

  /** Amount of print out:
  0 - none
  1 - just final
  2 - just factorizations
  3 - as 2 plus a bit more
  4 - verbose
  above that 8,16,32 etc just for selective debug
  */
  Cbc_setLogLevel(model, 1);
  if (VERBOSE>0) printf("%s Log Level %i\n", prefix, Cbc_logLevel(model));
  /* register callback */
  Cbc_registerCallBack(model,callBack);
  /* Keep names when reading an mps file */
  if (argc < 2) {
#if defined(SAMPLEDIR)
  /*
    SAMPLEDIR should be something like "path/to/mps/directory", including the
    quotes and excluding the final directory separator. Don't forget to properly escape
    '\' when using native Windows path syntax.
  */
    status=Cbc_readMps(model, SAMPLEDIR "/p0033.mps") ;
#else
    fprintf(stderr, "Please specify the full path to an MPS file on the command line\n");
    exit(1);
#endif
  }
  else
    status=Cbc_readMps(model,argv[1]);

  if (status) {
    fprintf(stderr,"Bad readMps %s\n",argv[1]);
    fprintf(stdout,"Bad readMps %s\n",argv[1]);
    exit(1);
  }
  Cbc_setOptimizationDirection(model, 1);
  if (1) Cbc_setIntegerTolerance(model,  1.0e-5);

  /* Solve initial LP relaxation */
  Cbc_initialSolve(model);

  if (VERBOSE>0) printf("%s return\n",prefix);
  return model;
} /* getDefaultModel() */

void printSolution(Cbc_Model *cbc_model) {
  
  /*  Now to print out solution.  The methods used return modifiable
      arrays while the alternative names return const pointers -
      which is of course much more virtuous.
     
      This version just does non-zero column and row values.
  */

  /* If we have not kept names (parameter to readMps) this will be 0 
      assert(Cbc_lengthNames(cbc_model));
  */
  {
    int  name_length = 256;
    char model_name[256];
    Cbc_problemName(cbc_model, name_length, model_name);
    printf("Model Name = '%s'\n", model_name);
  }
  printf("Iteration Count = %i\n",Cbc_getIterationCount(cbc_model));
  printf("Iteration Limit = %i\n",Cbc_maximumIterations(cbc_model));
  printf("Is Abandoned = %i\n",Cbc_isAbandoned(cbc_model));
  printf("Is Proven Optimal = %i\n",Cbc_isProvenOptimal(cbc_model));
  printf("Is Proven Infeasible = %i\n",Cbc_isProvenPrimalInfeasible(cbc_model));
  printf("Is Proven Dual Infeasible = %i\n",Cbc_isProvenDualInfeasible(cbc_model));
  printf("Is Proven Unbounded = %i\n",(Cbc_infeasibilityRay(cbc_model) == NULL) ? 0 : 1);
  printf("Is Primal Objective Limit Reached = %i\n",Cbc_isPrimalObjectiveLimitReached(cbc_model));
  printf("Is Dual Objective Limit Reached = %i\n",Cbc_isDualObjectiveLimitReached(cbc_model));
  printf("Is Iteration Limit Reached = %i\n",Cbc_isIterationLimitReached(cbc_model));
  printf("Objective Sense = %g\n",Cbc_getObjSense(cbc_model));  /* (1 - minimize, -1 - maximize, 0 - ignore) */
  printf("Primal Feasible = %i\n",Cbc_primalFeasible(cbc_model));
  printf("Dual Feasible = %i\n",Cbc_dualFeasible(cbc_model));
  printf("Dual Bound = %g\n",Cbc_dualBound(cbc_model));
  printf("Infeasibility Cost = %g\n",Cbc_infeasibilityCost(cbc_model));
  printf("Sum Dual Infeasibilities = %g\n",Cbc_sumDualInfeasibilities(cbc_model));
  printf("Number Dual Infeasibilities = %i\n",Cbc_numberDualInfeasibilities(cbc_model));
  printf("Sum Primal Infeasibilities = %g\n",Cbc_sumPrimalInfeasibilities(cbc_model));
  printf("Number Primal Infeasibilities = %i\n",Cbc_numberPrimalInfeasibilities(cbc_model));
  printf("Objective Value = %g\n",Cbc_objectiveValue(cbc_model)); 
  printf("Optimization Direction = %g\n", Cbc_optimizationDirection(cbc_model));
  printf("  (1 - minimize, -1 - maximize, 0 - ignore)\n");
  printf("--------------------------------------\n");

  /*  Rows */
  {
    int numberRows = Cbc_numberRows(cbc_model);
    int iRow;

    const double * rowPrimal = Cbc_getRowActivity(cbc_model);
    const double * rowDual   = Cbc_getRowPrice(cbc_model);
    const double * rowLower  = Cbc_getRowLower(cbc_model);
    const double * rowUpper  = Cbc_getRowUpper(cbc_model);

    assert(rowPrimal != NULL);
    assert(rowDual   != NULL);
    assert(rowLower  != NULL);
    assert(rowUpper  != NULL);

    printf("                       Primal          Dual         Lower         Upper\n");
    for (iRow=0;iRow<numberRows;iRow++) {
      double value;
      value = rowDual[iRow];
      if (value>1.0e-8||value<-1.0e-8) {
        char name[20];
        sprintf(name," Row%-4i",iRow);
        printf("%6d %8s",iRow,name);
        printf(" %13g",rowPrimal[iRow]);
        printf(" %13g",rowDual[iRow]);
        printf(" %13g",rowLower[iRow]);
        printf(" %13g",rowUpper[iRow]);
        printf("\n");
      }
    }
  }
  printf("--------------------------------------\n");
  /* Columns */
  {
    int numberColumns = Cbc_numberColumns(cbc_model);
    int iColumn;

    const double * columnPrimal    = Cbc_getColSolution(cbc_model);
    const double * columnDual      = Cbc_getReducedCost(cbc_model);
    const double * columnLower     = Cbc_getColLower(cbc_model);
    const double * columnUpper     = Cbc_getColUpper(cbc_model);
    const double * columnObjective = Cbc_getObjCoefficients(cbc_model);

    assert(columnPrimal    != NULL);
    assert(columnDual      != NULL);
    assert(columnLower     != NULL);
    assert(columnUpper     != NULL);
    assert(columnObjective != NULL);

    printf("                       Primal          Dual         Lower         Upper          Cost\n");
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      double value;
      value = columnPrimal[iColumn];
      if (value>1.0e-8||value<-1.0e-8) {
        char name[20];
        sprintf(name," Col%-4i",iColumn);
        /*    	Cbc_columnName(cbc_model,iColumn,name); */
        printf("%6d %8s",iColumn,name);
        printf(" %13g",columnPrimal[iColumn]);
        printf(" %13g",columnDual[iColumn]);
        printf(" %13g",columnLower[iColumn]);
        printf(" %13g",columnUpper[iColumn]);
        printf(" %13g",columnObjective[iColumn]);
        printf("\n");
      }
    }
  }
  printf("--------------------------------------\n");
} /*  printSolution() */

int main (int argc, const char *argv[])
{
  const char prefix[] = "cbc_driverC_sos.cpp:main(): ";
  const int  VERBOSE = 4;
  Cbc_Model * model, * model2;
  double time1;
  char modelName[80];
  int numberIntegers;
  int * integerVariable;

  if (VERBOSE>0) printf("%s begin\n",prefix);
  if (VERBOSE>0) printf("%s Version %g\n",prefix,Cbc_getVersion());

  /* set model using the local routine for reading an MPS file */
  model = getDefaultModel(argc, argv);
  model2 = NULL;  /* used to keep model around */

  /* This clause ought to set the initial basis, but does not yet work. */
  {
    int i;
    int row_length = Cbc_getNumRows(model);
    int col_length = Cbc_getNumCols(model);
    int rim_length = row_length + col_length;
    int elem_length = Cbc_getNumElements(model);

    int * cbc_rowStatus    = NULL;
    int * cbc_columnStatus = NULL;

    if (0) {
      fprintf(stdout,"%s row_length = %i\n", prefix, row_length);
      fprintf(stdout,"%s col_length = %i\n", prefix, col_length);
      fprintf(stdout,"%s rim_length = %i\n", prefix, rim_length);
      fprintf(stdout,"%s elem_length = %i\n", prefix, elem_length);
      fflush(stdout);
    }
    /* print solution status variables */
    if (0) {
      if (cbc_rowStatus) {
        for (i = 0; i < row_length; i++) {
          printf("%s cbc_rowStatus[%i] = %i\n", prefix, i, cbc_rowStatus[i]);
          fflush(stdout);
        }
      } else {
        fprintf(stdout,"%s cbc_rowStatus = %p\n", prefix, (void*)cbc_rowStatus);
        fflush(stdout);
      }
      if (cbc_columnStatus) {
        for (i = 0; i < row_length; i++) {
          fprintf(stdout,"%s cbc_rowStatus[%i] = %i\n", prefix, i, cbc_columnStatus[i]);
          fflush(stdout);
        }
      } else {
        fprintf(stdout,"%s cbc_columnStatus = %p\n", prefix, (void*)cbc_columnStatus);
        fflush(stdout);
      }
    }
  }

  /* Save model as a clone (does not work as of 2004). */
  model2 = Cbc_clone(model);

  /* Store which variables are integer as defined in the MPS file.
     They will be set to Continuous before adding the SOS constraints. */

  {
    int i;
    int numberColumns=Cbc_getNumCols(model);
    numberIntegers = 0;
    integerVariable = malloc(numberColumns*sizeof(int[1]));
    for ( i=0;i<numberColumns;i++) {
      if (Cbc_isInteger(model,i)) {
        integerVariable[numberIntegers++]=i;
        if (VERBOSE>3) printf("%s integerVariable[%i] = %i\n",prefix,
          numberIntegers-1, integerVariable[numberIntegers-1]);
      }
    }
  }
  Cbc_problemName(model,80,modelName);

  /* Solve the MPS version of the problem */
  if (1) {  
    if (VERBOSE>1) {
      printf("%s Solve MPS version of the problem\n",prefix);
      printf("%s Optimization Direction = %g (1 - minimize, -1 - maximize, 0 - ignore)\n",
        prefix, Cbc_optimizationDirection(model)); 
    }
    /*    Cbc_setLogLevel(model, VERBOSE); // 4 is verbose */
    time1 = Cbc_cpuTime(model) ;
    Cbc_branchAndBound(model);
    if (VERBOSE>1) printf("Model %s has %d rows and %d columns\n",
      modelName,Cbc_numberRows(model),Cbc_numberColumns(model));

    if (VERBOSE>1) printf("%s Solving model %s took %g seconds, %i nodes with objective %g\n",
      prefix, modelName, (Cbc_cpuTime(model)-time1), 
      Cbc_getNodeCount(model), Cbc_getObjValue(model));
    if (VERBOSE>0) (!Cbc_status(model)) ? printf(" Finished\n") : printf(" Not finished\n");
    if (VERBOSE>2) printSolution(model);
  }
  {
    
    /* SOS specification data
       specify numberSets, numPoints, and whichRanges explicitly
       NOTE: These need to be commented according to the MPS file.
         Example_1_1: Cbc0004I MPS reads 1081 only columns out of 1085
                      Then the 4th range goes out of bounds
         Example_2: Cbc0006I The LP relaxation is infeasible or too expensive
         Mod_RT_1: Cbc0006I The LP relaxation is infeasible or too expensive
    */
    const int numberSets = 
      1; /* dummy
    //    4; // Example_1_1
    //    2;  // Example_2
    //    2;  // Mod_RT_1
    //    2;  // SITH */
    const int numPoints = 
      1; /* dummy
    //    257; // Example_1_1
    //    256;  // Example_2
    //    257;  // Mod_RT_1
    //    257;  // SITH */

    const int whichRanges[1][2] = { /* counting from zero? */
      {0,0} /* dummy */
      /*   {56,312}, {313, 569}, {572, 828}, {829, 1085} // Example_1_1
           {48, 303}, {304, 559} // Example_2
           {45, 301}, {302, 558} // Mod_RT_1
           {45, 301}, {302, 558} // SITH
      */
    };
    /* the rest is determined parametrically */
    int *len = malloc(numberSets* sizeof(int[1]));
    int **which = malloc(numberSets* sizeof(int[1]));
    int setNum, pointNum;
    double *weights;
    int i, j;
    for (setNum = 0; setNum < numberSets; setNum++) {
      len[setNum] = whichRanges[setNum][1] - whichRanges[setNum][0] + 1;
      if (len[setNum] != numPoints) {
        printf("%s ERROR: len[%i] (%i) != numPoints (%i)\n",
          prefix, setNum, len[setNum], numPoints);
        return 1;
      }
      which[setNum] = malloc(numPoints*sizeof(int[1]));
      for (j = 0; j < len[setNum]; j++)
        which[setNum][j] = whichRanges[setNum][0] + j; /* Example_2 */
    }
    weights = malloc(numPoints*sizeof(double[1]));  
    for (pointNum = 0; pointNum < numPoints; pointNum++) weights[pointNum] = pointNum+1;

    /* Now use SOS2
       NOTE: Only enable this if known good SOS Specification (above)
    */
    if (1) {  
      if (VERBOSE>1) printf("%s Use SOS2\n",prefix);

      /* Restore model */
      Cbc_deleteModel(model);
      if (0) {
        model = Cbc_clone(model2);
      } else {
        model = getDefaultModel(argc, argv); 
      }
      /*    Cbc_setLogLevel(model, 4); // 4 is verbose */
      if (VERBOSE>1) printf("%s Model %s has %d rows and %d columns\n",
        prefix, modelName,Cbc_numberRows(model),Cbc_numberColumns(model));

      /* make SOS cuts */
      for (i=0;i<numberIntegers;i++) {
        int iColumn = integerVariable[i];
        /* Stop being integer */
        Cbc_setContinuous(model, iColumn);
      }
      /* add cut (0 - use dense, 1 - use sparse (TODO: test this)) */
      if (0) for (i=0;i<numberSets;i++) {
        printf(
          "%s calling Cbc_addSOS(), len[%i] = %i, which[%i][0] = %i, which[%i][%i] = %i,\n  weights[0] = %g, weights[%i] = %g\n",
          prefix, i, len[i], i, which[i][0], i, len[i]-1, which[i][len[i]-1], weights[0], len[i]-1, weights[len[i]-1]);
        /*      Cbc_addSOS_Sparse(model,len[i],which[i],weights,i,2); */
      } else {
        int numObjects = numberSets; /* cannot pass const int */
        Cbc_addSOS_Dense(model, numObjects, len, (const int* const *)which, (const double*)weights, 2);
      }
    }

    Cbc_setOptimizationDirection(model, 1); /* 1 minimize */
    if (VERBOSE>1) {
      printf("%s Solve MPS version of the problem\n",prefix);
      printf("%s Optimization Direction = %g (1 - minimize, -1 - maximize, 0 - ignore)\n",
        prefix, Cbc_optimizationDirection(model)); 
    }
    if (VERBOSE>1) printf("%s calling Cbc_scaling()\n",prefix);
    Cbc_scaling(model,1);
    time1 = Cbc_cpuTime(model) ;
    if (VERBOSE>1) printf("%s calling Cbc_initialSolve()\n",prefix);
    Cbc_initialSolve(model);
    if (VERBOSE>3) Cbc_printModel(model,prefix);
    if (VERBOSE>1) printf("%s calling Cbc_branchAndBound()\n",prefix);
    Cbc_branchAndBound(model);
    if (VERBOSE>1) printf("%s %s took %g seconds, %i nodes with objective %g\n",
      prefix, modelName, (Cbc_cpuTime(model)-time1), 
      Cbc_getNodeCount(model), Cbc_getObjValue(model));
    if (VERBOSE>0) (!Cbc_status(model)) ? printf(" Finished\n") : printf(" Not finished\n");
    if (VERBOSE>2) printSolution(model);
  }

  if (VERBOSE>1) printf("%s Log Level %i\n", prefix, Cbc_logLevel(model));
  if (VERBOSE>0) printf("%s return 0\n",prefix);
  return 0;
} /*  main() */
