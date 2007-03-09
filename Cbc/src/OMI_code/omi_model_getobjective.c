/**
  * @(#)omi_model_getobjective.c
  * @author Robert Entriken<br>
  * Copyright (C) 1999-2007 EPRI<br>
  * All rights reserved.
  * @version 07-03-07
  * @since OMI_1.0
  *
  * Revisions:
  *   06-02-20 bobe     added copyright
  *   06-01-28 jmm      new package layout
  *   06-04-10 jmm      minor adjustments to compile under Linux
  *   07-03-07 bobe     updated Copyright
  *
  */
#include <time.h>
#include "com_epri_omi_Model.h"
#include "omi_Model.h"
#include "omi_State.h"

/**
  *  Call back to get the objective value computed in Java.
  *
  */
void omi_model_getobjective__(integer *mode, integer *n, doublereal *x,
  doublereal *f, doublereal *g, integer* nstate, integer* nprob) {
  int VERBOSE = 0;
  
/*
  extern time_t time(int);
  extern char * asctime(int);
  extern int localtime(time_t *);
*/

  // local state variables
  JNIEnv         *env;
  jobject        jmodel;
  ClassOMI_Model *model;
  FILE           *out;
  
  // io variables
// unused  FILE *out2, *out_old;
  time_t     t;
  
  // local variables
  jclass    cls;  // holds pointer to omi.Model.class
  jmethodID mid;  // holds offset for omi.Model.getObjective()
  jfieldID  fid;
  int i;
  
  char prefix[] = "Model.c:omi_model_getobjective__():";
  char javaName[] = "getObjective";
//	char mySig[]  = "(II[DD[DII)V";
//  char javaSig[]  = "(IIII)V"; // for the minos version
  char javaSig[]  = "()V"; // for the omi version
  
  jsize start, len;

  // fill in state
  env    = omi_state.env;
  jmodel = omi_state.jmodel;
  model  = omi_state.model;
  out    = omi_state.out;

  /*
  out2 = fopen("minos.txt", "a");
  fprintf(out2,"begin\n");
  fprintf(out2,"%s begin\n",prefix);
  if (OMI_state_1.env    == NULL) fprintf(out2,"%s env is NULL\n",prefix);
  if (OMI_state_1.jmodel == NULL) fprintf(out2,"%s jmodel is NULL\n",prefix);
  if (OMI_state_1.out    == NULL) fprintf(out2,"%s out is NULL\n",prefix);
  fprintf(out2,"%s &env = %X, &jmodel = %X, &out = %X\n",prefix,&env,&jmodel,&out);
  fprintf(out2,"%s  env = %X,  jmodel = %X,  out = %X\n",prefix,env,jmodel,out);
  fprintf(out2,"%s &OMI_state_1.env = %X, &OMI_state_1.jmodel = %X, &OMI_state_1.out = %X\n",
    prefix,&OMI_state_1.env,&OMI_state_1.jmodel,&OMI_state_1.out);
  fclose(out2);
  */
  
  // Print entry time
  t = time(NULL);
  if (VERBOSE>0) fprintf(out,"\n%s started %s\n", prefix, asctime(localtime(&t)));
  
  // Print arguments
  if (VERBOSE>1) {
    fprintf(out, "%s *mode = %i, *n = %i, *nstate = %i, *nprob = %i\n",
      prefix, *mode, *n, *nstate, *nprob);
    for (i = 0; i < *n; i++) fprintf(out, "%s x[%i] = %g\n",prefix, i, x[i]);
  }
  
  // Check state.
  if (env == NULL) {
    fprintf(out, "%s !ERROR, NULL 'env' from OMI_state.\n", prefix);
    fprintf(out,"%s return\n",prefix);
    fclose(out);
    return;

  };  
  if (jmodel == NULL) {
    fprintf(out, "%s !ERROR, NULL 'jmodel' from OMI_state.\n", prefix);
    fprintf(out,"%s return\n",prefix);
    fclose(out);
    return;
  };
  if (VERBOSE>0) fprintf(out, "%s State OK.\n", prefix);

  // get the class pointer and method offset
  cls = (*env)->GetObjectClass(env, jmodel);
  if (VERBOSE>2) fprintf(out, "%s Got Class pointer.\n", prefix);
  
  mid = (*env)->GetMethodID(env, cls, javaName, javaSig);
  if (VERBOSE>2) fprintf(out, "%s Got Method ID.\n", prefix);

  if (mid == NULL) {
    fprintf(out,"%s ERROR: Cannot GetMethodID for '%s', '%s'\n",
      prefix,javaName,javaSig);
    fprintf(out,"%s return\n",prefix);
    fclose(out);
    return;
  }
  
  // put x in jgetObjectiveArg
  fid = (*env)->GetFieldID(env, cls, "getObjectiveArg", "[D");
  if (fid == 0) {
    printf("%s!ERROR, Cannot GetFieldID of 'getObjectiveArg'.\n",prefix);
    return;
  };
  if ((model->jgetObjectiveArg = (*env)->GetObjectField(env,jmodel,fid)) == NULL) {
      printf("%s!ERROR, Cannot GetDoubleField of 'getObjectiveArg'.\n",prefix);
    return;
  };
  start = 0;
  len   = *n;
  (*env)->SetDoubleArrayRegion(env, model->jgetObjectiveArg, start, len, x);

  // call back to Java
  if (VERBOSE>0) fprintf(out, "%s In C, calling omi_model_getobjective__\n", prefix);
  
//  (*env)->CallVoidMethod(env, jmodel, mid, *mode, *n, *nstate, *nprob);  // pass state
  (*env)->CallVoidMethod(env, jmodel, mid);  // get state from primal_solution
  
  if (VERBOSE>0) fprintf(out, "%s In C, back from omi_model_getobjective__\n", prefix);
  
  // get f from jmodel
  fid = (*env)->GetFieldID(env, cls, "objective_value", "D");
  if (fid == 0) {
    fprintf(out, "%s !ERROR: Cannot GetFieldID of 'objective_value'\n", prefix);
    fclose(out);
    return;
  };
  model->objective_value = (*env)->GetDoubleField(env, jmodel, fid);
  *f = model->objective_value;
  if (VERBOSE>1) {
    fprintf(out, "%s model->objective_value = %g\n", prefix, model->objective_value);
    fprintf(out, "%s f = %X, *f = %g\n", prefix, f, *f);
  }
  
  // transfer jmodel.jobjective to g 
  if ((model->objective=(*env)->GetDoubleArrayElements(env,model->jobjective,0))==NULL) {
    fprintf(out,"%s !ERROR, Cannot GetDoubleArrayElements of 'model->jobjective'.\n",prefix);
    fclose(out);
    return;
  };
  for (i = 0; i < *n; i++)  g[i] = model->objective[i];  // hack tbd pass a ref?
  
  if (VERBOSE>1) {
    fprintf(out, "%s &model->objective = %X, model->objective = %X, *model->objective = %g\n", 
      prefix, &model->objective, model->objective, *model->objective);
    fprintf(out, "%s &g = %X, g = %X, *g = %g\n", prefix, &g, g, *g);
    for (i = 0; i < *n; i++) 
      fprintf(out, "Model.c:omi_model_getobjective__(): &g[%i] = %X, g[%i] = %g\n",i,&g[i],i,g[i]);
  }

    if (VERBOSE>0) fprintf(out,"%s return\n",prefix);
  fclose(out);
  return;

}  //  omi_model_getobjective__()

