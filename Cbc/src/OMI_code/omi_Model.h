/**
  * @(#)omi_Model.h
  * @author Robert Entriken<br>
  * Copyright (C) 1999-2007 EPRI<br>
  * All rights reserved.
  * @version 07-03-07
  * @since OMI_1.0
  *
  * Based on attributes of Model.java
  *
  * Revisions:
  *   08mar00
  *   24mar02  bobe     added nprob
  *   04-05-22 bobe     added objective_direction
  *   06-02-20 bobe     added copyright
  *   06-04-10 bobe     ANSI standard comments
  *                     JNI calls now return jint instead of long
  *   07-03-07 bobe     updated Copyright
  *
  */

#ifndef _FN_Model_java
#define _FN_Model_java

/* from f2c.h */
#ifndef F2C_INCLUDE
#define F2C_INCLUDE
  typedef long int integer;
  typedef double   doublereal;
  typedef long int ftnlen;
#endif


#include "omi_Block.h"
typedef struct ClassOMI_Model {
    jobject               instance; /* pointer to this */
    struct ClassObject    *header;
    int                   nonlin_constr_length;
    int                   nonlin_obj_length;
    int                   nonlin_jac_length;
    int                   objective_row;
    int                   objective_direction; /* 1 for min, -1 for max, 0 for none */
    double                objective_constant;
    int                   model_names_length;
    int                   rim_names_length;
    int                   jacobian_length;
    int                   weightsSOS1_length;
    int                   weightsSOS2_length;
    int                   max_row_length;
    int                   max_col_length;
    int                   row_length;
    int                   col_length;
    int                   workspace_length;
    jcharArray            jmodel_names;
    jchar                 *model_names;
    jcharArray            jrim_names1;
    jchar                 *rim_names1;
    jcharArray            jrim_names2;
    jchar                 *rim_names2;
    jobjectArray          jJacobianArray;
    jobject               *jJacobian;        /* Array of Java Blocks */
    ClassOMI_Block        **jacobian;
    jobjectArray          jWeightsSOS1Array;
    jobject               *jWeightsSOS1;     /* Array of Java Blocks */
    ClassOMI_Block        **weightsSOS1;
    jobjectArray          jWeightsSOS2Array;
    jobject               *jWeightsSOS2;     /* Array of Java Blocks */
    ClassOMI_Block        **weightsSOS2;
    jdoubleArray          jobjective;
    doublereal            *objective;
    jdoubleArray          jrow_lower_bound;
    double                *row_lower_bound;
    jdoubleArray          jrow_upper_bound;
    double                *row_upper_bound;
    jdoubleArray          jcol_lower_bound;
    double                *col_lower_bound;
    jdoubleArray          jcol_upper_bound;
    double                *col_upper_bound;
    jintArray             jis_integer;
    jint                  *is_integer;
    int                   ispecs;
    int                   iprint;
    int                   isumm;
    int                   nout;
    int                   nprob;
    int                   iterations_limit;
    int                   number_superbasics;
    int                   number_infeasible;
    double                sum_infeasible;
    double                objective_value;
    jintArray             jbasis;
    jint                  *basis;
    jdoubleArray          jprimal_solution;
    double                *primal_solution;
    jdoubleArray          jdual_solution;
    double                *dual_solution;
    jdoubleArray          jreduced_costs;
    double                *reduced_costs;
    jdoubleArray          jworkspace;
    double                *workspace;
    jdoubleArray          jgetObjectiveArg;
    double                *getObjectiveArg;
} ClassOMI_Model;

/* function prototypes */
void OMI_Model_release (JNIEnv *env, jobject jmodel, ClassOMI_Model *model, jint mode);
ClassOMI_Model *OMI_Model_get (JNIEnv *env, jobject jmodel);
JNIEXPORT void JNICALL OMI_Model_oprint
  (JNIEnv *env, jobject jmodel, jstring jprefix, FILE *out);

JNIEXPORT void JNICALL Java_com_epri_omi_Model_pprint 
  (JNIEnv *env, jobject jmodel, jstring jprefix);
JNIEXPORT jint JNICALL Java_com_epri_omi_Model_pMPS
  (JNIEnv *env, jobject jmodel, jstring jfilename);
JNIEXPORT jint JNICALL Java_com_epri_omi_Model_pestimate
  (JNIEnv *env, jobject jmodel, jstring jsolver_name);
JNIEXPORT jint JNICALL Java_com_epri_omi_Model_poption
  (JNIEnv *env, jobject jmodel, jstring joption_name);
JNIEXPORT jint JNICALL Java_com_epri_omi_Model_poptioni
  (JNIEnv *env, jobject jmodel, jstring joption_name, jint value);
JNIEXPORT jint JNICALL Java_com_epri_omi_Model_poptiond
  (JNIEnv *env, jobject jmodel, jstring joption_name, jdouble value);
JNIEXPORT jint JNICALL Java_com_epri_omi_Model_psetup
  (JNIEnv *env, jobject jmodel, jstring jsolver_name);
JNIEXPORT jint JNICALL Java_com_epri_omi_Model_psolve
  (JNIEnv *env, jobject jmodel, jstring jsolver_name);
JNIEXPORT jint JNICALL Java_com_epri_omi_Model_pwriteMPS
  (JNIEnv *env, jobject jmodel, jstring jfile_name);

#endif
