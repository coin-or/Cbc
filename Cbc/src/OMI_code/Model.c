/**
  * @(#)Model.c
  * @author Robert Entriken<br>
  * Copyright (C) 1999-2007 EPRI<br>
  * All rights reserved.
  * @version 07-03-07
  * @since OMI_1.0
  *
  * Revisions:
  *   02-08-15  bobe    psolve returns 999 on ERROR
  *   02-08-16  bobe    added more consistent use of nprob for OMI callbacks
  *   04-10-18  bobe    turned off scaling in call to Cbc_scaling()
  *                     Cbc_setIntegerTolerance() set to 1.0e-5 
  *   06-02-20  bobe    aborting on both SOS1 and SOS2 in same model
  *                     using Cbc_addSOS_Sparse()
  *                     added lots of debugging prints
  *                     weights must have row dimension of the model or error
  *                     using Cbc_setLogLevel(1)
  *   06-03-30  jmm     returning jint instead of long in JNI methods
  *   07-03-07 bobe     updated Copyright
  *
  */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "omi_Model.h"
#include "omi_State.h"
#include "com_epri_omi_Block.h"
#include "com_epri_omi_Model.h"
#include "Cbc_C_Interface.h"

#ifndef VERBOSE_INCLUDE
#define VERBOSE_INCLUDE
const int VERBOSE = 0;
#endif

  // local function prototypes
  void index_2_rowIndex(const int    *from, jint *to, int start, int end);
  void rowIndex_2_index(const jint   *from, int  *to, int start, int end);
  void colIndex_2_index(const jint   *from, int  *to, int start, int end);
  void char_2_long     (const char   *from, jint *to, int start, int end);
  void long_2_char     (const jint   *from, char *to, int start, int end);

  void starts_2_colIndex(const CoinBigIndex *from, jint           *to, int start, int end);
  void index_2_starts   (const jint         *from, CoinBigIndex   *to, int start, int end);

  void status_2_basis(const unsigned char *from, long *to, int start, int end);
  void basis_2_status(const long *from, unsigned char  *to, int start, int end);

  void dcopy(const double *from, double *to, int start, int end);
  void jccopy(const jchar *from, jchar *to, int start, int end);
  void iadd(int val, int *array, int start, int end);

  void unsigned_char_print(FILE *out, const char *prefix, char* name, const unsigned char *array, int start, int end);
  void char_print        (FILE *out, const char *prefix, char* name, const char   *array, int start, int end);
  void int_print         (FILE *out, const char *prefix, char* name, const int    *array, int start, int end);
  void long_print        (FILE *out, const char *prefix, char* name, const long   *array, int start, int end);
  void CoinBigIndex_print(FILE *out, const char *prefix, char* name, CoinBigIndex *array, int start, int end);
  void double_print      (FILE *out, const char *prefix, char* name, const double *array, int start, int end);

  void printSolution(Cbc_Model *cbc_model);

/**
  *  Call back function - just says whenever it gets Clp0005 or Coin0002
  *
  */
static void callBack(Cbc_Model * cbc_model, int messageNumber,
         int nDouble, const double * vDouble,
         int nInt, const int * vInt,
         int nString, char ** vString) 
{
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
}

// the native methods

/**
  *  Copy model attributes from C back to Java and release any allocated memory
  *
  */
void OMI_Model_release
  (JNIEnv *env, jobject jmodel, ClassOMI_Model *model, jint mode)
{
  const char prefix[] = "Model.c:OMI_Model_release(cbc) ";
//	const int VERBOSE = 2;

  jclass         cls;
  jfieldID       fid;
  int            i;
    
  if (VERBOSE>0) printf("%sbegin\n",prefix);
  // Check arguments.
  if (env == NULL) {
    printf("%s!ERROR, NULL 'env' passed as argument.\n", prefix);
    //fflush(stdout);
    return;
  };  
  if (jmodel == NULL) {
    printf("%s!ERROR, NULL 'jmodel' passed as argument.\n", prefix);
    //fflush(stdout);
    return;
  };

  // Announce OMI_Model_release
  if (VERBOSE>0) printf("%sArguments OK.\n",prefix);

  if ((mode == NULL) || (mode == OMI_COMMIT) ) {
    // cls
    // In JDK1.1 use the class to get at fields.
    if ((cls = (*env)->GetObjectClass(env, jmodel)) == NULL) {
      printf("Model.c:OMI_Model_release "\
        "!ERROR, Cannot GetObjectClass of 'jmodel'.\n");
      //fflush(stdout);
      return;
    };
    
    // nonlin_constr_length
    if (VERBOSE>1) printf("%s nonlin_constr_length\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "nonlin_constr_length", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.nonlin_constr_length\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->nonlin_constr_length); 
    
    // nonlin_obj_length
    if (VERBOSE>1) printf("%s nonlin_obj_length\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "nonlin_obj_length", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.nonlin_obj_length\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->nonlin_obj_length); 

    // nonlin_jac_length
    if (VERBOSE>1) printf("%s nonlin_jac_length\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "nonlin_jac_length", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.nonlin_jac_length\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->nonlin_jac_length); 

    // objective_row
    if (VERBOSE>1) printf("%s objective_row\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "objective_row", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.objective_row\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->objective_row); 

    // objective_direction
    if (VERBOSE>1) printf("%s objective_direction\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "objective_direction", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.objective_direction\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->objective_direction); 

    // objective_constant
    if (VERBOSE>1) printf("%s objective_constant\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "objective_constant", "D");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.objective_constant\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetDoubleField(env, jmodel, fid, model->objective_constant); 

    // model_names_length
    if (VERBOSE>1) printf("%s model_names_length\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "model_names_length", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.model_names_length\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->model_names_length); 

    // rim_names_length
    if (VERBOSE>1) printf("%s rim_names_length\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "rim_names_length", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.rim_names_length\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->rim_names_length); 

    // jacobian_length
    if (VERBOSE>1) printf("%s jacobian_length\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "jacobian_length", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.jacobian_length\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->jacobian_length); 

    // weightsSOS1_length
    if (VERBOSE>1) printf("%s weightsSOS1_length\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "weightsSOS1_length", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.weightsSOS1_length\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->weightsSOS1_length); 

    // weightsSOS2_length
    if (VERBOSE>1) printf("%s weightsSOS2_length\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "weightsSOS2_length", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.weightsSOS2_length\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->weightsSOS2_length); 

    // max_row_length
    if (VERBOSE>1) printf("%s max_row_length\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "max_row_length", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.max_row_length\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->max_row_length); 

    // max_col_length
    if (VERBOSE>1) printf("%s max_col_length\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "max_col_length", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.max_col_length\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->max_col_length); 

    // row_length
    if (VERBOSE>1) printf("%s row_length\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "row_length", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.row_length\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->row_length); 

    // col_length
    if (VERBOSE>1) printf("%s col_length\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "col_length", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.col_length\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->col_length); 

    // workspace_length
    if (VERBOSE>1) printf("%s workspace_length\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "workspace_length", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.workspace_length\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->workspace_length); 

    // ispecs
    if (VERBOSE>1) printf("%s ispecs\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "ispecs", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.ispecs\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->ispecs); 

    // iprint
    if (VERBOSE>1) printf("%s iprint\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "iprint", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.iprint\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->iprint); 

    // isumm
    if (VERBOSE>1) printf("%s isumm\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "isumm", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.isumm\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->isumm); 

    // nout
    if (VERBOSE>1) printf("%s nout\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "nout", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.nout\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->nout); 

    // nprob
    if (VERBOSE>1) printf("%s nprob\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "nprob", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.nprob\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->nprob); 

    // iterations_limit
    if (VERBOSE>1) printf("%s iterations_limit\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "iterations_limit", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.iterations_limit\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->iterations_limit); 

    // number_superbasics
    if (VERBOSE>1) printf("%s number_superbasics\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "number_superbasics", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.number_superbasics\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->number_superbasics); 

    // number_infeasible
    if (VERBOSE>1) printf("%s number_infeasible\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "number_infeasible", "I");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.number_infeasible\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetIntField(env, jmodel, fid, model->number_infeasible); 

    // sum_infeasible
    if (VERBOSE>1) printf("%s sum_infeasible\n",prefix);
    fid = (*env)->GetFieldID(env, cls, "sum_infeasible", "D");
    if (fid == 0) {
      printf ("Model.c:OMI_Model_release "\
        "!ERROR: Cannot GetFieldID of model.sum_infeasible\n");
      //fflush(stdout);
      return;
    }
    (*env)->SetDoubleField(env, jmodel, fid, model->sum_infeasible); 

    // objective_value
    if (VERBOSE>1) printf("%s objective_value = %f\n",prefix,model->objective_value);
    fid = (*env)->GetFieldID(env, cls, "objective_value", "D");
    if (fid == 0) {
      printf ("%sERROR: Cannot GetFieldID of model.objective_value\n",prefix);
      //fflush(stdout);
      return;
    }
    (*env)->SetDoubleField(env, jmodel, fid, model->objective_value); 

  };  // if NULL or OMI_COMMIT
  
  if (VERBOSE>1) printf("%s Arrays\n",prefix);
  if ((mode == NULL) ||  (mode == OMI_ABORT)) {
    if (VERBOSE>1) printf("%s model_names\n",prefix);
    (*env)->ReleaseCharArrayElements(env, model->jmodel_names, model->model_names, mode);
    if (VERBOSE>1) printf("%s rim_names1\n",prefix);
    (*env)->ReleaseCharArrayElements(env, model->jrim_names1, model->rim_names1, mode);
    if (VERBOSE>1) printf("%s rim_names2\n",prefix);
    (*env)->ReleaseCharArrayElements(env, model->jrim_names2, model->rim_names2, mode);

    // release the arrays of jacobian blocks
    if (VERBOSE>1) printf("%s jacobian\n",prefix);
    for (i=0; i<model->jacobian_length; i++) 
      OMI_Block_release(env, model->jJacobian[i], model->jacobian[i], mode);  
    free(model->jJacobian); free(model->jacobian);
    
    // release the arrays of SOS1 weight blocks
    if (VERBOSE>1) printf("%s weightsSOS1\n",prefix);
    for (i=0; i<model->weightsSOS1_length; i++) 
      if (model->jWeightsSOS1[i] != NULL) // null weight blocks are allowed
        OMI_Block_release(env, model->jWeightsSOS1[i], model->weightsSOS1[i], mode);  
    free(model->jWeightsSOS1); free(model->weightsSOS1);
    
    // release the arrays of SOS2 weight blocks
    if (VERBOSE>1) printf("%s weightsSOS2\n",prefix);
    for (i=0; i<model->weightsSOS2_length; i++) 
      if (model->jWeightsSOS2[i] != NULL) // null weight blocks are allowed
        OMI_Block_release(env, model->jWeightsSOS2[i], model->weightsSOS2[i], mode);  
    free(model->jWeightsSOS2); free(model->weightsSOS2);
    
    if (VERBOSE>1) printf("%s objective\n",prefix);
    (*env)->ReleaseDoubleArrayElements(env, model->jobjective, model->objective, mode);
    if (VERBOSE>1) printf("%s row_lower_bound\n",prefix);
    (*env)->ReleaseDoubleArrayElements(env, model->jrow_lower_bound, model->row_lower_bound, mode);
    if (VERBOSE>1) printf("%s row_upper_bound\n",prefix);
    (*env)->ReleaseDoubleArrayElements(env, model->jrow_upper_bound, model->row_upper_bound, mode);
    if (VERBOSE>1) printf("%s col_lower_bound\n",prefix);
    (*env)->ReleaseDoubleArrayElements(env, model->jcol_lower_bound, model->col_lower_bound, mode);
    if (VERBOSE>1) printf("%s col_upper_bound\n",prefix);
    (*env)->ReleaseDoubleArrayElements(env, model->jcol_upper_bound, model->col_upper_bound, mode);
    if (VERBOSE>1) printf("%s is_integer\n",prefix);
    (*env)->ReleaseIntArrayElements(env, model->jis_integer, model->is_integer, mode);

    if (VERBOSE>1) printf("%s basis\n",prefix);
    (*env)->ReleaseIntArrayElements(env, model->jbasis, model->basis, mode);

    if (VERBOSE>1) printf("%s primal_solution\n",prefix);
    (*env)->ReleaseDoubleArrayElements(env, model->jprimal_solution, model->primal_solution, mode);
    if (VERBOSE>1) printf("%s dual_solution\n",prefix);
    (*env)->ReleaseDoubleArrayElements(env, model->jdual_solution, model->dual_solution, mode);
    if (VERBOSE>1) printf("%s reduced_costs\n",prefix);
    (*env)->ReleaseDoubleArrayElements(env, model->jreduced_costs, model->reduced_costs, mode);
    if (VERBOSE>1) printf("%s workspace\n",prefix);
    (*env)->ReleaseDoubleArrayElements(env, model->jworkspace, model->workspace, mode);

    free(model);
  }  // if NULL or OMI_ABORT
  if (VERBOSE > 0) printf("%sreturn\n",prefix);
} // OMI_Model_release 

/**
  *  Copy model attributes from Java into C for use here.  Result is NULL on error.
  *
  */
ClassOMI_Model *OMI_Model_get
  (JNIEnv *env, jobject jmodel)
{ 
  const char prefix[] = "Model.c:OMI_Model_get(cbc) ";  
//  const int VERBOSE = 2;
  
  jclass          cls;
  jfieldID        fid;
  jstring         jnewPrefix;
  ClassOMI_Model  *model;
  jthrowable      jexcept;
  
  integer i,j;
  
  char  newPrefix[OMI_MaxStrLen], 
        blockName[OMI_MaxStrLen];

  // external functions
  extern void       Java_com_epri_omi_block_convertToRowOrder (JNIEnv *, jobject);

  if (VERBOSE>0) printf("%sbegin\n",prefix);
  // Check arguments.
  if (env == NULL) {
    printf("%s!ERROR, NULL 'env' passed as argument.\n", prefix);
    //fflush(stdout);
    return NULL;
  };  
  if (jmodel == NULL) {
    printf("%s!ERROR, NULL 'jmodel' passed as argument.\n", prefix);
    //fflush(stdout);
    return NULL;
  };
    
  // Announce OMI_Model_get
  if (VERBOSE > 0) printf("%sArguments OK.\n",prefix);

  // malloc space for a model
  i = sizeof(ClassOMI_Model);
  if ((model = (ClassOMI_Model *) malloc(sizeof(ClassOMI_Model))) == NULL) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot malloc space for model\n");
    //fflush(stdout);
    return NULL;
  };

  // cls	
  // In JDK1.1 use the class to get at fields.
  if ((cls = (*env)->GetObjectClass(env, jmodel)) == NULL) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetObjectClass of 'jmodel'.\n");
    //fflush(stdout);
    return NULL;
  };
  
  // pointer to this
  model->instance = jmodel;

  // load local refs to attributes
  fid = (*env)->GetFieldID(env, cls, "nonlin_constr_length", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'nonlin_constr_length'\n");
    //fflush(stdout);
    return NULL;
  };
  model->nonlin_constr_length = (*env)->GetIntField(env,jmodel,fid);

  fid = (*env)->GetFieldID(env, cls, "nonlin_obj_length", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'nonlin_obj_length'\n");
    //fflush(stdout);
    return NULL;
  };
  model->nonlin_obj_length = (*env)->GetIntField(env,jmodel,fid);

  fid = (*env)->GetFieldID(env, cls, "nonlin_jac_length", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'nonlin_jac_length'\n");
    //fflush(stdout);
    return NULL;
  };
  model->nonlin_jac_length = (*env)->GetIntField(env,jmodel,fid);

  fid = (*env)->GetFieldID(env, cls, "objective_row", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'objective_row'\n");
    //fflush(stdout);
    return NULL;
  };
  model->objective_row = (*env)->GetIntField(env,jmodel,fid);

  fid = (*env)->GetFieldID(env, cls, "objective_direction", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'objective_direction'\n");
    //fflush(stdout);
    return NULL;
  };
  model->objective_direction = (*env)->GetIntField(env,jmodel,fid);

  fid = (*env)->GetFieldID(env, cls, "objective_constant", "D");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'objective_constant'\n");
    //fflush(stdout);
    return NULL;
  };
  model->objective_constant = (*env)->GetDoubleField(env,jmodel,fid);

  fid = (*env)->GetFieldID(env, cls, "model_names_length", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'model_names_length'\n");
    //fflush(stdout);
    return NULL;
  };
  model->model_names_length = (*env)->GetIntField(env,jmodel,fid);

  fid = (*env)->GetFieldID(env, cls, "rim_names_length", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'rim_names_length'\n");
    //fflush(stdout);
    return NULL;
  };
  model->rim_names_length = (*env)->GetIntField(env,jmodel,fid);

  // jacobian_length
  fid = (*env)->GetFieldID(env, cls, "jacobian_length", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'jacobian_length'\n");
    //fflush(stdout);
    return NULL;
  };
  model->jacobian_length = (*env)->GetIntField(env,jmodel,fid);
  
  // weightsSOS1_length
  fid = (*env)->GetFieldID(env, cls, "weightsSOS1_length", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'weightsSOS1_length'\n");
    //fflush(stdout);
    return NULL;
  };
  model->weightsSOS1_length = (*env)->GetIntField(env,jmodel,fid);

  // weightsSOS2_length
  fid = (*env)->GetFieldID(env, cls, "weightsSOS2_length", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'weightsSOS2_length'\n");
    //fflush(stdout);
    return NULL;
  };
  model->weightsSOS2_length = (*env)->GetIntField(env,jmodel,fid);

  // max_row_length
  fid = (*env)->GetFieldID(env, cls, "max_row_length", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'max_row_length'\n");
    //fflush(stdout);
    return NULL;
  };
  model->max_row_length = (*env)->GetIntField(env,jmodel,fid);

  fid = (*env)->GetFieldID(env, cls, "max_col_length", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'max_col_length'\n");
    //fflush(stdout);
    return NULL;
  };
  model->max_col_length = (*env)->GetIntField(env,jmodel,fid);

  fid = (*env)->GetFieldID(env, cls, "row_length", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'row_length'\n");
    //fflush(stdout);
    return NULL;
  };
  model->row_length = (*env)->GetIntField(env,jmodel,fid);

  fid = (*env)->GetFieldID(env, cls, "col_length", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'col_length'\n");
    //fflush(stdout);
    return NULL;
  };
  model->col_length = (*env)->GetIntField(env,jmodel,fid);

  fid = (*env)->GetFieldID(env, cls, "workspace_length", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'workspace_length'\n");
    //fflush(stdout);
    return NULL;
  };
  model->workspace_length = (*env)->GetIntField(env,jmodel,fid);
    
  // get model array attributes
  //
  // download model_names
  fid = (*env)->GetFieldID(env, cls, "model_names", "[C");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'model_names'\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->jmodel_names = (jcharArray) (*env)->GetObjectField(env,jmodel,fid)) == NULL) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetObjectField of 'model_names'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->model_names = (*env)->GetCharArrayElements(env, model->jmodel_names, 0)) == NULL) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetCharArrayElements of 'model->jmodel_names'.\n");
    //fflush(stdout);
    return NULL;
  };
  model->model_names_length = (*env)->GetArrayLength(env, model->jmodel_names);
  
  // download rim_names1
  fid = (*env)->GetFieldID(env, cls, "rim_names1", "[C");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'rim_names1'\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->jrim_names1 = (jcharArray) (*env)->GetObjectField(env,jmodel,fid)) == NULL) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetObjectField of 'rim_names1'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->rim_names1 = (*env)->GetCharArrayElements(env, model->jrim_names1, 0)) == NULL) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetCharArrayElements of 'model->jrim_names1'.\n");
    //fflush(stdout);
    return NULL;
  };
  // assuming lengths of rim_names1 and rim_names2 are identical
  model->rim_names_length = (*env)->GetArrayLength(env, model->jrim_names1);

  // download rim_names2
  fid = (*env)->GetFieldID(env, cls, "rim_names2", "[C");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'rim_names2'\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->jrim_names2 = (jcharArray) (*env)->GetObjectField(env,jmodel,fid)) == NULL) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetObjectField of 'rim_names2'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->rim_names2 = (*env)->GetCharArrayElements(env, model->jrim_names2, 0)) == NULL) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetCharArrayElements of 'model->jrim_names2'.\n");
    //fflush(stdout);
    return NULL;
  };
  // check that lengths of rim_names1 and rim_names2 are identical
  if ((*env)->GetArrayLength(env, model->jrim_names1)
    != (*env)->GetArrayLength(env, model->jrim_names2)) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, (GetArrayLength(rim_names1) = %d) != (GetArrayLength(rim_names2) = %d)\n",
        model->rim_names_length, (*env)->GetArrayLength(env, model->jrim_names2));
    //fflush(stdout);
    return NULL;
  };
  model->rim_names_length = (*env)->GetArrayLength(env, model->jrim_names2);

  // get the array of jacobian blocks
  fid = (*env)->GetFieldID(env, cls, "jacobian", "[Lcom/epri/omi/Block;");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'jacobian'\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->jJacobianArray = (*env)->GetObjectField(env, jmodel, fid)) == NULL) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetObjectField of 'jacobian'\n");
    //fflush(stdout);
    return NULL;
  };
  model->jacobian_length = (*env)->GetArrayLength(env, model->jJacobianArray);
  model->jJacobian = (jobject *        ) malloc(model->jacobian_length * sizeof(jobject         ));
  model->jacobian  = (ClassOMI_Block **) malloc(model->jacobian_length * sizeof(ClassOMI_Block *));
  for (i=0; i<model->jacobian_length; i++) {
    if ((model->jJacobian[i] = (*env)->GetObjectArrayElement(env, model->jJacobianArray, i)) == NULL) {
      printf("Model.c:OMI_Model_get "\
        "!ERROR, Cannot GetObjectArrayElements of 'model->jJacobianArray, %d'.\n", i);
      //fflush(stdout);
      return NULL;
    };
    // make sure the jacobian Blocks are in Column Order format
    Java_com_epri_omi_block_convertToColumnOrder(env, model->jJacobian[i]);
    model->jacobian[i] = (ClassOMI_Block *) OMI_Block_get(env, model->jJacobian[i]);
  };
  
  // get the array of SOS1 weight blocks
  fid = (*env)->GetFieldID(env, cls, "weightsSOS1", "[Lcom/epri/omi/Block;");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'weightsSOS1'\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->jWeightsSOS1Array = (*env)->GetObjectField(env, jmodel, fid)) == NULL) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetObjectField of 'weightsSOS1'\n");
    //fflush(stdout);
    return NULL;
  };
  model->weightsSOS1_length = (*env)->GetArrayLength(env, model->jWeightsSOS1Array);
  model->jWeightsSOS1 = (jobject *        ) malloc(model->weightsSOS1_length * sizeof(jobject         ));
  model->weightsSOS1  = (ClassOMI_Block **) malloc(model->weightsSOS1_length * sizeof(ClassOMI_Block *));
  for (i=0; i<model->weightsSOS1_length; i++) {
    if ((model->jWeightsSOS1[i] = (*env)->GetObjectArrayElement(env, model->jWeightsSOS1Array, i)) == NULL) {
      printf("Model.c:OMI_Model_get "\
        "!WARNING, Cannot GetObjectArrayElements of 'model->jWeightsSOS1Array, %d'.\n", i);
    } else {
      if (VERBOSE>3) {
        printf("%s i = %i\n",prefix,i);
        sprintf(blockName, "Model->weightsSOS1[%d] @ %X, ",i, 
          &model->jWeightsSOS1[i]);
        if (strlen(blockName) + strlen(prefix) > OMI_MaxStrLen) {
          fprintf(stdout, 
            "%s !ERROR: Cannot strcat,  strlen(newPrefix) (= %d) > %d\n",
            prefix, strlen(blockName) + strlen(prefix), OMI_MaxStrLen);
          fprintf(stdout, "  blockName = '%s'\n", blockName);
          fprintf(stdout, "  prefix    = '%s'\n", prefix   );
          OMI_Model_release(env, jmodel, model, OMI_ABORT);
          return NULL;
        };
        strcat(strcpy(newPrefix,prefix),blockName);  
        if ((jnewPrefix = (jstring) (*env)->NewStringUTF(env, newPrefix)) == NULL) {
          fprintf(stdout, "%s !ERROR: Cannot NewStringUTF of 'newPrefix'\n", prefix);
          OMI_Model_release(env, jmodel, model, OMI_ABORT);
          return NULL;
        };
        // pprint the ith SOS1 weight block
        OMI_Block_oprint(env, model->jWeightsSOS1[i], jnewPrefix, stdout); 
      }
      if (VERBOSE>1) {
        printf("%s calling Java_com_epri_omi_block_convertToRowOrder(jWeightsSOS1[%i])\n",
          prefix,i);
        fflush(stdout);
      }
      // make sure the weightsSOS1 are in Row Order format
      Java_com_epri_omi_block_convertToRowOrder(env, model->jWeightsSOS1[i]);
      if (VERBOSE>1) {
        printf("%s finished Java_com_epri_omi_block_convertToRowOrder(jWeightsSOS1[%i])\n",
          prefix,i);
      }
      model->weightsSOS1[i] = (ClassOMI_Block *) OMI_Block_get(env, model->jWeightsSOS1[i]);
    }
  }
  
  // get the array of SOS2 weight blocks
  fid = (*env)->GetFieldID(env, cls, "weightsSOS2", "[Lcom/epri/omi/Block;");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'weightsSOS2'\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->jWeightsSOS2Array = (*env)->GetObjectField(env, jmodel, fid)) == NULL) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetObjectField of 'weightsSOS2'\n");
    //fflush(stdout);
    return NULL;
  };
  model->weightsSOS2_length = (*env)->GetArrayLength(env, model->jWeightsSOS2Array);
  model->jWeightsSOS2 = (jobject *        ) malloc(model->weightsSOS2_length * sizeof(jobject         ));
  model->weightsSOS2  = (ClassOMI_Block **) malloc(model->weightsSOS2_length * sizeof(ClassOMI_Block *));
  for (i=0; i<model->weightsSOS2_length; i++) {
    if ((model->jWeightsSOS2[i] = (*env)->GetObjectArrayElement(env, model->jWeightsSOS2Array, i)) == NULL) {
      printf("Model.c:OMI_Model_get "\
        "!WARNING, Cannot GetObjectArrayElements of 'model->jWeightsSOS2Array, %d'.\n", i);
    } else {
      // make sure the weightsSOS2 are in Row Order format
      Java_com_epri_omi_block_convertToRowOrder(env, model->jWeightsSOS2[i]); 
      model->weightsSOS2[i] = (ClassOMI_Block *) OMI_Block_get(env, model->jWeightsSOS2[i]);
    }
  }
  
  // objective
  fid = (*env)->GetFieldID(env, cls, "objective", "[D");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetFieldID of 'objective'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->jobjective = (*env)->GetObjectField(env,jmodel,fid)) == NULL) {
      printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetObjectField of 'objective'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->objective = (*env)->GetDoubleArrayElements(env, model->jobjective, 0)) == NULL) {
      printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetDoubleArrayElements of 'model->jobjective'.\n");
    //fflush(stdout);
    return NULL;
  };
  if (0 && model->max_col_length != (*env)->GetArrayLength(env, model->jobjective)) {
      printf("Model.c:OMI_Model_get "\
      "!ERROR, model->max_col_length (= %d) != "\
         "GetArrayLength(env, model->jobjective) (= %d)\n",
         model->max_col_length, (*env)->GetArrayLength(env, model->jobjective));
    //fflush(stdout);
    return NULL;
  };

  // row_lower_bound
  fid = (*env)->GetFieldID(env, cls, "row_lower_bound", "[D");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetFieldID of 'row_lower_bound'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->jrow_lower_bound = (*env)->GetObjectField(env,jmodel,fid)) == NULL) {
      printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetObjectField of 'row_lower_bound'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->row_lower_bound = (*env)->GetDoubleArrayElements(env, model->jrow_lower_bound, 0)) == NULL) {
      printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetDoubleArrayElements of 'model->jrow_lower_bound'.\n");
    //fflush(stdout);
    return NULL;
  };
  if (model->max_row_length != (*env)->GetArrayLength(env, model->jrow_lower_bound)) {
      printf("Model.c:OMI_Model_get "\
      "!ERROR, model->max_row_length (= %d) != "\
         "GetArrayLength(env, model->jrow_lower_bound) (= %d)\n",
         model->max_row_length, (*env)->GetArrayLength(env, model->jrow_lower_bound));
    //fflush(stdout);
    return NULL;
  };

  // row_upper_bound
  fid = (*env)->GetFieldID(env, cls, "row_upper_bound", "[D");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetFieldID of 'row_upper_bound'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->jrow_upper_bound = (*env)->GetObjectField(env,jmodel,fid)) == NULL) {
      printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetObjectField of 'row_upper_bound'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->row_upper_bound = (*env)->GetDoubleArrayElements(env, model->jrow_upper_bound, 0)) == NULL) {
      printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetDoubleArrayElements of 'model->jrow_upper_bound'.\n");
    //fflush(stdout);
    return NULL;
  };
  if (model->max_row_length != (*env)->GetArrayLength(env, model->jrow_upper_bound)) {
      printf("Model.c:OMI_Model_get "\
      "!ERROR, model->max_row_length (= %d) != "\
         "GetArrayLength(env, model->jrow_upper_bound) (= %d)\n",
         model->max_row_length, (*env)->GetArrayLength(env, model->jrow_upper_bound));
    //fflush(stdout);
    return NULL;
  };

  // col_lower_bound
  fid = (*env)->GetFieldID(env, cls, "col_lower_bound", "[D");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetFieldID of 'col_lower_bound'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->jcol_lower_bound = (*env)->GetObjectField(env,jmodel,fid)) == NULL) {
      printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetObjectField of 'col_lower_bound'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->col_lower_bound = (*env)->GetDoubleArrayElements(env, model->jcol_lower_bound, 0)) == NULL) {
      printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetDoubleArrayElements of 'model->jcol_lower_bound'.\n");
    //fflush(stdout);
    return NULL;
  };
  if (model->max_col_length != (*env)->GetArrayLength(env, model->jcol_lower_bound)) {
      printf("Model.c:OMI_Model_get "\
      "!ERROR, model->max_col_length (= %d) != "\
         "GetArrayLength(env, model->jcol_lower_bound) (= %d)\n",
         model->max_col_length, (*env)->GetArrayLength(env, model->jcol_lower_bound));
    //fflush(stdout);
    return NULL;
  };

  // col_upper_bound
  fid = (*env)->GetFieldID(env, cls, "col_upper_bound", "[D");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetFieldID of 'col_upper_bound'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->jcol_upper_bound = (*env)->GetObjectField(env,jmodel,fid)) == NULL) {
      printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetObjectField of 'col_upper_bound'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->col_upper_bound = (*env)->GetDoubleArrayElements(env, model->jcol_upper_bound, 0)) == NULL) {
      printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetDoubleArrayElements of 'model->jcol_upper_bound'.\n");
    //fflush(stdout);
    return NULL;
  };
  if (model->max_col_length != (*env)->GetArrayLength(env, model->jcol_upper_bound)) {
      printf("Model.c:OMI_Model_get "\
      "!ERROR, model->max_col_length (= %d) != "\
         "GetArrayLength(env, model->jcol_upper_bound) (= %d)\n",
         model->max_col_length, (*env)->GetArrayLength(env, model->jcol_upper_bound));
    //fflush(stdout);
    return NULL;
  };

  // is_integer
  fid = (*env)->GetFieldID(env, cls, "is_integer", "[I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetFieldID of 'is_integer' with sig '[I'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->jis_integer = (*env)->GetObjectField(env,jmodel,fid)) == NULL) {
      printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetObjectField of 'is_integer'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->is_integer = (*env)->GetIntArrayElements(env, model->jis_integer, 0)) == NULL) {
      printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetDoubleArrayElements of 'model->jis_integer'.\n");
    //fflush(stdout);
    return NULL;
  };
  if (model->col_length > (*env)->GetArrayLength(env, model->jis_integer)) {
      printf("Model.c:OMI_Model_get "\
      "!ERROR, model->col_length (= %d) > "\
         "GetArrayLength(env, model->jis_integer) (= %d)\n",
         model->col_length, (*env)->GetArrayLength(env, model->jis_integer));
    //fflush(stdout);
    return NULL;
  };

  // load up on miscellaneous ints and doubles
  fid = (*env)->GetFieldID(env, cls, "ispecs", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'ispecs'\n");
    //fflush(stdout);
    return NULL;
  };
  model->ispecs = (*env)->GetIntField(env,jmodel,fid);

  fid = (*env)->GetFieldID(env, cls, "iprint", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'iprint'\n");
    //fflush(stdout);
    return NULL;
  };
  model->iprint = (*env)->GetIntField(env,jmodel,fid);

  fid = (*env)->GetFieldID(env, cls, "isumm", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'isumm'\n");
    //fflush(stdout);
    return NULL;
  };
  model->isumm = (*env)->GetIntField(env,jmodel,fid);

  fid = (*env)->GetFieldID(env, cls, "nout", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'nout'\n");
    //fflush(stdout);
    return NULL;
  };
  model->nout = (*env)->GetIntField(env,jmodel,fid);

  fid = (*env)->GetFieldID(env, cls, "nprob", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'nprob'\n");
    //fflush(stdout);
    return NULL;
  };
  model->nprob = (*env)->GetIntField(env,jmodel,fid);

  fid = (*env)->GetFieldID(env, cls, "iterations_limit", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'iterations_limit'\n");
    //fflush(stdout);
    return NULL;
  };
  model->iterations_limit = (*env)->GetIntField(env,jmodel,fid);

  fid = (*env)->GetFieldID(env, cls, "number_superbasics", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'number_superbasics'\n");
    //fflush(stdout);
    return NULL;
  };
  model->number_superbasics = (*env)->GetIntField(env,jmodel,fid);

  fid = (*env)->GetFieldID(env, cls, "number_infeasible", "I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'number_infeasible'\n");
    //fflush(stdout);
    return NULL;
  };
  model->number_infeasible = (*env)->GetIntField(env,jmodel,fid);

  fid = (*env)->GetFieldID(env, cls, "sum_infeasible", "D");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'sum_infeasible'\n");
    //fflush(stdout);
    return NULL;
  };
  model->sum_infeasible = (*env)->GetDoubleField(env,jmodel,fid);

  fid = (*env)->GetFieldID(env, cls, "objective_value", "D");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR: Cannot GetFieldID of 'objective_value'\n");
    //fflush(stdout);
    return NULL;
  };
  model->objective_value = (*env)->GetDoubleField(env,jmodel,fid);
    
  // basis
  fid = (*env)->GetFieldID(env, cls, "basis", "[I");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetFieldID of 'basis'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->jbasis = (*env)->GetObjectField(env,jmodel,fid)) == NULL) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetObjectField of 'basis'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->basis = (*env)->GetIntArrayElements(env, model->jbasis, 0)) == NULL) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetIntArrayElements of 'model->jbasis'.\n");
    //fflush(stdout);
    return NULL;
  };
  if (model->max_row_length + model->max_col_length 
    != (*env)->GetArrayLength(env, model->jbasis)) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, model->max_row_length (= %d) + model->max_col_length (= %d) != "\
         "GetArrayLength(env, model->jbasis) (= %d)\n",
         model->max_row_length, model->max_col_length, 
         (*env)->GetArrayLength(env, model->jbasis));
    //fflush(stdout);
    return NULL;
  };

  // primal_solution
  fid = (*env)->GetFieldID(env, cls, "primal_solution", "[D");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetFieldID of 'primal_solution'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->jprimal_solution = (*env)->GetObjectField(env,jmodel,fid)) == NULL) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetObjectField of 'primal_solution'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->primal_solution = (*env)->GetDoubleArrayElements(env, model->jprimal_solution, 0)) == NULL) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetDoubleArrayElements of 'model->jprimal_solution'.\n");
    //fflush(stdout);
    return NULL;
  };
  if (model->max_col_length != (*env)->GetArrayLength(env, model->jprimal_solution)) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, model->max_col_length (= %d) != "\
         "GetArrayLength(env, model->jprimal_solution) (= %d)\n",
         model->max_col_length, (*env)->GetArrayLength(env, model->jprimal_solution));
    //fflush(stdout);
    return NULL;
  };

  // dual_solution
  fid = (*env)->GetFieldID(env, cls, "dual_solution", "[D");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetFieldID of 'dual_solution'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->jdual_solution = (*env)->GetObjectField(env,jmodel,fid)) == NULL) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetObjectField of 'dual_solution'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->dual_solution = (*env)->GetDoubleArrayElements(env, model->jdual_solution, 0)) == NULL) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetDoubleArrayElements of 'model->jdual_solution'.\n");
    //fflush(stdout);
    return NULL;
  };
  if (model->max_row_length != (*env)->GetArrayLength(env, model->jdual_solution)) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, model->max_row_length (= %d) != ",
         "GetArrayLength(env, model->jdual_solution) (= %d)\n",
         model->max_row_length, (*env)->GetArrayLength(env, model->jdual_solution));
    //fflush(stdout);
    return NULL;
  };

  // reduced_costs
  fid = (*env)->GetFieldID(env, cls, "reduced_costs", "[D");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetFieldID of 'reduced_costs'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->jreduced_costs = (*env)->GetObjectField(env,jmodel,fid)) == NULL) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetObjectField of 'reduced_costs'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->reduced_costs = (*env)->GetDoubleArrayElements(env, model->jreduced_costs, 0)) == NULL) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetDoubleArrayElements of 'model->jreduced_costs'.\n");
    //fflush(stdout);
    return NULL;
  };
  if (model->max_col_length != (*env)->GetArrayLength(env, model->jreduced_costs)) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, model->max_col_length (= %d) != ",
         "GetArrayLength(env, model->jreduced_costs) (= %d)\n",
         model->max_col_length, (*env)->GetArrayLength(env, model->jreduced_costs));
    //fflush(stdout);
    return NULL;
  };

  // workspace
  fid = (*env)->GetFieldID(env, cls, "workspace", "[D");
  if (fid == 0) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetFieldID of 'workspace'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->jworkspace = (*env)->GetObjectField(env,jmodel,fid)) == NULL) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetObjectField of 'workspace'.\n");
    //fflush(stdout);
    return NULL;
  };
  if ((model->workspace = (*env)->GetDoubleArrayElements(env, model->jworkspace, 0)) == NULL) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, Cannot GetDoubleArrayElements of 'model->jworkspace'.\n");
    //fflush(stdout);
    return NULL;
  };
  if (model->workspace_length != (*env)->GetArrayLength(env, model->jworkspace)) {
    printf("Model.c:OMI_Model_get "\
      "!ERROR, model->max_col_length (= %d) != ",
         "GetArrayLength(env, model->jworkspace) (= %d)\n",
         model->max_col_length, (*env)->GetArrayLength(env, model->jworkspace));
    //fflush(stdout);
    return NULL;
  };
  
  // success!
  if (VERBOSE > 0) printf("%sreturn\n",prefix);
  return model;
}  // OMI_Model_get

/**
  *  Print model attributes to the FILE argument.
  *
  */
JNIEXPORT void JNICALL OMI_Model_oprint
  (JNIEnv *env, jobject jmodel, jstring jprefix, FILE *out)
{ 
  jstring         jnewPrefix;
  ClassOMI_Model  *model;
  const char      *prefix;
  const char thisPrefix[] = "Model.c:OMI_Model_oprint() ";
  
  integer i,j;
  char newPrefix[OMI_MaxStrLen], 
       blockName[OMI_MaxStrLen];

  // Check arguments.
  if (env == NULL) {
    fprintf(out, "%s!ERROR, NULL 'env' passed as argument.\n", thisPrefix);
    //fflush(out);
    return;
  };  
  if (jmodel == NULL) {
    fprintf(out, "%s!ERROR, NULL 'jmodel' passed as argument.\n", thisPrefix);
    //fflush(out);
    return;
  };
  
  // Announce oprint
  // In JDK 1.1 get a local reference to a string.
  // The UTF string is like a C string.
  // The alternative is a unicode string.
  if ((prefix = (*env)->GetStringUTFChars(env, jprefix, 0)) == NULL) {
    fprintf(out, 
      "%s !ERROR, Cannot GetStringUTFChars of 'jprefix'.\n", thisPrefix);
    //fflush(out);
    return;
  };
  fprintf(out, "%s Model->pprint('%s') Arguments OK.\n", thisPrefix, prefix);
  
  // get a C version of the model
  //
  if ((model = OMI_Model_get(env, jmodel)) == NULL) {
    fprintf(out, "%s!ERROR, Cannot OMI_Model_get of 'jmodel'.\n", thisPrefix);
    //fflush(out);
    return;
  };
  
  // print model dimensions
  //
  fprintf(out, "%sModel->max_row_length       @ %X = %d\n", 
    prefix, &model->max_row_length, model->max_row_length);
  fprintf(out, "%sModel->max_col_length       @ %X = %d\n", 
    prefix, &model->max_col_length, model->max_col_length);
  fprintf(out, "%sModel->row_length           @ %X = %d\n", 
    prefix, &model->row_length, model->row_length);
  fprintf(out, "%sModel->col_length           @ %X = %d\n", 
    prefix, &model->col_length, model->col_length);
  fprintf(out, "%sModel->workspace_length     @ %X = %d\n", 
    prefix, &model->workspace_length, model->workspace_length);
  fprintf(out, "%sModel->nonlin_constr_length @ %X = %d\n", 
    prefix, &model->nonlin_constr_length, model->nonlin_constr_length);
  fprintf(out, "%sModel->nonlin_obj_length    @ %X = %d\n", 
    prefix, &model->nonlin_obj_length, model->nonlin_obj_length); 
  fprintf(out, "%sModel->nonlin_jac_length    @ %X = %d\n", 
    prefix, &model->nonlin_jac_length, model->nonlin_jac_length);
  fprintf(out, "%sModel->objective_row        @ %X = %d\n", 
    prefix, &model->objective_row, model->objective_row); 	  
  fprintf(out, "%sModel->objective_direction  @ %X = %d\n", 
    prefix, &model->objective_direction, model->objective_direction); 	  
  fprintf(out, "%sModel->objective_constant   @ %X = %g\n", 
    prefix, &model->objective_constant, model->objective_constant);
  fprintf(out, "%sModel->model_names_length   @ %X = %d\n",   
    prefix, &model->model_names_length, model->model_names_length);
  fprintf(out, "%sModel->rim_names_length     @ %X = %d\n",   
    prefix, &model->rim_names_length, model->rim_names_length);
  fprintf(out, "%sModel->jacobian_length         @ %X = %d\n",   
    prefix, &model->jacobian_length, model->jacobian_length);
  fprintf(out, "%sModel->weightsSOS1_length         @ %X = %d\n",   
    prefix, &model->weightsSOS1_length, model->weightsSOS1_length);
  fprintf(out, "%sModel->weightsSOS2_length         @ %X = %d\n",   
    prefix, &model->weightsSOS2_length, model->weightsSOS2_length);
    
  // print model array attributes
  //
  
  fprintf(out, "%sModel->model_names @ %X =\n  ", prefix, model->model_names);

  if (0) { // print unicode second bytes
    for (i=0; i<model->model_names_length/8; i++) {
      for (j=0; j<8; j++) fprintf(out, "%c", model->model_names[(i*8+j)*2+1]);
      fprintf(out, ", ");
    };
    fprintf(out, "\n");
  };
  // print UTF8 of a C string
  for (i=0; i<model->model_names_length/8; i++) {
    for (j=0; j<8; j++) fprintf(out, "%c", model->model_names[i*8+j]);
    fprintf(out, ", ");
  };
  fprintf(out, "\n");
  
  fprintf(out, "%sModel->rim_names1  @ %X =\n  ", prefix, model->rim_names1         );
  for (i=0; i<model->rim_names_length/8; i++) {
    for (j=0; j<8; j++) fprintf(out, "%c", model->rim_names1[i*8+j]);
    fprintf(out, ", ");
  };
  fprintf(out, "\n");

  fprintf(out, "%sModel->rim_names2  @ %X =\n  ", prefix, model->rim_names2         );
  for (i=0; i<model->rim_names_length/8; i++) {
    for (j=0; j<8; j++) fprintf(out, "%c", model->rim_names2[i*8+j]);
    fprintf(out, ", ");
  };
  fprintf(out, "\n");

  // loop over and print jacobian blocks in Model
  for (i=0; i<model->jacobian_length; i++) {
    sprintf(blockName, "Model->jacobian[%d] @ %X, ",i, &model->jJacobian[i]);
    if (strlen(blockName) + strlen(prefix) > OMI_MaxStrLen) {
      fprintf(out, "%s !ERROR: Cannot strcat,  strlen(newPrefix) (= %d) > %d\n",
        thisPrefix, strlen(blockName) + strlen(prefix), OMI_MaxStrLen);
      fprintf(out, "  blockName = '%s'\n", blockName);
      fprintf(out, "  prefix    = '%s'\n", prefix   );
      OMI_Model_release(env, jmodel, model, OMI_ABORT);
      return;
    };
    strcat(strcpy(newPrefix,prefix),blockName);  
    if ((jnewPrefix = (jstring) (*env)->NewStringUTF(env, newPrefix)) == NULL) {
      fprintf(out, "%s !ERROR: Cannot NewStringUTF of 'newPrefix'\n", thisPrefix);
      OMI_Model_release(env, jmodel, model, OMI_ABORT);
      return;
    };
    // pprint the ith jacobian block
    OMI_Block_oprint(env, model->jJacobian[i], jnewPrefix, out); 
  };
  
  // loop over and print non-null SOS1 weight blocks in Model
  for (i=0; i<model->weightsSOS1_length; i++) {
    if (model->jWeightsSOS1[i] != NULL) {
      sprintf(blockName, "Model->weightsSOS1[%d] @ %X, ",i, &model->jWeightsSOS1[i]);
      if (strlen(blockName) + strlen(prefix) > OMI_MaxStrLen) {
        fprintf(out, "%s !ERROR: Cannot strcat,  strlen(newPrefix) (= %d) > %d\n",
          thisPrefix, strlen(blockName) + strlen(prefix), OMI_MaxStrLen);
        fprintf(out, "  blockName = '%s'\n", blockName);
        fprintf(out, "  prefix    = '%s'\n", prefix   );
        OMI_Model_release(env, jmodel, model, OMI_ABORT);
        return;
      };
      strcat(strcpy(newPrefix,prefix),blockName);  
      if ((jnewPrefix = (jstring) (*env)->NewStringUTF(env, newPrefix)) == NULL) {
        fprintf(out, "%s !ERROR: Cannot NewStringUTF of 'newPrefix'\n", thisPrefix);
        OMI_Model_release(env, jmodel, model, OMI_ABORT);
        return;
      };
      // pprint the ith SOS1 weight block
      OMI_Block_oprint(env, model->jWeightsSOS1[i], jnewPrefix, out); 
    }
  }
  
  // loop over and print non-null SOS2 weight blocks in Model
  for (i=0; i<model->weightsSOS2_length; i++) {
    if (model->jWeightsSOS2[i] != NULL) {
      sprintf(blockName, "Model->weightsSOS2[%d] @ %X, ",i, &model->jWeightsSOS2[i]);
      if (strlen(blockName) + strlen(prefix) > OMI_MaxStrLen) {
        fprintf(out, "%s !ERROR: Cannot strcat,  strlen(newPrefix) (= %d) > %d\n",
          thisPrefix, strlen(blockName) + strlen(prefix), OMI_MaxStrLen);
        fprintf(out, "  blockName = '%s'\n", blockName);
        fprintf(out, "  prefix    = '%s'\n", prefix   );
        OMI_Model_release(env, jmodel, model, OMI_ABORT);
        return;
      };
      strcat(strcpy(newPrefix,prefix),blockName);  
      if ((jnewPrefix = (jstring) (*env)->NewStringUTF(env, newPrefix)) == NULL) {
        fprintf(out, "%s !ERROR: Cannot NewStringUTF of 'newPrefix'\n", thisPrefix);
        OMI_Model_release(env, jmodel, model, OMI_ABORT);
        return;
      };
      // pprint the ith SOS2 weight block
      OMI_Block_oprint(env, model->jWeightsSOS2[i], jnewPrefix, out); 
    }
  }
  
  // Dump the double arrays using row and col dimensions
  
  fprintf(out, "%sModel->objective       @ %X =\n  ", prefix, model->objective);
  for (i=0; i<model->max_col_length; i++) fprintf(out, "%g ", model->objective[i]); 
  fprintf(out, "\n");

  fprintf(out, "%sModel->row_lower_bound @ %X =\n  ",  prefix, model->row_lower_bound);
  for (i=0; i<model->max_row_length; i++) fprintf(out, "%g ", model->row_lower_bound[i]); 
  fprintf(out, "\n");

  fprintf(out, "%sModel->row_upper_bound @ %X =\n  ",  prefix, model->row_upper_bound);
  for (i=0; i<model->max_row_length; i++) fprintf(out, "%g ", model->row_upper_bound[i]); 
  fprintf(out, "\n");

  fprintf(out, "%sModel->col_lower_bound @ %X =\n  ",  prefix, model->col_lower_bound);
  for (i=0; i<model->max_col_length; i++) fprintf(out, "%g ", model->col_lower_bound[i]); 
  fprintf(out, "\n");

  fprintf(out, "%sModel->col_upper_bound @ %X =\n  ",  prefix, model->col_upper_bound);
  for (i=0; i<model->max_col_length; i++) fprintf(out, "%g ", model->col_upper_bound[i]); 
  fprintf(out, "\n");
  
  fprintf(out, "%sModel->is_integer @ %X =\n  ",  prefix, model->is_integer);
  for (i=0; i<model->max_col_length; i++) fprintf(out, "%i ", model->is_integer[i]); 
  fprintf(out, "\n");
  
  // Print I/O and solution scalars.
  //
  fprintf(out, "%sModel->ispecs             @ %X = %d\n", 
    prefix, &model->ispecs, model->ispecs);
  fprintf(out, "%sModel->iprint             @ %X = %d\n", 
    prefix, &model->iprint, model->iprint);
  fprintf(out, "%sModel->isumm              @ %X = %d\n", 
    prefix, &model->isumm, model->isumm);
  fprintf(out, "%sModel->nout               @ %X = %d\n", 
    prefix, &model->nout, model->nout);
  fprintf(out, "%sModel->nprob              @ %X = %d\n", 
    prefix, &model->nprob, model->nprob);
  fprintf(out, "%sModel->iterations_limit   @ %X = %d\n", 
    prefix, &model->iterations_limit, model->iterations_limit);
  fprintf(out, "%sModel->number_superbasics @ %X = %d\n", 
    prefix, &model->number_superbasics, model->number_superbasics);
  fprintf(out, "%sModel->number_infeasible  @ %X = %d\n", 
    prefix, &model->number_infeasible, model->number_infeasible);
  fprintf(out, "%sModel->sum_infeasible     @ %X = %g\n", 
    prefix, &model->sum_infeasible, model->sum_infeasible);
  fprintf(out, "%sModel->objective_value    @ %X = %g\n", 
    prefix, &model->objective_value, model->objective_value);

  // Print solution arrays.
  //
  fprintf(out, "%sModel->basis           @ %X =\n  ", prefix, model->basis);
  for (i=0; i<model->max_row_length+model->max_col_length; i++) 
    fprintf(out, "%d ", model->basis[i]); 
  fprintf(out, "\n");

  fprintf(out, "%sModel->primal_solution @ %X =\n  ", prefix, model->primal_solution);
  for (i=0; i<model->max_col_length; i++) fprintf(out, "%g ", model->primal_solution[i]); 
  fprintf(out, "\n");

  fprintf(out, "%sModel->dual_solution @ %X =\n  ", prefix, model->dual_solution);
  for (i=0; i<model->max_row_length; i++) fprintf(out, "%g ", model->dual_solution[i]); 
  fprintf(out, "\n");

  fprintf(out, "%sModel->reduced_costs @ %X =\n  ", prefix, model->reduced_costs);
  for (i=0; i<model->max_col_length; i++) fprintf(out, "%g ", model->reduced_costs[i]); 
  fprintf(out, "\n");

  fprintf(out, "%sModel->workspace_length   @ %X = %d\n", 
    prefix, &model->workspace_length, model->workspace_length);
  fprintf(out, "Model->workspace @ %X =\n  ", prefix, model->workspace);
  if (0) for (i=0; 
       i<(model->workspace_length < 20) ? model->workspace_length : 20; 
       i++
      ) fprintf(out, "%X ", model->workspace[i]); 

  fprintf(out, "\n\n");
  fflush(out);

  // release booty gotten
  (*env)->ReleaseStringUTFChars(env, jprefix, prefix);
  OMI_Model_release(env, jmodel, model, OMI_ABORT);
  
  if (VERBOSE>0) printf("%s Model->pprint('%s') return\n", thisPrefix, prefix);
  return;
}  // OMI_Model_oprint
  
/**
  *  JNI method pprint(). Opens file and calls OMI_Model_oprint()
  *
  */
JNIEXPORT void JNICALL Java_com_epri_omi_Model_pprint 
  (JNIEnv *env, jobject jmodel, jstring jprefix)
{
  const char prefix[] = "Model.c:Java_com_epri_omi_Model_pprint(cbc): ";
//  const int VERBOSE = 0;
  FILE *out;
  time_t  t = time(NULL);
  int i;
  if (0) return;
  if (VERBOSE>0) printf("%sbegin\n",prefix);
  
  if ((out = fopen("Java_Native_Output", "a")) == NULL) {
    printf("%sERROR: Cannot fopen Java_Native_Output.\n",prefix);
    //fflush(stdout);
    return;
  };
  for (i=0; i<1; i++) {
    t = time(NULL);
    fprintf(out, "\nModel->pprint('') %s", asctime(localtime(&t)));
    OMI_Model_oprint(env, jmodel, jprefix, out);
  };
  fclose(out);
  
  if (VERBOSE>0) printf("%sreturn\n",prefix);
  return;
}  // Java_com_epri_omi_Model_pprint

  

/**
  *  JNI method pMPS() reads a model from an MPS file.
  *
  */
JNIEXPORT jint JNICALL Java_com_epri_omi_Model_pMPS
  (JNIEnv *env, jobject jmodel, jstring jfilename)
{ 
  const char prefix[] = "Model.c:Java_com_epri_omi_Model_pMPS(cbc): ";
//  const int VERBOSE = 1;
  
  ClassOMI_Model  *omi_model;
  Cbc_Model     *cbc_model;
  jint status = 0;
  const jint error_status = 999;

  integer    i,j;
  const char *filename;
  FILE       *out;
  time_t     t;

  jint return_code;
    
  // model dimensions
  integer num_rows, num_cols, num_coefs;

  // Open output file for append.
  if ((out = fopen("Java_Native_Output", "a")) == NULL) {
    printf("%s!ERROR, Cannot fopen Java_Native_Output.\n",prefix);
    //fflush(stdout);
    return error_status;
  };
  
  // Print entry time
  t = time(NULL);
  fprintf(out, "%sstarted %s", prefix, asctime(localtime(&t)));
  
  // Check arguments.
  if (env == NULL) {
    fprintf(out, "%s!ERROR, NULL 'env' passed as argument.\n", prefix);
    //fflush(out);
    return error_status;
  };  
  if (jmodel == NULL) {
    fprintf(out, "%s!ERROR, NULL 'jmodel' passed as argument.\n", prefix);
    //fflush(out);
    return error_status;
  };
  if (jfilename == NULL) {
    fprintf(out, "%s!ERROR, NULL 'jfilename' passed as argument.\n", prefix);
    //fflush(out);
    return error_status;
  };
  
  // Announce pMPS
  // In JDK 1.1 get a local reference to a string.
  // The UTF string is like a C string.
  // The alternative is a unicode string.
  if ((filename = (*env)->GetStringUTFChars(env, jfilename, 0)) == NULL) {
    fprintf(out, "%s!ERROR, Cannot GetStringUTFChars of 'jfilename'.\n");
    //fflush(out);
    return error_status;
  };
  fprintf(out, "%s('%s') Arguments OK.\n", prefix, filename);
  if (VERBOSE>0) printf("%s('%s') begin\n", prefix, filename);
  
  // get a C version of the model
  //
  if ((omi_model = OMI_Model_get(env, jmodel)) == NULL) {
    fprintf(out, "Model.c:Java_com_epri_omi_Model_pManne "\
      "!ERROR, Cannot OMI_Model_get of 'jmodel'.\n");
    //fflush(out);
    return error_status;
  };
  // save state
  omi_state.env    = env;
  omi_state.jmodel = jmodel;
  omi_state.model  = omi_model;
  omi_state.out    = out;

  if (omi_state.env    == NULL) fprintf(out,"pMPS(): env is NULL\n");
  if (omi_state.jmodel == NULL) fprintf(out,"pMPS(): jmodel is NULL\n");
  if (omi_state.out    == NULL) fprintf(out,"pMPS(): out is NULL\n");
  fprintf(out,"%s &env = %X, &jmodel = %X, &out = %X\n",prefix,&env,&jmodel,&out);
  fprintf(out,"%s  env = %X,  jmodel = %X,  out = %X\n",prefix,env,jmodel,out);
  fprintf(out,"%s &omi_state.env = %X, &omi_state.jmodel = %X, &omi_state.out = %X\n",
    prefix,&omi_state.env,&omi_state.jmodel,&omi_state.out);
  
  // Assume only one jacobian block
  if (omi_model->jacobian_length != 1) {
    fprintf(out,"%s !ERROR, omi_model->jacobian_length (= %d) != 1.\n", 
      prefix, omi_model->jacobian_length);
    OMI_Model_release(env, jmodel, omi_model, OMI_ABORT);
    //fflush(out);
    return error_status;
  };
  // Assume only one block in weightsSOS1
  if (omi_model->weightsSOS1_length != 1) {
    fprintf(out, "%s !ERROR, omi_model->weightsSOS1_length (= %d) != 1.\n", 
      prefix, omi_model->weightsSOS1_length);
    OMI_Model_release(env, jmodel, omi_model, OMI_ABORT);
    //fflush(out);
    return error_status;
  };
  // Assume only one block in weightsSOS2
  if (omi_model->weightsSOS2_length != 1) {
    fprintf(out, "%s !ERROR, omi_model->weightsSOS2_length (= %d) != 1.\n", 
      prefix, omi_model->weightsSOS2_length);
    OMI_Model_release(env, jmodel, omi_model, OMI_ABORT);
    //fflush(out);
    return error_status;
  };

  // ------------------------------------------------------------------ 
  {
    // Get default cbc_model 
    cbc_model = Cbc_newModel();
//    const char   filename[] = "Mps/Sample/p0033.mps\0";

    // Read an MPS file. 
    {
      // Cbc_readMPS(Cbc_Model * cbc_model,const char *filename, int keepNames, int ignoreErrors);
      status = Cbc_readMps(cbc_model, filename);
//      Cbc_setMaximumIterations    (cbc_model, -3840);
      Cbc_setOptimizationDirection(cbc_model, 1); // (1 - minimize, -1 - maximize, 0 - ignore)
      
      fprintf(out,"%s Return from Cbc_readMps, status = %i\n", prefix, status);
      if (status) {
        fprintf(stderr,"Bad readMps '%s'\n",filename);
        fprintf(out,"Bad readMps '%s'\n",filename);
        Cbc_deleteModel(cbc_model);
        return status;
      }
    }
    // initialSolve may be necessary in Cbc
    if (1) {
      printf("%s calling Cbc_scaling\n",prefix);
      Cbc_scaling(cbc_model,1);
//      Cbc_printModel(cbc_model, prefix);
      printf("%s calling Cbc_initialSolve\n",prefix);
      status = Cbc_initialSolve(cbc_model);
      printf("%s Finished Cbc_initialSolve, status = %i\n",prefix,status);
    }
    // solve it for fun
    if (0) {
      status = Cbc_branchAndBound(cbc_model);
      printf("%s Finished Cbc_branchAndBound, status = %i\n",prefix,status);
      Cbc_printSolution(cbc_model);
    }
    // Save the cbc_model to the omi_model
    fprintf(out,"%s Save the cbc_model to the omi_model\n",prefix);
    //fflush(out);
    {
      void * space;
      const int    * clp_indices, * clp_vectorLengths;
      CoinBigIndex * clp_vectorStarts;
      const double * clp_elements;
                         
      int row_length = Cbc_getNumRows(cbc_model);
      int col_length = Cbc_getNumCols(cbc_model);
      int rim_length = row_length + col_length;
      int elem_length = Cbc_getNumElements(cbc_model);
      
      if (VERBOSE>1) {
        fprintf(stdout,"%s row_length = %i\n", prefix, row_length);
        fprintf(stdout,"%s col_length = %i\n", prefix, col_length);
        fprintf(stdout,"%s rim_length = %i\n", prefix, rim_length);
        fprintf(stdout,"%s elem_length = %i\n", prefix, elem_length);
        //fflush(out);
      }
      
      // check that omi_model object is large enough
      if (omi_model->max_row_length < row_length) {
        fprintf(out,"%s ERROR: omi_model->max_row_length = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->max_row_length, row_length);
        printf("%s ERROR: omi_model->max_row_length = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->max_row_length, row_length);
        fprintf(stdout,"%s row_length = %i\n", prefix, row_length);
        OMI_Model_release(env, jmodel, omi_model, OMI_ABORT);
        //fflush(out);
        return error_status;
      }
      if (omi_model->max_col_length < col_length) {
        fprintf(out,"%s ERROR: omi_model->max_col_length = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->max_col_length, col_length);
        printf("%s ERROR: omi_model->max_col_length = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->max_col_length, col_length);
        OMI_Model_release(env, jmodel, omi_model, OMI_ABORT);
        //fflush(out);
        return error_status;
      }
      /* NOT USING WORKSPACE IN CBC
      if (omi_model->workspace_length < elem_length) {
        fprintf(out,"%s ERROR: omi_model->workspace_length = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->workspace_length, elem_length);
        printf("%s ERROR: omi_model->workspace_length = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->workspace_length, elem_length);
        OMI_Model_release(env, jmodel, omi_model, OMI_ABORT);
        //fflush(out);
        return error_status;
      }
      */
      // check jacobian lengths
      if (omi_model->jacobian[0]->rowIndexLength < row_length) {
        fprintf(out,"%s ERROR: omi_model->jacobian[0]->rowIndexLength = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->jacobian[0]->rowIndexLength, row_length);
        printf("%s ERROR: omi_model->jacobian[0]->rowIndexLength = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->jacobian[0]->rowIndexLength, row_length);
        OMI_Model_release(env, jmodel, omi_model, OMI_ABORT);
        //fflush(out);
        return error_status;
      }
      if (omi_model->jacobian[0]->columnIndexLength < col_length) {
        fprintf(out,"%s ERROR: omi_model->jacobian[0]->columnIndexLength = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->jacobian[0]->columnIndexLength, col_length);
        printf("%s ERROR: omi_model->jacobian[0]->columnIndexLength = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->jacobian[0]->columnIndexLength, col_length);
        OMI_Model_release(env, jmodel, omi_model, OMI_ABORT);
        //fflush(out);
        return error_status;
      }
      if (omi_model->jacobian[0]->valueLength < elem_length) {
        fprintf(out,"%s ERROR: omi_model->jacobian[0]->valueLength = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->jacobian[0]->valueLength, elem_length);
        printf("%s ERROR: omi_model->jacobian[0]->valueLength = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->jacobian[0]->valueLength, elem_length);
        OMI_Model_release(env, jmodel, omi_model, OMI_ABORT);
        //fflush(out);
        return error_status;
      }
      // check SOS1 weight lengths
      if (omi_model->weightsSOS1[0]->rowIndexLength < row_length) {
        fprintf(out,"%s ERROR: omi_model->weightsSOS1[0]->rowIndexLength = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->weightsSOS1[0]->rowIndexLength, row_length);
        printf("%s ERROR: omi_model->weightsSOS1[0]->rowIndexLength = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->weightsSOS1[0]->rowIndexLength, row_length);
        OMI_Model_release(env, jmodel, omi_model, OMI_ABORT);
        //fflush(out);
        return error_status;
      }
      if (omi_model->weightsSOS1[0]->columnIndexLength < col_length) {
        fprintf(out,"%s ERROR: omi_model->weightsSOS1[0]->columnIndexLength = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->weightsSOS1[0]->columnIndexLength, col_length);
        printf("%s ERROR: omi_model->weightsSOS1[0]->columnIndexLength = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->weightsSOS1[0]->columnIndexLength, col_length);
        OMI_Model_release(env, jmodel, omi_model, OMI_ABORT);
        //fflush(out);
        return error_status;
      }
      if (omi_model->weightsSOS1[0]->valueLength < elem_length) {
        fprintf(out,"%s ERROR: omi_model->weightsSOS1[0]->valueLength = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->weightsSOS1[0]->valueLength, elem_length);
        printf("%s ERROR: omi_model->weightsSOS1[0]->valueLength = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->weightsSOS1[0]->valueLength, elem_length);
        OMI_Model_release(env, jmodel, omi_model, OMI_ABORT);
        //fflush(out);
        return error_status;
      }
      // check SOS2 weight lengths
      if (omi_model->weightsSOS2[0]->rowIndexLength < row_length) {
        fprintf(out,"%s ERROR: omi_model->weightsSOS2[0]->rowIndexLength = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->weightsSOS2[0]->rowIndexLength, row_length);
        printf("%s ERROR: omi_model->weightsSOS2[0]->rowIndexLength = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->weightsSOS2[0]->rowIndexLength, row_length);
        OMI_Model_release(env, jmodel, omi_model, OMI_ABORT);
        //fflush(out);
        return error_status;
      }
      if (omi_model->weightsSOS2[0]->columnIndexLength < col_length) {
        fprintf(out,"%s ERROR: omi_model->weightsSOS2[0]->columnIndexLength = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->weightsSOS2[0]->columnIndexLength, col_length);
        printf("%s ERROR: omi_model->weightsSOS2[0]->columnIndexLength = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->weightsSOS2[0]->columnIndexLength, col_length);
        OMI_Model_release(env, jmodel, omi_model, OMI_ABORT);
        //fflush(out);
        return error_status;
      }
      if (omi_model->weightsSOS2[0]->valueLength < elem_length) {
        fprintf(out,"%s ERROR: omi_model->weightsSOS2[0]->valueLength = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->weightsSOS2[0]->valueLength, elem_length);
        printf("%s ERROR: omi_model->weightsSOS2[0]->valueLength = %i is too small. "\
          "Needs to be at least %i.\n", prefix,omi_model->weightsSOS2[0]->valueLength, elem_length);
        OMI_Model_release(env, jmodel, omi_model, OMI_ABORT);
        //fflush(out);
        return error_status;
      }
      
      // omi_model name
      {
        char model_name[OMI_MaxStrLen];
     
        int name_length = omi_model->model_names_length;
        if (VERBOSE>1) {
          printf("%s omi_model->model_names_length = %i\n",prefix, omi_model->model_names_length);
          printf("%s name_length = '%i',  OMI_MaxStrLen = '%i'\n", 
            prefix, name_length, OMI_MaxStrLen);
        }
        if (name_length < 0) {
          fprintf(out,"%s ERROR: omi_model->model_names_length = %i is too bad.\n",
            prefix, name_length);
          OMI_Model_release(env, jmodel, omi_model, OMI_ABORT);
          //fflush(out);
          return error_status;
        }
          
        if (name_length < OMI_MaxStrLen) {
          Cbc_problemName(cbc_model, name_length, model_name);
        } else {
          fprintf(out,"%s ERROR: omi_model->model_names_length = %i is too long. "\
            "Needs to be no more than %i.\n", prefix,name_length, OMI_MaxStrLen);
          OMI_Model_release(env, jmodel, omi_model, OMI_ABORT);
          //fflush(out);
          return error_status;
        }
//tbd        omi_model->identifier = model_name;
        if (VERBOSE>1) 
          printf("%s model_name = '%s'\n", prefix, model_name);
        if (model_name != NULL) 
          fprintf(out,"%s model_name = '%s'\n", prefix, model_name);
        
        // convert model_name to a jchar *
//bobe 06-02-15 this is printing junk
        if (0) { 
          jboolean is_Copy = NULL;
          // convert model_name to jstring using NewStringUTF()
          jstring jmodel_name, jmodel_name_save;
          const jchar *jmodel_name_copy;
          // NOTE: I do not think that jmodel_name needs to be released.
          if ((jmodel_name = (jstring) (*env)->NewStringUTF(env, model_name)) == NULL) {
            fprintf(out, "%s !ERROR: Cannot NewStringUTF of 'model_name'\n",prefix);
            OMI_Model_release(env, jmodel, omi_model, OMI_ABORT);
            return error_status;
          };
          
          // then get jchar * using GetStringChars()
          jmodel_name_copy = (*env)->GetStringChars(env, jmodel_name, &is_Copy);
          // copy over to omi_model->model_names and then dispose
          jccopy(jmodel_name_copy, omi_model->model_names, 0, name_length);
          if (VERBOSE>1) printf("%s ReleaseStringChars(jmodel_name, jmodel_name_copy)\n", prefix);
          (*env)->ReleaseStringChars(env, jmodel_name, jmodel_name_copy);
        }   
        if (VERBOSE>1) 
          printf("%s omi_model->model_names_length = %i\n", prefix, omi_model->model_names_length);
          for (i=0; i<omi_model->model_names_length/8; i++) {
            for (j=0; j<8; j++) fprintf(stdout, "%c", omi_model->model_names[i*8+j]);
          fprintf(stdout, ", ");
        };
        fprintf(stdout, "\n");
      }
      
      omi_model->nonlin_constr_length = 0;
      omi_model->nonlin_obj_length = 0;
      omi_model->nonlin_jac_length = 0;
      omi_model->objective_row = 0;
      
      // objective direction
      //fflush(out);
      omi_model->objective_direction = Cbc_optimizationDirection(cbc_model); // (-1 - minimize, 1 - maximize, 0 - ignore)
      omi_model->objective_constant  = Cbc_objectiveOffset(cbc_model);
      printf("%s Cbc_optimizationDirection = %g, Cbc_objectiveOffset = %g\n",
        prefix, Cbc_optimizationDirection(cbc_model), Cbc_objectiveOffset(cbc_model));
      printf("%s omi_model->objective_direction = %i, omi_model->objective_constant = %g\n",
        prefix, omi_model->objective_direction, omi_model->objective_constant);
      
      omi_model->rim_names_length = 0;
//      omi_model->rim_names1 = ? tbd 
//      use Cbc_lengthNames(cbc_model), 
//          Cbc_rowName(cbc_model, row, name)
//          Cbc_colName(cbc_model, col, name)

      // model dimensions
      fprintf(out,"%s model dimensions\n", prefix);
      //fflush(out);
      omi_model->row_length     = row_length;
      omi_model->col_length     = col_length;

      // model matrix
      fprintf(out,"%s model matrix\n", prefix);      
      //fflush(out);
      // we already checked that (omi_model->jacobian_length == 1)
      omi_model->jacobian[0]->numRows     = row_length;
      omi_model->jacobian[0]->numColumns  = col_length;
      omi_model->jacobian[0]->numElements = elem_length;
      omi_model->weightsSOS1[0]->numRows     = row_length;
      omi_model->weightsSOS1[0]->numColumns  = col_length;
      omi_model->weightsSOS1[0]->numElements = 0;  // HACK: assume MPS cannot get SOS
      omi_model->weightsSOS2[0]->numRows     = row_length;
      omi_model->weightsSOS2[0]->numColumns  = col_length;
      omi_model->weightsSOS2[0]->numElements = 0;  // HACK: assume MPS cannot get SOS

      // copy over the jacobian
      if (omi_model->jacobian[0]->columnIndexLength < col_length+1) {
        printf("%s ERROR: omi_model->jacobian[0]->columnIndexLength (%i) < col_length+1 (%i)\n",
          omi_model->jacobian[0]->columnIndexLength, col_length+1);
        return error_status;
      }
      starts_2_colIndex(Cbc_getVectorStarts(cbc_model), omi_model->jacobian[0]->colIndex, 0, col_length+1);      
      index_2_rowIndex(Cbc_getIndices(cbc_model), omi_model->jacobian[0]->rowIndex, 0, elem_length);      
      dcopy(Cbc_getElements(cbc_model), omi_model->jacobian[0]->value, 0, elem_length);

      // copy SOS1 weights from Cbc to OMI
      if (omi_model->weightsSOS1[0]->columnIndexLength < col_length+1) {
        printf("%s ERROR: omi_model->weightsSOS1[0]->columnIndexLength (%i) < col_length+1 (%i)\n",
          omi_model->weightsSOS1[0]->columnIndexLength, col_length+1);
        return error_status;
      }
/* TBD These may not need to be implemented.
      starts_2_colIndex(Cbc_getWeightStarts(cbc_model), omi_model->weightsSOS1[0]->colIndex, 0, col_length+1);      
      index_2_rowIndex(Cbc_getWeightIndices(cbc_model), omi_model->weightsSOS1[0]->rowIndex, 0, elem_length);      
      dcopy(Cbc_getWeightElements(cbc_model), omi_model->weightsSOS1[0]->value, 0, elem_length);
*/

      // copy SOS2 weights from Cbc to OMI
      if (omi_model->weightsSOS2[0]->columnIndexLength < col_length+1) {
        printf("%s ERROR: omi_model->weightsSOS2[0]->columnIndexLength (%i) < col_length+1 (%i)\n",
          omi_model->weightsSOS2[0]->columnIndexLength, col_length+1);
        return error_status;
      }
/* TBD These may not need to be implemented.
      starts_2_colIndex(Cbc_getWeightStarts(cbc_model), omi_model->weightsSOS2[0]->colIndex, 0, col_length+1);      
      index_2_rowIndex(Cbc_getWeightIndices(cbc_model), omi_model->weightsSOS2[0]->rowIndex, 0, elem_length);      
      dcopy(Cbc_getWeightElements(cbc_model), omi_model->weightsSOS2[0]->value, 0, elem_length);
*/

      // note also cbc_model->clpMatrix()->getVectorLengths();

      // objective and bounds
      fprintf(out,"%s objective and bounds\n",prefix);
      //fflush(out);
      dcopy(Cbc_getObjCoefficients(cbc_model),   omi_model->objective, 0, rim_length);
      dcopy(Cbc_getRowLower(cbc_model),    omi_model->row_lower_bound, 0, rim_length);
      dcopy(Cbc_getRowUpper(cbc_model),    omi_model->row_upper_bound, 0, rim_length);
      dcopy(Cbc_getColLower(cbc_model), omi_model->col_lower_bound, 0, rim_length);
      dcopy(Cbc_getColUpper(cbc_model), omi_model->col_upper_bound, 0, rim_length);
      char_2_long(Cbc_integerInformation(cbc_model), omi_model->is_integer, 0, col_length);
      
      
      // solution status
      fprintf(out,"%s solution status\n",prefix);
      //fflush(out);
      omi_model->iterations_limit   = Cbc_maximumIterations(cbc_model);
      omi_model->number_superbasics = 1; // always 1 for LP
      omi_model->number_infeasible  = Cbc_numberPrimalInfeasibilities(cbc_model);
      omi_model->sum_infeasible     = Cbc_sumPrimalInfeasibilities(cbc_model);
      omi_model->objective_value    = Cbc_objectiveValue(cbc_model);

      // solution
      fprintf(out,"%s solution\n", prefix);
      //fflush(out);
      
      // Cbc column status values
      // isFree = 0x00, basic = 0x01, atUpperBound = 0x02, atLowerBound = 0x03,
      //   superBasic = 0x04, isFixed = 0x05 
      // NOTE: That Cbc is row then column, whereas OMI convention is column then row
/*tbd set basis in Cbc
      fprintf(out,"%s basis\n", prefix); 
      if (0) {
        unsigned char * statusArray = Cbc_statusArray(cbc_model);
        status_2_basis(&statusArray[row_length], omi_model->basis, 0, col_length);
        status_2_basis(statusArray, &(omi_model->basis[col_length]), 0, row_length);
        for (i=0; i<rim_length; i++) 
          if (VERBOSE>10) printf("%s after status_2_basis, omi_model->basis[%i] = %i, statusArray[%i] = %i\n",
            prefix, i, omi_model->basis[i], i, statusArray[i]);
      }
*/      
      fprintf(out,"%s primal solution\n",prefix);
      dcopy(Cbc_getColSolution(cbc_model),omi_model->primal_solution, 0, col_length);
      dcopy(Cbc_getRowActivity(cbc_model),omi_model->primal_solution, col_length, rim_length);
      
      fprintf(out,"%s dual solution\n",prefix);
      dcopy(Cbc_getRowPrice(cbc_model),    omi_model->dual_solution, 0, row_length);
      dcopy(Cbc_getReducedCost(cbc_model), omi_model->reduced_costs, 0, col_length);
    }
    // quit here for debugging
    if (0) {
      // release booty gotten
      (*env)->ReleaseStringUTFChars(env, jfilename, filename);
      OMI_Model_release(env, jmodel, omi_model, NULL);
      // close output file
//      fflush(out);
      return error_status;
    }

    if (status >= 1) {
      fprintf(out,"%s Error while reading MPS file\n",prefix);
      printf("%s Error while reading MPS file\n",prefix);
      OMI_Model_release(env, jmodel, omi_model, OMI_ABORT);
      //fflush(out);
      return error_status;
    }
      
    // bobe debug
    fprintf(out,"%s * Finished reading MPS file\n",prefix);
    fprintf(stdout,"%s * Finished reading MPS file\n", prefix);
    fflush(out); fflush(stdout);

    // release booty gotten
    (*env)->ReleaseStringUTFChars(env, jfilename, filename);
    OMI_Model_release(env, jmodel, omi_model, NULL);
    
    // delete the Cbc_model
    fprintf(out,"%s delete the Cbc_model\n",prefix);
    fprintf(stdout,"%s delete the Cbc_model\n",prefix);
    fflush(out); fflush(stdout);
    Cbc_deleteModel(cbc_model);
  }

  // close output file
  fprintf (out,"%s return, closing output file\n",prefix);
  fprintf (stdout,"%s return, closing output file\n",prefix);
  fflush(out); fflush(stdout);
  fclose(out);

  if (VERBOSE>0) printf("%s('%s') return %i\n", prefix, filename, status);
  return status;
}  // Java_com_epri_omi_Model_pMPS

/**
  *  JNI method pestimate().  Estimates memory requirement for storing model data.
  *  Not used in Cbc.
  *
  */
JNIEXPORT jint JNICALL Java_com_epri_omi_Model_pestimate
  (JNIEnv *env, jobject jmodel, jstring jsolver_name)
{ 
  return (jint)NULL;
}  // Java_com_epri_omi_Model_pworkspace_estimate

/**
  *  JNI method poption() for setting string options.
  *  Not yet used in Cbc.
  *
  */
JNIEXPORT jint JNICALL Java_com_epri_omi_Model_poption
  (JNIEnv *env, jobject jmodel, jstring joption_name)
{ 
  const char *option_name;
  extern int miopt_(char *, integer *, integer *, integer *, ftnlen);
  FILE       *out;
  time_t     t;
  
  char* prefix = "Model.c:Java_com_epri_omi_Model_poption(cbc): ";
  const jint error_status = 999;
  jint return_code;
  ClassOMI_Model  *model;
  static integer ispecs, inform; 
  static ftnlen  len;
  char   temp[OMI_MaxStrLen];

  // Open output file for append.
  if ((out = fopen("Java_Native_Output", "a")) == NULL) {
    printf("%s\n\t!ERROR, Cannot fopen Java_Native_Output.\n",prefix);
    //fflush(stdout);
    return error_status;
  };
  // print entry time
  t = time(NULL);
  fprintf(out, "\n%s  %s", prefix, asctime(localtime(&t)));

  // Check arguments.
  if (env == NULL) {
    fprintf(out, "%s\n\t!ERROR, NULL 'env' passed as argument.\n",prefix);
    //fflush(out);
    return error_status;
  };  
  if (jmodel == NULL) {
    fprintf(out, "%s\n\t!ERROR, NULL 'jmodel' passed as argument.\n",prefix);
    //fflush(out);
    return error_status;
  };
  if (joption_name == NULL) {
    fprintf(out, "%s\n\t!ERROR, NULL 'joption_name' passed as argument.\n",prefix);
    //fflush(out);
    return error_status;
  };
  
  // announce poption
  if ((option_name = (*env)->GetStringUTFChars(env, joption_name, 0)) == NULL) {
    fprintf(out, "%s\n\t!ERROR, Cannot GetStringUTFChars of 'joption_name'.\n",prefix);
    //fflush(out);
    return error_status;
  };
  fprintf(out, "%s option_name = '%s'\n", prefix, option_name);

  // get a C version of the model
  //
  if ((model = OMI_Model_get(env, jmodel)) == NULL) {
    fprintf(out, "%s\n\t!ERROR, Cannot OMI_Model_get of 'jmodel'.\n",prefix);
    return_code = inform;
    return return_code;
  };	
  
  // actually set the option
  len = (ftnlen) sprintf(temp,"%s",option_name);
//tbd  miopt_(temp, &value, &m1file_1.iprint, &m1file_1.isumm, &inform, len);
  return_code = inform;
  
  // release booty gotten
  (*env)->ReleaseStringUTFChars(env, joption_name, option_name);
  OMI_Model_release(env, jmodel, model, NULL);

  // close output file
  fprintf(out, "%s, return %i, closing output file",prefix,return_code);
  fclose(out);
  return return_code;
} // Java_com_epri_omi_Model_poption

/**
  *  JNI method poptioni() for setting integer options.
  *  Not yet used in Cbc.
  *
  */
JNIEXPORT jint JNICALL Java_com_epri_omi_Model_poptioni
  (JNIEnv *env, jobject jmodel, jstring joption_name, jint value)
{ 
  const char *option_name;
  extern int miopti_(
    char *, integer *, integer *, integer *, integer *, ftnlen 
    );
  FILE       *out;
  time_t     t;
  
  const jint error_status = 999;
  jint return_code;
  ClassOMI_Model  *model;
  static integer ispecs, inform;
  static ftnlen len;
  char temp[OMI_MaxStrLen];
  char* prefix = "Model.c:Java_com_epri_omi_Model_poptioni(cbc): ";

  // Open output file for append.
  if ((out = fopen("Java_Native_Output", "a")) == NULL) {
    printf("%s\n\t!ERROR, Cannot fopen Java_Native_Output.\n",prefix);
    //fflush(stdout);
    return error_status;
  };
  // print entry time
  t = time(NULL);
  fprintf(out, "\n%s  %s", prefix, asctime(localtime(&t)));

  // Check arguments.
  if (env == NULL) {
    fprintf(out, "%s\n\t!ERROR, NULL 'env' passed as argument.\n",prefix);
    //fflush(out);
    return error_status;
  };  
  if (jmodel == NULL) {
    fprintf(out, "%s\n\t!ERROR, NULL 'jmodel' passed as argument.\n",prefix);
    //fflush(out);
    return error_status;
  };
  if (joption_name == NULL) {
    fprintf(out, "%s\n\t!ERROR, NULL 'joption_name' passed as argument.\n",prefix);
    //fflush(out);
    return error_status;
  };
  
  // announce poptioni
  if ((option_name = (*env)->GetStringUTFChars(env, joption_name, 0)) == NULL) {
    fprintf(out, "%s\n\t!ERROR, Cannot GetStringUTFChars of 'joption_name'.\n",prefix);
    //fflush(out);
    return error_status;
  };
  fprintf(out, "%s option_name = '%s', value = %i\n", prefix, option_name, value);

  // get a C version of the model
  //
  if ((model = OMI_Model_get(env, jmodel)) == NULL) {
    fprintf(out, "%s\n\t!ERROR, Cannot OMI_Model_get of 'jmodel'.\n",prefix);
    return_code = inform;
    return return_code;
  };
  
  // actually set the option
  len = sprintf(temp,"%s",option_name);
//tbd  miopti_(temp, &value, &m1file_1.iprint, &m1file_1.isumm, &inform, len);
  return_code = inform;
  
  // release booty gotten
  (*env)->ReleaseStringUTFChars(env, joption_name, option_name);
  OMI_Model_release(env, jmodel, model, NULL);

  // close output file
  fprintf(out, "%s, return %i, closing output file",prefix,return_code);
  fclose(out);
  return return_code;
}  // Java_com_epri_omi_Model_poptioni

/**
  *  JNI method poptiond() for setting double options.
  *  Not yet used in Cbc.
  *
  */
JNIEXPORT jint JNICALL Java_com_epri_omi_Model_poptiond
  (JNIEnv *env, jobject jmodel, jstring joption_name, jdouble value)
{ 
  const char *option_name;
  extern mioptr_(char *, doublereal *, integer *, integer *, integer *, ftnlen);
  FILE       *out;
  time_t     t;
  
  char* prefix = "Model.c:Java_com_epri_omi_Model_poptiond(cbc): ";
  
  const jint error_status = 999;
  jint return_code;
  ClassOMI_Model  *model;
  static integer ispecs, inform;
  static ftnlen  len;
  char   temp[OMI_MaxStrLen];

  // Open output file for append.
  if ((out = fopen("Java_Native_Output", "a")) == NULL) {
    printf("%s\n\t!ERROR, Cannot fopen Java_Native_Output.\n",prefix);
    //fflush(stdout);
    return error_status;
  };
  // print entry time
  t = time(NULL);
  fprintf(out, "\nModel.c:Java_com_epri_omi_Model_poptiond  %s", asctime(localtime(&t)));

  // Check arguments.
  if (env == NULL) {
    fprintf(out, "%s\n\t!ERROR, NULL 'env' passed as argument.\n",prefix);
    //fflush(out);
    return error_status;
  };  
  if (jmodel == NULL) {
    fprintf(out, "%s\n\t!ERROR, NULL 'jmodel' passed as argument.\n",prefix);
    //fflush(out);
    return error_status;
  };
  if (joption_name == NULL) {
    fprintf(out, "%s\n\t!ERROR, NULL 'joption_name' passed as argument.\n",prefix);
    //fflush(out);
    return error_status;
  };
  
  // announce poptiond
  if ((option_name = (*env)->GetStringUTFChars(env, joption_name, 0)) == NULL) {
    fprintf(out, "%s\n\t!ERROR, Cannot GetStringUTFChars of 'joption_name'.\n",prefix);
    //fflush(out);
    return error_status;
  };
  fprintf(out, "%s option_name = '%s', value = %g\n", prefix, option_name, value);

  // get a C version of the model
  //
  if ((model = OMI_Model_get(env, jmodel)) == NULL) {
    fprintf(out, "%s\n\t!ERROR, Cannot OMI_Model_get of 'jmodel'.\n",prefix);
    return_code = inform;
    return return_code;
  };
      
  // actually set the option
  len = sprintf(temp,"%s",option_name);
//tbd  mioptr_(temp, &value, &m1file_1.iprint, &m1file_1.isumm, &inform, len);
  return_code = inform;
  
  // release booty gotten
  (*env)->ReleaseStringUTFChars(env, joption_name, option_name);
  OMI_Model_release(env, jmodel, model, NULL);

  // close output file
  fprintf(out, "%s, return %i, closing output file",prefix,return_code);
  fclose(out);
  return return_code;
}  // Java_com_epri_omi_Model_poptiond

/**
  *  JNI method psetup() for setting up C to take a model.
  *  Not yet used in Cbc.
  *
  */
JNIEXPORT jint JNICALL Java_com_epri_omi_Model_psetup
  (JNIEnv *env, jobject jmodel, jstring jsolver_name)
{
  return (jint)NULL;
}  //  Java_com_epri_omi_Model_psetup()

/**
  *  JNI method psolve() for solving a model with Cbc.
  *  Return 0 if all OK.  
  *  Return positive value if something wrong with solving the problem.
  *  Return 999 on error in this glue.
  *
  */
JNIEXPORT jint JNICALL Java_com_epri_omi_Model_psolve
  (JNIEnv *env, jobject jmodel, jstring jsolver_name)
{ 
  const char prefix[] = "Model.c:Java_com_epri_omi_Model_psolve(cbc): ";
//  const int VERBOSE = 1;
 
  ClassOMI_Model  * omi_model;
  Cbc_Model       * cbc_model;  
  
  // local model properties
  char *names;
  char *names1;
  char *names2;

  integer    i,j;
  const char *solver_name;
  FILE       *out;
  time_t     t;
  
  const jint error_status = 999;
  jint status = 0;
  
  if (VERBOSE>0) printf("%sbegin\n",prefix);
  // Open output file for append.
  if ((out = fopen("Java_Native_Output", "a")) == NULL) {
    printf("%s!ERROR, Cannot fopen Java_Native_Output.\n",prefix);
    //fflush(stdout);
    return error_status;
  };
  
  // Print entry time
  t = time(NULL);
  if (VERBOSE>0) fprintf(out, "%sstarted %s\n", prefix, asctime(localtime(&t)));
  
  // Print version
  if (VERBOSE>0) fprintf(stdout, "%sVersion %g\n", prefix, Cbc_getVersion());
  
  // Check arguments.
  if (env == NULL) {
    fprintf(out, "%s!ERROR, NULL 'env' passed as argument.\n", prefix);
    //fflush(out);
    return error_status;
  };  
  if (jmodel == NULL) {
    fprintf(out, "%s!ERROR, NULL 'jmodel' passed as argument.\n", prefix);
    //fflush(out);
    return error_status;
  };
  if (jsolver_name == NULL) {
    fprintf(out, "%s!ERROR, NULL 'jsolver_name' passed as argument.\n", prefix);
    //fflush(out);
    return error_status;
  };
  
  // Announce psolve
  // In JDK 1.1 get a local reference to a string.
  // The UTF string is like a C string.
  // The alternative is a unicode string.
  if ((solver_name = (*env)->GetStringUTFChars(env, jsolver_name, 0)) == NULL) {
    fprintf(out, "%s!ERROR, Cannot GetStringUTFChars of 'jsolver_name'.\n");
    //fflush(out);
    return error_status;
  };
  if (VERBOSE>0) fprintf(out, "%s('%s') Arguments OK.\n", prefix, solver_name);
  if (VERBOSE>0) fprintf(stdout, "%s('%s') Arguments OK.\n", prefix, solver_name);
  fflush(out); fflush(stdout);
  (*env)->ReleaseStringUTFChars(env, jsolver_name, solver_name);
  
  // get a C version of the omi_model
  if (VERBOSE>1) printf("%s get a C version of the omi_model\n", prefix);
  if ((omi_model = OMI_Model_get(env, jmodel)) == NULL) {
    fprintf(out, "%s!ERROR, Cannot OMI_Model_get of 'jmodel'.\n", prefix);
    fflush(out);
    return error_status;
  };
  
  // save state
  if (VERBOSE>1) printf("%s save state\n", prefix);
  omi_state.env    = env;
  omi_state.jmodel = jmodel;
  omi_state.model  = omi_model;
  omi_state.out    = out;

  // debug prints
  if (VERBOSE>1) printf("%s debug prints\n", prefix);
  if (VERBOSE>1) {
    if (omi_state.env    == NULL) fprintf(out,"psolve(): env is NULL\n");
    if (omi_state.jmodel == NULL) fprintf(out,"psolve(): jmodel is NULL\n");
    if (omi_state.out    == NULL) fprintf(out,"psolve(): out is NULL\n");
    fprintf(out,"%s &env = %X, &jmodel = %X, &out = %X\n", prefix,&env,&jmodel,&out);
    fprintf(out,"%s  env = %X,  jmodel = %X,  out = %X\n", prefix,env,jmodel,out);
    fprintf(out,"%s &omi_state.env = %X, &omi_state.jmodel = %X, &omi_state.out = %X\n",
      prefix,&omi_state.env,&omi_state.jmodel,&omi_state.out);
    fflush(out);
  }
  
  // print the OMI model
  if (0) {
    jstring jprefix;
    if (VERBOSE>1) printf("%s print the omi_model\n", prefix);
    // NOTE: I do not think that jprefix needs to be released or freed.
    if ((jprefix = (jstring) (*env)->NewStringUTF(env, prefix)) == NULL) {
      fprintf(out, "Model.c:Java_com_epri_omi_Model_psolve "\
        "!ERROR: Cannot NewStringUTF of 'prefix'\n");
      OMI_Model_release(env, jmodel, omi_model, OMI_ABORT);
      return error_status;
    }
    OMI_Model_oprint(env, jmodel, jprefix, stdout);
  }
  
  // shift omi_model to cbc data
  if (VERBOSE>1) printf("%s shift omi_model to cbc data\n", prefix);
  
  // set up a default cbc model
  cbc_model = Cbc_newModel();
  if (VERBOSE>1) printf("%s set up a default cbc model\n", prefix);
  
  // set omi_model name
  //   TBD omi_model->model_names needs to be converted from jchar to char
//TBD  Cbc_setProblemName(cbc_model, omi_model->model_names_length, omi_model->model_names);

  // load the omi_model into cbc
  if (VERBOSE>1) printf("%s load the omi_model into cbc\n", prefix);
  {
    int i = 0;     
    const int numcols   = omi_model->jacobian[0]->numColumns;
    const int numrows   = omi_model->jacobian[0]->numRows;
    const int numelem   = omi_model->jacobian[0]->numElements;
    const int numrim    = numcols + numrows;
    
    const int numSOS1rows   = omi_model->weightsSOS1[0]->numRows;
    const int numSOS1elem   = omi_model->weightsSOS1[0]->numElements;
    
    const int numSOS2rows   = omi_model->weightsSOS2[0]->numRows;
    const int numSOS2elem   = omi_model->weightsSOS2[0]->numElements;
    
    // allocate local jacobian
    // column order
    CoinBigIndex  * start0  = (CoinBigIndex *)  malloc((numcols+1) * sizeof(CoinBigIndex ));
    int           * index0  = (int *)           malloc(numelem     * sizeof(int          )); 
    unsigned char * status0 = (unsigned char *) malloc(numrim      * sizeof(unsigned char));

    // allocate local integer information
    char          * is_integer0  = (char *)  malloc((numcols) * sizeof(char));
    
    // allocate local SOS1 weights
    CoinBigIndex  * weightSOS1Start0  = (CoinBigIndex *)  malloc((numrows+1) * sizeof(CoinBigIndex ));
    int           * weightSOS1Index0  = (int *)           malloc(numSOS1elem * sizeof(int          )); 
    
    // allocate local SOS2 weights
    CoinBigIndex  * weightSOS2Start0  = (CoinBigIndex *)  malloc((numrows+1) * sizeof(CoinBigIndex ));
    int           * weightSOS2Index0  = (int *)           malloc(numSOS2elem * sizeof(int          )); 

    if (VERBOSE>5) {
      printf("%s sizeof(CoinBigIndex ) = %i\n", prefix, sizeof(CoinBigIndex ));
      printf("%s sizeof(int          ) = %i\n", prefix, sizeof(int          ));
      printf("%s sizeof(unsigned char) = %i\n", prefix, sizeof(unsigned char));
    }
    
    if (VERBOSE>1) {
      printf("%s numcols = %i, numrows = %i, numelem = %i, numrim = %i\n",
        prefix, numcols, numrows, numelem, numrim);
      printf("%s numSOS1rows,  = %i, numSOS1elem = %i, numSOS2rows,  = %i, numSOS2elem = %i\n",
        prefix, numSOS1rows, numSOS1elem, numSOS2rows, numSOS2elem);
    }
    if (numSOS1elem > 0 && numSOS2elem > 0) {
        printf("%sERROR: This driver does yet not support SOS1 and SOS2 in the same model\n",prefix);
        return error_status;
    }

    // copy over Jacobian column starts and row indices
    index_2_starts(omi_model->jacobian[0]->colIndex, start0, 0, numcols+1);
    start0[numcols] = numelem;  // cbc has to have numcols element of starts set to numelem
    if (VERBOSE>3) {
      printf("%s after index_2_starts start0[%i] = %i\n", prefix, numcols, start0[numcols]); 
      fflush(stdout);
    }
    rowIndex_2_index (omi_model->jacobian[0]->rowIndex, index0, 0, numelem);
    if (VERBOSE>3) {
      printf("%s after rowIndex_2_index start0[%i] = %i\n", prefix, numcols, start0[numcols]); 
      fflush(stdout);
    }
    
    // copy over SOS1 Weights' column starts and row indices
    // row order
    index_2_starts(omi_model->weightsSOS1[0]->rowIndex, weightSOS1Start0, 0, numSOS1rows+1);
    // Cbc_C_Interface requires that numSOS1rows element of starts set to numSOS1elem
    weightSOS1Start0[numSOS1rows] = numSOS1elem;  
    if (VERBOSE>3) {
      printf("%s after index_2_starts weightStart0[%i] = %i\n", prefix, 
        numSOS1rows, weightSOS1Start0[numSOS1rows]); 
      fflush(stdout);
    }
    colIndex_2_index (omi_model->weightsSOS1[0]->colIndex, weightSOS1Index0, 0, numSOS1elem);
    if (VERBOSE>3) {
      printf("%s after colIndex_2_index weightStart0[%i] = %i\n", prefix, 
        numSOS1rows, weightSOS1Start0[numSOS1rows]); 
      fflush(stdout);
    }
    
    // copy over SOS2 Weights' column starts and row indices
    // row order
    index_2_starts(omi_model->weightsSOS2[0]->rowIndex, weightSOS2Start0, 0, numSOS2rows+1);
     // Cbc_C_Interface requires that numSOS2rows element of starts set to numSOS2elem
    weightSOS2Start0[numSOS2rows] = numSOS2elem; 
    if (VERBOSE>1) {
      printf("%s after index_2_starts weightSOS2Start0[%i] = %i\n", prefix, 
        numSOS2rows, weightSOS2Start0[numSOS2rows]); 
      fflush(stdout);
    }
    colIndex_2_index(omi_model->weightsSOS2[0]->colIndex, weightSOS2Index0, 
      0, numSOS2elem);
    if (VERBOSE>1) {
      printf("%s after colIndex_2_index weightSOS2Start0[%i] = %i\n", prefix, 
        numSOS2rows, weightSOS2Start0[numSOS2rows]); 
      fflush(stdout);
    }
    
    // set some omi_model parameters
    if (VERBOSE>1) printf("%s set some omi_model parameters\n", prefix);
    Cbc_setOptimizationDirection(cbc_model, omi_model->objective_direction); // (1 - minimize, -1 - maximize, 0 - ignore)
    Cbc_setObjectiveOffset      (cbc_model, omi_model->objective_constant);
    Cbc_setMaximumIterations    (cbc_model, omi_model->iterations_limit);
    
    {
      // jacobian in column order
      const CoinBigIndex * start = start0;
      const int          * index = index0;
      const double       * value = omi_model->jacobian[0]->value;
      // weightsSOS1 in row order
      const CoinBigIndex * weightSOS1Start = weightSOS1Start0;
      const int          * weightSOS1Index = weightSOS1Index0;
      const double       * weightsSOS1     = omi_model->weightsSOS1[0]->value;
      // weightsSOS2 in row order
      const CoinBigIndex * weightSOS2Start = weightSOS2Start0;
      const int          * weightSOS2Index = weightSOS2Index0;
      const double       * weightsSOS2     = omi_model->weightsSOS2[0]->value;
      // rim
      const double* obj   = omi_model->objective;
      const double* rowlb = omi_model->row_lower_bound;
      const double* rowub = omi_model->row_upper_bound;
      const double* collb = omi_model->col_lower_bound;
      const double* colub = omi_model->col_upper_bound;
  
      Cbc_setLogLevel(cbc_model, 0); // 4 is verbose
      if (VERBOSE>1) printf("%s Log Level %i\n", prefix, Cbc_logLevel(cbc_model));
      if (VERBOSE>1) printf("%s calling Cbc_loadProblem\n", prefix);
      Cbc_loadProblem (cbc_model, numcols, numrows, 
        start, index, value, collb, colub, obj, rowlb, rowub);
      if (VERBOSE>1) printf("%s finished Cbc_loadProblem\n", prefix);
        
      // set integer variables
      long_2_char(omi_model->is_integer, is_integer0, 0, numcols);
      Cbc_copyInIntegerInformation(cbc_model, is_integer0);
      
      // create SOS1 constraints
      if (VERBOSE>1) printf("%s create SOS1 constraints, numSOS1elem = %i\n",prefix,numSOS1elem);
      if (numSOS1elem > 0) {
        int                  row, col;
        int                  numSOS1Cols; // number of SOS columns
        const CoinBigIndex * colIndex;    // column indices
        const double       * colWeight;   // column weights
        assert(numSOS1elem == weightSOS1Start[numSOS1rows]);
        if (VERBOSE>2) {
          printf("%s numSOS1rows = %i\n",prefix,numSOS1rows);
            for (row=0;row<numSOS1rows;row++) {
            printf("%s row = %i\n",prefix,row);
            printf("%s weightSOS1Start[%i] = %i\n",prefix,row,weightSOS1Start[row]);
            printf("%s weightSOS1Start[%i+1] = %i\n",prefix,row,weightSOS1Start[row+1]);
            numSOS1Cols = weightSOS1Start[row+1]-weightSOS1Start[row];
            colIndex    = weightSOS1Index+weightSOS1Start[row];
            colWeight   = weightsSOS1+weightSOS1Start[row];
            printf("%s numSOS1Cols = %i\n",prefix,numSOS1Cols);
            printf("%s colIndex [%i] = %i\n",prefix,0,colIndex[0]);
            printf("%s colIndex [%i] = %i\n",prefix,numSOS1Cols-1,colIndex[numSOS1Cols-1]);
            if (VERBOSE>3 && numSOS1Cols>0) for (col=0;col<numSOS1Cols;col++) {
              printf("%s colIndex [%i] = %i\n",prefix,col,colIndex[col]);
              printf("%s colWeight[%i] = %f\n",prefix,col,colWeight[col]);
            }
            fflush(stdout);
          }
        }
        if (VERBOSE>1) printf("%s calling Cbc_addSOS_Sparse(type = 1)\n",prefix);
        Cbc_addSOS_Sparse(cbc_model, weightSOS1Start, weightSOS1Index, weightsSOS1, 1);
      } else {
        if (VERBOSE>1) printf("%s No SOS1 constraints.\n",prefix);
      }
      
      // create SOS2 constraints
      if (VERBOSE>1) printf("%s create SOS2 constraints, numSOS2elem = %i\n",prefix,numSOS2elem);
      if (numSOS2elem > 0) {
        int                  row, col;
        int                  numSOS2Cols; // number of SOS columns
        const CoinBigIndex * colIndex;    // column indices
        const double       * colWeight;   // column weights
        assert(numSOS2elem == weightSOS2Start[numSOS2rows]);
        if (VERBOSE>2) {
          printf("%s numSOS1rows = %i\n",prefix,numSOS1rows);
          for (row=0;row<numSOS2rows;row++) {
            printf("%s row = %i\n",prefix,row);
            printf("%s weightSOS2Start[%i] = %i\n",prefix,row,weightSOS2Start[row]);
            printf("%s weightSOS2Start[%i+1] = %i\n",prefix,row,weightSOS2Start[row+1]);
            numSOS2Cols = weightSOS2Start[row+1]-weightSOS2Start[row];
            colIndex    = weightSOS2Index+weightSOS2Start[row];
            colWeight   = weightsSOS2+weightSOS2Start[row];
            printf("%s numSOS2Cols = %i\n",prefix,numSOS2Cols);
            printf("%s colIndex [%i] = %i\n",prefix,0,colIndex[0]);
            printf("%s colIndex [%i] = %i\n",prefix,numSOS2Cols-1,colIndex[numSOS2Cols-1]);
            if (VERBOSE>3 && numSOS2Cols>0) for (col=0;col<numSOS2Cols;col++) {
              printf("%s colIndex [%i] = %i\n",prefix,col,colIndex[col]);
              printf("%s colWeight[%i] = %f\n",prefix,col,colWeight[col]);
            }
            fflush(stdout);
          }
        }
        if (VERBOSE>1) printf("%s calling Cbc_addSOS_Sparse(type = 2)\n",prefix);
        Cbc_addSOS_Sparse(cbc_model, weightSOS2Start, weightSOS2Index, weightsSOS2, 2);
      } else {
        if (VERBOSE>1) printf("%s No SOS2 constraints.\n",prefix);
      }
    }
    
    if (VERBOSE> 2) {
      printf("%s omi_model->iterations_limit = %i\n",prefix,omi_model->iterations_limit);
      printf("%s Cbc maximumIterations = %i\n",prefix, Cbc_maximumIterations(cbc_model));
    }
    // starting solution is in omi_model->basis 
    // convert it from long to unsigned char
    if (VERBOSE>1) printf("%s set starting solution from omi_model->basis\n", prefix);
    if (VERBOSE>1) printf("%s numrim = %i, sizeof(unsigned char) = %i\n", 
      prefix, numrim, sizeof(unsigned char));
    if (VERBOSE>1) printf("%s numcols = %i, numrows = %i, numrim = %i\n",
      prefix, numcols, numrows, numrim);
    // NOTE: That CLP is row then column, whereas OMI convention is column then row
    basis_2_status(omi_model->basis, &status0[numrows], 0, numcols);
    basis_2_status(&(omi_model->basis[numcols]), status0, 0, numrows);
    for (i=0; i<numrim; i++) 
      if (VERBOSE>10) printf("%s after basis_2_status, omi_model->basis[%i] = %i, status0[%i] = %i\n",
        prefix, i, omi_model->basis[i], i, status0[i]);
    if (VERBOSE>1) printf("%s calling Cbc_copyinStatus(cbc_model, status0)\n",prefix);
//tbd    Cbc_copyinStatus(cbc_model, status0);
    if (VERBOSE>1) printf("%s finished Cbc_copyinStatus(cbc_model, status0)\n",prefix);

    // print cbc_model
    if (VERBOSE>5) Cbc_printModel(cbc_model, prefix);
    
    // solve the problem
//      char statusMessage[][] = 
//        {"Unknown","Optimal","PrimalInfeasible","DualInfeasible","Stopped", "Errors"};
    if (VERBOSE>1) printf("%s solve the problem\n", prefix);
    
    {
      double tol = 1.0e-5;
      if (VERBOSE>1) printf("%s Calling Cbc_setIntegerTolerance %g\n",prefix,tol);
      Cbc_setIntegerTolerance(cbc_model,tol);
    }
    if (VERBOSE>1) printf("%s Calling Cbc_scaling\n",prefix);
    Cbc_scaling(cbc_model,1);
    if (VERBOSE>1) printf("%s Calling Cbc_initialSolve\n",prefix); fflush(stdout);
    status = Cbc_initialSolve(cbc_model);
    if (VERBOSE>1) printf("%s Calling Cbc_branchAndBound\n",prefix); fflush(stdout);
    status = Cbc_branchAndBound(cbc_model);
    if (VERBOSE>1) printf("%s Finished Cbc_branchAndBound, status = %i\n",prefix,status);

    // print the Cbc solution
    if (VERBOSE>2) Cbc_printSolution(cbc_model);
    
    // assign the cbc solution to the omi_model objects
    fprintf(out,"%s assign the cbc solution to the OMI omi_model objects\n",prefix);
    //fflush(out);
    if (VERBOSE>1) printf("%s assign the cbc solution to the omi_model objects\n", prefix);
    omi_model->iterations_limit    = Cbc_maximumIterations(cbc_model);
    omi_model->number_superbasics  = 1; // always 1 for LP
    omi_model->number_infeasible   = Cbc_numberPrimalInfeasibilities(cbc_model);
    omi_model->sum_infeasible      = Cbc_sumPrimalInfeasibilities(cbc_model);
    omi_model->objective_value     = Cbc_objectiveValue(cbc_model);
    omi_model->objective_direction = Cbc_optimizationDirection(cbc_model);
    if (VERBOSE>1) printf("%sobjective_value = %g\n",prefix,omi_model->objective_value);

      
    // Cbc column status values
    // isFree = 0x00, basic = 0x01, atUpperBound = 0x02, atLowerBound = 0x03,
    //   superBasic = 0x04, isFixed = 0x05 
    // NOTE: That Cbc is row then column, whereas OMI convention is column then row
/*tbd
    fprintf(out,"%s basis\n", prefix); //fflush(out); 
    {
      unsigned char * statusArray = Cbc_statusArray(cbc_model);
      status_2_basis(&statusArray[numrows], omi_model->basis, 0, numcols);
      status_2_basis(statusArray, &(omi_model->basis[numcols]), 0, numrows);
      for (i=0; i<numrim; i++) 
        if (VERBOSE>10) printf("%s after status_2_basis, omi_model->basis[%i] = %i, statusArray[%i] = %i\n",
          prefix, i, omi_model->basis[i], i, statusArray[i]);
    }
 */   
    fprintf(out,"%s primal solution\n",prefix);
    dcopy(Cbc_getColSolution(cbc_model), omi_model->primal_solution, 0, numcols);
    if (0)
      dcopy(Cbc_getRowActivity(cbc_model), omi_model->primal_solution, numcols, numrim);
    
    fprintf(out,"%s dual solution\n",prefix);
    dcopy(Cbc_getRowPrice(cbc_model),    omi_model->dual_solution, 0, numrows);
    if (0)
      dcopy(Cbc_getReducedCost(cbc_model), omi_model->dual_solution, numcols, numrim);
    dcopy(Cbc_getReducedCost(cbc_model), omi_model->reduced_costs, 0, numcols);

    // free up the sparse matrix compenents that we malloc'ed
    if (VERBOSE>1) fprintf(stdout,"%s free up the sparse matrix\n", prefix); fflush(stdout);
    free(start0); free(index0); free(status0); free(is_integer0);
    // SOS weights may have had zero allocation
    if (weightSOS1Start0) free(weightSOS1Start0); 
    if (weightSOS1Index0) free(weightSOS1Index0);
    if (weightSOS2Start0) free(weightSOS2Start0); 
    if (weightSOS2Index0) free(weightSOS2Index0);
  }
     
  if (VERBOSE>1) fprintf(out,"%s Cbc_branchAndBound returns status = %d\n", prefix, status);
  if (VERBOSE>1) fprintf(stdout,"%s Cbc_branchAndBound returns status = %d\n", prefix, status);

  // delete the Cbc_model
  fprintf(out,"%s delete the Cbc_model\n",prefix);
  if (VERBOSE>0) fprintf(stdout,"%s delete the Cbc_model\n",prefix);
  fflush(out); fflush(stdout);
  Cbc_deleteModel(cbc_model);

  // release booty gotten
  if (VERBOSE>1) fprintf(stdout,"%s release booty gotten\n", prefix);  fflush(stdout);
  OMI_Model_release(env, jmodel, omi_model, NULL);

  if (VERBOSE>0) fprintf(out,"%sreturn %d\n", prefix, status);
  if (VERBOSE>0) fprintf(stdout,"%sreturn %d\n", prefix, status);
  fflush(out); fflush(stdout);
  
  // close output file
  if (VERBOSE>0) fprintf (out,"%sreturn, closing output file",prefix);
  fclose(out);
  return status;

}  // Java_com_epri_omi_Model_psolve

/**
  *  JNI method pwriteMPS() for writing an MPS file.
  *
  */
JNIEXPORT jint JNICALL Java_com_epri_omi_Model_pwriteMPS
  (JNIEnv *env, jobject jmodel, jstring jfile_name)
{ 
  const char prefix[] = "Model.c:Java_com_epri_omi_Model_pwriteMPS(cbc): ";
//	const int VERBOSE = 1;
  
  ClassOMI_Model  * omi_model;
  Cbc_Model       * cbc_model;  
  
  // local model properties
  char *names;
  char *names1;
  char *names2;

  integer    i,j;
  const char *file_name;
  FILE       *out;
  time_t     t;
  
  const jint error_status = 999;
  jint status = 0;
  
  // Open output file for append.
  if ((out = fopen("Java_Native_Output", "a")) == NULL) {
    printf("%s!ERROR, Cannot fopen Java_Native_Output.\n",prefix);
    //fflush(stdout);
    return error_status;
  };
  
  // Print entry time
  t = time(NULL);
  if (VERBOSE>0) fprintf(out, "%sstarted %s", prefix, asctime(localtime(&t)));
  
  // Check arguments.
  if (env == NULL) {
    fprintf(out, "%s!ERROR, NULL 'env' passed as argument.\n", prefix);
    //fflush(out);
    return error_status;
  };  
  if (jmodel == NULL) {
    fprintf(out, "%s!ERROR, NULL 'jmodel' passed as argument.\n", prefix);
    //fflush(out);
    return error_status;
  };
  if (jfile_name == NULL) {
    fprintf(out, "%s!ERROR, NULL 'jfile_name' passed as argument.\n", prefix);
    //fflush(out);
    return error_status;
  };
  
  // Announce pwriteMPS
  // In JDK 1.1 get a local reference to a string.
  // The UTF string is like a C string.
  // The alternative is a unicode string.
  if ((file_name = (*env)->GetStringUTFChars(env, jfile_name, 0)) == NULL) {
    fprintf(out, "%s!ERROR, Cannot GetStringUTFChars of 'jfile_name'.\n");
    //fflush(out);
    return error_status;
  };
  if (VERBOSE>0) fprintf(out, "%s('%s') Arguments OK.\n", prefix, file_name);
  if (VERBOSE>0) fprintf(stdout, "%s('%s') Arguments OK.\n", prefix, file_name);
  fflush(out); fflush(stdout);
  
  // get a C version of the omi_model
  if (VERBOSE>1) printf("%s get a C version of the omi_model\n", prefix);
  if ((omi_model = OMI_Model_get(env, jmodel)) == NULL) {
    fprintf(out, "%s!ERROR, Cannot OMI_Model_get of 'jmodel'.\n", prefix);
    fflush(out);
    return error_status;
  };
  
  // save state
  if (VERBOSE>1) printf("%s save state\n", prefix);
  omi_state.env    = env;
  omi_state.jmodel = jmodel;
  omi_state.model  = omi_model;
  omi_state.out    = out;

  // debug prints
  if (VERBOSE>1) printf("%s debug prints\n", prefix);
  if (VERBOSE>1) {
    if (omi_state.env    == NULL) fprintf(out,"psolve(): env is NULL\n");
    if (omi_state.jmodel == NULL) fprintf(out,"psolve(): jmodel is NULL\n");
    if (omi_state.out    == NULL) fprintf(out,"psolve(): out is NULL\n");
    fprintf(out,"%s &env = %X, &jmodel = %X, &out = %X\n", prefix,&env,&jmodel,&out);
    fprintf(out,"%s  env = %X,  jmodel = %X,  out = %X\n", prefix,env,jmodel,out);
    fprintf(out,"%s &omi_state.env = %X, &omi_state.jmodel = %X, &omi_state.out = %X\n",
      prefix,&omi_state.env,&omi_state.jmodel,&omi_state.out);
    fflush(out);
  }
  
  // print the OMI model
  if (0) {
    jstring jprefix;
    if (VERBOSE>1) printf("%s print the omi_model\n", prefix);
    // NOTE: I do not think that jprefix needs to be released or freed.
    if ((jprefix = (jstring) (*env)->NewStringUTF(env, prefix)) == NULL) {
      fprintf(out, "Model.c:Java_com_epri_omi_Model_psolve "\
        "!ERROR: Cannot NewStringUTF of 'prefix'\n");
      OMI_Model_release(env, jmodel, omi_model, OMI_ABORT);
      return error_status;
    }
    OMI_Model_oprint(env, jmodel, jprefix, stdout);
  }
  
  // shift omi_model to cbc data
  if (VERBOSE>1) printf("%s shift omi_model to cbc data\n", prefix);
  
  // set up a default cbc model
  cbc_model = Cbc_newModel();
  if (VERBOSE>1) printf("%s set up a default cbc model\n", prefix);
  
  // set omi_model name
  //   TBD omi_model->model_names needs to be converted from jchar to char
//TBD  Cbc_setProblemName(cbc_model, omi_model->model_names_length, omi_model->model_names);

  // load the omi_model into cbc
  if (VERBOSE>1) printf("%s load the omi_model into cbc\n", prefix);
  {
    int i = 0;     
    const int numcols   = omi_model->jacobian[0]->numColumns;
    const int numrows   = omi_model->jacobian[0]->numRows;
    const int numelem   = omi_model->jacobian[0]->numElements;
    const int numrim    = numcols + numrows;
    
    const int numSOS1rows   = omi_model->weightsSOS1[0]->numRows;
    const int numSOS1elem   = omi_model->weightsSOS1[0]->numElements;
    
    const int numSOS2rows   = omi_model->weightsSOS2[0]->numRows;
    const int numSOS2elem   = omi_model->weightsSOS2[0]->numElements;
    
    // allocate local jacobian
    // column order
    CoinBigIndex  * start0  = (CoinBigIndex *)  malloc((numcols+1) * sizeof(CoinBigIndex ));
    int           * index0  = (int *)           malloc(numelem     * sizeof(int          )); 
    unsigned char * status0 = (unsigned char *) malloc(numrim      * sizeof(unsigned char));

    // allocate local integer information
    char          * is_integer0  = (char *)  malloc((numcols) * sizeof(char));
    
    // allocate local SOS1 weights
    CoinBigIndex  * weightSOS1Start0  = (CoinBigIndex *)  malloc((numrows+1) * sizeof(CoinBigIndex ));
    int           * weightSOS1Index0  = (int *)           malloc(numSOS1elem * sizeof(int          )); 
    
    // allocate local SOS2 weights
    CoinBigIndex  * weightSOS2Start0  = (CoinBigIndex *)  malloc((numrows+1) * sizeof(CoinBigIndex ));
    int           * weightSOS2Index0  = (int *)           malloc(numSOS2elem * sizeof(int          )); 

    if (VERBOSE>5) {
      printf("%s sizeof(CoinBigIndex ) = %i\n", prefix, sizeof(CoinBigIndex ));
      printf("%s sizeof(int          ) = %i\n", prefix, sizeof(int          ));
      printf("%s sizeof(unsigned char) = %i\n", prefix, sizeof(unsigned char));
    }
    
    if (VERBOSE>1) {
      printf("%s numcols = %i, numrows = %i, numelem = %i, numrim = %i\n",
        prefix, numcols, numrows, numelem, numrim);
      printf("%s numSOS1rows,  = %i, numSOS1elem = %i, numSOS2rows,  = %i, numSOS2elem = %i\n",
        prefix, numSOS1rows, numSOS1elem, numSOS2rows, numSOS2elem);
    }

    // copy over Jacobian column starts and row indices
    index_2_starts(omi_model->jacobian[0]->colIndex, start0, 0, numcols+1);
    start0[numcols] = numelem;  // cbc has to have numcols element of starts set to numelem
    if (VERBOSE>3) {
      printf("%s after index_2_starts start0[%i] = %i\n", prefix, numcols, start0[numcols]); 
      fflush(stdout);
    }
    rowIndex_2_index (omi_model->jacobian[0]->rowIndex, index0, 0, numelem);
    if (VERBOSE>3) {
      printf("%s after rowIndex_2_index start0[%i] = %i\n", prefix, numcols, start0[numcols]); 
      fflush(stdout);
    }
    
    // copy over SOS1 Weights' column starts and row indices
    // row order
    index_2_starts(omi_model->weightsSOS1[0]->rowIndex, weightSOS1Start0, 0, numSOS1rows+1);
    // Cbc_C_Interface requires that numSOS1rows element of starts set to numSOS1elem
    weightSOS1Start0[numSOS1rows] = numSOS1elem;  
    if (VERBOSE>3) {
      printf("%s after index_2_starts weightStart0[%i] = %i\n", prefix, 
        numSOS1rows, weightSOS1Start0[numSOS1rows]); 
      fflush(stdout);
    }
    colIndex_2_index (omi_model->weightsSOS1[0]->colIndex, weightSOS1Index0, 0, numSOS1elem);
    if (VERBOSE>3) {
      printf("%s after colIndex_2_index weightStart0[%i] = %i\n", prefix, 
        numSOS1rows, weightSOS1Start0[numSOS1rows]); 
      fflush(stdout);
    }
    
    // copy over SOS2 Weights' column starts and row indices
    // row order
    index_2_starts(omi_model->weightsSOS2[0]->rowIndex, weightSOS2Start0, 0, numSOS2rows+1);
     // Cbc_C_Interface requires that numSOS2rows element of starts set to numSOS2elem
    weightSOS2Start0[numSOS2rows] = numSOS2elem; 
    if (VERBOSE>1) {
      printf("%s after index_2_starts weightSOS2Start0[%i] = %i\n", prefix, 
        numSOS2rows, weightSOS2Start0[numSOS2rows]); 
      fflush(stdout);
    }
    colIndex_2_index(omi_model->weightsSOS2[0]->colIndex, weightSOS2Index0, 
      0, numSOS2elem);
    if (VERBOSE>1) {
      printf("%s after colIndex_2_index weightSOS2Start0[%i] = %i\n", prefix, 
        numSOS2rows, weightSOS2Start0[numSOS2rows]); 
      fflush(stdout);
    }
    
    // set some omi_model parameters
    if (VERBOSE>1) printf("%s set some omi_model parameters\n", prefix);
    Cbc_setOptimizationDirection(cbc_model, omi_model->objective_direction); // (1 - minimize, -1 - maximize, 0 - ignore)
    Cbc_setObjectiveOffset      (cbc_model, omi_model->objective_constant);
    Cbc_setMaximumIterations    (cbc_model, omi_model->iterations_limit);
    
    {
      // jacobian in column order
      const CoinBigIndex * start = start0;
      const int          * index = index0;
      const double       * value = omi_model->jacobian[0]->value;
      // weightsSOS1 in row order
      const CoinBigIndex * weightSOS1Start = weightSOS1Start0;
      const int          * weightSOS1Index = weightSOS1Index0;
      const double       * weightsSOS1     = omi_model->weightsSOS1[0]->value;
      // weightsSOS2 in row order
      const CoinBigIndex * weightSOS2Start = weightSOS2Start0;
      const int          * weightSOS2Index = weightSOS2Index0;
      const double       * weightsSOS2     = omi_model->weightsSOS2[0]->value;
      // rim
      const double* obj   = omi_model->objective;
      const double* rowlb = omi_model->row_lower_bound;
      const double* rowub = omi_model->row_upper_bound;
      const double* collb = omi_model->col_lower_bound;
      const double* colub = omi_model->col_upper_bound;
  
      Cbc_setLogLevel(cbc_model, 0); // 4 is verbose
      if (VERBOSE>1) printf("%s Log Level %i\n", prefix, Cbc_logLevel(cbc_model));
      if (VERBOSE>1) printf("%s calling Cbc_loadProblem\n", prefix);
      Cbc_loadProblem (cbc_model, numcols, numrows, 
        start, index, value, collb, colub, obj, rowlb, rowub);
      if (VERBOSE>1) printf("%s finished Cbc_loadProblem\n", prefix);

      // set integer variables
      long_2_char(omi_model->is_integer, is_integer0, 0, numcols);
      Cbc_copyInIntegerInformation(cbc_model, is_integer0);
      
      // create SOS1 constraints
      if (VERBOSE>1) printf("%s create SOS1 constraints, numSOS1elem = %i\n",prefix,numSOS1elem);
      if (numSOS1elem > 0) {
        int                  row, col;
        int                  numSOS1Cols; // number of SOS columns
        const CoinBigIndex * colIndex;    // column indices
        const double       * colWeight;   // column weights
        assert(numSOS1elem == weightSOS1Start[numSOS1rows]);
        if (VERBOSE>2) {
          printf("%s numSOS1rows = %i\n",prefix,numSOS1rows);
            for (row=0;row<numSOS1rows;row++) {
            printf("%s row = %i\n",prefix,row);
            printf("%s weightSOS1Start[%i] = %i\n",prefix,row,weightSOS1Start[row]);
            printf("%s weightSOS1Start[%i+1] = %i\n",prefix,row,weightSOS1Start[row+1]);
            numSOS1Cols = weightSOS1Start[row+1]-weightSOS1Start[row];
            colIndex    = weightSOS1Index+weightSOS1Start[row];
            colWeight   = weightsSOS1+weightSOS1Start[row];
            printf("%s numSOS1Cols = %i\n",prefix,numSOS1Cols);
            printf("%s colIndex [%i] = %i\n",prefix,0,colIndex[0]);
            printf("%s colIndex [%i] = %i\n",prefix,numSOS1Cols-1,colIndex[numSOS1Cols-1]);
            if (VERBOSE>3 && numSOS1Cols>0) for (col=0;col<numSOS1Cols;col++) {
              printf("%s colIndex [%i] = %i\n",prefix,col,colIndex[col]);
              printf("%s colWeight[%i] = %f\n",prefix,col,colWeight[col]);
            }
            fflush(stdout);
          }
        }
        if (VERBOSE>1) printf("%s calling Cbc_addSOS_Sparse(type = 1)\n",prefix);
        Cbc_addSOS_Sparse(cbc_model, weightSOS1Start, weightSOS1Index, weightsSOS1, 1);
      } else {
        if (VERBOSE>1) printf("%s No SOS1 constraints.\n",prefix);
      }
      // create SOS2 constraints
      if (VERBOSE>1) printf("%s create SOS2 constraints, numSOS2elem = %i\n",prefix,numSOS2elem);
      if (numSOS2elem > 0) {
        int                  row, col;
        int                  numSOS2Cols; // number of SOS columns
        const CoinBigIndex * colIndex;    // column indices
        const double       * colWeight;   // column weights
        assert(numSOS2elem == weightSOS2Start[numSOS2rows]);
        if (VERBOSE>2) {
          printf("%s numSOS1rows = %i\n",prefix,numSOS1rows);
          for (row=0;row<numSOS2rows;row++) {
            printf("%s row = %i\n",prefix,row);
            printf("%s weightSOS2Start[%i] = %i\n",prefix,row,weightSOS2Start[row]);
            printf("%s weightSOS2Start[%i+1] = %i\n",prefix,row,weightSOS2Start[row+1]);
            numSOS2Cols = weightSOS2Start[row+1]-weightSOS2Start[row];
            colIndex    = weightSOS2Index+weightSOS2Start[row];
            colWeight   = weightsSOS2+weightSOS2Start[row];
            printf("%s numSOS2Cols = %i\n",prefix,numSOS2Cols);
            printf("%s colIndex [%i] = %i\n",prefix,0,colIndex[0]);
            printf("%s colIndex [%i] = %i\n",prefix,numSOS2Cols-1,colIndex[numSOS2Cols-1]);
            if (VERBOSE>3 && numSOS2Cols>0) for (col=0;col<numSOS2Cols;col++) {
              printf("%s colIndex [%i] = %i\n",prefix,col,colIndex[col]);
              printf("%s colWeight[%i] = %f\n",prefix,col,colWeight[col]);
            }
            fflush(stdout);
          }
        }
        if (VERBOSE>1) printf("%s calling Cbc_addSOS_Sparse(type = 2)\n",prefix);
        Cbc_addSOS_Sparse(cbc_model, weightSOS2Start, weightSOS2Index, weightsSOS2, 2);
      } else {
        if (VERBOSE>1) printf("%s No SOS2 constraints.\n",prefix);
      }
    }
    
    if (VERBOSE> 2) {
      printf("%s omi_model->iterations_limit = %i\n",prefix,omi_model->iterations_limit);
      printf("%s Cbc maximumIterations = %i\n",prefix, Cbc_maximumIterations(cbc_model));
    }
    // starting solution is in omi_model->basis 
    // convert it from long to unsigned char
    if (VERBOSE>1) printf("%s set starting solution from omi_model->basis\n", prefix);
    if (VERBOSE>1) printf("%s numrim = %i, sizeof(unsigned char) = %i\n", 
      prefix, numrim, sizeof(unsigned char));
    if (VERBOSE>1) printf("%s numcols = %i, numrows = %i, numrim = %i\n",
      prefix, numcols, numrows, numrim);
    // NOTE: That CLP is row then column, whereas OMI convention is column then row
    basis_2_status(omi_model->basis, &status0[numrows], 0, numcols);
    basis_2_status(&(omi_model->basis[numcols]), status0, 0, numrows);
    for (i=0; i<numrim; i++) 
      if (VERBOSE>10) printf("%s after basis_2_status, omi_model->basis[%i] = %i, status0[%i] = %i\n",
        prefix, i, omi_model->basis[i], i, status0[i]);
    if (VERBOSE>1) printf("%s calling Cbc_copyinStatus(cbc_model, status0)\n",prefix);
//tbd    Cbc_copyinStatus(cbc_model, status0);
    if (VERBOSE>1) printf("%s finished Cbc_copyinStatus(cbc_model, status0)\n",prefix);

    // print cbc_model
    if (VERBOSE>5) Cbc_printModel(cbc_model, prefix);
  }

  // write the MPS file
  Cbc_writeMps(cbc_model, file_name);

  // delete the Cbc_model
  fprintf(out,"%s delete the Cbc_model\n",prefix);
  fprintf(stdout,"%s delete the Cbc_model\n",prefix);
  fflush(out); fflush(stdout);
  Cbc_deleteModel(cbc_model);

  // release booty gotten
  if (VERBOSE>1) fprintf(stdout,"%s release booty gotten\n", prefix);  fflush(stdout);
  (*env)->ReleaseStringUTFChars(env, jfile_name, file_name);
  OMI_Model_release(env, jmodel, omi_model, NULL);

  if (VERBOSE>0) fprintf(out,"%sreturn %d\n", prefix, status);
  if (VERBOSE>0) fprintf(stdout,"%sreturn %d\n", prefix, status);
  fflush(out); fflush(stdout);
  
  // close output file
  if (VERBOSE>0) fprintf (out,"%sreturn, closing output file",prefix);
  fclose(out);
  return status;

}  // Java_com_epri_omi_Model_pwriteMPS

/**
  * Converts row indices from cbc to omi format.
  * Cbc counts from zero && omi counts from one.
  */
void index_2_rowIndex(const int *from, jint *to, int start, int end){
  int i;
  if (VERBOSE>5) int_print(stdout, "Model.c:index_2_rowIndex(): ", "from", from, start, end);
  for (i=start; i<end; i++) to[i] = from[i]+1;
  if (VERBOSE>5) long_print(stdout, "Model.c:index_2_rowIndex(): ", "to", to, start, end);
  return;
}
/**
  * Converts row indices from omi to cbc format.
  * Cbc counts from zero && omi counts from one.
  */
void rowIndex_2_index(const jint *from, int *to, int start, int end){
  int i;
  if (VERBOSE>5) long_print(stdout, "Model.c:rowIndex_2_index(): ", "from", from, start, end);
  for (i=start; i<end; i++) to[i] = from[i]-1;
  if (VERBOSE>5) int_print(stdout, "Model.c:rowIndex_2_index(): ", "to", to, start, end);
  return;
}

/**
  * Converts column indices from omi to cbc format.
  * Used for integer information.
  */
void colIndex_2_index(const jint *from, int *to, int start, int end){
  int i;
  if (VERBOSE>5) long_print(stdout, "Model.c:colIndex_2_index(): ", "from", from, start, end);
  for (i=start; i<end; i++) to[i] = from[i]-1;
  if (VERBOSE>5) int_print(stdout, "Model.c:colIndex_2_index(): ", "to", to, start, end);
  return;
}
/**
  * Converts char to long.
  * Used for row major format in weights.
  * Cbc counts from zero && omi counts from one.
  */
void char_2_long(const char *from, jint *to, int start, int end){
//  const int VERBOSE = 6;
  int i;
  if (VERBOSE>5) char_print(stdout, "Model.c:char_2_long(): ", "from", from, start, end);
  for (i=start; i<end; i++) to[i] = from[i];
  if (VERBOSE>5) long_print(stdout, "Model.c:char_2_long(): ", "to", to, start, end);
  return;
}
/**
  * Converts char to long.
  * Used for row major format in weights.
  * Cbc counts from zero && omi counts from one.
  */
void long_2_char(const jint *from, char *to, int start, int end){
//  const int VERBOSE = 6;
  int i;
  if (VERBOSE>5) long_print(stdout, "Model.c:char_2_long(): ", "from", from, start, end);
  for (i=start; i<end; i++) to[i] = from[i];
  if (VERBOSE>5) char_print(stdout, "Model.c:char_2_long(): ", "to", to, start, end);
  return;
}

/**
  * Converts column starts from cbc to omi format.
  * Cbc counts from zero && omi counts from one.
  */
void starts_2_colIndex(const CoinBigIndex *from, jint *to, int start, int end){
  int i;
  if (VERBOSE>5) int_print(stdout, "Model.c:starts_2_colIndex(): ", "from", from, start, end);
  for (i=start; i<end; i++) to[i] = from[i]+1;
  if (VERBOSE>5) long_print(stdout, "Model.c:starts_2_colIndex(): ", "to", to, start, end);
  return;
}
/**
  * Converts column starts from omi to cbc format.
  * Cbc counts from zero && omi counts from one.
  */
void index_2_starts(const jint *from, CoinBigIndex *to, int start, int end){
  int i;
  if (VERBOSE>5) long_print(stdout, "Model.c:index_2_starts(): ", "from", from, start, end);
  for (i=start; i<end; i++) to[i] = from[i]-1;
  if (VERBOSE>5) int_print(stdout, "Model.c:index_2_starts(): ", "to", to, start, end);
  return;
}

/**
  * Converts status from cbc to omi basis format
  * this is mostly to convert from unsigned char to long
  */
void status_2_basis(const unsigned char *from, long *to, int start, int end){
  int i;
  if (VERBOSE>5) unsigned_char_print(stdout, "Model.c:status_2_basis(): ", "from", from, start, end);
//  for (i=start; i<end; i++) to[i] = from[i];
  if (VERBOSE>5) long_print(stdout, "Model.c:status_2_basis(): ", "to", to, start, end);
  return;
}

/**
  * Converts from omi basis to cbc status format
  * this is mostly to convert from long to unsigned char
  */
void basis_2_status(const long *from, unsigned char  *to, int start, int end){
  int i;
  if (VERBOSE>5) long_print(stdout, "Model.c:basis_2_status(): ", "from", from, start, end);
//  for (i=start; i<end; i++) to[i] = from[i];
  if (VERBOSE>5) unsigned_char_print(stdout, "Model.c:basis_2_status(): ", "to", to, start, end);
  return;
}

/**
  * Copies double values in an array
  */
void dcopy(const double *from, double *to, int start, int end) {
  int i;
  if (VERBOSE>5) double_print(stdout, "Model.c:dcopy(): ", "from", from, start, end);
  for (i=start; i<end; i++) to[i] = from[i];
  return;
  if (VERBOSE>5) double_print(stdout, "Model.c:dcopy(): ", "to", to, start, end);
}

/**
  * Copies int values in an array
  */
void iadd(int val, int *array, int start, int end) {
  int i;
  for (i=start; i<end; i++) array[i] += val;
  return;
}

/**
  * Copies jchar values in an array
  */
void jccopy(const jchar *from, jchar *to, int start, int end) {
  int i;
  for (i=start; i<end; i++) to[i] = from[i];
  return;
}

/**
  * Prints an array of unsigned char to the FILE argument
  */
void unsigned_char_print(FILE *out, const char *prefix, char* name, const unsigned char *array, int start, int end) {
  int i;
  for (i=start; i<end; i++) 
    fprintf(out, "%s %s[%i] = %i\n",prefix, name, i, array[i]);
  return;
}

/**
  * Prints an array of char to the FILE argument
  */
void char_print(FILE *out, const char *prefix, char* name, const char *array, int start, int end) {
  int i;
  for (i=start; i<end; i++) 
    fprintf(out, "%s %s[%i] = %i\n",prefix, name, i, array[i]);
  return;
}

/**
  * Prints an array of int to the FILE argument
  */
void int_print(FILE *out, const char *prefix, char* name, const int *array, int start, int end) {
  int i;
  for (i=start; i<end; i++) 
    fprintf(out, "%s %s[%i] = %i\n",prefix, name, i, array[i]);
  return;
}

/**
  * Prints an array of long to the FILE argument
  */
void long_print(FILE *out, const char *prefix, char* name, const long *array, int start, int end) {
  int i;
  for (i=start; i<end; i++) 
    fprintf(out, "%s %s[%i] = %i\n",prefix, name, i, array[i]);
  return;
}

/**
  * Prints an array of CoinBigIndex to the FILE argument
  */
void CoinBigIndex_print(FILE *out, const char *prefix, char* name, CoinBigIndex *array, int start, int end) {
  int i;
  for (i=start; i<end; i++) 
    fprintf(out, "%s %s[%i] = %i\n",prefix, name, i, array[i]);
  return;
}

/**
  * Prints an array of double to the FILE argument
  */
void double_print(FILE *out, const char *prefix, char* name, const double *array, int start, int end) {
  int i;
  for (i=start; i<end; i++) 
    fprintf(out, "%s %s[%i] = %g\n",prefix, name, i, array[i]);
  return;
}

/**
  * Prints the solution found by Cbc using its own methods
  * Good for debugging.
  */
void printSolution(Cbc_Model *cbc_model) {
  //*
  //  Now to print out solution.  The methods used return modifiable
  //  arrays while the alternative names return const pointers -
  //  which is of course much more virtuous.
  //  
  //  This version just does non-zero column and row values.
  //
  //*
  
  //* If we have not kept names (parameter to readMps) this will be 0 
//   assert(Cbc_lengthNames(cbc_model));
  {
    int  name_length = OMI_MaxStrLen;
    char model_name[OMI_MaxStrLen];
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
  printf("Objective Sense = %g\n",Cbc_getObjSense(cbc_model));  // (1 - minimize, -1 - maximize, 0 - ignore)
  printf("Primal Feasible = %i\n",Cbc_primalFeasible(cbc_model));
  printf("Dual Feasible = %i\n",Cbc_dualFeasible(cbc_model));
  printf("Dual Bound = %g\n",Cbc_dualBound(cbc_model));
  printf("Infeasibility Cost = %g\n",Cbc_infeasibilityCost(cbc_model));
  printf("Sum Dual Infeasibilities = %g\n",Cbc_sumDualInfeasibilities(cbc_model));
  printf("Number Dual Infeasibilities = %i\n",Cbc_numberDualInfeasibilities(cbc_model));
  printf("Sum Primal Infeasibilities = %g\n",Cbc_sumPrimalInfeasibilities(cbc_model));
  printf("Number Primal Infeasibilities = %i\n",Cbc_numberPrimalInfeasibilities(cbc_model));
  printf("Objective Value = %g\n",Cbc_objectiveValue(cbc_model)); 
  printf("--------------------------------------\n");

  //* Rows 
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
  //* Columns 
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
  //    	Cbc_columnName(cbc_model,iColumn,name);
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
} //  printSolution()
