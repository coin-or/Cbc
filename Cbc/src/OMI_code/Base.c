/**
  * @(#)Base.c
  * @author Robert Entriken<br>
  * Copyright (C) 1999-2007 EPRI<br>
  * All rights reserved.
  * @version 07-03-07
  * @since OMI_1.0
  *
  * Revisions:
  *   02-06-16 bobe     Switched to MW generated headers
  *   02-08-10 bobe     Got rid of OMI_Base_ppprint
  *   04-05-07 bobe     made to work with clp
  *   06-02-20 bobe     revised copyright
  *   06-01-28 jmm      new package layout + casts
  *   06-04-10 jmm      cast NULL to jint to avoid warnings
  *   07-03-07 bobe     updated Copyright
  *
  */

#include <stdio.h>
#include <time.h>
#include "omi_Base.h"
#include "com_epri_omi_Base.h" 

/* Begin Block Signatures
Compiled from Base.java
public class omi.Base extends java.lang.Object 
{
    row I
    col I
}
End Block Signatures */

/* native methods */
#pragma export on
/**
  * Copy scalars from C back to Java and free C version of base.
  * mode = 0,          copy and free
  * mode = OMI_COMMIT, copy w/o free
  * mode = OMI_ABORT,  free w/o copy
  *
  */
void OMI_Base_release
  (JNIEnv *env, jobject jbase, ClassOMI_Base  *base, int mode)
{
  extern void free(ClassOMI_Base *);

  int            VERBOSE = 0;
  jclass         cls;
  jfieldID       fid;
  
  /* Check arguments */
  if (env == NULL) {
    printf("Base.c:OMI_Base_release ",
      "!ERROR, NULL 'env' passed as argument.\n");
    return;
  }  
  if (jbase == NULL) {
    printf("Base.c:OMI_Base_release ",
      "!ERROR, NULL 'jbase' passed as argument.\n");
    return;
  }
  if (base == NULL) {
    printf("Base.c:OMI_Base_release ",
      "!ERROR, NULL 'base' passed as argument.\n");
    return;
  }
  if ((mode != (jint)NULL) && (mode != OMI_COMMIT) && (mode != OMI_ABORT)) {
    printf("Base.c:OMI_Base_release ",
      "!ERROR, 'mode' is not one of NULL, OMI_COMMIT, or OMI_ABORT.\n");
    return;
  }
  /* Announce OMI_Base_release */
  if (VERBOSE>0) printf ("OMI_Base_release() Arguments OK.\n");

  if ((mode == (jint)NULL) || (mode == OMI_COMMIT)) {
    /* cls
     * In JDK1.1 use the class to get at fields.
     */
    if ((cls = (*env)->GetObjectClass(env, jbase)) == NULL) {
      printf("Base.c:OMI_Base_release ",
        "!ERROR, Cannot GetObjectClass of 'jbase'.\n");
      return;
    }
    
    /* row */
    fid = (*env)->GetFieldID(env, cls, "row", "I");
    if (fid == 0) {
      printf ("Base.c:OMI_Base_release ",
        "!ERROR: Cannot GetFieldID of base.row\n");
      return;
    }
    (*env)->SetIntField(env, jbase, fid, base->row); 
    /* col */
    fid = (*env)->GetFieldID(env, cls, "col", "I");
    if (fid == 0) {
      printf ("Base.c:OMI_Base_release ",
        "!ERROR: Cannot GetFieldID of base.col\n");
      return;
    }
    (*env)->SetIntField(env, jbase, fid, base->col); 
  }  /* if NULL or OMI_COMMIT */
  
  if ((mode == (jint)NULL) ||  (mode == OMI_ABORT)) {
    free(base);
  }  /* if NULL or OMI_ABORT */
}  /* OMI_Base_release */

/**
  * Copy scalars from Java into C.
  *
  */
ClassOMI_Base *OMI_Base_get
  (JNIEnv *env, jobject jbase)
{
  extern malloc(int);
  
  int            VERBOSE = 0;
  jclass         cls;
  jfieldID       fid;
  ClassOMI_Base  *base;
  
  /* Check arguments */
  if (env == NULL) {
    printf("Base.c:OMI_Base_get ",
      "!ERROR, NULL 'env' passed as argument.\n");
    return NULL;
  }  
  if (jbase == NULL) {
    printf("Base.c:OMI_Base_get ",
      "!ERROR, NULL 'jbase' passed as argument.\n");
    return NULL;
  }

  /* Announce OMI_Base_get */
  if (VERBOSE>0) printf ("Base->OMI_Base_get() Attributes OK.\n");

  /* malloc space for a base */
  if ((base = (ClassOMI_Base *) malloc(sizeof(ClassOMI_Base))) == NULL) {
    printf("Model.c:OMI_Base_get ",
      "!ERROR, Cannot malloc space for base\n");
    return NULL;
  }
  
  /* cls
   * In JDK1.1 use the class to get at fields.
   */
  if ((cls = (*env)->GetObjectClass(env, jbase)) == NULL) {
    printf("Base.c:OMI_Base_get ",
      "!ERROR, Cannot GetObjectClass of 'jbase'.\n");
    OMI_Base_release(env, jbase, base, OMI_ABORT);
    return NULL;
  }
  
  /* In JDK1.1 work from the name of the field of the class to get the field id,
   * then with the field id get the value.
   * row
   */
  fid = (*env)->GetFieldID(env, cls, "row", "I");
  if (fid == 0) {
    printf ("Base.c:OMI_Base_get ",
      "!ERROR: Cannot GetFieldID of base.row\n");
    return NULL;
  }
  base->row = (*env)->GetIntField(env, jbase, fid);  /* Exceptions? TBD */
  /* col */
  fid = (*env)->GetFieldID(env, cls, "col", "I");
  if (fid == 0) {
    printf ("Base.c:OMI_Base_get ",
      "!ERROR: Cannot GetFieldID of base.col\n\n");
    return NULL;
  }
  base->col = (*env)->GetIntField(env, jbase, fid);  /* Exceptions? TBD */
  return base;
}  /* OMI_Base_get */

/**
  * Print the attributes of this Base to the given FILE argument.
  *
  */
void OMI_Base_oprint 
  (JNIEnv *env, jobject jbase, jstring jprefix, FILE *out)
{
  /* Declare C versions of base and prefix */
  ClassOMI_Base  *base;
  const char     *prefix;

  /* Check arguments */
  if (env == NULL) {
    fprintf(out,"Base.c:OMI_Base_oprint ",
      "!ERROR, NULL 'env' passed as argument.\n");
    return;
  }  
  if (jbase == NULL) {
    fprintf(out,"Base.c:OMI_Base_oprint ",
      "!ERROR, NULL 'jbase' passed as argument.\n");
    return;
  }
  if (jprefix == NULL) {
    fprintf(out,"Base.c:OMI_Base_oprint ",
      "!ERROR, NULL 'jprefix' passed as argument.\n");
    return;
  }
  /* Announce pprint */
  if ((prefix = (*env)->GetStringUTFChars(env, jprefix, 0)) == NULL) {
    fprintf(out,"Base.c:OMI_Base_oprint ",
      "!ERROR, Cannot GetStringUTFChars of 'jprefix'.\n");
    return;
  }
  fprintf(out,"%sBase->pprint() Attributes OK.\n", prefix);
  /* get base and print it */
  if ((base = OMI_Base_get(env, jbase)) == NULL) {
    fprintf(out,"Base.c:OMI_Base_oprint ",
      "!ERROR, Cannot OMI_Base_get of 'jbase'.\n");
    return;
  }
  fprintf(out,"%sBase->row @ %X = %i\n", prefix, &base->row, base->row);
  fprintf(out,"%sBase->col @ %X = %i\n", prefix, &base->col, base->col);
  /* release booty gotten */
  OMI_Base_release(env, jbase, base, (jint)NULL);
  (*env)->ReleaseStringUTFChars(env, jprefix, prefix);
}  /* OMI_Base_oprint */

/*
 * Class:     com_epri_omi_Base
 * Method:    pprint
 * Signature: (Ljava/lang/String;)V
 *
 *  JNI call from Java to pprint()
 */
JNIEXPORT void JNICALL Java_com_epri_omi_Base_pprint
  (JNIEnv *env, jobject jbase, jstring jprefix)
{
  FILE *out;
  time_t  t = time(NULL);
  if ((out = fopen("Java_Native_Output", "a")) == NULL) {
    printf("Base.c:Java_com_epri_omi_Base_pprint ",
      "!ERROR, Cannot fopen Java_Native_Output.\n");
    return;
  }
  fprintf(out,"\nBase->pprint('') %s\n", asctime(localtime(&t)));
  OMI_Base_oprint(env, jbase, jprefix, out);
  fclose(out);
}  /* Java_com_epri_omi_Base_pprint */
#pragma export off
