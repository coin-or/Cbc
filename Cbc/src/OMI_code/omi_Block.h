/**
  * @(#)omi_Block.h
  * @author Robert Entriken<br>
  * Copyright (C) 1999-2007 EPRI<br>
  * All rights reserved.
  * @version 07-03-07
  * @since OMI_1.0
  *
  * Revisions:
  *   06-02-20 bobe     added copyright
  *   06-04-10 bobe     ANSI standard comments
  *   07-03-07 bobe     updated Copyright
  *
  */
#ifndef _FN_Block_java
#define _FN_Block_java

#include <JavaVM/jni.h>
#include "OMI.h"
#include "omi_Base.h"

#ifndef F2C_INCLUDE
#define F2C_INCLUDE
  typedef long int integer;
  typedef double   doublereal;
  typedef long int ftnlen;
#endif

typedef struct ClassOMI_Block {
    jint                 rowIndexLength;
    jint                 columnIndexLength;
    jint                 valueLength;
    jint                 numRows;
    jint                 numColumns;
    jint                 numElements;
    jint                 blockType;
    jobject              joffset;
    struct ClassOMI_Base *offset;
  	jintArray            jrowIndex;
    jint                 *rowIndex;
	  jintArray            jcolIndex;
    jint                 *colIndex;
	  jdoubleArray         jvalue;
    jdouble              *value;
} ClassOMI_Block;

/* function prototypes for JNI methods */
void OMI_Block_release (JNIEnv *env, jobject jblock, ClassOMI_Block *block, int mode);
ClassOMI_Block *OMI_Block_get (JNIEnv *env, jobject jblock);
void OMI_Block_oprint (JNIEnv *env, jobject jblock, jstring jprefix, FILE *out);

/* function prototypes for callback methods */
integer Java_com_epri_omi_Block_get_rowIndex_length (JNIEnv *env, jobject jblock);
integer Java_com_epri_omi_Block_get_colIndex_length (JNIEnv *env, jobject jblock);
integer Java_com_epri_omi_Block_get_value_length (JNIEnv *env, jobject jblock);
integer Java_com_epri_omi_Block_get_num_elems (JNIEnv *env, jobject jblock);
integer *Java_com_epri_omi_Block_get_rowIndex (JNIEnv *env, jobject jblock);
integer *Java_com_epri_omi_Block_get_colIndex (JNIEnv *env, jobject jblock);
doublereal *Java_com_epri_omi_Block_get_value (JNIEnv *env, jobject jblock);
integer Java_com_epri_omi_Block_set_num_elems (JNIEnv *env, jobject jblock, jint len);
integer Java_com_epri_omi_Block_release_rowIndex (JNIEnv *env, jobject jblock, jint *rowIndex);
integer Java_com_epri_omi_Block_release_colIndex (JNIEnv *env, jobject jblock, jint *colIndex);
integer Java_com_epri_omi_Block_release_value (JNIEnv *env, jobject jblock, jdouble *value);

void Java_com_epri_omi_block_convertToRowOrder (JNIEnv *env, jobject jblock);
void Java_com_epri_omi_block_convertToColumnOrder (JNIEnv *env, jobject jblock);
void Java_com_epri_omi_block_convertToTriples(JNIEnv *env, jobject jblock);

#endif
