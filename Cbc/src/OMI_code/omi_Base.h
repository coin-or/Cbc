/**
  * @(#)omi_Base.h
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

#ifndef _FN_Base_java
#define _FN_Base_java

#ifndef F2C_INCLUDE
#define F2C_INCLUDE
  typedef long int integer;
  typedef double   doublereal;
  typedef long int ftnlen;
#endif

#include "OMI.h"
typedef struct ClassOMI_Base {
    int row;
    int col;
} ClassOMI_Base;

/* function prototypes */
void OMI_Base_release(JNIEnv *env, jobject jbase, ClassOMI_Base  *base, int mode);
ClassOMI_Base * OMI_Base_get(JNIEnv *env, jobject jbase);
void OMI_Base_oprint(JNIEnv *env, jobject jbase, jstring jprefix, FILE *out);

#endif
