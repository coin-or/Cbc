/**
  * @(#)omi_State.h
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
#include <JavaVM/jni.h>
#include "omi_Model.h"

#ifndef _OMI_State
#define _OMI_State

/* OMI Global State */
typedef struct OMI_state {
  JNIEnv         *env;    /* pointer to the JVM    */
  jobject        jmodel;  /* Java version of model */
  ClassOMI_Model *model;  /* C version of model    */
  FILE           *out;
} OMI_state;

OMI_state omi_state;

#endif