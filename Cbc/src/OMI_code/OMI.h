/**
  * @(#)OMI.h
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
#ifndef _FN_OMI_java
#define _FN_OMI_java

#include <JavaVM/jni.h>

#ifndef F2C_INCLUDE
#define F2C_INCLUDE
  typedef long int integer;
  typedef double   doublereal;
  typedef long int ftnlen;
#endif

/* OMI Global Constants */
#define OMI_MaxStrLen 256
/*                    NULL */     /* copy & free   */
#define OMI_COMMIT    JNI_COMMIT  /* copy w/o free */
#define OMI_ABORT     JNI_ABORT   /* free w/o copy */

/* OMI Callbacks */
void omi_model_getobjective__(integer *mode, integer *n, doublereal *x,
  doublereal *f, doublereal *g, integer* nstate, integer* nprob);

#endif