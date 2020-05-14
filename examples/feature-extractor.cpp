// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/*! \file feature-extractor.cpp
  \brief Example of using the class OsiFeatures to extract problem features from a Mixed-Integer Linear Program

*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <OsiClpSolverInterface.hpp>
#include <OsiFeatures.hpp>

char *basename(char *dest, const char *fileWithPath);

int main( int argc, char ** argv ) {
  if (argc<2) {
    fprintf(stderr, "Enter instance name or -header to print feature names.\n");
    exit(1);
  }

  if (strcmp(argv[1], "-header") == 0) {
    printf("%s", OsiFeatures::name(0));
    for ( int i=1 ; (i<OsiFeatures::n) ; ++i )
      printf(",%s", OsiFeatures::name(i));
    printf("\n"); fflush(stdout);
    exit(0);
  }

  OsiClpSolverInterface solver;
  solver.messageHandler()->setLogLevel(0);

  if (strstr(argv[1], ".lp") || strstr(argv[1], ".LP"))
    solver.readLp(argv[1]);
  else
    solver.readMps(argv[1]);
  char instance[256];
  basename(instance, argv[1]);

  double *features = new double[OsiFeatures::n];
  OsiFeatures::compute(features, &solver);

  printf("%s", instance);
  for ( int i=0 ; i<OsiFeatures::n ; ++i )
    printf(",%g", features[i]);
  printf("\n"); fflush(stdout);

  delete[] features;
}

char *basename(char *dest, const char *fileWithPath) {
  const char *s = fileWithPath;
  for ( const char *sl=fileWithPath ; *sl!='\0' ; ++sl ) {
    if (*sl == '/' || *sl == '\\')
      s = sl+1;
  }

  strcpy(dest, s);
  char *s2 = strstr(dest, ".mps");
  if (s2)
    *s2 = '\0';
  s2 = strstr(dest, ".MPS");
  if (s2)
    *s2 = '\0';
  s2 = strstr(dest, ".mps.gz");
  if (s2)
    *s2 = '\0';
  s2 = strstr(dest, ".MPS.GZ");
  if (s2)
    *s2 = '\0';
  s2 = strstr(dest, ".lp");
  if (s2)
    *s2 = '\0';
  s2 = strstr(dest, ".LP");
  if (s2)
    *s2 = '\0';
  s2 = strstr(dest, ".lp.gz");
  if (s2)
    *s2 = '\0';
  s2 = strstr(dest, ".LP.GZ");
  if (s2)
    *s2 = '\0';

  return dest;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
