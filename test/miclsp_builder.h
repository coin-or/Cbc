/* Shared MICLSP model builder, extracted from CInterfaceTest_miclsp.c.
 * Include this header in any file that needs to build MICLSP instances
 * via the C API.  All functions are static so the header is safe to include
 * from multiple translation units.
 */

#ifndef MICLSP_BUILDER_H
#define MICLSP_BUILDER_H

#include "Cbc_C_Interface.h"
#include "miclsp_fixtures.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Maximum dimensions used for stack-allocated arrays */
#define MAX_I 8
#define MAX_T 24

static Cbc_Model *miclsp_new_model(void)
{
  Cbc_Model *m = Cbc_newModel();
  Cbc_setObjSense(m, 1.0); /* minimise */
  Cbc_setLogLevel(m, 0);
  return m;
}

static void miclsp_add_row(Cbc_Model *m, int nz, int *idx, double *coef,
  char sense, double rhs, const char *name)
{
  Cbc_addRow(m, name, nz, idx, coef, sense, rhs);
}

/* Column index helpers — column layout per instance:
 *   [0 .. I*T)        : y[i,t]  binary
 *   [I*T .. 2*I*T)    : n[i,t]  integer
 *   [2*I*T .. 3*I*T)  : s[i,t]  continuous
 */
static int yidx(int T, int i, int t) { return i * T + t; }
static int nidx(int I, int T, int i, int t) { return I * T + i * T + t; }
static int sidx(int I, int T, int i, int t) { return 2 * I * T + i * T + t; }

/* Build and populate a MICLSP model from fixture data.
 * Returns a newly allocated model that the caller must free with Cbc_deleteModel. */
static Cbc_Model *build_miclsp(const MiclspInstance *inst)
{
  int I = inst->I, T = inst->T;
  Cbc_Model *m = miclsp_new_model();

  /* y[i,t]: binary, objective coefficient f[i] */
  for (int i = 0; i < I; i++) {
    for (int t = 0; t < T; t++) {
      char name[32];
      snprintf(name, sizeof(name), "y(%d,%d)", i, t);
      Cbc_addCol(m, name, 0.0, 1.0, (double)inst->f[i], 1, 0, NULL, NULL);
    }
  }

  /* n[i,t]: integer >= 0, no objective cost */
  for (int i = 0; i < I; i++) {
    for (int t = 0; t < T; t++) {
      char name[32];
      snprintf(name, sizeof(name), "n(%d,%d)", i, t);
      Cbc_addCol(m, name, 0.0, 1e30, 0.0, 2, 0, NULL, NULL);
    }
  }

  /* s[i,t]: continuous >= 0, objective coefficient h[i] */
  for (int i = 0; i < I; i++) {
    for (int t = 0; t < T; t++) {
      char name[32];
      snprintf(name, sizeof(name), "s(%d,%d)", i, t);
      Cbc_addCol(m, name, 0.0, 1e30, (double)inst->h[i], 0, 0, NULL, NULL);
    }
  }

  /* Inventory balance: B[i]*n[i,t] + s[i,t-1] - s[i,t] = d[i,t] */
  for (int i = 0; i < I; i++) {
    for (int t = 0; t < T; t++) {
      int idx[3];
      double coef[3];
      int nz = 0;
      idx[nz] = nidx(I, T, i, t);
      coef[nz++] = (double)inst->B[i];
      idx[nz] = sidx(I, T, i, t);
      coef[nz++] = -1.0;
      if (t > 0) {
        idx[nz] = sidx(I, T, i, t - 1);
        coef[nz++] = 1.0;
      }
      char name[32];
      snprintf(name, sizeof(name), "inv(%d,%d)", i, t);
      miclsp_add_row(m, nz, idx, coef, 'E', (double)inst->d[i * T + t], name);
    }
  }

  /* Machine capacity: sum_i B[i]*n[i,t] <= C[t] */
  for (int t = 0; t < T; t++) {
    int idx[MAX_I];
    double coef[MAX_I];
    for (int i = 0; i < I; i++) {
      idx[i] = nidx(I, T, i, t);
      coef[i] = (double)inst->B[i];
    }
    char name[32];
    snprintf(name, sizeof(name), "cap(%d)", t);
    miclsp_add_row(m, I, idx, coef, 'L', (double)inst->C[t], name);
  }

  /* Big-M setup linking: n[i,t] <= floor(C[t]/B[i]) * y[i,t] */
  for (int i = 0; i < I; i++) {
    for (int t = 0; t < T; t++) {
      int Mv = inst->C[t] / inst->B[i];
      if (Mv <= 0)
        continue;
      int idx[2] = { nidx(I, T, i, t), yidx(T, i, t) };
      double coef[2] = { 1.0, -(double)Mv };
      char name[32];
      snprintf(name, sizeof(name), "bigm(%d,%d)", i, t);
      miclsp_add_row(m, 2, idx, coef, 'L', 0.0, name);
    }
  }

  return m;
}

#endif /* MICLSP_BUILDER_H */
