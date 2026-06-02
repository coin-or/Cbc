/* mescn_builder.h — Model construction helpers for MESCN tests.
 *
 * Multi-Echelon Supply Chain Network (MESCN) problem:
 *   - Factories produce in batches, ship to warehouses
 *   - Warehouses store and distribute via trucks to retailers
 *   - Retailers satisfy demand
 *
 * Index functions map (facility, product, time) tuples to column indices.
 * All indices are 0-based.
 */

#ifndef MESCN_BUILDER_H
#define MESCN_BUILDER_H

#include "Cbc_C_Interface.h"
#include "mescn_fixtures.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ── Index calculation macros ─────────────────────────────────────────────
 * Variable layout in column order:
 *   x[f,p,t]     : factory production (continuous)
 *   y[f,w,p,t]   : factory-to-warehouse shipment (continuous)
 *   z[w,r,p,t]   : warehouse-to-retail shipment (continuous)
 *   u[f,p,t]     : factory setup (binary)
 *   v[w,p,t]     : warehouse setup (binary)
 *   n[f,p,t]     : production batches (integer)
 *   m[w,p,t]     : shipment trucks (integer)
 *   s_f[f,p,t]   : factory inventory (continuous)
 *   s_w[w,p,t]   : warehouse inventory (continuous)
 */

#define xidx(F, P, T, f, p, t) \
  ((f) * (P) * (T) + (p) * (T) + (t))

#define yidx(F, W, P, T, f, w, p, t) \
  ((F) * (P) * (T) + \
   (f) * (W) * (P) * (T) + (w) * (P) * (T) + (p) * (T) + (t))

#define zidx(F, W, R, P, T, w, r, p, t) \
  ((F) * (P) * (T) + (F) * (W) * (P) * (T) + \
   (w) * (R) * (P) * (T) + (r) * (P) * (T) + (p) * (T) + (t))

#define uidx(F, W, R, P, T, f, p, t) \
  ((F) * (P) * (T) + (F) * (W) * (P) * (T) + (W) * (R) * (P) * (T) + \
   (f) * (P) * (T) + (p) * (T) + (t))

#define vidx(F, W, R, P, T, w, p, t) \
  ((F) * (P) * (T) + (F) * (W) * (P) * (T) + (W) * (R) * (P) * (T) + \
   (F) * (P) * (T) + \
   (w) * (P) * (T) + (p) * (T) + (t))

#define nidx(F, W, R, P, T, f, p, t) \
  ((F) * (P) * (T) + (F) * (W) * (P) * (T) + (W) * (R) * (P) * (T) + \
   (F) * (P) * (T) + (W) * (P) * (T) + \
   (f) * (P) * (T) + (p) * (T) + (t))

#define midx(F, W, R, P, T, w, p, t) \
  ((F) * (P) * (T) + (F) * (W) * (P) * (T) + (W) * (R) * (P) * (T) + \
   (F) * (P) * (T) + (W) * (P) * (T) + (F) * (P) * (T) + \
   (w) * (P) * (T) + (p) * (T) + (t))

#define sf_idx(F, W, R, P, T, f, p, t) \
  ((F) * (P) * (T) + (F) * (W) * (P) * (T) + (W) * (R) * (P) * (T) + \
   (F) * (P) * (T) + (W) * (P) * (T) + (F) * (P) * (T) + (W) * (P) * (T) + \
   (f) * (P) * (T) + (p) * (T) + (t))

#define sw_idx(F, W, R, P, T, w, p, t) \
  ((F) * (P) * (T) + (F) * (W) * (P) * (T) + (W) * (R) * (P) * (T) + \
   (F) * (P) * (T) + (W) * (P) * (T) + (F) * (P) * (T) + (W) * (P) * (T) + \
   (F) * (P) * (T) + \
   (w) * (P) * (T) + (p) * (T) + (t))

/* Total number of variables */
static inline int mescn_num_vars(int F, int W, int R, int P, int T)
{
  return F * P * T +           /* x */
         F * W * P * T +       /* y */
         W * R * P * T +       /* z */
         F * P * T +           /* u */
         W * P * T +           /* v */
         F * P * T +           /* n */
         W * P * T +           /* m */
         F * P * T +           /* s_f */
         W * P * T;            /* s_w */
}

/* Check if factory-warehouse arc exists */
static inline int has_fw_arc(const MESCNInstance *inst, int f, int w)
{
  if (inst->fully_connected)
    return 1;
  for (int i = 0; i < inst->n_fw_arcs; i++) {
    if (inst->fw_arcs[2 * i] == f && inst->fw_arcs[2 * i + 1] == w)
      return 1;
  }
  return 0;
}

/* Check if warehouse-retail arc exists */
static inline int has_wr_arc(const MESCNInstance *inst, int w, int r)
{
  if (inst->fully_connected)
    return 1;
  for (int i = 0; i < inst->n_wr_arcs; i++) {
    if (inst->wr_arcs[2 * i] == w && inst->wr_arcs[2 * i + 1] == r)
      return 1;
  }
  return 0;
}

/* Build MESCN model from instance data.
 * Returns a Cbc_Model ready to solve (structure only, no time/node/cut params).
 */
static Cbc_Model *build_mescn(const MESCNInstance *inst)
{
  int F = inst->F, W = inst->W, R = inst->R, P = inst->P, T = inst->T;
  int n_vars = mescn_num_vars(F, W, R, P, T);

  Cbc_Model *model = Cbc_newModel();

  /* ── Add columns (variables) ──────────────────────────────────────────── */
  double *obj = calloc(n_vars, sizeof(double));
  double *lb = calloc(n_vars, sizeof(double));
  double *ub = calloc(n_vars, sizeof(double));
  char *is_int = calloc(n_vars, sizeof(char));

  /* x[f,p,t] — continuous production */
  for (int f = 0; f < F; f++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        int col = xidx(F, P, T, f, p, t);
        obj[col] = inst->prod_cost[f * P + p];
        lb[col] = 0.0;
        ub[col] = 1e20;
        is_int[col] = 0;
      }
    }
  }

  /* y[f,w,p,t] — continuous factory-to-warehouse shipment */
  for (int f = 0; f < F; f++) {
    for (int w = 0; w < W; w++) {
      if (!has_fw_arc(inst, f, w))
        continue;
      for (int p = 0; p < P; p++) {
        for (int t = 0; t < T; t++) {
          int col = yidx(F, W, P, T, f, w, p, t);
          obj[col] = inst->transport_fw[f * W * P + w * P + p];
          lb[col] = 0.0;
          ub[col] = 1e20;
          is_int[col] = 0;
        }
      }
    }
  }

  /* z[w,r,p,t] — continuous warehouse-to-retail shipment */
  for (int w = 0; w < W; w++) {
    for (int r = 0; r < R; r++) {
      if (!has_wr_arc(inst, w, r))
        continue;
      for (int p = 0; p < P; p++) {
        for (int t = 0; t < T; t++) {
          int col = zidx(F, W, R, P, T, w, r, p, t);
          obj[col] = inst->transport_wr[w * R * P + r * P + p];
          lb[col] = 0.0;
          ub[col] = 1e20;
          is_int[col] = 0;
        }
      }
    }
  }

  /* u[f,p,t] — binary factory setup */
  for (int f = 0; f < F; f++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        int col = uidx(F, W, R, P, T, f, p, t);
        obj[col] = inst->setup_cost_f[f * P + p];
        lb[col] = 0.0;
        ub[col] = 1.0;
        is_int[col] = 1;
      }
    }
  }

  /* v[w,p,t] — binary warehouse setup */
  for (int w = 0; w < W; w++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        int col = vidx(F, W, R, P, T, w, p, t);
        obj[col] = inst->setup_cost_w[w * P + p];
        lb[col] = 0.0;
        ub[col] = 1.0;
        is_int[col] = 1;
      }
    }
  }

  /* n[f,p,t] — integer production batches */
  for (int f = 0; f < F; f++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        int col = nidx(F, W, R, P, T, f, p, t);
        obj[col] = 0.0;
        lb[col] = 0.0;
        ub[col] = 1e20;
        is_int[col] = 1;
      }
    }
  }

  /* m[w,p,t] — integer shipment trucks */
  for (int w = 0; w < W; w++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        int col = midx(F, W, R, P, T, w, p, t);
        obj[col] = 0.0;
        lb[col] = 0.0;
        ub[col] = 1e20;
        is_int[col] = 1;
      }
    }
  }

  /* s_f[f,p,t] — continuous factory inventory */
  for (int f = 0; f < F; f++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        int col = sf_idx(F, W, R, P, T, f, p, t);
        obj[col] = inst->holding_f[f * P + p];
        lb[col] = 0.0;
        ub[col] = 1e20;
        is_int[col] = 0;
      }
    }
  }

  /* s_w[w,p,t] — continuous warehouse inventory */
  for (int w = 0; w < W; w++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        int col = sw_idx(F, W, R, P, T, w, p, t);
        obj[col] = inst->holding_w[w * P + p];
        lb[col] = 0.0;
        ub[col] = 1e20;
        is_int[col] = 0;
      }
    }
  }

  Cbc_loadProblem(model, n_vars, 0, NULL, NULL, NULL,
                  lb, ub, obj, NULL, NULL);

  for (int col = 0; col < n_vars; col++) {
    if (is_int[col])
      Cbc_setInteger(model, col);
  }

  free(obj);
  free(lb);
  free(ub);
  free(is_int);

  /* ── Add constraints ─────────────────────────────────────────────────── */

  /* We'll build constraints one at a time using Cbc_addRow */

  int *indices = malloc(n_vars * sizeof(int));
  double *values = malloc(n_vars * sizeof(double));

  /* Factory flow balance: s_f[f,p,t-1] + x[f,p,t] - s_f[f,p,t] = sum_w y[f,w,p,t] */
  for (int f = 0; f < F; f++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        int nnz = 0;
        /* s_f[f,p,t-1] */
        if (t > 0) {
          indices[nnz] = sf_idx(F, W, R, P, T, f, p, t - 1);
          values[nnz] = 1.0;
          nnz++;
        }
        /* x[f,p,t] */
        indices[nnz] = xidx(F, P, T, f, p, t);
        values[nnz] = 1.0;
        nnz++;
        /* -s_f[f,p,t] */
        indices[nnz] = sf_idx(F, W, R, P, T, f, p, t);
        values[nnz] = -1.0;
        nnz++;
        /* -sum_w y[f,w,p,t] */
        for (int w = 0; w < W; w++) {
          if (has_fw_arc(inst, f, w)) {
            indices[nnz] = yidx(F, W, P, T, f, w, p, t);
            values[nnz] = -1.0;
            nnz++;
          }
        }
        Cbc_addRow(model, "row", nnz, indices, values, 'E', 0.0);
      }
    }
  }

  /* Warehouse flow balance: s_w[w,p,t-1] + sum_f y[f,w,p,t] - s_w[w,p,t] = sum_r z[w,r,p,t] */
  for (int w = 0; w < W; w++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        int nnz = 0;
        /* s_w[w,p,t-1] */
        if (t > 0) {
          indices[nnz] = sw_idx(F, W, R, P, T, w, p, t - 1);
          values[nnz] = 1.0;
          nnz++;
        }
        /* sum_f y[f,w,p,t] */
        for (int f = 0; f < F; f++) {
          if (has_fw_arc(inst, f, w)) {
            indices[nnz] = yidx(F, W, P, T, f, w, p, t);
            values[nnz] = 1.0;
            nnz++;
          }
        }
        /* -s_w[w,p,t] */
        indices[nnz] = sw_idx(F, W, R, P, T, w, p, t);
        values[nnz] = -1.0;
        nnz++;
        /* -sum_r z[w,r,p,t] */
        for (int r = 0; r < R; r++) {
          if (has_wr_arc(inst, w, r)) {
            indices[nnz] = zidx(F, W, R, P, T, w, r, p, t);
            values[nnz] = -1.0;
            nnz++;
          }
        }
        Cbc_addRow(model, "row", nnz, indices, values, 'E', 0.0);
      }
    }
  }

  /* Demand satisfaction: sum_w z[w,r,p,t] = demand[r,p,t] */
  for (int r = 0; r < R; r++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        int nnz = 0;
        for (int w = 0; w < W; w++) {
          if (has_wr_arc(inst, w, r)) {
            indices[nnz] = zidx(F, W, R, P, T, w, r, p, t);
            values[nnz] = 1.0;
            nnz++;
          }
        }
        double rhs = inst->demand[(r * P + p) * T + t];
        Cbc_addRow(model, "row", nnz, indices, values, 'E', rhs);
      }
    }
  }

  /* Factory capacity: sum_p x[f,p,t] <= cap_factory[f] */
  for (int f = 0; f < F; f++) {
    for (int t = 0; t < T; t++) {
      int nnz = 0;
      for (int p = 0; p < P; p++) {
        indices[nnz] = xidx(F, P, T, f, p, t);
        values[nnz] = 1.0;
        nnz++;
      }
      Cbc_addRow(model, "row", nnz, indices, values, 'L', inst->cap_factory[f]);
    }
  }

  /* Warehouse capacity: sum_p s_w[w,p,t] <= cap_warehouse[w] */
  for (int w = 0; w < W; w++) {
    for (int t = 0; t < T; t++) {
      int nnz = 0;
      for (int p = 0; p < P; p++) {
        indices[nnz] = sw_idx(F, W, R, P, T, w, p, t);
        values[nnz] = 1.0;
        nnz++;
      }
      Cbc_addRow(model, "row", nnz, indices, values, 'L', inst->cap_warehouse[w]);
    }
  }

  /* Batch production: x[f,p,t] = batch_size[p] * n[f,p,t] */
  for (int f = 0; f < F; f++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        int nnz = 0;
        indices[nnz] = xidx(F, P, T, f, p, t);
        values[nnz] = 1.0;
        nnz++;
        indices[nnz] = nidx(F, W, R, P, T, f, p, t);
        values[nnz] = -(double)inst->batch_size[p];
        nnz++;
        Cbc_addRow(model, "row", nnz, indices, values, 'E', 0.0);
      }
    }
  }

  /* Setup link: n[f,p,t] <= M_f * u[f,p,t] */
  for (int f = 0; f < F; f++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        int M_f = inst->cap_factory[f] / inst->batch_size[p];
        int nnz = 0;
        indices[nnz] = nidx(F, W, R, P, T, f, p, t);
        values[nnz] = 1.0;
        nnz++;
        indices[nnz] = uidx(F, W, R, P, T, f, p, t);
        values[nnz] = -(double)M_f;
        nnz++;
        Cbc_addRow(model, "row", nnz, indices, values, 'L', 0.0);
      }
    }
  }

  /* Truck shipments: sum_r z[w,r,p,t] = truck_size * m[w,p,t] */
  for (int w = 0; w < W; w++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        int nnz = 0;
        for (int r = 0; r < R; r++) {
          if (has_wr_arc(inst, w, r)) {
            indices[nnz] = zidx(F, W, R, P, T, w, r, p, t);
            values[nnz] = 1.0;
            nnz++;
          }
        }
        indices[nnz] = midx(F, W, R, P, T, w, p, t);
        values[nnz] = -(double)inst->truck_size;
        nnz++;
        Cbc_addRow(model, "row", nnz, indices, values, 'E', 0.0);
      }
    }
  }

  /* Truck setup link: m[w,p,t] <= M_w * v[w,p,t] */
  for (int w = 0; w < W; w++) {
    for (int p = 0; p < P; p++) {
      for (int t = 0; t < T; t++) {
        int M_w = (inst->cap_warehouse[w] / inst->batch_size[p]) > 0
                      ? (inst->cap_warehouse[w] / inst->batch_size[p])
                      : 1;
        int nnz = 0;
        indices[nnz] = midx(F, W, R, P, T, w, p, t);
        values[nnz] = 1.0;
        nnz++;
        indices[nnz] = vidx(F, W, R, P, T, w, p, t);
        values[nnz] = -(double)M_w;
        nnz++;
        Cbc_addRow(model, "row", nnz, indices, values, 'L', 0.0);
      }
    }
  }

  free(indices);
  free(values);

  return model;
}

#endif /* MESCN_BUILDER_H */
