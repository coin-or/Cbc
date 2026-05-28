/* Multi-Item Capacitated Lot Sizing with Batch Production (MICLSP)
 * test fixtures.  All three variable types: binary (setup indicator),
 * integer (production batches), continuous (inventory).
 * Big-M links batch count to setup: n[i,t] <= floor(C[t]/B[i]) * y[i,t].
 * Capacity formula: C[t] = ceil(sum_i ceil(d[i,t]/B[i])*B[i] * rho).
 *
 * Do not edit manually; regenerate with test/gen_miclsp_fixtures.py.
 */

#ifndef MICLSP_FIXTURES_H
#define MICLSP_FIXTURES_H

typedef struct {
  int I;       /* number of items */
  int T;       /* number of periods */
  double rho;  /* capacity tightness factor */
  int opt;     /* certified MIP optimum */
  const int *B;  /* batch sizes [I] */
  const int *d;  /* demands [I*T], row-major: d[i*T+t] */
  const int *f;  /* setup costs [I] */
  const int *h;  /* holding costs [I] */
  const int *C;  /* machine capacity [T] */
} MiclspInstance;

/* miclsp_i04_t10_s042_r110: I=4 T=10 seed=42 rho=1.10 opt=6538 */
static const int B_miclsp_i04_t10_s042_r110[4] = {
    13, 4, 3, 14
};
static const int d_miclsp_i04_t10_s042_r110[40] = {
    9, 8, 8, 5, 4, 18, 3, 19, 14, 2, 1, 3,
    7, 8, 17, 20, 1, 18, 7, 18, 14, 8, 15, 19,
    9, 1, 6, 14, 11, 9, 5, 7, 11, 4, 3, 13,
    4, 12, 12, 20
};
static const int f_miclsp_i04_t10_s042_r110[4] = {
    235, 122, 473, 335
};
static const int h_miclsp_i04_t10_s042_r110[4] = {
    9, 2, 7, 2
};
static const int C_miclsp_i04_t10_s042_r110[10] = {
    51, 44, 56, 62, 62, 70, 41, 83, 66, 77
};

/* miclsp_i04_t10_s042_r114: I=4 T=10 seed=42 rho=1.15 opt=6490 */
static const int B_miclsp_i04_t10_s042_r114[4] = {
    13, 4, 3, 14
};
static const int d_miclsp_i04_t10_s042_r114[40] = {
    9, 8, 8, 5, 4, 18, 3, 19, 14, 2, 1, 3,
    7, 8, 17, 20, 1, 18, 7, 18, 14, 8, 15, 19,
    9, 1, 6, 14, 11, 9, 5, 7, 11, 4, 3, 13,
    4, 12, 12, 20
};
static const int f_miclsp_i04_t10_s042_r114[4] = {
    235, 122, 473, 335
};
static const int h_miclsp_i04_t10_s042_r114[4] = {
    9, 2, 7, 2
};
static const int C_miclsp_i04_t10_s042_r114[10] = {
    53, 46, 58, 65, 65, 73, 43, 87, 69, 81
};

/* miclsp_i04_t10_s042_r120: I=4 T=10 seed=42 rho=1.20 opt=6395 */
static const int B_miclsp_i04_t10_s042_r120[4] = {
    13, 4, 3, 14
};
static const int d_miclsp_i04_t10_s042_r120[40] = {
    9, 8, 8, 5, 4, 18, 3, 19, 14, 2, 1, 3,
    7, 8, 17, 20, 1, 18, 7, 18, 14, 8, 15, 19,
    9, 1, 6, 14, 11, 9, 5, 7, 11, 4, 3, 13,
    4, 12, 12, 20
};
static const int f_miclsp_i04_t10_s042_r120[4] = {
    235, 122, 473, 335
};
static const int h_miclsp_i04_t10_s042_r120[4] = {
    9, 2, 7, 2
};
static const int C_miclsp_i04_t10_s042_r120[10] = {
    56, 48, 60, 68, 68, 76, 45, 90, 72, 84
};

/* miclsp_i04_t10_s042_r130: I=4 T=10 seed=42 rho=1.30 opt=6238 */
static const int B_miclsp_i04_t10_s042_r130[4] = {
    13, 4, 3, 14
};
static const int d_miclsp_i04_t10_s042_r130[40] = {
    9, 8, 8, 5, 4, 18, 3, 19, 14, 2, 1, 3,
    7, 8, 17, 20, 1, 18, 7, 18, 14, 8, 15, 19,
    9, 1, 6, 14, 11, 9, 5, 7, 11, 4, 3, 13,
    4, 12, 12, 20
};
static const int f_miclsp_i04_t10_s042_r130[4] = {
    235, 122, 473, 335
};
static const int h_miclsp_i04_t10_s042_r130[4] = {
    9, 2, 7, 2
};
static const int C_miclsp_i04_t10_s042_r130[10] = {
    60, 52, 65, 73, 73, 82, 49, 98, 78, 91
};

/* miclsp_i04_t10_s137_r110: I=4 T=10 seed=137 rho=1.10 opt=8423 */
static const int B_miclsp_i04_t10_s137_r110[4] = {
    4, 5, 9, 8
};
static const int d_miclsp_i04_t10_s137_r110[40] = {
    7, 8, 16, 12, 1, 13, 17, 11, 13, 17, 16, 5,
    14, 17, 17, 11, 1, 20, 1, 5, 16, 12, 6, 12,
    8, 16, 13, 11, 11, 1, 13, 11, 3, 12, 17, 4,
    7, 2, 20, 11
};
static const int f_miclsp_i04_t10_s137_r110[4] = {
    378, 388, 381, 403
};
static const int h_miclsp_i04_t10_s137_r110[4] = {
    6, 1, 4, 10
};
static const int C_miclsp_i04_t10_s137_r110[10] = {
    69, 52, 53, 73, 63, 63, 57, 64, 70, 56
};

/* miclsp_i04_t10_s137_r114: I=4 T=10 seed=137 rho=1.15 opt=8302 */
static const int B_miclsp_i04_t10_s137_r114[4] = {
    4, 5, 9, 8
};
static const int d_miclsp_i04_t10_s137_r114[40] = {
    7, 8, 16, 12, 1, 13, 17, 11, 13, 17, 16, 5,
    14, 17, 17, 11, 1, 20, 1, 5, 16, 12, 6, 12,
    8, 16, 13, 11, 11, 1, 13, 11, 3, 12, 17, 4,
    7, 2, 20, 11
};
static const int f_miclsp_i04_t10_s137_r114[4] = {
    378, 388, 381, 403
};
static const int h_miclsp_i04_t10_s137_r114[4] = {
    6, 1, 4, 10
};
static const int C_miclsp_i04_t10_s137_r114[10] = {
    72, 55, 56, 76, 66, 66, 59, 67, 73, 58
};

/* miclsp_i04_t10_s137_r120: I=4 T=10 seed=137 rho=1.20 opt=7962 */
static const int B_miclsp_i04_t10_s137_r120[4] = {
    4, 5, 9, 8
};
static const int d_miclsp_i04_t10_s137_r120[40] = {
    7, 8, 16, 12, 1, 13, 17, 11, 13, 17, 16, 5,
    14, 17, 17, 11, 1, 20, 1, 5, 16, 12, 6, 12,
    8, 16, 13, 11, 11, 1, 13, 11, 3, 12, 17, 4,
    7, 2, 20, 11
};
static const int f_miclsp_i04_t10_s137_r120[4] = {
    378, 388, 381, 403
};
static const int h_miclsp_i04_t10_s137_r120[4] = {
    6, 1, 4, 10
};
static const int C_miclsp_i04_t10_s137_r120[10] = {
    75, 57, 58, 80, 69, 69, 62, 70, 76, 60
};

/* miclsp_i04_t12_s042_r110: I=4 T=12 seed=42 rho=1.10 opt=8414 */
static const int B_miclsp_i04_t12_s042_r110[4] = {
    13, 4, 3, 14
};
static const int d_miclsp_i04_t12_s042_r110[48] = {
    9, 8, 8, 5, 4, 18, 3, 19, 14, 2, 1, 3,
    7, 8, 17, 20, 1, 18, 7, 18, 14, 8, 15, 19,
    9, 1, 6, 14, 11, 9, 5, 7, 11, 4, 3, 13,
    4, 12, 12, 20, 9, 2, 15, 18, 4, 13, 3, 18,
};
static const int f_miclsp_i04_t12_s042_r110[4] = {
    250, 421, 416, 285
};
static const int h_miclsp_i04_t12_s042_r110[4] = {
    10, 4, 2, 1
};
static const int C_miclsp_i04_t12_s042_r110[12] = {
    49, 42, 59, 84, 48, 76, 61, 92, 75, 46, 51, 84,
};

/* miclsp_i04_t12_s042_r114: I=4 T=12 seed=42 rho=1.15 opt=8316 */
static const int B_miclsp_i04_t12_s042_r114[4] = {
    13, 4, 3, 14
};
static const int d_miclsp_i04_t12_s042_r114[48] = {
    9, 8, 8, 5, 4, 18, 3, 19, 14, 2, 1, 3,
    7, 8, 17, 20, 1, 18, 7, 18, 14, 8, 15, 19,
    9, 1, 6, 14, 11, 9, 5, 7, 11, 4, 3, 13,
    4, 12, 12, 20, 9, 2, 15, 18, 4, 13, 3, 18,
};
static const int f_miclsp_i04_t12_s042_r114[4] = {
    250, 421, 416, 285
};
static const int h_miclsp_i04_t12_s042_r114[4] = {
    10, 4, 2, 1
};
static const int C_miclsp_i04_t12_s042_r114[12] = {
    51, 44, 61, 88, 50, 80, 64, 96, 79, 48, 53, 88,
};

/* miclsp_i04_t12_s042_r120: I=4 T=12 seed=42 rho=1.20 opt=8181 */
static const int B_miclsp_i04_t12_s042_r120[4] = {
    13, 4, 3, 14
};
static const int d_miclsp_i04_t12_s042_r120[48] = {
    9, 8, 8, 5, 4, 18, 3, 19, 14, 2, 1, 3,
    7, 8, 17, 20, 1, 18, 7, 18, 14, 8, 15, 19,
    9, 1, 6, 14, 11, 9, 5, 7, 11, 4, 3, 13,
    4, 12, 12, 20, 9, 2, 15, 18, 4, 13, 3, 18,
};
static const int f_miclsp_i04_t12_s042_r120[4] = {
    250, 421, 416, 285
};
static const int h_miclsp_i04_t12_s042_r120[4] = {
    10, 4, 2, 1
};
static const int C_miclsp_i04_t12_s042_r120[12] = {
    53, 46, 64, 92, 52, 83, 66, 100, 82, 50, 56, 92,
};

/* miclsp_i04_t12_s137_r114: I=4 T=12 seed=137 rho=1.15 opt=8555 */
static const int B_miclsp_i04_t12_s137_r114[4] = {
    4, 5, 9, 8
};
static const int d_miclsp_i04_t12_s137_r114[48] = {
    7, 8, 16, 12, 1, 13, 17, 11, 13, 17, 16, 5,
    14, 17, 17, 11, 1, 20, 1, 5, 16, 12, 6, 12,
    8, 16, 13, 11, 11, 1, 13, 11, 3, 12, 17, 4,
    7, 2, 20, 11, 18, 19, 18, 19, 11, 2, 8, 19,
};
static const int f_miclsp_i04_t12_s137_r114[4] = {
    150, 296, 364, 439
};
static const int h_miclsp_i04_t12_s137_r114[4] = {
    1, 4, 6, 4
};
static const int C_miclsp_i04_t12_s137_r114[12] = {
    46, 63, 90, 71, 59, 80, 78, 68, 71, 71, 60, 65,
};

#define N_MICLSP_INSTANCES 11

static const MiclspInstance MICLSP_INSTANCES[] = {
  { 4, 10, 1.10, 6538,
    B_miclsp_i04_t10_s042_r110, d_miclsp_i04_t10_s042_r110, f_miclsp_i04_t10_s042_r110, h_miclsp_i04_t10_s042_r110, C_miclsp_i04_t10_s042_r110 },
  { 4, 10, 1.15, 6490,
    B_miclsp_i04_t10_s042_r114, d_miclsp_i04_t10_s042_r114, f_miclsp_i04_t10_s042_r114, h_miclsp_i04_t10_s042_r114, C_miclsp_i04_t10_s042_r114 },
  { 4, 10, 1.20, 6395,
    B_miclsp_i04_t10_s042_r120, d_miclsp_i04_t10_s042_r120, f_miclsp_i04_t10_s042_r120, h_miclsp_i04_t10_s042_r120, C_miclsp_i04_t10_s042_r120 },
  { 4, 10, 1.30, 6238,
    B_miclsp_i04_t10_s042_r130, d_miclsp_i04_t10_s042_r130, f_miclsp_i04_t10_s042_r130, h_miclsp_i04_t10_s042_r130, C_miclsp_i04_t10_s042_r130 },
  { 4, 10, 1.10, 8423,
    B_miclsp_i04_t10_s137_r110, d_miclsp_i04_t10_s137_r110, f_miclsp_i04_t10_s137_r110, h_miclsp_i04_t10_s137_r110, C_miclsp_i04_t10_s137_r110 },
  { 4, 10, 1.15, 8302,
    B_miclsp_i04_t10_s137_r114, d_miclsp_i04_t10_s137_r114, f_miclsp_i04_t10_s137_r114, h_miclsp_i04_t10_s137_r114, C_miclsp_i04_t10_s137_r114 },
  { 4, 10, 1.20, 7962,
    B_miclsp_i04_t10_s137_r120, d_miclsp_i04_t10_s137_r120, f_miclsp_i04_t10_s137_r120, h_miclsp_i04_t10_s137_r120, C_miclsp_i04_t10_s137_r120 },
  { 4, 12, 1.10, 8414,
    B_miclsp_i04_t12_s042_r110, d_miclsp_i04_t12_s042_r110, f_miclsp_i04_t12_s042_r110, h_miclsp_i04_t12_s042_r110, C_miclsp_i04_t12_s042_r110 },
  { 4, 12, 1.15, 8316,
    B_miclsp_i04_t12_s042_r114, d_miclsp_i04_t12_s042_r114, f_miclsp_i04_t12_s042_r114, h_miclsp_i04_t12_s042_r114, C_miclsp_i04_t12_s042_r114 },
  { 4, 12, 1.20, 8181,
    B_miclsp_i04_t12_s042_r120, d_miclsp_i04_t12_s042_r120, f_miclsp_i04_t12_s042_r120, h_miclsp_i04_t12_s042_r120, C_miclsp_i04_t12_s042_r120 },
  { 4, 12, 1.15, 8555,
    B_miclsp_i04_t12_s137_r114, d_miclsp_i04_t12_s137_r114, f_miclsp_i04_t12_s137_r114, h_miclsp_i04_t12_s137_r114, C_miclsp_i04_t12_s137_r114 },
};

#endif /* MICLSP_FIXTURES_H */
