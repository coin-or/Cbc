/* MESCN test fixtures with reference solutions.
 * 
 * Auto-generated from gen_mescn_fixtures_ci_final.py
 * DO NOT EDIT MANUALLY
 */

#ifndef MESCN_FIXTURES_H
#define MESCN_FIXTURES_H
#include <stddef.h>

typedef struct {
  const char *name;
  int certified_optimal;  /* -1 for infeasible, >= 0 for optimal */
  
  /* Dimensions */
  int F, W, R, P, T;
  
  /* Cost arrays */
  const int *setup_cost_f;    /* F*P */
  const int *setup_cost_w;    /* W*P */
  const int *prod_cost;       /* F*P */
  const int *transport_fw;    /* F*W*P */
  const int *transport_wr;    /* W*R*P */
  const int *holding_f;       /* F*P */
  const int *holding_w;       /* W*P */
  
  /* Capacities */
  const int *cap_factory;     /* F */
  const int *cap_warehouse;   /* W */
  
  /* Batch and truck sizes */
  const int *batch_size;      /* P */
  int truck_size;
  
  /* Demand */
  const int *demand;          /* R*P*T */
  
  /* Network structure */
  int fully_connected;
  int n_fw_arcs;
  const int *fw_arcs;         /* 2*n_fw_arcs: (f,w) pairs */
  int n_wr_arcs;
  const int *wr_arcs;         /* 2*n_wr_arcs: (w,r) pairs */
  
  /* Reference solution for debugging */
  const double *reference_solution;
  int reference_solution_size;
} MESCNInstance;


/* ========== tiny_feasible ========== */

static const int TINY_FEASIBLE_SETUP_COST_F[] = {
   100,
};

static const int TINY_FEASIBLE_SETUP_COST_W[] = {
    50,
};

static const int TINY_FEASIBLE_PROD_COST[] = {
     5,
};

static const int TINY_FEASIBLE_TRANSPORT_FW[] = {
     2,
};

static const int TINY_FEASIBLE_TRANSPORT_WR[] = {
     1,
};

static const int TINY_FEASIBLE_HOLDING_F[] = {
    10,
};

static const int TINY_FEASIBLE_HOLDING_W[] = {
    10,
};

static const int TINY_FEASIBLE_CAP_FACTORY[] = {
   500,
};

static const int TINY_FEASIBLE_CAP_WAREHOUSE[] = {
   500,
};

static const int TINY_FEASIBLE_BATCH_SIZE[] = {
    10,
};

static const int TINY_FEASIBLE_DEMAND[] = {
    10,   20,
};

static const double TINY_FEASIBLE_REFERENCE_SOL[] = {
      0.000000,     0.000000,     0.000000,     0.000000,    -9.000000,     1.000000,
    100.000000,   100.000000,    50.000000,    50.000000,   170.000000,    70.000000,
      0.000000,     0.000000,     0.000000,     8.000000,     0.000000,    10.000000,
    540.000000,     0.000000,     0.000000,     0.000000,     0.000000,    -9.000000,
      1.000000,   100.000000,   100.000000,    50.000000,    50.000000,   170.000000,
     70.000000,     0.000000,     0.000000,     0.000000,     8.000000,     0.000000,
     10.000000,
};


/* ========== infeas_demand_spike ========== */

static const int INFEAS_DEMAND_SPIKE_SETUP_COST_F[] = {
   100,  100,
};

static const int INFEAS_DEMAND_SPIKE_SETUP_COST_W[] = {
    50,   50,
};

static const int INFEAS_DEMAND_SPIKE_PROD_COST[] = {
     5,    5,
};

static const int INFEAS_DEMAND_SPIKE_TRANSPORT_FW[] = {
     2,    2,    2,    2,
};

static const int INFEAS_DEMAND_SPIKE_TRANSPORT_WR[] = {
     1,    1,    1,    1,
};

static const int INFEAS_DEMAND_SPIKE_HOLDING_F[] = {
     2,    2,
};

static const int INFEAS_DEMAND_SPIKE_HOLDING_W[] = {
     2,    2,
};

static const int INFEAS_DEMAND_SPIKE_CAP_FACTORY[] = {
   100,  100,
};

static const int INFEAS_DEMAND_SPIKE_CAP_WAREHOUSE[] = {
   200,  200,
};

static const int INFEAS_DEMAND_SPIKE_BATCH_SIZE[] = {
    10,
};

static const int INFEAS_DEMAND_SPIKE_DEMAND[] = {
    10,  150,   10,   10,  150,   10,
};

static const double INFEAS_DEMAND_SPIKE_REFERENCE_SOL[] = {
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,   100.000000,   100.000000,
    100.000000,   100.000000,    50.000000,    50.000000,    50.000000,    50.000000,
     70.000000,    50.000000,    70.000000,    50.000000,    20.000000,     0.000000,
     20.000000,    40.000000,     0.000000,     0.000000,     0.000000,  3600.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,   100.000000,   100.000000,
    100.000000,   100.000000,    50.000000,    50.000000,    50.000000,    50.000000,
     70.000000,    50.000000,    70.000000,    50.000000,    20.000000,     0.000000,
     20.000000,    40.000000,     0.000000,     0.000000,     0.000000,
};


/* ========== infeas_disconnected ========== */

static const int INFEAS_DISCONNECTED_SETUP_COST_F[] = {
   100,  100,
};

static const int INFEAS_DISCONNECTED_SETUP_COST_W[] = {
    50,   50,
};

static const int INFEAS_DISCONNECTED_PROD_COST[] = {
     5,    5,
};

static const int INFEAS_DISCONNECTED_TRANSPORT_FW[] = {
     2,    2,    2,    2,
};

static const int INFEAS_DISCONNECTED_TRANSPORT_WR[] = {
     1,    1,    1,    1,
};

static const int INFEAS_DISCONNECTED_HOLDING_F[] = {
     2,    2,
};

static const int INFEAS_DISCONNECTED_HOLDING_W[] = {
     2,    2,
};

static const int INFEAS_DISCONNECTED_CAP_FACTORY[] = {
   200,  200,
};

static const int INFEAS_DISCONNECTED_CAP_WAREHOUSE[] = {
   200,  200,
};

static const int INFEAS_DISCONNECTED_BATCH_SIZE[] = {
    10,
};

static const int INFEAS_DISCONNECTED_DEMAND[] = {
    30,   30,   30,   30,
};


/* ========== moderate_sparse_s7000 ========== */

static const int MODERATE_SPARSE_S7000_SETUP_COST_F[] = {
    70,   52,   71,   86,   64,   76,
};

static const int MODERATE_SPARSE_S7000_SETUP_COST_W[] = {
    68,   38,   70,   43,   52,   64,   37,   39,
};

static const int MODERATE_SPARSE_S7000_PROD_COST[] = {
    14,   18,   24,   15,   19,   26,
};

static const int MODERATE_SPARSE_S7000_TRANSPORT_FW[] = {
    33,   34,   17,   24,   36,   21,   30,   31,
    39,   16,   21,   15,   20,   19,   15,   40,
    19,   18,   33,   39,   32,   22,   20,   27,
};

static const int MODERATE_SPARSE_S7000_TRANSPORT_WR[] = {
    21,   14,   20,   10,   23,   11,   12,   22,
    22,   10,   26,   18,   18,   22,   13,   17,
    29,   13,   10,   27,   19,   18,   25,   14,
    12,   21,   25,   12,   23,   15,   18,   15,
    17,   26,   20,   30,   16,   21,   30,   18,
};

static const int MODERATE_SPARSE_S7000_HOLDING_F[] = {
     5,   10,    8,    5,    7,    4,
};

static const int MODERATE_SPARSE_S7000_HOLDING_W[] = {
     8,    9,    7,    5,    7,   10,    5,    8,
};

static const int MODERATE_SPARSE_S7000_CAP_FACTORY[] = {
   357,  306,  358,
};

static const int MODERATE_SPARSE_S7000_CAP_WAREHOUSE[] = {
   202,  222,  206,  264,
};

static const int MODERATE_SPARSE_S7000_BATCH_SIZE[] = {
     7,    7,
};

static const int MODERATE_SPARSE_S7000_DEMAND[] = {
    76,   57,   76,   57,  114,   95,   76,   76,
   114,   95,  133,  133,   95,  114,   95,   76,
    95,   76,   95,   57,   76,  114,   95,   57,
    95,   76,   57,   76,   95,   57,
};

static const int MODERATE_SPARSE_S7000_FW_ARCS[] = {
  0, 1,
  0, 3,
  1, 2,
  1, 1,
  2, 1,
  2, 2,
};

static const int MODERATE_SPARSE_S7000_WR_ARCS[] = {
  0, 1,
  0, 3,
  0, 4,
  1, 1,
  1, 0,
  1, 4,
  2, 0,
  2, 1,
  2, 3,
  3, 0,
  3, 3,
  3, 1,
  1, 2,
};

static const double MODERATE_SPARSE_S7000_REFERENCE_SOL[] = {
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,    -5.000000,
     -5.000000,    -5.000000,     9.000000,     9.000000,     9.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,    70.000000,
     70.000000,    70.000000,    86.000000,    86.000000,    86.000000,    64.000000,
     64.000000,    64.000000,    76.000000,    76.000000,    76.000000,    70.000000,
     70.000000,    70.000000,    43.000000,    43.000000,    43.000000,    52.000000,
     52.000000,    52.000000,    64.000000,    64.000000,    64.000000,    37.000000,
     37.000000,    37.000000,   105.000000,    70.000000,    35.000000,    28.000000,
     -7.000000,   112.000000,   252.000000,   217.000000,   182.000000,   126.000000,
    238.000000,   210.000000,   646.000000,   741.000000,   836.000000,   646.000000,
    741.000000,   418.000000,   494.000000,   589.000000,   684.000000,   570.000000,
    266.000000,   342.000000,   589.000000,   684.000000,   779.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
};


/* ========== moderate_sparse_s7001 ========== */

static const int MODERATE_SPARSE_S7001_SETUP_COST_F[] = {
    57,   95,   70,   92,   73,   94,
};

static const int MODERATE_SPARSE_S7001_SETUP_COST_W[] = {
    43,   58,   42,   57,   47,   66,   49,   44,
};

static const int MODERATE_SPARSE_S7001_PROD_COST[] = {
    21,   19,   21,   22,   28,   12,
};

static const int MODERATE_SPARSE_S7001_TRANSPORT_FW[] = {
    21,   33,   17,   19,   37,   33,   23,   15,
    15,   31,   36,   29,   27,   32,   17,   16,
    29,   16,   39,   19,   33,   24,   23,   31,
};

static const int MODERATE_SPARSE_S7001_TRANSPORT_WR[] = {
    22,   26,   18,   11,   26,   17,   19,   19,
    10,   29,   25,   14,   29,   16,   23,   18,
    19,   14,   16,   29,   24,   13,   22,   26,
    24,   28,   13,   17,   29,   15,   13,   20,
    10,   21,   20,   18,   25,   28,   30,   26,
};

static const int MODERATE_SPARSE_S7001_HOLDING_F[] = {
     7,   10,    4,    5,    9,    9,
};

static const int MODERATE_SPARSE_S7001_HOLDING_W[] = {
     7,    4,    4,    8,    9,    7,   10,    5,
};

static const int MODERATE_SPARSE_S7001_CAP_FACTORY[] = {
   416,  357,  418,
};

static const int MODERATE_SPARSE_S7001_CAP_WAREHOUSE[] = {
   261,  214,  299,  280,
};

static const int MODERATE_SPARSE_S7001_BATCH_SIZE[] = {
    11,   13,
};

static const int MODERATE_SPARSE_S7001_DEMAND[] = {
    95,   95,   95,  114,   76,  133,  133,   57,
   133,   57,   76,   57,  114,   76,  133,  133,
    95,   95,   95,  114,   76,   95,   76,   57,
    95,  114,   57,  133,   95,   95,
};

static const int MODERATE_SPARSE_S7001_FW_ARCS[] = {
  0, 3,
  0, 2,
  1, 3,
  1, 0,
  2, 1,
  2, 2,
};

static const int MODERATE_SPARSE_S7001_WR_ARCS[] = {
  0, 1,
  0, 2,
  0, 4,
  1, 0,
  1, 4,
  1, 2,
  2, 3,
  2, 4,
  2, 0,
  3, 2,
  3, 1,
  3, 4,
};

static const double MODERATE_SPARSE_S7001_REFERENCE_SOL[] = {
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,    14.000000,    14.000000,
     -6.000000,     2.000000,     2.000000,    31.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
     57.000000,    57.000000,    57.000000,    95.000000,    95.000000,    95.000000,
     70.000000,    70.000000,    70.000000,    94.000000,    94.000000,    94.000000,
     43.000000,    43.000000,    43.000000,    57.000000,    57.000000,    57.000000,
     47.000000,    47.000000,    47.000000,    66.000000,    66.000000,    66.000000,
     49.000000,    49.000000,    49.000000,    44.000000,    44.000000,    44.000000,
    396.000000,   352.000000,   308.000000,   507.000000,   442.000000,   377.000000,
    330.000000,   286.000000,   242.000000,   299.000000,   260.000000,   169.000000,
   -266.000000,  -190.000000,   266.000000,   646.000000,   513.000000,   342.000000,
    399.000000,   475.000000,   931.000000,   532.000000,   589.000000,   171.000000,
    152.000000,   228.000000,   190.000000,   399.000000,   304.000000,    38.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,
};


/* ========== moderate_sparse_s7002 ========== */

static const int MODERATE_SPARSE_S7002_SETUP_COST_F[] = {
    85,   83,   51,   54,   93,   87,
};

static const int MODERATE_SPARSE_S7002_SETUP_COST_W[] = {
    54,   63,   40,   35,   52,   38,   55,   67,
};

static const int MODERATE_SPARSE_S7002_PROD_COST[] = {
    22,   24,   20,   23,   20,   24,
};

static const int MODERATE_SPARSE_S7002_TRANSPORT_FW[] = {
    28,   17,   17,   38,   26,   28,   31,   20,
    34,   35,   20,   28,   36,   25,   34,   32,
    33,   26,   25,   39,   26,   32,   32,   23,
};

static const int MODERATE_SPARSE_S7002_TRANSPORT_WR[] = {
    22,   16,   27,   27,   24,   16,   21,   23,
    10,   23,   23,   12,   16,   12,   16,   28,
    23,   23,   24,   30,   30,   12,   18,   29,
    18,   30,   12,   29,   24,   29,   15,   19,
    29,   18,   28,   28,   19,   12,   11,   16,
};

static const int MODERATE_SPARSE_S7002_HOLDING_F[] = {
    10,   10,    4,    6,    9,   10,
};

static const int MODERATE_SPARSE_S7002_HOLDING_W[] = {
     8,    5,    4,    6,    4,    8,   10,    9,
};

static const int MODERATE_SPARSE_S7002_CAP_FACTORY[] = {
   379,  325,  381,
};

static const int MODERATE_SPARSE_S7002_CAP_WAREHOUSE[] = {
   173,  201,  218,  245,
};

static const int MODERATE_SPARSE_S7002_BATCH_SIZE[] = {
    13,   11,
};

static const int MODERATE_SPARSE_S7002_DEMAND[] = {
   133,  114,   95,   76,   76,   76,  133,   76,
    57,   76,  114,   76,   95,   76,  114,  133,
    57,  114,   57,   57,  114,   76,   95,   76,
   114,  133,   95,   76,  133,  114,
};

static const int MODERATE_SPARSE_S7002_FW_ARCS[] = {
  0, 3,
  0, 1,
  1, 1,
  1, 0,
  2, 3,
  2, 0,
};

static const int MODERATE_SPARSE_S7002_WR_ARCS[] = {
  0, 0,
  0, 2,
  0, 4,
  1, 2,
  1, 3,
  1, 1,
  2, 4,
  2, 0,
  2, 2,
  3, 2,
  3, 0,
  3, 3,
};

static const double MODERATE_SPARSE_S7002_REFERENCE_SOL[] = {
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,    12.000000,    12.000000,    12.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,    85.000000,    85.000000,    85.000000,
     83.000000,    83.000000,    83.000000,    51.000000,    51.000000,    51.000000,
     54.000000,    54.000000,    54.000000,    93.000000,    93.000000,    93.000000,
     87.000000,    87.000000,    87.000000,    54.000000,    54.000000,    54.000000,
     63.000000,    63.000000,    63.000000,    40.000000,    40.000000,    40.000000,
     35.000000,    35.000000,    35.000000,    55.000000,    55.000000,    55.000000,
     67.000000,    67.000000,    67.000000,   403.000000,   351.000000,   299.000000,
    462.000000,   396.000000,   374.000000,   416.000000,   364.000000,   312.000000,
    341.000000,   275.000000,   319.000000,   403.000000,   351.000000,   299.000000,
    484.000000,   429.000000,   374.000000,   418.000000,   494.000000,   570.000000,
   -114.000000,   -19.000000,    76.000000,   456.000000,   532.000000,   608.000000,
    380.000000,   494.000000,   418.000000,   399.000000,   475.000000,   551.000000,
     38.000000,   152.000000,   190.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,
};


/* ========== moderate_sparse_s7007 ========== */

static const int MODERATE_SPARSE_S7007_SETUP_COST_F[] = {
    88,   84,   73,   66,   84,   90,
};

static const int MODERATE_SPARSE_S7007_SETUP_COST_W[] = {
    40,   57,   44,   66,   62,   66,   45,   51,
};

static const int MODERATE_SPARSE_S7007_PROD_COST[] = {
    26,   20,   26,   21,   12,   26,
};

static const int MODERATE_SPARSE_S7007_TRANSPORT_FW[] = {
    39,   17,   31,   24,   38,   32,   23,   38,
    21,   37,   26,   30,   24,   36,   22,   33,
    30,   34,   19,   38,   21,   38,   29,   34,
};

static const int MODERATE_SPARSE_S7007_TRANSPORT_WR[] = {
    26,   12,   21,   19,   24,   14,   30,   19,
    26,   15,   25,   21,   27,   10,   10,   19,
    21,   21,   16,   28,   11,   25,   14,   10,
    13,   18,   30,   11,   28,   18,   28,   18,
    13,   16,   19,   11,   28,   22,   12,   16,
};

static const int MODERATE_SPARSE_S7007_HOLDING_F[] = {
     6,    8,    5,    6,    5,    4,
};

static const int MODERATE_SPARSE_S7007_HOLDING_W[] = {
     5,    9,    6,    5,    9,   10,   10,   10,
};

static const int MODERATE_SPARSE_S7007_CAP_FACTORY[] = {
   379,  325,  381,
};

static const int MODERATE_SPARSE_S7007_CAP_WAREHOUSE[] = {
   244,  195,  229,  217,
};

static const int MODERATE_SPARSE_S7007_BATCH_SIZE[] = {
     7,    7,
};

static const int MODERATE_SPARSE_S7007_DEMAND[] = {
    95,  133,  114,   76,   95,   76,   57,   57,
   114,   57,  114,   76,   57,   57,  114,   76,
   114,   57,  114,   76,   76,   57,  133,   95,
    76,   57,  114,  133,  133,   57,
};

static const int MODERATE_SPARSE_S7007_FW_ARCS[] = {
  0, 1,
  0, 0,
  1, 3,
  1, 1,
  2, 2,
  2, 0,
};

static const int MODERATE_SPARSE_S7007_WR_ARCS[] = {
  0, 4,
  0, 3,
  0, 1,
  1, 1,
  1, 3,
  1, 4,
  2, 2,
  2, 1,
  2, 0,
  3, 1,
  3, 2,
  3, 3,
};

static const double MODERATE_SPARSE_S7007_REFERENCE_SOL[] = {
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
     -3.000000,    -3.000000,    -3.000000,    25.000000,    25.000000,    25.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,    84.000000,    84.000000,    84.000000,
     73.000000,    73.000000,    73.000000,    66.000000,    66.000000,    66.000000,
     84.000000,    84.000000,    84.000000,    90.000000,    90.000000,    90.000000,
     40.000000,    57.000000,    57.000000,    57.000000,    44.000000,    44.000000,
     44.000000,    66.000000,    66.000000,    66.000000,    62.000000,    62.000000,
     62.000000,    66.000000,    66.000000,    66.000000,    51.000000,    51.000000,
     51.000000,   105.000000,    84.000000,    49.000000,   287.000000,   252.000000,
    217.000000,   154.000000,   133.000000,    91.000000,   189.000000,   154.000000,
    119.000000,   266.000000,   238.000000,   210.000000,   285.000000,   380.000000,
    437.000000,   532.000000,    38.000000,   304.000000,   399.000000,   551.000000,
    608.000000,   703.000000,   380.000000,   475.000000,   570.000000,   494.000000,
    570.000000,   646.000000,   608.000000,   665.000000,   779.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,     0.000000,     0.000000,
      0.000000,     0.000000,     0.000000,     0.000000,
};


/* ========== Instance array ========== */

static const MESCNInstance MESCN_FIXTURES[] = {
  { /* tiny_feasible */
    .name = "tiny_feasible",
    .certified_optimal = 540,
    .F = 1, .W = 1, .R = 1, .P = 1, .T = 2,
    .setup_cost_f = TINY_FEASIBLE_SETUP_COST_F,
    .setup_cost_w = TINY_FEASIBLE_SETUP_COST_W,
    .prod_cost = TINY_FEASIBLE_PROD_COST,
    .transport_fw = TINY_FEASIBLE_TRANSPORT_FW,
    .transport_wr = TINY_FEASIBLE_TRANSPORT_WR,
    .holding_f = TINY_FEASIBLE_HOLDING_F,
    .holding_w = TINY_FEASIBLE_HOLDING_W,
    .cap_factory = TINY_FEASIBLE_CAP_FACTORY,
    .cap_warehouse = TINY_FEASIBLE_CAP_WAREHOUSE,
    .batch_size = TINY_FEASIBLE_BATCH_SIZE,
    .truck_size = 10,
    .demand = TINY_FEASIBLE_DEMAND,
    .fully_connected = 1,
    .n_fw_arcs = 0,
    .fw_arcs = NULL,
    .n_wr_arcs = 0,
    .wr_arcs = NULL,
    .reference_solution = TINY_FEASIBLE_REFERENCE_SOL,
    .reference_solution_size = 37,
  },
  { /* infeas_demand_spike */
    .name = "infeas_demand_spike",
    .certified_optimal = 3600,
    .F = 2, .W = 2, .R = 2, .P = 1, .T = 3,
    .setup_cost_f = INFEAS_DEMAND_SPIKE_SETUP_COST_F,
    .setup_cost_w = INFEAS_DEMAND_SPIKE_SETUP_COST_W,
    .prod_cost = INFEAS_DEMAND_SPIKE_PROD_COST,
    .transport_fw = INFEAS_DEMAND_SPIKE_TRANSPORT_FW,
    .transport_wr = INFEAS_DEMAND_SPIKE_TRANSPORT_WR,
    .holding_f = INFEAS_DEMAND_SPIKE_HOLDING_F,
    .holding_w = INFEAS_DEMAND_SPIKE_HOLDING_W,
    .cap_factory = INFEAS_DEMAND_SPIKE_CAP_FACTORY,
    .cap_warehouse = INFEAS_DEMAND_SPIKE_CAP_WAREHOUSE,
    .batch_size = INFEAS_DEMAND_SPIKE_BATCH_SIZE,
    .truck_size = 10,
    .demand = INFEAS_DEMAND_SPIKE_DEMAND,
    .fully_connected = 1,
    .n_fw_arcs = 0,
    .fw_arcs = NULL,
    .n_wr_arcs = 0,
    .wr_arcs = NULL,
    .reference_solution = INFEAS_DEMAND_SPIKE_REFERENCE_SOL,
    .reference_solution_size = 71,
  },
  { /* infeas_disconnected */
    .name = "infeas_disconnected",
    .certified_optimal = -999,
    .F = 2, .W = 2, .R = 2, .P = 1, .T = 2,
    .setup_cost_f = INFEAS_DISCONNECTED_SETUP_COST_F,
    .setup_cost_w = INFEAS_DISCONNECTED_SETUP_COST_W,
    .prod_cost = INFEAS_DISCONNECTED_PROD_COST,
    .transport_fw = INFEAS_DISCONNECTED_TRANSPORT_FW,
    .transport_wr = INFEAS_DISCONNECTED_TRANSPORT_WR,
    .holding_f = INFEAS_DISCONNECTED_HOLDING_F,
    .holding_w = INFEAS_DISCONNECTED_HOLDING_W,
    .cap_factory = INFEAS_DISCONNECTED_CAP_FACTORY,
    .cap_warehouse = INFEAS_DISCONNECTED_CAP_WAREHOUSE,
    .batch_size = INFEAS_DISCONNECTED_BATCH_SIZE,
    .truck_size = 10,
    .demand = INFEAS_DISCONNECTED_DEMAND,
    .fully_connected = 1,
    .n_fw_arcs = 0,
    .fw_arcs = NULL,
    .n_wr_arcs = 0,
    .wr_arcs = NULL,
    .reference_solution = NULL,
    .reference_solution_size = 0,
  },
  { /* moderate_sparse_s7000 */
    .name = "moderate_sparse_s7000",
    .certified_optimal = 139005,
    .F = 3, .W = 4, .R = 5, .P = 2, .T = 3,
    .setup_cost_f = MODERATE_SPARSE_S7000_SETUP_COST_F,
    .setup_cost_w = MODERATE_SPARSE_S7000_SETUP_COST_W,
    .prod_cost = MODERATE_SPARSE_S7000_PROD_COST,
    .transport_fw = MODERATE_SPARSE_S7000_TRANSPORT_FW,
    .transport_wr = MODERATE_SPARSE_S7000_TRANSPORT_WR,
    .holding_f = MODERATE_SPARSE_S7000_HOLDING_F,
    .holding_w = MODERATE_SPARSE_S7000_HOLDING_W,
    .cap_factory = MODERATE_SPARSE_S7000_CAP_FACTORY,
    .cap_warehouse = MODERATE_SPARSE_S7000_CAP_WAREHOUSE,
    .batch_size = MODERATE_SPARSE_S7000_BATCH_SIZE,
    .truck_size = 19,
    .demand = MODERATE_SPARSE_S7000_DEMAND,
    .fully_connected = 0,
    .n_fw_arcs = 6,
    .fw_arcs = MODERATE_SPARSE_S7000_FW_ARCS,
    .n_wr_arcs = 13,
    .wr_arcs = MODERATE_SPARSE_S7000_WR_ARCS,
    .reference_solution = MODERATE_SPARSE_S7000_REFERENCE_SOL,
    .reference_solution_size = 125,
  },
  { /* moderate_sparse_s7001 */
    .name = "moderate_sparse_s7001",
    .certified_optimal = 161955,
    .F = 3, .W = 4, .R = 5, .P = 2, .T = 3,
    .setup_cost_f = MODERATE_SPARSE_S7001_SETUP_COST_F,
    .setup_cost_w = MODERATE_SPARSE_S7001_SETUP_COST_W,
    .prod_cost = MODERATE_SPARSE_S7001_PROD_COST,
    .transport_fw = MODERATE_SPARSE_S7001_TRANSPORT_FW,
    .transport_wr = MODERATE_SPARSE_S7001_TRANSPORT_WR,
    .holding_f = MODERATE_SPARSE_S7001_HOLDING_F,
    .holding_w = MODERATE_SPARSE_S7001_HOLDING_W,
    .cap_factory = MODERATE_SPARSE_S7001_CAP_FACTORY,
    .cap_warehouse = MODERATE_SPARSE_S7001_CAP_WAREHOUSE,
    .batch_size = MODERATE_SPARSE_S7001_BATCH_SIZE,
    .truck_size = 19,
    .demand = MODERATE_SPARSE_S7001_DEMAND,
    .fully_connected = 0,
    .n_fw_arcs = 6,
    .fw_arcs = MODERATE_SPARSE_S7001_FW_ARCS,
    .n_wr_arcs = 12,
    .wr_arcs = MODERATE_SPARSE_S7001_WR_ARCS,
    .reference_solution = MODERATE_SPARSE_S7001_REFERENCE_SOL,
    .reference_solution_size = 133,
  },
  { /* moderate_sparse_s7002 */
    .name = "moderate_sparse_s7002",
    .certified_optimal = 180576,
    .F = 3, .W = 4, .R = 5, .P = 2, .T = 3,
    .setup_cost_f = MODERATE_SPARSE_S7002_SETUP_COST_F,
    .setup_cost_w = MODERATE_SPARSE_S7002_SETUP_COST_W,
    .prod_cost = MODERATE_SPARSE_S7002_PROD_COST,
    .transport_fw = MODERATE_SPARSE_S7002_TRANSPORT_FW,
    .transport_wr = MODERATE_SPARSE_S7002_TRANSPORT_WR,
    .holding_f = MODERATE_SPARSE_S7002_HOLDING_F,
    .holding_w = MODERATE_SPARSE_S7002_HOLDING_W,
    .cap_factory = MODERATE_SPARSE_S7002_CAP_FACTORY,
    .cap_warehouse = MODERATE_SPARSE_S7002_CAP_WAREHOUSE,
    .batch_size = MODERATE_SPARSE_S7002_BATCH_SIZE,
    .truck_size = 19,
    .demand = MODERATE_SPARSE_S7002_DEMAND,
    .fully_connected = 0,
    .n_fw_arcs = 6,
    .fw_arcs = MODERATE_SPARSE_S7002_FW_ARCS,
    .n_wr_arcs = 12,
    .wr_arcs = MODERATE_SPARSE_S7002_WR_ARCS,
    .reference_solution = MODERATE_SPARSE_S7002_REFERENCE_SOL,
    .reference_solution_size = 154,
  },
  { /* moderate_sparse_s7007 */
    .name = "moderate_sparse_s7007",
    .certified_optimal = 159376,
    .F = 3, .W = 4, .R = 5, .P = 2, .T = 3,
    .setup_cost_f = MODERATE_SPARSE_S7007_SETUP_COST_F,
    .setup_cost_w = MODERATE_SPARSE_S7007_SETUP_COST_W,
    .prod_cost = MODERATE_SPARSE_S7007_PROD_COST,
    .transport_fw = MODERATE_SPARSE_S7007_TRANSPORT_FW,
    .transport_wr = MODERATE_SPARSE_S7007_TRANSPORT_WR,
    .holding_f = MODERATE_SPARSE_S7007_HOLDING_F,
    .holding_w = MODERATE_SPARSE_S7007_HOLDING_W,
    .cap_factory = MODERATE_SPARSE_S7007_CAP_FACTORY,
    .cap_warehouse = MODERATE_SPARSE_S7007_CAP_WAREHOUSE,
    .batch_size = MODERATE_SPARSE_S7007_BATCH_SIZE,
    .truck_size = 19,
    .demand = MODERATE_SPARSE_S7007_DEMAND,
    .fully_connected = 0,
    .n_fw_arcs = 6,
    .fw_arcs = MODERATE_SPARSE_S7007_FW_ARCS,
    .n_wr_arcs = 12,
    .wr_arcs = MODERATE_SPARSE_S7007_WR_ARCS,
    .reference_solution = MODERATE_SPARSE_S7007_REFERENCE_SOL,
    .reference_solution_size = 148,
  },
};

#define N_MESCN_FIXTURES 7

#endif /* MESCN_FIXTURES_H */
