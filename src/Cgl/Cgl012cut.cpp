// Copyright (C) 2010, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).
/** @file 012cut.c Definition file for C coded 0-1/2 separator */
#include "CoinFinite.hpp"
#include "CoinTime.hpp"
#include "Cgl012cut.hpp"
#include "CglZeroHalf.hpp"
#include <algorithm>
#include <cstdint>
#include <climits>
#include <cmath>
static const int MAX_CUTS = 10000000;
//#define CGL_ZH_ADVANCED_DEBUG_PRINT_TABU
//#define CGL_ZH_ADVANCED_DEBUG_PRINT_CUTS
//#define CGL_ZH_ADVANCED_DEBUG_PRINT_TIMING
//#define CGL_ZH_ADVANCED_DEBUG_TIMING

/* #define CGL_ZH_ADVANCED_DEBUG_TIMING  */
#undef CGL_ZH_ADVANCED_DEBUG_TIMING

#define TRUE 1
#define FALSE 0

#define ODD 1
#define EVEN 0

#define NONE -1
#define BOTH 2

#define IN 1
#define OUT 0
#define ADD 1
#define DEL 0

#define LOWER_BOUND 0
#define UPPER_BOUND 1

#define EPS 0.0001 /* small tolerance */
//#define EPS 0.000001 /* small tolerance */
#define ZERO 0.000001 /* estimated accuracy for doubles */
//#define ZERO 0.0001 /* estimated accuracy for doubles */
#define INF 1000000000.0
#define IINF 1000000000

#define MAX_SLACK 1.0
#define MAX_LOSS 1.0
#define MAX_CYCLE_WEIGHT 1.0
#define MIN_VIOLATION 0.001
#define MIN_SCORE_RANGE 10.0
#define MAX_SCORE_RANGE ZERO /* 1.0 */

#define ISCALE 10000

//#define MAX_CUTS  10

#define MAX_CUT_POOL 10000
#define MAX_CUT_COD 10000
#define MAX_ITER_POOL 100
#define CLEAN_THRESH 0.9
#define MANY_IT_ZERO 10

#define mod2(I) ( I % 2 == 0 ? 0 : 1 )

static const int MAX_SEP_GRAPH_ACTIVE_NODES = 8000;

static long long
cglZeroHalfSepGraphMaxEdges(int nnodes)
{
  return nnodes > 1 ? (static_cast<long long>(nnodes) * (static_cast<long long>(nnodes) - 1)) / 2 : 0;
}

static bool
cglZeroHalfSepGraphTooLarge(int nnodes, long long maxedges, double *estimatedGiB)
{
  const double bytes = static_cast<double>(maxedges) * 2.0 * static_cast<double>(sizeof(edge *));
  if (estimatedGiB)
    *estimatedGiB = bytes / (1024.0 * 1024.0 * 1024.0);

  if (nnodes > MAX_SEP_GRAPH_ACTIVE_NODES)
    return true;
  if (maxedges > static_cast<long long>(INT_MAX))
    return true;
  return false;
}

static bool
cglZeroHalfUseSparseSepGraph(int nnodes, long long maxedges, int sparseThreshold, double *estimatedGiB)
{
  if (cglZeroHalfSepGraphTooLarge(nnodes, maxedges, estimatedGiB))
    return true;
  if (sparseThreshold == 0)
    return true;
  if (sparseThreshold > 0 && nnodes > sparseThreshold)
    return true;
  return false;
}

static std::uint64_t
cglZeroHalfSparseEdgeKey(int j, int k, short int parity)
{
  const std::uint64_t endpoint1 = static_cast<std::uint64_t>(j < k ? j : k);
  const std::uint64_t endpoint2 = static_cast<std::uint64_t>(j < k ? k : j);
  return (endpoint1 << 33) | (endpoint2 << 1) | static_cast<std::uint64_t>(parity == ODD ? 1 : 0);
}

static edge *
cglZeroHalfGetAuxiliaryArcEdge(const auxiliary_graph *a_graph, int from, int to)
{
#ifndef CGL_NEW_SHORT
  for (arc *arcPtr = a_graph->nodes[from].first; arcPtr < a_graph->nodes[from + 1].first; ++arcPtr)
    if (arcPtr->head->index == to)
      return arcPtr->backEdge;
#else
  for (cgl_arc *arcPtr = a_graph->nodes[from].firstArc; arcPtr < a_graph->nodes[from + 1].firstArc; ++arcPtr)
    if (arcPtr->to == to)
      return arcPtr->backEdge;
#endif
  return NULL;
}

bool Cgl012Cut::checkTimeLimit(const char *, const char *)
{
  enum { TIME_CHECK_POLL_INTERVAL = 1024 };
  if (timeLimitReached_ || maxSeconds_ <= 0.0)
    return timeLimitReached_;
  if (--timeCheckCountdown_ > 0)
    return false;
  timeCheckCountdown_ = TIME_CHECK_POLL_INTERVAL;
  if (CoinGetTimeOfDay() - profileStartSeconds_ >= maxSeconds_)
    timeLimitReached_ = true;
  return timeLimitReached_;
}

bool Cgl012Cut::timeLimitReached() const
{
  return timeLimitReached_;
}

void Cgl012Cut::resetTimeCheckState()
{
  timeLimitReached_ = false;
  timeCheckCountdown_ = 1024;
}

void Cgl012Cut::ensureBestWeakeningBufferCapacity(int requiredSize)
{
  if (static_cast<int>(typeEvenWeakBuffer_.size()) < requiredSize) {
    typeEvenWeakBuffer_.resize(requiredSize);
    switchEvenWeakBuffer_.resize(requiredSize);
    typeOddWeakBuffer_.resize(requiredSize);
    switchOddWeakBuffer_.resize(requiredSize);
  }
}

void Cgl012Cut::ensureVarsToWeakBufferCapacity(int requiredSize)
{
  if (static_cast<int>(varsToWeakBuffer_.size()) < requiredSize)
    varsToWeakBuffer_.resize(requiredSize);
}

static void
cglZeroHalfCollectFractionality(const ilp *problem, int row, int histogram[10],
  int *fractionalCount, double *fractionalSum, int *veryFractionalCount)
{
  const int begin = problem->mtbeg[row];
  const int nz = problem->mtcnt[row];

  *fractionalCount = 0;
  *fractionalSum = 0.0;
  *veryFractionalCount = 0;
  for (int bucket = 0; bucket < 10; ++bucket)
    histogram[bucket] = 0;

  for (int offset = 0; offset < nz; ++offset) {
    const int column = problem->mtind[begin + offset];
    const double value = problem->xstar[column];
    double fractional = value - std::floor(value);
    if (fractional < ZERO || 1.0 - fractional < ZERO)
      continue;
    int bucket = static_cast<int>(std::floor(fractional * 10.0));
    if (bucket < 0)
      bucket = 0;
    else if (bucket > 9)
      bucket = 9;
    histogram[bucket]++;
    (*fractionalCount)++;
    *fractionalSum += fractional;
    if (fractional >= 0.3 - ZERO && fractional <= 0.7 + ZERO)
      (*veryFractionalCount)++;
  }
}

#ifdef CGL_ZH_ADVANCED_DEBUG_CONSTRAINT_PROFILING
struct CglZeroHalfConstraintProfileState {
  std::vector<double> rowElapsedSeconds;
  std::vector<int> rowPairCounts;
  std::vector<int> usefulCutCounts;
  std::vector<double> usefulViolationSums;
};

static std::unordered_map<const Cgl012Cut *, CglZeroHalfConstraintProfileState> &
cglZeroHalfConstraintProfileStates()
{
  static std::unordered_map<const Cgl012Cut *, CglZeroHalfConstraintProfileState> states;
  return states;
}

static int
cglZeroHalfProfileAbsInt(int value)
{
  return value < 0 ? -value : value;
}

static int
cglZeroHalfProfileGcd(int a, int b)
{
  a = cglZeroHalfProfileAbsInt(a);
  b = cglZeroHalfProfileAbsInt(b);
  while (b != 0) {
    const int next = a % b;
    a = b;
    b = next;
  }
  return a;
}

static bool
cglZeroHalfProfileIsBinaryVar(const ilp *problem, int column)
{
  return problem->vlb[column] == 0 && problem->vub[column] == 1;
}

static const char *
cglZeroHalfProfileRowType(const ilp *problem, int row)
{
  const int begin = problem->mtbeg[row];
  const int nz = problem->mtcnt[row];
  const char sense = problem->msense[row];
  const int rhs = problem->mrhs[row];

  int nBinary = 0;
  int nNegative = 0;
  int nPositive = 0;
  int minValue = 0;
  int maxValue = 0;
  bool allOnes = true;
  bool allIntegers = true;

  for (int offset = 0; offset < nz; ++offset) {
    const int column = problem->mtind[begin + offset];
    const int value = problem->mtval[begin + offset];
    if (offset == 0 || value < minValue)
      minValue = value;
    if (offset == 0 || value > maxValue)
      maxValue = value;
    if (value < 0)
      ++nNegative;
    else if (value > 0)
      ++nPositive;
    if (value != 1)
      allOnes = false;
    if (!cglZeroHalfProfileIsBinaryVar(problem, column))
      allIntegers = false;
    else
      ++nBinary;
  }

  if (nz == 1)
    return "singleton";
  if (nz == 2) {
    if (sense == 'E')
      return "aggregate";
    if (nBinary == 1)
      return "variable_bound";
    if ((nBinary % 2 == 0) && nNegative == 1 && nPositive == 1 && minValue == -maxValue)
      return "precedence";
  }

  if (nBinary == nz) {
    if (allOnes) {
      if (rhs == 1) {
        if (sense == 'E')
          return "partitioning";
        if (sense == 'L')
          return "packing";
        if (sense == 'G')
          return "covering";
      } else if (rhs >= 2) {
        if (sense == 'E')
          return "cardinality";
        if (sense == 'L')
          return "invariant_knapsack";
      }
    } else if (rhs >= 2) {
      if ((maxValue > minValue) && nNegative == 0)
        return allIntegers ? "integer_knapsack" : "knapsack";
      if (nNegative == 1 && nz >= 2)
        return "bin_packing";
    }

    if (nNegative >= 2 && nPositive >= 2 && sense == 'E')
      return "flow_binary";
  } else {
    if (nNegative >= 2 && nPositive >= 2 && sense == 'E')
      return "flow_mixed";
  }

  return "general";
}

static void
cglZeroHalfProfileWriteHeader(FILE *file)
{
  fprintf(file,
    "round,role,row,row_type,sense,rhs,row_nz,odd_nz,row_deleted,slack,row_elapsed_seconds,row_pair_count,"
    "useful_cut_count,useful_total_violation,positive_coeffs,negative_coeffs,min_coeff,max_coeff,coeff_gcd,"
    "fractional_count,fractional_mean,very_fractional_count,very_fractional_share,hist_0_0.1,hist_0.1_0.2,hist_0.2_0.3,hist_0.3_0.4,hist_0.4_0.5,"
    "hist_0.5_0.6,hist_0.6_0.7,hist_0.7_0.8,hist_0.8_0.9,hist_0.9_1.0,generated_cuts\n");
}

void Cgl012Cut::zhProfileStartRound()
{
  CglZeroHalfConstraintProfileState &state = cglZeroHalfConstraintProfileStates()[this];
  state.rowElapsedSeconds.assign(inp_ilp->mr, 0.0);
  state.rowPairCounts.assign(inp_ilp->mr, 0);
  state.usefulCutCounts.assign(inp_ilp->mr, 0);
  state.usefulViolationSums.assign(inp_ilp->mr, 0.0);
}

void Cgl012Cut::zhProfileRecordRow(int row, double elapsedSeconds, int pairCount)
{
  if (row < 0 || row >= inp_ilp->mr)
    return;
  CglZeroHalfConstraintProfileState &state = cglZeroHalfConstraintProfileStates()[this];
  state.rowElapsedSeconds[row] += elapsedSeconds;
  state.rowPairCounts[row] = pairCount;
}

void Cgl012Cut::zhProfileRecordCut(const cut *v_cut)
{
  if (v_cut == NULL)
    return;
  CglZeroHalfConstraintProfileState &state = cglZeroHalfConstraintProfileStates()[this];
  for (int i = 0; i < v_cut->n_of_constr; ++i) {
    const int row = v_cut->constr_list[i];
    if (row < 0 || row >= inp_ilp->mr)
      continue;
    state.usefulCutCounts[row]++;
    state.usefulViolationSums[row] += v_cut->violation;
  }
}

static void
cglZeroHalfProfileEmitRow(FILE *file, const ilp *problem, const parity_ilp *parity,
  const std::vector<double> &rowElapsedSeconds, const std::vector<int> &rowPairCounts,
  const std::vector<int> &usefulCutCounts, const std::vector<double> &usefulViolationSums,
  int round, const char *role, int row, int generatedCuts)
{
  const int begin = problem->mtbeg[row];
  const int nz = problem->mtcnt[row];

  int positiveCount = 0;
  int negativeCount = 0;
  int minCoeff = 0;
  int maxCoeff = 0;
  int coeffGcd = 0;
  for (int offset = 0; offset < nz; ++offset) {
    const int coeff = problem->mtval[begin + offset];
    if (offset == 0 || coeff < minCoeff)
      minCoeff = coeff;
    if (offset == 0 || coeff > maxCoeff)
      maxCoeff = coeff;
    if (coeff > 0)
      ++positiveCount;
    else if (coeff < 0)
      ++negativeCount;
    coeffGcd = (offset == 0) ? cglZeroHalfProfileAbsInt(coeff) : cglZeroHalfProfileGcd(coeffGcd, coeff);
  }

  int histogram[10];
  int fractionalCount = 0;
  int veryFractionalCount = 0;
  double fractionalSum = 0.0;
  cglZeroHalfCollectFractionality(problem, row, histogram, &fractionalCount, &fractionalSum,
    &veryFractionalCount);
  const double fractionalMean = fractionalCount ? fractionalSum / static_cast<double>(fractionalCount) : 0.0;
  const double veryFractionalShare = nz ? static_cast<double>(veryFractionalCount) / static_cast<double>(nz) : 0.0;

  fprintf(file,
    "%d,%s,%d,%s,%c,%d,%d,%d,%d,%.10g,%.10g,%d,%d,%.10g,%d,%d,%d,%d,%d,%d,%.10g,%d,%.10g,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
    round,
    role,
    row,
    cglZeroHalfProfileRowType(problem, row),
    problem->msense[row],
    problem->mrhs[row],
    nz,
    parity ? parity->mtcnt[row] : 0,
    parity ? parity->row_to_delete[row] : 0,
    parity ? parity->slack[row] : 0.0,
    rowElapsedSeconds[row],
    rowPairCounts[row],
    usefulCutCounts[row],
    usefulViolationSums[row],
    positiveCount,
    negativeCount,
    minCoeff,
    maxCoeff,
    coeffGcd,
    fractionalCount,
    fractionalMean,
    veryFractionalCount,
    veryFractionalShare,
    histogram[0],
    histogram[1],
    histogram[2],
    histogram[3],
    histogram[4],
    histogram[5],
    histogram[6],
    histogram[7],
    histogram[8],
    histogram[9],
    generatedCuts);
}

void Cgl012Cut::zhProfileFlushRound(int round, int generatedCuts)
{
  const char *path = std::getenv("CGL_ZH_CONSTRAINT_PROFILE_FILE");
  if (path == NULL || !path[0])
    return;

  const std::unordered_map<const Cgl012Cut *, CglZeroHalfConstraintProfileState>::const_iterator found =
    cglZeroHalfConstraintProfileStates().find(this);
  if (found == cglZeroHalfConstraintProfileStates().end())
    return;
  const CglZeroHalfConstraintProfileState &state = found->second;

  int expensiveRow = -1;
  int usefulRow = -1;
  for (int row = 0; row < inp_ilp->mr; ++row) {
    if (state.rowElapsedSeconds[row] > 0.0 &&
        (expensiveRow < 0 || state.rowElapsedSeconds[row] > state.rowElapsedSeconds[expensiveRow]))
      expensiveRow = row;

    if (state.usefulCutCounts[row] > 0 &&
        (usefulRow < 0 ||
         state.usefulCutCounts[row] > state.usefulCutCounts[usefulRow] ||
         (state.usefulCutCounts[row] == state.usefulCutCounts[usefulRow] &&
           state.usefulViolationSums[row] > state.usefulViolationSums[usefulRow])))
      usefulRow = row;
  }

  if (expensiveRow < 0 && usefulRow < 0)
    return;

  FILE *file = fopen(path, "a+");
  if (file == NULL)
    return;

  if (fseek(file, 0, SEEK_END) == 0 && ftell(file) == 0)
    cglZeroHalfProfileWriteHeader(file);

  if (expensiveRow >= 0)
    cglZeroHalfProfileEmitRow(file, inp_ilp, p_ilp, state.rowElapsedSeconds,
      state.rowPairCounts, state.usefulCutCounts, state.usefulViolationSums,
      round, "most_expensive", expensiveRow, generatedCuts);
  if (usefulRow >= 0)
    cglZeroHalfProfileEmitRow(file, inp_ilp, p_ilp, state.rowElapsedSeconds,
      state.rowPairCounts, state.usefulCutCounts, state.usefulViolationSums,
      round, "most_useful", usefulRow, generatedCuts);

  fclose(file);
}
#endif


#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING

static float tot_basic_sep_time = 0.0; /* total time spent for basic 
					  separation */
static float avg_basic_sep_time; /* average time per iteration spent for basic
				    separation */
static float total_time = 0.0; /* total time spent in the separation */
static float prep_time = 0.0; /* time spent for the definition of the
				 parity ILP data structure */
static float weak_time = 0.0; /* time spent for the construction of the
				 separation graph by weakening */
static float aux_time = 0.0; /* time spent for the definition of the 
				auxiliary graph */
static float path_time = 0.0; /* time spent in the computation of the
				 shortest paths */ 
static float cycle_time = 0.0; /* time spent in the determination of the
				  shortest cycles */ 
static float cut_time = 0.0; /* time spent in the determination of the
				violated cuts */
static float bw_time = 0.0; /* time spent in best_weakening */
static float coef_time = 0.0; /* time spent in the initial computation
				 of coef in get_cut */
static float pool_time = 0.0; /* time spent for the addition and
				 extraction of cuts from the pool */
static int cut_ncalls = 0; /* number of calls to get_cut */
static float tabu_time = 0.0; /* time spent within tabu search */
float ti, tf, td, tti, ttf, tsi, tsf, tii, tff, tpi, tpf, ttabi, ttabf;
void 
second_(float *t) {*t=CoinCpuTime();}
#endif

/* #endif */

/* global data structures */

#define CGGGGG
#ifndef CGGGGG
static ilp *inp_ilp; /* input ILP data structure */
static parity_ilp *p_ilp; /* parity ILP data structure */
#endif

#ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_CUTS
/* utility subroutines */

void print_int_vect(char *s,int *v,int n)
{
  int i;

  printf("integer vector %s:",s);
  for ( i = 0; i < n; i++ ) printf(" %d",v[i]);
  printf("\n");
}

void print_short_int_vect(char *s,short int *v,int n)
{
  int i;

  printf("short integer vector %s:",s);
  for ( i = 0; i < n; i++ ) printf(" %d",v[i]);
  printf("\n");
}

void print_double_vect(char *s,const double *v,int n)
{
  int i;

  printf("double vector %s:",s);
  for ( i = 0; i < n; i++ ) printf(" %f",v[i]);
  printf("\n");
}
#endif
void alloc_error(char *s)
{
  printf("\n Warning: Not enough memory to allocate %s\n",s);
  printf("\n Cannot proceed with 0-1/2 cut separation\n");
  exit(FALSE);
}

/* double2int: compute the integer equivalent of a double */

int double2int(double x)
{
  if ( x > IINF ) return (IINF);
  if ( x < - IINF ) return (- IINF);
  if ( x < ZERO && x > - ZERO ) return(0);
  if ( x > 0.0 ) return(static_cast<int> (x + ZERO));
  return(static_cast<int> (x - ZERO));
}

/* gcd: compute the greatest common divisor of two integers */

int gcd(int a,int b)
{
  int c;

  if ( a < 0 ) a = - a;
  if ( b < 0 ) b = - b;
  if ( a < b ) { c = a; a = b; b = c; }
  while ( b != 0 ) {
    c = a % b; a = b; b = c;
  }
  return(a);
}

/* ILP data structures subroutines */

/* ilp_load: load the input ILP into an internal data structure */

void Cgl012Cut::ilp_load(
	      int mr, /* number of rows in the ILP matrix */
	      int mc, /* number of columns in the ILP matrix */
	      int mnz, /* number of nonzero's in the ILP matrix */
	      int *mtbeg, /* starting position of each row in arrays mtind and mtval */
	      int *mtcnt, /* number of entries of each row in arrays mtind and mtval */
	      int *mtind, /* column indices of the nonzero entries of the ILP matrix */
	      int *mtval, /* values of the nonzero entries of the ILP matrix */
	      int *vlb, /* lower bounds on the variables */
	      int *vub, /* upper bounds on the variables */
	      int *mrhs, /* right hand sides of the constraints */
	      char *msense /* senses of the constraints: 'L', 'G' or 'E' */
	      )
{
  inp_ilp = reinterpret_cast<ilp *> (calloc(1,sizeof(ilp)));
  if ( inp_ilp == NULL ) alloc_error(const_cast<char*>("inp_ilp"));

  inp_ilp->mr = mr; inp_ilp->mc = mc; inp_ilp->mnz = mnz; 
  inp_ilp->mtbeg = mtbeg; inp_ilp->mtcnt = mtcnt; 
  inp_ilp->mtind = mtind; inp_ilp->mtval = mtval; 
  inp_ilp->vlb = vlb; inp_ilp->vub = vub; 
  inp_ilp->mrhs = mrhs; inp_ilp->msense = msense;
}
void Cgl012Cut::free_ilp()
{
  free(inp_ilp);
  inp_ilp=NULL;
}

#ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_CUTS
void Cgl012Cut::print_constr(int i /* constraint to be printed */)
{

  printf("\n content of constraint %d: nzcnt = %d, rhs = %d, sense = %c, slack = %f\n",
    i, inp_ilp->mtcnt[i], inp_ilp->mrhs[i], inp_ilp->msense[i], p_ilp->slack[i]);
  print_int_vect(const_cast<char*>("ind"),inp_ilp->mtind + inp_ilp->mtbeg[i],inp_ilp->mtcnt[i]);
  print_int_vect(const_cast<char*>("val"),inp_ilp->mtval + inp_ilp->mtbeg[i],inp_ilp->mtcnt[i]);
}
#endif

/* alloc_parity_ilp: allocate the memory for the parity ILP data structure */

void Cgl012Cut::alloc_parity_ilp(
		      int mr, /* number of rows in the ILP matrix */
		      int mc, /* number of columns in the ILP matrix */
		      int mnz /* number of nonzero's in the ILP matrix */
		      )
{
  p_ilp = reinterpret_cast<parity_ilp *> (calloc(1,sizeof(parity_ilp)));
  if ( p_ilp == NULL ) alloc_error(const_cast<char*>("p_ilp"));

  p_ilp->mtbeg = reinterpret_cast<int *> (calloc(mr,sizeof(int))); 
  if ( p_ilp->mtbeg == NULL ) alloc_error(const_cast<char*>("p_ilp->mtbeg"));
  p_ilp->mtcnt = reinterpret_cast<int *> (calloc(mr,sizeof(int)));
  if ( p_ilp->mtcnt == NULL ) alloc_error(const_cast<char*>("p_ilp->mtcnt"));
  p_ilp->mtind = reinterpret_cast<int *> (calloc(mnz,sizeof(int)));
  if ( p_ilp->mtind == NULL ) alloc_error(const_cast<char*>("p_ilp->mtind"));
  p_ilp->mrhs = reinterpret_cast<short int *> (calloc(mr,sizeof(short int)));
  if ( p_ilp->mrhs== NULL ) alloc_error(const_cast<char*>("p_ilp->mrhs"));
  p_ilp->slack = reinterpret_cast<double *> (calloc(mr,sizeof(double)));
  if ( p_ilp->slack == NULL ) alloc_error(const_cast<char*>("p_ilp->slack"));
  p_ilp->row_to_delete = reinterpret_cast<short int *> (calloc(mr,sizeof(short int)));
  if ( p_ilp->row_to_delete == NULL ) alloc_error(const_cast<char*>("p_ilp->row_to_delete"));
  p_ilp->col_to_delete = reinterpret_cast<short int *> (calloc(mc,sizeof(short int)));
  if ( p_ilp->col_to_delete == NULL ) alloc_error(const_cast<char*>("p_ilp->col_to_delete"));
  p_ilp->gcd = reinterpret_cast<int *> (calloc(mr,sizeof(int)));
  if ( p_ilp->gcd == NULL ) alloc_error(const_cast<char*>("p_ilp->gcd"));
  p_ilp->possible_weak = reinterpret_cast<short int *> (calloc(mc,sizeof(short int)));
  if ( p_ilp->possible_weak == NULL ) alloc_error(const_cast<char*>("p_ilp->possible_weak"));
  p_ilp->type_even_weak = reinterpret_cast<short int *> (calloc(mc,sizeof(short int)));
  if ( p_ilp->type_even_weak == NULL ) alloc_error(const_cast<char*>("p_ilp->type_even_weak"));
  p_ilp->type_odd_weak = reinterpret_cast<short int *> (calloc(mc,sizeof(short int)));
  if ( p_ilp->type_odd_weak == NULL ) alloc_error(const_cast<char*>("p_ilp->type_odd_weak"));
  p_ilp->loss_even_weak = reinterpret_cast<double *> (calloc(mc,sizeof(double)));
  if ( p_ilp->loss_even_weak == NULL ) alloc_error(const_cast<char*>("p_ilp->loss_even_weak"));
  p_ilp->loss_odd_weak = reinterpret_cast<double *> (calloc(mc,sizeof(double)));
  if ( p_ilp->loss_odd_weak == NULL ) alloc_error(const_cast<char*>("p_ilp->loss_odd_weak"));
  p_ilp->min_loss_by_weak = reinterpret_cast<double *> (calloc(mc,sizeof(double)));
  if ( p_ilp->min_loss_by_weak == NULL ) alloc_error(const_cast<char*>("p_ilp->min_loss_by_weak"));
  p_ilp->mr=mr;
  p_ilp->mc=mc;
  p_ilp->mnz=mnz;
}

#ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_CUTS
void Cgl012Cut::print_parity_ilp()
{
  printf("\n content of parity_ilp data structure: mc = %d, mr = %d, mnz = %d\n",
    p_ilp->mc,p_ilp->mr,p_ilp->mnz);
  print_int_vect(const_cast<char*>("mtbeg"),p_ilp->mtbeg,p_ilp->mr); 
  print_int_vect(const_cast<char*>("mtcnt"),p_ilp->mtcnt,p_ilp->mr); 
  print_int_vect(const_cast<char*>("mtind"),p_ilp->mtind,p_ilp->mnz); 
print_short_int_vect(const_cast<char*>("mrhs"),p_ilp->mrhs,p_ilp->mr); 
print_double_vect(const_cast<char*>("xstar"),p_ilp->xstar,p_ilp->mc); 
print_double_vect(const_cast<char*>("slack"),p_ilp->slack,p_ilp->mr); 
print_short_int_vect(const_cast<char*>("row_to_delete"),p_ilp->row_to_delete,p_ilp->mr); 
print_short_int_vect(const_cast<char*>("col_to_delete"),p_ilp->col_to_delete,p_ilp->mc); 
print_int_vect(const_cast<char*>("gcd"),p_ilp->gcd,p_ilp->mr); 
print_short_int_vect(const_cast<char*>("possible_weak"),p_ilp->possible_weak,p_ilp->mc);
print_short_int_vect(const_cast<char*>("type_even_weak"),p_ilp->type_even_weak,p_ilp->mc);
print_short_int_vect(const_cast<char*>("type_odd_weak"),p_ilp->type_odd_weak,p_ilp->mc);
print_double_vect(const_cast<char*>("loss_even_weak"),p_ilp->loss_even_weak,p_ilp->mc); 
print_double_vect(const_cast<char*>("loss_odd_weak"),p_ilp->loss_odd_weak,p_ilp->mc); 
print_double_vect(const_cast<char*>("min_loss_by_weak"),p_ilp->min_loss_by_weak,p_ilp->mc); 
}
#endif

void Cgl012Cut::free_parity_ilp()
{
  if (p_ilp) {
    free(p_ilp->mtbeg); 
    free(p_ilp->mtcnt); 
    free(p_ilp->mtind); 
    free(p_ilp->mrhs); 
    free(p_ilp->slack); 
    free(p_ilp->row_to_delete); 
    free(p_ilp->col_to_delete); 
    free(p_ilp->gcd); 
    free(p_ilp->possible_weak); 
    free(p_ilp->type_even_weak); 
    free(p_ilp->type_odd_weak); 
    free(p_ilp->loss_even_weak); 
    free(p_ilp->loss_odd_weak); 
    free(p_ilp->min_loss_by_weak); 
    free(p_ilp);
    p_ilp=NULL;
  }
}

/* alloc_info_weak: allocate memory for the weakening info data structure */

info_weak *alloc_info_weak(int nweak /* number of variables to be weakened */)
{
  info_weak *i_weak;

  i_weak = reinterpret_cast<info_weak *> (calloc(1,sizeof(info_weak)));
  if ( i_weak == NULL ) alloc_error(const_cast<char*>("i_weak"));
  if ( nweak > 0 ) {
    i_weak->var = reinterpret_cast<int *> (calloc(nweak,sizeof(int)));
    if ( i_weak->var == NULL ) alloc_error(const_cast<char*>("i_weak->var"));
    i_weak->type = reinterpret_cast<short int *> (calloc(nweak,sizeof(short int)));
    if ( i_weak->type == NULL ) alloc_error(const_cast<char*>("i_weak->type"));
  }

  return(i_weak);
}

#ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_CUTS
void print_info_weak(info_weak *i_weak)
{
  printf("\n content of info_weak: nweak = %d\n",i_weak->nweak); 
  if ( i_weak->nweak > 0 ) {
    print_int_vect(const_cast<char*>("var"),i_weak->var,i_weak->nweak);
    print_short_int_vect(const_cast<char*>("type"),i_weak->type,i_weak->nweak);
  }
}
#endif

void free_info_weak(info_weak *i_weak)
{
  if ( i_weak->nweak > 0 ) {
    free(i_weak->var);
    free(i_weak->type);
  }
  free(i_weak);
}

/* get_parity_ilp: construct an internal data structure containing all the 
   information which can be useful for  0-1/2 cut separation */

void Cgl012Cut::get_parity_ilp()
{
  int i, j, h, ij, aij, cnti, cnttot, begi, begh, ofsj, gcdi, ubj, lbj;
  double slacki, xstarj, loss_upper, loss_lower;
  short int equalih;

  /* allocate the memory for the parity ILP data structure */

  //alloc_parity_ilp(inp_ilp->mr,inp_ilp->mc,inp_ilp->mnz);

  p_ilp->mr = inp_ilp->mr;
  p_ilp->mc = inp_ilp->mc;
  p_ilp->xstar = inp_ilp->xstar;
  
  /* mark the variables equal to their lower/upper bound */

  for ( j = 0; j < inp_ilp->mc; j++ ) {
    if ((j & 4095) == 0 && checkTimeLimit("get_parity_ilp timeout", "during column scan")) {
      p_ilp->mnz = 0;
      return;
    }
    xstarj = inp_ilp->xstar[j];
    ubj = inp_ilp->vub[j]; lbj = inp_ilp->vlb[j];
    if ( xstarj > static_cast<double> (ubj - ZERO) ) {
      /* variable at its upper bound */
      p_ilp->col_to_delete[j] = TRUE;
      p_ilp->min_loss_by_weak[j] = 0.0;
      if ( mod2(ubj) == ODD ) {
	p_ilp->possible_weak[j] = ODD;
	p_ilp->type_odd_weak[j] = UPPER_BOUND;
	p_ilp->loss_odd_weak[j] = 0.0;
      }
      else {
	p_ilp->possible_weak[j] = EVEN;
	p_ilp->type_even_weak[j] = UPPER_BOUND;
	p_ilp->loss_even_weak[j] = 0.0;
      }
    }
    else if ( xstarj < static_cast<double> (lbj) + ZERO ) {
      /* variable at its lower bound */
      p_ilp->col_to_delete[j] = TRUE;
      p_ilp->min_loss_by_weak[j] = 0.0;
      if ( mod2(lbj) == ODD ) {
	p_ilp->possible_weak[j] = ODD;
	p_ilp->type_odd_weak[j] = LOWER_BOUND;
	p_ilp->loss_odd_weak[j] = 0.0;
      }
      else {
	p_ilp->possible_weak[j] = EVEN;
	p_ilp->type_even_weak[j] = LOWER_BOUND;
	p_ilp->loss_even_weak[j] = 0.0;
      }
    }
    else { 
      /* variable neither at its lower nor at its upper bound */
      p_ilp->col_to_delete[j] = FALSE;
      loss_upper = static_cast<double> (ubj) - xstarj; 
      loss_lower = xstarj - static_cast<double> (lbj); 
      if ( ( loss_upper > MAX_LOSS ) && ( loss_lower > MAX_LOSS ) ) 
	/* no weakening for the variable */
	p_ilp->possible_weak[j] = NONE;
      else if ( loss_upper > MAX_LOSS ) {
	/* lower weakening only */
	if ( mod2(lbj) == EVEN ) {
	  p_ilp->possible_weak[j] = EVEN;
	  p_ilp->type_even_weak[j] = LOWER_BOUND;
	  p_ilp->loss_even_weak[j] = loss_lower;
	}
	else {
	  p_ilp->possible_weak[j] = ODD;
	  p_ilp->type_odd_weak[j] = LOWER_BOUND;
	  p_ilp->loss_odd_weak[j] = loss_lower;
	}
      }
      else if ( loss_lower > MAX_LOSS ) {
	/* upper weakening only */
	if ( mod2(ubj) == EVEN ) {
	  p_ilp->possible_weak[j] = EVEN;
	  p_ilp->type_even_weak[j] = UPPER_BOUND;
	  p_ilp->loss_even_weak[j] = loss_upper;
	}
	else {
	  p_ilp->possible_weak[j] = ODD;
	  p_ilp->type_odd_weak[j] = UPPER_BOUND;
	  p_ilp->loss_odd_weak[j] = loss_upper;
	}
      }
      else if ( mod2(ubj) == mod2(lbj) ) {
	/* lower and upper bound have the same parity: 
	   choose the best weakening */
	if ( mod2(ubj) == EVEN ) {
	  p_ilp->possible_weak[j] = EVEN;
	  if ( loss_lower <= loss_upper ) {
	    p_ilp->type_even_weak[j] = LOWER_BOUND;
	    p_ilp->loss_even_weak[j] = loss_lower;
	  }
	  else {  
	    p_ilp->type_even_weak[j] = UPPER_BOUND;
	    p_ilp->loss_even_weak[j] = loss_upper;
	  }
	}
	else {
	  p_ilp->possible_weak[j] = ODD;
	  if ( loss_lower <= loss_upper ) {
	    p_ilp->type_odd_weak[j] = LOWER_BOUND;
	    p_ilp->loss_odd_weak[j] = loss_lower;
	  }
	  else {  
	    p_ilp->type_odd_weak[j] = UPPER_BOUND;
	    p_ilp->loss_odd_weak[j] = loss_upper;
	  }
	}
      }
      else {
	/* lower and upper bound have different parities: 
	   consider both weakenings */
	p_ilp->possible_weak[j] = BOTH;
	if ( mod2(ubj) == EVEN ) {
	  p_ilp->type_even_weak[j] = UPPER_BOUND;
	  p_ilp->loss_even_weak[j] = loss_upper;
	  p_ilp->type_odd_weak[j] = LOWER_BOUND;
	  p_ilp->loss_odd_weak[j] = loss_lower;
	}
	else {  
	  p_ilp->type_even_weak[j] = LOWER_BOUND;
	  p_ilp->loss_even_weak[j] = loss_lower;
	  p_ilp->type_odd_weak[j] = UPPER_BOUND;
	  p_ilp->loss_odd_weak[j] = loss_upper;
	}
      }
      if ( loss_upper > loss_lower ) p_ilp->min_loss_by_weak[j] = loss_lower;
      else p_ilp->min_loss_by_weak[j] = loss_upper;
    }
  }    
  
  /* scan the constraints and delete those which are trivially useless 
     in the 0-1/2 cut separation */

  cnttot = 0;
  for ( i = 0; i < inp_ilp->mr; i++ ) {
    if ((i & 1023) == 0 &&
        checkTimeLimit("get_parity_ilp timeout", "during row scan")) {
      p_ilp->mnz = 0;
      return;
    }
    begi = inp_ilp->mtbeg[i];
    
    /* compute the row slack and the GCD of the entries of the row */
    
    slacki = static_cast<double> (inp_ilp->mrhs[i]); 
    gcdi = inp_ilp->mrhs[i];
    for ( ofsj = 0; ofsj < inp_ilp->mtcnt[i]; ofsj++ ) {
      ij = begi + ofsj;
      j = inp_ilp->mtind[ij];
      aij = inp_ilp->mtval[ij];
      slacki -= static_cast<double> (aij ) * ( inp_ilp->xstar[j] );        
      gcdi = gcd(gcdi,aij);
    }
    if ( inp_ilp->msense[i] == 'G' ) slacki = -slacki;
    if ( slacki < -ZERO || ( inp_ilp->msense[i] == 'E' && slacki > ZERO ) ) {
#ifdef COIN_DEVELOP
      printf("\n Warning: constraint %d in the model is violated:\n",i);
      printf("\n 0-1/2 cut separation is not possible\n");
printf("\nnumber of nonzero's %d\n",inp_ilp->mtcnt[i]);
printf("nonzero's (col,coef,xstar) ");
for (ofsj=0;ofsj<inp_ilp->mtcnt[i];ofsj++) printf("(%d,%d,%f) ",
  inp_ilp->mtind[begi+ofsj], inp_ilp->mtval[begi+ofsj],
  inp_ilp->xstar[inp_ilp->mtind[begi+ofsj]]);
printf("\n");
printf("sense %c and rhs %d and slack %.5e\n",inp_ilp->msense[i],inp_ilp->mrhs[i], slacki);
#endif
      //exit(0);
	  slacki = INF;
    }
    p_ilp->slack[i] = slacki;
  
    /* mark the rows with slack greater than the maximum allowed */
    
    if ( slacki > MAX_SLACK - EPS ) p_ilp->row_to_delete[i] = TRUE;
    else p_ilp->row_to_delete[i] = FALSE;

    /* store the odd entries in the (possibly scaled) row i */
   
    //if ( gcdi != 1 ) 
    //printf("Warning: constraint %d with nonprime coefficients\n",i); 
    p_ilp->gcd[i] = gcdi;
    p_ilp->mrhs[i] = mod2(( inp_ilp->mrhs[i] / gcdi ));
    p_ilp->mtbeg[i] = cnttot;
    cnti = 0;
    for ( ofsj = 0; ofsj < inp_ilp->mtcnt[i]; ofsj++ ) {
      ij = begi + ofsj;
      j = inp_ilp->mtind[ij];
      aij = mod2(( inp_ilp->mtval[ij] / gcdi ));
      if ( aij == ODD ) { 
	if ( ! p_ilp->col_to_delete[j] ) {
	  p_ilp->mtind[cnttot] = j;
	  cnti++; cnttot++;       
	}
	else if ( p_ilp->possible_weak[j] == ODD ) {
	  if ( p_ilp->mrhs[i] == EVEN ) p_ilp->mrhs[i] = ODD; 
	  else p_ilp->mrhs[i] = EVEN;
	}
      }
    }
    p_ilp->mtcnt[i] = cnti;
    if ( cnti == 0 ) /* (scaled) row with even entries only */
      p_ilp->row_to_delete[i] = TRUE;
    else {
      if ( cnti == 1 && slacki < EPS ) { /* the row could be deleted */
#ifdef PRINT
	printf("get_parity_ilp: row %d could be deleted since it\n",i);
	printf("has only one odd entry and is tight, but it is not ...\n");
#endif
      }
    }
  }    
  p_ilp->mnz = cnttot;

#ifdef REDUCTION
  
  /* remove identical rows in the parity matrix */
  /* very trivial implementation */

  for ( i = 0; i < p_ilp->mr; i++ ) 
    for ( h = i+1; h < p_ilp->mr; h++ ) 
      if ( ( p_ilp->mrhs[i] == p_ilp->mrhs[h] ) &&
	   ( p_ilp->mtcnt[i] == p_ilp->mtcnt[h] ) && 
	   ( ! p_ilp->row_to_delete[i] ) && 
	   ( ! p_ilp->row_to_delete[h] ) ) {
	begi = p_ilp->mtbeg[i]; begh = p_ilp->mtbeg[h];
	equalih = TRUE;
	for ( ofsj = 0; ofsj < p_ilp->mtcnt[i]; ofsj++ )
	  /* the check assumes the indexes of the columns associated
	     with each row are ordered in p_ilp->mtind[] ... */
	  if ( p_ilp->mtind[begi+ofsj] != p_ilp->mtind[begh+ofsj] ) {
	    equalih = FALSE;
	    break;
	  }
	if ( equalih ) {
	  if ( p_ilp->slack[h] > p_ilp->slack[i] )
	    p_ilp->row_to_delete[h] = TRUE;
	  else 
	    p_ilp->row_to_delete[i] = TRUE;
	}
      }

  /* check for the existence of separate connected components in the 
     parity matrix row intersection graph */
  /* not implemented so far - if ever, the availability of the parity
     matrix in column form also would be really convenient */

#endif

}

/* separation graph subroutines */

#define SG_EDGE_INDEX(s_graph,J,K) ( ((J) < (K)) ? ( (s_graph->nnodes * (J)) - (((J)+1)*(J)/2) + (K) - (J) -1 ) : ( (s_graph->nnodes * (K)) - (((K)+1)*(K)/2) + (J) - (K) -1 ) ) 
  
/* initialize_sep_graph: allocate and initialize the data structure
   to contain the information associated with a separation graph */

separation_graph *Cgl012Cut::initialize_sep_graph()
{
  int maxnodes, nnodes, j;
  long long jk;
  long long maxedges;
  int *nodes, *ind;
  separation_graph *s_graph;
  double estimatedGiB;

  s_graph = reinterpret_cast<separation_graph *> (calloc(1,sizeof(separation_graph)));
  if ( s_graph == NULL ) alloc_error(const_cast<char*>("s_graph"));
  
  maxnodes = p_ilp->mc + 1;
  nnodes = 0;
  nodes = reinterpret_cast<int *> (calloc(maxnodes,sizeof(int)));
  if ( nodes == NULL ) alloc_error(const_cast<char*>("nodes"));
  ind = reinterpret_cast<int *> (calloc(maxnodes,sizeof(int)));
  if ( ind == NULL ) alloc_error(const_cast<char*>("ind"));
  for ( j = 0; j < p_ilp->mc; j++ )
    if ( ! p_ilp->col_to_delete[j] ) {
      /* variable not removed from the separation problem */
      nodes[nnodes] = j;
      ind[j] = nnodes;
      nnodes++;
    }

  /* take into account the special node */
  nodes[nnodes] = maxnodes - 1;
  ind[maxnodes-1] = nnodes;
  nnodes++; 
  
  s_graph->nnodes = nnodes;
  s_graph->nedges = 0;
  s_graph->nodes = reinterpret_cast<int *> (malloc(nnodes*sizeof(int)));
  if ( s_graph->nodes == NULL ) alloc_error(const_cast<char*>("s_graph->nodes"));
  for ( j = 0; j < nnodes; j++ ) s_graph->nodes[j] = nodes[j];
  free(nodes);
  s_graph->ind = reinterpret_cast<int *> (malloc(maxnodes*sizeof(int)));
  if ( s_graph->ind == NULL ) alloc_error(const_cast<char*>("s_graph->ind"));
  for ( j = 0; j < maxnodes; j++ ) s_graph->ind[j] = ind[j];
  free(ind);
  maxedges = ((long long)nnodes * ((long long)nnodes - 1)) / 2;
  if ( cglZeroHalfUseSparseSepGraph(nnodes,maxedges,sepGraphSparseThreshold_,&estimatedGiB) ) {
    s_graph->sparseMode = true;
    s_graph->sparseEdges = new std::unordered_map<std::uint64_t, edge *>();
    s_graph->sparseAdj = new std::vector<edge *>[nnodes];
    if ( cglZeroHalfSepGraphTooLarge(nnodes,maxedges,&estimatedGiB) ) {
      printf("Warning: using sparse 0-1/2 cut separation due to dense graph with %d active nodes\n", nnodes);
      printf("         dense graph would need %lld edge slots (about %.2f GiB for even/odd adjacency tables)\n",
        maxedges, estimatedGiB);
    }
    return(s_graph);
  }
  s_graph->even_adj_list = reinterpret_cast<edge **> (malloc(maxedges*sizeof(edge *)));
  if ( s_graph->even_adj_list == NULL ) alloc_error(const_cast<char*>("s_graph->even_adj_list"));
  s_graph->odd_adj_list = reinterpret_cast<edge **> (malloc(maxedges*sizeof(edge *)));
  if ( s_graph->odd_adj_list == NULL ) alloc_error(const_cast<char*>("s_graph->odd_adj_list"));
  for ( jk = 0; jk < maxedges; jk++ ) 
    s_graph->even_adj_list[jk] = s_graph->odd_adj_list[jk] = NULL;

  return(s_graph);
}

/* update_weight_sep_graph: consider a new edge obtained from the 
   (weakened) parity ILP and (possibly) add it to the separation graph */

separation_graph *update_weight_sep_graph(
					  int j, int k, /* endpoints of the new edge */
					  double weight, /* weight of the new edge */
					  short int parity, /* parity of the new edge */
					  int i, /* constraint associated with the new edge */
					  info_weak *i_weak, /* information associated with the weakening */
					  separation_graph *s_graph /* separation graph to be updated */
					  )
{
  int indj, indk, indjk;
  edge *old_edge, *new_edge;
  
  indj = s_graph->ind[j]; indk = s_graph->ind[k]; 
  if ( s_graph->sparseMode ) {
    const std::uint64_t key = cglZeroHalfSparseEdgeKey(indj,indk,parity);
    std::unordered_map<std::uint64_t, edge *>::iterator it = s_graph->sparseEdges->find(key);
    old_edge = (it == s_graph->sparseEdges->end()) ? NULL : it->second;
    if ( old_edge == NULL ) {
      new_edge = reinterpret_cast<edge *> (calloc(1,sizeof(edge)));
      if ( new_edge == NULL ) alloc_error(const_cast<char*>("new_edge"));
      new_edge->endpoint1 = indj; new_edge->endpoint2 = indk;
      new_edge->weight = weight; new_edge->parity = parity;
      new_edge->constr = i; new_edge->weak = i_weak;
      (*s_graph->sparseEdges)[key] = new_edge;
      s_graph->sparseAdj[indj].push_back(new_edge);
      s_graph->sparseAdj[indk].push_back(new_edge);
      (s_graph->nedges)++;
    }
    else {
      if ( old_edge->weight > weight ) {
        old_edge->weight = weight; old_edge->constr = i;
        free_info_weak(old_edge->weak);
        old_edge->weak = i_weak;
      }
      else {
        free_info_weak(i_weak);
      }
    }
    return(s_graph);
  }
  indjk = SG_EDGE_INDEX(s_graph,indj,indk);
  if ( parity == EVEN ) old_edge = s_graph->even_adj_list[indjk];
  else old_edge = s_graph->odd_adj_list[indjk];
  if ( old_edge == NULL ) {
    /* edge is not in the graph */
    new_edge = reinterpret_cast<edge *> (calloc(1,sizeof(edge)));
    if ( new_edge == NULL ) alloc_error(const_cast<char*>("new_edge"));
    new_edge->endpoint1 = indj; new_edge->endpoint2 = indk;
    new_edge->weight = weight; new_edge->parity = parity;
    new_edge->constr = i; new_edge->weak = i_weak;
    (s_graph->nedges)++;
    if ( parity == EVEN ) s_graph->even_adj_list[indjk] = new_edge;
    else s_graph->odd_adj_list[indjk] = new_edge;
  }
  else {
    /* edge is already in the graph */
    if ( old_edge->weight > weight ) {
      /* replace the old edge */
      old_edge->weight = weight; old_edge->constr = i; 
      free_info_weak(old_edge->weak); 
      old_edge->weak = i_weak;
    }   
    else {
      /* keep the old edge */
      free_info_weak(i_weak);
    }
  }

  return(s_graph);
}

#ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_CUTS
void print_edge(edge *e)
{
  printf("\n content of edge: endpoint1 = %d, endpoint2 = %d, weight = %f, parity = %d, constr = %d\n",
    e->endpoint1,e->endpoint2,e->weight,e->parity,e->constr);
  print_info_weak(e->weak);
}
#endif

void free_edge(edge *e)
{
  if ( e->weak != NULL )  
    free_info_weak(e->weak);
  free(e);
}

#ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_CUTS
void print_sep_graph(separation_graph *s_graph)
{
  int nnodes;
  long long maxedges, jk;
  
  nnodes = s_graph->nnodes;
  printf("\n content of separation_graph: nnodes = %d, nedges = %d\n",
    nnodes, s_graph->nedges);
  print_int_vect(const_cast<char*>("nodes"),s_graph->nodes,nnodes);
  print_int_vect(const_cast<char*>("ind"),s_graph->ind,nnodes);
  if ( s_graph->sparseMode ) {
    for (std::unordered_map<std::uint64_t, edge *>::const_iterator it = s_graph->sparseEdges->begin();
         it != s_graph->sparseEdges->end(); ++it)
      print_edge(it->second);
    return;
  }
  maxedges = (nnodes * (nnodes - 1)) / 2;
  for ( jk = 0; jk < maxedges; jk++ ) {
    if ( s_graph->even_adj_list[jk] != NULL )
      print_edge(s_graph->even_adj_list[jk]);
    if ( s_graph->odd_adj_list[jk] != NULL )
      print_edge(s_graph->odd_adj_list[jk]);
  }
}
#endif

void free_sep_graph(separation_graph *s_graph)
{
  int nnodes;
  long long maxedges, jk;

  if ( s_graph->sparseMode ) {
    if ( s_graph->sparseEdges != NULL ) {
      for (std::unordered_map<std::uint64_t, edge *>::const_iterator it = s_graph->sparseEdges->begin();
           it != s_graph->sparseEdges->end(); ++it)
        free_edge(it->second);
      delete s_graph->sparseEdges;
    }
    delete [] s_graph->sparseAdj;
    free(s_graph->nodes);
    free(s_graph->ind);
    free(s_graph);
    return;
  }

  nnodes = s_graph->nnodes;
  maxedges = (nnodes * (nnodes - 1)) / 2;
  for ( jk = 0; jk < maxedges; jk++ ) {
    if ( s_graph->even_adj_list[jk] != NULL )
      free_edge(s_graph->even_adj_list[jk]);
    if ( s_graph->odd_adj_list[jk] != NULL )
      free_edge(s_graph->odd_adj_list[jk]);
  }
  free(s_graph->nodes);
  free(s_graph->ind);
  free(s_graph->even_adj_list);
  free(s_graph->odd_adj_list);
  free(s_graph);
}

/* auxiliary graph subroutines - depend on the shortest path code used */

#ifndef CGL_NEW_SHORT
// will error if we get here
#include "Cgldikbd.c"
#endif

#define AG_TWIN1(J) 2 * J
#define AG_TWIN2(J) 2 * J + 1
#define AG_MATE(J) 2 * static_cast<int> ( J / 2 ) + ( J % 2 == EVEN ? 1 : 0 )
#define AG_TYPE(J,K) ( (J % 2) == (K % 2) ? EVEN : ODD )
#define SG_ORIG(J) static_cast<int> (J / 2)

/* define_aux_graph: construct the auxiliary graph for the shortest 
   path computation - the data structure is based on that used by
   Cherkassky, Goldberg and Radzik's shortest path codes */
     
auxiliary_graph *define_aux_graph(separation_graph *s_graph /* input separation graph */)
{
  int j, k, indjk, auxj1, auxj2, auxk1, auxk2, /*noutj,*/ totoutj, narcs;
  edge *s_edge;
  auxiliary_graph *a_graph;

  a_graph = reinterpret_cast<auxiliary_graph *> (calloc(1,sizeof(auxiliary_graph)));
  if ( a_graph == NULL ) alloc_error(const_cast<char*>("a_graph"));

  a_graph->nnodes = 2 * s_graph->nnodes;
  a_graph->narcs = 4 * s_graph->nedges;

#ifndef CGL_NEW_SHORT
  a_graph->nodes = reinterpret_cast<node *> (calloc((a_graph->nnodes + 1),sizeof(node)));
#else
  a_graph->nodes = reinterpret_cast<cgl_node *> (calloc((a_graph->nnodes + 1),sizeof(cgl_node)));
#endif
  if ( a_graph->nodes == NULL ) alloc_error(const_cast<char*>("a_graph->nodes"));
#ifndef CGL_NEW_SHORT
  a_graph->arcs = reinterpret_cast<arc *> (calloc(((a_graph->narcs) + 1),sizeof(arc)));
#else
  a_graph->arcs = reinterpret_cast<cgl_arc *> (calloc(((a_graph->narcs) + 1),sizeof(cgl_arc)));
#endif
  if ( a_graph->arcs == NULL ) alloc_error(const_cast<char*>("a_graph->arcs"));

  narcs = 0; 
  for ( j = 0; j < s_graph->nnodes; j++ ) {
    /* count the number of edges incident with j in the separation graph */
    if ( s_graph->sparseMode ) {
      std::vector<edge *> &incident = s_graph->sparseAdj[j];
      std::sort(incident.begin(), incident.end(),
        [j](const edge *lhs, const edge *rhs) {
          const int lhsOther = lhs->endpoint1 == j ? lhs->endpoint2 : lhs->endpoint1;
          const int rhsOther = rhs->endpoint1 == j ? rhs->endpoint2 : rhs->endpoint1;
          if (lhsOther != rhsOther)
            return lhsOther < rhsOther;
          if (lhs->parity != rhs->parity)
            return lhs->parity < rhs->parity;
          if (lhs->constr != rhs->constr)
            return lhs->constr < rhs->constr;
          return lhs < rhs;
        });
      totoutj = static_cast<int>(incident.size());
    }
    else {
      totoutj = 0;
      for ( k = 0; k < s_graph->nnodes; k++ ) 
        if ( k != j ) {
	  indjk = SG_EDGE_INDEX(s_graph,j,k);
	  if ( s_graph->even_adj_list[indjk] != NULL ) totoutj++;
	  if ( s_graph->odd_adj_list[indjk] != NULL ) totoutj++;
        }
    }
    auxj1 = AG_TWIN1(j); auxj2 = AG_TWIN2(j);
    a_graph->nodes[auxj1].index = auxj1;
    a_graph->nodes[auxj2].index = auxj2;
#ifndef CGL_NEW_SHORT
    a_graph->nodes[auxj1].first = &(a_graph->arcs[narcs]);
    a_graph->nodes[auxj2].first = &(a_graph->arcs[narcs+totoutj]);
#else
    a_graph->nodes[auxj1].firstArc = &(a_graph->arcs[narcs]);
    a_graph->nodes[auxj2].firstArc = &(a_graph->arcs[narcs+totoutj]);
#endif
    /* add the edges as arcs outgoing from j to the auxiliary graph */
    //noutj = 0;
    if ( s_graph->sparseMode ) {
      const std::vector<edge *> &incident = s_graph->sparseAdj[j];
      for (std::vector<edge *>::const_iterator edgeIt = incident.begin(); edgeIt != incident.end(); ++edgeIt) {
        s_edge = *edgeIt;
        k = s_edge->endpoint1 == j ? s_edge->endpoint2 : s_edge->endpoint1;
        auxk1 = AG_TWIN1(k); auxk2 = AG_TWIN2(k);
        if ( s_edge->parity == EVEN ) {
#ifndef CGL_NEW_SHORT
	  a_graph->arcs[narcs].len = a_graph->arcs[narcs+totoutj].len = 
	    (int) (s_edge->weight * ISCALE); 
	  a_graph->arcs[narcs].head = &(a_graph->nodes[auxk1]);
	  a_graph->arcs[narcs+totoutj].head = &(a_graph->nodes[auxk2]);
          a_graph->arcs[narcs].backEdge = s_edge;
          a_graph->arcs[narcs+totoutj].backEdge = s_edge;
#else
	  a_graph->arcs[narcs].length = a_graph->arcs[narcs+totoutj].length = 
	    static_cast<int> (s_edge->weight * ISCALE); 
	  a_graph->arcs[narcs].to = auxk1;
	  a_graph->arcs[narcs+totoutj].to = auxk2;
          a_graph->arcs[narcs].backEdge = s_edge;
          a_graph->arcs[narcs+totoutj].backEdge = s_edge;
#endif
	  narcs++; //noutj++;
	}
        else {
#ifndef CGL_NEW_SHORT
	  a_graph->arcs[narcs].len = a_graph->arcs[narcs+totoutj].len = 
	    (int) (s_edge->weight * ISCALE); 
	  a_graph->arcs[narcs].head = &(a_graph->nodes[auxk2]);
	  a_graph->arcs[narcs+totoutj].head = &(a_graph->nodes[auxk1]);
          a_graph->arcs[narcs].backEdge = s_edge;
          a_graph->arcs[narcs+totoutj].backEdge = s_edge;
#else
	  a_graph->arcs[narcs].length = a_graph->arcs[narcs+totoutj].length = 
	    static_cast<int> (s_edge->weight * ISCALE); 
	  a_graph->arcs[narcs].to = auxk2;
	  a_graph->arcs[narcs+totoutj].to = auxk1;
          a_graph->arcs[narcs].backEdge = s_edge;
          a_graph->arcs[narcs+totoutj].backEdge = s_edge;
#endif
	  narcs++; //noutj++;
	}
      }
    }
    else {
      for ( k = 0; k < s_graph->nnodes; k++ ) {
        if ( k != j ) {
	  auxk1 = AG_TWIN1(k); auxk2 = AG_TWIN2(k);
	  indjk = SG_EDGE_INDEX(s_graph,j,k);
	  s_edge = s_graph->even_adj_list[indjk];
	  if ( s_edge != NULL ) {
	    /* there is an even edge between j and k */        
#ifndef CGL_NEW_SHORT
	    a_graph->arcs[narcs].len = a_graph->arcs[narcs+totoutj].len = 
	      (int) (s_edge->weight * ISCALE); 
	    a_graph->arcs[narcs].head = &(a_graph->nodes[auxk1]);
	    a_graph->arcs[narcs+totoutj].head = &(a_graph->nodes[auxk2]);
            a_graph->arcs[narcs].backEdge = s_edge;
            a_graph->arcs[narcs+totoutj].backEdge = s_edge;
#else
	    a_graph->arcs[narcs].length = a_graph->arcs[narcs+totoutj].length = 
	      static_cast<int> (s_edge->weight * ISCALE); 
	    a_graph->arcs[narcs].to = auxk1;
	    a_graph->arcs[narcs+totoutj].to = auxk2;
            a_graph->arcs[narcs].backEdge = s_edge;
            a_graph->arcs[narcs+totoutj].backEdge = s_edge;
#endif
	    narcs++; //noutj++;
	  }
	  s_edge = s_graph->odd_adj_list[indjk];
	  if ( s_edge != NULL ) {
	    /* there is an odd edge between j and k */        
#ifndef CGL_NEW_SHORT
	    a_graph->arcs[narcs].len = a_graph->arcs[narcs+totoutj].len = 
	      (int) (s_edge->weight * ISCALE); 
	    a_graph->arcs[narcs].head = &(a_graph->nodes[auxk2]);
	    a_graph->arcs[narcs+totoutj].head = &(a_graph->nodes[auxk1]);
            a_graph->arcs[narcs].backEdge = s_edge;
            a_graph->arcs[narcs+totoutj].backEdge = s_edge;
#else
	    a_graph->arcs[narcs].length = a_graph->arcs[narcs+totoutj].length = 
	      static_cast<int> (s_edge->weight * ISCALE); 
	    a_graph->arcs[narcs].to = auxk2;
	    a_graph->arcs[narcs+totoutj].to = auxk1;
            a_graph->arcs[narcs].backEdge = s_edge;
            a_graph->arcs[narcs+totoutj].backEdge = s_edge;
#endif
	    narcs++; //noutj++;
	  }
        }
      }
    }
    narcs += totoutj;
  }
#ifndef CGL_NEW_SHORT
  a_graph->nodes[a_graph->nnodes].first = &(a_graph->arcs[narcs]);
#else
  a_graph->nodes[a_graph->nnodes].firstArc = &(a_graph->arcs[narcs]);
#endif
  
  return(a_graph);
}

/* cancel_node_aux_graph: remove the node j in the separation graph
  from the auxiliary graph - all the outgoing arc lengths are set to
  a large value */

auxiliary_graph *cancel_node_aux_graph(
				       int j, /* index of the node in the separation graph */
				       auxiliary_graph *a_graph /* auxiliary graph to be updated */
				       )
{
  int auxj1, auxj2;
#ifndef CGL_NEW_SHORT
  arc *arc_ptr;
#else
  cgl_arc *arc_ptr;
#endif

  auxj1 = AG_TWIN1(j); auxj2 = AG_TWIN2(j); 
#ifndef CGL_NEW_SHORT
  for ( arc_ptr = a_graph->nodes[auxj1].first; 
	arc_ptr < a_graph->nodes[auxj1+1].first; 
	arc_ptr++ )  
    (*arc_ptr).len = ISCALE;
  for ( arc_ptr = a_graph->nodes[auxj2].first; 
	arc_ptr < a_graph->nodes[auxj2+1].first; 
	arc_ptr++ )  
    (*arc_ptr).len = ISCALE;
#else
  for ( arc_ptr = a_graph->nodes[auxj1].firstArc; 
	arc_ptr < a_graph->nodes[auxj1+1].firstArc; 
	arc_ptr++ )  
    (*arc_ptr).length = ISCALE;
  for ( arc_ptr = a_graph->nodes[auxj2].firstArc; 
	arc_ptr < a_graph->nodes[auxj2+1].firstArc; 
	arc_ptr++ )  
    (*arc_ptr).length = ISCALE;
#endif

  return(a_graph);
}
  
#ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_CUTS
#ifndef CGL_NEW_SHORT
void print_node(node *n)
{
  printf("\n content of node (addr = %d): first = %d, dist = %d, parent = %d, next = %d, prev = %d, status = %d\n",
    (int)n,(int)(*n).first,(*n).dist,(int)(*n).parent,(int)(*n).next,(int)(*n).prev,
    (*n).status);
}

void print_arc(arc *a)
{
  printf("\n content of arc (addr = %d): len = %d, head = %d\n",
    (int)a,(*a).len,(int)(*a).head);
}

void print_node_vect(char *s,node *v,int n)
{
  int i;

  printf("node vector %s:",s);
  for ( i = 0; i < n; i++ ) print_node(&v[i]);
  printf("\n");
}

void print_arc_vect(char *s,arc *v,int n)
{
  int i;

  printf("arc vector %s:",s);
  for ( i = 0; i < n; i++ ) print_arc(&v[i]);
  printf("\n");
}
#else
void print_node(cgl_node *n)
{
  printf("\n content of node (addr = %p): first = %p, dist = %d, parent = %p\n",
	 reinterpret_cast<void *>(n),
	 reinterpret_cast<void *>((*n).firstArc),
	 (*n).distanceBack,
	 static_cast<int>((*n).parentNode));
}

void print_arc(cgl_arc *a)
{
  printf("\n content of arc (addr = %p): len = %d, head = %d\n",
	 reinterpret_cast<void *>(a),(*a).length,static_cast<int>((*a).to));
}

void print_node_vect(char *s,cgl_node *v,int n)
{
  int i;

  printf("node vector %s:",s);
  for ( i = 0; i < n; i++ ) print_node(&v[i]);
  printf("\n");
}

void print_arc_vect(char *s,cgl_arc *v,int n)
{
  int i;

  printf("arc vector %s:",s);
  for ( i = 0; i < n; i++ ) print_arc(&v[i]);
  printf("\n");
}
#endif

void print_aux_graph(auxiliary_graph *a_graph)
{
  printf("\n content of auxiliary graph: nnodes = %d, narcs = %d\n",
    a_graph->nnodes,a_graph->narcs);
  print_node_vect(const_cast<char*>("nodes"),a_graph->nodes,a_graph->nnodes);
  print_arc_vect(const_cast<char*>("nodes"),a_graph->arcs,a_graph->narcs);
}
#endif

void free_aux_graph(auxiliary_graph *a_graph)
{
  free(a_graph->nodes);
  free(a_graph->arcs);
  free(a_graph);
}

/* odd cycles management subroutines */

/* simple_cycle: check whether a given cycle is simple
   (and therefore may correspond to a non-dominated ineq.) */

short int simple_cycle(cycle *s_cyc /* cycle to be checked */)
{
  int i, e, maxnodes;
  int *cnt;
 
  maxnodes = 0; 
  for ( e = 0; e < s_cyc->length; e++ ) {
    if (!s_cyc->edge_list[e]) {
      // bad
      maxnodes=-1;
      abort();//break;
    }
    i = s_cyc->edge_list[e]->endpoint1;
    if ( i > maxnodes ) maxnodes = i;
    i = s_cyc->edge_list[e]->endpoint2;
    if ( i > maxnodes ) maxnodes = i;
  }
  if (maxnodes<0)
    return FALSE;
  cnt = reinterpret_cast<int *> (calloc(maxnodes+1,sizeof(int)));
  if ( cnt == NULL ) alloc_error(const_cast<char*>("cnt"));
  //for ( i = 0; i <= maxnodes; i++ ) cnt[i] = 0;

  for ( e = 0; e < s_cyc->length; e++ ) {
    i = s_cyc->edge_list[e]->endpoint1;
    cnt[i]++;
    if ( cnt[i] > 2 ) {
      free(cnt);
      return(FALSE);
    }
    i = s_cyc->edge_list[e]->endpoint2;
    cnt[i]++;
    if ( cnt[i] > 2 ) {
      free(cnt);
      return(FALSE);
    }
  }

  free(cnt);
  return(TRUE);
}
  
/* same_cycle: check whether two cycles are identical 
   (assumes the first nodes of the cycles coincide) */

short int same_cycle(cycle *s_cyc1, cycle *s_cyc2 /* cycles to be compared */)
{
  int e, eb;
  short int same;

  if ( s_cyc1->length != s_cyc2->length ) return(FALSE);
  /* check the cycles in the same direction ... */
  same = TRUE;
  for ( e = 0; e < s_cyc1->length; e++ ) {
    if ( s_cyc1->edge_list[e] != s_cyc2->edge_list[e] ) {
      same = FALSE;
      break;
    }
  }
  if ( same ) return(TRUE);
  /* ... and in reverse direction */
  same = TRUE;
  for ( e = 0, eb = s_cyc2->length - 1; e < s_cyc1->length; e++, eb-- ) {
    if ( s_cyc1->edge_list[e] != s_cyc2->edge_list[eb] ) {
      same = FALSE;
      break;
    }
  }
  if ( same ) return(TRUE);
  return(FALSE);
}

#ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_CUTS
void print_cycle(cycle *s_cycle)
{
  int e;

  printf("\n content of cycle: weight = %f, length = %d\n",
    s_cycle->weight,s_cycle->length);
  for ( e = 0; e < s_cycle->length; e++ )
    print_edge(s_cycle->edge_list[e]);
}
#endif

void free_cycle(cycle *s_cycle)
{
  free(s_cycle->edge_list);
  free(s_cycle);
}

/* initialize_cycle_list: allocate and initialize the cycle list data structure */

cycle_list *initialize_cycle_list(int max_cyc /* maximum number of cycles in the list */)
{
  cycle_list *s_cycle_list;

  s_cycle_list = reinterpret_cast<cycle_list *> (calloc(1,sizeof(cycle_list)));
  if ( s_cycle_list == NULL ) alloc_error(const_cast<char*>("s_cycle_list"));
  s_cycle_list->cnum = 0;
  s_cycle_list->list = reinterpret_cast<cycle **> (calloc(max_cyc,sizeof(cycle *)));
  if ( s_cycle_list->list == NULL ) alloc_error(const_cast<char*>("s_cycle_list->list"));
  return(s_cycle_list);
}

/* add_cycle_to_list: add a new cycle to the cycle list data structure
   (if not already in the list) */

cycle_list *add_cycle_to_list(
			      cycle *s_cycle, /* pointer to the cycle to be added to the list */
			      cycle_list *s_cycle_list /* input cycle list to be updated */
			      )
{
  int c;

  if ( ! simple_cycle(s_cycle) ) {
    free_cycle(s_cycle);
    return(s_cycle_list);
  }

  for ( c = 0; c < s_cycle_list->cnum; c++ )
    if ( same_cycle(s_cycle,s_cycle_list->list[c]) ) {
      free_cycle(s_cycle);
      return(s_cycle_list);
  }

  s_cycle_list->list[s_cycle_list->cnum] = s_cycle;
  (s_cycle_list->cnum)++;        

  return(s_cycle_list);
}

void free_cycle_list(cycle_list *s_cycle_list)
{
  int c;

  for ( c = 0; c < s_cycle_list->cnum; c++ ) 
    free_cycle(s_cycle_list->list[c]);
  free(s_cycle_list->list);
  free(s_cycle_list);
}

#ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_CUTS
void print_cycle_list(cycle_list *s_cycle_list)
{
  int c;

  printf("\n content of cycle_list: cnum = %d\n",s_cycle_list->cnum);
  for ( c = 0; c < s_cycle_list->cnum; c++ )
    print_cycle(s_cycle_list->list[c]);
}
#endif

/* get_shortest_odd_cycle_list: computation of the shortest odd cycles
   visiting a certain node in the separation graph, and each other 
   possible intermediate node, using the auxiliary graph data structure 
   for the shortest path computation - all the cycles in the list are
   different from each other */

cycle_list *Cgl012Cut::get_shortest_odd_cycle_list(
					int j, /* first node to be visited by the odd cycle */
					separation_graph *s_graph, /* current separation graph */
					auxiliary_graph *a_graph /* auxiliary graph for the shortest path computation */
					)
{
  int source, sink, curr, pred, totedges, k, t, kt;
  double weight;
#ifndef CGL_NEW_SHORT
  //node *source_ptr, *sink_ptr, *first_ptr; 
#else
  //cgl_node *source_ptr, *sink_ptr, *first_ptr; 
#endif
  edge *curr_edge;
  short_path_node *forw_arb, *backw_arb;
  cycle *s_cycle;
  cycle_list *s_cycle_list;

#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&tsi);
#endif
 
  s_cycle_list = initialize_cycle_list((a_graph->nnodes)-2);
  if (checkTimeLimit("shortest_odd_cycle timeout", "before shortest path"))
    return(s_cycle_list);

  source = AG_TWIN1(j); sink = AG_TWIN2(j);
  //source_ptr = &(a_graph->nodes[source]);
  //sink_ptr = &(a_graph->nodes[sink]);
  //first_ptr = &(a_graph->nodes[0]);

  /* compute the shortest path arborescence rooted at source and
     the shortest path arborescence rooted at sink (that comes for
     free due to symmetry) and store them (the path information is 
     hidden into aux_graph) */

#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&ti);
#endif
#ifndef CGL_NEW_SHORT
  {
    int nNodes = a_graph->nnodes;
    int nArcs = a_graph->narcs;
    cgl_arc * arcs = new cgl_arc [nArcs];
    for (int i=0;i<nArcs;i++) {
      arcs[i].length=a_graph->arcs[i].len;
      arcs[i].to=a_graph->arcs[i].head->index;
    }
    cgl_node * nodes = new cgl_node[nNodes+1];
    for (int i=0;i<nNodes;i++) {
      int iArc = a_graph->nodes[i].first-a_graph->arcs;
      nodes[i].firstArc=arcs+iArc;
      nodes[i].index=i;
    }
    int iArc = a_graph->nodes[nNodes].first-a_graph->arcs;
    nodes[nNodes].firstArc=arcs+iArc;
    cgl_graph graph;
    graph.nnodes=nNodes;
    graph.narcs=nArcs;
    graph.nodes=nodes;
    graph.arcs=arcs;
    cglShortestPath(&graph,source,ISCALE);
    dikbd(a_graph->nnodes,first_ptr,source_ptr,ISCALE);
    for ( k = 0; k < a_graph->nnodes; k++ ) { 
      if ( a_graph->nodes[k].parent != NULL ) {
	int distance1 = a_graph->nodes[k].dist;
	int distance2 = graph.nodes[k].distanceBack;
	assert (distance1==distance2);
      } else {
	//
	printf("null parent %d\n",k);
      }
    }
  }
#else
  cglShortestPath(a_graph,source,ISCALE);
#endif
  if (checkTimeLimit("shortest_odd_cycle timeout", "after shortest path"))
    return(s_cycle_list);
#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&tf);
  path_time += tf - ti;
#endif
  forw_arb = 
    reinterpret_cast<short_path_node *> (calloc(a_graph->nnodes,sizeof(short_path_node)));
  if ( forw_arb == NULL ) alloc_error(const_cast<char*>("forw_arb"));
  for ( k = 0; k < a_graph->nnodes; k++ ) { 
    if ((k & 255) == 0 &&
        checkTimeLimit("shortest_odd_cycle timeout", "while copying shortest-path trees")) {
      free(forw_arb);
      free_cycle_list(s_cycle_list);
      return(initialize_cycle_list(1));
    }
#ifndef CGL_NEW_SHORT
    if ( a_graph->nodes[k].parent != NULL ) {
      forw_arb[k].dist = a_graph->nodes[k].dist;
      forw_arb[k].pred = a_graph->nodes[k].parent->index;
    }
#else
    if ( a_graph->nodes[k].parentNode >=0 ) {
      forw_arb[k].dist = a_graph->nodes[k].distanceBack;
      forw_arb[k].pred = a_graph->nodes[k].parentNode;
    }
#endif
    else {
      forw_arb[k].dist = COIN_INT_MAX;
      forw_arb[k].pred = NONE;
    }
  }
  backw_arb = 
    reinterpret_cast<short_path_node *> (calloc(a_graph->nnodes,sizeof(short_path_node)));
  if ( backw_arb == NULL ) alloc_error(const_cast<char*>("backw_arb"));
  for ( k = 0; k < a_graph->nnodes; k++ ) { 
    if ((k & 255) == 0 &&
        checkTimeLimit("shortest_odd_cycle timeout", "while copying reverse shortest-path tree")) {
      free(forw_arb);
      free(backw_arb);
      free_cycle_list(s_cycle_list);
      return(initialize_cycle_list(1));
    }
#ifndef CGL_NEW_SHORT
    if ( a_graph->nodes[k].parent != NULL ) {
      backw_arb[AG_MATE(k)].dist = a_graph->nodes[k].dist;
      backw_arb[AG_MATE(k)].pred = AG_MATE(a_graph->nodes[k].parent->index);
    }
#else
    if ( a_graph->nodes[k].parentNode >=0) {
      backw_arb[AG_MATE(k)].dist = a_graph->nodes[k].distanceBack;
      backw_arb[AG_MATE(k)].pred = AG_MATE(a_graph->nodes[k].parentNode);
    }
#endif
    else {
      backw_arb[AG_MATE(k)].dist = COIN_INT_MAX;
      backw_arb[AG_MATE(k)].pred = NONE;
    }
  }

#ifdef USELESS
  /* compute second the shortest path anti-arborescence rooted at sink 
     (which coincides with the arborescence since aux_graph is 
     symmetrical) and store it */

#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&ti);
#endif
  cc = dikbd(a_graph->nnodes,first_ptr,sink_ptr,ISCALE);
#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&tf);
  path_time += tf - ti;
#endif
  backw_arb = 
    (short_path_node *) calloc(a_graph->nnodes,sizeof(short_path_node));
  if ( backw_arb == NULL ) alloc_error("backw_arb");
  for ( k = 0; k < a_graph->nnodes; k++ ) { 
    backw_arb[k].dist = a_graph->nodes[k].dist;
    backw_arb[k].pred = a_graph->nodes[k].parent->index;
  }
#endif

  /* consider each possible intermediate node in aux_graph */

  for ( k = 0; k < s_graph->nnodes; k++ ) {
    if ((k & 63) == 0 &&
        checkTimeLimit("shortest_odd_cycle timeout", "while enumerating intermediate nodes")) {
      free(forw_arb);
      free(backw_arb);
      return(s_cycle_list);
    }
    if ( k != j ) {
      for ( t = 1; t <= 2; t++ ) {
	if ( t == 1 ) kt = AG_TWIN1(k); 
	else kt = AG_TWIN2(k);
	weight = 
	  (static_cast<double> (forw_arb[kt].dist + backw_arb[kt].dist)) / 
	  (static_cast<double> (ISCALE));
	if ( weight < MAX_CYCLE_WEIGHT + EPS ) {
	  totedges = 0;
	  /* count how many edges are in the forward path from source ... */
	  curr = kt;
	  do {
	    if (curr<0) {
	      totedges=-1;
	      break;
	    }
	    curr = forw_arb[curr].pred; totedges++;
	  } while ( curr != source );
	  if (totedges>=0) {
	    /* ... and in the backward path to sink */
	    curr = kt;
	    do {
	      if (curr<0) {
		totedges=-1;
		break;
	      }
	      curr = backw_arb[curr].pred; totedges++;
	    } while ( curr != sink );
	  }
	  if (totedges>0) {
	    s_cycle = reinterpret_cast<cycle *> (calloc(1,sizeof(cycle)));
	    if ( s_cycle == NULL ) alloc_error(const_cast<char*>("s_cycle"));
	    s_cycle->weight = weight;
	    s_cycle->length = totedges;
	    s_cycle->edge_list = reinterpret_cast<edge **> (calloc(totedges,sizeof(edge *)));
	    if ( s_cycle->edge_list == NULL ) alloc_error(const_cast<char*>("s_cycle->edge_list"));
	    /* define the set of edges corresponding to the paths in sep_graph */
	    totedges = 0;
	    /* forward path from source ... */
	    curr = kt;
	    do {
	      pred = forw_arb[curr].pred;
	      curr_edge = cglZeroHalfGetAuxiliaryArcEdge(a_graph,pred,curr);
              if ( curr_edge == NULL ) abort();
	      s_cycle->edge_list[totedges] = curr_edge;
	      curr = pred; 
	      totedges++;
	    } while ( curr != source );
	    /* ... and backward path to sink */
	    curr = kt;
	    do {
	      pred = backw_arb[curr].pred;
	      curr_edge = cglZeroHalfGetAuxiliaryArcEdge(a_graph,pred,curr);
              if ( curr_edge == NULL ) abort();
	      s_cycle->edge_list[totedges] = curr_edge;
	      curr = pred; 
	      totedges++;
	    } while ( curr != sink );
	    /* insert the new cycle in the list */
	    s_cycle_list = add_cycle_to_list(s_cycle,s_cycle_list);
	  }
	}
      }
    }
  }
  free(forw_arb);
  free(backw_arb);

#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&tsf);
  cycle_time += tsf - tsi;
#endif

  return(s_cycle_list);
}

/* cut management subroutines */

/* initialize_cut_list: allocate and initialize the cut list data structure */

cut_list *initialize_cut_list(int max_cut /* maximum number of cuts in the list */)
{
  cut_list *cuts;

  cuts = reinterpret_cast<cut_list *> (calloc(1,sizeof(cut_list)));
  if ( cuts == NULL ) alloc_error(const_cast<char*>("cuts"));
  cuts->cnum = 0;
  cuts->list = reinterpret_cast<cut **> (calloc(max_cut,sizeof(cut *)));

  return(cuts);
}

#ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_CUTS
void Cgl012Cut::print_cut(cut *v_cut)
{
  printf("\n content of cut: n_of_constr = %d, cnzcnt = %d, crhs = %d, csense = %c, violation = %f\n",
    v_cut->n_of_constr,v_cut->cnzcnt,v_cut->crhs,v_cut->csense,v_cut->violation);
  print_int_vect(const_cast<char*>("cind"),v_cut->cind,v_cut->cnzcnt);
  print_int_vect(const_cast<char*>("cval"),v_cut->cval,v_cut->cnzcnt);
  if ( v_cut->constr_list != NULL ) 
    print_int_vect(const_cast<char*>("constr_list"),v_cut->constr_list,v_cut->n_of_constr);
  if ( v_cut->in_constr_list != NULL )
    print_short_int_vect(const_cast<char*>("in_constr_list"),v_cut->in_constr_list,inp_ilp->mr);
	;
}

void Cgl012Cut::print_cut_list(cut_list *cuts)
{
  int c;

  printf("\n content of cut_list: cnum = %d\n",cuts->cnum);
  for ( c = 0; c < cuts->cnum; c++ )
    print_cut(cuts->list[c]);
}
#endif

void free_cut(cut *v_cut)
{
  if ( v_cut->constr_list != NULL ) free(v_cut->constr_list);
  if ( v_cut->in_constr_list != NULL ) free(v_cut->in_constr_list);
  if ( v_cut->cind != NULL ) free(v_cut->cind);
  if ( v_cut->cval != NULL ) free(v_cut->cval);
  free(v_cut);
}

void free_cut_list(cut_list *cuts)
{
  int c;

  for ( c = 0; c < cuts->cnum; c++ ) 
    if ( cuts->list[c] != NULL ) free_cut(cuts->list[c]);
  free(cuts->list);
  free(cuts);
}

/* get_ori_cut_coef: get the coefficients of a cut, before dividing by 2 and
   rounding, starting from the list of the constraints combined to get 
   the cut */

short int Cgl012Cut::get_ori_cut_coef(
			   int n_of_constr, /* number of constraints combined */
			   int *constr_list, /* list of the constraints combined */
			   int *ccoef, /* cut left hand side coefficients */
			   int *crhs, /* cut right hand side */
			   short int only_viol /* flag which tells whether only an inequality of
			slack smaller than MAX_SLACK is of interest (TRUE)
			otherwise (FALSE) */
			   )
{
  int h, i, begi, gcdi, ofsj;
  double tot_slack;

  /* fast check of the possible violation of the cut */
  if ( only_viol ) {
    tot_slack = 0.0;
    for ( h = 0; h < n_of_constr; h++ ) {
      tot_slack += p_ilp->slack[constr_list[h]];
      if ( tot_slack > MAX_SLACK - EPS ) return(FALSE);
    }
  }
      
  //for ( j = 0; j < inp_ilp->mc; j++ )
  //ccoef[j] = 0;
  memset(ccoef,0,inp_ilp->mc*sizeof(int));
  (*crhs) = 0;

  for ( h = 0; h < n_of_constr; h++ ) {
    i = constr_list[h];
    begi = inp_ilp->mtbeg[i]; gcdi = p_ilp->gcd[i];
    if ( inp_ilp->msense[i] != 'G' ) {
      if ( gcdi == 1 ) {
	for ( ofsj = 0; ofsj < inp_ilp->mtcnt[i]; ofsj++ ) 
	  ccoef[inp_ilp->mtind[begi+ofsj]] += inp_ilp->mtval[begi+ofsj];
	(*crhs) += inp_ilp->mrhs[i];
      }
      else {
	for ( ofsj = 0; ofsj < inp_ilp->mtcnt[i]; ofsj++ ) 
	  ccoef[inp_ilp->mtind[begi+ofsj]] += inp_ilp->mtval[begi+ofsj] / gcdi;
	(*crhs) += inp_ilp->mrhs[i] / gcdi;
      }
    }
    else {
      if ( gcdi == 1 ) {
	for ( ofsj = 0; ofsj < inp_ilp->mtcnt[i]; ofsj++ ) 
	  ccoef[inp_ilp->mtind[begi+ofsj]] -= inp_ilp->mtval[begi+ofsj];
	(*crhs) -= inp_ilp->mrhs[i];
      }
      else {
	for ( ofsj = 0; ofsj < inp_ilp->mtcnt[i]; ofsj++ ) 
	  ccoef[inp_ilp->mtind[begi+ofsj]] -= inp_ilp->mtval[begi+ofsj] / gcdi;
	(*crhs) -= inp_ilp->mrhs[i] / gcdi;
      }
    }
  }
  
  return(TRUE);
}

/* best_cut: find the coefficients, the rhs and the violation of the
   best possible cut that can be obtained by weakening a given set of
   coefficients to even and a rhs to odd, dividing by 2 and rounding */

short int Cgl012Cut::best_cut(
		   int *ccoef, /* vector of the coefficients */
		   int *crhs, /* pointer to rhs value */
		   double *violation, /* violation of the cut */
		   short int update, /* TRUE/FALSE: if TRUE, the new ccoef and crhs are 
					given on output */ 
		   short int only_viol /* flag which tells whether only an inequality of
			slack smaller than MAX_SLACK is of interest (TRUE)
			otherwise (FALSE) */
		   )
{
  int j, n_to_weak;
  short int original_parity; 
  double original_slack, best_even_slack, best_odd_slack; 
  int *vars_to_weak;  
  info_weak *info_even_weak, *info_odd_weak; 

  /* choose the best weakening for the variables whose coefficient
     is not odd - this hopefully produces a stronger cut than that
     associated with the weakened inequalities to define the edges */

  ensureVarsToWeakBufferCapacity(inp_ilp->mc);
  vars_to_weak = varsToWeakBuffer_.data();
  
  n_to_weak = 0;
  original_slack = 0.0;
  for ( j = 0; j < inp_ilp->mc; j++ ) {
    if ((j & 4095) == 0 &&
        checkTimeLimit("best_cut timeout", "while scanning cut coefficients")) {
      return(FALSE);
    }
    if ( ccoef[j] != 0 ) {
      if ( mod2(ccoef[j]) == ODD ) {
	vars_to_weak[n_to_weak] = j;
	n_to_weak++;
      }
      original_slack -= inp_ilp->xstar[j] * static_cast<double> (ccoef[j]);
    }
  }
  original_slack += static_cast<double> (*crhs);
  if ( original_slack > MAX_SLACK - EPS ) {
    return(FALSE);
  }
  original_parity = mod2(*crhs);

  if ( best_weakening(n_to_weak,vars_to_weak,
		      original_parity,original_slack,
		      &best_even_slack,&best_odd_slack,
		      &info_even_weak,&info_odd_weak,
		      TRUE,only_viol) == ODD ) {
    *violation = ( 1.0 - best_odd_slack ) / 2.0;
    if ( ! update ) { 
      /* new ccoef and rhs are not required on output */
      free_info_weak(info_odd_weak);
      return(TRUE);
    }
    /* update ccoef and crhs according to the best odd weakening */
    for ( j = 0; j < n_to_weak; j++ )
      if ( info_odd_weak->type[j] == LOWER_BOUND ) {
        if ((j & 1023) == 0 &&
            checkTimeLimit("best_cut timeout", "while applying weakening")) {
          free_info_weak(info_odd_weak);
          return(FALSE);
        }
	ccoef[vars_to_weak[j]]--;
	*crhs -= inp_ilp->vlb[vars_to_weak[j]];
      }
      else {
        if ((j & 1023) == 0 &&
            checkTimeLimit("best_cut timeout", "while applying weakening")) {
          free_info_weak(info_odd_weak);
          return(FALSE);
        }
	ccoef[vars_to_weak[j]]++;
	*crhs += inp_ilp->vub[vars_to_weak[j]];
      }
    /* compute and check the correctness of the cut coefficients */
    for ( j = 0; j < inp_ilp->mc; j++ ) {
      if ((j & 4095) == 0 &&
          checkTimeLimit("best_cut timeout", "while normalizing coefficients")) {
        free_info_weak(info_odd_weak);
        return(FALSE);
      }
      if ( mod2(ccoef[j]) == ODD ) {
	printf("!!! Error 2 in weakening a cut !!!\n");
	exit(0);
      }
      if ( ccoef[j] != 0 ) ccoef[j] /= 2;
    }
    if ( mod2(*crhs) == EVEN ) {
      printf("!!! Error 1 in weakening a cut !!!\n");
      exit(0);
    }
    *crhs = (*crhs - 1) / 2;
    free_info_weak(info_odd_weak);
    return(TRUE);
  }
  else {
    return(FALSE);
  }
}

/* define_cut: construct a cut data structure from a vector of
   coefficients and a right-hand-side */

cut *Cgl012Cut::define_cut(
		int *ccoef, /* coefficients of the cut */
		int crhs /* right hand side of the cut */
		)
{
  int cnzcnt, j;
  cut *v_cut;

  v_cut = reinterpret_cast<cut *> (calloc(1,sizeof(cut)));
  if ( v_cut == NULL ) alloc_error(const_cast<char*>("v_cut"));
  v_cut->crhs = crhs;
  cnzcnt = 0;
  for ( j = 0; j < inp_ilp->mc; j++ ) {
    if ((j & 4095) == 0 &&
        checkTimeLimit("define_cut timeout", "while counting coefficients")) {
      free(v_cut);
      return(NULL);
    }
    if ( ccoef[j] != 0 ) cnzcnt++;
  }
  v_cut->cnzcnt = cnzcnt;
  v_cut->csense = 'L';
  v_cut->cind = reinterpret_cast<int *> (calloc(cnzcnt,sizeof(int)));
  if ( v_cut->cind == NULL ) alloc_error(const_cast<char*>("v_cut->cind"));
  v_cut->cval = reinterpret_cast<int *> (calloc(cnzcnt,sizeof(int)));
  if ( v_cut->cval == NULL ) alloc_error(const_cast<char*>("v_cut->cval"));
  cnzcnt = 0; v_cut->violation = 0.0;
  for ( j = 0; j < inp_ilp->mc; j++ ) {
    if ((j & 4095) == 0 &&
        checkTimeLimit("define_cut timeout", "while materializing coefficients")) {
      free_cut(v_cut);
      return(NULL);
    }
    if ( ccoef[j] != 0 ) {
      v_cut->cind[cnzcnt] = j;
      v_cut->cval[cnzcnt] = ccoef[j];
      v_cut->violation += inp_ilp->xstar[j] * static_cast<double> (ccoef[j]);
      cnzcnt++;
    }
  }
  v_cut->violation -= static_cast<double> (crhs); 
  return(v_cut);
}

/* get_cut: extract a hopefully violated cut from an odd cycle of the
   separation graph */

cut *Cgl012Cut::get_cut(
	     cycle *s_cyc /* shortest odd cycles identified in the separation graph */
	     )
{ 
  int i, e, crhs;
  short int ok;
  /* short int original_parity; */
  double violation;
  /* double original_slack, best_even_slack, best_odd_slack; */
  int *ccoef /*, *vars_to_weak */ ;  
  /* info_weak *info_even_weak, *info_odd_weak; */
  cut *v_cut;
  int ncomb;
  int *comb;
  short int *flag_comb;
#ifndef CGGGGG
  static int iter = 0;
  static double gap, maxgap = 0.0;
#endif

#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&tsi);
  cut_ncalls++;
#endif

  /* compute the cut obtained by adding-up all the constraints 
     corresponding to edges in the cycle, in their non-weak form */

#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&tii);
#endif
  
  ccoef = reinterpret_cast<int *> (calloc(inp_ilp->mc,sizeof(int)));
  if ( ccoef == NULL ) alloc_error(const_cast<char*>("ccoef"));
  ncomb = 0;
  comb = reinterpret_cast<int *> (calloc(inp_ilp->mr,sizeof(int)));
  if ( comb == NULL ) alloc_error(const_cast<char*>("comb"));
  flag_comb = reinterpret_cast<short int *> (calloc(inp_ilp->mr,sizeof(short int)));
  if ( flag_comb == NULL ) alloc_error(const_cast<char*>("flag_comb"));
#if 0
  // no need to as calloc used
  for ( i = 0; i < inp_ilp->mr; i++ ) flag_comb[i] = OUT;

  for ( j = 0; j < inp_ilp->mc; j++ )
    ccoef[j] = 0;
#endif
  crhs = 0;
  for ( e = 0; e < s_cyc->length; e++ ) {
    if ((e & 255) == 0 &&
        checkTimeLimit("get_cut timeout", "while gathering cycle constraints")) {
      free(ccoef);
      free(comb);
      free(flag_comb);
      return(NULL);
    }
    i = (s_cyc->edge_list[e])->constr; 
    if ( i >= 0 && flag_comb[i] != IN) {
      /* the edge is not associated with a bound constraint */
      assert (ncomb<inp_ilp->mr);
      comb[ncomb] = i; ncomb++; flag_comb[i] = IN;
    }
  }
  ok = get_ori_cut_coef(ncomb,comb,ccoef,&crhs,TRUE);

#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&tff);
  coef_time += tff - tii;
#endif
  
  ok = ok && best_cut(ccoef,&crhs,&violation,TRUE,TRUE); 
  if ( ! ok ) {
    free(ccoef);
    free(comb);
    free(flag_comb);
#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
    second_(&tsf);
    cut_time += tsf - tsi;
#endif
    return(NULL);
  }

  v_cut = define_cut(ccoef,crhs);
  if ( v_cut == NULL ) {
    free(ccoef);
    free(comb);
    free(flag_comb);
    return(NULL);
  }
iter++;
  if ( v_cut->violation > violation + EPS || 
       v_cut->violation < violation - EPS ) {
    //printf("Error in violation check\n");
    //printf("v_cut->violation %f  violation %f  gap %f  maxgap (previous) %f\n",
    //	   v_cut->violation,violation,v_cut->violation-violation,maxgap);
    //printf("iter %d\n",iter);
    //exit(0);
    free_cut(v_cut);
    free(ccoef);
    free(comb);
    free(flag_comb);
    errorNo=1;
    return(NULL);
  }
gap = v_cut->violation - violation;
if ( gap < 0.0 ) gap = -gap;
if ( gap > maxgap ) maxgap = gap;
  v_cut->n_of_constr = ncomb;
  v_cut->constr_list = comb;
  v_cut->in_constr_list = flag_comb;

  free(ccoef);

#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&tsf);
  cut_time += tsf - tsi;
#endif

  return(v_cut);
}

/* cut_score: define the score of a (violated) cut */

double Cgl012Cut::cut_score(
		 int *ccoef, /* cut left hand side coefficients */
		 int crhs, /* cut right hand side */
		 double viol, /* cut violation */
		 short int only_viol /* flag which tells whether only an inequality of
			slack smaller than MAX_SLACK is of interest (TRUE)
			otherwise (FALSE) */
		 )
{
  int j, norm;

  /* very simple score: violation divided/multiplied by the lhs norm */
  if ( only_viol && viol < MIN_VIOLATION ) return(-INF);
  norm = 0;
  for ( j = 0; j < p_ilp->mc; j++ ) {
    if ( ccoef[j] != 0 ) norm += ccoef[j] * ccoef[j];
  }
  if ( viol > 0.0 ) return (viol / sqrt(static_cast<double> (norm)));
  else return (viol * sqrt(static_cast<double> (norm)));
}

/* same_cut: check whether two cuts are identical - not too clever
   (assumes the sparse coefficients are sorted by column index) */

short int same_cut(cut *cut1, cut *cut2 /* cuts to be compared */)
{
  int j;

  if ( cut1->cnzcnt != cut2->cnzcnt ) return(FALSE);
  if ( cut1->crhs != cut2->crhs ) return(FALSE);
  if ( cut1->csense != cut2->csense ) return(FALSE);
  for ( j = 0; j < cut1->cnzcnt; j++ ) {
    if ( cut1->cind[j] != cut2->cind[j] ) return(FALSE);
    if ( cut1->cval[j] != cut2->cval[j] ) return(FALSE);
  }
  return(TRUE);
}

/* add_cut_to_list: adds a cut to a list after checking that a copy of
   the same cut is not already in the list - no checking is made about 
   cuts dominating each other or implied by other cuts in the list plus
   the constraints of the original problem */

cut_list *add_cut_to_list(
			  cut *v_cut, /* pointer to the violated cut to be added to the list */
			  cut_list *cuts /* input cut list to be updated */
			  )
{
  int c;

  for ( c = 0; c < cuts->cnum; c++ ) {
    if ( same_cut(v_cut,cuts->list[c]) ) {
      free_cut(v_cut);
      return(cuts);
    }
  }
  cuts->list[cuts->cnum] = v_cut;
  cuts->cnum++;
  return(cuts);
}

/* getcuts: pick the 0-1/2 cuts in the list and give them on output */

void getcuts(
	     cut_list *cuts, /* input cut list */
	     int *cnum, /* number of violated 0-1/2 cuts identified by the procedure */
	     int *cnzcnt, /* overall number of nonzero's in the cuts */
	     int **cbeg, /* starting position of each cut in arrays cind and cval */
	     int **ccnt, /* number of entries of each cut in arrays cind and cval */
	     int **cind, /* column indices of the nonzero entries of the cuts */
	     int **cval, /* values of the nonzero entries of the cuts */
	     int **crhs, /* right hand sides of the cuts */
	     char **csense /* senses of the cuts: 'L', 'G' or 'E' */
	     )
{
  int i, ofsj, count;
  cut *cut_ptr;

  /* allocate the memory for the output vectors */

  (*cnum) = cuts->cnum;
  (*cnzcnt) = 0;
  for ( i = 0; i < cuts->cnum; i++ )
    (*cnzcnt) += (cuts->list[i])->cnzcnt;
  /* if ( (*cbeg) != NULL ) free(*cbeg); */
  (*cbeg) = reinterpret_cast<int *> (calloc((*cnum),sizeof(int)));
  if ( (*cbeg) == NULL ) alloc_error(const_cast<char*>("*cbeg"));
  /* if ( (*ccnt) != NULL ) free(*ccnt); */
  (*ccnt) = reinterpret_cast<int *> (calloc((*cnum),sizeof(int)));
  if ( (*ccnt) == NULL ) alloc_error(const_cast<char*>("*ccnt"));
  /* if ( (*crhs) != NULL ) free(*crhs); */
  (*crhs) = reinterpret_cast<int *> (calloc((*cnum),sizeof(int)));
  if ( (*crhs) == NULL ) alloc_error(const_cast<char*>("*crhs"));
  /* if ( (*csense) != NULL ) free(*csense); */
  (*csense) = reinterpret_cast<char *> (calloc((*cnum),sizeof(char)));
  if ( (*csense) == NULL ) alloc_error(const_cast<char*>("*csense"));
  /* if ( (*cind) != NULL ) free(*cind); */
  (*cind) = reinterpret_cast<int *> (calloc((*cnzcnt),sizeof(int)));
  if ( (*cind) == NULL ) alloc_error(const_cast<char*>("*cind")); 
  /* if ( (*cval) != NULL ) free(*cval); */
  (*cval) = reinterpret_cast<int *> (calloc((*cnzcnt),sizeof(int)));
  if ( (*cval) == NULL ) alloc_error(const_cast<char*>("*cval"));

  /* transfer the cuts information into the output data structures */

  count = 0;
  for ( i = 0; i < cuts->cnum; i++ ) {
    cut_ptr = cuts->list[i];
    (*cbeg)[i] = count; 
    (*ccnt)[i] = cut_ptr->cnzcnt;
    (*crhs)[i] = cut_ptr->crhs;
    (*csense)[i] = cut_ptr->csense;
    for ( ofsj = 0; ofsj < cut_ptr->cnzcnt; ofsj++ ) {
      (*cind)[count] = cut_ptr->cind[ofsj];
      (*cval)[count] = cut_ptr->cval[ofsj];
      count++;
    }
  }
}

/* actual separation subroutines */

/* best_weakening: find the best upper/lower bound weakening of a set
   of variables */

int Cgl012Cut::best_weakening(
		   int n_to_weak, /* number of variables to weaken */
int *vars_to_weak, /* indices of the variables to weaken */
short int original_parity, /* original parity of the constraint to weaken */
double original_slack, /* original slack of the constraint to weaken */
double *best_even_slack, /* best possible slack of a weakened constraint 
			   with even right-hand-side */
double *best_odd_slack, /* best possible slack of a weakened constraint 
			  with odd right-hand-side */
info_weak **info_even_weak, /* weakening information about the best possible
			       even weakened constraint */ 
info_weak **info_odd_weak, /* weakening information about the best possible
			      odd weakened constraint */ 
short int only_odd, /* flag which tells whether only an odd weakening is of
		       interest (TRUE) or both weakenings are (FALSE) */
short int only_viol /* flag which tells whether only an inequality of
			slack smaller than MAX_SLACK is of interest (TRUE)
			otherwise (FALSE) */
		   )
{
  int nweak, cntweak, ofsl, l;
  short int flag_even, flag_odd, ok_even, ok_odd;
  double best_even_e, best_even_o, best_odd_e, best_odd_o;
  short int *type_even_weak, *type_odd_weak,   
	    *switch_even_weak, *switch_odd_weak;

#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&tii);
#endif

  ensureBestWeakeningBufferCapacity(p_ilp->mc);
  type_even_weak = typeEvenWeakBuffer_.data();
  switch_even_weak = switchEvenWeakBuffer_.data();
  type_odd_weak = typeOddWeakBuffer_.data();
  switch_odd_weak = switchOddWeakBuffer_.data();

  if ( original_parity == EVEN ) {
    (*best_even_slack) = original_slack;
    (*best_odd_slack) = INF;
  }
  else {
    (*best_odd_slack) = original_slack;
    (*best_even_slack) = INF;
  }
  nweak = 0;
  for ( ofsl = 0; ofsl < n_to_weak; ofsl++ ) {
    if ((ofsl & 1023) == 0 &&
        checkTimeLimit("best_weakening timeout", "while evaluating weakening options")) {
      return(NONE);
    }
    l = vars_to_weak[nweak];
    if ( p_ilp->possible_weak[l] == NONE ) {
#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
      second_(&tff);
      bw_time += tff - tii;
#endif
      return(NONE);
    }
    else if ( p_ilp->possible_weak[l] == EVEN ) { 
      /* only even weakening of l is possible */
      (*best_even_slack) += p_ilp->loss_even_weak[l];
      type_even_weak[nweak] = p_ilp->type_even_weak[l];
      switch_even_weak[nweak] = FALSE;
      (*best_odd_slack) += p_ilp->loss_even_weak[l];
      type_odd_weak[nweak] = p_ilp->type_even_weak[l];
      switch_odd_weak[nweak] = FALSE;
    }
    else if ( p_ilp->possible_weak[l] == ODD ) { 
      /* only odd weakening of l is possible */
      best_even_e = (*best_even_slack);
      best_odd_o = (*best_odd_slack);
      (*best_even_slack) = best_odd_o + p_ilp->loss_odd_weak[l];
      type_even_weak[nweak] = p_ilp->type_odd_weak[l];
      switch_even_weak[nweak] = TRUE;
      (*best_odd_slack) = best_even_e + p_ilp->loss_odd_weak[l];
      type_odd_weak[nweak] = p_ilp->type_odd_weak[l];
      switch_odd_weak[nweak] = TRUE;
    }
    else {
      /* both weakenings of l are possible */
      best_even_e = (*best_even_slack) + p_ilp->loss_even_weak[l];
      best_even_o = (*best_odd_slack) + p_ilp->loss_odd_weak[l];
      best_odd_e = (*best_odd_slack) + p_ilp->loss_even_weak[l];
      best_odd_o = (*best_even_slack) + p_ilp->loss_odd_weak[l];
      if ( best_even_e <= best_even_o ) {
	(*best_even_slack) = best_even_e;
	type_even_weak[nweak] = p_ilp->type_even_weak[l];
	switch_even_weak[nweak] = FALSE;
      }
      else {
	(*best_even_slack) = best_even_o;
	type_even_weak[nweak] = p_ilp->type_odd_weak[l];
	switch_even_weak[nweak] = TRUE;
      }
      if ( best_odd_e <= best_odd_o ) {
	(*best_odd_slack) = best_odd_e;
	type_odd_weak[nweak] = p_ilp->type_even_weak[l];
	switch_odd_weak[nweak] = FALSE;
      }
      else {
	(*best_odd_slack) = best_odd_o;
	type_odd_weak[nweak] = p_ilp->type_odd_weak[l];
	switch_odd_weak[nweak] = TRUE;
      }
    }
    if ( ( only_viol ) &&
	 ( (*best_even_slack) > MAX_SLACK - EPS ) &&
	 ( (*best_odd_slack) > MAX_SLACK - EPS ) ) {
#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
      second_(&tff);
      bw_time += tff - tii;
#endif
      return(NONE);
    }
    nweak++;
  }

  /* construct the weakening vectors associated with the best
     even and odd pairs (if the associated slack is not too big) */

  if ( (! only_odd) && 
       ( ( (*best_even_slack) <= MAX_SLACK - EPS ) ||
       ( (! only_viol) && (*best_even_slack) <= INF - EPS ) ) ) {
    ok_even = TRUE;
    (*info_even_weak) = alloc_info_weak(nweak);
    (*info_even_weak)->nweak = nweak;
    flag_even = EVEN;
    cntweak = nweak; 
    for ( ofsl = n_to_weak - 1; ofsl >= 0; ofsl-- ) {
      if (((n_to_weak - 1 - ofsl) & 1023) == 0 &&
          checkTimeLimit("best_weakening timeout", "while building even weakening")) {
        free_info_weak(*info_even_weak);
        *info_even_weak = NULL;
        return(NONE);
      }
      cntweak--;
      (*info_even_weak)->var[cntweak] = vars_to_weak[ofsl];
      if ( flag_even == EVEN ) {
	(*info_even_weak)->type[cntweak] = type_even_weak[cntweak];
	if ( switch_even_weak[cntweak] ) flag_even = ODD;
      }
      else {
	(*info_even_weak)->type[cntweak] = type_odd_weak[cntweak];
	if ( switch_odd_weak[cntweak] ) flag_even = EVEN;
      }
    }
  }
  else ok_even = FALSE;

  if ( ( (*best_odd_slack) <= MAX_SLACK - EPS ) || 
       ( (! only_viol) && (*best_odd_slack) <= INF - EPS ) ) {
    ok_odd = TRUE;
    (*info_odd_weak) = alloc_info_weak(nweak);
    (*info_odd_weak)->nweak = nweak;
    flag_odd = ODD;
    cntweak = nweak; 
    for ( ofsl = n_to_weak - 1; ofsl >= 0; ofsl-- ) {
      if (((n_to_weak - 1 - ofsl) & 1023) == 0 &&
          checkTimeLimit("best_weakening timeout", "while building odd weakening")) {
        if (ok_even)
          free_info_weak(*info_even_weak);
        *info_even_weak = NULL;
        free_info_weak(*info_odd_weak);
        *info_odd_weak = NULL;
        return(NONE);
      }
      cntweak--;
      (*info_odd_weak)->var[cntweak] = vars_to_weak[ofsl];
      if ( flag_odd == EVEN ) {
	(*info_odd_weak)->type[cntweak] = type_even_weak[cntweak];
	if ( switch_even_weak[cntweak] ) flag_odd = ODD;
      }
      else {
	(*info_odd_weak)->type[cntweak] = type_odd_weak[cntweak];
	if ( switch_odd_weak[cntweak] ) flag_odd = EVEN;
      }
    }
  }
  else ok_odd = FALSE;

#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&tff);        
  bw_time += tff - tii;
#endif

  if ( ok_odd && ok_even ) return(BOTH);
  if ( ok_even ) return(EVEN);
  if ( ok_odd ) return(ODD);
  return(NONE);
}

/* basic_separation: try to identify violated 0-1/2 cuts by using the 
   original procedure described in Caprara and Fischetti's MP paper */

cut_list *Cgl012Cut::basic_separation()
{
  int i, j, k, l, begi, special, ofsj, ofsk, ofsl, n_to_weak, c;
  short int parity, original_parity, ok_weak;
  double weight, original_slack, best_even_slack, best_odd_slack;
  int *vars_to_weak;
  info_weak *info_even_weak, *info_odd_weak, *i_weak;
  separation_graph *sep_graph;
  auxiliary_graph *aux_graph;
  cycle_list *short_cycle_list;
  cut *violated_cut;
  cut_list *out_cuts;

  /* construct the separation graph by the standard weakening procedure */
  
#ifdef PRINT
  print_parity_ilp();
#endif

#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&td);
#endif
#ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_TIMING
  printf("... time elapsed at the beginning of basic_separation: %f\n",td - tti);
#endif
#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&td);
#endif
  #ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_TIMING
  printf("... time elapsed before initialize_sep_graph: %f\n",td - tti);
  #endif
  sep_graph = initialize_sep_graph();
  if ( sep_graph == NULL ) return(initialize_cut_list(1));
  special = p_ilp->mc;
#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&td);
#endif
#ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_TIMING
  printf("... time elapsed before weakening: %f\n",td - tti);
#endif

  /* edges associated with actual constraints in the ILP */

  for ( i = 0; i < p_ilp->mr; i++ ) {
    if ((i & 255) == 0 &&
        checkTimeLimit("basic_separation timeout", "during separation-graph construction")) {
      free_sep_graph(sep_graph);
      return(initialize_cut_list(1));
    }
    if ( ! p_ilp->row_to_delete[i] ) {
#ifdef CGL_ZH_ADVANCED_DEBUG_CONSTRAINT_PROFILING
      const double zhRowStartSeconds = CoinGetTimeOfDay();
#endif
      begi = p_ilp->mtbeg[i];
      int zhRowPairCount = 0;
      if ( p_ilp->mtcnt[i] <= 2 )
        zhRowPairCount = (p_ilp->mtcnt[i] >= 1) ? 1 : 0;
      else
        zhRowPairCount = (p_ilp->mtcnt[i] * (p_ilp->mtcnt[i] - 1)) / 2;
      if ( rowMaxPairCount_ >= 0 && zhRowPairCount > rowMaxPairCount_ )
        continue;
      if ( rowMaxFractionalCount_ >= 0 ) {
        int zhFractionalHistogram[10];
        int zhFractionalCount = 0;
        int zhVeryFractionalCount = 0;
        double zhFractionalSum = 0.0;
        cglZeroHalfCollectFractionality(inp_ilp, i, zhFractionalHistogram,
          &zhFractionalCount, &zhFractionalSum, &zhVeryFractionalCount);
        if ( zhFractionalCount > rowMaxFractionalCount_ )
          continue;
      }
      if ( p_ilp->mtcnt[i] == 1 ) {
	/* row with one odd entry only: edge j -- special */
	weight = p_ilp->slack[i];
	if ( weight < MAX_SLACK - EPS ) {
	  j = p_ilp->mtind[begi];
	  parity = p_ilp->mrhs[i];
	  i_weak = alloc_info_weak(0);
	  sep_graph = update_weight_sep_graph
	    (j,special,weight,parity,i,i_weak,sep_graph);
	}
      }
      else if ( p_ilp->mtcnt[i] == 2 ) {
	/* row with two odd entries only: edge j -- k */
	weight = p_ilp->slack[i];
	if ( weight < MAX_SLACK - EPS ) {
	  j = p_ilp->mtind[begi];
	  k = p_ilp->mtind[begi+1];
	  parity = p_ilp->mrhs[i];
	  i_weak = alloc_info_weak(0);
	  sep_graph = update_weight_sep_graph
	    (j,k,weight,parity,i,i_weak,sep_graph);
	}
      }
      else {
	/* row with three or more odd entries: weakening for all 1's pairs */
	for ( ofsj = 0; ofsj < p_ilp->mtcnt[i]; ofsj++ ) {
	  for ( ofsk = ofsj + 1; ofsk < p_ilp->mtcnt[i]; ofsk++ ) {
            if (((ofsj * p_ilp->mtcnt[i]) + ofsk) % 16384 == 0 &&
                checkTimeLimit("basic_separation timeout", "during row weakening")) {
#ifdef CGL_ZH_ADVANCED_DEBUG_CONSTRAINT_PROFILING
              zhProfileRecordRow(i, CoinGetTimeOfDay() - zhRowStartSeconds, zhRowPairCount);
#endif
              free_sep_graph(sep_graph);
              return(initialize_cut_list(1));
            }
	    /* edge(s) j -- k */
	    j = p_ilp->mtind[begi+ofsj];
	    k = p_ilp->mtind[begi+ofsk];
	    original_slack = p_ilp->slack[i];
	    original_parity = p_ilp->mrhs[i];
	    n_to_weak = 0;
	    ensureVarsToWeakBufferCapacity(inp_ilp->mc);
	    vars_to_weak = varsToWeakBuffer_.data();
	    for ( ofsl = 0; ofsl < p_ilp->mtcnt[i]; ofsl++ ) 
	      if ( ofsl != ofsj && ofsl != ofsk ) {
		l = p_ilp->mtind[begi+ofsl];
		vars_to_weak[n_to_weak] = l;
		n_to_weak++;
	      }
	    ok_weak = best_weakening(n_to_weak,vars_to_weak,
				     original_parity,original_slack,
				     &best_even_slack,&best_odd_slack,
				     &info_even_weak,&info_odd_weak,
				     FALSE,TRUE); 
            if ( timeLimitReached() ) {
#ifdef CGL_ZH_ADVANCED_DEBUG_CONSTRAINT_PROFILING
              zhProfileRecordRow(i, CoinGetTimeOfDay() - zhRowStartSeconds, zhRowPairCount);
#endif
              free_sep_graph(sep_graph);
              return(initialize_cut_list(1));
            }
	    if ( ok_weak == NONE ) goto EXITJK;
	    if ( ok_weak == BOTH || ok_weak == EVEN ) {
	      if ( best_even_slack < MAX_SLACK - EPS ) {
		weight = best_even_slack; parity = EVEN;
		sep_graph = update_weight_sep_graph
		  (j,k,weight,parity,i,info_even_weak,sep_graph);
	      }
	    }
	    if ( ok_weak == BOTH || ok_weak == ODD ) {
	      if ( best_odd_slack < MAX_SLACK - EPS ) {
		weight = best_odd_slack; parity = ODD;
		sep_graph = update_weight_sep_graph
		  (j,k,weight,parity,i,info_odd_weak,sep_graph);
	      }
	    }
 EXITJK:;  }
	}
      }
#ifdef CGL_ZH_ADVANCED_DEBUG_CONSTRAINT_PROFILING
      zhProfileRecordRow(i, CoinGetTimeOfDay() - zhRowStartSeconds, zhRowPairCount);
#endif
    }
  }

  /* edges associated with the bound constraints (probably useless
     but necessary in some cases */

  for ( j = 0; j < p_ilp->mc; j++ ) {
    if ( ! p_ilp->col_to_delete[j] ) {
      weight = p_ilp->xstar[j] - inp_ilp->vlb[j];
      if ( weight < MAX_SLACK - EPS ) {
	parity = mod2(inp_ilp->vlb[j]);
	i_weak = alloc_info_weak(0);
	sep_graph = update_weight_sep_graph
	  (j,special,weight,parity,NONE,i_weak,sep_graph);
      }
      weight = inp_ilp->vub[j] - p_ilp->xstar[j];
      if ( weight < MAX_SLACK - EPS ) {
	parity = mod2(inp_ilp->vub[j]);
	i_weak = alloc_info_weak(0);
	sep_graph = update_weight_sep_graph
	  (j,special,weight,parity,NONE,i_weak,sep_graph);
      }
    }
  }   
#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&tf);
  weak_time += tf - ti;
#endif

  /* construct the auxiliary graph for the shortest path computation 
     and compute the smallest cost odd cycle visiting each node -
     this part is strongly dependent on the data structure used by
     the shortest path subroutine used */ 
     
#ifdef PRINT
  print_sep_graph(sep_graph);
#endif

#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&ti);
#endif
#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&td);
#endif
  #ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_TIMING
  printf("... time elapsed before define_aux_graph: %f\n",td - tti);
  #endif
  aux_graph = define_aux_graph(sep_graph);
#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&tf);
  aux_time += tf - ti;
#endif
#ifdef PRINT
  print_aux_graph(aux_graph);
#endif
  
/* exit(1); */
  
#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&td);
#endif
  #ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_TIMING
  printf("... time elapsed before cycles and cuts: %f\n",td - tti);
  printf("%d nodes on list\n",sep_graph->nnodes);
  #endif
  out_cuts = initialize_cut_list(MAX_CUTS);
  for ( j = 0; j < sep_graph->nnodes; j++ ) {
    if ((j & 15) == 0 &&
        checkTimeLimit("basic_separation timeout", "during odd-cycle search"))
      goto EXIT_CUTS;
    short_cycle_list = get_shortest_odd_cycle_list(j,sep_graph,aux_graph);
    if ( timeLimitReached() ) {
      free_cycle_list(short_cycle_list);
      goto EXIT_CUTS;
    }
    if ( short_cycle_list == NULL ) goto EXIT_NODE;
#ifdef PRINT
    print_cycle_list(short_cycle_list);
#endif
    for ( c = 0; c < short_cycle_list->cnum; c++ ) {
      violated_cut = get_cut(short_cycle_list->list[c]);
      if ( timeLimitReached() ) {
        free_cycle_list(short_cycle_list);
        goto EXIT_CUTS;
      }
      if ( violated_cut == NULL ) {
	if (!errorNo)
	  continue;
	else
	  break;
      }
#ifdef PRINT
      print_cut(violated_cut);
#endif
      if ( violated_cut->violation > MIN_VIOLATION + EPS ) {
	/* violated 0-1/2 cut found */  
#ifdef CGL_ZH_ADVANCED_DEBUG_CONSTRAINT_PROFILING
        const int previousCutCount = out_cuts->cnum;
#endif
	out_cuts = add_cut_to_list(violated_cut,out_cuts);  
#ifdef CGL_ZH_ADVANCED_DEBUG_CONSTRAINT_PROFILING
        if (out_cuts->cnum > previousCutCount)
          zhProfileRecordCut(violated_cut);
#endif
	if ( out_cuts->cnum >= MAX_CUTS ) {
	  free_cycle_list(short_cycle_list);
	  goto EXIT_CUTS; 
	}
      }
      else free_cut(violated_cut);
    }
    /* remove the current node from the auxiliary graph */
EXIT_NODE:
    aux_graph = cancel_node_aux_graph(j,aux_graph); 
    free_cycle_list(short_cycle_list);
  }

EXIT_CUTS:
  free_sep_graph(sep_graph);
  free_aux_graph(aux_graph);

#ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_CUTS
  print_cut_list(out_cuts);
#endif
#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&td);
#endif
  #ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_TIMING
  printf("... time elapsed at the end of basic_separation: %f\n",td - tti);
  #endif

  return(out_cuts);
}

/*
  012cut: main procedure for 0-1/2 cut separation
  first release: Aug 12 1996
  last revision: Jun 10 1997
*/  

/* static data structures for log information about separation */

#ifndef CGGGGG
static int sep_iter = 0; /* number of the current separation iteration */
//#define POOL 
#ifdef POOL 
static pool_cut_list *pool = NULL; /* information about the cuts separated
				      so far, used to decide when they should 
				      be added to the current LP */
#endif
static log_var **vlog = NULL; /* information about the value attained
				  by the variables in the last iterations,
				  used to possibly set to 0 some coefficient
				  > 0 in a cut to be added */ 
static bool aggr; /* flag saying whether as many cuts as possible are required
		   from the separation procedure (TRUE) or not (FALSE) */
#endif
/* include the reactive local search heuristic */

//was #include "Cgltabu_012.c"
//start include "Cgltabu_012.c"

#define MAX_TABU_ITER 100
#define NUM_HASH_ENTRIES MAX_CUT_POOL
/* initial length of the tabu list */
#define IN_PROHIB_PERIOD 3
#define MAX_TIME_FACTOR 3

/* data structure for the current local search solution */

typedef struct {
int n_of_constr; /* number of constraints in the current cut */
short int *in_constr_list; /* flag saying whether a given constraint is
			      in the list of constraints of the cut (IN)
			      or not (OUT) */
int *non_weak_coef; /* coefficients of the cut before weakening */
int non_weak_rhs; /* coefficient of the rhs before weakening */
double slack_sum; /* sum of the slacks of the constraints in the cut */
double min_weak_loss; /* minimum loss by weakening the non even 
			 coefficients */
int one_norm; /* 1-norm of the lhs, i.e. sum of the absolute values of
		 the coefficients */
short int ok; /* logical flag telling whether the cut could be weakened
		 to a 0-1/2 cut or not - if false the two fields below
		 have no meaning */
int *coef; /* actual coefficients of the cut */
int rhs; /* actual rhs of the cut */
double violation; /* violation of the cut */
} tabu_cut;

/* data structure for the hash table used in memory reaction */

typedef struct h_e {
int n_of_el; /* number of components to be considered */
short int *flag_vect; /* vector of flags for the components */
int last_vis; /* last iteration when this element was visited */
struct h_e *next; /* pointer to the next element in the hash chain */
} hash_element;

typedef hash_element **hash_table;

/* global variables for local search */

static int n; /* number of variables in the ILP */
static int m; /* number of constraints in the ILP */
static int it; /* number of tabu search iterations so far */
static tabu_cut *cur_cut; /* information about the current cut in local search */
static int *last_moved; /* last iteration when a given constraint was added/
		    deleted from the list of constraints of the cut */
static int last_it_add; /* last iteration when a cut was added to the list */
static int last_it_restart; /* last iteration when a restart was performed */
static int prohib_period; /* current prohibition period */
static int last_prohib_period_mod; /* last iteration where prohibition period was modified */
static hash_table hash_tab; /* hash table */
static int A; /* parameter A in Battiti and Protasi */
static int B; /* parameter B in Battiti and Protasi */
static float elapsed_time; /* time elapsed since the beginning of the current
		       tabu search call */

/* clear_cur_cut: clear the current solution (no constraint in the cut) */

void clear_cur_cut()
{
  int i, j;

  cur_cut->n_of_constr = 0;
  cur_cut->rhs = 0;
  cur_cut->non_weak_rhs = 0;
  cur_cut->violation = 0.0;
  cur_cut->slack_sum = 0.0;
  cur_cut->min_weak_loss = 0.0;
  cur_cut->one_norm = 0;
  for ( j = 0; j < n; j++ ) {
    cur_cut->coef[j] = 0;
    cur_cut->non_weak_coef[j] = 0;
  }
  for ( i = 0; i < m; i++ ) {
    cur_cut->in_constr_list[i] = OUT;
  }
  cur_cut->ok = FALSE;
}

/* initialize_cur_cut: allocate the memory for cur_cut */

void initialize_cur_cut() 
{
  cur_cut = reinterpret_cast<tabu_cut *> (calloc(1,sizeof(tabu_cut)));
  if ( cur_cut == NULL ) alloc_error(const_cast<char*>("cur_cut"));
  cur_cut->coef = reinterpret_cast<int *> (calloc(n,sizeof(int)));
  if ( cur_cut->coef == NULL ) alloc_error(const_cast<char*>("cur_cut->coef"));
  cur_cut->non_weak_coef = reinterpret_cast<int *> (calloc(n,sizeof(int)));
  if ( cur_cut->non_weak_coef == NULL ) alloc_error(const_cast<char*>("cur_cut->non_weak_coef"));
  cur_cut->in_constr_list = reinterpret_cast<short int *> (calloc(m,sizeof(short int)));
  if ( cur_cut->in_constr_list == NULL ) alloc_error(const_cast<char*>("cur_cut->in_constr_list"));
  clear_cur_cut();
}

/* free_cur_cut: free the memory for cur_cut */

void free_cur_cut()
{
  free(cur_cut->coef);
  free(cur_cut->non_weak_coef);
  free(cur_cut->in_constr_list);
  free(cur_cut);
}

#ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_TABU
/* print_cur_cut: display cur_cut on output */

void Cgl012Cut::print_cur_cut()
{
  int i, j; 

  printf("iteration %d  prohib_period %d\n",it,prohib_period);
  printf("\n content of cur_cut data structure: n_of_constr = %d, ok = %d\n", cur_cut->n_of_constr, cur_cut->ok);
  for ( i = 0; i < m; i++ ) 
    if ( cur_cut->in_constr_list[i] == IN ) 
      printf("constr. %d\n",i);
  /*
  printf(" list of constraints:\n");
  for ( i = 0; i < m; i++ ) 
    if ( cur_cut->in_constr_list[i] == IN ) 
      print_constr(i);
  print_int_vect("non_weak_coef",cur_cut->non_weak_coef,n);
  printf(" non_weak_rhs = %d\n",cur_cut->non_weak_rhs);
  print_int_vect("coef",cur_cut->coef,n);
  printf(" rhs = %d\n",cur_cut->rhs);
  */
  for ( j = 0 /* , viol = - (double) cur_cut->rhs */ ; j < n; j++ ) 
    if ( ( p_ilp->xstar[j] > ZERO || p_ilp->xstar[j] < -ZERO ) && cur_cut->non_weak_coef[j] != 0 ) {
      printf("var. %d  xstar %f  non_weak_coef %d  coef %d\n", j, p_ilp->xstar[j], cur_cut->non_weak_coef[j], cur_cut->coef[j]);
      /* viol += p_ilp->xstar[j] * cur_cut->coef[j]; */
    }
  printf("rhs %d  viol %f  slack_sum %f  min_weak_loss %f  one_norm %d\n", 
	 cur_cut->rhs, cur_cut->violation, cur_cut->slack_sum,
	 cur_cut->min_weak_loss, cur_cut->one_norm);
}  
#endif  
/* same_short_vect: check whether two short int vectors have the same content */

short int same_short_vect(
			  int n_of_el, /* number of components in the vectors */
			  short int *vec_1, 
			  short int *vec_2 /* vectors to be checked */
			  )
{
  int i;
  for ( i = 0; i < n_of_el; i++ ) 
    if ( vec_1[i] != vec_2[i] ) return(FALSE);
  return(TRUE);
}

/* initialize_hash_table: allocate the memory for the hash table */

void initialize_hash_table()
{
  int i;
  hash_tab = reinterpret_cast<hash_element **> (calloc(NUM_HASH_ENTRIES,sizeof(hash_element *)));
  if ( hash_tab == NULL ) alloc_error(const_cast<char*>("hash_tab"));
  for ( i = 0; i < NUM_HASH_ENTRIES; i++ ) hash_tab[i] = NULL;
}

/* clear_hash_table: clear the current hash table */

void clear_hash_table()
{
  int i; 
  hash_element *hash_ptr, *hash_el;

  for ( i = 0; i < NUM_HASH_ENTRIES; i++ ) {
    if ( hash_tab[i] != NULL ) {
      hash_ptr = hash_tab[i];
      do {
	hash_el = hash_ptr->next;
	free(hash_ptr->flag_vect);
	free(hash_ptr);
	hash_ptr = hash_el;
      } while ( hash_ptr != NULL );
      hash_tab[i] = NULL;
    } 
  }
}

/* free_hash_table: deallocate the memory for the hash table */

void free_hash_table()
{
  clear_hash_table();
  free(hash_tab);
}

/* hash_addr: compute the hash address associated with the current cut */

int hash_addr(
	      int n_of_el, /* number of elements to be considered */
	      short int *flag_vect /* vector of flags for the elements */
	      )
{
  int i, addr;

  /* very simple algorithm: just add-up the squared indices of the IN elements */
  addr = 0;
  for ( i = 0; i < n_of_el; i++ ) 
    if ( flag_vect[i] == IN ) addr += i * i;
  return(addr % NUM_HASH_ENTRIES);
}

/* hash_search: search for the current cut in the hash list of all cuts -
   if found return TRUE and update the last iteration the cut was found */

short int hash_search(int *cyc_len /* length of the cycle if the current cut is found */)
{ 
  int addr;
  hash_element *hash_el;

  addr = hash_addr(m,cur_cut->in_constr_list);
  hash_el = hash_tab[addr];
  while ( hash_el != NULL ) {
    if ( same_short_vect(m,cur_cut->in_constr_list,hash_el->flag_vect) ) {
      *cyc_len = it - hash_el->last_vis;
      hash_el->last_vis = it;
      return(TRUE);
    }
    hash_el = hash_el->next;
  }
  return(FALSE);
}

/* hash_insert: insert a new cut in the hash list of all cuts */

void hash_insert()
{
  int addr, i;
  hash_element *hash_el, *hash_ptr;

  addr = hash_addr(m,cur_cut->in_constr_list);
  hash_el = reinterpret_cast<hash_element *> (calloc(1,sizeof(hash_element)));
  if ( hash_el == NULL ) alloc_error(const_cast<char*>("hash_el"));
  hash_el->n_of_el = m;
  hash_el->last_vis = it;
  hash_el->next = NULL;
  hash_el->flag_vect = reinterpret_cast<short int *> (calloc(m,sizeof(short int)));
  if ( hash_el->flag_vect == NULL ) alloc_error(const_cast<char*>("hash_el->flag_vect"));
  for ( i = 0; i < m; i++ )
    hash_el->flag_vect[i] = cur_cut->in_constr_list[i];
  if ( hash_tab[addr] == NULL ) 
    hash_tab[addr] = hash_el;
  else {
    hash_ptr = hash_tab[addr];
    while ( hash_ptr->next != NULL ) {
#if 0
      /* this check can be omitted to save time */
      if ( same_short_vect(m,cur_cut->in_constr_list,hash_ptr->flag_vect) ) {
	printf("attempt to insert in the hash an already present cut\n");
	exit(0);
      }
#endif
      hash_ptr = hash_ptr->next;
    }
    hash_ptr->next = hash_el;
  }
}        

/* increase_prohib_period: implemented as in Battiti and Protasi */

void increase_prohib_period()
{
  if ( prohib_period * 1.1 > prohib_period + 1 ) 
    if ( prohib_period * 1.1 < m - 2 ) prohib_period = 
					 static_cast<int> (prohib_period*1.1);
    else prohib_period = m - 2;
  else
    if ( prohib_period + 1 < m - 2 ) prohib_period += 1;
    else prohib_period = m - 2;
  last_prohib_period_mod = it;
}

/* decrease_prohib_period: implemented as in Battiti and Protasi */

void decrease_prohib_period()
{
  if ( prohib_period * 0.9 < prohib_period - 1 ) 
    if ( prohib_period * 0.9 > IN_PROHIB_PERIOD ) prohib_period = static_cast<int> (prohib_period* 0.9);
    else prohib_period = IN_PROHIB_PERIOD;
  else
    if ( prohib_period - 1 > IN_PROHIB_PERIOD ) prohib_period -= 1;
    else prohib_period = IN_PROHIB_PERIOD;
  last_prohib_period_mod = it;
}

/* allowed: check if moving (adding/deleting) a given constraint 
   is not a tabu move */

short int allowed(int i /* constraint to be checked */)
{
  if ( last_moved[i]  < it - prohib_period ) {
    if ( cur_cut->in_constr_list[i] == IN ) {
      if ( cur_cut->n_of_constr > 1 ) return(TRUE);
      else return(FALSE);
    }
    else {
      if ( cur_cut->n_of_constr < m - 1 ) return(TRUE);
      else return(FALSE);
    }
  }
  else return(FALSE);
}

/* in_cur_cut: check whether a given constraint is in the list of
   constraints defining the current cut */

short int in_cur_cut(int i /* constraint to be checked */)
{
  if ( cur_cut->in_constr_list[i] == OUT ) return(FALSE);
  else return(TRUE);
}

/* tabu_score: define the score of a potential new cut */

double tabu_score(
		  int *ccoef, /* cut left hand side coefficients */
		  int crhs, /* cut right hand side */
		  double viol, /* cut violation */
double norm /* cut norm - 1-norm is used below */
		  )
{
  /* very simple score: violation divided/multiplied by the lhs 1-norm */
   
  if ( norm == 0 ) norm = 1;
  if ( viol > 0.0 ) return (viol / norm);
  else return (viol * norm);
}

/* score_by_moving: compute the score of the best cut obtainable from 
   the current local search solution by inserting/deleting a constraint */

double Cgl012Cut::score_by_moving(
		       int i, /* constraint to be moved */
		       short int itype, /* type of move - ADD or DEL */
		       double thresh /* minimum value of an interesting score */
		       )
{
#define PENALTY_NON_WEAKABLE -1.0
#define FAST_SCORE_EVAL 1
  int j, begi, gcdi, ofsj, ij, crhs, one_norm, support_inter; 
  short int flag_gt, ok;
  double slack_sum, weak_loss, score, viol;
  int *new_coef, *ccoef;
  begi = inp_ilp->mtbeg[i];
  gcdi = p_ilp->gcd[i];
  
  /* fast check - optimistic evaluation of the score */

  slack_sum = cur_cut->slack_sum;

  if ( itype == ADD ) slack_sum += p_ilp->slack[i] / static_cast<double> (gcdi);
  else slack_sum -= p_ilp->slack[i] / static_cast<double> (gcdi);

  viol = ( 1.0 - slack_sum ) / 2.0;

  score = tabu_score(NULL,0,viol,1.0);
/*
printf("Score estimate 1 %f   Threshold %f\n",score,thresh);
*/
  if ( score < thresh + ZERO ) {
    return (score);
  }

  /* discard the cuts that have empty support intersection with the 
     current one */

  support_inter = 0;

  for ( ofsj = 0, ij = begi; ofsj < inp_ilp->mtcnt[i]; ofsj++, ij++ ) 
    if ( cur_cut->non_weak_coef[inp_ilp->mtind[ij]] !=0 ) support_inter++;

  if ( support_inter == 0 ) return(-INF);

  /* compute the new cut coefficients and rhs (before weakening) */

  new_coef = reinterpret_cast<int *> (calloc(inp_ilp->mtcnt[i],sizeof(int)));
  if ( new_coef == NULL ) alloc_error(const_cast<char*>("new_coef"));

  if ( ( itype == ADD && inp_ilp->msense[i] != 'G' ) ||
    ( itype == DEL && inp_ilp->msense[i] == 'G' ) ) flag_gt = 1;
  else flag_gt = -1;
  if ( flag_gt == 1 ) {
    if ( gcdi == 1 ) {
      for ( ofsj = 0, ij = begi; ofsj < inp_ilp->mtcnt[i]; ofsj++, ij++ ) {
	j = inp_ilp->mtind[ij];
	new_coef[ofsj] = cur_cut->non_weak_coef[j] + inp_ilp->mtval[ij];
      }
      crhs = cur_cut->non_weak_rhs + inp_ilp->mrhs[i];
    }
    else {
      for ( ofsj = 0, ij = begi; ofsj < inp_ilp->mtcnt[i]; ofsj++, ij++ ) {
	j = inp_ilp->mtind[ij];
	new_coef[ofsj] = cur_cut->non_weak_coef[j] + inp_ilp->mtval[ij] / gcdi;
      }
      crhs = cur_cut->non_weak_rhs + inp_ilp->mrhs[i] / gcdi;
    }
  }
  else {
    if ( gcdi == 1 ) {
      for ( ofsj = 0, ij= begi; ofsj < inp_ilp->mtcnt[i]; ofsj++, ij++ ) {
	j = inp_ilp->mtind[ij];
	new_coef[ofsj] = cur_cut->non_weak_coef[j] - inp_ilp->mtval[ij];
      }
      crhs = cur_cut->non_weak_rhs - inp_ilp->mrhs[i];
    }
    else {
      for ( ofsj = 0, ij = begi; ofsj < inp_ilp->mtcnt[i]; ofsj++, ij++ ) {
	j = inp_ilp->mtind[ij];
	new_coef[ofsj] = cur_cut->non_weak_coef[j] - inp_ilp->mtval[ij] / gcdi;
      }
      crhs = cur_cut->non_weak_rhs - inp_ilp->mrhs[i] / gcdi;
    }
  }

  /* other - relatively fast - check by optimistic evaluation of the 
     cut score */

  weak_loss = cur_cut->min_weak_loss;
  one_norm = cur_cut->one_norm;
  for ( ofsj = 0, ij = begi; ofsj < inp_ilp->mtcnt[i]; ofsj++, ij++ ) {
    j = inp_ilp->mtind[ij];
    if ( cur_cut->coef[j] > 0 ) one_norm -= cur_cut->coef[j];
    else one_norm += cur_cut->coef[j];
    if ( new_coef[ofsj] >= 2 ) one_norm += new_coef[ofsj] / 2;
    else one_norm -= new_coef[ofsj] / 2;
    if ( mod2(cur_cut->non_weak_coef[j]) == ODD ) {
      if ( mod2(new_coef[ofsj]) == EVEN ) 
	weak_loss -= p_ilp->min_loss_by_weak[j];
    }
    else {
      if ( mod2(new_coef[ofsj]) == ODD ) 
	weak_loss += p_ilp->min_loss_by_weak[j];
    }
  }
    
  viol = ( 1.0 - slack_sum - weak_loss ) / 2.0;

  score = tabu_score(NULL,0,viol,static_cast<double> (one_norm));
/*
printf("Score estimate 2 %f   Threshold %f\n",score,thresh);
*/
  if ( score < thresh + ZERO || FAST_SCORE_EVAL ) {

    free(new_coef);
    return (score);
  }
  /* get the actual cut coefficients and the violation of the 
     best cut obtainable trough weakening */

  ccoef = reinterpret_cast<int *> (calloc(n,sizeof(int)));
  if ( ccoef == NULL ) alloc_error(const_cast<char*>("ccoef"));
  for ( j = 0; j < n; j++ ) ccoef[j] = cur_cut->non_weak_coef[j];
  for ( ofsj = 0; ofsj < inp_ilp->mtcnt[i]; ofsj++ ) {
    ij = begi + ofsj;
    j = inp_ilp->mtind[ij];
    ccoef[j] = new_coef[ofsj];
  }

  ok = best_cut(ccoef,&crhs,&viol,FALSE,FALSE);
  if ( ok ) score = tabu_score(ccoef,crhs,viol,static_cast<double> (one_norm));
  else {
    viol = ( - slack_sum - weak_loss ) / 2.0;
    score = PENALTY_NON_WEAKABLE + tabu_score(ccoef,crhs,viol,static_cast<double> (one_norm));
  }

/*
printf("Score actual %f   Threshold %f\n",score,thresh);
*/

  free(new_coef);
  free(ccoef);
  return(score);
}

/* modify_current: update the current local search solution by inserting/
   deleting a constraint */

void Cgl012Cut::modify_current(
		    int i, /* constraint to be moved */
		    short int itype /* type of move - ADD or DEL */
		    )
{
  int j, begi, gcdi, ofsj, ij; 
  short int flag_gt;

  if ( itype == ADD ) {
    cur_cut->n_of_constr++;
    cur_cut->in_constr_list[i] = IN;
  }
  else {
    cur_cut->n_of_constr--;
    cur_cut->in_constr_list[i] = OUT;
  }
  last_moved[i] = it;

  /* compute the new cut coefficients and rhs (before weakening) */

  if ( ( itype == ADD && inp_ilp->msense[i] != 'G' ) ||
    ( itype == DEL && inp_ilp->msense[i] == 'G' ) ) flag_gt = 1;
  else flag_gt = -1;
  begi = inp_ilp->mtbeg[i];
  gcdi = p_ilp->gcd[i];
  for ( ofsj = 0; ofsj < inp_ilp->mtcnt[i]; ofsj++ ) {
    ij = begi + ofsj;
    j = inp_ilp->mtind[ij];
    /* the '*' and '/' operations can be saved by writing some more code ... */
    cur_cut->non_weak_coef[j] += flag_gt * (inp_ilp->mtval[ij] / gcdi);
  }
  cur_cut->non_weak_rhs += flag_gt * (inp_ilp->mrhs[i] / gcdi);

  if ( itype == ADD ) 
    cur_cut->slack_sum += p_ilp->slack[i] / static_cast<double> (gcdi);
  else 
    cur_cut->slack_sum -= p_ilp->slack[i] / static_cast<double> (gcdi);

  /* get the best possible cut */

  cur_cut->min_weak_loss = 0.0;
  for ( j = 0; j < n; j++ ) {
    cur_cut->coef[j] = cur_cut->non_weak_coef[j];
    if ( mod2(cur_cut->coef[j]) == ODD ) 
      cur_cut->min_weak_loss += p_ilp->min_loss_by_weak[j];
  }
  cur_cut->rhs = cur_cut->non_weak_rhs;
  cur_cut->ok = 
    best_cut(cur_cut->coef,&cur_cut->rhs,&cur_cut->violation,TRUE,FALSE);
  cur_cut->one_norm = 0;
  for ( j = 0; j < n; j++ ) {
    if ( cur_cut->coef[j] > 0 ) cur_cut->one_norm += cur_cut->coef[j];
    else cur_cut->one_norm -= cur_cut->coef[j];
  }
}

/* get_current_cut: return a cut data type with the information about
   the current cut of the search procedure */

cut *Cgl012Cut::get_current_cut()
{
  int i, j, nz;
  /*double viol;*/
  cut *cut_ptr;

  cut_ptr = reinterpret_cast<cut *> (calloc(1,sizeof(cut)));
  if ( cut_ptr == NULL ) alloc_error(const_cast<char*>("cut_ptr"));  
  cut_ptr->crhs = cur_cut->rhs;
  cut_ptr->csense = 'L';
  /* count the number of nonzeroes in the cut */
  for ( j = 0, nz = 0; j < n; j++ ) if ( cur_cut->coef[j] != 0 ) nz++;
  cut_ptr->cnzcnt = nz; 
  cut_ptr->cind = reinterpret_cast<int *> (calloc(nz,sizeof(int)));
  if ( cut_ptr->cind == NULL ) alloc_error(const_cast<char*>("cut_ptr->cind"));
  cut_ptr->cval = reinterpret_cast<int *> (calloc(nz,sizeof(int)));
  if ( cut_ptr->cval == NULL ) alloc_error(const_cast<char*>("cut_ptr->cval"));
  nz = 0; /*viol = 0.0;*/
  for ( j = 0; j < n; j++ ) {
    if ( cur_cut->coef[j] != 0 ) {
      cut_ptr->cind[nz] = j;
      cut_ptr->cval[nz] = cur_cut->coef[j];
      nz++;
      /* viol += p_ilp->xstar[j] * (double) cur_cut->coef[j]; */
    }
  }
  /* viol -= (double) cur_cut->rhs; */
  /* cut_ptr->violation = viol; */
  cut_ptr->violation = cur_cut->violation;
  cut_ptr->n_of_constr = 0;
  cut_ptr->constr_list = reinterpret_cast<int *> (calloc(inp_ilp->mr,sizeof(int)));
  if ( cut_ptr->constr_list == NULL ) alloc_error(const_cast<char*>("cut_ptr->constr_list"));
  cut_ptr->in_constr_list = reinterpret_cast<short int *> (calloc(inp_ilp->mr,sizeof(short int)));
  if ( cut_ptr->in_constr_list == NULL ) alloc_error(const_cast<char*>("cut_ptr->in_constr_list"));
  for ( i = 0; i < m; i++ ) {
    if ( cur_cut->in_constr_list[i] == IN ) {
      cut_ptr->in_constr_list[i] = IN;
      cut_ptr->constr_list[cut_ptr->n_of_constr] = i;
      (cut_ptr->n_of_constr)++;
    }
    else cut_ptr->in_constr_list[i] = OUT;
  }
  return(cut_ptr);
}

/* best neighbour: find the cut to be added/deleted from the current
   solution among those allowed by the tabu rules */

short int Cgl012Cut::best_neighbour(cut_list *out_cuts /* list of the violated cuts found */)
{
  int i, ibest;
  short int itype, itypebest=-1;
  double score, max_score;
  cut *new_cut;

  /* cycle through all the constraints in your problem ... */

  max_score = -INF;  
  ibest = NONE;
  for ( i = 0; i < m; i++ ) {
    if ( ! p_ilp->row_to_delete[i] && allowed(i) ) {
      if ( in_cur_cut(i) ) {
	/* constraint i is in the current set of constraints  
	   (those defining the current cut) */
	itype = DEL;
      }
      else {
	/* constraint i is not in the current set of constraints  
	   (those defining the current cut) */
	itype = ADD;
      } 
      score = score_by_moving(i,itype,max_score); 
      if ( score > max_score ) {
	/* best cut found in this iteration: store it */
	ibest = i;
	itypebest = itype;
	max_score = score;
      }
    }
  } /* for ( i = 0; i < m; i++ ) */
  
  if ( ibest == NONE ) {
#ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_TABU
    printf("No move could be performed by best_neighbour\n");
#endif
    return(TRUE);
  }
  modify_current(ibest,itypebest);  
  if ( cur_cut->violation > MIN_VIOLATION + EPS ) {
#ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_TABU
    printf("... adding the current cut to the output list - it = %d viol = %f\n",it, cur_cut->violation);
#endif
    new_cut = get_current_cut();
    out_cuts = add_cut_to_list(new_cut,out_cuts);
    last_it_add = it;
  }
  return(FALSE);
}

/* memory_reaction: perform the long term reaction by cheching whether the
   current solution has already been visited or the best solution has not 
   been updated for too many iterations */

void memory_reaction()
{
  int cycle_length;

  if ( hash_search(&cycle_length) ) {
    if ( cycle_length < 2 * ( m - 1 ) ) {
      increase_prohib_period();
      return;
    }
  }
  else hash_insert();
  if ( it - last_prohib_period_mod > B ) 
    decrease_prohib_period();
}

/* add_tight_constraint: initialize the current cut by adding a tight 
   constraint to it */
       
void Cgl012Cut::add_tight_constraint()
{
  int i, ntight;
  double smin=COIN_DBL_MAX;
  abort();
  int *tight;
    
  ntight = 0;
  tight = reinterpret_cast<int *> (calloc(m,sizeof(int)));
  if ( tight == NULL ) alloc_error(const_cast<char*>("tight"));
  for ( i = 0; i < m; i++ ) {
    /* search for the tightest constraint never added to cut */
    if ( last_moved[i] < 0 && p_ilp->slack[i] < smin ) {
      if ( p_ilp->slack[i] < ZERO ) {
	/* tight constraint */
	smin = ZERO;
	tight[ntight] = i;   
	ntight++;
      }
      else {
	/* best constraint so far */
	smin = p_ilp->slack[i];
	tight[0] = i;
	ntight = 1;
      }
    }
  }
  if ( ntight > 0 ) i = tight[rand() % ntight];
  /* if all constraints have already been in cur_cut choose first at random */
  else i = rand() % m;
  free(tight);
  modify_current(i,ADD);
}    

/* initialize: initialize the data structures for local search */

void Cgl012Cut::initialize()
{
  int i;

  m = inp_ilp->mr;
  n = inp_ilp->mc;
  it = 0;
  last_it_add = 0;
  last_it_restart = 0;
  last_prohib_period_mod = 0;
  prohib_period = IN_PROHIB_PERIOD; 
  initialize_cur_cut();
  last_moved = reinterpret_cast<int *> (calloc(m,sizeof(int)));
  if ( last_moved == NULL ) alloc_error(const_cast<char*>("last_moved"));
  for ( i = 0; i < m; i++ ) {
    last_moved[i] = -COIN_INT_MAX;
  }
  initialize_hash_table();
  add_tight_constraint();
  A = m;
  B = 10 * m;
}

/* restart: perform a restart of the search - IMPORTANT: in the current
   implementation vector last_moved is not cleared at restart */
       
void Cgl012Cut::restart(short int failure /* flag forcing the restart if some trouble occurred */)
{
  if ( failure || ( it - last_it_add > A && it - last_it_restart > A ) ) {
    /* perform restart */
    last_it_restart = it;
    prohib_period = IN_PROHIB_PERIOD;
    last_prohib_period_mod = it;
    clear_hash_table();
    clear_cur_cut();
    add_tight_constraint();
  }
}

/* free_memory: free the memory used by local search */

void free_memory()
{
  free_cur_cut();
  free(last_moved);
  free_hash_table();
}

/* tabu_012: try to identify violated 0-1/2 cuts by a simple tabu search
   procedure adapted from that used by Battiti and Protasi for finding
   large cliques */

cut_list *Cgl012Cut::tabu_012()
{
  short int failure;
  cut_list *out_cuts;

  out_cuts = initialize_cut_list(MAX_CUTS);
  initialize();
 
  it = 0; 
  do {
    memory_reaction();
    failure = best_neighbour(out_cuts);
#ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_TABU
    print_cur_cut();
#endif
    it++;
    restart(failure);

  }
  while ( out_cuts->cnum < MAX_CUTS && it < MAX_TABU_ITER );
  free_memory();
#ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_TABU
    printf("Number of violated cuts found by Tabu Search %d\n",out_cuts->cnum);
    printf("Tabu Search timings: best_neighbour %f  score_by_moving %f coefficient %f  best_cut %f\n",
	   time_best_neigh, time_scor_by_mov, time_coef, time_best_cut);
#endif
  return(out_cuts);
}
//end include "Cgltabu_012.c"
#ifdef POOL 

/* same_pool_cut: check whether two pool cuts are in fact the same cut */

short int same_pool_cut(pool_cut *p_cut1, pool_cut *p_cut2)
{
  int c;

  if ( p_cut1->n_of_constr != p_cut2->n_of_constr ) return(FALSE);
  if ( p_cut1->code != p_cut2->code ) return(FALSE);
  /* assumes constr_list is sorted for increasing/decreasing constraint index */
  for ( c = 1; c < p_cut1->n_of_constr; c++ )
    if ( p_cut1->constr_list[c] != p_cut2->constr_list[c] ) return(FALSE);
  return(TRUE);
}

/* free_pool_cut: free the memory for a non-empty pool cut */

void free_pool_cut(pool_cut *p_cut)
{
  if ( p_cut == NULL ) return;
  if ( p_cut->constr_list != NULL ) free(p_cut->constr_list);
  free(p_cut);
}

/* initialize_pool: initialize the pool data structure */

void initialize_pool()
{
  pool = (pool_cut_list *) calloc(1,sizeof(pool_cut_list));
  if ( pool == NULL ) alloc_error(const_cast<char*>("pool"));
  pool->cnum = 0;
  pool->list = (pool_cut **) calloc(MAX_CUT_POOL,sizeof(pool_cut *));
  if ( pool->list == NULL ) alloc_error(const_cast<char*>("pool->list"));
  pool->ncod = (int *) calloc(MAX_CUT_COD,sizeof(int));
  if ( pool->ncod == NULL ) alloc_error(const_cast<char*>("pool->ncod"));
}

/* free_pool: free the memory used by the pool */

void free_pool()
{
  int c;

  if ( pool == NULL ) return;
  for ( c = 0; c < pool->cnum; c++ ) 
    free_pool_cut(pool->list[c]);
  free(pool);
}

/* clean_pool: remove form the pool the cuts which are inactive since a
   large number of iterations */

void clean_pool()
{
  int c, d;

  /* the pool compression could be implemented more efficiently */

  for ( c = 0; c < pool->cnum; c++ ) {
    if ( pool->list[c]->it_found > MAX_ITER_POOL ) {
      free_pool_cut(pool->list[c]);
      pool->list[c] = NULL;
      for ( d = c ; d < pool->cnum; d++ ) 
	pool->list[d] = pool->list[d + 1];
      (pool->cnum)--;
    }
  }
}

/* insert_cut_in_pool: add a cut to the pool if there is space */

void insert_cut_in_pool(pool_cut *p_cut)
{

  if ( pool->cnum == MAX_CUT_POOL ) {
#ifdef COIN_DEVELOP
    printf("Warning: pool is full and separated cuts cannot be added\n");
#endif
    return;
  }
  pool->list[pool->cnum] = p_cut;
  (pool->cnum)++;
  (pool->ncod[p_cut->code])++;
}

/* cut_is_in_pool: check whether a given cut is already in the pool */

short int cut_is_in_pool(cut *v_cut)
{
  int c, i, cod;
  short int equal;
  short int *flag_v;
  pool_cut *p_cut;

/*
print_cut(v_cut);
printf("checking for a cut in the pool ...\n");
*/
  cod = hash_addr(inp_ilp->mr,v_cut->in_constr_list);
  if ( pool->ncod[cod] == 0 ) return(FALSE);
  /* trivial sequential search */
/*
printf("... sequential search needed ...\n");
*/
  flag_v = (short int *) calloc(inp_ilp->mr,sizeof(short int));
  if ( flag_v == NULL ) alloc_error(const_cast<char*>("flag_v"));
  for ( c = 0; c < pool->cnum; c++ ) {
    if ( pool->list[c]->code != cod ) continue;
    p_cut = pool->list[c];
    equal = TRUE;
    for ( i = 0; i < inp_ilp->mr; i++ ) flag_v[i] = OUT;
    for ( i = 0; i < p_cut->n_of_constr; i++ ) 
      flag_v[p_cut->constr_list[i]] = IN;
    for ( i = 0; i < inp_ilp->mr; i++ ) {
      if ( v_cut->in_constr_list[i] != flag_v[i] ) {
	equal = FALSE;
	break;
      }
    }
    if ( equal ) {
      free(flag_v);
/*
printf("... cut was in the pool!\n");
*/
      return(TRUE);
    }
  }
  return(FALSE);
}

/* good_pool_cut: check whether the current cut is worth adding to the pool */

short int good_pool_cut(cut *v_cut)
{
  /* no check performed */
  return(TRUE);
}

/* add_cuts_to_pool: add the cuts separated to the pool structure */

void add_cuts_to_pool(cut_list *out_cuts)
{
  int i, c;
  cut *v_cut;
  pool_cut *p_cut;
  /* float ti, tf, tti, ttf; */
  /* static float tcutis = 0.0, taddcut = 0.0; */

  /* second_(&tti); */

  if ( pool == NULL ) initialize_pool();
  if ( pool->cnum >= MAX_CUT_POOL * CLEAN_THRESH ) clean_pool();

  for ( i = 0; i < out_cuts->cnum; i++ ) {
    v_cut = out_cuts->list[i];
    /* second_(&ti); */
    if ( cut_is_in_pool(v_cut) ) {
      /* second_(&tf); */
      /* tcutis += tf - ti; */
      continue;
    }
    else {
      /* second_(&tf); */
      /* tcutis += tf - ti; */
    }
    if ( good_pool_cut(v_cut) ) {
      /* add the cut to the pool list */
      p_cut = (pool_cut *) calloc(1,sizeof(pool_cut));
      if ( p_cut == NULL ) alloc_error(const_cast<char*>("p_cut"));
      p_cut->n_of_constr = v_cut->n_of_constr;
      p_cut->constr_list = (int *) calloc(v_cut->n_of_constr,sizeof(int));
      if ( p_cut->constr_list == NULL ) alloc_error(const_cast<char*>("p_cut->constr_list"));
      for ( c = 0; c < v_cut->n_of_constr; c++ ) 
	p_cut->constr_list[c] = v_cut->constr_list[c];
      p_cut->code = hash_addr(inp_ilp->mr,v_cut->in_constr_list);
      p_cut->n_it_violated = 0;
      p_cut->it_found = sep_iter;
      insert_cut_in_pool(p_cut);
    }
  }

  /* second_(&ttf); */
  /* taddcut += ttf - tti; */
  /* printf("add_cuts_to_pool: tcutis %f  taddcut %f\n",tcutis,taddcut); */

}

/* interesting_var: decides whether a variable is relevant in the
   separation or not */

short int interesting_var(int j /* variable to be evaluated */)
{
  /* return ( vlog[j]->n_it_zero < MANY_IT_ZERO ); */
  /* if ( aggr ) return (TRUE); */
  return( ! p_ilp->col_to_delete[j] );
}

/* get_cuts_from_pool: select from the pool a convenient set of violated
   constraints to be added to the current LP */

static double max_score_ever = ZERO; /* maximum score of a violated cut during
					the whole cutting plane procedure */
#define MIN_CUT_SCORE ( max_score_ever / MIN_SCORE_RANGE )
#define MAX_CUT_SCORE ( max_score_ever / MAX_SCORE_RANGE )
#define MIN_IT_VIOL 2

cut_list *get_cuts_from_pool(
short int after_sep /* flag telling whether the pool is searched after
			a new separation in which case only new cuts are
			checked */
			     )
{
  int c, crhs, j, k, l, maxc, n_interest_var;
  double viol, score, maxscore, min_cut_score, max_cut_score;
  short int ok;
  int *interest_var, *ccoef;
  short int *added;
  pool_cut *p_cut;
  select_cut **best_var_cut;
  cut *a_cut;
  cut_list *add_cuts;
  /* float ti, tf, tti, ttf; */
  /* static float tgetcut = 0.0, talloc = 0.0, tgetcoef = 0.0, tbestcut = 0.0, tscore = 0.0, tupdbest = 0.0, tadd = 0.0; */

  /* second_(&tti); */

  if ( pool == NULL ) {
    add_cuts = initialize_cut_list(1);  
    return(add_cuts);
  }

  /* in the current implementation, the cut with best score with nonzero
     coefficient for each "interesting" variable is added to the current LP, 
     provided the cut satisfies some requirements (violation, depth, etc.) */

  /* define the set of the interesting variables */
  interest_var = (int *) calloc(p_ilp->mc,sizeof(int));
  if ( interest_var == NULL ) alloc_error(const_cast<char*>("interest_var"));
  n_interest_var = 0;
  for ( j = 0; j < p_ilp->mc; j++ ) {
    if ( interesting_var(j) ) {
      interest_var[n_interest_var] = j;
      n_interest_var++;
    }
  }
  best_var_cut = (select_cut **) calloc(n_interest_var,sizeof(select_cut *));
  if ( best_var_cut == NULL ) alloc_error(const_cast<char*>("best_var_cut"));
  for ( k = 0; k < n_interest_var; k++ ) {
    best_var_cut[k] = (select_cut *) calloc(1,sizeof(select_cut));
    if ( best_var_cut[k] == NULL ) alloc_error(const_cast<char*>("best_var_cut[k]"));
    best_var_cut[k]->ccoef = (int *) calloc(p_ilp->mc,sizeof(int));
    if ( best_var_cut[k]->ccoef == NULL )
      alloc_error(const_cast<char*>("best_var_cut[k]->ccoef"));
    best_var_cut[k]->score = -INF;
  }

  /* find the cuts with the best scores and the list of the cuts
     violated by the current fractional point */
  maxscore = -INF;
  ccoef = (int *) calloc(p_ilp->mc,sizeof(int));
  if ( ccoef == NULL ) alloc_error(const_cast<char*>("ccoef"));

  /* second_(&tf); */
  /* talloc += tf - tti; */

  min_cut_score = MIN_CUT_SCORE;
  max_cut_score = MAX_CUT_SCORE;

  for ( c = 0; c < pool->cnum; c++ ) {
    p_cut = pool->list[c];
    /* if a new separation was made before the call check only new cuts */
    if ( after_sep && p_cut->it_found != sep_iter ) continue; 
    /* determine the actual coefficients of the cut */
    /* second_(&ti); */
    ok = get_ori_cut_coef(p_cut->n_of_constr,p_cut->constr_list,
			  ccoef,&crhs,TRUE);
    /* second_(&tf); */
    /* tgetcoef += tf - ti; */
    /* second_(&ti); */
    ok = ok && best_cut(ccoef,&crhs,&viol,TRUE,TRUE);
    /* second_(&tf); */
    /* tbestcut += tf - ti; */
    if ( ok && viol > MIN_VIOLATION ) {
      (p_cut->n_it_violated)++; 
      /* second_(&ti); */
      score = cut_score(ccoef,crhs,viol,TRUE);
      if ( score > maxscore ) {
	maxscore = score;
	maxc = c;
      }
      /* second_(&tf); */
      /* tscore += tf - ti; */
      if ( ! aggr ) {
	if ( score < min_cut_score ) continue;
	if ( score < max_cut_score && p_cut->n_it_violated < MIN_IT_VIOL ) 
	  continue;
      }
      /* second_(&ti); */
      for ( k = 0; k < n_interest_var; k++ ) {
	j = interest_var[k];
	if ( ccoef[j] != 0 && best_var_cut[k]->score < score ) {
	  best_var_cut[k]->score = score;
	  for ( l = 0; l < p_ilp->mc; l++ ) 
	    best_var_cut[k]->ccoef[l] = ccoef[l];
	  best_var_cut[k]->crhs = crhs;
	  best_var_cut[k]->pool_index = c;
	}
      }
      /* second_(&tf); */
      /* tupdbest += tf - ti; */
    }
    else p_cut->n_it_violated = 0;
  }
  free(ccoef);

/* printf("maxscore of a cut : %f  ever: %f\n",maxscore,max_score_ever); */
  if ( maxscore > max_score_ever ) max_score_ever = maxscore;

  /* second_(&ti); */

  add_cuts = initialize_cut_list(n_interest_var);
  added = (short int *) calloc(pool->cnum,sizeof(short int));
  if ( added == NULL ) alloc_error(const_cast<char*>("added"));
  for ( c = 0; c < pool->cnum; c++ ) added[c] = FALSE;

  for ( k = 0; k < n_interest_var; k++ ) {
    j = interest_var[k];  
    if ( ! added[best_var_cut[k]->pool_index] &&
	 best_var_cut[k]->score >= MIN_CUT_SCORE ) {
      /* add the cut to the cut list on output */
      a_cut = define_cut(best_var_cut[k]->ccoef,best_var_cut[k]->crhs);
      add_cuts = add_cut_to_list(a_cut,add_cuts);
      added[best_var_cut[k]->pool_index] = TRUE;
    }
    free(best_var_cut[k]->ccoef);
    free(best_var_cut[k]);
  }
  free(best_var_cut);
  free(interest_var);
  free(added);

  /* second_(&tf); */
  /* tadd += tf - ti; */

  /* second_(&ttf); */
  /* tgetcut += ttf - tti; */
  /* printf("get_cuts_from_pool: talloc %f  tgetcoef %f  tbestcut %f  tscore %f  tupdbest %f  tadd %f  tgetcut %f\n",talloc,tgetcoef,tbestcut,tscore,tupdbest,tadd,tgetcut); */

  return(add_cuts);
}
#endif
/* initialize_log_var: initialize the log information for the problem variables */  

void Cgl012Cut::initialize_log_var()
{
  int j;
  if (!vlog) {
    if (p_ilp->mc) {
      vlog = reinterpret_cast<log_var **> (calloc(p_ilp->mc,sizeof(log_var *)));
      if ( vlog == NULL ) alloc_error(const_cast<char*>("vlog"));
      for ( j = 0; j < p_ilp->mc; j++ ) {
	vlog[j] = reinterpret_cast<log_var *> (calloc(1,sizeof(log_var)));
	if ( vlog[j] == NULL ) alloc_error(const_cast<char*>("vlog[j]"));
	vlog[j]->n_it_zero = 0;
      }
    }
  } else {
    // just initialize counts
    for ( j = 0; j < p_ilp->mc; j++ ) {
      vlog[j]->n_it_zero = 0;
    }
  }
}

/* free_log_var */

void Cgl012Cut::free_log_var()
{
  if (vlog) {
    int j;
    for ( j = 0; j < p_ilp->mc; j++ ) free(vlog[j]);
    free(vlog);
    vlog=NULL;
  }
}

/* update_log_var: update the log information for the problem variables */

void Cgl012Cut::update_log_var()
{
  int j;

  /* so far one counts only the number of consecutive iterations with
     0 value for each variable */

  if ( vlog == NULL ) initialize_log_var();
    
  for ( j = 0; j < p_ilp->mc; j++ ) {
    if ( p_ilp->xstar[j] < ZERO && p_ilp->xstar[j] > - ZERO ) 
      vlog[j]->n_it_zero++;
    else  vlog[j]->n_it_zero = 0;
  }
}

/* the final implementation should use the following additional functions:
   init_sep_012_cut: defines the inp_ilp and p_ilp data structures and 
		     initializes pool and vlog
   sep_012_cut: updates inp_ilp and p_ilp, performs separation and pool
		management
   kill_sep_012_cut: frees all the permanent memory (inp_ilp, p_ilp,
		     pool, vlog, etc.)
*/

int Cgl012Cut::sep_012_cut(
/*
  INPUT parameters:
*/
int mr, /* number of rows in the ILP matrix */
int mc, /* number of columns in the ILP matrix */
int mnz, /* number of nonzero's in the ILP matrix */
int *mtbeg, /* starting position of each row in arrays mtind and mtval */
int *mtcnt, /* number of entries of each row in arrays mtind and mtval */
int *mtind, /* column indices of the nonzero entries of the ILP matrix */
int *mtval, /* values of the nonzero entries of the ILP matrix */
int *vlb, /* lower bounds on the variables */
int *vub, /* upper bounds on the variables */
int *mrhs, /* right hand sides of the constraints */
char *msense, /* senses of the constraints: 'L', 'G' or 'E' */
const double *xstar, /* current optimal solution of the LP relaxation */
bool aggressive, /* flag asking whether as many cuts as possible are
			 required on output (TRUE) or not (FALSE) */
/*
  OUTPUT parameters (the memory for the vectors is allocated INTERNALLY
  by the procedure: if some memory is already allocated, it is FREED):
*/
int *cnum, /* number of violated 0-1/2 cuts identified by the procedure */
int *cnzcnt, /* overall number of nonzero's in the cuts */
int **cbeg, /* starting position of each cut in arrays cind and cval */
int **ccnt, /* number of entries of each cut in arrays cind and cval */
int **cind, /* column indices of the nonzero entries of the cuts */
int **cval, /* values of the nonzero entries of the cuts */
int **crhs, /* right hand sides of the cuts */
char **csense /* senses of the cuts: 'L', 'G' or 'E' */
/* 
  NOTE that all the numerical input/output vectors are INTEGER (with
  the exception of xstar), since the procedure is intended to work
  with pure ILP's, and that the ILP matrix has to be given on input
  in ROW format.
*/
		)
{
#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  float tbasi, tbasf;
#endif
  errorNo=0;
  cut_list *out_cuts, *add_cuts;
#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&tti);
#endif

  aggr = aggressive; 
  profileStartSeconds_ = CoinGetTimeOfDay();
  resetTimeCheckState();

  /* load the input ILP into an internal data structure */
  
  //ilp_load(mr,mc,mnz,mtbeg,mtcnt,mtind,mtval,
  //   vlb,vub,mrhs,msense);
  if (!inp_ilp)
      return FALSE;
  inp_ilp->xstar = xstar;
  

  /* construct an internal data structure containing all the information
     which can be useful for 0-1/2 cut separation  - this may in fact be
     done only at the first call of the separation procedure */

#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&ti);
#endif

  get_parity_ilp();

/*
print_double_vect("xstar",p_ilp->xstar,p_ilp->mc); 
print_parity_ilp();
*/

#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&tf);
  prep_time += tf - ti;
#endif

  if ( p_ilp->mnz == 0 ) {
#ifdef COIN_DEVELOP
    printf("Warning: no significant constraint for 0-1/2 cut separation\n");
    printf("... end separation\n");
#endif
    //free_ilp();
    //free_parity_ilp();
#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&ttf);
  total_time += ttf - tti;
#endif 
    return(FALSE);
  }

/* print_double_vect("xstar",p_ilp->xstar,p_ilp->mc); */

  sep_iter++;
  update_log_var();
#ifdef CGL_ZH_ADVANCED_DEBUG_CONSTRAINT_PROFILING
  zhProfileStartRound();
#endif

#ifdef POOL

  /* search for possible violated cuts in the pool */

#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&tpi);
#endif
  add_cuts = get_cuts_from_pool(FALSE);
#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&tpf);
  pool_time += tpf - tpi;
#ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_TIMING
  printf("... time elapsed at the end of get_cuts_from_pool: %f\n",tpf - tti);
#endif
#endif

  if ( add_cuts->cnum == 0 ) {
    free_cut_list(add_cuts);
  }
  else {
/* printf("Violated cuts found in the pool - no separation procedure used\n"); */
    goto free_memory;
  }

#endif

  /* try to identify violated 0-1/2 cuts by using the original procedure 
     described in Caprara and Fischetti's MP paper */

#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&tbasi);
#endif
  
  out_cuts = basic_separation();
#ifdef CGL_ZH_ADVANCED_DEBUG_CONSTRAINT_PROFILING
  zhProfileFlushRound(sep_iter, out_cuts ? out_cuts->cnum : 0);
#endif

#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&tbasf);
  tot_basic_sep_time += tbasf - tbasi;
  avg_basic_sep_time = tot_basic_sep_time / (float) sep_iter;
#endif

//#define TABU_SEARCH 
#ifdef TABU_SEARCH

  /* try to identify violated cuts by tabu search if none was found */

  if ( out_cuts->cnum == 0 ) {
    free_cut_list(out_cuts); 

#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&ttabi);
#endif

    out_cuts = tabu_012(); 

#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&ttabf);
  tabu_time += ttabf - ttabi;
#endif

  }

#endif

#ifdef POOL

  /* add the cuts separated to the pool */

#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&tpi);
#endif
  add_cuts_to_pool(out_cuts);
  free_cut_list(out_cuts); 

  /* select from the pool a convenient set of violated constraints
     to be added to the current LP */

  add_cuts = get_cuts_from_pool(TRUE);
#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&tpf);
  pool_time += tpf - tpi;
#ifdef CGL_ZH_ADVANCED_DEBUG_PRINT_TIMING
  printf("... time elapsed at the end of get_cuts_from_pool: %f\n",tpf - tti);
#endif
#endif

#else

  /* give on output the cuts separated */

  add_cuts = out_cuts;

#endif

  //free_ilp();
  //free_parity_ilp();
  
#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
  second_(&ttf);
  total_time += ttf - tti;
#endif 

  if ( add_cuts->cnum > 0 ) {
    getcuts(add_cuts,cnum,cnzcnt,cbeg,ccnt,cind,cval,crhs,csense);
/* print_cut_list(add_cuts); */
    free_cut_list(add_cuts); 
    return(TRUE); 
  }
  else {
    free_cut_list(add_cuts); 
    return(FALSE); 
  }
}

#ifdef CGL_ZH_ADVANCED_DEBUG_TIMING
/* print_times: print the timings of the separation procedure */

void print_times()
{
  printf("... separation timings \n");
  printf("times  total: %f  prep: %f  weak: %f  aux: %f  path: %f cycle: %f  cut: %f (%d calls)  bw: %f  coef: %f  pool: %f  tabu: %f\n",
    total_time, prep_time, weak_time, aux_time, path_time, cycle_time, cut_time, cut_ncalls, bw_time, coef_time, pool_time, tabu_time);
}
#endif
//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
Cgl012Cut::Cgl012Cut () :
  inp_ilp(NULL),
  p_ilp(NULL),
  iter(0),
  gap(0.0),
  maxgap(0.0),
  errorNo(0),
  sep_iter(0),
  vlog(NULL),
  aggr(true),
  maxSeconds_(0.0),
  profileStartSeconds_(0.0),
  timeLimitReached_(false),
  timeCheckCountdown_(1),
  sepGraphSparseThreshold_(MAX_SEP_GRAPH_ACTIVE_NODES),
  rowMaxPairCount_(-1),
  rowMaxFractionalCount_(-1)
{
  // nothing to do here
}
//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
Cgl012Cut::Cgl012Cut (const Cgl012Cut & rhs) :
  inp_ilp(NULL),
  p_ilp(NULL),
  iter(rhs.iter),
  gap(rhs.gap),
  maxgap(rhs.maxgap),
  errorNo(rhs.errorNo),
  sep_iter(rhs.sep_iter),
  vlog(NULL),
  aggr(rhs.aggr),
  maxSeconds_(rhs.maxSeconds_),
  profileStartSeconds_(0.0),
  timeLimitReached_(false),
  timeCheckCountdown_(1),
  sepGraphSparseThreshold_(rhs.sepGraphSparseThreshold_),
  rowMaxPairCount_(rhs.rowMaxPairCount_),
  rowMaxFractionalCount_(rhs.rowMaxFractionalCount_)
{
  if (rhs.p_ilp||rhs.vlog||inp_ilp)
    abort();  
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
Cgl012Cut::~Cgl012Cut ()
{
#ifdef CGL_ZH_ADVANCED_DEBUG_CONSTRAINT_PROFILING
  cglZeroHalfConstraintProfileStates().erase(this);
#endif
  free_log_var();
  free_parity_ilp();
  free_ilp();
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
Cgl012Cut &
Cgl012Cut::operator=(
                   const Cgl012Cut& rhs)
{
  if (this != &rhs) {
#ifdef CGL_ZH_ADVANCED_DEBUG_CONSTRAINT_PROFILING
    cglZeroHalfConstraintProfileStates().erase(this);
#endif
    if (rhs.p_ilp||rhs.vlog||inp_ilp)
      abort();  
    free_log_var();
    free_parity_ilp();
    free_ilp();
#if 0
    inp_ilp = reinterpret_cast<ilp *> (calloc(1,sizeof(ilp)));
    if ( inp_ilp == NULL ) alloc_error(const_cast<char*>("inp_ilp"));

    inp_ilp->mr = rhs.inp_ilp->mr; inp_ilp->mc = rhs.inp_ilp->mc;
    inp_ilp->mnz = rhs.inp_ilp->mnz; 
    inp_ilp->mtbeg = rhs.inp_ilp->mtbeg; inp_ilp->mtcnt = rhs.inp_ilp->mtcnt; 
    inp_ilp->mtind = rhs.inp_ilp->mtind; inp_ilp->mtval = rhs.inp_ilp->mtval; 
    inp_ilp->vlb = rhs.inp_ilp->vlb; inp_ilp->vub = rhs.inp_ilp->vub; 
    inp_ilp->mrhs = rhs.inp_ilp->mrhs; inp_ilp->msense = rhs.inp_ilp->msense;
#endif
    iter = rhs.iter;
    gap = rhs.gap;
    maxgap = rhs.maxgap;
    errorNo = rhs.errorNo;
    sep_iter = rhs.sep_iter;
    aggr = rhs.aggr;
    maxSeconds_ = rhs.maxSeconds_;
    profileStartSeconds_ = 0.0;
    resetTimeCheckState();
    sepGraphSparseThreshold_ = rhs.sepGraphSparseThreshold_;
    rowMaxPairCount_ = rhs.rowMaxPairCount_;
    rowMaxFractionalCount_ = rhs.rowMaxFractionalCount_;
  }
  return *this;
}

void Cgl012Cut::setMaxSeconds(double value)
{
  maxSeconds_ = value;
}

double Cgl012Cut::getMaxSeconds() const
{
  return maxSeconds_;
}

void Cgl012Cut::setSepGraphSparseThreshold(int value)
{
  sepGraphSparseThreshold_ = value;
}

int Cgl012Cut::getSepGraphSparseThreshold() const
{
  return sepGraphSparseThreshold_;
}

void Cgl012Cut::setRowMaxPairCount(int value)
{
  rowMaxPairCount_ = value;
}

int Cgl012Cut::getRowMaxPairCount() const
{
  return rowMaxPairCount_;
}

void Cgl012Cut::setRowMaxFractionalCount(int value)
{
  rowMaxFractionalCount_ = value;
}

int Cgl012Cut::getRowMaxFractionalCount() const
{
  return rowMaxFractionalCount_;
}
