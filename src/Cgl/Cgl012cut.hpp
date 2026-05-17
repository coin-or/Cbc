// Copyright (C) 2010, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).
/** @file 012cut.h Include file for C coded 0-1/2 separator */
#ifndef CGL012CUT
#define CGL012CUT
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "CglConfig.h"

struct edge_st;
typedef struct edge_st edge;

#define CGL_NEW_SHORT
#ifndef CGL_NEW_SHORT
typedef  /* arc */
   struct arc_st
{
   int              len;            /* length of the arc */
   struct node_st   *head;           /* head node */
   edge             *backEdge;       /* originating separation-graph edge */
}
  arc;

typedef  /* node */
   struct node_st
{
   arc              *first;           /* first outgoing arc */
   int              dist;	      /* tentative shortest path length */
   struct node_st   *parent;          /* parent pointer */
   struct node_st   *next;            /* next node in queue */
   struct node_st   *prev;            /* previous node in queue */
   int               status;          /* status of node */
   int               temp;            /* for temporary labels */
   int               index;           /* index of the node in the graph */
} node;
#endif
typedef struct 
{
  int length; // Length of arc
  int to; // To node
  edge * backEdge; // Originating separation-graph edge
} cgl_arc;

typedef struct
{
  cgl_arc * firstArc; // First outgoing arc 
  int parentNode; // Parent node in shortest path 
  int index; // Which node I am
  int distanceBack; // Distance back to source
} cgl_node;

typedef struct 
{
  int nnodes; // Number of nodes in graph
  int narcs; // Number of arcs in graph
  cgl_node * nodes;
  cgl_arc * arcs;
} cgl_graph;
/* #define PRINT */
/* #define CGL_ZH_ADVANCED_DEBUG_PRINT_CUTS */
#define REDUCTION

typedef struct {
int mr; /* number of rows in the ILP matrix */
int mc; /* number of columns in the ILP matrix */
int mnz; /* number of nonzero's in the ILP matrix */
int *mtbeg; /* starting position of each row in arrays mtind and mtval */
int *mtcnt; /* number of entries of each row in arrays mtind and mtval */
int *mtind; /* column indices of the nonzero entries of the ILP matrix */
int *mtval; /* values of the nonzero entries of the ILP matrix */
int *vlb; /* lower bounds on the variables */
int *vub; /* upper bounds on the variables */
int *mrhs; /* right hand sides of the constraints */
char *msense; /* senses of the constraints: 'L', 'G' or 'E' */
const double *xstar; /* current optimal solution of the LP relaxation */
} ilp; 

typedef struct {
int mr; /* number of rows in the parity ILP matrix */
int mc; /* number of columns in the parity ILP matrix */
int mnz; /* number of 1's in the parity ILP matrix */
int *mtbeg; /* starting position of each row in arrays mtind and mtval */
int *mtcnt; /* number of entries of each row in arrays mtind and mtval */
int *mtind; /* column indices of the 1's of the parity ILP matrix */
short int *mrhs; /* right hand side parity of the constraints */
 const double *xstar; /* current optimal solution of the LP relaxation */
double *slack; /* slack of the constraints w.r.t. xstar */
short int *row_to_delete; /* flag for marking rows not to be considered */
short int *col_to_delete; /* flag for marking columns not to be considered */
int *gcd; /* greatest common divisor of each row in the input ILP matrix */
short int *possible_weak; /* possible weakening types of each column */
short int *type_even_weak; /* type of even weakening of each column 
                              (lower or upper bound weakening) */
short int *type_odd_weak; /* type of odd weakening of each column 
                              (lower or upper bound weakening) */
double *loss_even_weak; /* loss for the even weakening of each column */
double *loss_odd_weak; /* loss for the odd weakening of each column */
double *min_loss_by_weak; /* minimum loss for the weakening of each column */
} parity_ilp; 

typedef struct {
int nweak; /* number of variables weakened */
int *var; /* list of variables weakened */
short int *type; /* type of weakening (lower or upper bound weakening) */
} info_weak;

typedef struct edge_st {
  int endpoint1, endpoint2; /* endpoints of the edge */
  double weight; /* edge weight */
  short int parity; /* edge parity (even or odd) */
  int constr; /* constraint associated with the edge */
  info_weak *weak; /* weakening information */
} edge;

typedef struct {
  int nnodes; /* number of nodes */
  int nedges; /* number of edges */
  int *nodes; /* indexes of the ILP columns corresponding to the nodes */
  int *ind; /* indexes of the nodes corresponding to the ILP columns */
  bool sparseMode; /* storage mode for separation graph */
  edge **even_adj_list; /* pointers to the even edges */ 
  edge **odd_adj_list; /* pointers to the odd edges */ 
  std::unordered_map<std::uint64_t, edge *> *sparseEdges; /* sparse edge lookup */
  std::vector<edge *> *sparseAdj; /* sparse incident-edge lists per node */
} separation_graph;

#ifndef CGL_NEW_SHORT
typedef struct {
int nnodes; /* number of nodes */
int narcs; /* number of arcs */
node *nodes; /* array of the nodes - see "types_db.h" */ 
arc *arcs; /* array of the arcs - see "types_db.h" */ 
} auxiliary_graph;
#else
typedef struct {
int nnodes; /* number of nodes */
int narcs; /* number of arcs */
cgl_node *nodes; /* array of the nodes - see "types_db.h" */ 
cgl_arc *arcs; /* array of the arcs - see "types_db.h" */ 
} auxiliary_graph;
#endif

typedef struct {
long dist; /* distance from/to root */
int pred; /* index of the predecessor */
} short_path_node;

typedef struct {
double weight; /* overall weight of the cycle */
int length; /* number of edges in the cycle */
edge **edge_list; /* list of edges in the cycle */
} cycle;

typedef struct {
int cnum; /* overall number of cycles */
cycle **list; /* pointers to the cycles in the list */
} cycle_list;

typedef struct {
int n_of_constr; /* number of constraints combined to get the cut */
int *constr_list; /* list of the constraints combined */
short int *in_constr_list; /* flag saying whether a given constraint is
                              in the list of constraints of the cut (IN)
                              or not (OUT) */
int cnzcnt; /* overall number of nonzero's in the cut */
int *cind; /* column indices of the nonzero entries of the cut */
int *cval; /* values of the nonzero entries of the cut */
int crhs; /* right hand side of the cut */
char csense; /* sense of the cut: 'L', 'G' or 'E' */
double violation; /* violation of the cut w.r.t. the current LP solution */
} cut;

typedef struct {
int cnum; /* overall number of cuts */
cut **list; /* pointers to the cuts in the list */
} cut_list;

typedef struct {
int n_of_constr; /* number of constraints combined to get the cut */
int *constr_list; /* list of the constraints combined */
int code; /* identifier of the cut */
int n_it_violated; /* number of consecutive iterations (starting from the
                      last and going backward) in which the cut was
                      violated by the LP solution */
int it_found; /* iteration in which the cut was separated */
double score; /* score of the cut, used to choose wich cut should be 
                 added to the current LP (if any) */
} pool_cut;
 
typedef struct {
int cnum; /* overall number of cuts */
pool_cut **list; /* pointers to the cuts in the list */
int *ncod; /* number of cuts with a given code in the pool */
} pool_cut_list;

typedef struct {
int *ccoef; /* coefficients of the cut */
int crhs; /* right hand side of the cut */
int pool_index; /* index of the cut in the pool */
double score; /* cut score (to be maximized) */
} select_cut;

typedef struct {
int n_it_zero; /* number of consecutive iterations (starting from the
                   last and going backward) in which each variable took
                   the value 0 in the LP solution */
} log_var;
/** 012Cut Generator Class

 This class is to make Cgl01cut thread safe etc
*/

class CGLLIB_EXPORT Cgl012Cut {
 
public:

  /**@name Generate Cuts */
  //@{
int sep_012_cut(
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
		);
void ilp_load(
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
	      );
void free_ilp();
/* alloc_parity_ilp: allocate the memory for the parity ILP data structure */

void alloc_parity_ilp(
		      int mr, /* number of rows in the ILP matrix */
		      int mc, /* number of columns in the ILP matrix */
		      int mnz /* number of nonzero's in the ILP matrix */
		      );
void free_parity_ilp();
  void initialize_log_var();
/* free_log_var */
  void free_log_var();
private:
  bool checkTimeLimit(const char *phase, const char *detail = NULL);
  bool timeLimitReached() const;
  void resetTimeCheckState();
  void ensureBestWeakeningBufferCapacity(int requiredSize);
  void ensureVarsToWeakBufferCapacity(int requiredSize);
#ifdef CGL_ZH_ADVANCED_DEBUG_CONSTRAINT_PROFILING
  void zhProfileStartRound();
  void zhProfileRecordRow(int row, double elapsedSeconds, int pairCount);
  void zhProfileRecordCut(const cut *v_cut);
  void zhProfileFlushRound(int round, int generatedCuts);
#endif
/* best_weakening: find the best upper/lower bound weakening of a set
   of variables */

int best_weakening(
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
		   );

/* best_cut: find the coefficients, the rhs and the violation of the
   best possible cut that can be obtained by weakening a given set of
   coefficients to even and a rhs to odd, dividing by 2 and rounding */

short int best_cut(
		   int *ccoef, /* vector of the coefficients */
		   int *crhs, /* pointer to rhs value */
		   double *violation, /* violation of the cut */
		   short int update, /* TRUE/FALSE: if TRUE, the new ccoef and crhs are 
					given on output */ 
		   short int only_viol /* flag which tells whether only an inequality of
			slack smaller than MAX_SLACK is of interest (TRUE)
			otherwise (FALSE) */
			      );
/* get_cut: extract a hopefully violated cut from an odd cycle of the
   separation graph */

cut *get_cut(
	     cycle *s_cyc /* shortest odd cycles identified in the separation graph */
	     );

/* update_log_var: update the log information for the problem variables */
  void update_log_var();

/* basic_separation: try to identify violated 0-1/2 cuts by using the 
   original procedure described in Caprara and Fischetti's MP paper */

  cut_list *basic_separation();

/* score_by_moving: compute the score of the best cut obtainable from 
   the current local search solution by inserting/deleting a constraint */

double score_by_moving(
		       int i, /* constraint to be moved */
		       short int itype, /* type of move - ADD or DEL */
		       double thresh /* minimum value of an interesting score */
		       );
/* modify_current: update the current local search solution by inserting/
   deleting a constraint */

void modify_current(
		    int i, /* constraint to be moved */
		    short int itype /* type of move - ADD or DEL */
		    );

/* best neighbour: find the cut to be added/deleted from the current
   solution among those allowed by the tabu rules */

  short int best_neighbour(cut_list *out_cuts /* list of the violated cuts found */);

/* add_tight_constraint: initialize the current cut by adding a tight 
   constraint to it */
       
  void add_tight_constraint();

/* tabu_012: try to identify violated 0-1/2 cuts by a simple tabu search
   procedure adapted from that used by Battiti and Protasi for finding
   large cliques */

  cut_list *tabu_012();
/* initialize: initialize the data structures for local search */

  void initialize();
/* restart: perform a restart of the search - IMPORTANT: in the current
   implementation vector last_moved is not cleared at restart */
       
  void restart(short int failure /* flag forcing the restart if some trouble occurred */);
  void print_constr(int i /* constraint to be printed */);
  void print_parity_ilp();

/* get_parity_ilp: construct an internal data structure containing all the 
   information which can be useful for  0-1/2 cut separation */

  void get_parity_ilp();
/* initialize_sep_graph: allocate and initialize the data structure
   to contain the information associated with a separation graph */
 
  separation_graph *initialize_sep_graph();
  cycle_list *get_shortest_odd_cycle_list(
    int j,
    separation_graph *s_graph,
    auxiliary_graph *a_graph);
  void print_cut(cut *v_cut);
/* get_ori_cut_coef: get the coefficients of a cut, before dividing by 2 and
   rounding, starting from the list of the constraints combined to get 
   the cut */

short int get_ori_cut_coef(
			   int n_of_constr, /* number of constraints combined */
			   int *constr_list, /* list of the constraints combined */
			   int *ccoef, /* cut left hand side coefficients */
			   int *crhs, /* cut right hand side */
			   short int only_viol /* flag which tells whether only an inequality of
			slack smaller than MAX_SLACK is of interest (TRUE)
			otherwise (FALSE) */
			   );
/* define_cut: construct a cut data structure from a vector of
   coefficients and a right-hand-side */

cut *define_cut(
		int *ccoef, /* coefficients of the cut */
		int crhs /* right hand side of the cut */
		);

/* cut_score: define the score of a (violated) cut */

double cut_score(
		 int *ccoef, /* cut left hand side coefficients */
		 int crhs, /* cut right hand side */
		 double viol, /* cut violation */
		 short int only_viol /* flag which tells whether only an inequality of
			slack smaller than MAX_SLACK is of interest (TRUE)
			otherwise (FALSE) */
		 );
/* get_current_cut: return a cut data type with the information about
   the current cut of the search procedure */

  cut *get_current_cut();
/* print_cur_cut: display cur_cut on output */

  void print_cur_cut();
  void print_cut_list(cut_list *cuts);
  //@}
public:
  void setMaxSeconds(double value);
  double getMaxSeconds() const;
  void setSepGraphSparseThreshold(int value);
  int getSepGraphSparseThreshold() const;
  void setRowMaxPairCount(int value);
  int getRowMaxPairCount() const;
  void setRowMaxFractionalCount(int value);
  int getRowMaxFractionalCount() const;
  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  Cgl012Cut ();
 
  /// Copy constructor 
  Cgl012Cut (
    const Cgl012Cut &);

  /// Assignment operator 
  Cgl012Cut &
    operator=(
    const Cgl012Cut& rhs);
  
  /// Destructor 
  virtual ~Cgl012Cut ();
  //@}

private:
  
  // Private member methods
   
  /**@name Private methods */
  //@{
  //@}
  
  
  /**@name Private member data */
  //@{

ilp *inp_ilp; /* input ILP data structure */
parity_ilp *p_ilp; /* parity ILP data structure */
int iter;
double gap;
double maxgap;
int errorNo;
  int sep_iter; /* number of the current separation iteration */
  log_var **vlog; /* information about the value attained
				  by the variables in the last iterations,
				  used to possibly set to 0 some coefficient
				  > 0 in a cut to be added */ 
  bool aggr; /* flag saying whether as many cuts as possible are required
		   from the separation procedure (TRUE) or not (FALSE) */
  double maxSeconds_; /* best-effort wall-clock budget for one separation call */
  double profileStartSeconds_;
  bool timeLimitReached_;
  int timeCheckCountdown_;
  std::vector<short int> typeEvenWeakBuffer_;
  std::vector<short int> switchEvenWeakBuffer_;
  std::vector<short int> typeOddWeakBuffer_;
  std::vector<short int> switchOddWeakBuffer_;
  std::vector<int> varsToWeakBuffer_;
  int sepGraphSparseThreshold_; /* active-node threshold for sparse separation graph */
  int rowMaxPairCount_; /* skip rows whose pair count exceeds this threshold; negative disables */
  int rowMaxFractionalCount_; /* skip rows whose fractional count exceeds this threshold; negative disables */
  //@}
};
#endif
