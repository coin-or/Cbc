// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdio>
#include <cassert>

#include "CoinHelperFunctions.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "CglClique.hpp"

/* to prevent the creation of very
 * large incidence matrixes */
#ifndef MAX_CGLCLIQUE_COLS
#define MAX_CGLCLIQUE_COLS 10000
#endif

/*****************************************************************************/

CglClique::CglClique(bool setPacking, bool justOriginalRows) :
  CglCutGenerator(),
   setPacking_(setPacking),
  justOriginalRows_(justOriginalRows),
   sp_numrows(0),
   sp_orig_row_ind(0),
   sp_numcols(0),
   sp_orig_col_ind(0),
   sp_colsol(0),
   sp_col_start(0),
   sp_col_ind(0),
   sp_row_start(0),
   sp_row_ind(0),
   node_node(0),
   petol(-1.0),
   maxNumber_(5000),
   do_row_clique(true),
   do_star_clique(true),
   scl_next_node_rule(SCL_MAX_XJ_MAX_DEG),
   scl_candidate_length_threshold(12),
   scl_report_result(true),
   rcl_candidate_length_threshold(12),
   rcl_report_result(true),
   cl_perm_indices(0),
   cl_perm_length(0),
   cl_indices(0),
   cl_length(0),
   cl_del_indices(0),
   cl_del_length(0)
{}
// Copy constructor
CglClique::CglClique(const CglClique& rhs)
  : CglCutGenerator(rhs),
    setPacking_(rhs.setPacking_),
    justOriginalRows_(rhs.justOriginalRows_),
    sp_numrows(rhs.sp_numrows),
    sp_orig_row_ind(rhs.sp_orig_row_ind),
    sp_numcols(rhs.sp_numcols),
    sp_orig_col_ind(rhs.sp_orig_col_ind),
    sp_colsol(rhs.sp_colsol),
    sp_col_start(rhs.sp_col_start),
    sp_col_ind(rhs.sp_col_ind),
    sp_row_start(rhs.sp_row_start),
    sp_row_ind(rhs.sp_row_ind),
    node_node(rhs.node_node),
    petol(rhs.petol),
    maxNumber_(rhs.maxNumber_),
    do_row_clique(rhs.do_row_clique),
    do_star_clique(rhs.do_star_clique),
    scl_next_node_rule(rhs.scl_next_node_rule),
    scl_candidate_length_threshold(rhs.scl_candidate_length_threshold),
    scl_report_result(rhs.scl_report_result),
    rcl_candidate_length_threshold(rhs.rcl_candidate_length_threshold),
    rcl_report_result(rhs.rcl_report_result),
    cl_perm_indices(rhs.cl_perm_indices),
    cl_perm_length(rhs.cl_perm_length),
    cl_indices(rhs.cl_indices),
    cl_length(rhs.cl_length),
    cl_del_indices(rhs.cl_del_indices),
    cl_del_length(rhs.cl_del_length)
{
}

/*****************************************************************************/

/*****************************************************************************/

void
CglClique::generateCuts(const OsiSolverInterface& si, OsiCuts & cs,
			const CglTreeInfo info)
{
   int i;
   bool has_petol_set = petol != -1.0;
   si.getColType(true); // update integer information
   if (! has_petol_set)
      si.getDblParam(OsiPrimalTolerance, petol);
   int numberOriginalRows = si.getNumRows();
   if (info.inTree&&justOriginalRows_)
     numberOriginalRows = info.formulation_rows;
   int numberRowCutsBefore = cs.sizeRowCuts();
   // First select which rows/columns we are interested in.
   if (!setPacking_) {
      selectFractionalBinaries(si);
      if (!sp_orig_row_ind) {
	 selectRowCliques(si,numberOriginalRows);
      }
   } else {
      selectFractionals(si);
      delete[] sp_orig_row_ind;
      sp_numrows = numberOriginalRows;
      //sp_numcols = si.getNumCols();
      sp_orig_row_ind = new int[sp_numrows];
      for (i = 0; i < sp_numrows; ++i)
	 sp_orig_row_ind[i] = i;
   }
   // Just original rows
   if (justOriginalRows_&&info.inTree) 
     sp_numrows = std::min(info.formulation_rows,sp_numrows);
   
#ifndef MAX_CGLCLIQUE_ROWS
#define MAX_CGLCLIQUE_ROWS 100000
#endif
   if (sp_numrows > MAX_CGLCLIQUE_ROWS || sp_numcols < 2 || sp_numcols>MAX_CGLCLIQUE_COLS) {
     //printf("sp_numrows is %d\n",sp_numrows);
     deleteSetPackingSubMatrix();
     return; // too many rows or too few columns!
   }

   createSetPackingSubMatrix(si);
   fgraph.edgenum = createNodeNode();
   createFractionalGraph();

   if (sp_numcols>1) {
     cl_indices = new int[sp_numcols];
     cl_del_indices = new int[sp_numcols];
     
     if (do_row_clique)
       find_rcl(cs);
     if (do_star_clique)
       find_scl(cs);
     if (!info.inTree&&((info.options&4)==4||((info.options&8)&&!info.pass))) {
       int numberRowCutsAfter = cs.sizeRowCuts();
       for (int i=numberRowCutsBefore;i<numberRowCutsAfter;i++)
	 cs.rowCutPtr(i)->setGloballyValid();
     }
   }

   delete[] cl_indices;     cl_indices = 0;
   delete[] cl_del_indices; cl_del_indices = 0;

   deleteFractionalGraph();
   delete[] node_node;      node_node = 0;
   deleteSetPackingSubMatrix();

   if (! has_petol_set)
      petol = -1;
}

/*****************************************************************************/
/*===========================================================================*
 * Find violated row cliques.
 *
 * Algorithm: For each row of the matrix collect all variables not in the
 * row's support that are non-orthogonal to all variables in the row. These
 * variables form a candidate set in which we look for maximal cliques
 * (enumerate all or try to find one greedily, depending on the size of the
 * candidate set). When a violated maximal clique is found, it is recorded in
 * the cut set.
 *===========================================================================*/

void
CglClique::find_rcl(OsiCuts& cs)
{
   const int nodenum = fgraph.nodenum;
   const fnode *nodes = fgraph.nodes;

   /* A flag for each column that might be used to extend the current row
      clique */
   bool *cand = new bool[nodenum];
   /* In cl_indices we'll list the indices of the 'true' entries in cand */
   /* The degree of each candidate (those listed in cl_indices) */
   int *degrees = new int[nodenum];

   /** An array used in the recursive complete enumeration of maximal cliques.
       The first cl_length entries are used. */
   bool* label = new bool[nodenum];

   int i, j, k;

   /* initialize global variables */
   cl_del_length = 0;
   cl_length = 0;

   int clique_count = 0;
   int largest_length = 0;

   /* for each row of the matrix */
   for (j = 0; j < sp_numrows; j++) {

      /* if the row is of zero length, take the next row */
      const int len = sp_row_start[j+1] - sp_row_start[j];
      if (!len)
	 continue;

      /* the beginning of the row to be considered */
      const int *row = sp_row_ind + sp_row_start[j];

      /* copy the row of node_node corresponding to the first column in 'row'
	 into cand, and take the AND of this vector with every row of
	 node_node corresponding to the rest of the columns in 'row' to
	 determine those columns that are non-orthog to every column in row */
      std::copy(node_node + row[0]*nodenum, node_node + (row[0]+1)*nodenum,
		cand);
      for (i = 1; i < len; i++) {
	 const bool* node_node_col = node_node + row[i] * nodenum;
	 for (k = 0; k < nodenum; k++)
	    cand[k] &= node_node_col[k];
      }
      cl_length = 0;
      for (k = 0; k < nodenum; k++)
	 if (cand[k])
	    cl_indices[cl_length++] = k;
      largest_length = std::max(cl_length, largest_length);

      /* if there is anything in indices, enumerate (or greedily find)
	 maximal cliques */
      if (cl_length > 0) {
	 cl_perm_length = len;
	 cl_perm_indices = row;
	 if (cl_length <= rcl_candidate_length_threshold) {
	    for (i = 0; i < cl_length; i++)
	       label[i] = false;
	    int pos = 0;
	    clique_count += enumerate_maximal_cliques(pos, label, cs);
	 } else {
	    /* order cl_indices into decreasing order of their degrees */
	    for (i = 0; i < cl_length; i++)
	       degrees[i] = nodes[cl_indices[i]].degree;
	    CoinSort_2(degrees, degrees + cl_length, cl_indices,
		       CoinFirstGreater_2<int,int>());
	    clique_count += greedy_maximal_clique(cs);
	 }
      }
   }
   if (rcl_report_result) {
      printf("\nrcl Found %i new violated cliques with the row-clique method",
	     clique_count);
      printf("\nrcl The largest admissible number was %i (threshold %i)\n",
	     largest_length, rcl_candidate_length_threshold);
      if (largest_length < rcl_candidate_length_threshold)
	 printf("rcl    all row cliques have been enumerated\n");
      else
	 printf("rcl    not all row cliques have been eliminated\n");
   }

   delete[] degrees;
   delete[] cand;
   delete[] label;
}

/*****************************************************************************/

/*===========================================================================*
 * Find violated star cliques, a la Hoffman-Padberg.
 *
 * Algorithm: Take min degree node. Check for violated cuts in the subgraph
 * consisting of this node and its neighbors (the "star" of this node). Then
 * delete the node and continue with the now min degree node.
 * 
 * Implementation: Two arrays are defined, one contains the indices the other
 * the degrees of all the nodes still in the graph. If the min degree is 0
 * or 1 then the min degree node can be deleted at once.
 * All cliques are enumerated in v U star(v) if the min degree is smaller
 * than the threshold  scl_candidate_length_threshold, otherwise attemp to
 * find maximal clique greedily.
 * 
 * Note: Indices in current_indices are always kept in increasing order.
 *===========================================================================*/

void
CglClique::find_scl(OsiCuts& cs)
{
   const int nodenum = fgraph.nodenum;
   const fnode *nodes = fgraph.nodes;

   // Return at once if no nodes - otherwise we get invalid reads
   assert (nodenum>1);
   //if (!nodenum)
   //return;
   int *current_indices = new int[nodenum];
   int *current_degrees = new int[nodenum];
   double *current_values = new double[nodenum];

   int *star = cl_indices;
   int *star_deg = new int[nodenum];

   /** An array used in the recursive complete enumeration of maximal cliques.
       The first cl_length entries are used. */
   bool* label = new bool[nodenum];

   int i, cnt1 = 0, cnt2 = 0, cnt3 = 0;
   int clique_cnt_e = 0, clique_cnt_g = 0;
   int largest_star_size = 0;

   /* initialize global variables */
   cl_del_length = 0;

   /* initialize current_nodes, current_degrees and current_values */
   int current_nodenum = nodenum;
   for (i = 0; i < nodenum; i++) {
      current_indices[i] = i;
      current_degrees[i] = nodes[i].degree;
      current_values[i] = nodes[i].val;
   }
   
   /* find first node to be checked */
   int best_ind = scl_choose_next_node(current_nodenum, current_indices,
				       current_degrees, current_values);

   int v = current_indices[best_ind];
   int v_deg = current_degrees[best_ind];
   double v_val = current_values[best_ind];
      
   /* while there are nodes left in the graph ... (more precisely, while
      there are at least 3 nodes in the graph) */
   while (current_nodenum > 2) {

      /* if the best node is of degree < 2 then it can be deleted */
      if (v_deg < 2) {
	 cl_del_indices[cl_del_length++] = v;
	 scl_delete_node(best_ind, current_nodenum,
			 current_indices, current_degrees, current_values);
	 best_ind = scl_choose_next_node(current_nodenum, current_indices,
					 current_degrees, current_values);
	 v = current_indices[best_ind];
	 v_deg = current_degrees[best_ind];
	 v_val = current_values[best_ind];
	 largest_star_size = std::max(largest_star_size, v_deg);
	 continue;
      }

      /* star will contain the indices of v's neighbors (but not v's index) */
      const bool* node_node_start = node_node + nodenum * v;
      int& star_length = cl_length;
      star_length = 0;
      double star_val = v_val;
      for (i = 0; i < current_nodenum; i++) {
	 const int other_node = current_indices[i];
	 if (node_node_start[other_node]) {
	    star[star_length] = other_node;
	    star_deg[star_length++] = current_degrees[i];
	    star_val += current_values[i];
	 }
      }

      /* quick check: if sum of values for the star does not exceed 1 then
	 there won't be a violated clique in the star */
      if (star_val >= 1 + petol) {
	 /* node whose star we're evaluating is always in */
	 cl_perm_length = 1;
	 cl_perm_indices = &v;

	 /* find maximal violated cliques in star. cliques found here might not
	    be maximal wrt to entire fractional graph, only for the current
	    subset of it (some nodes might be already deleted...)
	    Note that star is the same as cl_indices and start_length is
	    cl_length... 
	 */
	 if (v_deg < scl_candidate_length_threshold) { // par
	    /* enumerate if v_deg is small enough */
	    for (i = 0; i < star_length; i++)
	       label[i] = false;
	    int pos = 0;
	    clique_cnt_e += enumerate_maximal_cliques(pos, label, cs);
	    cnt1++;
	 } else {
	    /* greedily find if v_deg is too big */
	    /* order nodes in *decreasing* order of their degrees in star */
	    CoinSort_2(star_deg, star_deg + star_length, star,
		       CoinFirstGreater_2<int,int>());
	    /* find maxl clique greedily, including v */
	    clique_cnt_g += greedy_maximal_clique(cs);
	    cnt2++;
	 }
      } else {
	 cnt3++;
      }
      /* delete v from current_indices */
      cl_del_indices[cl_del_length++] = v;
      scl_delete_node(best_ind, current_nodenum,
		      current_indices, current_degrees, current_values);
      best_ind = scl_choose_next_node(current_nodenum, current_indices,
				      current_degrees, current_values);
      v = current_indices[best_ind];
      v_deg = current_degrees[best_ind];
      v_val = current_values[best_ind];
      largest_star_size = std::max(largest_star_size, v_deg);
   }

   const int clique_cnt = clique_cnt_e + clique_cnt_g;
   assert (nodenum>1||!clique_cnt);

   if (scl_report_result) {
      printf("\nscl Found %i new violated cliques with the star-clique method",
	     clique_cnt);
      printf("\nscl The largest star size was %i (threshold %i)\n",
	     largest_star_size, scl_candidate_length_threshold); // par
      printf("scl Enumeration %i times, found %i maxl cliques\n",
	     cnt1, clique_cnt_e);
      printf("scl Greedy %i times, found %i maxl cliques\n",
	     cnt2, clique_cnt_g);
      printf("scl Skipped a star b/c of small solution value %i times\n",
	     cnt3);

      if (cnt2 == 0)
	 printf("scl    all cliques have been enumerated\n");
      else
	 printf("scl    not all cliques have been eliminated\n");
   }

   delete[] current_indices;
   delete[] current_degrees;
   delete[] current_values;
   delete[] star_deg;
   delete[] label;
}

/*****************************************************************************/

/*===========================================================================*
 * returns the index of the "best" node wrt current_indices, etc.
 *===========================================================================*/

int
CglClique::scl_choose_next_node(const int current_nodenum,
				const int * /* current_indices */,
				const int *current_degrees,
				const double *current_values)
{
   int best = 0;
   int best_deg = current_degrees[0];
   double best_val = current_values[0];
   int i;

   switch (scl_next_node_rule) { // p->par.scl_which_node
   case SCL_MIN_DEGREE: // NOTE: could use stl::min_element
      for (i = 1; i < current_nodenum; i++)
	 if (current_degrees[i] < best_deg) {
	    best = i;
	    best_deg = current_degrees[i];
	 }
      break;
   case SCL_MAX_DEGREE: // NOTE: could use stl::max_element
      for (i = 1; i < current_nodenum; i++)
	 if (current_degrees[i] > best_deg) {
	    best = i;
	    best_deg = current_degrees[i];
	 }
      break;
   case SCL_MAX_XJ_MAX_DEG:
      for (i = 1; i < current_nodenum; i++) {
	 if (current_values[i] > best_val) {
	    best = i;
	    best_val = current_values[i];
	    best_deg = current_degrees[i];
	 } else if (current_values[i] == best_val &&
		    current_degrees[i] > best_deg) {
	    best = i;
	    best_deg = current_degrees[i];
	 }
      }
      break;
   default:
      printf("ERROR: bad starcl_which_node (in scl_choose_next_node\n");
      break;
   }
   return(best);
}

/*****************************************************************************/

/*===========================================================================*
 * Delete the node of index del_ind (this index is wrt current_indices) from
 * the list current_indices (current_degrees and current_values) and based on
 * the graph stored in fgraph decrease the degrees of its neighbors.
 *
 * There are at least 3 nodes in the graph when this function is invoked.
 *
 * Note that the node indices in current_indices are in increasing order,
 * and that this ordering is maintained here.
 *
 * fgraph: IN, pointer to the fractional graph
 * node_node: IN, node_node incidence matrix
 * del_ind: IN, the index of the node to be deleted (wrt to current_indices)
 * pcurrent_nodenum: INOUT, pointer to the current number of nodes
 * current_indices: INOUT, array of current node indices
 * current_degrees: INOUT, array of current node degrees, in the dame order
 *                  as in current_indices
 * current_values: INOUT, array of solution values
 *===========================================================================*/

void
CglClique::scl_delete_node(const int del_ind, int& current_nodenum,
			   int *current_indices, int *current_degrees,
			   double *current_values)
{
   const int v = current_indices[del_ind];

   /* delete the entry corresponding to del_ind from current_indices, 
      current_degrees and current_values */
   memmove(reinterpret_cast<char *>(current_indices + del_ind),
	   reinterpret_cast<char *>(current_indices + (del_ind+1)),
	   (current_nodenum-del_ind-1) * sizeof(int));
   memmove(reinterpret_cast<char *>(current_degrees + del_ind),
	   reinterpret_cast<char *>(current_degrees + (del_ind+1)),
	   (current_nodenum-del_ind-1) * sizeof(int));
   memmove(reinterpret_cast<char *>(current_values + del_ind),
	   reinterpret_cast<char *>(current_values + (del_ind+1)),
	   (current_nodenum-del_ind-1) * sizeof(double));
   current_nodenum--;
   
   /* decrease the degrees of v's neighbors by 1 */
   const bool* node_node_start = node_node + (fgraph.nodenum * v);
   for (int i = 0; i < current_nodenum; ++i)
      if (node_node_start[current_indices[i]])
	 current_degrees[i]--;
}

/*****************************************************************************/

/*===========================================================================*
 * Enumerate all maximal cliques on the nodes in scl_indices. Maximal cliques
 * that are violated are added to the cut list. Returns the number of maximal
 * violated cliques found. The algorithm is recursive.
 *
 * Data members used from CglClique
 *   fgraph: IN, the description of the intersection graph
 *   cl_perm_length: IN, the length of cl_perm_indices
 *   cl_perm_indices: IN, indices of nodes that MUST be in the clique, these
 *                    nodes are supposed to be connected to all nodes in
 *                    cl_indices
 *   cl_length: IN, length of cl_indices and label
 *   cl_indices: IN, indices of nodes on which maximal cliques are sought
 *   cl_del_length: IN, length of cl_del_indices
 *   cl_del_indices: IN, indices of nodes that are already deleted. these
 *                   nodes are tested whether they can be added to a max
 *                   clique discovered in scl_indices. if any of them can be
 *                   added then the clique is not maximal after all...
 *
 * Arguments:
 *   label: INOUT, indicates which nodes are in the clique at the moment
 *   pos: INOUT, position within cl_indices (and label), nodes up to
 *        position pos in cl_indices are permanently labeled (backtrack cannot
 *        change labels)
 *===========================================================================*/

int
CglClique::enumerate_maximal_cliques(int& pos, bool* label, OsiCuts& cs)
{
   const fnode *nodes = fgraph.nodes;
   const int nodenum = fgraph.nodenum;

   int i, j, k, cnt;

   /* starting from position pos, find the first node in cl_indices that
      can be added to the clique, and label it with true */
   while (pos < cl_length) {
      label[pos] = true;
      const bool* node_node_start = node_node + cl_indices[pos] * nodenum;
      for (j = 0; j < pos; j++)
	 if (label[j] && ! node_node_start[cl_indices[j]]) {
	    label[pos] = false;
	    break;
	 }
      if (label[pos++] == true)
	 break;
   }

   /* found counts the number of maximal violated cliques that have been sent
      to the lp under the current level of recursion */
   int found = 0;

   /* if not all nodes are labeled: recurse by setting the last node
      labeled true once to true and once to false;
      otherwise check whether the clique found is maximal and violated */
   if (pos < cl_length) {
      found += enumerate_maximal_cliques(pos, label, cs);
      label[pos-1] = false;
      found += enumerate_maximal_cliques(pos, label, cs);
   } else {
      /* check if the clique can be extended on cl_indices */

      /* copy indices of the clique into coef (not user inds, coef is a tmp) */
      int* coef = new int[cl_length + cl_perm_length];
      for (j = cl_length - 1, cnt = 0; j >= 0; j--)
	 if (label[j])
	    coef[cnt++] = cl_indices[j];
      if (!cnt) {
	 delete[] coef;
	 return(found);
      }
      
      /* check if the clique can be extended on cl_indices */
      for (k = cl_length - 1; k >= 0; k--) {
	 if (!label[k]) {
	    const bool* node_node_start = node_node + cl_indices[k] * nodenum;
	    for (i = cnt - 1; i >= 0; i--)
	       if (!node_node_start[coef[i]])
		  break;
	    /* if k can be added to the clique, return (the clique is not
	       maximal, so it will be or was recorded) */
	    if (i < 0) {
	       delete[] coef;
	       return(found);
	    }
	 }
      }

      /* now the clique is maximal on cl_indices.
	 fill relative indices into coef */
      for (j = 0; j < cl_perm_length; j++)
	 coef[cnt++] = cl_perm_indices[j];
      
      /* check if clique is violated */
      double lhs = 0;
      for (j = 0; j < cnt; j++)
	 lhs += nodes[coef[j]].val;
      if (lhs < 1 + petol) {
	 delete[] coef;
	 return(found);
      }
      
      /* if clique can be extended on cl_del_indices then it can be
	 discarded (was already counted) */
      for (i = 0; i < cl_del_length; i++) {
	 const bool* node_node_start = node_node + cl_del_indices[i]*nodenum;
	 for (j = cnt - 1; j >= 0; j--)
	    if (!node_node_start[coef[j]])
	       break;
	 /* if cl_del_indices[i] can be added to the clique, return */
	 if (j < 0) {
	    delete[] coef;
	    return(found);
	 }
      }

      recordClique(cnt, coef, cs);
      delete[] coef;
      
      ++found;
   }

   return(found);
}

/*****************************************************************************/

/*===========================================================================*
 * Find a violated clique greedily in the given array of indices, starting
 * from pos. Return the number of violated cliques found (1 or 0).
 * This routine overwrites the array indices: those variables in the
 * clique will be shuffled to the beginning of the array (after pos).
 *===========================================================================*/

int
CglClique::greedy_maximal_clique(OsiCuts& cs)
{
   assert(cl_length > 0);
   const fnode *nodes = fgraph.nodes;
   const int nodenum = fgraph.nodenum;
   int i, j;

   int * coef = new int[cl_length + cl_perm_length];
   coef[0] = cl_indices[0];
   int cnt = 1;
   for (j = 1; j < cl_length; j++) {
      const int var = cl_indices[j];
      const bool* node_node_start = node_node + var * nodenum;
      for (i = cnt-1; i >= 0; i--)
	 if (!node_node_start[coef[i]])
	    break;
      if (i < 0)
	 coef[cnt++] = var;
   }

   for (j = 0; j < cl_perm_length; j++)
      coef[cnt++] = cl_perm_indices[j];

   /* now coef contains the clique */
   /* only cliques of size at least 3 are interesting */
   if (cnt < 3) {
      delete[] coef;
      return(0);
   }

   /* compute lhs */
   double lhs = 0;
   for (j = 0; j < cnt; j++)
      lhs += nodes[coef[j]].val;

   if (lhs > 1 + petol) {
      recordClique(cnt, coef, cs);
      delete[] coef;
      return(1);
   }

   delete[] coef;
   return(0);
}

/*****************************************************************************/

/*===========================================================================*
 *===========================================================================*/

void
CglClique::recordClique(const int len, int* indices, OsiCuts& cs)
{
   /* transform relative indices into user indices and order them */
   for (int j = len - 1; j >= 0; j--)
      indices[j] = sp_orig_col_ind[indices[j]];
   std::sort(indices, indices + len);
   OsiRowCut rowcut;
   double* coef = new double[len];
   std::fill(coef, coef + len, 1.0);
   rowcut.setRow(len, indices, coef);
   rowcut.setUb(1.0);
   CoinAbsFltEq equal(1.0e-12);
   cs.insertIfNotDuplicate(rowcut,equal);
   delete[] coef;
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglClique::clone() const
{
  return new CglClique(*this);
}

/*****************************************************************************/
// Create C++ lines to get to current state
std::string
CglClique::generateCpp( FILE * fp) 
{
  CglClique other;
  fprintf(fp,"0#include \"CglClique.hpp\"\n");
  fprintf(fp,"3  CglClique clique;\n");
  std::string types[] = {"SCL_MIN_DEGREE","SCL_MAX_DEGREE",
			 "SCL_MAX_XJ_MAX_DEG"};
  if (scl_next_node_rule!=other.scl_next_node_rule)
    fprintf(fp,"3  clique.setStarCliqueNextNodeMethod(CglClique::%s);\n",
	    types[scl_next_node_rule].c_str());
  else
    fprintf(fp,"4  clique.setStarCliqueNextNodeMethod(CglClique::%s);\n",
	    types[scl_next_node_rule].c_str());
  if (scl_candidate_length_threshold!=other.scl_candidate_length_threshold)
    fprintf(fp,"3  clique.setStarCliqueCandidateLengthThreshold(%d);\n",
	    scl_candidate_length_threshold);
  else
    fprintf(fp,"4  clique.setStarCliqueCandidateLengthThreshold(%d);\n",
	    scl_candidate_length_threshold);
  if (rcl_candidate_length_threshold!=other.rcl_candidate_length_threshold)
    fprintf(fp,"3  clique.setRowCliqueCandidateLengthThreshold(%d);\n",
	    rcl_candidate_length_threshold);
  else
    fprintf(fp,"4  clique.setRowCliqueCandidateLengthThreshold(%d);\n",
	    rcl_candidate_length_threshold);
  if (scl_report_result!=other.scl_report_result)
    fprintf(fp,"3  clique.setStarCliqueReport(%s);\n",
	    scl_report_result ? "true" : "false");
  else
    fprintf(fp,"4  clique.setStarCliqueReport(%s);\n",
	    scl_report_result ? "true" : "false");
  if (rcl_report_result!=other.rcl_report_result)
    fprintf(fp,"3  clique.setRowCliqueReport(%s);\n",
	    rcl_report_result ? "true" : "false");
  else
    fprintf(fp,"4  clique.setRowCliqueReport(%s);\n",
	    rcl_report_result ? "true" : "false");
  if (do_star_clique!=other.do_star_clique)
    fprintf(fp,"3  clique.setDoStarClique(%s);\n",
	    do_star_clique ? "true" : "false");
  else
    fprintf(fp,"4  clique.setDoStarClique(%s);\n",
	    do_star_clique ? "true" : "false");
  if (do_row_clique!=other.do_row_clique)
    fprintf(fp,"3  clique.setDoRowClique(%s);\n",
	    do_row_clique ? "true" : "false");
  else
    fprintf(fp,"4  clique.setDoRowClique(%s);\n",
	    do_row_clique ? "true" : "false");
  if (petol!=other.petol)
    fprintf(fp,"3  clique.setMinViolation(%g);\n",petol);
  else
    fprintf(fp,"4  clique.setMinViolation(%g);\n",petol);
  if (getAggressiveness()!=other.getAggressiveness())
    fprintf(fp,"3  clique.setAggressiveness(%d);\n",getAggressiveness());
  else
    fprintf(fp,"4  clique.setAggressiveness(%d);\n",getAggressiveness());
  return "clique";
}
/*****************************************************************************/

#include "CglProbing.hpp"
CglFakeClique::CglFakeClique(OsiSolverInterface * solver, bool setPacking) :
  CglClique(setPacking,true)
{
  if (solver)
    fakeSolver_ = solver->clone();
  else
    fakeSolver_ = NULL;
  if (fakeSolver_) {
    probing_ = new CglProbing();
    probing_->refreshSolver(fakeSolver_);
  } else {
    probing_ = NULL;
  }
			      
}
// Copy constructor
CglFakeClique::CglFakeClique(const CglFakeClique& rhs)
  : CglClique(rhs)
{
  if (rhs.fakeSolver_) {
    fakeSolver_ = rhs.fakeSolver_->clone();
    probing_ = new CglProbing(*rhs.probing_);
    probing_->refreshSolver(fakeSolver_);
  } else {
    fakeSolver_ = NULL;
    probing_ = NULL;
  }
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglFakeClique::clone() const
{
  return new CglFakeClique(*this);
}

// Destructor
CglFakeClique::~CglFakeClique()
{
  delete fakeSolver_;
  delete probing_;
}
// Assign solver (generator takes over ownership)
void 
CglFakeClique::assignSolver(OsiSolverInterface * fakeSolver)
{
  delete fakeSolver_;
  fakeSolver_ = fakeSolver;
  if (fakeSolver_) {
    delete [] sp_orig_row_ind;
    sp_orig_row_ind=NULL;
  }
  if (probing_)
    probing_->refreshSolver(fakeSolver_);
}
#ifdef CGL_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#endif
// Generate cuts
void
CglFakeClique::generateCuts(const OsiSolverInterface& si, OsiCuts & cs,
			const CglTreeInfo info)
{
  if (fakeSolver_) {
    assert (si.getNumCols()==fakeSolver_->getNumCols());
    fakeSolver_->setColLower(si.getColLower());
    const double * solution = si.getColSolution();
    fakeSolver_->setColSolution(solution);
    fakeSolver_->setColUpper(si.getColUpper());
    // get and set branch and bound cutoff
    double cutoff;
    si.getDblParam(OsiDualObjectiveLimit,cutoff);
    fakeSolver_->setDblParam(OsiDualObjectiveLimit,COIN_DBL_MAX);
#ifdef CGL_HAS_CLP
    OsiClpSolverInterface * clpSolver
      = getClpSolver(fakeSolver_);
    if (CBC_SKIP_CLP_TEST||clpSolver) {
      // fix up fake solver
      const ClpSimplex * siSimplex = clpSolver->getModelPtr();
      // need to set djs
      memcpy(siSimplex->primalColumnSolution(),
	     si.getReducedCost(),si.getNumCols()*sizeof(double));
      fakeSolver_->setDblParam(OsiDualObjectiveLimit,cutoff);
    }
#endif
    const CoinPackedMatrix * matrixByRow = si.getMatrixByRow();
    const double * elementByRow = matrixByRow->getElements();
    const int * column = matrixByRow->getIndices();
    const CoinBigIndex * rowStart = matrixByRow->getVectorStarts();
    const int * rowLength = matrixByRow->getVectorLengths();
    const double * rowUpper = si.getRowUpper();
    const double * rowLower = si.getRowLower();
    
    // Scan all rows looking for possibles
    int numberRows = si.getNumRows();
    double tolerance = 1.0e-3;
    for (int iRow=0;iRow<numberRows;iRow++) {
      CoinBigIndex start = rowStart[iRow];
      CoinBigIndex end = start + rowLength[iRow];
      double upRhs = rowUpper[iRow]; 
      double loRhs = rowLower[iRow]; 
      double sum = 0.0;
      for (CoinBigIndex j=start;j<end;j++) {
	int iColumn=column[j];
	double value = elementByRow[j];
	sum += solution[iColumn]*value;
      }
      if (sum<loRhs-tolerance||sum>upRhs+tolerance) {
	// add as cut
	OsiRowCut rc;
	rc.setLb(loRhs);
	rc.setUb(upRhs);
	rc.setRow(end-start,column+start,elementByRow+start,false);
	CoinAbsFltEq equal(1.0e-12);
	cs.insertIfNotDuplicate(rc,equal);
      }
    }
    CglClique::generateCuts(*fakeSolver_,cs,info);
    if (probing_) {
      probing_->generateCuts(*fakeSolver_,cs,info);
    }
  } else {
    // just use real solver
    CglClique::generateCuts(si,cs,info);
  }
}
