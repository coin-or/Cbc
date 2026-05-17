// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <numeric>
#include <cassert>

#include "CoinPackedMatrix.hpp"
#include "OsiSolverInterface.hpp"
#include "CglClique.hpp"

/*****************************************************************************/

/*===========================================================================*
  Scan through the variables and select those that are binary and are at a
  fractional level.
 *===========================================================================*/
void
CglClique::selectFractionalBinaries(const OsiSolverInterface& si)
{
   // extract the primal tolerance from the solver
   double lclPetol = 0.0;
   si.getDblParam(OsiPrimalTolerance, lclPetol);
   const int numcols = si.getNumCols();
   const char * intVar = si.getColType();
   if (petol<0.0) {
     // do all if not too many
     int n=0;
     for (int i = 0; i < numcols; ++i) {
       if (intVar[i]==1)
	 n++;
     }
     if (n<maxNumber_)
       lclPetol=-1.0e-5;
   }
   const double* x = si.getColSolution();
   std::vector<int> fracind;
   int i;
   for (i = 0; i < numcols; ++i) {
      if (intVar[i]==1 && x[i] > lclPetol && x[i] < 1-petol)
	 fracind.push_back(i);
   }
   sp_numcols = static_cast<int>(fracind.size());
   sp_orig_col_ind = new int[sp_numcols];
   sp_colsol = new double[sp_numcols];
   for (i = 0; i < sp_numcols; ++i) {
      sp_orig_col_ind[i] = fracind[i];
      sp_colsol[i] = x[fracind[i]];
   }
}

/*****************************************************************************/

/*===========================================================================*
  Scan through the variables and select those that are at a fractional
  level. We already know that everything is binary.
 *===========================================================================*/

void
CglClique::selectFractionals(const OsiSolverInterface& si)
{
   // extract the primal tolerance from the solver
   double lclPetol = 0.0;
   si.getDblParam(OsiPrimalTolerance, lclPetol);

   const int numcols = si.getNumCols();
   const double* x = si.getColSolution();
   std::vector<int> fracind;
   int i;
   for (i = 0; i < numcols; ++i) {
      if (x[i] > lclPetol && x[i] < 1-lclPetol)
	 fracind.push_back(i);
   }
   sp_numcols = static_cast<int>(fracind.size());
   sp_orig_col_ind = new int[sp_numcols];
   sp_colsol = new double[sp_numcols];
   for (i = 0; i < sp_numcols; ++i) {
      sp_orig_col_ind[i] = fracind[i];
      sp_colsol[i] = x[fracind[i]];
   }
}

/*****************************************************************************/

/*===========================================================================*
 *===========================================================================*/

void
CglClique::selectRowCliques(const OsiSolverInterface& si,int numOriginalRows)
{
   const int numrows = si.getNumRows();
#ifndef INTEL_COMPILER
   std::vector<int> clique(numrows, 1);
#else
   std::vector<int> clique(numrows);
#endif

   int i, j, k;
   
   // First scan through the binary fractional variables and see where do they
   // have a 1 coefficient
   const CoinPackedMatrix& mcol = *si.getMatrixByCol();
   for (j = 0; j < sp_numcols; ++j) {
      const CoinShallowPackedVector& vec = mcol.getVector(sp_orig_col_ind[j]);
      const int* ind = vec.getIndices();
      const double* elem = vec.getElements();
      for (i = vec.getNumElements() - 1; i >= 0; --i) {
	 if (elem[i] != 1.0) {
	    clique[ind[i]] = 0;
	 }
      }
   }

   // Now check the sense and rhs (by checking rowupper) and the rest of the
   // coefficients 
   const CoinPackedMatrix& mrow = *si.getMatrixByRow();
   const double* rub = si.getRowUpper();
   for (i = 0; i < numrows; ++i) {
      if (rub[i] != 1.0||i>=numOriginalRows) {
	 clique[i] = 0;
	 continue;
      }
      if (clique[i] == 1) {
	 const CoinShallowPackedVector& vec = mrow.getVector(i);
	 const double* elem = vec.getElements();
	 for (j = vec.getNumElements() - 1; j >= 0; --j) {
	    if (elem[j] < 0) {
	       clique[i] = 0;
	       break;
	    }
	 }
      }
   }

   // Finally collect the still standing rows into sp_orig_row_ind
   sp_numrows = std::accumulate(clique.begin(), clique.end(), 0);
   sp_orig_row_ind = new int[sp_numrows];
   for (i = 0, k = 0; i < numrows; ++i) {
      if (clique[i] == 1) {
	 sp_orig_row_ind[k++] = i;
      }
   }
}

/*****************************************************************************/

/*===========================================================================*
  Create the set packing submatrix
 *===========================================================================*/
void
CglClique::createSetPackingSubMatrix(const OsiSolverInterface& si)
{
   sp_col_start = new int[sp_numcols+1];
   sp_row_start = new int[sp_numrows+1];
   std::fill(sp_col_start, sp_col_start + (sp_numcols+1), 0);
   std::fill(sp_row_start, sp_row_start + (sp_numrows+1), 0);

   int i, j;

   const CoinPackedMatrix& mcol = *si.getMatrixByCol();
   const int numrows = si.getNumRows();
   int* clique = new int[numrows];
   std::fill(clique, clique+numrows, -1);
   for (i = 0; i < sp_numrows; ++i)
      clique[sp_orig_row_ind[i]] = i;

   for (j = 0; j < sp_numcols; ++j) {
      const CoinShallowPackedVector& vec = mcol.getVector(sp_orig_col_ind[j]);
      const int* ind = vec.getIndices();
      for (i = vec.getNumElements() - 1; i >= 0; --i) {
	 if (clique[ind[i]] >= 0) {
	    ++sp_col_start[j];
	    ++sp_row_start[clique[ind[i]]];
	 }
      }
   }

   std::partial_sum(sp_col_start, sp_col_start+sp_numcols, sp_col_start);
   std::rotate(sp_col_start, sp_col_start+sp_numcols,
	       sp_col_start + (sp_numcols+1));
   std::partial_sum(sp_row_start, sp_row_start+sp_numrows, sp_row_start);
   std::rotate(sp_row_start, sp_row_start+sp_numrows,
	       sp_row_start + (sp_numrows+1));
   const int nzcnt = sp_col_start[sp_numcols];
   assert(nzcnt == sp_row_start[sp_numrows]);
/*
  Now create the vectors with row indices for each column (sp_col_ind) and
  column indices for each row (sp_row_ind). It turns out that
  CoinIsOrthogonal assumes that the row indices for a given column are listed
  in ascending order. This is *not* a solver-independent assumption! At best,
  one can hope that the underlying solver will produce an index vector that's
  either ascending or descending. Under that assumption, compare the first
  and last entries and proceed accordingly. Eventually some solver will come
  along that hands back an index vector in random order, and CoinIsOrthogonal
  will break.  Until then, try and avoid the cost of a sort.
*/
   sp_col_ind = new int[nzcnt];
   sp_row_ind = new int[nzcnt];
   int last=0;
   for (j = 0; j < sp_numcols; ++j) {
      const CoinShallowPackedVector& vec = mcol.getVector(sp_orig_col_ind[j]);
      const int len = vec.getNumElements();
      const int* ind = vec.getIndices();
      if (ind[0] < ind[len-1]) {
	for (i = 0; i < len; ++i) {
	   const int sp_row = clique[ind[i]];
	   if (sp_row >= 0) {
	      sp_col_ind[sp_col_start[j]++] = sp_row;
	      sp_row_ind[sp_row_start[sp_row]++] = j;
	   }
	}
      }
      else {
	for (i = len-1; i >= 0; --i) {
	   const int sp_row = clique[ind[i]];
	   if (sp_row >= 0) {
	      sp_col_ind[sp_col_start[j]++] = sp_row;
	      sp_row_ind[sp_row_start[sp_row]++] = j;
	   }
	}
      }
      // sort
      std::sort(sp_col_ind+last,sp_col_ind+sp_col_start[j]);
      last=sp_col_start[j];
   }
   std::rotate(sp_col_start, sp_col_start+sp_numcols,
	       sp_col_start + (sp_numcols+1));
   sp_col_start[0] = 0;
   std::rotate(sp_row_start, sp_row_start+sp_numrows,
	       sp_row_start + (sp_numrows+1));
   sp_row_start[0] = 0;

   delete[] clique;
}

/*****************************************************************************/

static inline bool
CoinIsOrthogonal(const int* first0, const int* last0,
		 const int* first1, const int* last1)
{
   while (first0 != last0 && first1 != last1) {
      if (*first0 == *first1)
	 return false;
      if (*first0 < *first1)
	 ++first0;
      else
	 ++first1;
   }
   return true;
}

/*===========================================================================*
  Build up the fractional graph
 *===========================================================================*/

void
CglClique::createFractionalGraph()
{
   // fgraph.edgenum is filled when createNodeNode is invoked
   fgraph.nodenum = sp_numcols;
   fgraph.all_nbr = new int[2*fgraph.edgenum];
   fgraph.nodes = new fnode[sp_numcols+1];

   int *all_nbr = fgraph.all_nbr;
   fnode *nodes = fgraph.nodes;
   int min_degree, max_degree, min_deg_node, max_deg_node;

#  ifdef ZEROFAULT
   // May be read below even if sp_numcols == 0
   nodes[0].degree = 0 ;
#  endif

   int i, j, total_deg, old_total;

   /*========================================================================*
      Construct the adjacency lists (neighbors) of the nodes in fgraph.
      Two nodes are adjacent iff the columns corresponding to them are
      non-orthogonal.
    *========================================================================*/

   for ( i = 0, total_deg = 0; i < sp_numcols; i++ ) {
      old_total = total_deg;
      const bool* node_node_i = node_node + i * sp_numcols;
      for ( j = 0; j < sp_numcols; j++ ) {
	 if ( node_node_i[j] ) {
	    all_nbr[total_deg++] = j;
	 }
      }
      nodes[i].val = sp_colsol[i];
      nodes[i].degree = total_deg - old_total;
      nodes[i].nbrs = all_nbr + old_total;
   }

   // not used fgraph.density = static_cast<double> (total_deg) / (sp_numcols * (sp_numcols-1));

   /*========================================================================*
     Compute the min and max degree.
    *========================================================================*/
   min_deg_node = 0; max_deg_node = 0;
   min_degree = max_degree = nodes[0].degree;
   for ( i = 0; i < sp_numcols; i++ ) {
      if ( nodes[i].degree < min_degree ) {
	 min_deg_node = i;
	 min_degree = nodes[i].degree;
      }
      if ( nodes[i].degree > max_degree ) {
	 max_deg_node = i;
	 max_degree = nodes[i].degree;
      }
   }
   fgraph.min_degree = min_degree;
   fgraph.max_degree = max_degree;
   fgraph.min_deg_node = min_deg_node;
   fgraph.max_deg_node = max_deg_node;
}

/*****************************************************************************/

/*===========================================================================*
 * Construct the node-node incidence matrix from the fractional graph.
 *===========================================================================*/
int
CglClique::createNodeNode()
{
   node_node = new bool[sp_numcols * sp_numcols];
   std::fill(node_node, node_node + sp_numcols * sp_numcols, false);

   int i, j;
   int edgenum = 0;
   for (i = 0; i < sp_numcols; ++i) {
      for (j = i+1; j < sp_numcols; ++j) {
	 if (! CoinIsOrthogonal(sp_col_ind + sp_col_start[i],
				sp_col_ind + sp_col_start[i+1],
				sp_col_ind + sp_col_start[j],
				sp_col_ind + sp_col_start[j+1]) ) {
	    node_node[i * sp_numcols + j] = true;
	    node_node[j * sp_numcols + i] = true;
	    ++edgenum;
	 }
      }
   }
   return edgenum;
}

/*****************************************************************************/

/*===========================================================================*
 * Cleanup routines
 *===========================================================================*/
void
CglClique::deleteSetPackingSubMatrix()
{
   delete[] sp_orig_row_ind; sp_orig_row_ind = 0; 
   delete[] sp_orig_col_ind; sp_orig_col_ind = 0; 
   delete[] sp_colsol;	     sp_colsol = 0;       
   delete[] sp_col_start;    sp_col_start = 0;    
   delete[] sp_col_ind;	     sp_col_ind = 0;      
   delete[] sp_row_start;    sp_row_start = 0;    
   delete[] sp_row_ind;      sp_row_ind = 0;
}

void
CglClique::deleteFractionalGraph()
{
   fgraph.nodenum = 0;
   fgraph.edgenum = 0;
   fgraph.density = 0;
   fgraph.min_deg_node = 0;
   fgraph.min_degree = 0;
   fgraph.max_deg_node = 0;
   fgraph.max_degree = 0;
   delete[] fgraph.all_nbr;      fgraph.all_nbr = 0;
   delete[] fgraph.nodes;        fgraph.nodes = 0;
   delete[] fgraph.all_edgecost; fgraph.all_edgecost = 0;
}
