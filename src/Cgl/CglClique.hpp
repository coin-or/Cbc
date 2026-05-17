// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef _CglClique_h_
#define _CglClique_h_

#include "CglCutGenerator.hpp"

//class OsiCuts;
//class OsiSolverInterface;

class CGLLIB_EXPORT CglClique : public CglCutGenerator {

    friend CGLLIB_EXPORT void CglCliqueUnitTest(const OsiSolverInterface * siP,
				  const std::string mpdDir );
public:
    /// Copy constructor
    CglClique(const CglClique& rhs);
    /// Clone
    virtual CglCutGenerator * clone() const;

    /// Assignment operator
    CglClique& operator=(const CglClique& rhs);
   
public:
    
    virtual void
    generateCuts(const OsiSolverInterface& si, OsiCuts & cs,
		 const CglTreeInfo info = CglTreeInfo());
   
    /**@name Constructors and destructors */
    //@{
    /** Default constructor.
	If the setPacking argument is set to true then CglClique will assume that the
	problem in the solverinterface passed to the generateCuts() method
	describes a set packing problem, i.e.,
	- all variables are binary
	- the matrix is a 0-1 matrix
	- all constraints are '= 1' or '<= 1'
       
	Otherwise the user can use the considerRows() method to set the list of
	clique rows, that is,
	- all coeffs corresponding to binary variables at fractional level is 1
	- all other coeffs are non-negative
	- the constraint is '= 1' or '<= 1'.

	If the user does not set the list of clique rows then CglClique will
	start the generateCuts() methods by scanning the matrix for them.
	Also justOriginalRows can be set to true to limit clique creation
    */
    CglClique(bool setPacking = false, bool justOriginalRows = false);
    /// Destructor
    virtual ~CglClique() {}
    /// Create C++ lines to get to current state
    virtual std::string generateCpp( FILE * fp);

    void considerRows(const int numRows, const int* rowInd);

public:
    /** possible choices for selecting the next node in the star clique search
     */
    enum scl_next_node_method {
	SCL_MIN_DEGREE,
	SCL_MAX_DEGREE,
	SCL_MAX_XJ_MAX_DEG
    };

    void setStarCliqueNextNodeMethod(scl_next_node_method method) {
	scl_next_node_rule = method;
    }

    void setStarCliqueCandidateLengthThreshold(int maxlen) {
	scl_candidate_length_threshold = maxlen;
    }
    void setRowCliqueCandidateLengthThreshold(int maxlen) {
	rcl_candidate_length_threshold = maxlen;
    }

    void setStarCliqueReport(bool yesno = true) { scl_report_result = yesno; }
    void setRowCliqueReport(bool yesno = true) { rcl_report_result = yesno; }

    void setDoStarClique(bool yesno = true) { do_star_clique = yesno; }
    void setDoRowClique(bool yesno = true) { do_row_clique = yesno; }

    void setMinViolation(double minviol) { petol = minviol; }
    double getMinViolation() const { return petol; }
    /// Maximum number of binaries for looking at all
    inline void setMaxNumber(int value) { maxNumber_ = value; }

private:

    struct frac_graph ;
    friend struct frac_graph ;

    /** A node of the fractional graph. There is a node for every variable at
	fractional level. */
    struct fnode {
	/** pointer into all_nbr */
	int          *nbrs;
	/** 1-x_i-x_j, needed for odd holes, in the same order as the adj list,
	    pointer into all_edgecost */
	double       *edgecosts;
	/** degree of the node */
	int           degree;
	/** the fractional value of the variable corresponding to this node */
	double        val;
    };

    /** A graph corresponding to a fractional solution of an LP. Two nodes are
	adjacent iff their columns are non-orthogonal. */
    struct frac_graph {
	/** # of nodes = # of fractional values in the LP solution */
	int    nodenum;
	/** # of edges in the graph */
	int    edgenum;
	/** density= edgenum/(nodenum choose 2) */
	double density;
	int    min_deg_node;
	int    min_degree;
	int    max_deg_node;
	int    max_degree;
	/** The array of the nodes in the graph */
	fnode  *nodes;
	/** The array of all the neighbors. First the indices of the nodes
	    adjacent to node 0 are listed, then those adjacent to node 1, etc. */
	int    *all_nbr;
	/** The array of the costs of the edges going to the neighbors */
	double *all_edgecost;

	frac_graph() :
	    nodenum(0), edgenum(0), density(0),
	    min_deg_node(0), min_degree(0), max_deg_node(0), max_degree(0),
	    nodes(0), all_nbr(0), all_edgecost(0) {}
    };

protected:
    /** An indicator showing whether the whole matrix in the solverinterface is
	a set packing problem or not */
    bool setPacking_;
    /// True if just look at original rows
    bool justOriginalRows_;
    /** pieces of the set packing part of the solverinterface */
    int sp_numrows;
    int* sp_orig_row_ind;
    int sp_numcols;
    int* sp_orig_col_ind;
    double* sp_colsol;
    int* sp_col_start;
    int* sp_col_ind;
    int* sp_row_start;
    int* sp_row_ind;

    /** the intersection graph corresponding to the set packing problem */
    frac_graph fgraph;
    /** the node-node incidence matrix of the intersection graph. */
    bool* node_node;

    /** The primal tolerance in the solverinterface. */
    double petol;
    /// Maximum number of binaries for looking at all
    int maxNumber_; 

    /** data for the star clique algorithm */

    /** Parameters */
    /**@{*/
    /** whether to do the row clique algorithm or not. */
    bool do_row_clique;
    /** whether to do the star clique algorithm or not. */
    bool do_star_clique;

    /** How the next node to be added to the star clique should be selected */
    scl_next_node_method scl_next_node_rule;
    /** In the star clique method the maximal length of the candidate list
	(those nodes that are in a star, i.e., connected to the center of the
	star) to allow complete enumeration of maximal cliques. Otherwise a
	greedy algorithm is used. */
    int scl_candidate_length_threshold;
    /** whether to give a detailed statistics on the star clique method */
    bool scl_report_result;

    /** In the row clique method the maximal length of the candidate list
	(those nodes that can extend the row clique, i.e., connected to all
	nodes in the row clique) to allow complete enumeration of maximal
	cliques. Otherwise a greedy algorithm is used. */
    int rcl_candidate_length_threshold;
    /** whether to give a detailed statistics on the row clique method */
    bool rcl_report_result;
    /**@}*/

    /** variables/arrays that are used across many methods */
    /**@{*/
    /** List of indices that must be in the to be created clique. This is just
	a pointer, it is never new'd and therefore does not need to be
	delete[]'d either. */
    const int* cl_perm_indices;
    /** The length of cl_perm_indices */
    int cl_perm_length;

    /** List of indices that should be considered for extending the ones listed
	in cl_perm_indices. */
    int* cl_indices;
    /** The length of cl_indices */
    int cl_length;

    /** An array of nodes discarded from the candidate list. These are
	rechecked when a maximal clique is found just to make sure that the
	clique is really maximal. */
    int* cl_del_indices;
    /** The length of cl_del_indices */
    int cl_del_length;

    /**@}*/

private:
    /** Scan through the variables and select those that are binary and are at
	a fractional level. */
    void selectFractionalBinaries(const OsiSolverInterface& si);
    /** Scan through the variables and select those that are at a fractional
	level. We already know that everything is binary. */
    void selectFractionals(const OsiSolverInterface& si);
    /**  */
    void selectRowCliques(const OsiSolverInterface& si,int numOriginalRows);
    /**  */
    void createSetPackingSubMatrix(const OsiSolverInterface& si);
    /**  */
    void createFractionalGraph();
    /**  */
    int createNodeNode();
    /**  */
    void deleteSetPackingSubMatrix();
    /**  */
    void deleteFractionalGraph();
    /**  */
    void find_scl(OsiCuts& cs);
    /**  */
    void find_rcl(OsiCuts& cs);
    /**  */
    int scl_choose_next_node(const int current_nodenum,
			     const int *current_indices,
			     const int *current_degrees,
			     const double *current_values);
    /**  */
    void scl_delete_node(const int del_ind, int& current_nodenum,
			 int *current_indices, int *current_degrees,
			 double *current_values);
    /**  */
    int enumerate_maximal_cliques(int& pos, bool* scl_label, OsiCuts& cs);
    /**  */
    int greedy_maximal_clique(OsiCuts& cs);
    /**  */
    void recordClique(const int len, int* indices, OsiCuts& cs);
};
//#############################################################################
/** A function that tests the methods in the CglClique class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
    
CGLLIB_EXPORT
void CglCliqueUnitTest(const OsiSolverInterface * siP,
		       const std::string mpdDir);

/// This works on a fake solver i.e. invented rows
class CglProbing;
class CGLLIB_EXPORT CglFakeClique : public CglClique {
  
public:
  /// Copy constructor
  CglFakeClique(const CglFakeClique& rhs);
  /// Clone
  virtual CglCutGenerator * clone() const;
  
  /// Assignment operator
  CglFakeClique& operator=(const CglFakeClique& rhs);
  
  virtual void
  generateCuts(const OsiSolverInterface& si, OsiCuts & cs,
	       const CglTreeInfo info = CglTreeInfo());
  
  /**@name Constructors and destructors */
  //@{
  /** Default constructor.
      If the setPacking argument is set to true then CglFakeClique will assume that the
      problem in the solverinterface passed to the generateCuts() method
      describes a set packing problem, i.e.,
      - all variables are binary
      - the matrix is a 0-1 matrix
      - all constraints are '= 1' or '<= 1'
      
      Otherwise the user can use the considerRows() method to set the list of
      clique rows, that is,
      - all coeffs corresponding to binary variables at fractional level is 1
      - all other coeffs are non-negative
      - the constraint is '= 1' or '<= 1'.
      
      If the user does not set the list of clique rows then CglFakeClique will
      start the generateCuts() methods by scanning the matrix for them.
  */
  CglFakeClique(OsiSolverInterface * solver=NULL,bool setPacking = false);
  /// Destructor
  virtual ~CglFakeClique();
  /// Assign solver (generator takes over ownership)
  void assignSolver(OsiSolverInterface * fakeSolver);
protected:
  /// fake solver to use
  OsiSolverInterface * fakeSolver_;
  /// Probing object
  CglProbing * probing_;
};

#endif
