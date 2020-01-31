/* $Id: CbcSymmetry.hpp 1033 2013-12-14 19:34:28Z pbelotti $
 *
 * Name:    Hacked from CouenneProblem.hpp
 * Author:  Pietro Belotti, Lehigh University
 *          Andreas Waechter, IBM
 * Purpose: define the class CouenneProblem
 *
 * (C) Carnegie-Mellon University, 2006-11.
 * This file is licensed under the Eclipse Public License (EPL)
 */
/*
  If this is much used then we could improve build experience
  Download nauty - say to /disk/nauty25r9
  In that directory ./configure --enable-tls --enable-wordsize=32
  make
  copy nauty.a to libnauty.a

  In Cbc's configure 
  add -DCOIN_HAS_NTY to CXXDEFS
  add -I/disk/nauty25r9 to CXXDEFS or ADD_CXXFLAGS
  add -L/disk/nauty25r9 -lnauty to LDFLAGS

  If you wish to use Traces rather than nauty then add -DNTY_TRACES

  To use it is -orbit on

  Nauty has this -
*   Permission 
*   is hereby given for use and/or distribution with the exception of        *
*   sale for profit or application with nontrivial military significance.    *
  so you can use it internally even if you are a company.
 */
#ifndef CBC_SYMMETRY_HPP
#define CBC_SYMMETRY_HPP
extern "C" {
#include "nauty.h"
#include "nausparse.h"
#ifdef NTY_TRACES
#include "traces.h"
#endif
}

#include <vector>
#include <map>
#include <string.h>

#include "CbcModel.hpp"

class OsiObject;
// when to give up (depth since last success)
#ifndef NTY_BAD_DEPTH
#define NTY_BAD_DEPTH 10
#endif
class CbcNauty;

#define COUENNE_HACKED_EPS 1.e-07
#define COUENNE_HACKED_EPS_SYMM 1e-8
#define COUENNE_HACKED_EXPRGROUP 8

/** Class to deal with symmetry
 *
 *  Hacked from Couenne
 *  Thanks, but it had been nice to make sure that there are no symbol collisions when building Couenne with this Cbc.
 */

class CbcSymmetry {
  class Node {
    int index;
    double coeff;
    double lb;
    double ub;
    int color;
    int code;
    int sign;

  public:
    void node(int, double, double, double, int, int);
    inline void color_vertex(register int k) { color = k; }
    inline int get_index() const { return index; }
    inline double get_coeff() const { return coeff; }
    inline double get_lb() const { return lb; }
    inline double get_ub() const { return ub; }
    inline int get_color() const { return color; }
    inline int get_code() const { return code; }
    inline int get_sign() const { return sign; }
    inline void bounds(register double a, register double b)
    {
      lb = a;
      ub = b;
    }
  };

  struct myclass0 {
    inline bool operator()(register const Node &a, register const Node &b)
    {

      return ((a.get_code() < b.get_code()) || ((a.get_code() == b.get_code() && ((a.get_coeff() < b.get_coeff() - COUENNE_HACKED_EPS_SYMM) || ((fabs(a.get_coeff() - b.get_coeff()) < COUENNE_HACKED_EPS_SYMM) && ((a.get_lb() < b.get_lb() - COUENNE_HACKED_EPS_SYMM) || ((fabs(a.get_lb() - b.get_lb()) < COUENNE_HACKED_EPS_SYMM) && ((a.get_ub() < b.get_ub() - COUENNE_HACKED_EPS_SYMM) || ((fabs(a.get_ub() - b.get_ub()) < COUENNE_HACKED_EPS_SYMM) && ((a.get_index() < b.get_index())))))))))));
    }
  };

  struct myclass {
    inline bool operator()(register const Node &a, register const Node &b)
    {
      return (a.get_index() < b.get_index());
    }
  };

  struct less_than_str {
    inline bool operator()(register const char *a, register const char *b) const
    {
      return strcmp(a, b) < 0;
    }
  };

public:
  /**@name Constructors and destructors */
  //@{
  /// Default constructor
  CbcSymmetry();

  /// Copy constructor
  CbcSymmetry(const CbcSymmetry &);

  /// Assignment operator
  CbcSymmetry &operator=(const CbcSymmetry &rhs);

  /// Destructor
  ~CbcSymmetry();
  //@}

  // Symmetry Info

  std::vector< int > *Find_Orbit(int) const;

  myclass0 node_sort;
  myclass index_sort;

  void Compute_Symmetry() const;
  int statsOrbits(CbcModel *model, int type) const;
  //double timeNauty () const;
  void Print_Orbits() const;
  void fillOrbits();
  /// Fixes variables using orbits (returns number fixed)
  int orbitalFixing(OsiSolverInterface *solver);
  inline int *whichOrbit()
  {
    return numberUsefulOrbits_ ? whichOrbit_ : NULL;
  }
  inline int numberUsefulOrbits() const
  {
    return numberUsefulOrbits_;
  }
  inline int numberUsefulObjects() const
  {
    return numberUsefulObjects_;
  }
  int largestOrbit(const double *lower, const double *upper) const;
  void ChangeBounds(const double *lower, const double *upper,
    int numberColumns, bool justFixedAtOne) const;
  inline bool compare(register Node &a, register Node &b) const;
  CbcNauty *getNtyInfo() { return nauty_info_; }

  // bool node_sort (  Node  a, Node  b);
  // bool index_sort (  Node  a, Node  b);

  /// empty if no NTY, symmetry data structure setup otherwise
  void setupSymmetry(CbcModel * model);

private:
  mutable std::vector< Node > node_info_;
  mutable CbcNauty *nauty_info_;
  int numberColumns_;
  int numberUsefulOrbits_;
  int numberUsefulObjects_;
  int *whichOrbit_;
};

class CbcNauty {

public:
  enum VarStatus { FIX_AT_ZERO,
    FIX_AT_ONE,
    FREE };

  /**@name Constructors and destructors */
  //@{
private:
  /// Default constructor
  CbcNauty();

public:
  /// Normal constructor (if dense - NULLS)
  CbcNauty(int n, const size_t *v, const int *d, const int *e);

  /// Copy constructor
  CbcNauty(const CbcNauty &);

  /// Assignment operator
  CbcNauty &operator=(const CbcNauty &rhs);

  /// Destructor
  ~CbcNauty();
  //@}

  void addElement(int ix, int jx);
  void clearPartitions();
  void computeAuto();
  void deleteElement(int ix, int jx);
  void color_node(int ix, int color) { vstat_[ix] = color; }
  void insertRHS(int rhs, int cons) { constr_rhs.insert(std::pair< int, int >(rhs, cons)); }

  double getGroupSize() const;
  //int getNautyCalls() const { return nautyCalls_; }
  //double getNautyTime() const { return nautyTime_; }

  int getN() const { return n_; }

  int getNumGenerators() const;
  int getNumOrbits() const;

  /// Returns the orbits in a "convenient" form
  std::vector< std::vector< int > > *getOrbits() const;

  void getVstat(double *v, int nv);
  inline bool isSparse() const
  {
    return GSparse_ != NULL;
  }
  inline int errorStatus() const
#ifndef NTY_TRACES
  {
    return stats_->errstatus;
  }
#else
  {
    return 0;
  }
#endif
  /// Pointer to options
  inline optionblk *options() const
  {
    return options_;
  }
  /**
   * Methods to classify orbits.  Not horribly efficient, but gets the job done
   */
  //  bool isAllFixOneOrbit(const std::vector<int> &orbit) const;
  // bool isAllFreeOrbit(const std::vector<int> &orbit) const;
  //bool isAutoComputed() const { return autoComputed_; }
  //bool isConstraintOrbit(const std::vector<int> &orbit) const;
  //bool isMixedFreeZeroOrbit(const std::vector<int> &orbit) const;
  //void makeFree(int ix) { vstat_[ix] = FREE; }

  void setWriteAutoms(const std::string &afilename);
  void unsetWriteAutoms();

private:
  // The base nauty stuff
  graph *G_;
  sparsegraph *GSparse_;
  int *lab_;
  int *ptn_;
  set *active_;
  int *orbits_;
#ifndef NTY_TRACES
  optionblk *options_;
  statsblk *stats_;
#else
  TracesOptions *options_;
  TracesStats *stats_;
#endif
  setword *workspace_;
  int worksize_;
  int m_;
  int n_;
  size_t nel_;
  graph *canonG_;

  bool autoComputed_;

  int *vstat_;

  //static int nautyCalls_;
  //static double nautyTime_;

  std::multimap< int, int > constr_rhs;
  std::multimap< int, int >::iterator it;

  std::pair< std::multimap< int, int >::iterator,
    std::multimap< int, int >::iterator >
    ret;

  // File pointer for automorphism group
  FILE *afp_;
};

/** Branching object for Orbital branching

    Variable_ is the set id number (redundant, as the object also holds a
    pointer to the set. 
 */
class CbcOrbitalBranchingObject : public CbcBranchingObject {

public:
  // Default Constructor
  CbcOrbitalBranchingObject();

  // Useful constructor
  CbcOrbitalBranchingObject(CbcModel *model, int column,
    int way,
    int numberExtra, const int *extraToZero);

  // Copy constructor
  CbcOrbitalBranchingObject(const CbcOrbitalBranchingObject &);

  // Assignment operator
  CbcOrbitalBranchingObject &operator=(const CbcOrbitalBranchingObject &rhs);

  /// Clone
  virtual CbcBranchingObject *clone() const;

  // Destructor
  virtual ~CbcOrbitalBranchingObject();

  using CbcBranchingObject::branch;
  /// Does next branch and updates state
  virtual double branch();
  /** Update bounds in solver as in 'branch' and update given bounds.
        branchState is -1 for 'down' +1 for 'up' */
  virtual void fix(OsiSolverInterface *solver,
    double *lower, double *upper,
    int branchState) const;

  /** Reset every information so that the branching object appears to point to
        the previous child. This method does not need to modify anything in any
        solver. */
  virtual void previousBranch()
  {
    CbcBranchingObject::previousBranch();
  }

  using CbcBranchingObject::print;
  /** \brief Print something about branch - only if log level high
    */
  virtual void print();

  /** Return the type (an integer identifier) of \c this */
  virtual CbcBranchObjType type() const
  {
    return SoSBranchObj;
  }

  /** Compare the original object of \c this with the original object of \c
        brObj. Assumes that there is an ordering of the original objects.
        This method should be invoked only if \c this and brObj are of the same
        type.
        Return negative/0/positive depending on whether \c this is
        smaller/same/larger than the argument.
    */
  virtual int compareOriginalObject(const CbcBranchingObject *brObj) const;

  /** Compare the \c this with \c brObj. \c this and \c brObj must be os the
        same type and must have the same original object, but they may have
        different feasible regions.
        Return the appropriate CbcRangeCompare value (first argument being the
        sub/superset if that's the case). In case of overlap (and if \c
        replaceIfOverlap is true) replace the current branching object with one
        whose feasible region is the overlap.
     */
  virtual CbcRangeCompare compareBranchingObject(const CbcBranchingObject *brObj, const bool replaceIfOverlap = false);

private:
  /// Column to go to 1
  int column_;
  /// Number (without column) going to zero on down branch
  int numberOther_;
  /// Number extra
  int numberExtra_;
  /// Fix to zero
  int *fixToZero_;
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
