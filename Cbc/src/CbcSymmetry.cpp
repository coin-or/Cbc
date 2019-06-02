/* $Id: CbcSymmetry.cpp 925 2012-11-27 19:11:04Z stefan $
 *
 * hacked from CouenneSymmetry.cpp
 * Name:    CbcSymmetry.cpp
 * Author:  Jim Ostrowski (the good bits - rest JJHF)
 * Purpose: methods for exploiting symmetry
 * Date:    October 13, 2010
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */
//#define PRINT_MORE 1

#include "CbcConfig.h"

#ifdef COIN_HAS_NTY

extern "C" {
#include "nauty.h"
#include "nausparse.h"
#ifdef NTY_TRACES
#include "traces.h"
#endif
}

#include <stdio.h>
#include <cassert>
#include <vector>
#include <algorithm>
#include <ostream>
#include <iterator>

#include "CbcSymmetry.hpp"
#include "CbcBranchingObject.hpp"
#include "CoinTime.hpp"
#define NAUTY_MAX_LEVEL 0
#if NAUTY_MAX_LEVEL
extern int nauty_maxalllevel;
#endif
/* Deliberately not threadsafe to save effort
   Just for statistics
   and not worth gathering across threads
   can redo later
 */
static int nautyBranchCalls_ = 0;
static int lastNautyBranchCalls_ = 0;
static int nautyBranchSucceeded_ = 0;
static int nautyFixCalls_ = 0;
static int lastNautyFixCalls_ = 0;
static int nautyFixSucceeded_ = 0;
static double nautyTime_ = 0.0;
static double nautyFixes_ = 0.0;
static double nautyOtherBranches_ = 0.0;

void CbcSymmetry::Node::node(int i, double c, double l, double u, int cod, int s)
{
  index = i;
  coeff = c;
  lb = l;
  ub = u;
  color = -1;
  code = cod;
  sign = s;
}

inline bool CbcSymmetry::compare(register Node &a, register Node &b) const
{
  if (a.get_code() == b.get_code())
    if (a.get_coeff() == b.get_coeff())
      if (a.get_sign() == b.get_sign())
        if (fabs(a.get_lb() - b.get_lb()) <= COUENNE_HACKED_EPS)
          if (fabs(a.get_ub() - b.get_ub()) <= COUENNE_HACKED_EPS)
            return 1;
  return 0;
}

// simple nauty definitely not thread safe
static int calls = 0;
static int maxLevel = 0;
static void
userlevelproc(int *lab, int *ptn, int level, int *orbits, statsblk *stats,
  int tv, int index, int tcellsize,
  int numcells, int childcount, int n)
{
  calls++;
  if (level > maxLevel) {
    printf("Level %d after %d calls\n", level, calls);
    fprintf(stderr, "Level %d after %d calls\n", level, calls);
    maxLevel = level;
  }
  if (level > 1500) {
    throw CoinError("May take too long", "", "CbcSymmetry");
  }
  //}
  return;
}
void CbcSymmetry::Compute_Symmetry() const
{

  nauty_info_->options()->userlevelproc = userlevelproc;
  std::sort(node_info_.begin(), node_info_.end(), node_sort);

  for (std::vector< Node >::iterator i = node_info_.begin(); i != node_info_.end(); ++i)
    (*i).color_vertex(-1);

  int color = 1;
  for (std::vector< Node >::iterator i = node_info_.begin(); i != node_info_.end(); ++i) {
    if ((*i).get_color() == -1) {
      (*i).color_vertex(color);
#ifdef PRINT_MORE
      printf("Graph vertex %d is given color %d\n", (*i).get_index(), color);
#endif
      nauty_info_->color_node((*i).get_index(), color);
      for (std::vector< Node >::iterator j = i + 1; j != node_info_.end(); ++j)
        if (compare((*i), (*j)) == 1) {
          (*j).color_vertex(color);
          nauty_info_->color_node((*j).get_index(), color);
#ifdef PRINT_MORE
          printf("Graph vertex %d is given color %d, the same as vertex %d\n", (*j).get_index(), color, (*i).get_index());
#endif
        }
      //       else
      // j = node_info_. end();
      color++;
    }
  }

  //Print_Orbits ();
  nauty_info_->computeAuto();
  //Print_Orbits ();
}

int CbcSymmetry::statsOrbits(CbcModel *model, int type) const
{
  char general[200];
  int returnCode = 0;
  bool printSomething = true;
  if (type) {
    double branchSuccess = 0.0;
    if (nautyBranchSucceeded_)
      branchSuccess = nautyOtherBranches_ / nautyBranchSucceeded_;
    double fixSuccess = 0.0;
    if (nautyFixSucceeded_)
      fixSuccess = nautyFixes_ / nautyFixSucceeded_;
    if (nautyBranchCalls_ > lastNautyBranchCalls_ || nautyFixCalls_ > lastNautyFixCalls_) {
      sprintf(general, "Orbital branching tried %d times, succeeded %d times - average extra %7.3f, fixing %d times (%d, %7.3f)",
        nautyBranchCalls_, nautyBranchSucceeded_, branchSuccess,
        nautyFixCalls_, nautyFixSucceeded_, fixSuccess);
      lastNautyBranchCalls_ = nautyBranchCalls_;
      lastNautyFixCalls_ = nautyFixCalls_;
    } else {
      printSomething = false;
    }
  } else {
    returnCode = nauty_info_->getNumGenerators();
    if (!nauty_info_->errorStatus()) {
      if (returnCode && numberUsefulOrbits_) {
        sprintf(general, "Nauty: %d orbits (%d useful covering %d variables), %d generators, group size: %g - dense size %d, sparse %d - took %g seconds",
          nauty_info_->getNumOrbits(), numberUsefulOrbits_, numberUsefulObjects_,
          nauty_info_->getNumGenerators(),
          nauty_info_->getGroupSize(),
          whichOrbit_[0], whichOrbit_[1], nautyTime_);
      } else {
        if ((model->moreSpecialOptions2() & (128 | 256)) != (128 | 256))
          sprintf(general, "Nauty did not find any useful orbits in time %g", nautyTime_);
        else
          sprintf(general, "Nauty did not find any useful orbits - but keeping Nauty on");
      }
    } else {
      // error
      sprintf(general, "Nauty failed with error code %d (%g seconds)",
        nauty_info_->errorStatus(), nautyTime_);
      model->setMoreSpecialOptions2(model->moreSpecialOptions2() & ~(128 | 256));
    }
  }
  if (printSomething)
    model->messageHandler()->message(CBC_GENERAL,
      model->messages())
      << general << CoinMessageEol;
  return returnCode;
}

void CbcSymmetry::Print_Orbits() const
{

  //printf ("num gens = %d, num orbits = %d \n", nauty_info_ -> getNumGenerators(), nauty_info_ -> getNumOrbits() );

  std::vector< std::vector< int > > *new_orbits = nauty_info_->getOrbits();

  printf("Nauty: %d generators, group size: %.0g",
    //  nauty_info_->getNumOrbits(),
    nauty_info_->getNumGenerators(),
    nauty_info_->getGroupSize());

  int nNonTrivialOrbits = 0;

  for (unsigned int i = 0; i < new_orbits->size(); i++) {
    if ((*new_orbits)[i].size() > 1)
      nNonTrivialOrbits++;
    else
      continue;

    // int orbsize = (*new_orbits)[i].size();
    // printf( "Orbit %d [size: %d] [", i, orbsize);
    // copy ((*new_orbits)[i].begin(), (*new_orbits)[i].end(),
    // 	  std::ostream_iterator<int>(std::cout, " "));
    // printf("] \n");
  }

  printf(" (%d non-trivial orbits).\n", nNonTrivialOrbits);

#if 1
  if (nNonTrivialOrbits) {

    int orbCnt = 0;

    std::vector< std::vector< int > > *orbits = nauty_info_->getOrbits();

    for (std::vector< std::vector< int > >::iterator i = orbits->begin(); i != orbits->end(); ++i) {

      printf("Orbit %d: ", orbCnt++);

      for (std::vector< int >::iterator j = i->begin(); j != i->end(); ++j)
        printf(" %d", *j);

      printf("\n");
    }
  }
#endif

#if 0
  if (nNonTrivialOrbits)
    for (int i=0; i< nVars (); i++) {

      std::vector< int > *branch_orbit = Find_Orbit (i);

      if (branch_orbit -> size () > 1) {
  	printf ("x%04d: ", i);

  	for (std::vector<int>::iterator it = branch_orbit -> begin (); it != branch_orbit -> end (); ++it) 
  	  printf ("%d ", *it);
  	printf ("\n");
      }
    }
#endif
  delete new_orbits;
}
void CbcSymmetry::fillOrbits()
{
  for (int i = 0; i < numberColumns_; i++)
    whichOrbit_[i] = -1;
  numberUsefulOrbits_ = 0;
  numberUsefulObjects_ = 0;

  std::vector< std::vector< int > > *orbits = nauty_info_->getOrbits();

  for (std::vector< std::vector< int > >::iterator i = orbits->begin(); i != orbits->end(); ++i) {
    int nUseful = 0;
    int jColumn = -2;
    for (std::vector< int >::iterator j = i->begin(); j != i->end(); ++j) {
      int iColumn = *j;
      if (iColumn < numberColumns_) {
        whichOrbit_[iColumn] = numberUsefulOrbits_;
        nUseful++;
        jColumn = iColumn;
      }
    }
    if (nUseful > 1) {
      numberUsefulOrbits_++;
      numberUsefulObjects_ += nUseful;
    } else if (jColumn >= 0) {
      assert(nUseful);
      whichOrbit_[jColumn] = -2;
    }
  }
  delete orbits;
}
int CbcSymmetry::largestOrbit(const double *lower, const double *upper) const
{
  int *counts = new int[numberUsefulOrbits_];
  memset(counts, 0, numberUsefulOrbits_ * sizeof(int));
  for (int i = 0; i < numberColumns_; i++) {
    int iOrbit = whichOrbit_[i];
    if (iOrbit >= 0) {
      if (lower[i] == 0.0 && upper[i] == 1.0)
        counts[iOrbit]++;
    }
  }
  int iOrbit = -1;
  int maxOrbit = 0;
  for (int i = 0; i < numberUsefulOrbits_; i++) {
    if (counts[i] > maxOrbit) {
      maxOrbit = counts[i];
      iOrbit = i;
    }
  }
  delete[] counts;
  return iOrbit;
}

std::vector< int > *CbcSymmetry::Find_Orbit(int index) const
{

  std::vector< int > *orbit = new std::vector< int >;
  int which_orbit = -1;
  std::vector< std::vector< int > > *new_orbits = nauty_info_->getOrbits();

  for (unsigned int i = 0; i < new_orbits->size(); i++) {
    for (unsigned int j = 0; j < (*new_orbits)[i].size(); j++) {
      //   for (std::vector <int>:: iterator j = new_orbits[i].begin(); new_orbits[i].end(); ++j){
      if ((*new_orbits)[i][j] == index)
        which_orbit = i;
    }
  }

  //  for (std::vector <int>:: iterator j = new_orbits[which_orbit].begin(); new_orbits[which_orbit].end(), ++j)
  for (unsigned int j = 0; j < (*new_orbits)[which_orbit].size(); j++)
    orbit->push_back((*new_orbits)[which_orbit][j]);

  delete new_orbits;

  return orbit;
}

void CbcSymmetry::ChangeBounds(const double *new_lb, const double *new_ub,
  int num_cols, bool justFixedAtOne) const
{
  if (justFixedAtOne)
    nautyFixCalls_++;
  else
    nautyBranchCalls_++;
  std::sort(node_info_.begin(), node_info_.end(), index_sort);

  for (int i = 0; i < num_cols; i++) {
    //   printf("Var %d  lower bound: %f   upper bound %f \n", i, new_lb[i], new_ub[i]);

    assert(node_info_[i].get_index() == i);
    double newLower = new_lb[i];
    double newUpper = new_ub[i];
    if (justFixedAtOne) {
      // free up all fixed at zero
      if (!newLower)
        newUpper = 1.0;
    }
    node_info_[i].bounds(newLower, newUpper);
    //printf("Var %d  INPUT lower bound: %f   upper bound %f \n", i, node_info_[i].get_lb(), node_info_[i].get_ub());
  }
}
void CbcSymmetry::setupSymmetry(CbcModel * model)
{
  OsiSolverInterface * solver = model->continuousSolver();
  double startCPU = CoinCpuTime();
  const double *objective = solver->getObjCoefficients();
  const double *columnLower = solver->getColLower();
  const double *columnUpper = solver->getColUpper();
  int numberColumns = solver->getNumCols();
  int numberRows = solver->getNumRows();
  int iRow, iColumn;

  // Row copy
  CoinPackedMatrix matrixByRow(*solver->getMatrixByRow());
  const double *elementByRow = matrixByRow.getElements();
  const int *column = matrixByRow.getIndices();
  const CoinBigIndex *rowStart = matrixByRow.getVectorStarts();
  const int *rowLength = matrixByRow.getVectorLengths();

  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();
  //  // Find Coefficients

  /// initialize nauty

  int num_affine = 0;

  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    if (objective[iColumn] && objective[iColumn] != 1.0)
      num_affine++;
  }
  for (iRow = 0; iRow < numberRows; iRow++) {
    for (CoinBigIndex j = rowStart[iRow];
         j < rowStart[iRow] + rowLength[iRow]; j++) {
      int jColumn = column[j];
      double value = elementByRow[j];
      if (value != 1.0)
        num_affine++;
    }
  }

  // Create Nauty object

  int coef_count = numberRows + numberColumns + 1;
  int nc = num_affine + coef_count;
  // create graph (part 1)

  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    Node var_vertex;
    var_vertex.node(iColumn, 0.0, columnLower[iColumn], columnUpper[iColumn], -1, -1);
    node_info_.push_back(var_vertex);
  }
  // objective
  int index = numberColumns;
  {
    Node vertex;
    vertex.node(index, 0.0, -COIN_DBL_MAX, COIN_DBL_MAX,
      COUENNE_HACKED_EXPRGROUP, 0);
    node_info_.push_back(vertex);
  }

  // compute space for sparse
  size_t *v = NULL;
  int *d = NULL;
  int *e = NULL;
  bool sparse = false;
  double spaceDense = nc + WORDSIZE - 1;
  spaceDense *= nc + WORDSIZE - 1;
  spaceDense /= WORDSIZE;
  int spaceSparse = 0;
  {
    size_t numberElements = 0;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value = objective[iColumn];
      if (value) {
        if (value == 1.0) {
          numberElements += 2;
        } else {
          numberElements += 4;
          coef_count++;
        }
      }
    }
    for (iRow = 0; iRow < numberRows; iRow++) {
      for (CoinBigIndex j = rowStart[iRow];
           j < rowStart[iRow] + rowLength[iRow]; j++) {
        int jColumn = column[j];
        double value = elementByRow[j];
        if (value == 1.0) {
          numberElements += 2;
        } else {
          numberElements += 4;
          coef_count++;
        }
      }
    }
    spaceSparse = 2 * nc + numberElements;
    //printf("Space for sparse is %d for dense %g\n",
    //	   spaceSparse,spaceDense);
#ifdef NTY_TRACES
    bool goSparse = true;
#else
    bool goSparse = (spaceSparse < spaceDense);
#endif
    // for now always sparse
    goSparse = true;
    if (goSparse) {
      sparse = true;
      v = new size_t[nc + 1];
      d = new int[nc];
      e = new int[numberElements];
      size_t *counts = new size_t[coef_count + 1];
      memset(counts, 0, coef_count * sizeof(size_t));
      coef_count = numberRows + numberColumns + 1;
      // create graph (part 2)
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        double value = objective[iColumn];
        if (value) {
          if (value == 1.0) {
            counts[index]++;
            counts[iColumn]++;
          } else {
            counts[index]++;
            counts[coef_count] += 2;
            counts[iColumn]++;
            coef_count++;
          }
        }
      }
      index++;
      for (iRow = 0; iRow < numberRows; iRow++) {
        for (CoinBigIndex j = rowStart[iRow];
             j < rowStart[iRow] + rowLength[iRow]; j++) {
          int jColumn = column[j];
          double value = elementByRow[j];
          if (value == 1.0) {
            counts[index]++;
            counts[jColumn]++;
          } else {
            counts[index]++;
            counts[coef_count] += 2;
            counts[jColumn]++;
            coef_count++;
          }
        }
        index++;
      }
      // create graph (part 3)
      assert(nc == coef_count);
      numberElements = 0;
      v[0] = 0;
      for (int i = 0; i < nc; i++) {
        int count = counts[i];
        d[i] = count;
        counts[i] = v[i];
        numberElements += count;
        ;
        v[i + 1] = numberElements;
      }
      index = numberColumns;
      coef_count = numberRows + numberColumns + 1;
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        double value = objective[iColumn];
        if (value) {
          int where;
          if (value == 1.0) {
            where = counts[index];
            counts[index]++;
            e[where] = iColumn;
            where = counts[iColumn];
            counts[iColumn]++;
            e[where] = index;
          } else {
            Node coef_vertex;
            coef_vertex.node(coef_count, value, value, value, -2, 0);
            node_info_.push_back(coef_vertex);
            where = counts[index];
            counts[index]++;
            e[where] = coef_count;
            where = counts[coef_count];
            counts[coef_count]++;
            e[where] = index;
            where = counts[coef_count];
            counts[coef_count]++;
            e[where] = iColumn;
            where = counts[iColumn];
            counts[iColumn]++;
            e[where] = coef_count;
            coef_count++;
          }
        }
      }
      index++;
      for (iRow = 0; iRow < numberRows; iRow++) {
        Node vertex;
        vertex.node(index, 0.0, rowLower[iRow], rowUpper[iRow],
          COUENNE_HACKED_EXPRGROUP, 0);
        node_info_.push_back(vertex);
        for (CoinBigIndex j = rowStart[iRow];
             j < rowStart[iRow] + rowLength[iRow]; j++) {
          int jColumn = column[j];
          double value = elementByRow[j];
          int where;
          if (value == 1.0) {
            where = counts[index];
            counts[index]++;
            e[where] = jColumn;
            where = counts[jColumn];
            counts[jColumn]++;
            e[where] = index;
          } else {
            Node coef_vertex;
            coef_vertex.node(coef_count, value, value, value, -2, 0);
            node_info_.push_back(coef_vertex);
            where = counts[index];
            counts[index]++;
            e[where] = coef_count;
            where = counts[coef_count];
            counts[coef_count]++;
            e[where] = index;
            where = counts[coef_count];
            counts[coef_count]++;
            e[where] = jColumn;
            where = counts[jColumn];
            counts[jColumn]++;
            e[where] = coef_count;
            coef_count++;
          }
        }
        index++;
      }
      delete[] counts;
    }
  }

  nauty_info_ = new CbcNauty(nc, v, d, e);
  delete[] v;
  delete[] d;
  delete[] e;
  if (!sparse) {
    // create graph (part 2)
    coef_count = numberRows + numberColumns + 1;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value = objective[iColumn];
      if (value) {
        if (value == 1.0) {
          nauty_info_->addElement(index, iColumn);
          nauty_info_->addElement(iColumn, index);
        } else {
          Node coef_vertex;
          coef_vertex.node(coef_count, value, value, value, -2, 0);
          node_info_.push_back(coef_vertex);
          nauty_info_->addElement(index, coef_count);
          nauty_info_->addElement(coef_count, index);
          nauty_info_->addElement(coef_count, iColumn);
          nauty_info_->addElement(iColumn, coef_count);
          coef_count++;
        }
      }
    }
    index++;
    for (iRow = 0; iRow < numberRows; iRow++) {
      Node vertex;
      vertex.node(index, 0.0, rowLower[iRow], rowUpper[iRow],
        COUENNE_HACKED_EXPRGROUP, 0);
      node_info_.push_back(vertex);
      for (CoinBigIndex j = rowStart[iRow];
           j < rowStart[iRow] + rowLength[iRow]; j++) {
        int jColumn = column[j];
        double value = elementByRow[j];
        if (value == 1.0) {
          nauty_info_->addElement(index, jColumn);
          nauty_info_->addElement(jColumn, index);
        } else {
          Node coef_vertex;
          coef_vertex.node(coef_count, value, value, value, -2, 0);
          node_info_.push_back(coef_vertex);
          nauty_info_->addElement(index, coef_count);
          nauty_info_->addElement(coef_count, index);
          nauty_info_->addElement(coef_count, jColumn);
          nauty_info_->addElement(jColumn, coef_count);
          coef_count++;
        }
      }
      index++;
    }
  }
  numberColumns_ = numberColumns;
  whichOrbit_ = new int[2 * numberColumns_];
  nautyBranchCalls_ = 0;
  nautyBranchSucceeded_ = 0;
  nautyFixCalls_ = 0;
  nautyFixSucceeded_ = 0;
  nautyTime_ = 0.0;
  nautyFixes_ = 0.0;
  nautyOtherBranches_ = 0.0;
  try {
    Compute_Symmetry();
  } catch (CoinError &e) {
    char general[200];
    sprintf(general, "Nauty - initial level %d - will probably take too long",
        maxLevel);
    model->messageHandler()->message(CBC_GENERAL,model->messages())
      <<general <<CoinMessageEol;
  }
  fillOrbits();
  //whichOrbit_[2]=numberUsefulOrbits_;
  //Print_Orbits ();
  // stats in array
  if (spaceDense < COIN_INT_MAX)
    whichOrbit_[0] = spaceDense;
  else
    whichOrbit_[0] = COIN_INT_MAX;
  whichOrbit_[1] = spaceSparse;
  double endCPU = CoinCpuTime();
  nautyTime_ = endCPU - startCPU;
}
// Fixes variables using orbits (returns number fixed)
int CbcSymmetry::orbitalFixing(OsiSolverInterface *solver)
{
  int numberColumns = solver->getNumCols();
  char *status = new char[numberColumns];
  ChangeBounds(solver->getColLower(),
    solver->getColUpper(),
    solver->getNumCols(), true);
  Compute_Symmetry();
  fillOrbits();
  int n = 0;
  //#define PRINT_MORE 1
  const int *alternativeOrbits = whichOrbit();
  if (alternativeOrbits) {
    for (int i = 0; i < numberColumns; i++) {
      char type = '0';
      if (solver->getColUpper()[i]) {
        if (solver->getColLower()[i]) {
          type = '1';
        } else {
          double value = solver->getColSolution()[i];
          if (value < 0.0001)
            type = 'L';
          else if (value > 0.9999)
            type = 'U';
          else
            type = 'X';
        }
      }
      status[i] = type;
    }
    n = 0;
    for (int i = 0; i < numberColumns; i++) {
      if (status[i] != '0' && status[i] != '1') {
        int iOrbit = alternativeOrbits[i];
        if (iOrbit < 0)
          continue;
        for (int j = i + 1; j < numberColumns; j++) {
          if (status[j] == '0' && alternativeOrbits[j] == iOrbit) {
#if PRINT_MORE > 1
            printf("In alternative orbit %d - %d free (%c), %d fixed to 0\n",
              iOrbit, i, status[i], j);
#endif
            status[i] = '0'; // can fix on both branches
            solver->setColUpper(i, 0.0);
            n++;
            break;
          }
        }
      }
    }
  }
  delete[] status;
  if (n) {
    nautyFixSucceeded_++;
    nautyFixes_ += n;
#if PRINT_MORE
    printf("%d orbital fixes\n", n);
#endif
  }
  return n;
}
// Default Constructor
CbcSymmetry::CbcSymmetry()
  : nauty_info_(NULL)
  , numberColumns_(0)
  , numberUsefulOrbits_(0)
  , numberUsefulObjects_(0)
  , whichOrbit_(NULL)
{
}
// Copy constructor
CbcSymmetry::CbcSymmetry(const CbcSymmetry &rhs)
{
  node_info_ = rhs.node_info_;
  nauty_info_ = new CbcNauty(*rhs.nauty_info_);
  numberUsefulOrbits_ = rhs.numberUsefulOrbits_;
  numberUsefulObjects_ = rhs.numberUsefulObjects_;
  numberColumns_ = rhs.numberColumns_;
  if (rhs.whichOrbit_)
    whichOrbit_ = CoinCopyOfArray(rhs.whichOrbit_, numberColumns_);
  else
    whichOrbit_ = NULL;
}

// Assignment operator
CbcSymmetry &
CbcSymmetry::operator=(const CbcSymmetry &rhs)
{
  if (this != &rhs) {
    delete nauty_info_;
    node_info_ = rhs.node_info_;
    nauty_info_ = new CbcNauty(*rhs.nauty_info_);
    delete[] whichOrbit_;
    numberColumns_ = rhs.numberColumns_;
    numberUsefulOrbits_ = rhs.numberUsefulOrbits_;
    numberUsefulObjects_ = rhs.numberUsefulObjects_;
    if (rhs.whichOrbit_)
      whichOrbit_ = CoinCopyOfArray(rhs.whichOrbit_, numberColumns_);
    else
      whichOrbit_ = NULL;
  }
  return *this;
}

// Destructor
CbcSymmetry::~CbcSymmetry()
{
  delete nauty_info_;
  delete[] whichOrbit_;
}

CbcNauty::CbcNauty(int vertices, const size_t *v, const int *d, const int *e)
{
  //printf("Need sparse nauty - wordsize %d\n",WORDSIZE);
  n_ = vertices;
  m_ = (n_ + WORDSIZE - 1) / WORDSIZE;
  if (v)
    nel_ = v[n_];
  else
    nel_ = 0;

  //printf ("size of long = %d (%d)\nwordsize = %d\nn,m = %d,%d\n",
  //          SIZEOF_LONG, sizeof (long), WORDSIZE, n_, m_);

  nauty_check(WORDSIZE, m_, n_, NAUTYVERSIONID);

  /// Apparently sizes are skewed on 64bit machines

#define MULTIPLIER 1

  if (!nel_) {
    G_ = (graph *)malloc(MULTIPLIER * m_ * n_ * sizeof(int));
    GSparse_ = NULL;
  } else {
    G_ = NULL;
    GSparse_ = (sparsegraph *)malloc(sizeof(sparsegraph));
    SG_INIT(*GSparse_);
    SG_ALLOC(*GSparse_, n_, nel_, "malloc");
    GSparse_->nv = n_; /* Number of vertices */
    GSparse_->nde = nel_;
  }
  lab_ = (int *)malloc(MULTIPLIER * n_ * sizeof(int));
  ptn_ = (int *)malloc(MULTIPLIER * n_ * sizeof(int));
  active_ = NULL;
  orbits_ = (int *)malloc(MULTIPLIER * n_ * sizeof(int));
#ifndef NTY_TRACES
  options_ = (optionblk *)malloc(MULTIPLIER * sizeof(optionblk));
  stats_ = (statsblk *)malloc(MULTIPLIER * sizeof(statsblk));
#else
  options_ = (TracesOptions *)malloc(MULTIPLIER * sizeof(TracesOptions));
  stats_ = (TracesStats *)malloc(MULTIPLIER * sizeof(TracesStats));
#endif
  worksize_ = 100 * m_;
  workspace_ = (setword *)malloc(MULTIPLIER * worksize_ * sizeof(setword));
  canonG_ = NULL;
  if ((G_ == 0 && GSparse_ == 0) || lab_ == 0 || ptn_ == 0 || orbits_ == 0 || options_ == 0 || stats_ == 0 || workspace_ == 0)
    assert(0);

  // Zero allocated memory
  if (G_) {
    memset(G_, 0, m_ * n_ * sizeof(int));
  } else {
    //for (int i=0;i<n_;i++) {
    //GSparse_->v[i]=v[i];
    //}
    memcpy(GSparse_->v, v, n_ * sizeof(size_t));
    memcpy(GSparse_->d, d, n_ * sizeof(int));
    memcpy(GSparse_->e, e, nel_ * sizeof(int));
  }
  memset(lab_, 0, n_ * sizeof(int));
  memset(ptn_, 0, n_ * sizeof(int));
  memset(orbits_, 0, n_ * sizeof(int));
  memset(workspace_, 0, worksize_ * sizeof(setword));
#ifndef NTY_TRACES
  memset(options_, 0, MULTIPLIER * sizeof(optionblk));
#else
  memset(options_, 0, MULTIPLIER * sizeof(TracesOptions));
#endif

  // Set the options you want
#ifndef NTY_TRACES
  options_->getcanon = FALSE;
  options_->digraph = FALSE;
  options_->writeautoms = FALSE;
  options_->writemarkers = FALSE;
  options_->defaultptn = TRUE;
  options_->cartesian = FALSE;
  options_->linelength = 78;
  options_->outfile = NULL;
  options_->userrefproc = NULL;
  options_->userautomproc = NULL;
  options_->userlevelproc = NULL;
  options_->usernodeproc = NULL;
  //  options_->usertcellproc = NULL;
  options_->invarproc = NULL;
  options_->tc_level = 100;
  options_->mininvarlevel = 0;
  options_->maxinvarlevel = 1;
  options_->invararg = 0;
  options_->dispatch = &dispatch_graph;
#else
  options_->getcanon = FALSE;
  options_->writeautoms = FALSE;
  options_->cartesian = FALSE;
  options_->digraph = FALSE;
  options_->defaultptn = TRUE;
  options_->linelength = 78;
#endif
  if (G_) {
    // Make an empty graph
    for (int j = 0; j < n_; j++) {
      set *gv = GRAPHROW(G_, j, m_);
      EMPTYSET(gv, m_);
    }
  }

  vstat_ = new int[n_];
  clearPartitions();
  afp_ = NULL;
}

CbcNauty::~CbcNauty()
{
  if (G_)
    free(G_);
  if (GSparse_) {
    SG_FREE(*GSparse_);
    free(GSparse_);
  }
  if (lab_)
    free(lab_);
  if (ptn_)
    free(ptn_);
  if (active_)
    free(active_);
  if (orbits_)
    free(orbits_);
  if (options_)
    free(options_);
  if (stats_)
    free(stats_);
  if (workspace_)
    free(workspace_);
  if (canonG_)
    free(canonG_);
  if (vstat_)
    delete[] vstat_;
}
// Copy constructor
CbcNauty::CbcNauty(const CbcNauty &rhs)
{
  n_ = rhs.n_;
  m_ = rhs.m_;
  nel_ = rhs.nel_;
  G_ = NULL;
  GSparse_ = NULL;
  if (!nel_) {
    G_ = (graph *)malloc(MULTIPLIER * m_ * n_ * sizeof(int));
  } else {
    GSparse_ = (sparsegraph *)malloc(sizeof(sparsegraph));
    SG_INIT(*GSparse_);
    SG_ALLOC(*GSparse_, n_, nel_, "malloc");
    GSparse_->nv = n_; /* Number of vertices */
    GSparse_->nde = nel_;
  }
  lab_ = (int *)malloc(MULTIPLIER * n_ * sizeof(int));
  ptn_ = (int *)malloc(MULTIPLIER * n_ * sizeof(int));
  orbits_ = (int *)malloc(MULTIPLIER * n_ * sizeof(int));
#ifndef NTY_TRACES
  options_ = (optionblk *)malloc(MULTIPLIER * sizeof(optionblk));
  stats_ = (statsblk *)malloc(MULTIPLIER * sizeof(statsblk));
#else
  options_ = (TracesOptions *)malloc(MULTIPLIER * sizeof(TracesOptions));
  stats_ = (TracesStats *)malloc(MULTIPLIER * sizeof(TracesStats));
#endif
  worksize_ = 100 * m_;
  workspace_ = (setword *)malloc(MULTIPLIER * worksize_ * sizeof(setword));
  vstat_ = new int[n_];
  canonG_ = NULL;
  if ((G_ == 0 && GSparse_ == 0) || lab_ == 0 || ptn_ == 0 || orbits_ == 0 || options_ == 0 || stats_ == 0 || workspace_ == 0)
    assert(0);

  // Copy allocated memory
  if (G_) {
    memcpy(G_, rhs.G_, m_ * n_ * sizeof(int));
  } else {
    memcpy(GSparse_->v, rhs.GSparse_->v, n_ * sizeof(size_t));
    memcpy(GSparse_->d, rhs.GSparse_->d, n_ * sizeof(int));
    memcpy(GSparse_->e, rhs.GSparse_->e, nel_ * sizeof(int));
  }
  memcpy(lab_, rhs.lab_, n_ * sizeof(int));
  memcpy(ptn_, rhs.ptn_, n_ * sizeof(int));
  memcpy(orbits_, rhs.orbits_, n_ * sizeof(int));
  memcpy(workspace_, rhs.workspace_, worksize_ * sizeof(setword));
#ifndef NTY_TRACES
  memcpy(options_, rhs.options_, MULTIPLIER * sizeof(optionblk));
  memcpy(stats_, rhs.stats_, MULTIPLIER * sizeof(statsblk));
#else
  memcpy(options_, rhs.options_, MULTIPLIER * sizeof(TracesOptions));
  memcpy(stats_, rhs.stats_, MULTIPLIER * sizeof(TracesStats));
#endif
  memcpy(vstat_, rhs.vstat_, n_ * sizeof(int));

  // ? clearPartitions();
  active_ = NULL;
  afp_ = rhs.afp_; // ? no copy ?
}

// Assignment operator
CbcNauty &
CbcNauty::operator=(const CbcNauty &rhs)
{
  if (this != &rhs) {
    if (G_)
      free(G_);
    if (GSparse_) {
      SG_FREE(*GSparse_);
      free(GSparse_);
    }
    if (lab_)
      free(lab_);
    if (ptn_)
      free(ptn_);
    if (active_)
      free(active_);
    if (orbits_)
      free(orbits_);
    if (options_)
      free(options_);
    if (stats_)
      free(stats_);
    if (workspace_)
      free(workspace_);
    if (canonG_)
      free(canonG_);
    if (vstat_)
      delete[] vstat_;
    {
      n_ = rhs.n_;
      m_ = rhs.m_;
      nel_ = rhs.nel_;
      G_ = NULL;
      GSparse_ = NULL;
      if (!nel_) {
        G_ = (graph *)malloc(MULTIPLIER * m_ * n_ * sizeof(int));
      } else {
        GSparse_ = (sparsegraph *)malloc(sizeof(sparsegraph));
        SG_INIT(*GSparse_);
        SG_ALLOC(*GSparse_, n_, nel_, "malloc");
        GSparse_->nv = n_; /* Number of vertices */
        GSparse_->nde = nel_;
      }
      lab_ = (int *)malloc(MULTIPLIER * n_ * sizeof(int));
      ptn_ = (int *)malloc(MULTIPLIER * n_ * sizeof(int));
      orbits_ = (int *)malloc(MULTIPLIER * n_ * sizeof(int));
#ifndef NTY_TRACES
      options_ = (optionblk *)malloc(MULTIPLIER * sizeof(optionblk));
      stats_ = (statsblk *)malloc(MULTIPLIER * sizeof(statsblk));
#else
      options_ = (TracesOptions *)malloc(MULTIPLIER * sizeof(TracesOptions));
      stats_ = (TracesStats *)malloc(MULTIPLIER * sizeof(TracesStats));
#endif
      worksize_ = 100 * m_;
      workspace_ = (setword *)malloc(MULTIPLIER * worksize_ * sizeof(setword));
      vstat_ = new int[n_];
      canonG_ = NULL;
      if ((G_ == 0 && GSparse_ == 0) || lab_ == 0 || ptn_ == 0 || orbits_ == 0 || options_ == 0 || stats_ == 0 || workspace_ == 0)
        assert(0);

      // Copy allocated memory
      if (!nel_) {
        memcpy(G_, rhs.G_, m_ * n_ * sizeof(int));
      } else {
        memcpy(GSparse_->v, rhs.GSparse_->v, n_ * sizeof(size_t));
        memcpy(GSparse_->d, rhs.GSparse_->d, n_ * sizeof(int));
        memcpy(GSparse_->e, rhs.GSparse_->e, nel_ * sizeof(int));
      }
      memcpy(lab_, rhs.lab_, n_ * sizeof(int));
      memcpy(ptn_, rhs.ptn_, n_ * sizeof(int));
      memcpy(orbits_, rhs.orbits_, n_ * sizeof(int));
      memcpy(workspace_, rhs.workspace_, worksize_ * sizeof(setword));
#ifndef NTY_TRACES
      memcpy(options_, rhs.options_, MULTIPLIER * sizeof(optionblk));
      memcpy(stats_, rhs.stats_, MULTIPLIER * sizeof(statsblk));
#else
      memcpy(options_, rhs.options_, MULTIPLIER * sizeof(TracesOptions));
      memcpy(stats_, rhs.stats_, MULTIPLIER * sizeof(TracesStats));
#endif
      memcpy(vstat_, rhs.vstat_, n_ * sizeof(int));

      // ? clearPartitions();
      active_ = NULL;
      afp_ = rhs.afp_; // ? no copy ?
    }
  }
  return *this;
}

void CbcNauty::addElement(int ix, int jx)
{
  // Right now die if bad index.  Can throw exception later
  //printf("addelement %d %d \n", ix, jx);
  assert(ix < n_ && jx < n_);
  if (ix != jx) { //No Loops
    set *gv = GRAPHROW(G_, ix, m_);
    ADDELEMENT(gv, jx);
    set *gv2 = GRAPHROW(G_, jx, m_);
    ADDELEMENT(gv2, ix);
    autoComputed_ = false;
  }
}

void CbcNauty::clearPartitions()
{
  for (int j = 0; j < n_; j++) {
    vstat_[j] = 1;
    //printf("vstat %d = %d", j, vstat_[j]);
  }

  autoComputed_ = false;
}

void CbcNauty::computeAuto()
{

  //  if (autoComputed_) return;

  //double startCPU = CoinCpuTime ();

  options_->defaultptn = FALSE;

  // Here we only implement the partitions
  // [ fix1 | fix0 (union) free | constraints ]
  int ix = 0;

  for (int color = 1; color <= n_; color++) {
    for (int j = 0; j < n_; j++) {
      if (vstat_[j] == color) {
        lab_[ix] = j;
        ptn_[ix] = color;
        ix++;
      }
    }
    if (ix > 0)
      ptn_[ix - 1] = 0;
  }

  /*
  for (int j = 0; j < n_; j++)
    printf("ptn %d = %d      lab = %d \n", j, ptn_[j], lab_[j]);
  */

  // Should be number of columns
  assert(ix == n_);
  // Now the constraints if needed

  // Compute Partition

  if (G_) {
#ifndef NTY_TRACES
    nauty(G_, lab_, ptn_, active_, orbits_, options_,
      stats_, workspace_, worksize_, m_, n_, canonG_);
#else
    abort();
#endif
  } else {
#if NAUTY_MAX_LEVEL
    nauty_maxalllevel = NAUTY_MAX_LEVEL;
#endif
#ifndef NTY_TRACES
    options_->dispatch = &dispatch_sparse;
    sparsenauty(GSparse_, lab_, ptn_, orbits_, options_,
      stats_, NULL);
#else
    //options_->dispatch = &dispatch_sparse;
    Traces(GSparse_, lab_, ptn_, orbits_, options_,
      stats_, NULL);
#endif
  }
  autoComputed_ = true;

  //double endCPU = CoinCpuTime ();

  //nautyTime_ += endCPU - startCPU;
  // Need to make sure all generators are written
  if (afp_)
    fflush(afp_);
}

void CbcNauty::deleteElement(int ix, int jx)
{
  // Right now die if bad index.  Can throw exception later
  assert(ix < n_ && jx < n_);
  set *gv = GRAPHROW(G_, ix, m_);
  if (ISELEMENT(gv, jx)) {
    DELELEMENT(gv, jx);
  }
  autoComputed_ = false;
}

double
CbcNauty::getGroupSize() const
{
  if (!autoComputed_)
    return -1.0;
  return (stats_->grpsize1 * pow(10.0, (double)stats_->grpsize2));
}

int CbcNauty::getNumGenerators() const
{
  if (!autoComputed_)
    return -1;
  return (stats_->numgenerators);
}

int CbcNauty::getNumOrbits() const
{
  if (!autoComputed_)
    return -1;
  return (stats_->numorbits);
}

std::vector< std::vector< int > >
  *CbcNauty::getOrbits() const
{
  std::vector< std::vector< int > > *orb = new std::vector< std::vector< int > >;
  if (!autoComputed_)
    return orb;
  orb->resize(getNumOrbits());
  std::multimap< int, int > orbmap;
  std::set< int > orbkeys;
  for (int j = 0; j < n_; j++) {
    orbkeys.insert(orbits_[j]);
    orbmap.insert(std::make_pair(orbits_[j], j));
  }

  int orbix = 0;
  for (std::set< int >::iterator it = orbkeys.begin();
       it != orbkeys.end(); ++it) {
    std::multimap< int, int >::iterator pos;
    for (pos = orbmap.lower_bound(*it);
         pos != orbmap.upper_bound(*it); ++pos) {
      (*orb)[orbix].push_back(pos->second);
    }
    orbix++;
  }

  assert(orbix == getNumOrbits());
  return orb;
}

void CbcNauty::getVstat(double *v, int nv)
{
  assert(nv == n_);
  memcpy(v, vstat_, nv * sizeof(VarStatus));
}

/*
bool
CbcNauty::isAllFixOneOrbit(const std::vector<int> &orbit) const
{

  for(std::vector<int>::const_iterator it = orbit.begin();
      it != orbit.end(); ++it) {
    if (*it >= n_) return false;
    if (vstat_[*it] != FIX_AT_ONE) return false;
  }
  return true;
}

bool
CbcNauty::isAllFreeOrbit(const std::vector<int> &orbit) const
{
  for(std::vector<int>::const_iterator it = orbit.begin();
      it != orbit.end(); ++it) {
    if (*it >= n_) return false;
    if (vstat_[*it] != FREE) return false;
  }
  return true;
}

bool
CbcNauty::isConstraintOrbit(const std::vector<int> &orbit) const
{
  for(std::vector<int>::const_iterator it = orbit.begin();
      it != orbit.end(); ++it) {
    if (*it >= n_) return true;
  }
  return false;
  
}

bool
CbcNauty::isMixedFreeZeroOrbit(const std::vector<int> &orbit) const
{
  bool containsFree = false;
  bool containsZero = false;

  for(std::vector<int>::const_iterator it = orbit.begin();
      it != orbit.end(); ++it) {
    if (*it >= n_) return false;    
    if (vstat_[*it] == FREE) containsFree = true;
    if (vstat_[*it] == FIX_AT_ZERO) containsZero = true;    
    if (containsFree && containsZero) break;
  }  
  return (containsFree && containsZero);
}
*/

void CbcNauty::setWriteAutoms(const std::string &fname)
{
  afp_ = fopen(fname.c_str(), "w");
  options_->writeautoms = TRUE;
#ifndef NTY_TRACES
  options_->writemarkers = FALSE;
#endif
  options_->outfile = afp_;
}

void CbcNauty::unsetWriteAutoms()
{
  fclose(afp_);
  options_->writeautoms = FALSE;
}

// Default Constructor
CbcOrbitalBranchingObject::CbcOrbitalBranchingObject()
  : CbcBranchingObject()
  , column_(-1)
  , numberOther_(0)
  , numberExtra_(0)
  , fixToZero_(NULL)
{
}

// Useful constructor
CbcOrbitalBranchingObject::CbcOrbitalBranchingObject(CbcModel *model, int column,
  int way,
  int numberExtra,
  const int *extraToZero)
  : CbcBranchingObject(model, -1, way, 0.5)
  , column_(column)
  , numberOther_(0)
  , numberExtra_(0)
  , fixToZero_(NULL)
{
  CbcSymmetry *symmetryInfo = model->symmetryInfo();
  assert(symmetryInfo);
  // Filled in (hopefully)
  const int *orbit = symmetryInfo->whichOrbit();
  int iOrbit = orbit[column];
  assert(iOrbit >= 0);
  int numberColumns = model->getNumCols();
  numberOther_ = -1;
  for (int i = 0; i < numberColumns; i++) {
    if (orbit[i] == iOrbit)
      numberOther_++;
  }
  assert(numberOther_ > 0);
  nautyBranchSucceeded_++;
  nautyOtherBranches_ += numberOther_;
  numberExtra_ = numberExtra;
  fixToZero_ = new int[numberOther_ + numberExtra_];
  int n = 0;
  for (int i = 0; i < numberColumns; i++) {
    if (orbit[i] == iOrbit && i != column)
      fixToZero_[n++] = i;
  }
  for (int i = 0; i < numberExtra; i++) {
    fixToZero_[n++] = extraToZero[i];
  }
}

// Copy constructor
CbcOrbitalBranchingObject::CbcOrbitalBranchingObject(const CbcOrbitalBranchingObject &rhs)
  : CbcBranchingObject(rhs)
  , column_(rhs.column_)
  , numberOther_(rhs.numberOther_)
  , numberExtra_(rhs.numberExtra_)
{
  fixToZero_ = CoinCopyOfArray(rhs.fixToZero_, numberOther_ + numberExtra_);
}

// Assignment operator
CbcOrbitalBranchingObject &
CbcOrbitalBranchingObject::operator=(const CbcOrbitalBranchingObject &rhs)
{
  if (this != &rhs) {
    CbcBranchingObject::operator=(rhs);
    delete[] fixToZero_;
    column_ = rhs.column_;
    numberOther_ = rhs.numberOther_;
    numberExtra_ = rhs.numberExtra_;
    fixToZero_ = CoinCopyOfArray(rhs.fixToZero_, numberOther_ + numberExtra_);
  }
  return *this;
}
CbcBranchingObject *
CbcOrbitalBranchingObject::clone() const
{
  return (new CbcOrbitalBranchingObject(*this));
}

// Destructor
CbcOrbitalBranchingObject::~CbcOrbitalBranchingObject()
{
  delete[] fixToZero_;
}

double
CbcOrbitalBranchingObject::branch()
{
  decrementNumberBranchesLeft();
  if (model_->logLevel() > 1)
    print();
  OsiSolverInterface *solver = model_->solver();
  if (way_ < 0) {
    solver->setColUpper(column_, 0.0);
    for (int i = 0; i < numberOther_ + numberExtra_; i++) {
      solver->setColUpper(fixToZero_[i], 0.0);
    }
    way_ = 1; // Swap direction
  } else {
    solver->setColLower(column_, 1.0);
    for (int i = numberOther_; i < numberOther_ + numberExtra_; i++) {
      solver->setColUpper(fixToZero_[i], 0.0);
    }
    way_ = -1; // Swap direction
  }
  return 0.0;
}
/* Update bounds in solver as in 'branch' and update given bounds.
   branchState is -1 for 'down' +1 for 'up' */
void CbcOrbitalBranchingObject::fix(OsiSolverInterface *solver,
  double *lower, double *upper,
  int branchState) const
{
  if (branchState < 0) {
    upper[column_] = 0.0;
    for (int i = 0; i < numberOther_ + numberExtra_; i++) {
      upper[fixToZero_[i]] = 0.0;
      ;
    }
  } else {
    lower[column_] = 1.0;
    for (int i = numberOther_; i < numberOther_ + numberExtra_; i++) {
      upper[fixToZero_[i]] = 0.0;
      ;
    }
  }
}
// Print what would happen
void CbcOrbitalBranchingObject::print()
{
  if (way_ < 0) {
    printf("Orbital Down - to zero %d", column_);
    for (int i = 0; i < numberOther_ + numberExtra_; i++) {
      printf(" %d", fixToZero_[i]);
    }
  } else {
    printf("Orbital Up - to one %d, to zero", column_);
    for (int i = numberOther_; i < numberOther_ + numberExtra_; i++) {
      printf(" %d", fixToZero_[i]);
    }
  }
  printf("\n");
}

/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type.
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int CbcOrbitalBranchingObject::compareOriginalObject(const CbcBranchingObject *brObj) const
{
  const CbcOrbitalBranchingObject *br = dynamic_cast< const CbcOrbitalBranchingObject * >(brObj);
  assert(!br);
  abort();
  return 0;
}

/** Compare the \c this with \c brObj. \c this and \c brObj must be os the
    same type and must have the same original object, but they may have
    different feasible regions.
    Return the appropriate CbcRangeCompare value (first argument being the
    sub/superset if that's the case). In case of overlap (and if \c
    replaceIfOverlap is true) replace the current branching object with one
    whose feasible region is the overlap.
*/
CbcRangeCompare
CbcOrbitalBranchingObject::compareBranchingObject(const CbcBranchingObject *brObj, const bool replaceIfOverlap)
{
  const CbcOrbitalBranchingObject *br = dynamic_cast< const CbcOrbitalBranchingObject * >(brObj);
  assert(!br);
  abort();
  return CbcRangeDisjoint;
}
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
