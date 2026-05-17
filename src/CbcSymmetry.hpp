/*
 * Stub symmetry interface kept for source compatibility after nauty removal.
 * This file is licensed under the Eclipse Public License (EPL).
 */
#ifndef CBC_SYMMETRY_HPP
#define CBC_SYMMETRY_HPP

#include <vector>

class CbcModel;
class OsiSolverInterface;

typedef struct {
  int numberInPerm;
  int numberPerms;
  int *orbits;
} cbc_permute;

class CbcSymmetry {
public:
  CbcSymmetry() = default;
  CbcSymmetry(const CbcSymmetry &) = default;
  CbcSymmetry &operator=(const CbcSymmetry &) = default;
  ~CbcSymmetry() = default;

  std::vector<int> *Find_Orbit(int) const { return new std::vector<int>; }
  void Compute_Symmetry() const {}
  int statsOrbits(CbcModel *, int) const { return 0; }
  void Print_Orbits(int = 0) const {}
  void fillOrbits() {}
  int orbitalFixing(OsiSolverInterface *) { return 0; }
  int orbitalFixing2(OsiSolverInterface *) { return 0; }
  int *whichOrbit() { return nullptr; }
  int *fixedToZero() const { return nullptr; }
  int numberUsefulOrbits() const { return 0; }
  int numberUsefulObjects() const { return 0; }
  int largestOrbit(const double *, const double *) const { return -1; }
  void ChangeBounds(const double *, const double *, int, bool) const {}
  int changeBounds(int, double *, double *, OsiSolverInterface *, int) const { return 0; }
  int changeBounds(double *, double *, OsiSolverInterface *) const { return 0; }
  int changeBounds2(double *, double *, OsiSolverInterface *) const { return 0; }
  int fixSome(int, double *, double *) const { return 0; }
  int worthBranching(const double *, const double *, int, int &numberCouldFix) const
  {
    numberCouldFix = 0;
    return 0;
  }
  void fixSuccess(int) {}
  void adjustStats(const CbcSymmetry *) {}
  int numberColumns() const { return 0; }
  void *getNtyInfo() { return nullptr; }
  void setupSymmetry(CbcModel *) {}
  void addPermutation(cbc_permute permutation) { delete[] permutation.orbits; }
  int numberPermutations() const { return 0; }
  int *permutation(int) const { return nullptr; }
  int numberInPermutation(int) const { return 0; }
  void incrementNautyBranches(int) {}
  void incrementBranchSucceeded() {}
  double getNautyTime() const { return 0.0; }
};

#endif
