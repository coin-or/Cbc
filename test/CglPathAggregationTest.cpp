#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "CoinPackedVector.hpp"
#include "CglPathAggregation.hpp"
#include "OsiClpSolverInterface.hpp"

namespace {
double rowActivity(const OsiRowCut &cut, const double *solution)
{
  CoinPackedVector row = cut.row();
  double activity = 0.0;
  const int *indices = row.getIndices();
  const double *elements = row.getElements();
  for (int i = 0; i < row.getNumElements(); ++i)
    activity += elements[i] * solution[indices[i]];
  return activity;
}

void assertCutsKeepIntegerFeasiblePoints(const OsiCuts &cuts)
{
  for (int y = 0; y <= 1; ++y) {
    for (int z = 0; z <= 1; ++z) {
      const double x = 3.0 * y;
      if (x + 3.0 * z > 4.5 + 1.0e-9)
        continue;
      const double solution[] = { x, static_cast< double >(y), static_cast< double >(z) };
      for (int i = 0; i < cuts.sizeRowCuts(); ++i)
        assert(rowActivity(cuts.rowCut(i), solution) <= cuts.rowCut(i).ub() + 1.0e-8);
    }
  }
}

OsiClpSolverInterface buildPathAggregationFixture()
{
  OsiClpSolverInterface si;
  si.messageHandler()->setLogLevel(0);

  CoinPackedVector empty;
  si.addCol(empty, 0.0, 3.0, 0.0); // x
  si.addCol(empty, 0.0, 1.0, -1.0); // y
  si.addCol(empty, 0.0, 1.0, -1.1); // z
  si.setInteger(1);
  si.setInteger(2);

  CoinPackedVector link;
  link.insert(0, 1.0);
  link.insert(1, -3.0);
  si.addRow(link, 0.0, 0.0);

  CoinPackedVector capacity;
  capacity.insert(0, 1.0);
  capacity.insert(2, 3.0);
  si.addRow(capacity, -si.getInfinity(), 4.5);

  return si;
}
}

void CglPathAggregationUnitTest()
{
  {
    CglPathAggregation generator;
    generator.setMinViolation(1.0e-8);
    OsiClpSolverInterface si = buildPathAggregationFixture();
    si.initialSolve();

    const double *solution = si.getColSolution();
    assert(fabs(solution[0] - 1.5) < 1.0e-7);
    assert(fabs(solution[1] - 0.5) < 1.0e-7);
    assert(fabs(solution[2] - 1.0) < 1.0e-7);

    OsiCuts cuts;
    generator.generateCuts(si, cuts);
    std::cout << "Generated " << cuts.sizeRowCuts() << " cuts:" << std::endl;
    for (int i = 0; i < cuts.sizeRowCuts(); ++i) {
      const OsiRowCut &rc = cuts.rowCut(i);
      std::cout << "  Cut " << i << ": ub=" << rc.ub() << " elements: ";
      for (int j = 0; j < rc.row().getNumElements(); ++j) {
        std::cout << rc.row().getIndices()[j] << ":" << rc.row().getElements()[j] << " ";
      }
      std::cout << std::endl;
    }
    assert(cuts.sizeRowCuts() > 0);
    assertCutsKeepIntegerFeasiblePoints(cuts);

    const double lpValue = si.getObjValue();
    si.applyCuts(cuts);
    si.resolve();
    assert(si.getObjValue() > lpValue + 1.0e-7);
  }

  {
    CglPathAggregation generator;
    generator.setMinViolation(1.0e-8);
    OsiClpSolverInterface si;
    si.messageHandler()->setLogLevel(0);
    CoinPackedVector empty;
    si.addCol(empty, 0.0, 1.0, 0.0);
    si.addCol(empty, 0.0, 2.0, 0.0);
    si.setInteger(0);
    si.setInteger(1);

    std::vector< int > active;
    active.push_back(0);
    active.push_back(1);
    std::vector< double > dense(2, 0.0);
    dense[0] = -3.0;
    dense[1] = 3.0;
    const double lpSolution[] = { 0.8, 1.2 };

    OsiRowCut cut;
    bool generated = generator.generateMirCutFromRow(si, active, dense, 1.5, lpSolution, cut);
    assert(generated);
    assert(fabs(cut.ub() - 0.0) < 1.0e-9);
    assert(rowActivity(cut, lpSolution) > cut.ub() + 1.0e-8);

    const double feasibleA[] = { 0.0, 0.0 };
    const double feasibleB[] = { 1.0, 1.0 };
    assert(rowActivity(cut, feasibleA) <= cut.ub() + 1.0e-9);
    assert(rowActivity(cut, feasibleB) <= cut.ub() + 1.0e-9);
  }
}

int main()
{
  CglPathAggregationUnitTest();
  std::cout << "CglPathAggregationTest passed" << std::endl;
  return 0;
}
