/*
 * cbcmain1-test — checks what CbcMain1() leaves in the caller's CbcModel.
 *
 * CbcMain1() solves inside a CbcSolver that has its own copy of the
 * caller's model, and then copies the results back. Before CbcSolver,
 * CbcMain1() worked on the caller's model itself and ended with moveInfo()
 * from the branch-and-bound model, so programs that call CbcMain1() (OS,
 * for example) query the caller's model afterwards: isNodeLimitReached()
 * and isSolutionLimitReached(), which compare the node and solution counts
 * with the limits set by -maxNodes and -maxSolutions, getSolutionCount(),
 * getBestPossibleObjValue() and the other parameters set on the command
 * line. The copy-back used to leave these out, so after -maxNodes,
 * isNodeLimitReached() was false (the caller's model still had no node
 * limit), the solution count was 0, and the best possible objective was
 * the incumbent's objective.
 *
 * The fixture is a multidimensional knapsack (40 binaries, 3 rows, fixed
 * pseudo-random data) that takes many nodes with preprocessing, cuts and
 * heuristics off, so both limits stop the search before it is complete.
 *
 * Usage: cbcmain1-test
 * Exit code: 0 = all checks passed, 1 = one or more checks failed.
 */

#include "CbcCompareDefault.hpp"
#include "CbcEventHandler.hpp"
#include "CbcFullNodeInfo.hpp"
#include "CbcModel.hpp"
#include "CbcSolver.hpp"
#include "CoinPackedMatrix.hpp"
#include "OsiClpSolverInterface.hpp"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

int g_failures = 0;

void check(bool cond, const std::string &msg)
{
  if (!cond) {
    ++g_failures;
    printf("  FAIL: %s\n", msg.c_str());
  } else {
    printf("  ok:   %s\n", msg.c_str());
  }
}

/* max sum p_j x_j s.t. sum w_ij x_j <= c_i (i = 0, 1, 2), x binary */
void loadKnapsack(OsiClpSolverInterface &solver)
{
  const int numberColumns = 40;
  const int numberRows = 3;
  unsigned int seed = 12345;
  std::vector< double > profit(numberColumns);
  std::vector< double > capacity(numberRows, 0.0);
  CoinPackedMatrix matrix(false, 0, 0);
  matrix.setDimensions(0, numberColumns);
  for (int j = 0; j < numberColumns; j++) {
    seed = 1103515245u * seed + 12345u;
    profit[j] = 10 + (seed >> 16) % 51;
  }
  for (int i = 0; i < numberRows; i++) {
    std::vector< int > index(numberColumns);
    std::vector< double > weight(numberColumns);
    for (int j = 0; j < numberColumns; j++) {
      seed = 1103515245u * seed + 12345u;
      index[j] = j;
      weight[j] = 10 + (seed >> 16) % 51;
      capacity[i] += weight[j];
    }
    capacity[i] = floor(capacity[i] / 2);
    matrix.appendRow(numberColumns, &index[0], &weight[0]);
  }
  std::vector< double > columnLower(numberColumns, 0.0);
  std::vector< double > columnUpper(numberColumns, 1.0);
  std::vector< double > rowLower(numberRows, -COIN_DBL_MAX);
  solver.loadProblem(matrix, &columnLower[0], &columnUpper[0], &profit[0],
    &rowLower[0], &capacity[0]);
  for (int j = 0; j < numberColumns; j++)
    solver.setInteger(j);
  solver.setObjSense(-1.0);
  solver.messageHandler()->setLogLevel(0);
}

/* Solves the knapsack with CbcMain1() and the given options, which come
 * before -solve, and returns the caller's model through model. */
void solve(CbcModel &model, const std::vector< std::string > &options)
{
  CbcParameters parameters;
  CbcMain0(model, parameters);
  std::vector< const char * > argv;
  argv.push_back("cbcmain1-test");
  argv.push_back("-log=0");
  for (size_t i = 0; i < options.size(); i++)
    argv.push_back(options[i].c_str());
  argv.push_back("-solve");
  argv.push_back("-quit");
  CbcMain1(static_cast< int >(argv.size()), &argv[0], model, parameters);
}

const char *const searchOff[] = { "-preprocess=off", "-cutsOnOff=off",
  "-heuristicsOnOff=off" };

std::vector< std::string > optionsWith(const char *option)
{
  std::vector< std::string > options(searchOff, searchOff + 3);
  options.push_back(option);
  return options;
}

void testNodeLimit()
{
  printf("node limit (-maxNodes=1):\n");
  OsiClpSolverInterface solver;
  loadKnapsack(solver);
  CbcModel model(solver);
  solve(model, optionsWith("-maxNodes=1"));
  check(model.status() == 1 && model.secondaryStatus() == 3,
    "status is stopped on nodes");
  check(model.getMaximumNodes() == 1, "getMaximumNodes() is 1");
  check(model.isNodeLimitReached(), "isNodeLimitReached()");
  check(!model.isSolutionLimitReached(), "!isSolutionLimitReached()");
}

void testSolutionLimit()
{
  printf("solution limit (-maxSolutions=1):\n");
  OsiClpSolverInterface solver;
  loadKnapsack(solver);
  CbcModel model(solver);
  solve(model, optionsWith("-maxSolutions=1"));
  check(model.status() == 1 && model.secondaryStatus() == 6,
    "status is stopped on solutions");
  check(model.getMaximumSolutions() == 1, "getMaximumSolutions() is 1");
  check(model.getSolutionCount() >= 1, "getSolutionCount() >= 1");
  check(model.isSolutionLimitReached(), "isSolutionLimitReached()");
  check(!model.isNodeLimitReached(), "!isNodeLimitReached()");
  check(model.bestSolution() != NULL, "bestSolution() is set");
  // The search is incomplete, so the best possible objective is a bound
  // better than the incumbent's (larger, as the problem is a maximization)
  printf("  objective %g, best possible %g\n", model.getObjValue(),
    model.getBestPossibleObjValue());
  check(model.getBestPossibleObjValue() > model.getObjValue() + 1.0e-6,
    "getBestPossibleObjValue() is a bound better than the incumbent");
}

void testParameters()
{
  printf("parameters (-integerTolerance=1e-5 -allowableGap=0.5):\n");
  OsiClpSolverInterface solver;
  loadKnapsack(solver);
  CbcModel model(solver);
  std::vector< std::string > options;
  options.push_back("-integerTolerance=1e-5");
  options.push_back("-allowableGap=0.5");
  solve(model, options);
  check(model.getIntegerTolerance() == 1.0e-5, "getIntegerTolerance() is 1e-5");
  check(model.getAllowableGap() == 0.5, "getAllowableGap() is 0.5");
}

void testOptimal()
{
  printf("complete search:\n");
  OsiClpSolverInterface solver;
  loadKnapsack(solver);
  CbcModel model(solver);
  solve(model, std::vector< std::string >());
  check(model.isProvenOptimal(), "isProvenOptimal()");
  check(!model.isNodeLimitReached(), "!isNodeLimitReached()");
  check(!model.isSolutionLimitReached(), "!isSolutionLimitReached()");
  check(model.getSolutionCount() >= 1, "getSolutionCount() >= 1");
  check(fabs(model.getBestPossibleObjValue() - model.getObjValue()) <= 1.0e-6,
    "getBestPossibleObjValue() equals the optimal objective");
}

void testLegacyNodeOrdering()
{
  CbcCompareDefault depth(-3.0);
  depth.sayThreaded();
  CbcNode nodes[3];
  for (int i = 0; i < 3; ++i) {
    nodes[i].setNodeNumber(i);
    nodes[i].setDepth(i + 1);
    nodes[i].setNumberUnsatisfied(6000 * i);
  }
  check(depth.test(&nodes[0], &nodes[1]) && depth.test(&nodes[1], &nodes[2])
      && depth.test(&nodes[0], &nodes[2]) && !depth.test(&nodes[2], &nodes[0]),
    "depth-first order is transitive across large infeasibility differences");

  OsiClpSolverInterface solver;
  loadKnapsack(solver);
  CbcModel model(solver);
  CbcCompareDefault fewest(-2.0);
  check(fewest.every1000Nodes(&model, 1002) && fewest.getWeight() == -1.0,
    "fewest-first transitions even if the node count skips 1000");
}

void testLegacyDiveSchedule()
{
  class DiveProbe : public CbcCompareDefault {
  public:
    bool diving() const { return startNodeNumber_ >= 0; }
  } compare;
  OsiClpSolverInterface solver;
  loadKnapsack(solver);
  CbcModel model(solver);
  model.tree()->setComparison(compare);
  CbcNode nodes[2];
  for (int i = 0; i < 2; ++i) {
    CbcFullNodeInfo *info = new CbcFullNodeInfo();
    info->setNumberBranchesLeft(1);
    info->setNumberPointingToThis(1);
    info->setNodeNumber(i);
    nodes[i].setNodeInfo(info);
    nodes[i].setObjectiveValue(i);
    model.tree()->push(&nodes[i]);
  }
  model.setNodeCount(1002);
  compare.updateDive(&model);
  check(compare.diving() && !compare.test(model.tree()->top(), model.tree()->top()),
    "skipped entry starts a dive without violating self-comparison");
  model.setNodeCount(2001);
  compare.updateDive(&model);
  check(!compare.diving(), "leaving the window ends the dive");
  model.setNodeCount(5002);
  compare.updateDive(&model);
  check(compare.diving(), "a subsequent window starts another dive");
  model.tree()->setComparison(compare);
  model.setNodeCount(5003);
  compare.updateDive(&model);
  check(!compare.diving(), "an early heap rebuild does not restart the same dive");
  while (!model.tree()->empty())
    model.tree()->pop();
}

void testDirectionalFewest()
{
  class StrategyProbe : public CbcEventHandler {
  public:
    StrategyProbe(int way, bool &observed)
      : way_(way)
      , observed_(observed)
    {
    }
    CbcEventHandler *clone() const { return new StrategyProbe(*this); }
    CbcAction event(CbcEvent event)
    {
      if (event == startUp && !model_->parentModel()) {
        const CbcCompareDefault *compare = dynamic_cast< const CbcCompareDefault * >(model_->nodeComparison());
        observed_ = compare && compare->getWeight() == -2.0 && model_->getPreferredWay() == way_;
      }
      return noAction;
    }

  private:
    int way_;
    bool &observed_;
  };
  for (int way : { 1, -1 }) {
    OsiClpSolverInterface solver;
    loadKnapsack(solver);
    CbcModel model(solver);
    bool observed = false;
    StrategyProbe probe(way, observed);
    model.passInEventHandler(&probe);
    auto options = optionsWith(way > 0 ? "-nodeStrategy=upfewest" : "-nodeStrategy=downfewest");
    options.push_back("-maxNodes=1");
    solve(model, options);
    check(observed, way > 0 ? "upfewest sets both traversal and up direction" : "downfewest sets both traversal and down direction");
  }
}

} // namespace

int main(int argc, char **argv)
{
  if (argc > 1) {
    printf("Usage: %s\nChecks what CbcMain1() leaves in the caller's CbcModel.\n",
      argv[0]);
    return 0;
  }
  testNodeLimit();
  testSolutionLimit();
  testParameters();
  testOptimal();
  testLegacyNodeOrdering();
  testLegacyDiveSchedule();
  testDirectionalFewest();
  if (g_failures) {
    printf("%d check(s) failed\n", g_failures);
    return 1;
  }
  printf("all checks passed\n");
  return 0;
}
