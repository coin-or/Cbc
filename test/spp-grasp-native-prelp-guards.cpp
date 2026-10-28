// Copyright (C) 2026, COIN-OR Foundation and others.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CbcSolver.hpp"
#include "CbcEventHandler.hpp"
#include "CbcSOS.hpp"
#include "CglProbing.hpp"
#include "OsiClpSolverInterface.hpp"

#include <cmath>
#include <chrono>
#include <deque>
#include <iostream>
#include <string>
#include <thread>

namespace {
int callbackCalls = 0;
bool callbackSawNative = false;
bool preprocessingArmed = false;
bool preprocessingTimedOut = false;

class PreprocessingTimeoutSolver : public OsiClpSolverInterface {
public:
  OsiSolverInterface *clone(bool copyData = true) const override
  {
    return copyData ? new PreprocessingTimeoutSolver(*this)
                    : new PreprocessingTimeoutSolver;
  }

  void initialSolve() override
  {
    bool inBranchAndCut = false;
    getHintParam(OsiDoInBranchAndCut, inBranchAndCut);
    if (preprocessingArmed && inBranchAndCut) {
      // Inject an LP time limit only after CBC has entered preprocessing.
      std::this_thread::sleep_for(std::chrono::milliseconds(1100));
      getModelPtr()->setProblemStatus(3);
      getModelPtr()->setSecondaryStatus(9);
      preprocessingTimedOut = true;
      return;
    }
    OsiClpSolverInterface::initialSolve();
  }
};

int preprocessingCallback(CbcModel *, int where)
{
  if (where == 1)
    preprocessingArmed = true;
  return 0;
}

int callback(CbcModel *model, int where)
{
  ++callbackCalls;
  if (where == 1 && model->bestSolution())
    callbackSawNative = true;
  return 0;
}
}

int main(int argc, char **argv)
{
  if (argc != 4)
    return 2;
  const std::string mode(argv[2]);
  PreprocessingTimeoutSolver osi;
  if (osi.readMps(argv[1], ""))
    return 2;
  CbcModel model(osi);
  if (mode.find("preprocess-") == 0) {
    if (mode == "preprocess-existing") {
      const double incumbent[] = { 1, 1, 1, 0, 0, 0 };
      model.setBestSolution(incumbent, 6, 6.0, false);
    }
    CbcSolver solver(model);
    solver.initialize();
    const std::string solution = std::string(argv[3]) + ".sol";
    std::deque< std::string > arguments = {
      "-log", "2", "-timeMode", "elapsed", "-threads", "0",
      "-presolve", "off", "-preprocess", "on",
      "-preRootLPStrenghtening", "off", "-sppGrasp",
      mode == "preprocess-original" ? "on" : "off",
      "-sppGraspWork", ".001", "-sppGraspNoImprove", "1",
      "-sec", "1", "-progressCsv", argv[3],
      "-solve", "-solu", solution, "-quit"
    };
    const int result = solver.run(arguments, preprocessingCallback);
    if (result || !preprocessingTimedOut) {
      std::cerr << "did not reach the preprocessing LP time limit\n";
      return 1;
    }
    return 0;
  }
  const double better[] = { 0.0, 0.0, 1.0 };
  const double worse[] = { 1.0, 1.0, 0.0 };
  if (mode == "better" || mode == "worse") {
    model.setBestSolution(mode == "better" ? better : worse, 3,
      mode == "better" ? 1.0 : 5.0, false);
    model.setCutoff(100.0);
  } else if (mode == "event") {
    CbcEventHandler handler(&model);
    model.passInEventHandler(&handler);
  } else if (mode == "cuts") {
    CglProbing probing;
    model.addCutGenerator(&probing, 1, "custom");
  } else if (mode == "sos") {
    int members[] = { 0, 1 };
    double weights[] = { 1.0, 2.0 };
    CbcSOS sos(&model, 2, members, weights, 0, 1);
    CbcObject *object = &sos;
    model.addObjects(1, &object);
  } else if (mode == "hotstart") {
    model.setHotstartSolution(better);
  }
  CbcSolver solver(model);
  solver.initialize();
  std::deque< std::string > arguments = {
    "-log", "2", "-preprocess", "off", "-presolve", "off",
    "-preRootLPStrenghtening", "off", "-sppGraspNativePreLp", "on",
    "-sppGraspWork", ".001", "-sppGraspNoImprove", "1",
    "-sec", "5", "-progressCsv", argv[3], "-solve", "-quit"
  };
  int result = solver.run(arguments, mode == "callback" ? callback : nullptr);
  if (result || callbackSawNative || (mode == "callback" && !callbackCalls)) {
    std::cerr << "native completion bypassed callback guard\n";
    return 1;
  }
  if (mode == "sos")
    return 0;
  if (!solver.hasSolution() || std::fabs(solver.objectiveValue() - 1.0) > 1.0e-8) {
    std::cerr << "incumbent or final solution was not preserved\n";
    return 1;
  }
  return 0;
}
