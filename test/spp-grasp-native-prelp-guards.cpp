// Copyright (C) 2026, COIN-OR Foundation and others.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CbcSolver.hpp"
#include "CbcEventHandler.hpp"
#include "CbcSOS.hpp"
#include "CglProbing.hpp"
#include "OsiClpSolverInterface.hpp"

#include <cmath>
#include <deque>
#include <iostream>
#include <string>

namespace {
int callbackCalls = 0;
bool callbackSawNative = false;

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
  OsiClpSolverInterface osi;
  if (osi.readMps(argv[1], ""))
    return 2;
  CbcModel model(osi);
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
