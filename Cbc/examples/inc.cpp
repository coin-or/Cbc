// inc.cpp: example of event handler to save
// every incumbent solution to a new file

#include <cstdio>
#include <cstdlib>
#include <OsiCbcSolverInterface.hpp>
#include <CbcModel.hpp>
#include <CglPreProcess.hpp>
#include <CbcSolver.hpp>
#include <Cbc_C_Interface.h>
#include "CbcEventHandler.hpp"

static int callBack(CbcModel *model, int whereFrom)
{
  return 0;
}

class SolHandler : public CbcEventHandler {
public:
  virtual CbcAction event(CbcEvent whichEvent);
  SolHandler();
  SolHandler(CbcModel *model);
  virtual ~SolHandler();
  SolHandler(const SolHandler &rhs);
  SolHandler &operator=(const SolHandler &rhs);
  virtual CbcEventHandler *clone() const;

  double bestCost;
};

SolHandler::SolHandler()
  : CbcEventHandler()
  , bestCost(COIN_DBL_MAX)
{
}

SolHandler::SolHandler(const SolHandler &rhs)
  : CbcEventHandler(rhs)
  , bestCost(rhs.bestCost)
{
}

SolHandler::SolHandler(CbcModel *model)
  : CbcEventHandler(model)
  , bestCost(COIN_DBL_MAX)
{
}

SolHandler::~SolHandler()
{
}

SolHandler &SolHandler::operator=(const SolHandler &rhs)
{
  if (this != &rhs) {
    CbcEventHandler::operator=(rhs);
    this->bestCost = rhs.bestCost;
  }
  return *this;
}

CbcEventHandler *SolHandler::clone() const
{
  return new SolHandler(*this);
}

CbcEventHandler::CbcAction SolHandler::event(CbcEvent whichEvent)
{
  // If in sub tree carry on
  if ((model_->specialOptions() & 2048) == 0) {
    if ((whichEvent == solution || whichEvent == heuristicSolution)) {
      OsiSolverInterface *origSolver = model_->solver();
      const OsiSolverInterface *pps = model_->postProcessedSolver(1);

      const OsiSolverInterface *solver = pps ? pps : origSolver;

      if (bestCost > solver->getObjValue() + 1e-6) {
        bestCost = solver->getObjValue();
        char outFName[256] = "", strCost[256] = "";
        sprintf(strCost, "%.12e", bestCost);
        char *s = strstr(strCost, ".");
        if (s)
          *s = 'p';
        sprintf(outFName, "inc%s.sol", strCost);
        printf("** solution improved to %g, saving it to file %s ** \n", solver->getObjValue(), outFName);

        // saving .sol
        FILE *f = fopen(outFName, "w");
        assert(f);
        fprintf(f, "Stopped on iterations - objective value %.8f\n", solver->getObjValue());
        const double *x = solver->getColSolution();
        const double *obj = solver->getObjCoefficients();
        for (int i = 0; (i < solver->getNumCols()); ++i) {
          if (fabs(x[i]) < 1e-6)
            continue;
          fprintf(f, "%-7d %22s %22g %g\n", i,
            solver->getColName(i).c_str(), x[i], obj[i]);
        }
        fclose(f);
      } // improved cost
    } // solution found
  } // not in subtree

  return noAction;
}

int main(int argc, char **argv)
{
  if (argc < 2) {
    fprintf(stderr, "enter instance name\n");
    exit(EXIT_FAILURE);
  }

  CbcSolverUsefulData cbcData;
  cbcData.noPrinting_ = false;
  OsiClpSolverInterface lp;

  lp.readMps(argv[1]);

  CbcModel model(lp);

  SolHandler sh;
  model.passInEventHandler(&sh);

  CbcMain0(model, cbcData);
  CbcMain1(argc - 1, (const char **)(argv + 1), model, callBack, cbcData);

  return 0;
}
