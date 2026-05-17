// Copyright (C) 2010, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).
#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "CoinPragma.hpp"
#include "CglZeroHalf.hpp" 
//#include "CglKnapsackCover.hpp" 
#include <stdio.h>

namespace {

std::string canonicalCut(const OsiRowCut &cut)
{
  const CoinPackedVector &row = cut.row();
  std::vector< std::pair<int, int> > entries;
  entries.reserve(row.getNumElements());
  for (int i = 0; i < row.getNumElements(); ++i)
    entries.push_back(std::make_pair(row.getIndices()[i], static_cast<int>(std::lround(row.getElements()[i]))));
  std::sort(entries.begin(), entries.end());

  std::ostringstream out;
  if (cut.lb() <= -1.0e20)
    out << "lb=-inf;";
  else
    out << "lb=" << static_cast<int>(std::lround(cut.lb())) << ';';
  if (cut.ub() >= 1.0e20)
    out << "ub=inf;";
  else
    out << "ub=" << static_cast<int>(std::lround(cut.ub())) << ';';
  for (std::vector< std::pair<int, int> >::const_iterator it = entries.begin(); it != entries.end(); ++it)
    out << it->first << ':' << it->second << ',';
  return out.str();
}

std::vector<std::string> canonicalCuts(const OsiCuts &cuts)
{
  std::vector<std::string> result;
  result.reserve(cuts.sizeRowCuts());
  for (int i = 0; i < cuts.sizeRowCuts(); ++i)
    result.push_back(canonicalCut(cuts.rowCut(i)));
  std::sort(result.begin(), result.end());
  return result;
}

}

//--------------------------------------------------------------------------
// test the zero half cut generators methods.
void
CglZeroHalfUnitTest(
  const OsiSolverInterface * baseSiP,
  const std::string mpsDir )
{

  // Test default constructor
  {
    CglZeroHalf cg;
  }

  // Test copy & assignment
  {
    CglZeroHalf rhs;
    {
      CglZeroHalf cg;
      CglZeroHalf cgC(cg);
      rhs=cg;
    }
  }



  // Test generate cuts method on lseu
  {
    CglZeroHalf cg;
    
    OsiSolverInterface * siP = baseSiP->clone();
    std::string fn = mpsDir+"lseu.mps";
    siP->readMps(fn.c_str(),"");
    // test if there
    if (!siP->getNumRows()) {
      printf("** Unable to find lseu in %s\n",
	     mpsDir.c_str());
      return;
    }
    siP->initialSolve();
    cg.refreshSolver(siP);
    OsiCuts cuts;
    cg.generateCuts(*siP,cuts);
    CglZeroHalf sparseCg;
    sparseCg.setSepGraphSparseThreshold(0);
    sparseCg.refreshSolver(siP);
    OsiCuts sparseCuts;
    sparseCg.generateCuts(*siP,sparseCuts);
    assert(canonicalCuts(cuts) == canonicalCuts(sparseCuts));

    // lseu is the optimal solution to lseu
    // Optimal IP solution to lseu    
    int objIndices[13]={0,1,6,13,26,33,38,43,50,52,63,65,85};
    CoinPackedVector lseu(13,objIndices,1.0);

    // test that none of the generated cuts
    // chops off the optimal solution
    int nRowCuts = cuts.sizeRowCuts();
    OsiRowCut rcut;
    CoinPackedVector rpv;
    int i;
    for (i=0; i<nRowCuts; i++){
      rcut = cuts.rowCut(i);
      rpv = rcut.row();
      double lseuSum = (rpv*lseu).sum();
      double rcutub = rcut.ub();
      assert (lseuSum <= rcutub);
    }

    // test that the cuts improve the 
    // lp objective function value
    double lpRelaxBefore=siP->getObjValue();
    OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cuts);
    siP->resolve();
    double lpRelaxAfter=siP->getObjValue(); 
#ifdef CGL_DEBUG
    printf("\n\nOrig LP min=%f\n",lpRelaxBefore);
    printf("Final LP min=%f\n\n",lpRelaxAfter);
#endif
    printf("Zero cuts %d\n",nRowCuts);
    if (!(lpRelaxBefore < lpRelaxAfter)){
       printf("***Warning: Bound did not improve after addition of cut.\n");
       printf("***This can happen, but is generally not expected\n");
    }else{
       printf("Good zero %s\n",fn.c_str());
    }

    delete siP;

  }

}
