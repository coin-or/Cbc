/* $Id$ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif

#include "CbcConfig.h"

#include <string>
//#define CBC_DEBUG 1
//#define CHECK_CUT_COUNTS
//#define CHECK_NODE
//#define CHECK_NODE_FULL
//#define NODE_LOG
//#define GLOBAL_CUTS_JUST_POINTERS
#ifdef CGL_DEBUG_GOMORY
extern int gomory_try;
#endif
#include <cassert>
#include <cmath>
#include <cfloat>
#ifdef COIN_HAS_CLP
// include Presolve from Clp
#include "ClpPresolve.hpp"
#include "OsiClpSolverInterface.hpp"
#include "ClpNode.hpp"
#include "ClpDualRowDantzig.hpp"
#include "ClpSimplexPrimal.hpp"
#endif

#include "CbcEventHandler.hpp"

#include "OsiSolverInterface.hpp"
#include "OsiAuxInfo.hpp"
#include "OsiSolverBranch.hpp"
#include "OsiChooseVariable.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinHelperFunctions.hpp"
#include "CbcBranchActual.hpp"
#include "CbcBranchDynamic.hpp"
#include "CbcHeuristic.hpp"
#include "CbcHeuristicFPump.hpp"
#include "CbcHeuristicRINS.hpp"
#include "CbcHeuristicDive.hpp"
#include "CbcModel.hpp"
#include "CbcTreeLocal.hpp"
#include "CbcStatistics.hpp"
#include "CbcStrategy.hpp"
#include "CbcMessage.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiCuts.hpp"
#include "CbcCountRowCut.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcFeasibilityBase.hpp"
#include "CbcFathom.hpp"
#include "CbcFullNodeInfo.hpp"
#ifdef COIN_HAS_NTY
#include "CbcSymmetry.hpp"
#endif
// include Probing
#include "CglProbing.hpp"
#include "CglGomory.hpp"
#include "CglTwomir.hpp"
// include preprocessing
#include "CglPreProcess.hpp"
#include "CglDuplicateRow.hpp"
#include "CglStored.hpp"
#include "CglClique.hpp"
#include "CglKnapsackCover.hpp"

#include "CoinTime.hpp"
#include "CoinMpsIO.hpp"

#include "CbcCompareActual.hpp"
#include "CbcTree.hpp"
// This may be dummy
#include "CbcThread.hpp"
/* Various functions local to CbcModel.cpp */

typedef struct {
  double useCutoff;
  CbcModel *model;
  int switches;
} rootBundle;
static void *doRootCbcThread(void *voidInfo);

namespace {

//-------------------------------------------------------------------
// Returns the greatest common denominator of two
// positive integers, a and b, found using Euclid's algorithm
//-------------------------------------------------------------------
static int gcd(int a, int b)
{
  int remainder = -1;
  // make sure a<=b (will always remain so)
  if (a > b) {
    // Swap a and b
    int temp = a;
    a = b;
    b = temp;
  }
  // if zero then gcd is nonzero (zero may occur in rhs of packed)
  if (!a) {
    if (b) {
      return b;
    } else {
      printf("**** gcd given two zeros!!\n");
      abort();
    }
  }
  while (remainder) {
    remainder = b % a;
    b = a;
    a = remainder;
  }
  return b;
}

#ifdef CHECK_NODE_FULL

/*
  Routine to verify that tree linkage is correct. The invariant that is tested
  is

  reference count = (number of actual references) + (number of branches left)

  The routine builds a set of paired arrays, info and count, by traversing the
  tree. Each CbcNodeInfo is recorded in info, and the number of times it is
  referenced (via the parent field) is recorded in count. Then a final check is
  made to see if the numberPointingToThis_ field agrees.
*/

void verifyTreeNodes(const CbcTree *branchingTree, const CbcModel &model)

{
  if (model.getNodeCount() == 661)
    return;
  printf("*** CHECKING tree after %d nodes\n", model.getNodeCount());

  int j;
  int nNodes = branchingTree->size();
#define MAXINFO 1000
  int *count = new int[MAXINFO];
  CbcNodeInfo **info = new CbcNodeInfo *[MAXINFO];
  int nInfo = 0;
  /*
      Collect all CbcNodeInfo objects in info, by starting from each live node and
      traversing back to the root. Nodes in the live set should have unexplored
      branches remaining.

      TODO: The `while (nodeInfo)' loop could be made to break on reaching a
    	common ancester (nodeInfo is found in info[k]). Alternatively, the
    	check could change to signal an error if nodeInfo is not found above a
    	common ancestor.
    */
  for (j = 0; j < nNodes; j++) {
    CbcNode *node = branchingTree->nodePointer(j);
    if (!node)
      continue;
    CbcNodeInfo *nodeInfo = node->nodeInfo();
    int change = node->nodeInfo()->numberBranchesLeft();
    assert(change);
    while (nodeInfo) {
      int k;
      for (k = 0; k < nInfo; k++) {
        if (nodeInfo == info[k])
          break;
      }
      if (k == nInfo) {
        assert(nInfo < MAXINFO);
        nInfo++;
        info[k] = nodeInfo;
        count[k] = 0;
      }
      nodeInfo = nodeInfo->parent();
    }
  }
  /*
      Walk the info array. For each nodeInfo, look up its parent in info and
      increment the corresponding count.
    */
  for (j = 0; j < nInfo; j++) {
    CbcNodeInfo *nodeInfo = info[j];
    nodeInfo = nodeInfo->parent();
    if (nodeInfo) {
      int k;
      for (k = 0; k < nInfo; k++) {
        if (nodeInfo == info[k])
          break;
      }
      assert(k < nInfo);
      count[k]++;
    }
  }
  /*
      Walk the info array one more time and check that the invariant holds. The
      number of references (numberPointingToThis()) should equal the sum of the
      number of actual references (held in count[]) plus the number of potential
      references (unexplored branches, numberBranchesLeft()).
    */
  for (j = 0; j < nInfo; j++) {
    CbcNodeInfo *nodeInfo = info[j];
    if (nodeInfo) {
      int k;
      for (k = 0; k < nInfo; k++)
        if (nodeInfo == info[k])
          break;
      printf("Nodeinfo %x - %d left, %d count\n",
        nodeInfo,
        nodeInfo->numberBranchesLeft(),
        nodeInfo->numberPointingToThis());
      assert(nodeInfo->numberPointingToThis() == count[k] + nodeInfo->numberBranchesLeft());
    }
  }

  delete[] count;
  delete[] info;

  return;
}

#endif /* CHECK_NODE_FULL */

#ifdef CHECK_CUT_COUNTS

/*
  Routine to verify that cut reference counts are correct.
*/
void verifyCutCounts(const CbcTree *branchingTree, CbcModel &model)

{
  printf("*** CHECKING cuts after %d nodes\n", model.getNodeCount());

  int j;
  int nNodes = branchingTree->size();

  /*
      cut.tempNumber_ exists for the purpose of doing this verification. Clear it
      in all cuts. We traverse the tree by starting from each live node and working
      back to the root. At each CbcNodeInfo, check for cuts.
    */
  for (j = 0; j < nNodes; j++) {
    CbcNode *node = branchingTree->nodePointer(j);
    CbcNodeInfo *nodeInfo = node->nodeInfo();
    assert(node->nodeInfo()->numberBranchesLeft());
    while (nodeInfo) {
      int k;
      for (k = 0; k < nodeInfo->numberCuts(); k++) {
        CbcCountRowCut *cut = nodeInfo->cuts()[k];
        if (cut)
          cut->tempNumber_ = 0;
      }
      nodeInfo = nodeInfo->parent();
    }
  }
  /*
      Walk the live set again, this time collecting the list of cuts in use at each
      node. addCuts1 will collect the cuts in model.addedCuts_. Take into account
      that when we recreate the basis for a node, we compress out the slack cuts.
    */
  for (j = 0; j < nNodes; j++) {
    CoinWarmStartBasis *debugws = model.getEmptyBasis();
    CbcNode *node = branchingTree->nodePointer(j);
    CbcNodeInfo *nodeInfo = node->nodeInfo();
    int change = node->nodeInfo()->numberBranchesLeft();
    printf("Node %d %x (info %x) var %d way %d obj %g", j, node,
      node->nodeInfo(), node->columnNumber(), node->way(),
      node->objectiveValue());

    model.addCuts1(node, debugws);

    int i;
    int numberRowsAtContinuous = model.numberRowsAtContinuous();
    CbcCountRowCut **addedCuts = model.addedCuts();
    for (i = 0; i < model.currentNumberCuts(); i++) {
      CoinWarmStartBasis::Status status = debugws->getArtifStatus(i + numberRowsAtContinuous);
      if (status != CoinWarmStartBasis::basic && addedCuts[i]) {
        addedCuts[i]->tempNumber_ += change;
      }
    }

    while (nodeInfo) {
      nodeInfo = nodeInfo->parent();
      if (nodeInfo)
        printf(" -> %x", nodeInfo);
    }
    printf("\n");
    delete debugws;
  }
  /*
      The moment of truth: We've tallied up the references by direct scan of the  search tree. Check for agreement with the count in the cut.

      TODO: Rewrite to check and print mismatch only when tempNumber_ == 0?
    */
  for (j = 0; j < nNodes; j++) {
    CbcNode *node = branchingTree->nodePointer(j);
    CbcNodeInfo *nodeInfo = node->nodeInfo();
    while (nodeInfo) {
      int k;
      for (k = 0; k < nodeInfo->numberCuts(); k++) {
        CbcCountRowCut *cut = nodeInfo->cuts()[k];
        if (cut && cut->tempNumber_ >= 0) {
          if (cut->tempNumber_ != cut->numberPointingToThis())
            printf("mismatch %x %d %x %d %d\n", nodeInfo, k,
              cut, cut->tempNumber_, cut->numberPointingToThis());
          else
            printf("   match %x %d %x %d %d\n", nodeInfo, k,
              cut, cut->tempNumber_, cut->numberPointingToThis());
          cut->tempNumber_ = -1;
        }
      }
      nodeInfo = nodeInfo->parent();
    }
  }

  return;
}

#endif /* CHECK_CUT_COUNTS */

#ifdef CHECK_CUT_SIZE

/*
  Routine to verify that cut reference counts are correct.
*/
void verifyCutSize(const CbcTree *branchingTree, CbcModel &model)
{

  int j;
  int nNodes = branchingTree->size();
  int totalCuts = 0;

  /*
      cut.tempNumber_ exists for the purpose of doing this verification. Clear it
      in all cuts. We traverse the tree by starting from each live node and working
      back to the root. At each CbcNodeInfo, check for cuts.
    */
  for (j = 0; j < nNodes; j++) {
    CbcNode *node = branchingTree->nodePointer(j);
    CbcNodeInfo *nodeInfo = node->nodeInfo();
    assert(node->nodeInfo()->numberBranchesLeft());
    while (nodeInfo) {
      totalCuts += nodeInfo->numberCuts();
      nodeInfo = nodeInfo->parent();
    }
  }
  printf("*** CHECKING cuts (size) after %d nodes - %d cuts\n", model.getNodeCount(), totalCuts);
  return;
}

#endif /* CHECK_CUT_SIZE */

}

/* End unnamed namespace for CbcModel.cpp */

void CbcModel::analyzeObjective()
/*
  Try to find a minimum change in the objective function. The first scan
  checks that there are no continuous variables with non-zero coefficients,
  and grabs the largest objective coefficient associated with an unfixed
  integer variable. The second scan attempts to scale up the objective
  coefficients to a point where they are sufficiently close to integer that
  we can pretend they are integer, and calculate a gcd over the coefficients
  of interest. This will be the minimum increment for the scaled coefficients.
  The final action is to scale the increment back for the original coefficients
  and install it, if it's better than the existing value.

  John's note: We could do better than this.

  John's second note - apologies for changing s to z
*/
{
  const double *objective = getObjCoefficients();
  const double *lower = getColLower();
  const double *upper = getColUpper();
  /*
      Scan continuous and integer variables to see if continuous
      are cover or network with integral rhs.
    */
  double continuousMultiplier = 1.0;
  double *coeffMultiplier = NULL;
  double largestObj = 0.0;
  double smallestObj = COIN_DBL_MAX;
  {
    const double *rowLower = getRowLower();
    const double *rowUpper = getRowUpper();
    int numberRows = solver_->getNumRows();
    double *rhs = new double[numberRows];
    memset(rhs, 0, numberRows * sizeof(double));
    int iColumn;
    int numberColumns = solver_->getNumCols();
    // Column copy of matrix
    int problemType = -1;
    const double *element = solver_->getMatrixByCol()->getElements();
    const int *row = solver_->getMatrixByCol()->getIndices();
    const CoinBigIndex *columnStart = solver_->getMatrixByCol()->getVectorStarts();
    const int *columnLength = solver_->getMatrixByCol()->getVectorLengths();
    int numberInteger = 0;
    int numberIntegerObj = 0;
    int numberGeneralIntegerObj = 0;
    int numberIntegerWeight = 0;
    int numberContinuousObj = 0;
    double cost = COIN_DBL_MAX;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (upper[iColumn] == lower[iColumn]) {
        CoinBigIndex start = columnStart[iColumn];
        CoinBigIndex end = start + columnLength[iColumn];
        for (CoinBigIndex j = start; j < end; j++) {
          int iRow = row[j];
          rhs[iRow] += lower[iColumn] * element[j];
        }
      } else {
        double objValue = objective[iColumn];
        if (solver_->isInteger(iColumn))
          numberInteger++;
        if (objValue) {
          if (!solver_->isInteger(iColumn)) {
            numberContinuousObj++;
          } else {
            largestObj = CoinMax(largestObj, fabs(objValue));
            smallestObj = CoinMin(smallestObj, fabs(objValue));
            numberIntegerObj++;
            if (cost == COIN_DBL_MAX)
              cost = objValue;
            else if (cost != objValue)
              cost = -COIN_DBL_MAX;
            int gap = static_cast< int >(upper[iColumn] - lower[iColumn]);
            if (gap > 1) {
              numberGeneralIntegerObj++;
              numberIntegerWeight += gap;
            }
          }
        }
      }
    }
    int iType = 0;
    if (!numberContinuousObj && numberIntegerObj <= 5 && numberIntegerWeight <= 100 && numberIntegerObj * 3 < numberObjects_ && !parentModel_ && solver_->getNumRows() > 100)
      iType = 1 + 4 + (((moreSpecialOptions_ & 536870912) == 0) ? 2 : 0);
    else if (!numberContinuousObj && numberIntegerObj <= 100 && numberIntegerObj * 5 < numberObjects_ && numberIntegerWeight <= 100 && !parentModel_ && solver_->getNumRows() > 100 && cost != -COIN_DBL_MAX)
      iType = 4 + (((moreSpecialOptions_ & 536870912) == 0) ? 2 : 0);
    else if (!numberContinuousObj && numberIntegerObj <= 100 && numberIntegerObj * 5 < numberObjects_ && !parentModel_ && solver_->getNumRows() > 100 && cost != -COIN_DBL_MAX)
      iType = 8;
    int iTest = getMaximumNodes();
    if (iTest >= 987654320 && iTest < 987654330 && numberObjects_ && !parentModel_) {
      iType = iTest - 987654320;
      printf("Testing %d integer variables out of %d objects (%d integer) have cost of %g - %d continuous\n",
        numberIntegerObj, numberObjects_, numberInteger, cost, numberContinuousObj);
      if (iType == 9)
        exit(77);
      if (numberContinuousObj)
        iType = 0;
    }

    //if (!numberContinuousObj&&(numberIntegerObj<=5||cost!=-COIN_DBL_MAX)&&
    //numberIntegerObj*3<numberObjects_&&!parentModel_&&solver_->getNumRows()>100) {
    if (iType) {
      /*
            A) put high priority on (if none)
            B) create artificial objective (if clp)
            */
      int iPriority = -1;
      for (int i = 0; i < numberObjects_; i++) {
        int k = object_[i]->priority();
        if (iPriority == -1)
          iPriority = k;
        else if (iPriority != k)
          iPriority = -2;
      }
      bool branchOnSatisfied = ((iType & 1) != 0);
      bool createFake = ((iType & 2) != 0);
      bool randomCost = ((iType & 4) != 0);
      if (iPriority >= 0) {
        char general[200];
        if (cost == -COIN_DBL_MAX) {
          sprintf(general, "%d integer variables out of %d objects (%d integer) have costs - high priority",
            numberIntegerObj, numberObjects_, numberInteger);
        } else if (cost == COIN_DBL_MAX) {
          sprintf(general, "No integer variables out of %d objects (%d integer) have costs",
            numberObjects_, numberInteger);
          branchOnSatisfied = false;
        } else {
          sprintf(general, "%d integer variables out of %d objects (%d integer) have cost of %g - high priority",
            numberIntegerObj, numberObjects_, numberInteger, cost);
        }
        messageHandler()->message(CBC_GENERAL,
          messages())
          << general << CoinMessageEol;
        sprintf(general, "branch on satisfied %c create fake objective %c random cost %c",
          branchOnSatisfied ? 'Y' : 'N',
          createFake ? 'Y' : 'N',
          randomCost ? 'Y' : 'N');
        messageHandler()->message(CBC_GENERAL,
          messages())
          << general << CoinMessageEol;
        // switch off clp type branching
        // no ? fastNodeDepth_ = -1;
        int highPriority = (branchOnSatisfied) ? -999 : 100;
        for (int i = 0; i < numberObjects_; i++) {
          CbcSimpleInteger *thisOne = dynamic_cast< CbcSimpleInteger * >(object_[i]);
          object_[i]->setPriority(1000);
          if (thisOne) {
            int iColumn = thisOne->columnNumber();
            if (objective[iColumn])
              thisOne->setPriority(highPriority);
          }
        }
      }
#ifdef COIN_HAS_CLP
      OsiClpSolverInterface *clpSolver
        = dynamic_cast< OsiClpSolverInterface * >(solver_);
      if (clpSolver && createFake) {
        // Create artificial objective to be used when all else fixed
        int numberColumns = clpSolver->getNumCols();
        double *fakeObj = new double[numberColumns];
        // Column copy
        const CoinPackedMatrix *matrixByCol = clpSolver->getMatrixByCol();
        //const double * element = matrixByCol.getElements();
        //const int * row = matrixByCol.getIndices();
        //const CoinBigIndex * columnStart = matrixByCol.getVectorStarts();
        const int *columnLength = matrixByCol->getVectorLengths();
        const double *solution = clpSolver->getColSolution();
#ifdef JJF_ZERO
        int nAtBound = 0;
        for (int i = 0; i < numberColumns; i++) {
          double lowerValue = lower[i];
          double upperValue = upper[i];
          if (clpSolver->isInteger(i)) {
            double lowerValue = lower[i];
            double upperValue = upper[i];
            double value = solution[i];
            if (value < lowerValue + 1.0e-6 || value > upperValue - 1.0e-6)
              nAtBound++;
          }
        }
#endif
        /*
                  Generate a random objective function for problems where the given objective
                  function is not terribly useful. (Nearly feasible, single integer variable,
                  that sort of thing.
                */
        CoinDrand48(true, 1234567);
        for (int i = 0; i < numberColumns; i++) {
          double lowerValue = lower[i];
          double upperValue = upper[i];
          double value = (randomCost) ? ceil((CoinDrand48() + 0.5) * 1000)
                                      : i + 1 + columnLength[i] * 1000;
          value *= 0.001;
          //value += columnLength[i];
          if (lowerValue > -1.0e5 || upperValue < 1.0e5) {
            if (fabs(lowerValue) > fabs(upperValue))
              value = -value;
            if (clpSolver->isInteger(i)) {
              double solValue = solution[i];
              // Better to add in 0.5 or 1.0??
              if (solValue < lowerValue + 1.0e-6)
                value = fabs(value) + 0.5; //fabs(value*1.5);
              else if (solValue > upperValue - 1.0e-6)
                value = -fabs(value) - 0.5; //-fabs(value*1.5);
            }
          } else {
            value = 0.0;
          }
          fakeObj[i] = value;
        }
        // pass to solver
        clpSolver->setFakeObjective(fakeObj);
        delete[] fakeObj;
      }
#endif
    } else if (largestObj < smallestObj * 5.0 && !parentModel_ && !numberContinuousObj && !numberGeneralIntegerObj && numberIntegerObj * 2 < CoinMin(numberColumns, 20)) {
      // up priorities on costed
      int iPriority = -1;
      for (int i = 0; i < numberObjects_; i++) {
        int k = object_[i]->priority();
        if (iPriority == -1)
          iPriority = k;
        else if (iPriority != k)
          iPriority = -2;
      }
      if (iPriority >= 100) {
#if CBC_USEFUL_PRINTING > 1
        printf("Setting variables with obj to high priority\n");
#endif
        for (int i = 0; i < numberObjects_; i++) {
          CbcSimpleInteger *obj = dynamic_cast< CbcSimpleInteger * >(object_[i]);
          if (obj) {
            int iColumn = obj->columnNumber();
            if (objective[iColumn])
              object_[i]->setPriority(iPriority - 1);
          }
        }
      }
    }
    int iRow;
    for (iRow = 0; iRow < numberRows; iRow++) {
      if (rowLower[iRow] > -1.0e20 && fabs(rowLower[iRow] - rhs[iRow] - floor(rowLower[iRow] - rhs[iRow] + 0.5)) > 1.0e-10) {
        continuousMultiplier = 0.0;
        break;
      }
      if (rowUpper[iRow] < 1.0e20 && fabs(rowUpper[iRow] - rhs[iRow] - floor(rowUpper[iRow] - rhs[iRow] + 0.5)) > 1.0e-10) {
        continuousMultiplier = 0.0;
        break;
      }
      // set rhs to limiting value
      if (rowLower[iRow] != rowUpper[iRow]) {
        if (rowLower[iRow] > -1.0e20) {
          if (rowUpper[iRow] < 1.0e20) {
            // no good
            continuousMultiplier = 0.0;
            break;
          } else {
            rhs[iRow] = rowLower[iRow] - rhs[iRow];
            if (problemType < 0)
              problemType = 3; // set cover
            else if (problemType != 3)
              problemType = 4;
          }
        } else {
          rhs[iRow] = rowUpper[iRow] - rhs[iRow];
          if (problemType < 0)
            problemType = 1; // set partitioning <=
          else if (problemType != 1)
            problemType = 4;
        }
      } else {
        rhs[iRow] = rowUpper[iRow] - rhs[iRow];
        if (problemType < 0)
          problemType = 3; // set partitioning ==
        else if (problemType != 2)
          problemType = 2;
      }
      if (fabs(rhs[iRow] - 1.0) > 1.0e-12)
        problemType = 4;
    }
    if (continuousMultiplier) {
      // 1 network, 2 cover, 4 negative cover
      int possible = 7;
      bool unitRhs = true;
      // See which rows could be set cover
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        if (upper[iColumn] > lower[iColumn] + 1.0e-8) {
          CoinBigIndex start = columnStart[iColumn];
          CoinBigIndex end = start + columnLength[iColumn];
          for (CoinBigIndex j = start; j < end; j++) {
            double value = element[j];
            if (value == 1.0) {
            } else if (value == -1.0) {
              rhs[row[j]] = -0.5;
            } else {
              rhs[row[j]] = -COIN_DBL_MAX;
            }
          }
        }
      }
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        if (upper[iColumn] > lower[iColumn] + 1.0e-8) {
          if (!isInteger(iColumn)) {
            CoinBigIndex start = columnStart[iColumn];
            CoinBigIndex end = start + columnLength[iColumn];
            double rhsValue = 0.0;
            // 1 all ones, -1 all -1s, 2 all +- 1, 3 no good
            int type = 0;
            for (CoinBigIndex j = start; j < end; j++) {
              double value = element[j];
              if (fabs(value) != 1.0) {
                type = 3;
                break;
              } else if (value == 1.0) {
                if (!type)
                  type = 1;
                else if (type != 1)
                  type = 2;
              } else {
                if (!type)
                  type = -1;
                else if (type != -1)
                  type = 2;
              }
              int iRow = row[j];
              if (rhs[iRow] == -COIN_DBL_MAX) {
                type = 3;
                break;
              } else if (rhs[iRow] == -0.5) {
                // different values
                unitRhs = false;
              } else if (rhsValue) {
                if (rhsValue != rhs[iRow])
                  unitRhs = false;
              } else {
                rhsValue = rhs[iRow];
              }
            }
            // if no elements OK
            if (type == 3) {
              // no good
              possible = 0;
              break;
            } else if (type == 2) {
              if (end - start > 2) {
                // no good
                possible = 0;
                break;
              } else {
                // only network
                possible &= 1;
                if (!possible)
                  break;
              }
            } else if (type == 1) {
              // only cover
              possible &= 2;
              if (!possible)
                break;
            } else if (type == -1) {
              // only negative cover
              possible &= 4;
              if (!possible)
                break;
            }
          }
        }
      }
      if ((possible == 2 || possible == 4) && !unitRhs) {
#if COIN_DEVELOP > 1
        printf("XXXXXX Continuous all +1 but different rhs\n");
#endif
        possible = 0;
      }
      // may be all integer
      if (possible != 7) {
        if (!possible)
          continuousMultiplier = 0.0;
        else if (possible == 1)
          continuousMultiplier = 1.0;
        else
          continuousMultiplier = 0.0; // 0.5 was incorrect;
#if COIN_DEVELOP > 1
        if (continuousMultiplier)
          printf("XXXXXX multiplier of %g\n", continuousMultiplier);
#endif
        if (continuousMultiplier == 0.5) {
          coeffMultiplier = new double[numberColumns];
          bool allOne = true;
          for (iColumn = 0; iColumn < numberColumns; iColumn++) {
            coeffMultiplier[iColumn] = 1.0;
            if (upper[iColumn] > lower[iColumn] + 1.0e-8) {
              if (!isInteger(iColumn)) {
                CoinBigIndex start = columnStart[iColumn];
                int iRow = row[start];
                double value = rhs[iRow];
                assert(value >= 0.0);
                if (value != 0.0 && value != 1.0)
                  allOne = false;
                coeffMultiplier[iColumn] = 0.5 * value;
              }
            }
          }
          if (allOne) {
            // back to old way
            delete[] coeffMultiplier;
            coeffMultiplier = NULL;
          }
        }
      } else {
        // all integer
        problemType_ = problemType;
#if COIN_DEVELOP > 1
        printf("Problem type is %d\n", problemType_);
#endif
      }
    }

    // But try again
    if (continuousMultiplier < 1.0) {
      memset(rhs, 0, numberRows * sizeof(double));
      int *count = new int[numberRows];
      memset(count, 0, numberRows * sizeof(int));
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        CoinBigIndex start = columnStart[iColumn];
        CoinBigIndex end = start + columnLength[iColumn];
        if (upper[iColumn] == lower[iColumn]) {
          for (CoinBigIndex j = start; j < end; j++) {
            int iRow = row[j];
            rhs[iRow] += lower[iColumn] * element[j];
          }
        } else if (solver_->isInteger(iColumn)) {
          for (CoinBigIndex j = start; j < end; j++) {
            int iRow = row[j];
            if (fabs(element[j] - floor(element[j] + 0.5)) > 1.0e-10)
              rhs[iRow] = COIN_DBL_MAX;
          }
        } else {
          for (CoinBigIndex j = start; j < end; j++) {
            int iRow = row[j];
            count[iRow]++;
            if (fabs(element[j]) != 1.0)
              rhs[iRow] = COIN_DBL_MAX;
          }
        }
      }
      // now look at continuous
      bool allGood = true;
      double direction = solver_->getObjSense();
      int numberObj = 0;
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        if (upper[iColumn] > lower[iColumn]) {
          double objValue = objective[iColumn] * direction;
          if (objValue && !solver_->isInteger(iColumn)) {
            numberObj++;
            CoinBigIndex start = columnStart[iColumn];
            CoinBigIndex end = start + columnLength[iColumn];
            if (objValue > 0.0) {
              // wants to be as low as possible
              if (lower[iColumn] < -1.0e10 || fabs(lower[iColumn] - floor(lower[iColumn] + 0.5)) > 1.0e-10) {
                allGood = false;
                break;
              } else if (upper[iColumn] < 1.0e10 && fabs(upper[iColumn] - floor(upper[iColumn] + 0.5)) > 1.0e-10) {
                allGood = false;
                break;
              }
              bool singletonRow = true;
              bool equality = false;
              for (CoinBigIndex j = start; j < end; j++) {
                int iRow = row[j];
                if (count[iRow] > 1)
                  singletonRow = false;
                else if (rowLower[iRow] == rowUpper[iRow])
                  equality = true;
                double rhsValue = rhs[iRow];
                double lowerValue = rowLower[iRow];
                double upperValue = rowUpper[iRow];
                if (rhsValue < 1.0e20) {
                  if (lowerValue > -1.0e20)
                    lowerValue -= rhsValue;
                  if (upperValue < 1.0e20)
                    upperValue -= rhsValue;
                }
                if (fabs(rhsValue) > 1.0e20 || fabs(rhsValue - floor(rhsValue + 0.5)) > 1.0e-10
                  || fabs(element[j]) != 1.0) {
                  // no good
                  allGood = false;
                  break;
                }
                if (element[j] > 0.0) {
                  if (lowerValue > -1.0e20 && fabs(lowerValue - floor(lowerValue + 0.5)) > 1.0e-10) {
                    // no good
                    allGood = false;
                    break;
                  }
                } else {
                  if (upperValue < 1.0e20 && fabs(upperValue - floor(upperValue + 0.5)) > 1.0e-10) {
                    // no good
                    allGood = false;
                    break;
                  }
                }
              }
              if (!singletonRow && end > start + 1 && !equality)
                allGood = false;
              if (!allGood)
                break;
            } else {
              // wants to be as high as possible
              if (upper[iColumn] > 1.0e10 || fabs(upper[iColumn] - floor(upper[iColumn] + 0.5)) > 1.0e-10) {
                allGood = false;
                break;
              } else if (lower[iColumn] > -1.0e10 && fabs(lower[iColumn] - floor(lower[iColumn] + 0.5)) > 1.0e-10) {
                allGood = false;
                break;
              }
              bool singletonRow = true;
              bool equality = false;
              for (CoinBigIndex j = start; j < end; j++) {
                int iRow = row[j];
                if (count[iRow] > 1)
                  singletonRow = false;
                else if (rowLower[iRow] == rowUpper[iRow])
                  equality = true;
                double rhsValue = rhs[iRow];
                double lowerValue = rowLower[iRow];
                double upperValue = rowUpper[iRow];
                if (rhsValue < 1.0e20) {
                  if (lowerValue > -1.0e20)
                    lowerValue -= rhsValue;
                  if (upperValue < 1.0e20)
                    upperValue -= rhsValue;
                }
                if (fabs(rhsValue) > 1.0e20 || fabs(rhsValue - floor(rhsValue + 0.5)) > 1.0e-10
                  || fabs(element[j]) != 1.0) {
                  // no good
                  allGood = false;
                  break;
                }
                if (element[j] < 0.0) {
                  if (lowerValue > -1.0e20 && fabs(lowerValue - floor(lowerValue + 0.5)) > 1.0e-10) {
                    // no good
                    allGood = false;
                    break;
                  }
                } else {
                  if (upperValue < 1.0e20 && fabs(upperValue - floor(upperValue + 0.5)) > 1.0e-10) {
                    // no good
                    allGood = false;
                    break;
                  }
                }
              }
              if (!singletonRow && end > start + 1 && !equality)
                allGood = false;
              if (!allGood)
                break;
            }
          }
        }
      }
      delete[] count;
      if (allGood) {
#if COIN_DEVELOP > 1
        if (numberObj)
          printf("YYYY analysis says all continuous with costs will be integer\n");
#endif
        continuousMultiplier = 1.0;
      }
    }
    delete[] rhs;
  }
  /*
      Take a first scan to see if there are unfixed continuous variables in the
      objective.  If so, the minimum objective change could be arbitrarily small.
      Also pick off the maximum coefficient of an unfixed integer variable.

      If the objective is found to contain only integer variables, set the
      fathoming discipline to strict.
    */
  double maximumCost = 0.0;
  //double trueIncrement=0.0;
  int iColumn;
  int numberColumns = getNumCols();
  double scaleFactor = 1.0; // due to rhs etc
  /*
      Original model did not have integer bounds.
    */
  if ((specialOptions_ & 65536) == 0) {
    /* be on safe side (later look carefully as may be able to
           to get 0.5 say if bounds are multiples of 0.5 */
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (upper[iColumn] > lower[iColumn] + 1.0e-8) {
        double value;
        value = fabs(lower[iColumn]);
        if (floor(value + 0.5) != value) {
          scaleFactor = CoinMin(scaleFactor, 0.5);
          if (floor(2.0 * value + 0.5) != 2.0 * value) {
            scaleFactor = CoinMin(scaleFactor, 0.25);
            if (floor(4.0 * value + 0.5) != 4.0 * value) {
              scaleFactor = 0.0;
            }
          }
        }
        value = fabs(upper[iColumn]);
        if (floor(value + 0.5) != value) {
          scaleFactor = CoinMin(scaleFactor, 0.5);
          if (floor(2.0 * value + 0.5) != 2.0 * value) {
            scaleFactor = CoinMin(scaleFactor, 0.25);
            if (floor(4.0 * value + 0.5) != 4.0 * value) {
              scaleFactor = 0.0;
            }
          }
        }
      }
    }
  }
  bool possibleMultiple = continuousMultiplier != 0.0 && scaleFactor != 0.0;
  if (possibleMultiple) {
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (upper[iColumn] > lower[iColumn] + 1.0e-8) {
        maximumCost = CoinMax(maximumCost, fabs(objective[iColumn]));
      }
    }
  }
  setIntParam(CbcModel::CbcFathomDiscipline, possibleMultiple);
  /*
      If a nontrivial increment is possible, try and figure it out. We're looking
      for gcd(c<j>) for all c<j> that are coefficients of unfixed integer
      variables. Since the c<j> might not be integers, try and inflate them
      sufficiently that they look like integers (and we'll deflate the gcd
      later).

      2520.0 is used as it is a nice multiple of 2,3,5,7
    */
  if (possibleMultiple && maximumCost) {
    int increment = 0;
    double multiplier = 2520.0;
    while (10.0 * multiplier * maximumCost < 1.0e8)
      multiplier *= 10.0;
    int bigIntegers = 0; // Count of large costs which are integer
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (upper[iColumn] > lower[iColumn] + 1.0e-8) {
        double objValue = fabs(objective[iColumn]);
        if (!isInteger(iColumn)) {
          if (!coeffMultiplier)
            objValue *= continuousMultiplier;
          else
            objValue *= coeffMultiplier[iColumn];
        }
        if (objValue) {
          double value = objValue * multiplier;
          if (value < 2.1e9) {
            int nearest = static_cast< int >(floor(value + 0.5));
            if (fabs(value - floor(value + 0.5)) > 1.0e-8) {
              increment = 0;
              break;
            } else if (!increment) {
              increment = nearest;
            } else {
              increment = gcd(increment, nearest);
            }
          } else {
            // large value - may still be multiple of 1.0
            if (fabs(objValue - floor(objValue + 0.5)) > 1.0e-8) {
              increment = 0;
              break;
            } else {
              bigIntegers++;
            }
          }
        }
      }
    }
    if (coeffMultiplier) {
        delete[] coeffMultiplier;
        coeffMultiplier = NULL;
    }
    /*
          If the increment beats the current value for objective change, install it.
        */
    if (increment) {
      double value = increment;
      double cutoff = getDblParam(CbcModel::CbcCutoffIncrement);
      if (bigIntegers) {
        // allow for 1.0
        increment = gcd(increment, static_cast< int >(multiplier));
        value = increment;
      }
      value /= multiplier;
      value *= scaleFactor;
      //trueIncrement=CoinMax(cutoff,value);;
      if (value * 0.999 > cutoff) {
        messageHandler()->message(CBC_INTEGERINCREMENT,
          messages())
          << value << CoinMessageEol;
        setDblParam(CbcModel::CbcCutoffIncrement, CoinMax(value * 0.999, value - 1.0e-4));
      }
    }
  }

  if (coeffMultiplier)
      delete[] coeffMultiplier;

  return;
}

/*
saveModel called (carved out of) BranchandBound
*/
void CbcModel::saveModel(OsiSolverInterface *saveSolver, double *checkCutoffForRestart, bool *feasible)
{
  if (saveSolver && (specialOptions_ & 32768) != 0) {
    // See if worth trying reduction
    *checkCutoffForRestart = getCutoff();
    bool tryNewSearch = solverCharacteristics_->reducedCostsAccurate() && (*checkCutoffForRestart < 1.0e20);
    int numberColumns = getNumCols();
    if (tryNewSearch) {
#if CBC_USEFUL_PRINTING > 1
      printf("after %d nodes, cutoff %g - looking\n",
        numberNodes_, getCutoff());
#endif
      saveSolver->resolve();
      double direction = saveSolver->getObjSense();
      double gap = *checkCutoffForRestart - saveSolver->getObjValue() * direction;
      double tolerance;
      saveSolver->getDblParam(OsiDualTolerance, tolerance);
      if (gap <= 0.0)
        gap = tolerance;
      gap += 100.0 * tolerance;
      double integerTolerance = getDblParam(CbcIntegerTolerance);

      const double *lower = saveSolver->getColLower();
      const double *upper = saveSolver->getColUpper();
      const double *solution = saveSolver->getColSolution();
      const double *reducedCost = saveSolver->getReducedCost();

      int numberFixed = 0;
      int numberFixed2 = 0;
      for (int i = 0; i < numberIntegers_; i++) {
        int iColumn = integerVariable_[i];
        double djValue = direction * reducedCost[iColumn];
        if (upper[iColumn] - lower[iColumn] > integerTolerance) {
          if (solution[iColumn] < lower[iColumn] + integerTolerance && djValue > gap) {
            saveSolver->setColUpper(iColumn, lower[iColumn]);
            numberFixed++;
          } else if (solution[iColumn] > upper[iColumn] - integerTolerance && -djValue > gap) {
            saveSolver->setColLower(iColumn, upper[iColumn]);
            numberFixed++;
          }
        } else {
          numberFixed2++;
        }
      }
#ifdef COIN_DEVELOP
      /*
              We're debugging. (specialOptions 1)
            */
      if ((specialOptions_ & 1) != 0) {
        const OsiRowCutDebugger *debugger = saveSolver->getRowCutDebugger();
        if (debugger) {
          printf("Contains optimal\n");
          OsiSolverInterface *temp = saveSolver->clone();
          const double *solution = debugger->optimalSolution();
          const double *lower = temp->getColLower();
          const double *upper = temp->getColUpper();
          int n = temp->getNumCols();
          for (int i = 0; i < n; i++) {
            if (temp->isInteger(i)) {
              double value = floor(solution[i] + 0.5);
              assert(value >= lower[i] && value <= upper[i]);
              temp->setColLower(i, value);
              temp->setColUpper(i, value);
            }
          }
          temp->writeMps("reduced_fix");
          delete temp;
          saveSolver->writeMps("reduced");
        } else {
          abort();
        }
      }
      printf("Restart could fix %d integers (%d already fixed)\n",
        numberFixed + numberFixed2, numberFixed2);
#endif
      numberFixed += numberFixed2;
      if (numberFixed * 20 < numberColumns)
        tryNewSearch = false;
    }
    if (tryNewSearch) {
      // back to solver without cuts?
      OsiSolverInterface *solver2 = continuousSolver_->clone();
      const double *lower = saveSolver->getColLower();
      const double *upper = saveSolver->getColUpper();
      for (int i = 0; i < numberIntegers_; i++) {
        int iColumn = integerVariable_[i];
        solver2->setColLower(iColumn, lower[iColumn]);
        solver2->setColUpper(iColumn, upper[iColumn]);
      }
      // swap
      delete saveSolver;
      saveSolver = solver2;
      double *newSolution = new double[numberColumns];
      double objectiveValue = *checkCutoffForRestart;
      CbcSerendipity heuristic(*this);
      if (bestSolution_)
        heuristic.setInputSolution(bestSolution_, bestObjective_);
      heuristic.setFractionSmall(0.9);
      heuristic.setFeasibilityPumpOptions(1008013);
      // Use numberNodes to say how many are original rows
      heuristic.setNumberNodes(continuousSolver_->getNumRows());
#ifdef COIN_DEVELOP
      if (continuousSolver_->getNumRows() < saveSolver->getNumRows())
        printf("%d rows added ZZZZZ\n",
          solver_->getNumRows() - continuousSolver_->getNumRows());
#endif
      int returnCode = heuristic.smallBranchAndBound(saveSolver,
        -1, newSolution,
        objectiveValue,
        *checkCutoffForRestart, "Reduce");
      if (returnCode < 0) {
#ifdef COIN_DEVELOP
        printf("Restart - not small enough to do search after fixing\n");
#endif
        delete[] newSolution;
      } else {
        if ((returnCode & 1) != 0) {
          // increment number of solutions so other heuristics can test
          numberSolutions_++;
          numberHeuristicSolutions_++;
          lastHeuristic_ = NULL;
          setBestSolution(CBC_ROUNDING, objectiveValue, newSolution);
        }
        delete[] newSolution;
        *feasible = false; // stop search
      }
#if 0 // probably not needed def CBC_THREAD
            if (master_) {
                lockThread();
                if (parallelMode() > 0) {
                    while (master_->waitForThreadsInTree(0)) {
                        lockThread();
                        double dummyBest;
                        tree_->cleanTree(this, -COIN_DBL_MAX, dummyBest) ;
                        //unlockThread();
                    }
                }
                master_->waitForThreadsInTree(2);
                delete master_;
                master_ = NULL;
                masterThread_ = NULL;
            }
#endif
    }
  }
}
/*
Adds integers, called from BranchandBound()
*/
void CbcModel::AddIntegers()
{
  int numberColumns = continuousSolver_->getNumCols();
  int numberRows = continuousSolver_->getNumRows();
  int numberOriginalIntegers = numberIntegers_;
  int *del = new int[CoinMax(numberColumns, numberRows)];
  int *original = new int[numberColumns];
  char *possibleRow = new char[numberRows];
  {
    const CoinPackedMatrix *rowCopy = continuousSolver_->getMatrixByRow();
    const int *column = rowCopy->getIndices();
    const int *rowLength = rowCopy->getVectorLengths();
    const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
    const double *rowLower = continuousSolver_->getRowLower();
    const double *rowUpper = continuousSolver_->getRowUpper();
    const double *element = rowCopy->getElements();
    for (int i = 0; i < numberRows; i++) {
      int nLeft = 0;
      bool possible = false;
      if (rowLower[i] < -1.0e20) {
        double value = rowUpper[i];
        if (fabs(value - floor(value + 0.5)) < 1.0e-8)
          possible = true;
      } else if (rowUpper[i] > 1.0e20) {
        double value = rowLower[i];
        if (fabs(value - floor(value + 0.5)) < 1.0e-8)
          possible = true;
      } else {
        double value = rowUpper[i];
        if (rowLower[i] == rowUpper[i] && fabs(value - floor(value + 0.5)) < 1.0e-8)
          possible = true;
      }
      double allSame = (possible) ? 0.0 : -1.0;
      for (CoinBigIndex j = rowStart[i];
           j < rowStart[i] + rowLength[i]; j++) {
        int iColumn = column[j];
        if (continuousSolver_->isInteger(iColumn)) {
          if (fabs(element[j]) != 1.0)
            possible = false;
        } else {
          nLeft++;
          if (!allSame) {
            allSame = fabs(element[j]);
          } else if (allSame > 0.0) {
            if (allSame != fabs(element[j]))
              allSame = -1.0;
          }
        }
      }
      if (nLeft == rowLength[i] && allSame > 0.0)
        possibleRow[i] = 2;
      else if (possible || !nLeft)
        possibleRow[i] = 1;
      else
        possibleRow[i] = 0;
    }
  }
  int nDel = 0;
  for (int i = 0; i < numberColumns; i++) {
    original[i] = i;
    if (continuousSolver_->isInteger(i))
      del[nDel++] = i;
  }
  {
    // we must not exclude current best solution (rounding errors)
    // also not if large values
    const int *row = continuousSolver_->getMatrixByCol()->getIndices();
    const CoinBigIndex *columnStart = continuousSolver_->getMatrixByCol()->getVectorStarts();
    const int *columnLength = continuousSolver_->getMatrixByCol()->getVectorLengths();
    const double *solution = continuousSolver_->getColSolution();
    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (!continuousSolver_->isInteger(iColumn)) {
        double value = bestSolution_ ? bestSolution_[iColumn] : 0.0;
        double value2 = solution[iColumn];
        if (fabs(value - floor(value + 0.5)) > 1.0e-8 || fabs(value2) > 1.0e3) {
          CoinBigIndex start = columnStart[iColumn];
          CoinBigIndex end = start + columnLength[iColumn];
          for (CoinBigIndex j = start; j < end; j++) {
            int iRow = row[j];
            possibleRow[iRow] = 0;
          }
        }
      }
    }
  }
  int nExtra = 0;
  OsiSolverInterface *copy1 = continuousSolver_->clone();
  int nPass = 0;
  while (nDel && nPass < 10) {
    nPass++;
    OsiSolverInterface *copy2 = copy1->clone();
    int nLeft = 0;
    for (int i = 0; i < nDel; i++)
      original[del[i]] = -1;
    for (int i = 0; i < numberColumns; i++) {
      int kOrig = original[i];
      if (kOrig >= 0)
        original[nLeft++] = kOrig;
    }
    assert(nLeft == numberColumns - nDel);
    copy2->deleteCols(nDel, del);
    numberColumns = copy2->getNumCols();
    const CoinPackedMatrix *rowCopy = copy2->getMatrixByRow();
    numberRows = rowCopy->getNumRows();
    const int *column = rowCopy->getIndices();
    const int *rowLength = rowCopy->getVectorLengths();
    const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
    const double *rowLower = copy2->getRowLower();
    const double *rowUpper = copy2->getRowUpper();
    const double *element = rowCopy->getElements();
    const CoinPackedMatrix *columnCopy = copy2->getMatrixByCol();
    const int *columnLength = columnCopy->getVectorLengths();
    nDel = 0;
    // Could do gcd stuff on ones with costs
    for (int i = 0; i < numberRows; i++) {
      if (!rowLength[i]) {
        del[nDel++] = i;
      } else if (possibleRow[i]) {
        if (rowLength[i] == 1) {
          CoinBigIndex k = rowStart[i];
          int iColumn = column[k];
          if (!copy2->isInteger(iColumn)) {
            double mult = 1.0 / fabs(element[k]);
            if (rowLower[i] < -1.0e20) {
              // treat rhs as multiple of 1 unless elements all same
              double value = ((possibleRow[i] == 2) ? rowUpper[i] : 1.0) * mult;
              if (fabs(value - floor(value + 0.5)) < 1.0e-8) {
                del[nDel++] = i;
                if (columnLength[iColumn] == 1) {
                  copy2->setInteger(iColumn);
                  int kOrig = original[iColumn];
                  setOptionalInteger(kOrig);
                }
              }
            } else if (rowUpper[i] > 1.0e20) {
              // treat rhs as multiple of 1 unless elements all same
              double value = ((possibleRow[i] == 2) ? rowLower[i] : 1.0) * mult;
              if (fabs(value - floor(value + 0.5)) < 1.0e-8) {
                del[nDel++] = i;
                if (columnLength[iColumn] == 1) {
                  copy2->setInteger(iColumn);
                  int kOrig = original[iColumn];
                  setOptionalInteger(kOrig);
                }
              }
            } else {
              // treat rhs as multiple of 1 unless elements all same
              double value = ((possibleRow[i] == 2) ? rowUpper[i] : 1.0) * mult;
              if (rowLower[i] == rowUpper[i] && fabs(value - floor(value + 0.5)) < 1.0e-8) {
                del[nDel++] = i;
                copy2->setInteger(iColumn);
                int kOrig = original[iColumn];
                setOptionalInteger(kOrig);
              }
            }
          }
        } else {
          // only if all singletons
          bool possible = false;
          if (rowLower[i] < -1.0e20) {
            double value = rowUpper[i];
            if (fabs(value - floor(value + 0.5)) < 1.0e-8)
              possible = true;
          } else if (rowUpper[i] > 1.0e20) {
            double value = rowLower[i];
            if (fabs(value - floor(value + 0.5)) < 1.0e-8)
              possible = true;
          } else {
            double value = rowUpper[i];
            if (rowLower[i] == rowUpper[i] && fabs(value - floor(value + 0.5)) < 1.0e-8)
              possible = true;
          }
          if (possible) {
            for (CoinBigIndex j = rowStart[i];
                 j < rowStart[i] + rowLength[i]; j++) {
              int iColumn = column[j];
              if (columnLength[iColumn] != 1 || fabs(element[j]) != 1.0) {
                possible = false;
                break;
              }
            }
            if (possible) {
              for (CoinBigIndex j = rowStart[i];
                   j < rowStart[i] + rowLength[i]; j++) {
                int iColumn = column[j];
                if (!copy2->isInteger(iColumn)) {
                  copy2->setInteger(iColumn);
                  int kOrig = original[iColumn];
                  setOptionalInteger(kOrig);
                }
              }
              del[nDel++] = i;
            }
          }
        }
      }
    }
    if (nDel) {
      copy2->deleteRows(nDel, del);
      // pack down possible
      int n = 0;
      for (int i = 0; i < nDel; i++)
        possibleRow[del[i]] = -1;
      for (int i = 0; i < numberRows; i++) {
        if (possibleRow[i] >= 0)
          possibleRow[n++] = possibleRow[i];
      }
    }
    if (nDel != numberRows) {
      nDel = 0;
      for (int i = 0; i < numberColumns; i++) {
        if (copy2->isInteger(i)) {
          del[nDel++] = i;
          nExtra++;
        }
      }
    } else {
      nDel = 0;
    }
    delete copy1;
    copy1 = copy2->clone();
    delete copy2;
  }
  // See if what's left is a network
  bool couldBeNetwork = false;
  if (copy1->getNumRows() && copy1->getNumCols()) {
#ifdef COIN_HAS_CLP
    OsiClpSolverInterface *clpSolver
      = dynamic_cast< OsiClpSolverInterface * >(copy1);
    if (false && clpSolver) {
      numberRows = clpSolver->getNumRows();
      char *rotate = new char[numberRows];
      int n = clpSolver->getModelPtr()->findNetwork(rotate, 1.0);
      delete[] rotate;
#if CBC_USEFUL_PRINTING > 1
      printf("INTA network %d rows out of %d\n", n, numberRows);
#endif
      if (CoinAbs(n) == numberRows) {
        couldBeNetwork = true;
        for (int i = 0; i < numberRows; i++) {
          if (!possibleRow[i]) {
            couldBeNetwork = false;
#if CBC_USEFUL_PRINTING > 1
            printf("but row %d is bad\n", i);
#endif
            break;
          }
        }
      }
    } else
#endif
    {
      numberColumns = copy1->getNumCols();
      numberRows = copy1->getNumRows();
      const double *rowLower = copy1->getRowLower();
      const double *rowUpper = copy1->getRowUpper();
      couldBeNetwork = true;
      for (int i = 0; i < numberRows; i++) {
        if (rowLower[i] > -1.0e20 && fabs(rowLower[i] - floor(rowLower[i] + 0.5)) > 1.0e-12) {
          couldBeNetwork = false;
          break;
        }
        if (rowUpper[i] < 1.0e20 && fabs(rowUpper[i] - floor(rowUpper[i] + 0.5)) > 1.0e-12) {
          couldBeNetwork = false;
          break;
        }
        if (possibleRow[i] == 0) {
          couldBeNetwork = false;
          break;
        }
      }
      if (couldBeNetwork) {
        const CoinPackedMatrix *matrixByCol = copy1->getMatrixByCol();
        const double *element = matrixByCol->getElements();
        //const int * row = matrixByCol->getIndices();
        const CoinBigIndex *columnStart = matrixByCol->getVectorStarts();
        const int *columnLength = matrixByCol->getVectorLengths();
        for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
          CoinBigIndex start = columnStart[iColumn];
          CoinBigIndex end = start + columnLength[iColumn];
          if (end > start + 2) {
            couldBeNetwork = false;
            break;
          }
          int type = 0;
          for (CoinBigIndex j = start; j < end; j++) {
            double value = element[j];
            if (fabs(value) != 1.0) {
              couldBeNetwork = false;
              break;
            } else if (value == 1.0) {
              if ((type & 1) == 0)
                type |= 1;
              else
                type = 7;
            } else if (value == -1.0) {
              if ((type & 2) == 0)
                type |= 2;
              else
                type = 7;
            }
          }
          if (type > 3) {
            couldBeNetwork = false;
            break;
          }
        }
      }
    }
  }
  if (couldBeNetwork) {
    for (int i = 0; i < numberColumns; i++)
      setOptionalInteger(original[i]);
  }
  if (nExtra || couldBeNetwork) {
    numberColumns = copy1->getNumCols();
    numberRows = copy1->getNumRows();
    if (!numberColumns || !numberRows) {
      int numberColumns = solver_->getNumCols();
      for (int i = 0; i < numberColumns; i++)
        assert(solver_->isInteger(i));
    }
#if CBC_USEFUL_PRINTING > 1
    if (couldBeNetwork || nExtra)
      printf("INTA %d extra integers, %d left%s\n", nExtra,
        numberColumns,
        couldBeNetwork ? ", all network" : "");
#endif
    findIntegers(true, 2);
    convertToDynamic();
  }
#if CBC_USEFUL_PRINTING > 1
  if (!couldBeNetwork && copy1->getNumCols() && copy1->getNumRows()) {
    printf("INTA %d rows and %d columns remain\n",
      copy1->getNumRows(), copy1->getNumCols());
    if (copy1->getNumCols() < 200) {
      copy1->writeMps("moreint");
      printf("INTA Written remainder to moreint.mps.gz %d rows %d cols\n",
        copy1->getNumRows(), copy1->getNumCols());
    }
  }
#endif
  delete copy1;
  delete[] del;
  delete[] original;
  delete[] possibleRow;
  // double check increment
  analyzeObjective();
  // If any changes - tell code
  if (numberOriginalIntegers < numberIntegers_)
    synchronizeModel();
}
/**
  \todo
  Normally, it looks like we enter here from command dispatch in the main
  routine, after calling the solver for an initial solution
  (CbcModel::initialSolve, which simply calls the solver's initialSolve
  routine.) The first thing we do is call resolve. Presumably there are
  circumstances where this is nontrivial? There's also a call from
  CbcModel::originalModel (tied up with integer presolve), which should be
  checked.

*/

#ifdef CONFLICT_CUTS
#if PRINT_CONFLICT == 1
static int numberConflictCuts = 0;
static int lastNumberConflictCuts = 0;
static double lengthConflictCuts = 0.0;
#endif
#endif
/*
  The overall flow can be divided into three stages:
    * Prep: Check that the lp relaxation remains feasible at the root. If so,
      do all the setup for B&C.
    * Process the root node: Generate cuts, apply heuristics, and in general do
      the best we can to resolve the problem without B&C.
    * Do B&C search until we hit a limit or exhaust the search tree.

  Keep in mind that in general there is no node in the search tree that
  corresponds to the active subproblem. The active subproblem is represented
  by the current state of the model,  of the solver, and of the constraint
  system held by the solver.
*/
#ifdef COIN_HAS_CPX
#include "OsiCpxSolverInterface.hpp"
#include "cplex.h"
#endif
void CbcModel::branchAndBound(int doStatistics)

{
  if (!parentModel_) {
    /*
	Capture a time stamp before we start (unless set).
      */
    if (!dblParam_[CbcStartSeconds]) {
      if (!useElapsedTime())
        dblParam_[CbcStartSeconds] = CoinCpuTime();
      else
        dblParam_[CbcStartSeconds] = CoinGetTimeOfDay();
    }
  }
  dblParam_[CbcSmallestChange] = COIN_DBL_MAX;
  dblParam_[CbcSumChange] = 0.0;
  dblParam_[CbcLargestChange] = 0.0;
  intParam_[CbcNumberBranches] = 0;
  double lastBestPossibleObjective = -COIN_DBL_MAX;
  // when to check for restart
  int nextCheckRestart = 50;
  // Force minimization !!!!
  bool flipObjective = (solver_->getObjSense() < 0.0);
  if (flipObjective)
    flipModel();
  dblParam_[CbcOptimizationDirection] = 1.0; // was solver_->getObjSense();
  strongInfo_[0] = 0;
  strongInfo_[1] = 0;
  strongInfo_[2] = 0;
  strongInfo_[3] = 0;
  strongInfo_[4] = 0;
  strongInfo_[5] = 0;
  strongInfo_[6] = 0;
  numberStrongIterations_ = 0;
  currentNode_ = NULL;
  // See if should do cuts old way
  if (parallelMode() < 0) {
    specialOptions_ |= 4096 + 8192;
  } else if (parallelMode() > 0) {
    specialOptions_ |= 4096;
  }
  int saveMoreSpecialOptions = moreSpecialOptions_;
  if (dynamic_cast< CbcTreeLocal * >(tree_))
    specialOptions_ |= 4096 + 8192;
#ifdef COIN_HAS_CLP
  {
    OsiClpSolverInterface *clpSolver
      = dynamic_cast< OsiClpSolverInterface * >(solver_);
    if (clpSolver) {
      // pass in disaster handler
      CbcDisasterHandler handler(this);
      clpSolver->passInDisasterHandler(&handler);
      // Initialise solvers seed (unless users says not)
      if ((specialOptions_ & 4194304) == 0)
        clpSolver->getModelPtr()->setRandomSeed(1234567);
#ifdef JJF_ZERO
      // reduce factorization frequency
      int frequency = clpSolver->getModelPtr()->factorizationFrequency();
      clpSolver->getModelPtr()->setFactorizationFrequency(CoinMin(frequency, 120));
#endif
    }
  }
#endif
  // original solver (only set if pre-processing)
  OsiSolverInterface *originalSolver = NULL;
  int numberOriginalObjects = numberObjects_;
  OsiObject **originalObject = NULL;
  // Save whether there were any objects
  bool noObjects = (numberObjects_ == 0);
  // Set up strategies
  /*
      See if the user has supplied a strategy object and deal with it if present.
      The call to setupOther will set numberStrong_ and numberBeforeTrust_, and
      perform integer preprocessing, if requested.

      We need to hang on to a pointer to solver_. setupOther will assign a
      preprocessed solver to model, but will instruct assignSolver not to trash the
      existing one.
    */
  if (strategy_) {
    // May do preprocessing
    originalSolver = solver_;
    strategy_->setupOther(*this);
    if (strategy_->preProcessState()) {
      // pre-processing done
      if (strategy_->preProcessState() < 0) {
        // infeasible (or unbounded)
        status_ = 0;
        if (!solver_->isProvenDualInfeasible()) {
          handler_->message(CBC_INFEAS, messages_) << CoinMessageEol;
          secondaryStatus_ = 1;
        } else {
          handler_->message(CBC_UNBOUNDED,
            messages_)
            << CoinMessageEol;
          secondaryStatus_ = 7;
        }
        originalContinuousObjective_ = COIN_DBL_MAX;
        if (flipObjective)
          flipModel();
        return;
      } else if (numberObjects_ && object_) {
        numberOriginalObjects = numberObjects_;
        // redo sequence
        numberIntegers_ = 0;
        int numberColumns = getNumCols();
        int nOrig = originalSolver->getNumCols();
        CglPreProcess *process = strategy_->process();
        assert(process);
        const int *originalColumns = process->originalColumns();
        // allow for cliques etc
        nOrig = CoinMax(nOrig, originalColumns[numberColumns - 1] + 1);
        // try and redo debugger
        OsiRowCutDebugger *debugger = const_cast< OsiRowCutDebugger * >(solver_->getRowCutDebuggerAlways());
        if (debugger) {
          if (numberColumns <= debugger->numberColumns())
            debugger->redoSolution(numberColumns, originalColumns);
          else
            debugger = NULL; // no idea how to handle (SOS?)
        }
        // User-provided solution might have been best. Synchronise.
        if (bestSolution_) {
          // need to redo - in case no better found in BAB
          // just get integer part right
          for (int i = 0; i < numberColumns; i++) {
            int jColumn = originalColumns[i];
            bestSolution_[i] = bestSolution_[jColumn];
          }
        }
        originalObject = object_;
        // object number or -1
        int *temp = new int[nOrig];
        int iColumn;
        for (iColumn = 0; iColumn < nOrig; iColumn++)
          temp[iColumn] = -1;
        int iObject;
        int nNonInt = 0;
        for (iObject = 0; iObject < numberOriginalObjects; iObject++) {
          iColumn = originalObject[iObject]->columnNumber();
          if (iColumn < 0) {
            nNonInt++;
          } else {
            temp[iColumn] = iObject;
          }
        }
        int numberNewIntegers = 0;
        int numberOldIntegers = 0;
        int numberOldOther = 0;
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          int jColumn = originalColumns[iColumn];
          if (temp[jColumn] >= 0) {
            int iObject = temp[jColumn];
            CbcSimpleInteger *obj = dynamic_cast< CbcSimpleInteger * >(originalObject[iObject]);
            if (obj)
              numberOldIntegers++;
            else
              numberOldOther++;
          } else if (isInteger(iColumn)) {
            numberNewIntegers++;
          }
        }
        /*
                  Allocate an array to hold the indices of the integer variables.
                  Make a large enough array for all objects
                */
        numberObjects_ = numberNewIntegers + numberOldIntegers + numberOldOther + nNonInt;
        object_ = new OsiObject *[numberObjects_];
        delete[] integerVariable_;
        integerVariable_ = new int[numberNewIntegers + numberOldIntegers];
        /*
                  Walk the variables again, filling in the indices and creating objects for
                  the integer variables. Initially, the objects hold the index and upper &
                  lower bounds.
                */
        numberIntegers_ = 0;
        int n = originalColumns[numberColumns - 1] + 1;
        int *backward = new int[n];
        int i;
        for (i = 0; i < n; i++)
          backward[i] = -1;
        for (i = 0; i < numberColumns; i++)
          backward[originalColumns[i]] = i;
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          int jColumn = originalColumns[iColumn];
          if (temp[jColumn] >= 0) {
            int iObject = temp[jColumn];
            CbcSimpleInteger *obj = dynamic_cast< CbcSimpleInteger * >(originalObject[iObject]);
            if (obj) {
              object_[numberIntegers_] = originalObject[iObject]->clone();
              // redo ids etc
              //object_[numberIntegers_]->resetSequenceEtc(numberColumns,originalColumns);
              object_[numberIntegers_]->resetSequenceEtc(numberColumns, backward);
              integerVariable_[numberIntegers_++] = iColumn;
            }
          } else if (isInteger(iColumn)) {
            object_[numberIntegers_] = new CbcSimpleInteger(this, iColumn);
            integerVariable_[numberIntegers_++] = iColumn;
          }
        }
        delete[] backward;
        numberObjects_ = numberIntegers_;
        // Now append other column stuff
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          int jColumn = originalColumns[iColumn];
          if (temp[jColumn] >= 0) {
            int iObject = temp[jColumn];
            CbcSimpleInteger *obj = dynamic_cast< CbcSimpleInteger * >(originalObject[iObject]);
            if (!obj) {
              object_[numberObjects_] = originalObject[iObject]->clone();
              // redo ids etc
              CbcObject *obj = dynamic_cast< CbcObject * >(object_[numberObjects_]);
              assert(obj);
              obj->redoSequenceEtc(this, numberColumns, originalColumns);
              numberObjects_++;
            }
          }
        }
        // now append non column stuff
        for (iObject = 0; iObject < numberOriginalObjects; iObject++) {
          iColumn = originalObject[iObject]->columnNumber();
          if (iColumn < 0) {
            // already has column numbers changed
            object_[numberObjects_] = originalObject[iObject]->clone();
#ifdef JJF_ZERO
            // redo ids etc
            CbcObject *obj = dynamic_cast< CbcObject * >(object_[numberObjects_]);
            assert(obj);
            obj->redoSequenceEtc(this, numberColumns, originalColumns);
#endif
            numberObjects_++;
          }
        }
        delete[] temp;
        if (!numberObjects_)
          handler_->message(CBC_NOINT, messages_) << CoinMessageEol;
      } else {
        int numberColumns = getNumCols();
        CglPreProcess *process = strategy_->process();
        assert(process);
        const int *originalColumns = process->originalColumns();
        // try and redo debugger
        OsiRowCutDebugger *debugger = const_cast< OsiRowCutDebugger * >(solver_->getRowCutDebuggerAlways());
        if (debugger)
          debugger->redoSolution(numberColumns, originalColumns);
      }
    } else {
      //no preprocessing
      originalSolver = NULL;
    }
    strategy_->setupCutGenerators(*this);
    strategy_->setupHeuristics(*this);
    // Set strategy print level to models
    strategy_->setupPrinting(*this, handler_->logLevel());
  }
  eventHappened_ = false;
  CbcEventHandler *eventHandler = getEventHandler();
  if (eventHandler)
    eventHandler->setModel(this);
#define CLIQUE_ANALYSIS
#ifdef CLIQUE_ANALYSIS
  // set up for probing
  // If we're doing clever stuff with cliques, additional info here.
  if (!parentModel_)
    probingInfo_ = new CglTreeProbingInfo(solver_);
  else
    probingInfo_ = NULL;
#else
  probingInfo_ = NULL;
#endif

  // Try for dominated columns
  if ((specialOptions_ & 64) != 0) {
    CglDuplicateRow dupcuts(solver_);
    dupcuts.setMode(2);
    CglStored *storedCuts = dupcuts.outDuplicates(solver_);
    if (storedCuts) {
      COIN_DETAIL_PRINT(printf("adding dup cuts\n"));
      addCutGenerator(storedCuts, 1, "StoredCuts from dominated",
        true, false, false, -200);
    }
  }
  if (!nodeCompare_)
    nodeCompare_ = new CbcCompareDefault();
  ;
  // See if hot start wanted
  CbcCompareBase *saveCompare = NULL;
  // User supplied hotstart. Adapt for preprocessing.
  if (hotstartSolution_) {
    if (strategy_ && strategy_->preProcessState() > 0) {
      CglPreProcess *process = strategy_->process();
      assert(process);
      int n = solver_->getNumCols();
      const int *originalColumns = process->originalColumns();
      // columns should be in order ... but
      double *tempS = new double[n];
      for (int i = 0; i < n; i++) {
        int iColumn = originalColumns[i];
        tempS[i] = hotstartSolution_[iColumn];
      }
      delete[] hotstartSolution_;
      hotstartSolution_ = tempS;
      if (hotstartPriorities_) {
        int *tempP = new int[n];
        for (int i = 0; i < n; i++) {
          int iColumn = originalColumns[i];
          tempP[i] = hotstartPriorities_[iColumn];
        }
        delete[] hotstartPriorities_;
        hotstartPriorities_ = tempP;
      }
    }
    saveCompare = nodeCompare_;
    // depth first
    nodeCompare_ = new CbcCompareDepth();
  }
  if (!problemFeasibility_)
    problemFeasibility_ = new CbcFeasibilityBase();
#ifdef CBC_DEBUG
  std::string problemName;
  solver_->getStrParam(OsiProbName, problemName);
  printf("Problem name - %s\n", problemName.c_str());
  solver_->setHintParam(OsiDoReducePrint, false, OsiHintDo, 0);
#endif
  /*
      Assume we're done, and see if we're proven wrong.
    */
  status_ = 0;
  secondaryStatus_ = 0;
  phase_ = 0;
  /*
      Scan the variables, noting the integer variables. Create an
      CbcSimpleInteger object for each integer variable.
    */
  findIntegers(false);
  // Say not dynamic pseudo costs
  ownership_ &= ~0x40000000;
  // If dynamic pseudo costs then do
  if (numberBeforeTrust_)
    convertToDynamic();
  // Set up char array to say if integer (speed)
  delete[] integerInfo_;
  {
    int n = solver_->getNumCols();
    integerInfo_ = new char[n];
    for (int i = 0; i < n; i++) {
      if (solver_->isInteger(i))
        integerInfo_[i] = 1;
      else
        integerInfo_[i] = 0;
    }
  }
  if (preferredWay_) {
    // set all unset ones
    for (int iObject = 0; iObject < numberObjects_; iObject++) {
      CbcObject *obj = dynamic_cast< CbcObject * >(object_[iObject]);
      if (obj && !obj->preferredWay())
        obj->setPreferredWay(preferredWay_);
    }
  }
  /*
      Ensure that objects on the lists of OsiObjects, heuristics, and cut
      generators attached to this model all refer to this model.
    */
  synchronizeModel();
  if (!solverCharacteristics_) {
    OsiBabSolver *solverCharacteristics = dynamic_cast< OsiBabSolver * >(solver_->getAuxiliaryInfo());
    if (solverCharacteristics) {
      solverCharacteristics_ = solverCharacteristics;
    } else {
      // replace in solver
      OsiBabSolver defaultC;
      solver_->setAuxiliaryInfo(&defaultC);
      solverCharacteristics_ = dynamic_cast< OsiBabSolver * >(solver_->getAuxiliaryInfo());
    }
  }

  solverCharacteristics_->setSolver(solver_);
  // Set so we can tell we are in initial phase in resolve
  continuousObjective_ = -COIN_DBL_MAX;
  /*
      Solve the relaxation.

      Apparently there are circumstances where this will be non-trivial --- i.e.,
      we've done something since initialSolve that's trashed the solution to the
      continuous relaxation.
    */
  /* Tell solver we are in Branch and Cut
       Could use last parameter for subtle differences */
  solver_->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, NULL);
#ifdef COIN_HAS_CLP
  {
    OsiClpSolverInterface *clpSolver
      = dynamic_cast< OsiClpSolverInterface * >(solver_);
    if (clpSolver) {
      if ( this->keepNamesPreproc == false )
        solver_->setIntParam( OsiNameDiscipline, 0 );
      ClpSimplex *clpSimplex = clpSolver->getModelPtr();
      if ((specialOptions_ & 32) == 0) {
        // take off names (unless going to be saving)
        int nameDisc; solver_->getIntParam( OsiNameDiscipline, nameDisc );
        if ( (numberAnalyzeIterations_ >= 0 || (-numberAnalyzeIterations_ & 64) == 0) && (!nameDisc) )
          clpSimplex->dropNames();
      }
      // no crunch if mostly continuous
      if ((clpSolver->specialOptions() & (1 + 8)) != (1 + 8)) {
        int numberColumns = solver_->getNumCols();
        if (numberColumns > 1000 && numberIntegers_ * 4 < numberColumns)
          clpSolver->setSpecialOptions(clpSolver->specialOptions() & (~1));
      }
      //#define NO_CRUNCH
#ifdef NO_CRUNCH
      printf("TEMP switching off crunch\n");
      int iOpt = clpSolver->specialOptions();
      iOpt &= ~1;
      iOpt |= 65536;
      clpSolver->setSpecialOptions(iOpt);
#endif
    }
  }
#endif
  bool feasible;
  numberSolves_ = 0;
  {
    // check
    int numberOdd = 0;
    for (int i = 0; i < numberObjects_; i++) {
      CbcSimpleInteger *obj = dynamic_cast< CbcSimpleInteger * >(object_[i]);
      if (!obj)
        numberOdd++;
    }
    if (numberOdd) {
      moreSpecialOptions_ |= 1073741824;
      // also switch off checking for restart (as preprocessing may be odd)
      specialOptions_ &= ~(512|32768);
    }
  }
  // If NLP then we assume already solved outside branchAndbound
  if (!solverCharacteristics_->solverType() || solverCharacteristics_->solverType() == 4) {
    feasible = resolve(NULL, 0) != 0;
  } else {
    // pick up given status
    feasible = (solver_->isProvenOptimal() && !solver_->isDualObjectiveLimitReached());
  }
  if (problemFeasibility_->feasible(this, 0) < 0) {
    feasible = false; // pretend infeasible
  }
  numberSavedSolutions_ = 0;
  int saveNumberStrong = numberStrong_;
  int saveNumberBeforeTrust = numberBeforeTrust_;
  /*
      If the linear relaxation of the root is infeasible, bail out now. Otherwise,
      continue with processing the root node.
    */
  if (!feasible) {
    status_ = 0;
    if (!solver_->isProvenDualInfeasible()) {
      handler_->message(CBC_INFEAS, messages_) << CoinMessageEol;
      secondaryStatus_ = 1;
    } else {
      handler_->message(CBC_UNBOUNDED, messages_) << CoinMessageEol;
      secondaryStatus_ = 7;
    }
    originalContinuousObjective_ = COIN_DBL_MAX;
    if (bestSolution_ && ((specialOptions_ & 8388608) == 0 || (specialOptions_ & 2048) != 0)) {
      // best solution found by various heuristics - set solution
      char general[200];
      sprintf(general, "Solution of %g already found by heuristic",
        bestObjective_);
      messageHandler()->message(CBC_GENERAL,
        messages())
        << general << CoinMessageEol;
      setCutoff(1.0e50); // As best solution should be worse than cutoff
      // change cutoff as constraint if wanted
      if (cutoffRowNumber_ >= 0) {
        if (solver_->getNumRows() > cutoffRowNumber_)
          solver_->setRowUpper(cutoffRowNumber_, 1.0e50);
      }
      // also in continuousSolver_
      if (continuousSolver_) {
        // Solvers know about direction
        double direction = solver_->getObjSense();
        continuousSolver_->setDblParam(OsiDualObjectiveLimit, 1.0e50 * direction);
      } else {
        continuousSolver_ = solver_->clone();
      }
      phase_ = 5;
      double increment = getDblParam(CbcModel::CbcCutoffIncrement);
      if ((specialOptions_ & 4) == 0)
        bestObjective_ += 100.0 * increment + 1.0e-3; // only set if we are going to solve
      setBestSolution(CBC_END_SOLUTION, bestObjective_, bestSolution_, 1);
      continuousSolver_->resolve();
      if (!continuousSolver_->isProvenOptimal()) {
        continuousSolver_->messageHandler()->setLogLevel(2);
        continuousSolver_->initialSolve();
      }
      delete solver_;
      solverCharacteristics_ = NULL;
      solver_ = continuousSolver_;
      setPointers(solver_);
      continuousSolver_ = NULL;
    }
    solverCharacteristics_ = NULL;
    if (flipObjective)
      flipModel();
    return;
  } else if (!numberObjects_) {
    // nothing to do
    // Undo preprocessing performed during BaB.
    if (strategy_ && strategy_->preProcessState() > 0) {
      // undo preprocessing
      CglPreProcess *process = strategy_->process();
      assert(process);
      int n = originalSolver->getNumCols();
      if (bestSolution_) {
        delete[] bestSolution_;
        bestSolution_ = new double[n];
        process->postProcess(*solver_);
      }
      strategy_->deletePreProcess();
      // Solution now back in originalSolver
      delete solver_;
      solver_ = originalSolver;
      if (bestSolution_) {
        bestObjective_ = solver_->getObjValue() * solver_->getObjSense();
        memcpy(bestSolution_, solver_->getColSolution(), n * sizeof(double));
      }
      // put back original objects if there were any
      if (originalObject) {
        int iColumn;
        assert(ownObjects_);
        for (iColumn = 0; iColumn < numberObjects_; iColumn++)
          delete object_[iColumn];
        delete[] object_;
        numberObjects_ = numberOriginalObjects;
        object_ = originalObject;
        delete[] integerVariable_;
        numberIntegers_ = 0;
        for (iColumn = 0; iColumn < n; iColumn++) {
          if (solver_->isInteger(iColumn))
            numberIntegers_++;
        }
        integerVariable_ = new int[numberIntegers_];
        numberIntegers_ = 0;
        for (iColumn = 0; iColumn < n; iColumn++) {
          if (solver_->isInteger(iColumn))
            integerVariable_[numberIntegers_++] = iColumn;
        }
      }
    }
    if (flipObjective)
      flipModel();
    solverCharacteristics_ = NULL;
    bestObjective_ = solver_->getObjValue() * solver_->getObjSense();
    int numberColumns = solver_->getNumCols();
    delete[] bestSolution_;
    bestSolution_ = new double[numberColumns];
    CoinCopyN(solver_->getColSolution(), numberColumns, bestSolution_);
    return;
  }
  /*
      See if we're using the Osi side of the branching hierarchy. If so, either
      convert existing CbcObjects to OsiObjects, or generate them fresh. In the
      first case, CbcModel owns the objects on the object_ list. In the second
      case, the solver holds the objects and object_ simply points to the
      solver's list.

      080417 The conversion code here (the block protected by `if (obj)') cannot
      possibly be correct. On the Osi side, descent is OsiObject -> OsiObject2 ->
      all other Osi object classes. On the Cbc side, it's OsiObject -> CbcObject
      -> all other Cbc object classes. It's structurally impossible for any Osi
      object to descend from CbcObject. The only thing I can see is that this is
      really dead code, and object detection is now handled from the Osi side.
    */
  // Convert to Osi if wanted
  //OsiBranchingInformation * persistentInfo = NULL;
  if (branchingMethod_ && branchingMethod_->chooseMethod()) {
    //persistentInfo = new OsiBranchingInformation(solver_);
    if (numberOriginalObjects) {
      for (int iObject = 0; iObject < numberObjects_; iObject++) {
        CbcObject *obj = dynamic_cast< CbcObject * >(object_[iObject]);
        if (obj) {
          CbcSimpleInteger *obj2 = dynamic_cast< CbcSimpleInteger * >(obj);
          if (obj2) {
            // back to Osi land
            object_[iObject] = obj2->osiObject();
            delete obj;
          } else {
            OsiSimpleInteger *obj3 = dynamic_cast< OsiSimpleInteger * >(obj);
            if (!obj3) {
              OsiSOS *obj4 = dynamic_cast< OsiSOS * >(obj);
              if (!obj4) {
                CbcSOS *obj5 = dynamic_cast< CbcSOS * >(obj);
                if (obj5) {
                  // back to Osi land
                  object_[iObject] = obj5->osiObject(solver_);
                } else {
                  printf("Code up CbcObject type in Osi land\n");
                  abort();
                }
              }
            }
          }
        }
      }
      // and add to solver
      //if (!solver_->numberObjects()) {
      solver_->addObjects(numberObjects_, object_);
      //} else {
      //if (solver_->numberObjects()!=numberOriginalObjects) {
      //printf("should have trapped that solver has objects before\n");
      //abort();
      //}
      //}
    } else {
      /*
              As of 080104, findIntegersAndSOS is misleading --- the default OSI
              implementation finds only integers.
            */
      // do from solver
      deleteObjects(false);
      solver_->findIntegersAndSOS(false);
      numberObjects_ = solver_->numberObjects();
      object_ = solver_->objects();
      ownObjects_ = false;
    }
    branchingMethod_->chooseMethod()->setSolver(solver_);
  }
  // take off heuristics if have to (some do not work with SOS, for example)
  // object should know what's safe.
  {
    int numberOdd = 0;
    int numberSOS = 0;
    for (int i = 0; i < numberObjects_; i++) {
      if (!object_[i]->canDoHeuristics())
        numberOdd++;
      CbcSOS *obj = dynamic_cast< CbcSOS * >(object_[i]);
      if (obj)
        numberSOS++;
    }
    if (numberOdd) {
      if (numberHeuristics_ && (specialOptions_ & 1024) == 0) {
        int k = 0;
        for (int i = 0; i < numberHeuristics_; i++) {
          if (!heuristic_[i]->canDealWithOdd())
            delete heuristic_[i];
          else
            heuristic_[k++] = heuristic_[i];
        }
        if (!k) {
          delete[] heuristic_;
          heuristic_ = NULL;
        }
        numberHeuristics_ = k;
        handler_->message(CBC_HEURISTICS_OFF, messages_) << numberOdd << CoinMessageEol;
      }
      // If odd switch off AddIntegers
      specialOptions_ &= ~65536;
      // switch off fast nodes for now
      fastNodeDepth_ = -1;
      moreSpecialOptions_ &= ~33554432; // no diving
    } else if (numberSOS) {
      specialOptions_ |= 128; // say can do SOS in dynamic mode
      // switch off fast nodes for now
      fastNodeDepth_ = -1;
      moreSpecialOptions_ &= ~33554432; // no diving
    }
    if (numberThreads_ > 0) {
      /* switch off fast nodes for now
	       Trouble is that by time mini bab finishes code is
	       looking at a different node
	     */
      fastNodeDepth_ = -1;
    }
  }
  // Save objective (just so user can access it)
  originalContinuousObjective_ = solver_->getObjValue() * solver_->getObjSense();
  bestPossibleObjective_ = originalContinuousObjective_;
  sumChangeObjective1_ = 0.0;
  sumChangeObjective2_ = 0.0;
  /*
      OsiRowCutDebugger knows an optimal answer for a subset of MIP problems.
      Assuming it recognises the problem, when called upon it will check a cut to
      see if it cuts off the optimal answer.
    */
  // If debugger exists set specialOptions_ bit
  if (solver_->getRowCutDebuggerAlways()) {
    specialOptions_ |= 1;
  }

#ifdef CBC_DEBUG
  if ((specialOptions_ & 1) == 0)
    solver_->activateRowCutDebugger(problemName.c_str());
  if (solver_->getRowCutDebuggerAlways())
    specialOptions_ |= 1;
#endif

  /*
      Begin setup to process a feasible root node.
    */
  bestObjective_ = CoinMin(bestObjective_, 1.0e50);
  if (!bestSolution_) {
    numberSolutions_ = 0;
    numberHeuristicSolutions_ = 0;
  }
  stateOfSearch_ = 0;
  // Everything is minimization
  {
    // needed to sync cutoffs
    double value;
    solver_->getDblParam(OsiDualObjectiveLimit, value);
    dblParam_[CbcCurrentCutoff] = value * solver_->getObjSense();
  }
  double cutoff = getCutoff();
  double direction = solver_->getObjSense();
  dblParam_[CbcOptimizationDirection] = direction;
  if (cutoff < 1.0e20 && direction < 0.0)
    messageHandler()->message(CBC_CUTOFF_WARNING1,
      messages())
      << cutoff << -cutoff << CoinMessageEol;
  if (cutoff > bestObjective_)
    cutoff = bestObjective_;
  setCutoff(cutoff);
  /*
      We probably already have a current solution, but just in case ...
    */
  int numberColumns = getNumCols();
  if (!currentSolution_)
    currentSolution_ = new double[numberColumns];
  testSolution_ = currentSolution_;
  /*
      Create a copy of the solver, thus capturing the original (root node)
      constraint system (aka the continuous system).
    */
  delete continuousSolver_;
  continuousSolver_ = solver_->clone();
#ifdef CONFLICT_CUTS
  if ((moreSpecialOptions_ & 4194304) != 0) {
#ifdef COIN_HAS_CLP
    OsiClpSolverInterface *clpSolver
      = dynamic_cast< OsiClpSolverInterface * >(solver_);
    if (clpSolver) {
      int specialOptions = clpSolver->getModelPtr()->specialOptions();
      // 2097152 switches on rays in crunch
      if (!parentModel_)
        clpSolver->getModelPtr()->setSpecialOptions(specialOptions | 32 | 2097152);
      else
        clpSolver->getModelPtr()->setSpecialOptions(specialOptions & ~(32 | 2097152));
    }
  }
#endif
#endif
#ifdef COIN_HAS_NTY
  // maybe allow on fix and restart later
  if ((moreSpecialOptions2_ & (128 | 256)) != 0 && !parentModel_) {
    symmetryInfo_ = new CbcSymmetry();
    symmetryInfo_->setupSymmetry(this);
    int numberGenerators = symmetryInfo_->statsOrbits(this, 0);
    if (!symmetryInfo_->numberUsefulOrbits() && (moreSpecialOptions2_ & (128 | 256)) != (128 | 256)) {
      delete symmetryInfo_;
      symmetryInfo_ = NULL;
      moreSpecialOptions2_ &= ~(128 | 256);
    }
    if ((moreSpecialOptions2_ & (128 | 256)) == (128 | 256)) {
      //moreSpecialOptions2_ &= ~256;
    }
  }
#endif

  // add cutoff as constraint if wanted
  if (cutoffRowNumber_ == -2) {
    if (!parentModel_) {
      int numberColumns = solver_->getNumCols();
      double *obj = CoinCopyOfArray(solver_->getObjCoefficients(), numberColumns);
      int *indices = new int[numberColumns];
      int n = 0;
      for (int i = 0; i < numberColumns; i++) {
        if (obj[i]) {
          indices[n] = i;
          obj[n++] = obj[i];
        }
      }
      if (n) {
        double cutoff = getCutoff();
        // relax a little bit
        cutoff += 1.0e-4;
        double offset;
        solver_->getDblParam(OsiObjOffset, offset);
        cutoffRowNumber_ = solver_->getNumRows();
        solver_->addRow(n, indices, obj, -COIN_DBL_MAX, CoinMin(cutoff, 1.0e25) + offset);
      } else {
        // no objective!
        cutoffRowNumber_ = -1;
      }
      delete[] indices;
      delete[] obj;
    } else {
      // switch off
      cutoffRowNumber_ = -1;
    }
  }
  numberRowsAtContinuous_ = getNumRows();
  solver_->saveBaseModel();
  /*
      Check the objective to see if we can deduce a nontrivial increment. If
      it's better than the current value for CbcCutoffIncrement, it'll be
      installed.
    */
  if (solverCharacteristics_->reducedCostsAccurate())
    analyzeObjective();
  {
    // may be able to change cutoff now
    double cutoff = getCutoff();
    double increment = getDblParam(CbcModel::CbcCutoffIncrement);
    if (cutoff > bestObjective_ - increment) {
      cutoff = bestObjective_ - increment;
      setCutoff(cutoff);
    }
  }
#ifdef COIN_HAS_CLP
  // Possible save of pivot method
  ClpDualRowPivot *savePivotMethod = NULL;
  {
    // pass tolerance and increment to solver
    OsiClpSolverInterface *clpSolver
      = dynamic_cast< OsiClpSolverInterface * >(solver_);
    if (clpSolver)
      clpSolver->setStuff(getIntegerTolerance(), getCutoffIncrement());
#ifdef CLP_RESOLVE
    if ((moreSpecialOptions_ & 1048576) != 0 && !parentModel_ && clpSolver) {
      resolveClp(clpSolver, 0);
    }
#endif
  }
#endif
  /*
      Set up for cut generation. addedCuts_ holds the cuts which are relevant for
      the active subproblem. whichGenerator will be used to record the generator
      that produced a given cut.
    */
#define INITIAL_MAXIMUM_WHICH 1000
  maximumWhich_ = INITIAL_MAXIMUM_WHICH;
  delete[] whichGenerator_;
  whichGenerator_ = new int[maximumWhich_];
  memset(whichGenerator_, 0, maximumWhich_ * sizeof(int));
  maximumNumberCuts_ = 0;
  currentNumberCuts_ = 0;
  delete[] addedCuts_;
  addedCuts_ = NULL;
  OsiObject **saveObjects = NULL;
  maximumRows_ = numberRowsAtContinuous_;
  currentDepth_ = 0;
  workingBasis_.resize(maximumRows_, numberColumns);
  /*
      Set up an empty heap and associated data structures to hold the live set
      (problems which require further exploration).
    */
  CbcCompareDefault *compareActual
    = dynamic_cast< CbcCompareDefault * >(nodeCompare_);
  if (compareActual) {
    compareActual->setBestPossible(direction * solver_->getObjValue());
    compareActual->setCutoff(getCutoff());
#ifdef JJF_ZERO
    if (false && !numberThreads_ && !parentModel_) {
      printf("CbcTreeArray ? threads ? parentArray\n");
      // Setup new style tree
      delete tree_;
      tree_ = new CbcTreeArray();
    }
#endif
  }
  tree_->setComparison(*nodeCompare_);
  /*
      Used to record the path from a node to the root of the search tree, so that
      we can then traverse from the root to the node when restoring a subproblem.
    */
  maximumDepth_ = 10;
  delete[] walkback_;
  walkback_ = new CbcNodeInfo *[maximumDepth_];
  lastDepth_ = 0;
  delete[] lastNodeInfo_;
  lastNodeInfo_ = new CbcNodeInfo *[maximumDepth_];
  delete[] lastNumberCuts_;
  lastNumberCuts_ = new int[maximumDepth_];
  maximumCuts_ = 100;
  lastNumberCuts2_ = 0;
  delete[] lastCut_;
  lastCut_ = new const OsiRowCut *[maximumCuts_];
  /*
      Used to generate bound edits for CbcPartialNodeInfo.
    */
  double *lowerBefore = new double[numberColumns];
  double *upperBefore = new double[numberColumns];
  /*
    Set up to run heuristics and generate cuts at the root node. The heavy
    lifting is hidden inside the calls to doHeuristicsAtRoot and solveWithCuts.

    To start, tell cut generators they can be a bit more aggressive at the
    root node.

    QUESTION: phase_ = 0 is documented as `initial solve', phase = 1 as `solve
        with cuts at root'. Is phase_ = 1 the correct indication when
        doHeurisiticsAtRoot is called to run heuristics outside of the main
        cut / heurisitc / reoptimise loop in solveWithCuts?

      Generate cuts at the root node and reoptimise. solveWithCuts does the heavy
      lifting. It will iterate a generate/reoptimise loop (including reduced cost
      fixing) until no cuts are generated, the change in objective falls off,  or
      the limit on the number of rounds of cut generation is exceeded.

      At the end of all this, any cuts will be recorded in cuts and also
      installed in the solver's constraint system. We'll have reoptimised, and
      removed any slack cuts (numberOldActiveCuts_ and numberNewCuts_ have been
      adjusted accordingly).

      Tell cut generators they can be a bit more aggressive at root node

      TODO: Why don't we make a copy of the solution after solveWithCuts?
      TODO: If numberUnsatisfied == 0, don't we have a solution?
    */
  phase_ = 1;
  int iCutGenerator;
  for (iCutGenerator = 0; iCutGenerator < numberCutGenerators_; iCutGenerator++) {
    // If parallel switch off global cuts
    if (numberThreads_) {
      generator_[iCutGenerator]->setGlobalCuts(false);
      generator_[iCutGenerator]->setGlobalCutsAtRoot(false);
    }
    CglCutGenerator *generator = generator_[iCutGenerator]->generator();
    generator->setAggressiveness(generator->getAggressiveness() + 100);
    if (!generator->canDoGlobalCuts())
      generator->setGlobalCuts(false);
  }
  OsiCuts cuts;
  int anyAction = -1;
  numberOldActiveCuts_ = 0;
  numberNewCuts_ = 0;
  // Array to mark solution
  delete[] usedInSolution_;
  usedInSolution_ = new int[numberColumns];
  CoinZeroN(usedInSolution_, numberColumns);
  /*
      For printing totals and for CbcNode (numberNodes_)
    */
  numberIterations_ = 0;
  numberNodes_ = 0;
  numberNodes2_ = 0;
  maximumStatistics_ = 0;
  maximumDepthActual_ = 0;
  numberDJFixed_ = 0.0;
  if (!parentModel_) {
    if ((specialOptions_ & 262144) != 0) {
      // create empty stored cuts
      //storedRowCuts_ = new CglStored(solver_->getNumCols());
    } else if ((specialOptions_ & 524288) != 0 && storedRowCuts_) {
      // tighten and set best solution
      // A) tight bounds on integer variables
      /*
                storedRowCuts_ are coming in from outside, probably for nonlinear.
              John was unsure about origin.
            */
      const double *lower = solver_->getColLower();
      const double *upper = solver_->getColUpper();
      const double *tightLower = storedRowCuts_->tightLower();
      const double *tightUpper = storedRowCuts_->tightUpper();
      int nTightened = 0;
      for (int i = 0; i < numberIntegers_; i++) {
        int iColumn = integerVariable_[i];
        if (tightLower[iColumn] > lower[iColumn]) {
          nTightened++;
          solver_->setColLower(iColumn, tightLower[iColumn]);
        }
        if (tightUpper[iColumn] < upper[iColumn]) {
          nTightened++;
          solver_->setColUpper(iColumn, tightUpper[iColumn]);
        }
      }
      if (nTightened)
        COIN_DETAIL_PRINT(printf("%d tightened by alternate cuts\n", nTightened));
      if (storedRowCuts_->bestObjective() < bestObjective_) {
        // B) best solution
        double objValue = storedRowCuts_->bestObjective();
        setBestSolution(CBC_SOLUTION, objValue,
          storedRowCuts_->bestSolution());
        // Do heuristics
        // Allow RINS
        for (int i = 0; i < numberHeuristics_; i++) {
          CbcHeuristicRINS *rins
            = dynamic_cast< CbcHeuristicRINS * >(heuristic_[i]);
          if (rins) {
            rins->setLastNode(-100);
          }
        }
      }
    }
  }
#ifdef SWITCH_VARIABLES
  // see if any switching variables
  if (numberIntegers_ < solver_->getNumCols())
    findSwitching();
#endif
  /*
      Run heuristics at the root. This is the only opportunity to run FPump; it
      will be removed from the heuristics list by doHeuristicsAtRoot.
    */
  // See if multiple runs wanted
  CbcModel **rootModels = NULL;
  if (!parentModel_ && multipleRootTries_ % 100) {
    double rootTimeCpu = CoinCpuTime();
    double startTimeRoot = CoinGetTimeOfDay();
    int numberRootThreads = 1;
    /* undocumented fine tuning
	 aabbcc where cc is number of tries
	 bb if nonzero is number of threads
	 aa if nonzero just do heuristics
      */
    int numberModels = multipleRootTries_ % 100;
#ifdef CBC_THREAD
    numberRootThreads = (multipleRootTries_ / 100) % 100;
    if (!numberRootThreads) {
      if (numberThreads_ < 2)
        numberRootThreads = numberModels;
      else
        numberRootThreads = CoinMin(numberThreads_, numberModels);
    }
#endif
    int otherOptions = (multipleRootTries_ / 10000) % 100;
    rootModels = new CbcModel *[numberModels];
    int newSeed = randomSeed_;
    if (newSeed == 0) {
      double time = fabs(CoinGetTimeOfDay());
      while (time >= COIN_INT_MAX)
        time *= 0.5;
      newSeed = static_cast< unsigned int >(time);
    } else if (newSeed < 0) {
      newSeed = 123456789;
#ifdef COIN_HAS_CLP
      OsiClpSolverInterface *clpSolver
        = dynamic_cast< OsiClpSolverInterface * >(solver_);
      if (clpSolver) {
        newSeed += clpSolver->getModelPtr()->randomNumberGenerator()->getSeed();
      }
#endif
    }
    CoinWarmStartBasis *basis = dynamic_cast< CoinWarmStartBasis * >(solver_->getEmptyWarmStart());
    for (int i = 0; i < numberModels; i++) {
      rootModels[i] = new CbcModel(*this);
      rootModels[i]->setNumberThreads(0);
      rootModels[i]->setMaximumNodes(otherOptions ? -1 : 0);
      rootModels[i]->setRandomSeed(newSeed + 10000000 * i);
      rootModels[i]->randomNumberGenerator()->setSeed(newSeed + 50000000 * i);
      rootModels[i]->setMultipleRootTries(0);
#ifdef COIN_HAS_NTY
      rootModels[i]->zapSymmetry();
      rootModels[i]->moreSpecialOptions2_ &= ~(128 | 256); // off nauty
#endif
      // use seed
      rootModels[i]->setSpecialOptions(specialOptions_ | (4194304 | 8388608));
      rootModels[i]->setMoreSpecialOptions(moreSpecialOptions_ & (~(134217728 | 4194304)));
      rootModels[i]->setMoreSpecialOptions2(moreSpecialOptions2_ & (~(128 | 256)));
      rootModels[i]->solver_->setWarmStart(basis);
#ifdef COIN_HAS_CLP
      OsiClpSolverInterface *clpSolver
        = dynamic_cast< OsiClpSolverInterface * >(rootModels[i]->solver_);
#define NEW_RANDOM_BASIS
#ifdef NEW_RANDOM_BASIS
      if (i == 0)
        continue;
#endif
      if (clpSolver) {
        ClpSimplex *simplex = clpSolver->getModelPtr();
        if (defaultHandler_)
          simplex->setDefaultMessageHandler();
        simplex->setRandomSeed(newSeed + 20000000 * i);
        simplex->allSlackBasis();
        int logLevel = simplex->logLevel();
        if (logLevel == 1)
          simplex->setLogLevel(0);
        if (i != 0) {
#ifdef NEW_RANDOM_BASIS
          int numberRows = simplex->numberRows();
          int throwOut = 20; //2+numberRows/100;
          for (int iThrow = 0; iThrow < throwOut; iThrow++) {
            double random = simplex->randomNumberGenerator()->randomDouble();
            int iStart = static_cast< int >(random * numberRows);
            for (int j = iStart; j < numberRows; j++) {
              if (simplex->getRowStatus(j) != ClpSimplex::basic) {
                simplex->setRowStatus(j, ClpSimplex::basic);
                break;
              }
            }
          }
          clpSolver->setWarmStart(NULL);
#else
            double random = simplex->randomNumberGenerator()->randomDouble();
            int bias = static_cast< int >(random * (numberIterations / 4));
            simplex->setMaximumIterations(numberIterations / 2 + bias);
            simplex->primal();
            simplex->setMaximumIterations(COIN_INT_MAX);
            simplex->dual();
#endif
        } else {
#ifndef NEW_RANDOM_BASIS
          simplex->primal();
#endif
#endif
        }
#ifdef NEW_RANDOM_BASIS
        simplex->setLogLevel(logLevel);
        clpSolver->setWarmStart(NULL);
#endif
      }
      for (int j = 0; j < numberHeuristics_; j++)
        rootModels[i]->heuristic_[j]->setSeed(rootModels[i]->heuristic_[j]->getSeed() + 100000000 * i);
      for (int j = 0; j < numberCutGenerators_; j++)
        rootModels[i]->generator_[j]->generator()->refreshSolver(rootModels[i]->solver_);
    }
    delete basis;
#ifdef CBC_THREAD
    if (numberRootThreads == 1) {
#endif
      for (int iModel = 0; iModel < numberModels; iModel++) {
        doRootCbcThread(rootModels[iModel]);
        // see if solved at root node
        if (rootModels[iModel]->getMaximumNodes()) {
          feasible = false;
          break;
        }
      }
#ifdef CBC_THREAD
    } else {
      Coin_pthread_t *threadId = new Coin_pthread_t[numberRootThreads];
      for (int kModel = 0; kModel < numberModels; kModel += numberRootThreads) {
        bool finished = false;
        for (int iModel = kModel; iModel < CoinMin(numberModels, kModel + numberRootThreads); iModel++) {
          pthread_create(&(threadId[iModel - kModel].thr), NULL,
            doRootCbcThread,
            rootModels[iModel]);
        }
        // wait
        for (int iModel = kModel; iModel < CoinMin(numberModels, kModel + numberRootThreads); iModel++) {
          pthread_join(threadId[iModel - kModel].thr, NULL);
        }
        // see if solved at root node
        for (int iModel = kModel; iModel < CoinMin(numberModels, kModel + numberRootThreads); iModel++) {
          if (rootModels[iModel]->getMaximumNodes())
            finished = true;
        }
        if (finished) {
          feasible = false;
          break;
        }
      }
      delete[] threadId;
    }
#endif
    // sort solutions
    int *which = new int[numberModels];
    double *value = new double[numberModels];
    int numberSolutions = 0;
    for (int iModel = 0; iModel < numberModels; iModel++) {
      if (rootModels[iModel]->bestSolution()) {
        which[numberSolutions] = iModel;
        value[numberSolutions++] = -rootModels[iModel]->getMinimizationObjValue();
      }
    }
    char general[100];
    rootTimeCpu = CoinCpuTime() - rootTimeCpu;
    if (numberRootThreads == 1)
      sprintf(general, "Multiple root solvers took a total of %.2f seconds\n",
        rootTimeCpu);
    else
      sprintf(general, "Multiple root solvers took a total of %.2f seconds (%.2f elapsed)\n",
        rootTimeCpu, CoinGetTimeOfDay() - startTimeRoot);
    messageHandler()->message(CBC_GENERAL,
      messages())
      << general << CoinMessageEol;
    CoinSort_2(value, value + numberSolutions, which);
    // to get name
    CbcHeuristicRINS dummyHeuristic;
    dummyHeuristic.setHeuristicName("Multiple root solvers");
    lastHeuristic_ = &dummyHeuristic;
    for (int i = 0; i < numberSolutions; i++) {
      double objValue = -value[i];
      if (objValue < getCutoff()) {
        int iModel = which[i];
        setBestSolution(CBC_ROUNDING, objValue,
          rootModels[iModel]->bestSolution());
      }
    }
    lastHeuristic_ = NULL;
    delete[] which;
    delete[] value;
  }
  // Do heuristics
  if (numberObjects_ && !rootModels)
    doHeuristicsAtRoot();
  if (solverCharacteristics_->solutionAddsCuts()) {
    // With some heuristics solver needs a resolve here
    solver_->resolve();
    if (!isProvenOptimal()) {
      solver_->initialSolve();
    }
  }
  /*
      Grepping through the code, it would appear that this is a command line
      debugging hook.  There's no obvious place in the code where this is set to
      a negative value.

      User hook, says John.
    */
  if (intParam_[CbcMaxNumNode] < 0
    || numberSolutions_ >= getMaximumSolutions())
    eventHappened_ = true; // stop as fast as possible
  stoppedOnGap_ = false;
  // See if can stop on gap
  bestPossibleObjective_ = solver_->getObjValue() * solver_->getObjSense();
  if (canStopOnGap()) {
    if (bestPossibleObjective_ < getCutoff())
      stoppedOnGap_ = true;
    feasible = false;
    //eventHappened_=true; // stop as fast as possible
  }
  /*
      Set up for statistics collection, if requested. Standard values are
      documented in CbcModel.hpp. The magic number 100 will trigger a dump of
      CbcSimpleIntegerDynamicPseudoCost objects (no others). Looks like another
      command line debugging hook.
    */
  statistics_ = NULL;
  // Do on switch
  if (doStatistics > 0 && doStatistics <= 100) {
    maximumStatistics_ = 10000;
    statistics_ = new CbcStatistics *[maximumStatistics_];
    memset(statistics_, 0, maximumStatistics_ * sizeof(CbcStatistics *));
  }
  // See if we can add integers
  if (noObjects && numberIntegers_ < solver_->getNumCols() && (specialOptions_ & 65536) != 0 && !parentModel_ && false) {
    int numberIntegers1 = 0;
    int numberColumns = solver_->getNumCols();
    for (int i = 0; i < numberColumns; i++) {
      if (solver_->isInteger(i))
        numberIntegers1++;
    }
    AddIntegers();
    // make sure in sync
    int numberIntegers2 = 0;
    for (int i = 0; i < numberColumns; i++) {
      if (solver_->isInteger(i))
        numberIntegers2++;
    }
    if (numberIntegers1 < numberIntegers2) {
      findIntegers(true, 2);
      convertToDynamic();
    }
  }

  /*
      Do an initial round of cut generation for the root node. Depending on the
      type of underlying solver, we may want to do this even if the initial query
      to the objects indicates they're satisfied.

      solveWithCuts does the heavy lifting. It will iterate a generate/reoptimise
      loop (including reduced cost fixing) until no cuts are generated, the
      change in objective falls off,  or the limit on the number of rounds of cut
      generation is exceeded.

      At the end of all this, any cuts will be recorded in cuts and also
      installed in the solver's constraint system. We'll have reoptimised, and
      removed any slack cuts (numberOldActiveCuts_ and numberNewCuts_ have been
      adjusted accordingly).
    */
  int iObject;
  int numberUnsatisfied = 0;
  delete[] currentSolution_;
  currentSolution_ = new double[numberColumns];
  testSolution_ = currentSolution_;
  memcpy(currentSolution_, solver_->getColSolution(),
    numberColumns * sizeof(double));
  // point to useful information
  OsiBranchingInformation usefulInfo = usefulInformation();

  for (iObject = 0; iObject < numberObjects_; iObject++) {
    double infeasibility = object_[iObject]->checkInfeasibility(&usefulInfo);
    if (infeasibility)
      numberUnsatisfied++;
  }
  // replace solverType
  double *tightBounds = NULL;
  if (solverCharacteristics_->tryCuts()) {

    if (numberUnsatisfied) {
      // User event
      if (!eventHappened_ && feasible) {
        if (rootModels) {
          // for fixings
          int numberColumns = solver_->getNumCols();
          tightBounds = new double[2 * numberColumns];
          {
            const double *lower = solver_->getColLower();
            const double *upper = solver_->getColUpper();
            for (int i = 0; i < numberColumns; i++) {
              tightBounds[2 * i + 0] = lower[i];
              tightBounds[2 * i + 1] = upper[i];
            }
          }
          int numberModels = multipleRootTries_ % 100;
          const OsiSolverInterface **solvers = new const OsiSolverInterface *[numberModels];
          int numberRows = continuousSolver_->getNumRows();
          int maxCuts = 0;
          for (int i = 0; i < numberModels; i++) {
            solvers[i] = rootModels[i]->solver();
            const double *lower = solvers[i]->getColLower();
            const double *upper = solvers[i]->getColUpper();
            for (int j = 0; j < numberColumns; j++) {
              tightBounds[2 * j + 0] = CoinMax(lower[j], tightBounds[2 * j + 0]);
              tightBounds[2 * j + 1] = CoinMin(upper[j], tightBounds[2 * j + 1]);
            }
            int numberRows2 = solvers[i]->getNumRows();
            assert(numberRows2 >= numberRows);
            maxCuts += numberRows2 - numberRows;
            // accumulate statistics
            for (int j = 0; j < numberCutGenerators_; j++) {
              generator_[j]->addStatistics(rootModels[i]->cutGenerator(j));
            }
          }
          for (int j = 0; j < numberCutGenerators_; j++) {
            generator_[j]->scaleBackStatistics(numberModels);
          }
          //CbcRowCuts rowCut(maxCuts);
          const OsiRowCutDebugger *debugger = NULL;
          if ((specialOptions_ & 1) != 0)
            debugger = solver_->getRowCutDebugger();
          for (int iModel = 0; iModel < numberModels; iModel++) {
            int numberRows2 = solvers[iModel]->getNumRows();
            const CoinPackedMatrix *rowCopy = solvers[iModel]->getMatrixByRow();
            const int *rowLength = rowCopy->getVectorLengths();
            const double *elements = rowCopy->getElements();
            const int *column = rowCopy->getIndices();
            const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
            const double *rowLower = solvers[iModel]->getRowLower();
            const double *rowUpper = solvers[iModel]->getRowUpper();
            for (int iRow = numberRows; iRow < numberRows2; iRow++) {
              OsiRowCut rc;
              rc.setLb(rowLower[iRow]);
              rc.setUb(rowUpper[iRow]);
              CoinBigIndex start = rowStart[iRow];
              rc.setRow(rowLength[iRow], column + start, elements + start, false);
              if (debugger)
                CoinAssert(!debugger->invalidCut(rc));
              globalCuts_.addCutIfNotDuplicate(rc);
            }
            //int cutsAdded=globalCuts_.numberCuts()-numberCuts;
            //numberCuts += cutsAdded;
            //printf("Model %d gave %d cuts (out of %d possible)\n",
            //	   iModel,cutsAdded,numberRows2-numberRows);
          }
          // normally replace global cuts
          //if (!globalCuts_.())
          //globalCuts_=rowCutrowCut.addCuts(globalCuts_);
          //rowCut.addCuts(globalCuts_);
          int nTightened = 0;
          assert(feasible);
          {
            double tolerance = 1.0e-5;
            const double *lower = solver_->getColLower();
            const double *upper = solver_->getColUpper();
            for (int i = 0; i < numberColumns; i++) {
              if (tightBounds[2 * i + 0] > tightBounds[2 * i + 1] + 1.0e-9) {
                feasible = false;
                char general[200];
                sprintf(general, "Solvers give infeasible bounds on %d %g,%g was %g,%g - search finished\n",
                  i, tightBounds[2 * i + 0], tightBounds[2 * i + 1], lower[i], upper[i]);
                messageHandler()->message(CBC_GENERAL, messages())
                  << general << CoinMessageEol;
                break;
              }
              double oldLower = lower[i];
              double oldUpper = upper[i];
              if (tightBounds[2 * i + 0] > oldLower + tolerance) {
                nTightened++;
                solver_->setColLower(i, tightBounds[2 * i + 0]);
              }
              if (tightBounds[2 * i + 1] < oldUpper - tolerance) {
                nTightened++;
                solver_->setColUpper(i, tightBounds[2 * i + 1]);
              }
            }
          }
          delete[] tightBounds;
          tightBounds = NULL;
          char printBuffer[200];
          sprintf(printBuffer, "%d solvers added %d different cuts out of pool of %d",
            numberModels, globalCuts_.sizeRowCuts(), maxCuts);
          messageHandler()->message(CBC_GENERAL, messages())
            << printBuffer << CoinMessageEol;
          if (nTightened) {
            sprintf(printBuffer, "%d bounds were tightened",
              nTightened);
            messageHandler()->message(CBC_GENERAL, messages())
              << printBuffer << CoinMessageEol;
          }
          delete[] solvers;
        }
        if (!parentModel_ && (moreSpecialOptions_ & 67108864) != 0) {
          // load cuts from file
          FILE *fp = fopen("global.cuts", "rb");
          if (fp) {
            size_t nRead;
            int numberColumns = solver_->getNumCols();
            int numCols;
            nRead = fread(&numCols, sizeof(int), 1, fp);
            if (nRead != 1)
              throw("Error in fread");
            if (numberColumns != numCols) {
              printf("Mismatch on columns %d %d\n", numberColumns, numCols);
              fclose(fp);
            } else {
              // If rootModel just do some
              double threshold = -1.0;
              if (!multipleRootTries_)
                threshold = 0.5;
              int initialCuts = 0;
              int initialGlobal = globalCuts_.sizeRowCuts();
              double *elements = new double[numberColumns + 2];
              int *indices = new int[numberColumns];
              int numberEntries = 1;
              while (numberEntries > 0) {
                nRead = fread(&numberEntries, sizeof(int), 1, fp);
                if (nRead != 1)
                  throw("Error in fread");
                double randomNumber = randomNumberGenerator_.randomDouble();
                if (numberEntries > 0) {
                  initialCuts++;
                  nRead = fread(elements, sizeof(double), numberEntries + 2, fp);
                  if (nRead != static_cast< size_t >(numberEntries + 2))
                    throw("Error in fread");
                  nRead = fread(indices, sizeof(int), numberEntries, fp);
                  if (nRead != static_cast< size_t >(numberEntries))
                    throw("Error in fread");
                  if (randomNumber > threshold) {
                    OsiRowCut rc;
                    rc.setLb(elements[numberEntries]);
                    rc.setUb(elements[numberEntries + 1]);
                    rc.setRow(numberEntries, indices, elements,
                      false);
                    rc.setGloballyValidAsInteger(2);
                    globalCuts_.addCutIfNotDuplicate(rc);
                  }
                }
              }
              fclose(fp);
              // fixes
              int nTightened = 0;
              fp = fopen("global.fix", "rb");
              if (fp) {
                nRead = fread(indices, sizeof(int), 2, fp);
                if (nRead != 2)
                  throw("Error in fread");
                if (numberColumns != indices[0]) {
                  printf("Mismatch on columns %d %d\n", numberColumns,
                    indices[0]);
                } else {
                  indices[0] = 1;
                  while (indices[0] >= 0) {
                    nRead = fread(indices, sizeof(int), 2, fp);
                    if (nRead != 2)
                      throw("Error in fread");
                    int iColumn = indices[0];
                    if (iColumn >= 0) {
                      nTightened++;
                      nRead = fread(elements, sizeof(double), 4, fp);
                      if (nRead != 4)
                        throw("Error in fread");
                      solver_->setColLower(iColumn, elements[0]);
                      solver_->setColUpper(iColumn, elements[1]);
                    }
                  }
                }
              }
              if (fp)
                fclose(fp);
              char printBuffer[200];
              sprintf(printBuffer, "%d cuts read in of which %d were unique, %d bounds tightened",
                initialCuts,
                globalCuts_.sizeRowCuts() - initialGlobal, nTightened);
              messageHandler()->message(CBC_GENERAL, messages())
                << printBuffer << CoinMessageEol;
              delete[] elements;
              delete[] indices;
            }
          }
        }
        if (feasible)
          feasible = solveWithCuts(cuts, maximumCutPassesAtRoot_,
            NULL);
        if (multipleRootTries_ && (moreSpecialOptions_ & 134217728) != 0) {
          FILE *fp = NULL;
          size_t nRead;
          int numberColumns = solver_->getNumCols();
          int initialCuts = 0;
          if ((moreSpecialOptions_ & 134217728) != 0) {
            // append so go down to end
            fp = fopen("global.cuts", "r+b");
            if (fp) {
              int numCols;
              nRead = fread(&numCols, sizeof(int), 1, fp);
              if (nRead != 1)
                throw("Error in fread");
              if (numberColumns != numCols) {
                printf("Mismatch on columns %d %d\n", numberColumns, numCols);
                fclose(fp);
                fp = NULL;
              }
            }
          }
          double *elements = new double[numberColumns + 2];
          int *indices = new int[numberColumns];
          if (fp) {
            int numberEntries = 1;
            while (numberEntries > 0) {
              fpos_t position;
              fgetpos(fp, &position);
              nRead = fread(&numberEntries, sizeof(int), 1, fp);
              if (nRead != 1)
                throw("Error in fread");
              if (numberEntries > 0) {
                initialCuts++;
                nRead = fread(elements, sizeof(double), numberEntries + 2, fp);
                if (nRead != static_cast< size_t >(numberEntries + 2))
                  throw("Error in fread");
                nRead = fread(indices, sizeof(int), numberEntries, fp);
                if (nRead != static_cast< size_t >(numberEntries))
                  throw("Error in fread");
              } else {
                // end
                fsetpos(fp, &position);
              }
            }
          } else {
            fp = fopen("global.cuts", "wb");
            size_t nWrite;
            nWrite = fwrite(&numberColumns, sizeof(int), 1, fp);
            if (nWrite != 1)
              throw("Error in fwrite");
          }
          size_t nWrite;
          // now append binding cuts
          int numberC = continuousSolver_->getNumRows();
          int numberRows = solver_->getNumRows();
          printf("Saving %d cuts (up from %d)\n",
            initialCuts + numberRows - numberC, initialCuts);
          const double *rowLower = solver_->getRowLower();
          const double *rowUpper = solver_->getRowUpper();
          // Row copy
          CoinPackedMatrix matrixByRow(*solver_->getMatrixByRow());
          const double *elementByRow = matrixByRow.getElements();
          const int *column = matrixByRow.getIndices();
          const CoinBigIndex *rowStart = matrixByRow.getVectorStarts();
          const int *rowLength = matrixByRow.getVectorLengths();
          for (int iRow = numberC; iRow < numberRows; iRow++) {
            int n = rowLength[iRow];
            assert(n);
            CoinBigIndex start = rowStart[iRow];
            memcpy(elements, elementByRow + start, n * sizeof(double));
            memcpy(indices, column + start, n * sizeof(int));
            elements[n] = rowLower[iRow];
            elements[n + 1] = rowUpper[iRow];
            nWrite = fwrite(&n, sizeof(int), 1, fp);
            if (nWrite != 1)
              throw("Error in fwrite");
            nWrite = fwrite(elements, sizeof(double), n + 2, fp);
            if (nWrite != static_cast< size_t >(n + 2))
              throw("Error in fwrite");
            nWrite = fwrite(indices, sizeof(int), n, fp);
            if (nWrite != static_cast< size_t >(n))
              throw("Error in fwrite");
          }
          // eof marker
          int eofMarker = -1;
          nWrite = fwrite(&eofMarker, sizeof(int), 1, fp);
          if (nWrite != 1)
            throw("Error in fwrite");
          fclose(fp);
          // do tighter bounds (? later extra to original columns)
          int nTightened = 0;
          const double *lower = solver_->getColLower();
          const double *upper = solver_->getColUpper();
          const double *originalLower = continuousSolver_->getColLower();
          const double *originalUpper = continuousSolver_->getColUpper();
          double tolerance = 1.0e-5;
          for (int i = 0; i < numberColumns; i++) {
            if (lower[i] > originalLower[i] + tolerance) {
              nTightened++;
            }
            if (upper[i] < originalUpper[i] - tolerance) {
              nTightened++;
            }
          }
          if (nTightened) {
            fp = fopen("global.fix", "wb");
            size_t nWrite;
            indices[0] = numberColumns;
            if (originalColumns_)
              indices[1] = COIN_INT_MAX;
            else
              indices[1] = -1;
            nWrite = fwrite(indices, sizeof(int), 2, fp);
            if (nWrite != 2)
              throw("Error in fwrite");
            for (int i = 0; i < numberColumns; i++) {
              int nTightened = 0;
              if (lower[i] > originalLower[i] + tolerance) {
                nTightened++;
              }
              if (upper[i] < originalUpper[i] - tolerance) {
                nTightened++;
              }
              if (nTightened) {
                indices[0] = i;
                if (originalColumns_)
                  indices[1] = originalColumns_[i];
                elements[0] = lower[i];
                elements[1] = upper[i];
                elements[2] = originalLower[i];
                elements[3] = originalUpper[i];
                nWrite = fwrite(indices, sizeof(int), 2, fp);
                if (nWrite != 2)
                  throw("Error in fwrite");
                nWrite = fwrite(elements, sizeof(double), 4, fp);
                if (nWrite != 4)
                  throw("Error in fwrite");
              }
            }
            // eof marker
            indices[0] = -1;
            nWrite = fwrite(indices, sizeof(int), 2, fp);
            if (nWrite != 2)
              throw("Error in fwrite");
            fclose(fp);
          }
          delete[] elements;
          delete[] indices;
        }
        if ((specialOptions_ & 524288) != 0 && !parentModel_
          && storedRowCuts_) {
          if (feasible) {
            /* pick up stuff and try again
                        add cuts, maybe keep around
                        do best solution and if so new heuristics
                        obviously tighten bounds
                        */
            // A and B probably done on entry
            // A) tight bounds on integer variables
            const double *lower = solver_->getColLower();
            const double *upper = solver_->getColUpper();
            const double *tightLower = storedRowCuts_->tightLower();
            const double *tightUpper = storedRowCuts_->tightUpper();
            int nTightened = 0;
            for (int i = 0; i < numberIntegers_; i++) {
              int iColumn = integerVariable_[i];
              if (tightLower[iColumn] > lower[iColumn]) {
                nTightened++;
                solver_->setColLower(iColumn, tightLower[iColumn]);
              }
              if (tightUpper[iColumn] < upper[iColumn]) {
                nTightened++;
                solver_->setColUpper(iColumn, tightUpper[iColumn]);
              }
            }
            if (nTightened)
              COIN_DETAIL_PRINT(printf("%d tightened by alternate cuts\n", nTightened));
            if (storedRowCuts_->bestObjective() < bestObjective_) {
              // B) best solution
              double objValue = storedRowCuts_->bestObjective();
              setBestSolution(CBC_SOLUTION, objValue,
                storedRowCuts_->bestSolution());
              // Do heuristics
              // Allow RINS
              for (int i = 0; i < numberHeuristics_; i++) {
                CbcHeuristicRINS *rins
                  = dynamic_cast< CbcHeuristicRINS * >(heuristic_[i]);
                if (rins) {
                  rins->setLastNode(-100);
                }
              }
              doHeuristicsAtRoot();
            }
#ifdef JJF_ZERO
            int nCuts = storedRowCuts_->sizeRowCuts();
            // add to global list
            for (int i = 0; i < nCuts; i++) {
              OsiRowCut newCut(*storedRowCuts_->rowCutPointer(i));
              newCut.setGloballyValidAsInteger(2);
              newCut.mutableRow().setTestForDuplicateIndex(false);
              globalCuts_.insert(newCut);
            }
#else
          addCutGenerator(storedRowCuts_, -99, "Stored from previous run",
            true, false, false, -200);
#endif
            // Set cuts as active
            delete[] addedCuts_;
            maximumNumberCuts_ = cuts.sizeRowCuts();
            if (maximumNumberCuts_) {
              addedCuts_ = new CbcCountRowCut *[maximumNumberCuts_];
            } else {
              addedCuts_ = NULL;
            }
            for (int i = 0; i < maximumNumberCuts_; i++)
              addedCuts_[i] = new CbcCountRowCut(*cuts.rowCutPtr(i),
                NULL, -1, -1, 2);
            COIN_DETAIL_PRINT(printf("size %d\n", cuts.sizeRowCuts()));
            cuts = OsiCuts();
            currentNumberCuts_ = maximumNumberCuts_;
            feasible = solveWithCuts(cuts, maximumCutPassesAtRoot_,
              NULL);
            for (int i = 0; i < maximumNumberCuts_; i++)
              delete addedCuts_[i];
          }
          delete storedRowCuts_;
          storedRowCuts_ = NULL;
        }
      } else {
        feasible = false;
      }
    } else if (solverCharacteristics_->solutionAddsCuts() || solverCharacteristics_->alwaysTryCutsAtRootNode()) {
      // may generate cuts and turn the solution
      //to an infeasible one
      feasible = solveWithCuts(cuts, 2,
        NULL);
    }
  }
  if (rootModels) {
    int numberModels = multipleRootTries_ % 100;
    for (int i = 0; i < numberModels; i++)
      delete rootModels[i];
    delete[] rootModels;
  }
  // check extra info on feasibility
  if (!solverCharacteristics_->mipFeasible())
    feasible = false;
  // If max nodes==0 - don't do strong branching
  if (!getMaximumNodes()) {
    if (feasible)
      feasible = false;
    else
      setMaximumNodes(1); //allow to stop on success
  }
  topOfTree_ = NULL;
#ifdef CLP_RESOLVE
  if ((moreSpecialOptions_ & 2097152) != 0 && !parentModel_ && feasible) {
    OsiClpSolverInterface *clpSolver
      = dynamic_cast< OsiClpSolverInterface * >(solver_);
    if (clpSolver)
      resolveClp(clpSolver, 0);
  }
#endif
  // make cut generators less aggressive
  for (iCutGenerator = 0; iCutGenerator < numberCutGenerators_; iCutGenerator++) {
    CglCutGenerator *generator = generator_[iCutGenerator]->generator();
    generator->setAggressiveness(generator->getAggressiveness() - 100);
  }
  currentNumberCuts_ = numberNewCuts_;
  if (solverCharacteristics_->solutionAddsCuts()) {
    // With some heuristics solver needs a resolve here (don't know if this is bug in heuristics)
    solver_->resolve();
    if (!isProvenOptimal()) {
      solver_->initialSolve();
    }
  }
  // See if can stop on gap
  bestPossibleObjective_ = solver_->getObjValue() * solver_->getObjSense();
  if (canStopOnGap()) {
    if (bestPossibleObjective_ < getCutoff())
      stoppedOnGap_ = true;
    feasible = false;
  }
  // User event
  if (eventHappened_)
    feasible = false;
#if defined(COIN_HAS_CLP) && defined(COIN_HAS_CPX)
  /*
      This is the notion of using Cbc stuff to get going, then calling cplex to
      finish off.
    */
  if (feasible && (specialOptions_ & 16384) != 0 && fastNodeDepth_ == -2 && !parentModel_) {
    // Use Cplex to do search!
    double time1 = CoinCpuTime();
    OsiClpSolverInterface *clpSolver
      = dynamic_cast< OsiClpSolverInterface * >(solver_);
    OsiCpxSolverInterface cpxSolver;
    double direction = clpSolver->getObjSense();
    cpxSolver.setObjSense(direction);
    // load up cplex
    const CoinPackedMatrix *matrix = continuousSolver_->getMatrixByCol();
    const double *rowLower = continuousSolver_->getRowLower();
    const double *rowUpper = continuousSolver_->getRowUpper();
    const double *columnLower = continuousSolver_->getColLower();
    const double *columnUpper = continuousSolver_->getColUpper();
    const double *objective = continuousSolver_->getObjCoefficients();
    cpxSolver.loadProblem(*matrix, columnLower, columnUpper,
      objective, rowLower, rowUpper);
    double *setSol = new double[numberIntegers_];
    int *setVar = new int[numberIntegers_];
    // cplex doesn't know about objective offset
    double offset = clpSolver->getModelPtr()->objectiveOffset();
    for (int i = 0; i < numberIntegers_; i++) {
      int iColumn = integerVariable_[i];
      cpxSolver.setInteger(iColumn);
      if (bestSolution_) {
        setSol[i] = bestSolution_[iColumn];
        setVar[i] = iColumn;
      }
    }
    CPXENVptr env = cpxSolver.getEnvironmentPtr();
    CPXLPptr lpPtr = cpxSolver.getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL);
    cpxSolver.switchToMIP();
    if (bestSolution_) {
#if 0
            CPXcopymipstart(env, lpPtr, numberIntegers_, setVar, setSol);
#else
        int zero = 0;
        CPXaddmipstarts(env, lpPtr, 1, numberIntegers_, &zero, setVar, setSol, NULL, NULL);
#endif
    }
    if (clpSolver->getNumRows() > continuousSolver_->getNumRows() && false) {
      // add cuts
      const CoinPackedMatrix *matrix = clpSolver->getMatrixByRow();
      const double *rhs = clpSolver->getRightHandSide();
      const char *rowSense = clpSolver->getRowSense();
      const double *elementByRow = matrix->getElements();
      const int *column = matrix->getIndices();
      const CoinBigIndex *rowStart = matrix->getVectorStarts();
      const int *rowLength = matrix->getVectorLengths();
      int nStart = continuousSolver_->getNumRows();
      int nRows = clpSolver->getNumRows();
      int size = rowStart[nRows - 1] + rowLength[nRows - 1] - rowStart[nStart];
      int nAdd = 0;
      double *rmatval = new double[size];
      int *rmatind = new int[size];
      int *rmatbeg = new int[nRows - nStart + 1];
      size = 0;
      rmatbeg[0] = 0;
      for (int i = nStart; i < nRows; i++) {
        for (int k = rowStart[i]; k < rowStart[i] + rowLength[i]; k++) {
          rmatind[size] = column[k];
          rmatval[size++] = elementByRow[k];
        }
        nAdd++;
        rmatbeg[nAdd] = size;
      }
      CPXaddlazyconstraints(env, lpPtr, nAdd, size,
        rhs, rowSense, rmatbeg,
        rmatind, rmatval, NULL);
      CPXsetintparam(env, CPX_PARAM_REDUCE,
        // CPX_PREREDUCE_NOPRIMALORDUAL (0)
        CPX_PREREDUCE_PRIMALONLY);
    }
    if (getCutoff() < 1.0e50) {
      double useCutoff = getCutoff() + offset;
      if (bestObjective_ < 1.0e50)
        useCutoff = bestObjective_ + offset + 1.0e-7;
      cpxSolver.setDblParam(OsiDualObjectiveLimit, useCutoff * direction);
      if (direction > 0.0)
        CPXsetdblparam(env, CPX_PARAM_CUTUP, useCutoff); // min
      else
        CPXsetdblparam(env, CPX_PARAM_CUTLO, useCutoff); // max
    }
    CPXsetdblparam(env, CPX_PARAM_EPGAP, dblParam_[CbcAllowableFractionGap]);
    delete[] setSol;
    delete[] setVar;
    char printBuffer[200];
    if (offset) {
      sprintf(printBuffer, "Add %g to all Cplex messages for true objective",
        -offset);
      messageHandler()->message(CBC_GENERAL, messages())
        << printBuffer << CoinMessageEol;
      cpxSolver.setDblParam(OsiObjOffset, offset);
    }
    cpxSolver.branchAndBound();
    double timeTaken = CoinCpuTime() - time1;
    sprintf(printBuffer, "Cplex took %g seconds",
      timeTaken);
    messageHandler()->message(CBC_GENERAL, messages())
      << printBuffer << CoinMessageEol;
    numberExtraNodes_ = CPXgetnodecnt(env, lpPtr);
    numberExtraIterations_ = CPXgetmipitcnt(env, lpPtr);
    double value = cpxSolver.getObjValue() * direction;
    if (cpxSolver.isProvenOptimal() && value <= getCutoff()) {
      feasible = true;
      clpSolver->setWarmStart(NULL);
      // try and do solution
      double *newSolution = CoinCopyOfArray(cpxSolver.getColSolution(),
        getNumCols());
      setBestSolution(CBC_STRONGSOL, value, newSolution);
      delete[] newSolution;
    }
    feasible = false;
  }
#endif
  if (!parentModel_ && (moreSpecialOptions_ & 268435456) != 0) {
    // try idiotic idea
    CbcObject *obj = new CbcIdiotBranch(this);
    obj->setPriority(1); // temp
    addObjects(1, &obj);
    delete obj;
  }

  /*
      A hook to use clp to quickly explore some part of the tree.
    */
  if (fastNodeDepth_ == 1000 && /*!parentModel_*/ (specialOptions_ & 2048) == 0) {
    fastNodeDepth_ = -1;
    CbcObject *obj = new CbcFollowOn(this);
    obj->setPriority(1);
    addObjects(1, &obj);
    delete obj;
  }
  int saveNumberSolves = numberSolves_;
  int saveNumberIterations = numberIterations_;
  if ((fastNodeDepth_ >= 0 || (moreSpecialOptions_ & 33554432) != 0)
    && /*!parentModel_*/ (specialOptions_ & 2048) == 0) {
    // add in a general depth object doClp
    int type = (fastNodeDepth_ <= 100) ? fastNodeDepth_ : -(fastNodeDepth_ - 100);
    if ((moreSpecialOptions_ & 33554432) != 0)
      type = 12;
    else
      fastNodeDepth_ += 1000000; // mark as done
    CbcObject *obj = new CbcGeneralDepth(this, type);
    addObjects(1, &obj);
    delete obj;
    // fake number of objects
    numberObjects_--;
    if (parallelMode() < -1) {
      // But make sure position is correct
      OsiObject *obj2 = object_[numberObjects_];
      obj = dynamic_cast< CbcObject * >(obj2);
      assert(obj);
      obj->setPosition(numberObjects_);
    }
  }
#ifdef COIN_HAS_CLP
#ifdef NO_CRUNCH
  if (true) {
    OsiClpSolverInterface *clpSolver
      = dynamic_cast< OsiClpSolverInterface * >(solver_);
    if (clpSolver && !parentModel_) {
      ClpSimplex *clpSimplex = clpSolver->getModelPtr();
      clpSimplex->setSpecialOptions(clpSimplex->specialOptions() | 131072);
      //clpSimplex->startPermanentArrays();
      clpSimplex->setPersistenceFlag(2);
    }
  }
#endif
#endif
  // Save copy of solver
  OsiSolverInterface *saveSolver = NULL;
  if (!parentModel_ && (specialOptions_ & (512 + 32768)) != 0)
    saveSolver = solver_->clone();
  double checkCutoffForRestart = 1.0e100;
  saveModel(saveSolver, &checkCutoffForRestart, &feasible);
  if ((specialOptions_ & 262144) != 0 && !parentModel_) {
    // Save stuff and return!
    storedRowCuts_->saveStuff(bestObjective_, bestSolution_,
      solver_->getColLower(),
      solver_->getColUpper());
    delete[] lowerBefore;
    delete[] upperBefore;
    delete saveSolver;
    if (flipObjective)
      flipModel();
    return;
  }
  /*
      We've taken the continuous relaxation as far as we can. Time to branch.
      The first order of business is to actually create a node. chooseBranch
      currently uses strong branching to evaluate branch object candidates,
      unless forced back to simple branching. If chooseBranch concludes that a
      branching candidate is monotone (anyAction == -1) or infeasible (anyAction
      == -2) when forced to integer values, it returns here immediately.

      Monotone variables trigger a call to resolve(). If the problem remains
      feasible, try again to choose a branching variable. At the end of the loop,
      resolved == true indicates that some variables were fixed.

      Loss of feasibility will result in the deletion of newNode.
    */

  bool resolved = false;
  CbcNode *newNode = NULL;
  numberFixedAtRoot_ = 0;
  numberFixedNow_ = 0;
  if (!parentModel_ && (moreSpecialOptions2_ & 2) != 0) {
#ifdef COIN_HAS_CLP
    OsiClpSolverInterface *clpSolver
      = dynamic_cast< OsiClpSolverInterface * >(solver_);
    if (clpSolver) {
      if (getCutoff() > 1.0e20) {
        printf("Zapping costs\n");
        int numberColumns = solver_->getNumCols();
        double *zeroCost = new double[numberColumns];
        // could make small random
        memset(zeroCost, 0, numberColumns * sizeof(double));
        solver_->setObjective(zeroCost);
        double objValue = solver_->getObjValue();
        solver_->setDblParam(OsiObjOffset, -objValue);
        clpSolver->getModelPtr()->setObjectiveValue(objValue);
        delete[] zeroCost;
      } else {
        moreSpecialOptions2_ &= ~2;
      }
    } else {
#endif
      moreSpecialOptions2_ &= ~2;
#ifdef COIN_HAS_CLP
    }
#endif
  }
  int numberIterationsAtContinuous = numberIterations_;
  //solverCharacteristics_->setSolver(solver_);
  if (feasible) {
    // mark all cuts as globally valid
    int numberCuts = cuts.sizeRowCuts();
    resizeWhichGenerator(0, numberCuts);
    for (int i = 0; i < numberCuts; i++) {
      cuts.rowCutPtr(i)->setGloballyValid();
      whichGenerator_[i] = 20000 + (whichGenerator_[i] % 10000);
    }
#define HOTSTART -1
#if HOTSTART < 0
    if (bestSolution_ && !parentModel_ && !hotstartSolution_ && (moreSpecialOptions_ & 1024) != 0 && (specialOptions_ & 2048) == 0) {
      // Set priorities so only branch on ones we need to
      // use djs and see if only few branches needed
#ifndef NDEBUG
      double integerTolerance = getIntegerTolerance();
#endif
      bool possible = true;
      const double *saveLower = continuousSolver_->getColLower();
      const double *saveUpper = continuousSolver_->getColUpper();
      for (int i = 0; i < numberObjects_; i++) {
        const CbcSimpleInteger *thisOne = dynamic_cast< const CbcSimpleInteger * >(object_[i]);
        if (thisOne) {
          int iColumn = thisOne->columnNumber();
          if (saveUpper[iColumn] > saveLower[iColumn] + 1.5) {
            possible = false;
            break;
          }
        } else {
          possible = false;
          break;
        }
      }
      if (possible) {
        OsiSolverInterface *solver = continuousSolver_->clone();
        int numberColumns = solver->getNumCols();
        for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
          double value = bestSolution_[iColumn];
          value = CoinMax(value, saveLower[iColumn]);
          value = CoinMin(value, saveUpper[iColumn]);
          value = floor(value + 0.5);
          if (solver->isInteger(iColumn)) {
            solver->setColLower(iColumn, value);
            solver->setColUpper(iColumn, value);
          }
        }
        solver->setHintParam(OsiDoDualInResolve, false, OsiHintTry);
        // objlim and all slack
        double direction = solver->getObjSense();
        solver->setDblParam(OsiDualObjectiveLimit, 1.0e50 * direction);
        CoinWarmStartBasis *basis = dynamic_cast< CoinWarmStartBasis * >(solver->getEmptyWarmStart());
        solver->setWarmStart(basis);
        delete basis;
        bool changed = true;
        hotstartPriorities_ = new int[numberColumns];
        for (int iColumn = 0; iColumn < numberColumns; iColumn++)
          hotstartPriorities_[iColumn] = 1;
        while (changed) {
          changed = false;
          solver->resolve();
          if (!solver->isProvenOptimal()) {
            possible = false;
            break;
          }
          const double *dj = solver->getReducedCost();
          const double *colLower = solver->getColLower();
          const double *colUpper = solver->getColUpper();
          const double *solution = solver->getColSolution();
          int nAtLbNatural = 0;
          int nAtUbNatural = 0;
          int nZeroDj = 0;
          int nForced = 0;
          for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
            double value = solution[iColumn];
            value = CoinMax(value, saveLower[iColumn]);
            value = CoinMin(value, saveUpper[iColumn]);
            if (solver->isInteger(iColumn)) {
              assert(fabs(value - solution[iColumn]) <= integerTolerance);
              if (hotstartPriorities_[iColumn] == 1) {
                if (dj[iColumn] < -1.0e-6) {
                  // negative dj
                  if (saveUpper[iColumn] == colUpper[iColumn]) {
                    nAtUbNatural++;
                    hotstartPriorities_[iColumn] = 2;
                    solver->setColLower(iColumn, saveLower[iColumn]);
                    solver->setColUpper(iColumn, saveUpper[iColumn]);
                  } else {
                    nForced++;
                  }
                } else if (dj[iColumn] > 1.0e-6) {
                  // positive dj
                  if (saveLower[iColumn] == colLower[iColumn]) {
                    nAtLbNatural++;
                    hotstartPriorities_[iColumn] = 2;
                    solver->setColLower(iColumn, saveLower[iColumn]);
                    solver->setColUpper(iColumn, saveUpper[iColumn]);
                  } else {
                    nForced++;
                  }
                } else {
                  // zero dj
                  nZeroDj++;
                }
              }
            }
          }
#if CBC_USEFUL_PRINTING > 1
          printf("%d forced, %d naturally at lower, %d at upper - %d zero dj\n",
            nForced, nAtLbNatural, nAtUbNatural, nZeroDj);
#endif
          if (nAtLbNatural || nAtUbNatural) {
            changed = true;
          } else {
            if (nForced + nZeroDj > 5000 || (nForced + nZeroDj) * 2 > numberIntegers_)
              possible = false;
          }
        }
        delete solver;
      }
      if (possible) {
        setHotstartSolution(bestSolution_);
        if (!saveCompare) {
          // create depth first comparison
          saveCompare = nodeCompare_;
          // depth first
          nodeCompare_ = new CbcCompareDepth();
          tree_->setComparison(*nodeCompare_);
        }
      } else {
        delete[] hotstartPriorities_;
        hotstartPriorities_ = NULL;
      }
    }
#endif
#if HOTSTART > 0
    if (hotstartSolution_ && !hotstartPriorities_) {
      // Set up hot start
      OsiSolverInterface *solver = solver_->clone();
      double direction = solver_->getObjSense();
      int numberColumns = solver->getNumCols();
      double *saveLower = CoinCopyOfArray(solver->getColLower(), numberColumns);
      double *saveUpper = CoinCopyOfArray(solver->getColUpper(), numberColumns);
      // move solution
      solver->setColSolution(hotstartSolution_);
      // point to useful information
      const double *saveSolution = testSolution_;
      testSolution_ = solver->getColSolution();
      OsiBranchingInformation usefulInfo = usefulInformation();
      testSolution_ = saveSolution;
      /*
            Run through the objects and use feasibleRegion() to set variable bounds
            so as to fix the variables specified in the objects at their value in this
            solution. Since the object list contains (at least) one object for every
            integer variable, this has the effect of fixing all integer variables.
            */
      for (int i = 0; i < numberObjects_; i++)
        object_[i]->feasibleRegion(solver, &usefulInfo);
      solver->resolve();
      assert(solver->isProvenOptimal());
      double gap = CoinMax((solver->getObjValue() - solver_->getObjValue()) * direction, 0.0);
      const double *dj = solver->getReducedCost();
      const double *colLower = solver->getColLower();
      const double *colUpper = solver->getColUpper();
      const double *solution = solver->getColSolution();
      int nAtLbNatural = 0;
      int nAtUbNatural = 0;
      int nAtLbNaturalZero = 0;
      int nAtUbNaturalZero = 0;
      int nAtLbFixed = 0;
      int nAtUbFixed = 0;
      int nAtOther = 0;
      int nAtOtherNatural = 0;
      int nNotNeeded = 0;
      delete[] hotstartSolution_;
      hotstartSolution_ = new double[numberColumns];
      delete[] hotstartPriorities_;
      hotstartPriorities_ = new int[numberColumns];
      int *order = (int *)saveUpper;
      int nFix = 0;
      double bestRatio = COIN_DBL_MAX;
      for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
        double value = solution[iColumn];
        value = CoinMax(value, saveLower[iColumn]);
        value = CoinMin(value, saveUpper[iColumn]);
        double sortValue = COIN_DBL_MAX;
        if (solver->isInteger(iColumn)) {
          assert(fabs(value - solution[iColumn]) <= 1.0e-5);
          double value2 = floor(value + 0.5);
          if (dj[iColumn] < -1.0e-6) {
            // negative dj
            //assert (value2==colUpper[iColumn]);
            if (saveUpper[iColumn] == colUpper[iColumn]) {
              nAtUbNatural++;
              sortValue = 0.0;
              double value = -dj[iColumn];
              if (value > gap)
                nFix++;
              else if (gap < value * bestRatio)
                bestRatio = gap / value;
              if (saveLower[iColumn] != colLower[iColumn]) {
                nNotNeeded++;
                sortValue = 1.0e20;
              }
            } else if (saveLower[iColumn] == colUpper[iColumn]) {
              nAtLbFixed++;
              sortValue = dj[iColumn];
            } else {
              nAtOther++;
              sortValue = 0.0;
              if (saveLower[iColumn] != colLower[iColumn] && saveUpper[iColumn] != colUpper[iColumn]) {
                nNotNeeded++;
                sortValue = 1.0e20;
              }
            }
          } else if (dj[iColumn] > 1.0e-6) {
            // positive dj
            //assert (value2==colLower[iColumn]);
            if (saveLower[iColumn] == colLower[iColumn]) {
              nAtLbNatural++;
              sortValue = 0.0;
              double value = dj[iColumn];
              if (value > gap)
                nFix++;
              else if (gap < value * bestRatio)
                bestRatio = gap / value;
              if (saveUpper[iColumn] != colUpper[iColumn]) {
                nNotNeeded++;
                sortValue = 1.0e20;
              }
            } else if (saveUpper[iColumn] == colLower[iColumn]) {
              nAtUbFixed++;
              sortValue = -dj[iColumn];
            } else {
              nAtOther++;
              sortValue = 0.0;
              if (saveLower[iColumn] != colLower[iColumn] && saveUpper[iColumn] != colUpper[iColumn]) {
                nNotNeeded++;
                sortValue = 1.0e20;
              }
            }
          } else {
            // zero dj
            if (value2 == saveUpper[iColumn]) {
              nAtUbNaturalZero++;
              sortValue = 0.0;
              if (saveLower[iColumn] != colLower[iColumn]) {
                nNotNeeded++;
                sortValue = 1.0e20;
              }
            } else if (value2 == saveLower[iColumn]) {
              nAtLbNaturalZero++;
              sortValue = 0.0;
            } else {
              nAtOtherNatural++;
              sortValue = 0.0;
              if (saveLower[iColumn] != colLower[iColumn] && saveUpper[iColumn] != colUpper[iColumn]) {
                nNotNeeded++;
                sortValue = 1.0e20;
              }
            }
          }
#if HOTSTART == 3
          sortValue = -fabs(dj[iColumn]);
#endif
        }
        hotstartSolution_[iColumn] = value;
        saveLower[iColumn] = sortValue;
        order[iColumn] = iColumn;
      }
      COIN_DETAIL_PRINT(printf("** can fix %d columns - best ratio for others is %g on gap of %g\n",
        nFix, bestRatio, gap));
      int nNeg = 0;
      CoinSort_2(saveLower, saveLower + numberColumns, order);
      for (int i = 0; i < numberColumns; i++) {
        if (saveLower[i] < 0.0) {
          nNeg++;
#if HOTSTART == 2 || HOTSTART == 3
          // swap sign ?
          saveLower[i] = -saveLower[i];
#endif
        }
      }
      CoinSort_2(saveLower, saveLower + nNeg, order);
      for (int i = 0; i < numberColumns; i++) {
#if HOTSTART == 1
        hotstartPriorities_[order[i]] = 100;
#else
          hotstartPriorities_[order[i]] = -(i + 1);
#endif
      }
      COIN_DETAIL_PRINT(printf("nAtLbNat %d,nAtUbNat %d,nAtLbNatZero %d,nAtUbNatZero %d,nAtLbFixed %d,nAtUbFixed %d,nAtOther %d,nAtOtherNat %d, useless %d %d\n",
        nAtLbNatural,
        nAtUbNatural,
        nAtLbNaturalZero,
        nAtUbNaturalZero,
        nAtLbFixed,
        nAtUbFixed,
        nAtOther,
        nAtOtherNatural, nNotNeeded, nNeg));
      delete[] saveLower;
      delete[] saveUpper;
      if (!saveCompare) {
        // create depth first comparison
        saveCompare = nodeCompare_;
        // depth first
        nodeCompare_ = new CbcCompareDepth();
        tree_->setComparison(*nodeCompare_);
      }
    }
#endif
    newNode = new CbcNode;
    // Set objective value (not so obvious if NLP etc)
    setObjectiveValue(newNode, NULL);
    anyAction = -1;
    // To make depth available we may need a fake node
    CbcNode fakeNode;
    phase_ = 3;
    // only allow 1000 passes
    int numberPassesLeft = 1000;
    // This is first crude step
    if (numberAnalyzeIterations_ && !parentModel_) {
      delete[] analyzeResults_;
      //int numberColumns = solver_->getNumCols();
      analyzeResults_ = new double[5 * numberIntegers_];
      numberFixedAtRoot_ = newNode->analyze(this, analyzeResults_);
      if (numberFixedAtRoot_ > 0) {
        COIN_DETAIL_PRINT(printf("%d fixed by analysis\n", numberFixedAtRoot_));
        setPointers(solver_);
        numberFixedNow_ = numberFixedAtRoot_;
      } else if (numberFixedAtRoot_ < 0) {
        COIN_DETAIL_PRINT(printf("analysis found to be infeasible\n"));
        anyAction = -2;
        delete newNode;
        newNode = NULL;
        feasible = false;
      }
    }
    OsiSolverBranch *branches = NULL;
    if (feasible)
      anyAction = chooseBranch(newNode, numberPassesLeft, NULL, cuts, resolved,
        NULL, NULL, NULL, branches);
    if (anyAction == -2 || newNode->objectiveValue() >= cutoff) {
      if (anyAction != -2) {
        // zap parent nodeInfo
#ifdef COIN_DEVELOP
        printf("zapping CbcNodeInfo %x\n", newNode->nodeInfo()->parent());
#endif
        if (newNode->nodeInfo())
          newNode->nodeInfo()->nullParent();
      }
      deleteNode(newNode);
      newNode = NULL;
      feasible = false;
    }
  }
  if (newNode && probingInfo_) {
    int number01 = probingInfo_->numberIntegers();
    //const fixEntry * entry = probingInfo_->fixEntries();
    const int *toZero = probingInfo_->toZero();
    //const int * toOne = probingInfo_->toOne();
    //const int * integerVariable = probingInfo_->integerVariable();
    if (toZero[number01]) {
      CglTreeProbingInfo info(*probingInfo_);
      if ((moreSpecialOptions2_ & 64) != 0 && !parentModel_) {
        /*
		Marginal idea. Further exploration probably good. Build some extra
		cliques from probing info. Not quite worth the effort?
	      */
        CglProbing generator1;
        generator1.setUsingObjective(false);
        generator1.setMaxPass(1);
        generator1.setMaxPassRoot(1);
        generator1.setMaxLook(100);
        generator1.setRowCuts(3);
        generator1.setMaxElements(300);
        generator1.setMaxProbeRoot(solver_->getNumCols());
        CoinThreadRandom randomGenerator;
        //CglTreeProbingInfo info(solver_);
        info.level = 0;
        info.formulation_rows = solver_->getNumRows();
        info.inTree = false;
        if (parentModel_) {
          info.parentSolver = parentModel_->continuousSolver();
          // indicate if doing full search
          info.hasParent = ((specialOptions_ & 67108864) == 0) ? 1 : 2;
        } else {
          info.hasParent = 0;
          info.parentSolver = NULL;
        }
        info.originalColumns = originalColumns();
        info.randomNumberGenerator = &randomGenerator;
        info.pass = 4;
        generator1.setMode(8);
        OsiCuts cs;
        generator1.generateCutsAndModify(*solver_, cs, &info);
        // very clunky
        OsiSolverInterface *temp = generator1.cliqueModel(solver_, 2);
        CglPreProcess dummy;
        OsiSolverInterface *newSolver = dummy.cliqueIt(*temp, 0.0001);
        delete temp;
        OsiSolverInterface *fake = NULL;
        if (newSolver) {
#if 0
		int numberCliques = generator1.numberCliques();
		cliqueEntry * entry = generator1.cliqueEntry();
		cliqueType * type = new cliqueType [numberCliques];
		int * start = new int [numberCliques+1];
		start[numberCliques]=2*numberCliques;
		int n=0;
		for (int i=0;i<numberCliques;i++) {
		  start[i]=2*i;
		  setOneFixesInCliqueEntry(entry[2*i],true);
		  setOneFixesInCliqueEntry(entry[2*i+1],true);
		  type[i]=0;
		}
		fake = info.analyze(*solver_, 1,numberCliques,start,
				    entry,type);
		delete [] type;
		delete [] entry;
#else
        fake = info.analyze(*newSolver, 1, -1);
#endif
          delete newSolver;
        } else {
          fake = info.analyze(*solver_, 1);
        }
        if (fake) {
          //fake->writeMps("fake");
          CglFakeClique cliqueGen(fake);
          cliqueGen.setStarCliqueReport(false);
          cliqueGen.setRowCliqueReport(false);
          cliqueGen.setMinViolation(0.1);
          addCutGenerator(&cliqueGen, 1, "Fake cliques", true, false, false, -200);
          generator_[numberCutGenerators_ - 1]->setTiming(true);
          for (int i = 0; i < numberCutGenerators_; i++) {
            CglKnapsackCover *cutGen = dynamic_cast< CglKnapsackCover * >(generator_[i]->generator());
            if (cutGen) {
              cutGen->createCliques(*fake, 2, 200, false);
            }
          }
        }
      }
      if (probingInfo_->packDown()) {
#if CBC_USEFUL_PRINTING > 1
        printf("%d implications on %d 0-1\n", toZero[number01], number01);
#endif
        // Create a cut generator that remembers implications discovered at root.
        CglImplication implication(probingInfo_);
        addCutGenerator(&implication, 1, "ImplicationCuts", true, false, false, -200);
        generator_[numberCutGenerators_ - 1]->setGlobalCuts(true);
        generator_[numberCutGenerators_ - 1]->setTiming(true);
      } else {
        delete probingInfo_;
        probingInfo_ = NULL;
      }
    } else {
      delete probingInfo_;

      probingInfo_ = NULL;
    }
  }
  /*
      At this point, the root subproblem is infeasible or fathomed by bound
      (newNode == NULL), or we're live with an objective value that satisfies the
      current objective cutoff.
    */
  assert(!newNode || newNode->objectiveValue() <= cutoff);
  // Save address of root node as we don't want to delete it
  /*
      The common case is that the lp relaxation is feasible but doesn't satisfy
      integrality (i.e., newNode->branchingObject(), indicating we've been able to
      select a branching variable). Remove any cuts that have gone slack due to
      forcing monotone variables. Then tack on an CbcFullNodeInfo object and full
      basis (via createInfo()) and stash the new cuts in the nodeInfo (via
      addCuts()). If, by some miracle, we have an integral solution at the root
      (newNode->branchingObject() is NULL), takeOffCuts() will ensure that the solver holds
      a valid solution for use by setBestSolution().
    */
  CoinWarmStartBasis *lastws = NULL;
  if (feasible && newNode->branchingObject()) {
    if (resolved) {
      takeOffCuts(cuts, false, NULL);
#ifdef CHECK_CUT_COUNTS
      {
        printf("Number of rows after chooseBranch fix (root)"
               "(active only) %d\n",
          numberRowsAtContinuous_ + numberNewCuts_ + numberOldActiveCuts_);
        const CoinWarmStartBasis *debugws = dynamic_cast< const CoinWarmStartBasis * >(solver_->getWarmStart());
        debugws->print();
        delete debugws;
      }
#endif
    }
    //newNode->createInfo(this,NULL,NULL,NULL,NULL,0,0) ;
    //newNode->nodeInfo()->addCuts(cuts,
    //			 newNode->numberBranches(),whichGenerator_) ;
    if (lastws)
      delete lastws;
    lastws = dynamic_cast< CoinWarmStartBasis * >(solver_->getWarmStart());
  }
  /*
      Continuous data to be used later
    */
  continuousObjective_ = solver_->getObjValue() * solver_->getObjSense();
  continuousInfeasibilities_ = 0;
  if (newNode) {
    continuousObjective_ = newNode->objectiveValue();
    delete[] continuousSolution_;
    continuousSolution_ = CoinCopyOfArray(solver_->getColSolution(),
      numberColumns);
    continuousInfeasibilities_ = newNode->numberUnsatisfied();
  }
  /*
      Bound may have changed so reset in objects
    */
  {
    int i;
    for (i = 0; i < numberObjects_; i++)
      object_[i]->resetBounds(solver_);
  }
  /*
      Feasible? Then we should have either a live node prepped for future
      expansion (indicated by variable() >= 0), or (miracle of miracles) an
      integral solution at the root node.

      initializeInfo sets the reference counts in the nodeInfo object.  Since
      this node is still live, push it onto the heap that holds the live set.
    */
  if (newNode) {
    if (newNode->branchingObject()) {
      newNode->initializeInfo();
      tree_->push(newNode);
      // save pointer to root node - so can pick up bounds
      if (!topOfTree_)
        topOfTree_ = dynamic_cast< CbcFullNodeInfo * >(newNode->nodeInfo());
      if (statistics_) {
        if (numberNodes2_ == maximumStatistics_) {
          maximumStatistics_ = 2 * maximumStatistics_;
          CbcStatistics **temp = new CbcStatistics *[maximumStatistics_];
          memset(temp, 0, maximumStatistics_ * sizeof(CbcStatistics *));
          memcpy(temp, statistics_, numberNodes2_ * sizeof(CbcStatistics *));
          delete[] statistics_;
          statistics_ = temp;
        }
        assert(!statistics_[numberNodes2_]);
        statistics_[numberNodes2_] = new CbcStatistics(newNode, this);
      }
      numberNodes2_++;
#ifdef CHECK_NODE
      printf("Node %x on tree\n", newNode);
#endif
    } else {
      // continuous is integer
      double objectiveValue = newNode->objectiveValue();
      setBestSolution(CBC_SOLUTION, objectiveValue,
        solver_->getColSolution());
      if (eventHandler) {
        // we are stopping anyway so no need to test return code
        eventHandler->event(CbcEventHandler::solution);
      }
      delete newNode;
      newNode = NULL;
    }
  }

  if (printFrequency_ <= 0) {
    printFrequency_ = 1000;
    if (getNumCols() > 2000)
      printFrequency_ = 100;
  }
  /*
      It is possible that strong branching fixes one variable and then the code goes round
      again and again.  This can take too long.  So we need to warn user - just once.
    */
  numberLongStrong_ = 0;
  CbcNode *createdNode = NULL;
#ifdef CBC_THREAD
  if ((specialOptions_ & 2048) != 0)
    numberThreads_ = 0;
  if (numberThreads_) {
    nodeCompare_->sayThreaded(); // need to use addresses
    master_ = new CbcBaseModel(*this,
      (parallelMode() < -1) ? 1 : 0);
    masterThread_ = master_->masterThread();
  }
#endif
#ifdef COIN_HAS_CLP
  {
    OsiClpSolverInterface *clpSolver
      = dynamic_cast< OsiClpSolverInterface * >(solver_);
    if (clpSolver && !parentModel_) {
      clpSolver->computeLargestAway();
    }
  }
#endif
  /*
      At last, the actual branch-and-cut search loop, which will iterate until
      the live set is empty or we hit some limit (integrality gap, time, node
      count, etc.). The overall flow is to rebuild a subproblem, reoptimise using
      solveWithCuts(), choose a branching pattern with chooseBranch(), and finally
      add the node to the live set.

      The first action is to winnow the live set to remove nodes which are worse
      than the current objective cutoff.
    */
  if (solver_->getRowCutDebuggerAlways()) {
    OsiRowCutDebugger *debuggerX = const_cast< OsiRowCutDebugger * >(solver_->getRowCutDebuggerAlways());
    const OsiRowCutDebugger *debugger = solver_->getRowCutDebugger();
    if (!debugger) {
      // infeasible!!
      printf("before search\n");
      const double *lower = solver_->getColLower();
      const double *upper = solver_->getColUpper();
      const double *solution = debuggerX->optimalSolution();
      int numberColumns = solver_->getNumCols();
      for (int i = 0; i < numberColumns; i++) {
        if (solver_->isInteger(i)) {
          if (solution[i] < lower[i] - 1.0e-6 || solution[i] > upper[i] + 1.0e-6)
            printf("**** ");
          printf("%d %g <= %g <= %g\n",
            i, lower[i], solution[i], upper[i]);
        }
      }
      //abort();
    }
  }
  {
    // may be able to change cutoff now
    double cutoff = getCutoff();
    double increment = getDblParam(CbcModel::CbcCutoffIncrement);
    if (cutoff > bestObjective_ - increment) {
      cutoff = bestObjective_ - increment;
      setCutoff(cutoff);
    }
  }
#ifdef CBC_THREAD
  bool goneParallel = false;
#endif
#define MAX_DEL_NODE 1
  CbcNode *delNode[MAX_DEL_NODE + 1];
  int nDeleteNode = 0;
  // For Printing etc when parallel
  int lastEvery1000 = 0;
  int lastPrintEvery = 0;
  int numberConsecutiveInfeasible = 0;
#define PERTURB_IN_FATHOM
#ifdef PERTURB_IN_FATHOM
  // allow in fathom
  if ((moreSpecialOptions_ & 262144) != 0)
    specialOptions_ |= 131072;
#endif
  while (true) {
    lockThread();
#ifdef COIN_HAS_CLP
    // See if we want dantzig row choice
    goToDantzig(100, savePivotMethod);
#endif
    //#define REPORT_DYNAMIC 2
#if REPORT_DYNAMIC
    if (numberNodes_ && !parentModel_ && (tree_->empty() || (numberNodes_ % 10000) == 0)) {
      // Get average up and down costs
      double averageUp = 0.0;
      double averageDown = 0.0;
      int numberUp = 0;
      int numberDown = 0;
      int minTimesDown = COIN_INT_MAX;
      int maxTimesDown = 0;
      int neverBranchedDown = 0;
      int infeasibleTimesDown = 0;
      int minTimesUp = COIN_INT_MAX;
      int maxTimesUp = 0;
      int infeasibleTimesUp = 0;
      int neverBranchedUp = 0;
      int neverBranched = 0;
      int i;
      int numberInts = 0;
      bool endOfSearch = tree_->empty();
      int numberUp2 = 0;
      int numberDown2 = 0;
      for (i = 0; i < numberObjects_; i++) {
        OsiObject *object = object_[i];
        CbcSimpleIntegerDynamicPseudoCost *dynamicObject = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object);
        if (dynamicObject) {
          numberInts++;
          if (dynamicObject->numberTimesUp() || dynamicObject->numberTimesDown()) {
            int nUp = 0;
            int nDown = 0;
            double up = 0.0;
            double down = 0.0;
            if (dynamicObject->numberTimesUp()) {
              numberUp++;
              nUp = dynamicObject->numberTimesUp();
              minTimesUp = CoinMin(minTimesUp, nUp);
              maxTimesUp = CoinMax(maxTimesUp, nUp);
              up = dynamicObject->upDynamicPseudoCost();
              averageUp += up;
              numberUp2 += nUp;
              infeasibleTimesUp += dynamicObject->numberTimesUpInfeasible();
            } else {
              neverBranchedUp++;
            }
            if (dynamicObject->numberTimesDown()) {
              numberDown++;
              nDown = dynamicObject->numberTimesDown();
              minTimesDown = CoinMin(minTimesDown, nDown);
              maxTimesDown = CoinMax(maxTimesDown, nDown);
              down = dynamicObject->downDynamicPseudoCost();
              averageDown += down;
              numberDown2 += dynamicObject->numberTimesDown();
              infeasibleTimesDown += dynamicObject->numberTimesDownInfeasible();
            } else {
              neverBranchedDown++;
            }
#if REPORT_DYNAMIC > 1
#if REPORT_DYNAMIC == 2
            if (endOfSearch && numberIntegers_ < 400) {
#elif REPORT_DYNAMIC == 3
            if (endOfSearch) {
#else
            {
#endif
              dynamicObject->print(0, 0.0);
            }
#endif
          } else {
            neverBranched++;
#if REPORT_DYNAMIC > 2
#if REPORT_DYNAMIC == 3
            if (endOfSearch && numberIntegers_ < 400) {
#elif REPORT_DYNAMIC == 4
            if (endOfSearch) {
#else
            {
#endif
              printf("col %d - never branched on\n", dynamicObject->columnNumber());
            }
#endif
          }
        }
      }
      if (numberUp)
        averageUp /= static_cast< double >(numberUp);
      else
        averageUp = 0.0;
      if (numberDown)
        averageDown /= static_cast< double >(numberDown);
      else
        averageDown = 0.0;
      printf("Report for %d variables (%d never branched on) after %d nodes - total solves down %d up %d\n",
        numberInts, neverBranched, numberNodes_, numberDown2, numberUp2);
      if ((neverBranchedDown || neverBranchedUp) && endOfSearch)
        printf("odd %d never branched down and %d never branched up\n",
          neverBranchedDown, neverBranchedUp);
      printf("down average %g times (%d infeasible) average increase %g min/max times (%d,%d)\n",
        static_cast< double >(numberDown2) / numberDown, infeasibleTimesDown, averageDown,
        minTimesDown, maxTimesDown);
      printf("up average %g times (%d infeasible) average increase %g min/max times (%d,%d)\n",
        static_cast< double >(numberUp2) / numberUp, infeasibleTimesUp, averageUp,
        minTimesUp, maxTimesUp);
    }
#endif
    if (tree_->empty()) {
#ifdef CBC_THREAD
      if (parallelMode() > 0 && master_) {
        int anyLeft = master_->waitForThreadsInTree(0);
        if (!anyLeft) {
          master_->stopThreads(-1);
          break;
        }
      } else {
        break;
      }
#else
    break;
#endif
    } else {
      unlockThread();
    }
    // If done 50/100 nodes see if worth trying reduction
    if (numberNodes_ >= nextCheckRestart) {
      if (nextCheckRestart < 100)
        nextCheckRestart = 100;
      else
        nextCheckRestart = COIN_INT_MAX;
#ifdef COIN_HAS_CLP
      OsiClpSolverInterface *clpSolver
        = dynamic_cast< OsiClpSolverInterface * >(solver_);
      if (clpSolver && ((specialOptions_ & 131072) == 0) && true) {
        ClpSimplex *simplex = clpSolver->getModelPtr();
        int perturbation = simplex->perturbation();
#if CBC_USEFUL_PRINTING > 1
        printf("Testing its n,s %d %d solves n,s %d %d - pert %d\n",
          numberIterations_, saveNumberIterations,
          numberSolves_, saveNumberSolves, perturbation);
#endif
        if (perturbation == 50 && (numberIterations_ - saveNumberIterations) < 8 * (numberSolves_ - saveNumberSolves)) {
          // switch off perturbation
          simplex->setPerturbation(100);
#if CBC_USEFUL_PRINTING > 1
          printf("Perturbation switched off\n");
#endif
        }
      }
#endif
      /*
              Decide if we want to do a restart.
            */
      if (saveSolver && (specialOptions_ & (512 + 32768)) != 0) {
        bool tryNewSearch = solverCharacteristics_->reducedCostsAccurate() && (getCutoff() < 1.0e20 && getCutoff() < checkCutoffForRestart);
        int numberColumns = getNumCols();
        if (tryNewSearch) {
          // adding increment back allows current best - tiny bit weaker
          checkCutoffForRestart = getCutoff() + getCutoffIncrement();
#if CBC_USEFUL_PRINTING > 1
          printf("after %d nodes, cutoff %g - looking\n",
            numberNodes_, getCutoff());
#endif
          saveSolver->resolve();
          double direction = saveSolver->getObjSense();
          double gap = checkCutoffForRestart - saveSolver->getObjValue() * direction;
          double tolerance;
          saveSolver->getDblParam(OsiDualTolerance, tolerance);
          if (gap <= 0.0)
            gap = tolerance;
          gap += 100.0 * tolerance;
          double integerTolerance = getDblParam(CbcIntegerTolerance);

          const double *lower = saveSolver->getColLower();
          const double *upper = saveSolver->getColUpper();
          const double *solution = saveSolver->getColSolution();
          const double *reducedCost = saveSolver->getReducedCost();

          int numberFixed = 0;
          int numberFixed2 = 0;
#ifdef COIN_DEVELOP
          printf("gap %g\n", gap);
#endif
          for (int i = 0; i < numberIntegers_; i++) {
            int iColumn = integerVariable_[i];
            double djValue = direction * reducedCost[iColumn];
            if (upper[iColumn] - lower[iColumn] > integerTolerance) {
              if (solution[iColumn] < lower[iColumn] + integerTolerance && djValue > gap) {
                //printf("%d to lb on dj of %g - bounds %g %g\n",
                //     iColumn,djValue,lower[iColumn],upper[iColumn]);
                saveSolver->setColUpper(iColumn, lower[iColumn]);
                numberFixed++;
              } else if (solution[iColumn] > upper[iColumn] - integerTolerance && -djValue > gap) {
                //printf("%d to ub on dj of %g - bounds %g %g\n",
                //     iColumn,djValue,lower[iColumn],upper[iColumn]);
                saveSolver->setColLower(iColumn, upper[iColumn]);
                numberFixed++;
              }
            } else {
              //printf("%d has dj of %g - already fixed to %g\n",
              //     iColumn,djValue,lower[iColumn]);
              numberFixed2++;
            }
          }
#ifdef COIN_DEVELOP
          if ((specialOptions_ & 1) != 0) {
            const OsiRowCutDebugger *debugger = saveSolver->getRowCutDebugger();
            if (debugger) {
              printf("Contains optimal\n");
              OsiSolverInterface *temp = saveSolver->clone();
              const double *solution = debugger->optimalSolution();
              const double *lower = temp->getColLower();
              const double *upper = temp->getColUpper();
              int n = temp->getNumCols();
              for (int i = 0; i < n; i++) {
                if (temp->isInteger(i)) {
                  double value = floor(solution[i] + 0.5);
                  assert(value >= lower[i] && value <= upper[i]);
                  temp->setColLower(i, value);
                  temp->setColUpper(i, value);
                }
              }
              temp->writeMps("reduced_fix");
              delete temp;
              saveSolver->writeMps("reduced");
            } else {
              abort();
            }
          }
          printf("Restart could fix %d integers (%d already fixed)\n",
            numberFixed + numberFixed2, numberFixed2);
#endif
          numberFixed += numberFixed2;
          if (numberFixed * 10 < numberColumns && numberFixed * 4 < numberIntegers_)
            tryNewSearch = false;
        }
#ifdef CONFLICT_CUTS
        // temporary
        //if ((moreSpecialOptions_&4194304)!=0)
        //tryNewSearch=false;
#endif
        if (tryNewSearch) {
          // back to solver without cuts?
          OsiSolverInterface *solver2 = saveSolver->clone();
          const double *lower = saveSolver->getColLower();
          const double *upper = saveSolver->getColUpper();
          for (int i = 0; i < numberIntegers_; i++) {
            int iColumn = integerVariable_[i];
            solver2->setColLower(iColumn, lower[iColumn]);
            solver2->setColUpper(iColumn, upper[iColumn]);
          }
          // swap
          delete saveSolver;
          saveSolver = solver2;
          double *newSolution = new double[numberColumns];
          double objectiveValue = checkCutoffForRestart;
          // Save the best solution so far.
          CbcSerendipity heuristic(*this);
          if (bestSolution_)
            heuristic.setInputSolution(bestSolution_, bestObjective_);
          // Magic number
          heuristic.setFractionSmall(0.8);
          // `pumpTune' to stand-alone solver for explanations.
          heuristic.setFeasibilityPumpOptions(1008013);
          // Use numberNodes to say how many are original rows
          heuristic.setNumberNodes(continuousSolver_->getNumRows());
#ifdef COIN_DEVELOP
          if (continuousSolver_->getNumRows() < solver_->getNumRows())
            printf("%d rows added ZZZZZ\n",
              solver_->getNumRows() - continuousSolver_->getNumRows());
#endif
          int returnCode = heuristic.smallBranchAndBound(saveSolver,
            -1, newSolution,
            objectiveValue,
            checkCutoffForRestart, "Reduce");
          if (returnCode < 0) {
#ifdef COIN_DEVELOP
            printf("Restart - not small enough to do search after fixing\n");
#endif
            delete[] newSolution;
          } else {
            // 1 for sol'n, 2 for finished, 3 for both
            if ((returnCode & 1) != 0) {
              // increment number of solutions so other heuristics can test
              numberSolutions_++;
              numberHeuristicSolutions_++;
              lastHeuristic_ = NULL;
              setBestSolution(CBC_ROUNDING, objectiveValue, newSolution);
            }
            delete[] newSolution;
#ifdef CBC_THREAD
            if (master_) {
              lockThread();
              if (parallelMode() > 0) {
                while (master_->waitForThreadsInTree(0)) {
                  lockThread();
                  double dummyBest;
                  tree_->cleanTree(this, -COIN_DBL_MAX, dummyBest);
                  //unlockThread();
                }
              } else {
                double dummyBest;
                tree_->cleanTree(this, -COIN_DBL_MAX, dummyBest);
              }
              master_->waitForThreadsInTree(2);
              delete master_;
              master_ = NULL;
              masterThread_ = NULL;
            }
#endif
            if (tree_->size()) {
              double dummyBest;
              tree_->cleanTree(this, -COIN_DBL_MAX, dummyBest);
            }
            break;
          }
        }
        delete saveSolver;
        saveSolver = NULL;
      }
    }
    /*
          Check for abort on limits: node count, solution count, time, integrality gap.
        */
    if (!(numberNodes_ < intParam_[CbcMaxNumNode] && numberSolutions_ < intParam_[CbcMaxNumSol] && !maximumSecondsReached() && !stoppedOnGap_ && !eventHappened_ && (maximumNumberIterations_ < 0 || numberIterations_ < maximumNumberIterations_))) {
      // out of loop
      break;
    }
#ifdef BONMIN
    assert(!solverCharacteristics_->solutionAddsCuts() || solverCharacteristics_->mipFeasible());
#endif
// Sets percentage of time when we try diving. Diving requires a bit of heap reorganisation, because
// we need to replace the comparison function to dive, and that requires reordering to retain the
// heap property.
#define DIVE_WHEN 1000
#define DIVE_STOP 2000
    int kNode = numberNodes_ % 4000;
    if (numberNodes_ < 100000 && kNode > DIVE_WHEN && kNode <= DIVE_STOP) {
      if (!parallelMode()) {
        if (kNode == DIVE_WHEN + 1 || numberConsecutiveInfeasible > 1) {
          CbcCompareDefault *compare = dynamic_cast< CbcCompareDefault * >(nodeCompare_);
          // Don't interfere if user has replaced the compare function.
          if (compare) {
            //printf("Redoing tree\n");
            compare->startDive(this);
            numberConsecutiveInfeasible = 0;
          }
        }
      }
    }
    // replace current cutoff?
    if (cutoff > getCutoff()) {
      double newCutoff = getCutoff();
      if (analyzeResults_) {
        // see if we could fix any (more)
        int n = 0;
        double *newLower = analyzeResults_;
        double *objLower = newLower + numberIntegers_;
        double *newUpper = objLower + numberIntegers_;
        double *objUpper = newUpper + numberIntegers_;
        for (int i = 0; i < numberIntegers_; i++) {
          if (objLower[i] > newCutoff) {
            n++;
            if (objUpper[i] > newCutoff) {
              newCutoff = -COIN_DBL_MAX;
              break;
            }
            // add as global cut
            objLower[i] = -COIN_DBL_MAX;
            OsiRowCut rc;
            rc.setLb(newLower[i]);
            rc.setUb(COIN_DBL_MAX);
            double one = 1.0;
            rc.setRow(1, integerVariable_ + i, &one, false);
            rc.setGloballyValidAsInteger(2);
            globalCuts_.addCutIfNotDuplicate(rc);
          } else if (objUpper[i] > newCutoff) {
            n++;
            // add as global cut
            objUpper[i] = -COIN_DBL_MAX;
            OsiRowCut rc;
            rc.setLb(-COIN_DBL_MAX);
            rc.setUb(newUpper[i]);
            double one = 1.0;
            rc.setRow(1, integerVariable_ + i, &one, false);
            rc.setGloballyValidAsInteger(2);
            globalCuts_.addCutIfNotDuplicate(rc);
          }
        }
        if (newCutoff == -COIN_DBL_MAX) {
          COIN_DETAIL_PRINT(printf("Root analysis says finished\n"));
        } else if (n > numberFixedNow_) {
          COIN_DETAIL_PRINT(printf("%d more fixed by analysis - now %d\n", n - numberFixedNow_, n));
          numberFixedNow_ = n;
        }
      }
      if (eventHandler) {
        if (!eventHandler->event(CbcEventHandler::solution)) {
          eventHappened_ = true; // exit
        }
        newCutoff = getCutoff();
      }
      lockThread();
      /*
              Clean the tree to reflect the new solution, then see if the
              node comparison predicate wants to make any changes. If so,
              call setComparison for the side effect of rebuilding the heap.
            */
      tree_->cleanTree(this, newCutoff, bestPossibleObjective_);
      if (nodeCompare_->newSolution(this) || nodeCompare_->newSolution(this, continuousObjective_, continuousInfeasibilities_)) {
        tree_->setComparison(*nodeCompare_);
      }
      if (tree_->empty()) {
        continue;
      }
      unlockThread();
    }
    cutoff = getCutoff();
    /*
            Periodic activities: Opportunities to
            + tweak the nodeCompare criteria,
            + check if we've closed the integrality gap enough to quit,
            + print a summary line to let the user know we're working
        */
    if (numberNodes_ >= lastEvery1000) {
      lockThread();
#ifdef COIN_HAS_CLP
      // See if we want dantzig row choice
      goToDantzig(1000, savePivotMethod);
#endif
      lastEvery1000 = numberNodes_ + 1000;
      bool redoTree = nodeCompare_->every1000Nodes(this, numberNodes_);
#ifdef CHECK_CUT_SIZE
      verifyCutSize(tree_, *this);
#endif
      // redo tree if requested
      if (redoTree)
        tree_->setComparison(*nodeCompare_);
      unlockThread();
    }
    // Had hotstart before, now switched off
    if (saveCompare && !hotstartSolution_) {
      // hotstart switched off
      delete nodeCompare_; // off depth first
      nodeCompare_ = saveCompare;
      saveCompare = NULL;
      // redo tree
      lockThread();
      tree_->setComparison(*nodeCompare_);
      unlockThread();
    }
    if (numberNodes_ >= lastPrintEvery) {
      lastPrintEvery = numberNodes_ + printFrequency_;
      lockThread();
      int nNodes = tree_->size();

      //MODIF PIERRE
      bestPossibleObjective_ = tree_->getBestPossibleObjective();
#ifdef CBC_THREAD
      if (parallelMode() > 0 && master_) {
        // need to adjust for ones not on tree
        int numberThreads = master_->numberThreads();
        for (int i = 0; i < numberThreads; i++) {
          CbcThread *child = master_->child(i);
          if (child->node()) {
            // adjust
            double value = child->node()->objectiveValue();
            bestPossibleObjective_ = CoinMin(bestPossibleObjective_, value);
          }
        }
      }
#endif
      unlockThread();
#if CBC_USEFUL_PRINTING > 1
      if (getCutoff() < 1.0e20) {
        if (fabs(getCutoff() - (bestObjective_ - getCutoffIncrement())) > 1.0e-6 && !parentModel_)
          printf("model cutoff in status %g, best %g, increment %g\n",
            getCutoff(), bestObjective_, getCutoffIncrement());
        assert(getCutoff() < bestObjective_ - getCutoffIncrement() + 1.0e-6 + 1.0e-10 * fabs(bestObjective_));
      }
#endif
      if (!intParam_[CbcPrinting]) {
        // Parallel may not have any nodes
        if (!nNodes)
          bestPossibleObjective_ = lastBestPossibleObjective;
        else
          lastBestPossibleObjective = bestPossibleObjective_;
        messageHandler()->message(CBC_STATUS, messages())
          << numberNodes_ << CoinMax(nNodes, 1) << bestObjective_ << bestPossibleObjective_
          << getCurrentSeconds()
          << CoinMessageEol;
      } else if (intParam_[CbcPrinting] == 1) {
        messageHandler()->message(CBC_STATUS2, messages())
          << numberNodes_ << nNodes << bestObjective_ << bestPossibleObjective_
          << tree_->lastDepth() << tree_->lastUnsatisfied()
          << tree_->lastObjective() << numberIterations_
          << getCurrentSeconds()
          << CoinMessageEol;
      } else if (!numberExtraIterations_) {
        messageHandler()->message(CBC_STATUS2, messages())
          << numberNodes_ << nNodes << bestObjective_ << bestPossibleObjective_
          << tree_->lastDepth() << tree_->lastUnsatisfied() << numberIterations_
          << getCurrentSeconds()
          << CoinMessageEol;
      } else {
        messageHandler()->message(CBC_STATUS3, messages())
          << numberNodes_ << numberFathoms_ << numberExtraNodes_ << nNodes
          << bestObjective_ << bestPossibleObjective_
          << tree_->lastDepth() << tree_->lastUnsatisfied() << numberIterations_ << numberExtraIterations_
          << getCurrentSeconds()
          << CoinMessageEol;
      }
#ifdef COIN_HAS_NTY
      if (symmetryInfo_)
        symmetryInfo_->statsOrbits(this, 1);
#endif
#if PRINT_CONFLICT == 1
      if (numberConflictCuts > lastNumberConflictCuts) {
        double length = lengthConflictCuts / numberConflictCuts;
        printf("%d new conflict cuts - total %d - average length %g\n",
          numberConflictCuts - lastNumberConflictCuts,
          numberConflictCuts, length);
        lastNumberConflictCuts = numberConflictCuts;
      }
#endif
      if (eventHandler && !eventHandler->event(CbcEventHandler::treeStatus)) {
        eventHappened_ = true; // exit
      }
    }
    // See if can stop on gap
    if (canStopOnGap()) {
      stoppedOnGap_ = true;
    }

#ifdef CHECK_NODE_FULL
    verifyTreeNodes(tree_, *this);
#endif
#ifdef CHECK_CUT_COUNTS
    verifyCutCounts(tree_, *this);
#endif
    /*
          Now we come to the meat of the loop. To create the active subproblem, we'll
          pop the most promising node in the live set, rebuild the subproblem it
          represents, and then execute the current arm of the branch to create the
          active subproblem.
        */
    CbcNode *node = NULL;
#ifdef CBC_THREAD
    if (!parallelMode() || parallelMode() == -1) {
#endif
      node = tree_->bestNode(cutoff);
      // Possible one on tree worse than cutoff
      // Weird comparison function can leave ineligible nodes on tree
      if (!node || node->objectiveValue() > cutoff)
        continue;
      // Do main work of solving node here
      doOneNode(this, node, createdNode);
#ifdef JJF_ZERO
      if (node) {
        if (createdNode) {
          printf("Node %d depth %d, created %d depth %d\n",
            node->nodeNumber(), node->depth(),
            createdNode->nodeNumber(), createdNode->depth());
        } else {
          printf("Node %d depth %d,  no created node\n",
            node->nodeNumber(), node->depth());
        }
      } else if (createdNode) {
        printf("Node exhausted, created %d depth %d\n",
          createdNode->nodeNumber(), createdNode->depth());
      } else {
        printf("Node exhausted,  no created node\n");
        numberConsecutiveInfeasible = 2;
      }
#endif
      //if (createdNode)
      //numberConsecutiveInfeasible=0;
      //else
      //numberConsecutiveInfeasible++;
#ifdef CBC_THREAD
    } else if (parallelMode() > 0) {
      //lockThread();
      //node = tree_->bestNode(cutoff) ;
      // Possible one on tree worse than cutoff
      if (true || !node || node->objectiveValue() > cutoff) {
        assert(master_);
        if (master_) {
          int anyLeft = master_->waitForThreadsInTree(1);
          // may need to go round again
          if (anyLeft) {
            continue;
          } else {
            master_->stopThreads(-1);
          }
        }
      }
      //unlockThread();
    } else {
      // Deterministic parallel
      if ((tree_->size() < CoinMax(numberThreads_, 8) || hotstartSolution_) && !goneParallel) {
        node = tree_->bestNode(cutoff);
        // Possible one on tree worse than cutoff
        if (!node || node->objectiveValue() > cutoff)
          continue;
        // Do main work of solving node here
        doOneNode(this, node, createdNode);
        assert(createdNode);
        if (!createdNode->active()) {
          delete createdNode;
          createdNode = NULL;
        } else {
          // Say one more pointing to this
          node->nodeInfo()->increment();
          tree_->push(createdNode);
        }
        if (node->active()) {
          assert(node->nodeInfo());
          if (node->nodeInfo()->numberBranchesLeft()) {
            tree_->push(node);
          } else {
            node->setActive(false);
          }
        } else {
          if (node->nodeInfo()) {
            if (!node->nodeInfo()->numberBranchesLeft())
              node->nodeInfo()->allBranchesGone(); // can clean up
            // So will delete underlying stuff
            node->setActive(true);
          }
          delNode[nDeleteNode++] = node;
          node = NULL;
        }
        if (nDeleteNode >= MAX_DEL_NODE) {
          for (int i = 0; i < nDeleteNode; i++) {
            //printf("trying to del %d %x\n",i,delNode[i]);
            delete delNode[i];
            //printf("done to del %d %x\n",i,delNode[i]);
          }
          nDeleteNode = 0;
        }
      } else {
        // Split and solve
        master_->deterministicParallel();
        goneParallel = true;
      }
    }
#endif
  }
  if (nDeleteNode) {
    for (int i = 0; i < nDeleteNode; i++) {
      delete delNode[i];
    }
    nDeleteNode = 0;
  }
#ifdef CBC_THREAD
  if (master_) {
    master_->stopThreads(-1);
    master_->waitForThreadsInTree(2);
    // adjust time to allow for children on some systems
    //dblParam_[CbcStartSeconds] -= CoinCpuTimeJustChildren();
  }
#endif
  /*
      End of the non-abort actions. The next block of code is executed if we've
      aborted because we hit one of the limits. Clean up by deleting the live set
      and break out of the node processing loop. Note that on an abort, node may
      have been pushed back onto the tree for further processing, in which case
      it'll be deleted in cleanTree. We need to check.
    */
  if (!(numberNodes_ < intParam_[CbcMaxNumNode] && numberSolutions_ < intParam_[CbcMaxNumSol] && !maximumSecondsReached() && !stoppedOnGap_ && !eventHappened_ && (maximumNumberIterations_ < 0 || numberIterations_ < maximumNumberIterations_))) {
    if (tree_->size()) {
      double dummyBest;
      tree_->cleanTree(this, -COIN_DBL_MAX, dummyBest);
#if 0 // Does not seem to be needed def CBC_THREAD
	    if (parallelMode() > 0 && master_) {
	      // see if any dangling nodes
	      int numberThreads = master_->numberThreads();
	      for (int i=0;i<numberThreads;i++) {
		CbcThread * child = master_->child(i);
		//if (child->createdNode())
		//printf("CHILD_NODE %p\n",child->createdNode());
		delete child->createdNode();
	      }
	    }
#endif
    }
    delete nextRowCut_;
    /* order is important here:
         * maximumSecondsReached() should be checked before eventHappened_ and
         * isNodeLimitReached() should be checked after eventHappened_
         * reason is, that at timelimit, eventHappened_ is set to true to make Cbc stop fast
         *   and if Ctrl+C is hit, then the nodelimit is set to -1 to make Cbc stop
         */
    if (stoppedOnGap_) {
      messageHandler()->message(CBC_GAP, messages())
        << bestObjective_ - bestPossibleObjective_
        << dblParam_[CbcAllowableGap]
        << dblParam_[CbcAllowableFractionGap] * 100.0
        << CoinMessageEol;
      secondaryStatus_ = 2;
      status_ = 0;
    } else if (maximumSecondsReached()) {
      handler_->message(CBC_MAXTIME, messages_) << CoinMessageEol;
      secondaryStatus_ = 4;
      status_ = 1;
    } else if (numberSolutions_ >= intParam_[CbcMaxNumSol]) {
      handler_->message(CBC_MAXSOLS, messages_) << CoinMessageEol;
      secondaryStatus_ = 6;
      status_ = 1;
    } else if (isNodeLimitReached()) {
      handler_->message(CBC_MAXNODES, messages_) << CoinMessageEol;
      secondaryStatus_ = 3;
      status_ = 1;
    } else if (maximumNumberIterations_ >= 0 && numberIterations_ >= maximumNumberIterations_) {
      handler_->message(CBC_MAXITERS, messages_) << CoinMessageEol;
      secondaryStatus_ = 8;
      status_ = 1;
    } else {
      handler_->message(CBC_EVENT, messages_) << CoinMessageEol;
      secondaryStatus_ = 5;
      status_ = 5;
    }
  }
#ifdef CBC_THREAD
  if (master_) {
    delete master_;
    master_ = NULL;
    masterThread_ = NULL;
  }
#endif
  /*
      That's it, we've exhausted the search tree, or broken out of the loop because
      we hit some limit on evaluation.

      We may have got an intelligent tree so give it one more chance
    */
  // Tell solver we are not in Branch and Cut
  solver_->setHintParam(OsiDoInBranchAndCut, false, OsiHintDo, NULL);
  tree_->endSearch();
  //  If we did any sub trees - did we give up on any?
  if (numberStoppedSubTrees_)
    status_ = 1;
  numberNodes_ += numberExtraNodes_;
  numberIterations_ += numberExtraIterations_;
  if (eventHandler) {
    eventHandler->event(CbcEventHandler::endSearch);
  }
  if (!status_) {
    // Set best possible unless stopped on gap
    if (secondaryStatus_ != 2)
      bestPossibleObjective_ = bestObjective_;
    handler_->message(CBC_END_GOOD, messages_)
      << bestObjective_ << numberIterations_ << numberNodes_ << getCurrentSeconds()
      << CoinMessageEol;
  } else {
    handler_->message(CBC_END, messages_)
      << bestObjective_ << bestPossibleObjective_
      << numberIterations_ << numberNodes_ << getCurrentSeconds()
      << CoinMessageEol;
  }
  if ((moreSpecialOptions_ & 4194304) != 0) {
    // Conflict cuts
    int numberCuts = globalCuts_.sizeRowCuts();
    int nConflict = 0;
    double sizeConflict = 0.0;
    for (int i = 0; i < numberCuts; i++) {
      OsiRowCut2 *cut = globalCuts_.cut(i);
      if (cut->whichRow() == 1) {
        nConflict++;
        sizeConflict += cut->row().getNumElements();
      }
    }
    if (nConflict) {
      sizeConflict /= nConflict;
      char general[200];
      sprintf(general, "%d conflict cuts generated - average length %g",
        nConflict, sizeConflict);
      messageHandler()->message(CBC_GENERAL,
        messages())
        << general << CoinMessageEol;
    }
  }
  if (numberStrongIterations_)
    handler_->message(CBC_STRONG_STATS, messages_)
      << strongInfo_[0] << numberStrongIterations_ << strongInfo_[2]
      << strongInfo_[1] << CoinMessageEol;
  if (!numberExtraNodes_)
    handler_->message(CBC_OTHER_STATS, messages_)
      << maximumDepthActual_
      << numberDJFixed_ << CoinMessageEol;
  else
    handler_->message(CBC_OTHER_STATS2, messages_)
      << maximumDepthActual_
      << numberDJFixed_ << numberFathoms_ << numberExtraNodes_ << numberExtraIterations_
      << CoinMessageEol;
#ifdef COIN_HAS_NTY
  if (symmetryInfo_)
    symmetryInfo_->statsOrbits(this, 1);
#endif
  if (doStatistics == 100) {
    for (int i = 0; i < numberObjects_; i++) {
      CbcSimpleIntegerDynamicPseudoCost *obj = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object_[i]);
      if (obj)
        obj->print();
    }
  }
  if (statistics_) {
    // report in some way
    int *lookup = new int[numberObjects_];
    int i;
    for (i = 0; i < numberObjects_; i++)
      lookup[i] = -1;
    bool goodIds = false; //true;
    for (i = 0; i < numberObjects_; i++) {
      int iColumn = object_[i]->columnNumber();
      if (iColumn >= 0 && iColumn < numberColumns) {
        if (lookup[i] == -1) {
          lookup[i] = iColumn;
        } else {
          goodIds = false;
          break;
        }
      } else {
        goodIds = false;
        break;
      }
    }
    if (!goodIds) {
      delete[] lookup;
      lookup = NULL;
    }
    if (doStatistics >= 3) {
      printf("  node parent depth column   value                    obj      inf\n");
      for (i = 0; i < numberNodes2_; i++) {
        statistics_[i]->print(lookup);
      }
    }
    if (doStatistics > 1) {
      // Find last solution
      int k;
      for (k = numberNodes2_ - 1; k >= 0; k--) {
        if (statistics_[k]->endingObjective() != COIN_DBL_MAX && !statistics_[k]->endingInfeasibility())
          break;
      }
      if (k >= 0) {
        int depth = statistics_[k]->depth();
        int *which = new int[depth + 1];
        for (i = depth; i >= 0; i--) {
          which[i] = k;
          k = statistics_[k]->parentNode();
        }
        printf("  node parent depth column   value                    obj      inf\n");
        for (i = 0; i <= depth; i++) {
          statistics_[which[i]]->print(lookup);
        }
        delete[] which;
      }
    }
    // now summary
    int maxDepth = 0;
    double averageSolutionDepth = 0.0;
    int numberSolutions = 0;
    double averageCutoffDepth = 0.0;
    double averageSolvedDepth = 0.0;
    int numberCutoff = 0;
    int numberDown = 0;
    int numberFirstDown = 0;
    double averageInfDown = 0.0;
    double averageObjDown = 0.0;
    int numberCutoffDown = 0;
    int numberUp = 0;
    int numberFirstUp = 0;
    double averageInfUp = 0.0;
    double averageObjUp = 0.0;
    int numberCutoffUp = 0;
    double averageNumberIterations1 = 0.0;
    double averageValue = 0.0;
    for (i = 0; i < numberNodes2_; i++) {
      int depth = statistics_[i]->depth();
      int way = statistics_[i]->way();
      double value = statistics_[i]->value();
      double startingObjective = statistics_[i]->startingObjective();
      int startingInfeasibility = statistics_[i]->startingInfeasibility();
      double endingObjective = statistics_[i]->endingObjective();
      int endingInfeasibility = statistics_[i]->endingInfeasibility();
      maxDepth = CoinMax(depth, maxDepth);
      // Only for completed
      averageNumberIterations1 += statistics_[i]->numberIterations();
      averageValue += value;
      if (endingObjective != COIN_DBL_MAX && !endingInfeasibility) {
        numberSolutions++;
        averageSolutionDepth += depth;
      }
      if (endingObjective == COIN_DBL_MAX) {
        numberCutoff++;
        averageCutoffDepth += depth;
        if (way < 0) {
          numberDown++;
          numberCutoffDown++;
          if (way == -1)
            numberFirstDown++;
        } else {
          numberUp++;
          numberCutoffUp++;
          if (way == 1)
            numberFirstUp++;
        }
      } else {
        averageSolvedDepth += depth;
        if (way < 0) {
          numberDown++;
          averageInfDown += startingInfeasibility - endingInfeasibility;
          averageObjDown += endingObjective - startingObjective;
          if (way == -1)
            numberFirstDown++;
        } else {
          numberUp++;
          averageInfUp += startingInfeasibility - endingInfeasibility;
          averageObjUp += endingObjective - startingObjective;
          if (way == 1)
            numberFirstUp++;
        }
      }
    }
    // Now print
    if (numberSolutions)
      averageSolutionDepth /= static_cast< double >(numberSolutions);
    int numberSolved = numberNodes2_ - numberCutoff;
    double averageNumberIterations2 = numberIterations_ - averageNumberIterations1
      - numberIterationsAtContinuous;
    if (numberCutoff) {
      averageCutoffDepth /= static_cast< double >(numberCutoff);
      averageNumberIterations2 /= static_cast< double >(numberCutoff);
    }
    if (numberNodes2_)
      averageValue /= static_cast< double >(numberNodes2_);
    if (numberSolved) {
      averageNumberIterations1 /= static_cast< double >(numberSolved);
      averageSolvedDepth /= static_cast< double >(numberSolved);
    }
    printf("%d solution(s) were found (by branching) at an average depth of %g\n",
      numberSolutions, averageSolutionDepth);
    printf("average value of variable being branched on was %g\n",
      averageValue);
    printf("%d nodes were cutoff at an average depth of %g with iteration count of %g\n",
      numberCutoff, averageCutoffDepth, averageNumberIterations2);
    printf("%d nodes were solved at an average depth of %g with iteration count of %g\n",
      numberSolved, averageSolvedDepth, averageNumberIterations1);
    if (numberDown) {
      averageInfDown /= static_cast< double >(numberDown);
      averageObjDown /= static_cast< double >(numberDown);
    }
    printf("Down %d nodes (%d first, %d second) - %d cutoff, rest decrease numinf %g increase obj %g\n",
      numberDown, numberFirstDown, numberDown - numberFirstDown, numberCutoffDown,
      averageInfDown, averageObjDown);
    if (numberUp) {
      averageInfUp /= static_cast< double >(numberUp);
      averageObjUp /= static_cast< double >(numberUp);
    }
    printf("Up %d nodes (%d first, %d second) - %d cutoff, rest decrease numinf %g increase obj %g\n",
      numberUp, numberFirstUp, numberUp - numberFirstUp, numberCutoffUp,
      averageInfUp, averageObjUp);
    for (i = 0; i < numberNodes2_; i++)
      delete statistics_[i];
    delete[] statistics_;
    statistics_ = NULL;
    maximumStatistics_ = 0;
    delete[] lookup;
  }
  /*
      If we think we have a solution, restore and confirm it with a call to
      setBestSolution().  We need to reset the cutoff value so as not to fathom
      the solution on bounds.  Note that calling setBestSolution( ..., true)
      leaves the continuousSolver_ bounds vectors fixed at the solution value.

      Running resolve() here is a failsafe --- setBestSolution has already
      reoptimised using the continuousSolver_. If for some reason we fail to
      prove optimality, run the problem again after instructing the solver to
      tell us more.

      If all looks good, replace solver_ with continuousSolver_, so that the
      outside world will be able to obtain information about the solution using
      public methods.

      Don't replace if we are trying to save cuts
    */
  if (bestSolution_ && (solverCharacteristics_->solverType() < 2 || solverCharacteristics_->solverType() == 4) && ((specialOptions_ & 8388608) == 0 || (specialOptions_ & 2048) != 0)) {
    setCutoff(1.0e50); // As best solution should be worse than cutoff
    // change cutoff as constraint if wanted
    if (cutoffRowNumber_ >= 0) {
      if (solver_->getNumRows() > cutoffRowNumber_)
        solver_->setRowUpper(cutoffRowNumber_, 1.0e50);
    }
    // also in continuousSolver_
    if (continuousSolver_) {
      // Solvers know about direction
      double direction = solver_->getObjSense();
      continuousSolver_->setDblParam(OsiDualObjectiveLimit, 1.0e50 * direction);
    }
    phase_ = 5;
    double increment = getDblParam(CbcModel::CbcCutoffIncrement);
    if ((specialOptions_ & 4) == 0)
      bestObjective_ += 100.0 * increment + 1.0e-3; // only set if we are going to solve
    setBestSolution(CBC_END_SOLUTION, bestObjective_, bestSolution_, 1);
    currentNode_ = NULL;
    continuousSolver_->resolve();
    // Deal with funny variables
    if ((moreSpecialOptions2_ & 32768) != 0)
      cleanBounds(continuousSolver_, NULL);
    if (!continuousSolver_->isProvenOptimal()) {
      continuousSolver_->messageHandler()->setLogLevel(2);
      continuousSolver_->initialSolve();
    }
    delete solver_;
    // above deletes solverCharacteristics_
    solverCharacteristics_ = NULL;
    solver_ = continuousSolver_;
    setPointers(solver_);
    continuousSolver_ = NULL;
  }
  /*
      Clean up dangling objects. continuousSolver_ may already be toast.
    */
  delete lastws;
  if (saveObjects) {
    for (int i = 0; i < numberObjects_; i++)
      delete saveObjects[i];
    delete[] saveObjects;
  }
  numberStrong_ = saveNumberStrong;
  numberBeforeTrust_ = saveNumberBeforeTrust;
  delete[] whichGenerator_;
  whichGenerator_ = NULL;
  delete[] lowerBefore;
  delete[] upperBefore;
  delete[] walkback_;
  walkback_ = NULL;
  delete[] lastNodeInfo_;
  lastNodeInfo_ = NULL;
  delete[] lastNumberCuts_;
  lastNumberCuts_ = NULL;
  delete[] lastCut_;
  lastCut_ = NULL;
  delete[] addedCuts_;
  addedCuts_ = NULL;
  //delete persistentInfo;
  // Get rid of characteristics
  solverCharacteristics_ = NULL;
  if (continuousSolver_) {
    delete continuousSolver_;
    continuousSolver_ = NULL;
  }
  /*
      Destroy global cuts by replacing with an empty OsiCuts object.
    */
  globalCuts_ = CbcRowCuts();
  delete globalConflictCuts_;
  globalConflictCuts_ = NULL;
  if (!bestSolution_ && (specialOptions_ & 8388608) == 0 && false) {
    // make sure lp solver is infeasible
    int numberColumns = solver_->getNumCols();
    const double *columnLower = solver_->getColLower();
    int iColumn;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (solver_->isInteger(iColumn))
        solver_->setColUpper(iColumn, columnLower[iColumn]);
    }
    solver_->initialSolve();
  }
#ifdef COIN_HAS_CLP
  {
    OsiClpSolverInterface *clpSolver
      = dynamic_cast< OsiClpSolverInterface * >(solver_);
    if (clpSolver) {
      // Possible restore of pivot method
      if (savePivotMethod) {
        // model may have changed
        savePivotMethod->setModel(NULL);
        clpSolver->getModelPtr()->setDualRowPivotAlgorithm(*savePivotMethod);
        delete savePivotMethod;
      }
      clpSolver->setLargestAway(-1.0);
    }
  }
#endif
  if ((fastNodeDepth_ >= 1000000 || (moreSpecialOptions_ & 33554432) != 0)
    && !parentModel_) {
    // delete object off end
    delete object_[numberObjects_];
    if ((moreSpecialOptions_ & 33554432) == 0)
      fastNodeDepth_ -= 1000000;
  }
  delete saveSolver;
  // Undo preprocessing performed during BaB.
  if (strategy_ && strategy_->preProcessState() > 0) {
    // undo preprocessing
    CglPreProcess *process = strategy_->process();
    assert(process);
    int n = originalSolver->getNumCols();
    if (bestSolution_) {
      delete[] bestSolution_;
      bestSolution_ = new double[n];
      process->postProcess(*solver_);
    }
    strategy_->deletePreProcess();
    // Solution now back in originalSolver
    delete solver_;
    solver_ = originalSolver;
    if (bestSolution_) {
      bestObjective_ = solver_->getObjValue() * solver_->getObjSense();
      memcpy(bestSolution_, solver_->getColSolution(), n * sizeof(double));
    }
    // put back original objects if there were any
    if (originalObject) {
      int iColumn;
      assert(ownObjects_);
      for (iColumn = 0; iColumn < numberObjects_; iColumn++)
        delete object_[iColumn];
      delete[] object_;
      numberObjects_ = numberOriginalObjects;
      object_ = originalObject;
      delete[] integerVariable_;
      numberIntegers_ = 0;
      for (iColumn = 0; iColumn < n; iColumn++) {
        if (solver_->isInteger(iColumn))
          numberIntegers_++;
      }
      integerVariable_ = new int[numberIntegers_];
      numberIntegers_ = 0;
      for (iColumn = 0; iColumn < n; iColumn++) {
        if (solver_->isInteger(iColumn))
          integerVariable_[numberIntegers_++] = iColumn;
      }
    }
  }
  if (flipObjective)
    flipModel();
#ifdef COIN_HAS_CLP
  {
    OsiClpSolverInterface *clpSolver
      = dynamic_cast< OsiClpSolverInterface * >(solver_);
    if (clpSolver)
      clpSolver->setFakeObjective(reinterpret_cast< double * >(NULL));
  }
#endif
  moreSpecialOptions_ = saveMoreSpecialOptions;
  return;
}

// Solve the initial LP relaxation
void CbcModel::initialSolve()
{
  assert(solver_);
  // Double check optimization directions line up
  dblParam_[CbcOptimizationDirection] = solver_->getObjSense();
  // Check if bounds are all integral (as may get messed up later)
  checkModel();
  if (!solverCharacteristics_) {
    OsiBabSolver *solverCharacteristics = dynamic_cast< OsiBabSolver * >(solver_->getAuxiliaryInfo());
    if (solverCharacteristics) {
      solverCharacteristics_ = solverCharacteristics;
    } else {
      // replace in solver
      OsiBabSolver defaultC;
      solver_->setAuxiliaryInfo(&defaultC);
      solverCharacteristics_ = dynamic_cast< OsiBabSolver * >(solver_->getAuxiliaryInfo());
    }
  }
  solverCharacteristics_->setSolver(solver_);
  solver_->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, NULL);
  // doesn't seem to be uniform time limit
#ifdef COIN_HAS_CLP
  OsiClpSolverInterface *clpSolver
    = dynamic_cast< OsiClpSolverInterface * >(solver_);
  if (clpSolver) {
    double maxTime = dblParam_[CbcMaximumSeconds]-dblParam_[CbcStartSeconds];
    if ((moreSpecialOptions_&131072)==0)
      clpSolver->getModelPtr()->setMaximumSeconds(maxTime);
    else
      clpSolver->getModelPtr()->setMaximumWallSeconds(maxTime);
  }
#endif
  solver_->initialSolve();
  solver_->setHintParam(OsiDoInBranchAndCut, false, OsiHintDo, NULL);
  if (!solver_->isProvenOptimal())
    solver_->resolve();
  // But set up so Jon Lee will be happy
  status_ = -1;
  secondaryStatus_ = -1;
  originalContinuousObjective_ = solver_->getObjValue() * solver_->getObjSense();
  bestPossibleObjective_ = originalContinuousObjective_;
  if (solver_->isProvenDualInfeasible())
    originalContinuousObjective_ = -COIN_DBL_MAX;
  delete[] continuousSolution_;
  continuousSolution_ = CoinCopyOfArray(solver_->getColSolution(),
    solver_->getNumCols());
  setPointers(solver_);
  solverCharacteristics_ = NULL;
}

/*! \brief Get an empty basis object

  Return an empty CoinWarmStartBasis object with the requested capacity,
  appropriate for the current solver. The object is cloned from the object
  cached as emptyWarmStart_. If there is no cached object, the routine
  queries the solver for a warm start object, empties it, and caches the
  result.
*/

CoinWarmStartBasis *CbcModel::getEmptyBasis(int ns, int na) const

{
  CoinWarmStartBasis *emptyBasis;
  /*
      Acquire an empty basis object, if we don't yet have one.
    */
  if (emptyWarmStart_ == 0) {
    if (solver_ == 0) {
      throw CoinError("Cannot construct basis without solver!",
        "getEmptyBasis", "CbcModel");
    }
    emptyBasis = dynamic_cast< CoinWarmStartBasis * >(solver_->getEmptyWarmStart());
    if (emptyBasis == 0) {
      throw CoinError(
        "Solver does not appear to use a basis-oriented warm start.",
        "getEmptyBasis", "CbcModel");
    }
    emptyBasis->setSize(0, 0);
    emptyWarmStart_ = dynamic_cast< CoinWarmStart * >(emptyBasis);
  }
  /*
      Clone the empty basis object, resize it as requested, and return.
    */
  emptyBasis = dynamic_cast< CoinWarmStartBasis * >(emptyWarmStart_->clone());
  assert(emptyBasis);
  if (ns != 0 || na != 0)
    emptyBasis->setSize(ns, na);

  return (emptyBasis);
}

/** Default Constructor

  Creates an empty model without an associated solver.
*/
CbcModel::CbcModel()

  : solver_(NULL)
  , ownership_(0x80000000)
  , continuousSolver_(NULL)
  , referenceSolver_(NULL)
  , defaultHandler_(true)
  , emptyWarmStart_(NULL)
  , bestObjective_(COIN_DBL_MAX)
  , bestPossibleObjective_(COIN_DBL_MAX)
  , sumChangeObjective1_(0.0)
  , sumChangeObjective2_(0.0)
  , bestSolution_(NULL)
  , savedSolutions_(NULL)
  , currentSolution_(NULL)
  , testSolution_(NULL)
  , globalConflictCuts_(NULL)
  , minimumDrop_(1.0e-7)
  , numberSolutions_(0)
  , numberSavedSolutions_(0)
  , maximumSavedSolutions_(0)
  , stateOfSearch_(0)
  , whenCuts_(-1)
  , hotstartSolution_(NULL)
  , hotstartPriorities_(NULL)
  , numberHeuristicSolutions_(0)
  , numberNodes_(0)
  , numberNodes2_(0)
  , numberIterations_(0)
  , numberSolves_(0)
  , status_(-1)
  , secondaryStatus_(-1)
  , numberIntegers_(0)
  , numberRowsAtContinuous_(0)
  , cutoffRowNumber_(-1)
  , maximumNumberCuts_(0)
  , phase_(0)
  , currentNumberCuts_(0)
  , maximumDepth_(0)
  , walkback_(NULL)
  , preProcess_(NULL)
  , lastNodeInfo_(NULL)
  , lastCut_(NULL)
  , lastDepth_(0)
  , lastNumberCuts2_(0)
  , maximumCuts_(0)
  , lastNumberCuts_(NULL)
  , addedCuts_(NULL)
  , nextRowCut_(NULL)
  , currentNode_(NULL)
  , integerVariable_(NULL)
  , integerInfo_(NULL)
  , continuousSolution_(NULL)
  , usedInSolution_(NULL)
  , specialOptions_(0)
  , moreSpecialOptions_(0)
  , moreSpecialOptions2_(0)
  , topOfTree_(NULL)
  , subTreeModel_(NULL)
  , heuristicModel_(NULL)
  , numberStoppedSubTrees_(0)
  , presolve_(0)
  , numberStrong_(5)
  , numberBeforeTrust_(10)
  , numberPenalties_(20)
  , stopNumberIterations_(-1)
  , penaltyScaleFactor_(3.0)
  , numberAnalyzeIterations_(0)
  , analyzeResults_(NULL)
  , numberInfeasibleNodes_(0)
  , problemType_(0)
  , printFrequency_(0)
  , numberCutGenerators_(0)
  , generator_(NULL)
  , virginGenerator_(NULL)
  , numberHeuristics_(0)
  , heuristic_(NULL)
  , lastHeuristic_(NULL)
  , fastNodeDepth_(-1)
  , eventHandler_(NULL)
#ifdef COIN_HAS_NTY
  , symmetryInfo_(NULL)
#endif
  , numberObjects_(0)
  , object_(NULL)
  , ownObjects_(true)
  , originalColumns_(NULL)
  , howOftenGlobalScan_(3)
  , numberGlobalViolations_(0)
  , numberExtraIterations_(0)
  , numberExtraNodes_(0)
  , numberFathoms_(0)
  , continuousObjective_(COIN_DBL_MAX)
  , originalContinuousObjective_(COIN_DBL_MAX)
  , continuousInfeasibilities_(COIN_INT_MAX)
  , maximumCutPassesAtRoot_(20)
  , maximumCutPasses_(10)
  , preferredWay_(0)
  , currentPassNumber_(0)
  , maximumWhich_(INITIAL_MAXIMUM_WHICH)
  , maximumRows_(0)
  , randomSeed_(-1)
  , multipleRootTries_(0)
  , currentDepth_(0)
  , whichGenerator_(NULL)
  , maximumStatistics_(0)
  , statistics_(NULL)
  , maximumDepthActual_(0)
  , numberDJFixed_(0.0)
  , probingInfo_(NULL)
  , numberFixedAtRoot_(0)
  , numberFixedNow_(0)
  , stoppedOnGap_(false)
  , eventHappened_(false)
  , numberLongStrong_(0)
  , numberOldActiveCuts_(0)
  , numberNewCuts_(0)
  , searchStrategy_(-1)
  , strongStrategy_(0)
  , numberStrongIterations_(0)
  , resolveAfterTakeOffCuts_(true)
  , maximumNumberIterations_(-1)
  , continuousPriority_(COIN_INT_MAX)
  , numberUpdateItems_(0)
  , maximumNumberUpdateItems_(0)
  , updateItems_(NULL)
  , storedRowCuts_(NULL)
  , numberThreads_(0)
  , threadMode_(0)
  , numberGlobalCutsIn_(0)
  , master_(NULL)
  , masterThread_(NULL)
{
  memset(intParam_, 0, sizeof(intParam_));
  intParam_[CbcMaxNumNode] = COIN_INT_MAX;
  intParam_[CbcMaxNumSol] = COIN_INT_MAX;

  memset(dblParam_, 0, sizeof(dblParam_));
  dblParam_[CbcIntegerTolerance] = 1e-6;
  dblParam_[CbcCutoffIncrement] = 1e-5;
  dblParam_[CbcAllowableGap] = 1.0e-10;
  dblParam_[CbcMaximumSeconds] = 1.0e100;
  dblParam_[CbcCurrentCutoff] = 1.0e100;
  dblParam_[CbcOptimizationDirection] = 1.0;
  dblParam_[CbcCurrentObjectiveValue] = 1.0e100;
  dblParam_[CbcCurrentMinimizationObjectiveValue] = 1.0e100;
  strongInfo_[0] = 0;
  strongInfo_[1] = 0;
  strongInfo_[2] = 0;
  strongInfo_[3] = 0;
  strongInfo_[4] = 0;
  strongInfo_[5] = 0;
  strongInfo_[6] = 0;
  keepNamesPreproc = false;
  solverCharacteristics_ = NULL;
  nodeCompare_ = new CbcCompareDefault();
  problemFeasibility_ = new CbcFeasibilityBase();
  tree_ = new CbcTree();
  branchingMethod_ = NULL;
  cutModifier_ = NULL;
  strategy_ = NULL;
  parentModel_ = NULL;
  cbcColLower_ = NULL;
  cbcColUpper_ = NULL;
  cbcRowLower_ = NULL;
  cbcRowUpper_ = NULL;
  cbcColSolution_ = NULL;
  cbcRowPrice_ = NULL;
  cbcReducedCost_ = NULL;
  cbcRowActivity_ = NULL;
  appData_ = NULL;
  handler_ = new CoinMessageHandler();
  handler_->setLogLevel(2);
  messages_ = CbcMessage();
  //eventHandler_ = new CbcEventHandler() ;
}

/** Constructor from solver.

  Creates a model complete with a clone of the solver passed as a parameter.
*/

CbcModel::CbcModel(const OsiSolverInterface &rhs)
  : ownership_(0x80000000)
  , continuousSolver_(NULL)
  , referenceSolver_(NULL)
  , defaultHandler_(true)
  , emptyWarmStart_(NULL)
  , bestObjective_(COIN_DBL_MAX)
  , bestPossibleObjective_(COIN_DBL_MAX)
  , sumChangeObjective1_(0.0)
  , sumChangeObjective2_(0.0)
  , globalConflictCuts_(NULL)
  , minimumDrop_(1.0e-7)
  , numberSolutions_(0)
  , numberSavedSolutions_(0)
  , maximumSavedSolutions_(0)
  , stateOfSearch_(0)
  , whenCuts_(-1)
  , hotstartSolution_(NULL)
  , hotstartPriorities_(NULL)
  , numberHeuristicSolutions_(0)
  , numberNodes_(0)
  , numberNodes2_(0)
  , numberIterations_(0)
  , numberSolves_(0)
  , status_(-1)
  , secondaryStatus_(-1)
  , numberRowsAtContinuous_(0)
  , cutoffRowNumber_(-1)
  , maximumNumberCuts_(0)
  , phase_(0)
  , currentNumberCuts_(0)
  , maximumDepth_(0)
  , walkback_(NULL)
  , preProcess_(NULL)
  , lastNodeInfo_(NULL)
  , lastCut_(NULL)
  , lastDepth_(0)
  , lastNumberCuts2_(0)
  , maximumCuts_(0)
  , lastNumberCuts_(NULL)
  , addedCuts_(NULL)
  , nextRowCut_(NULL)
  , currentNode_(NULL)
  , integerInfo_(NULL)
  , specialOptions_(0)
  , moreSpecialOptions_(0)
  , moreSpecialOptions2_(0)
  , topOfTree_(NULL)
  , subTreeModel_(NULL)
  , heuristicModel_(NULL)
  , numberStoppedSubTrees_(0)
  , presolve_(0)
  , numberStrong_(5)
  , numberBeforeTrust_(10)
  , numberPenalties_(20)
  , stopNumberIterations_(-1)
  , penaltyScaleFactor_(3.0)
  , numberAnalyzeIterations_(0)
  , analyzeResults_(NULL)
  , numberInfeasibleNodes_(0)
  , problemType_(0)
  , printFrequency_(0)
  , numberCutGenerators_(0)
  , generator_(NULL)
  , virginGenerator_(NULL)
  , numberHeuristics_(0)
  , heuristic_(NULL)
  , lastHeuristic_(NULL)
  , fastNodeDepth_(-1)
  , eventHandler_(NULL)
#ifdef COIN_HAS_NTY
  , symmetryInfo_(NULL)
#endif
  , numberObjects_(0)
  , object_(NULL)
  , ownObjects_(true)
  , originalColumns_(NULL)
  , howOftenGlobalScan_(3)
  , numberGlobalViolations_(0)
  , numberExtraIterations_(0)
  , numberExtraNodes_(0)
  , numberFathoms_(0)
  , continuousObjective_(COIN_DBL_MAX)
  , originalContinuousObjective_(COIN_DBL_MAX)
  , continuousInfeasibilities_(COIN_INT_MAX)
  , maximumCutPassesAtRoot_(20)
  , maximumCutPasses_(10)
  , preferredWay_(0)
  , currentPassNumber_(0)
  , maximumWhich_(INITIAL_MAXIMUM_WHICH)
  , maximumRows_(0)
  , randomSeed_(-1)
  , multipleRootTries_(0)
  , currentDepth_(0)
  , whichGenerator_(NULL)
  , maximumStatistics_(0)
  , statistics_(NULL)
  , maximumDepthActual_(0)
  , numberDJFixed_(0.0)
  , probingInfo_(NULL)
  , numberFixedAtRoot_(0)
  , numberFixedNow_(0)
  , stoppedOnGap_(false)
  , eventHappened_(false)
  , numberLongStrong_(0)
  , numberOldActiveCuts_(0)
  , numberNewCuts_(0)
  , searchStrategy_(-1)
  , strongStrategy_(0)
  , numberStrongIterations_(0)
  , resolveAfterTakeOffCuts_(true)
  , maximumNumberIterations_(-1)
  , continuousPriority_(COIN_INT_MAX)
  , numberUpdateItems_(0)
  , maximumNumberUpdateItems_(0)
  , updateItems_(NULL)
  , storedRowCuts_(NULL)
  , numberThreads_(0)
  , threadMode_(0)
  , numberGlobalCutsIn_(0)
  , master_(NULL)
  , masterThread_(NULL)
{
  memset(intParam_, 0, sizeof(intParam_));
  intParam_[CbcMaxNumNode] = COIN_INT_MAX;
  intParam_[CbcMaxNumSol] = COIN_INT_MAX;

  memset(dblParam_, 0, sizeof(dblParam_));
  dblParam_[CbcIntegerTolerance] = 1e-6;
  dblParam_[CbcCutoffIncrement] = 1e-5;
  dblParam_[CbcAllowableGap] = 1.0e-10;
  dblParam_[CbcMaximumSeconds] = 1.0e100;
  dblParam_[CbcCurrentCutoff] = 1.0e100;
  dblParam_[CbcOptimizationDirection] = 1.0;
  dblParam_[CbcCurrentObjectiveValue] = 1.0e100;
  dblParam_[CbcCurrentMinimizationObjectiveValue] = 1.0e100;
  strongInfo_[0] = 0;
  strongInfo_[1] = 0;
  strongInfo_[2] = 0;
  strongInfo_[3] = 0;
  strongInfo_[4] = 0;
  strongInfo_[5] = 0;
  strongInfo_[6] = 0;
  solverCharacteristics_ = NULL;
  keepNamesPreproc = false;
  nodeCompare_ = new CbcCompareDefault();
  problemFeasibility_ = new CbcFeasibilityBase();
  tree_ = new CbcTree();
  branchingMethod_ = NULL;
  cutModifier_ = NULL;
  strategy_ = NULL;
  parentModel_ = NULL;
  appData_ = NULL;
  solver_ = rhs.clone();
  ownership_ |= 0x80000000; // model now owns solver
  handler_ = new CoinMessageHandler();
  if (!solver_->defaultHandler() && solver_->messageHandler()->logLevel(0) != -1000)
    passInMessageHandler(solver_->messageHandler());
  handler_->setLogLevel(2);
  messages_ = CbcMessage();
  //eventHandler_ = new CbcEventHandler() ;
  referenceSolver_ = solver_->clone();
  ownership_ = 0x80000000;
  cbcColLower_ = NULL;
  cbcColUpper_ = NULL;
  cbcRowLower_ = NULL;
  cbcRowUpper_ = NULL;
  cbcColSolution_ = NULL;
  cbcRowPrice_ = NULL;
  cbcReducedCost_ = NULL;
  cbcRowActivity_ = NULL;

  // Initialize solution and integer variable vectors
  bestSolution_ = NULL; // to say no solution found
  savedSolutions_ = NULL;
  numberIntegers_ = 0;
  int numberColumns = solver_->getNumCols();
  int iColumn;
  if (numberColumns) {
    // Space for current solution
    currentSolution_ = new double[numberColumns];
    continuousSolution_ = CoinCopyOfArray(solver_->getColSolution(), numberColumns);
    usedInSolution_ = new int[numberColumns];
    CoinZeroN(usedInSolution_, numberColumns);
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (solver_->isInteger(iColumn))
        numberIntegers_++;
    }
  } else {
    // empty model
    currentSolution_ = NULL;
    continuousSolution_ = NULL;
    usedInSolution_ = NULL;
  }
  testSolution_ = currentSolution_;
  if (numberIntegers_) {
    integerVariable_ = new int[numberIntegers_];
    numberIntegers_ = 0;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (solver_->isInteger(iColumn))
        integerVariable_[numberIntegers_++] = iColumn;
    }
  } else {
    integerVariable_ = NULL;
  }
}

static int *resizeInt(int *array, int oldLength, int newLength)
{
  if (!array)
    return NULL;
  assert(newLength > oldLength);
  int *newArray = new int[newLength];
  memcpy(newArray, array, oldLength * sizeof(int));
  delete[] array;
  memset(newArray + oldLength, 0, (newLength - oldLength) * sizeof(int));
  return newArray;
}
static double *resizeDouble(double *array, int oldLength, int newLength)
{
  if (!array)
    return NULL;
  assert(newLength > oldLength);
  double *newArray = new double[newLength];
  memcpy(newArray, array, oldLength * sizeof(double));
  delete[] array;
  memset(newArray + oldLength, 0, (newLength - oldLength) * sizeof(double));
  return newArray;
}
/*
  Assign a solver to the model (model assumes ownership)

  The integer variable vector is initialized if it's not already present.
  If deleteSolver then current solver deleted (if model owned)

  Assuming ownership matches usage in OsiSolverInterface
  (cf. assignProblem, loadProblem).

  TODO: What to do about solver parameters? A simple copy likely won't do it,
	because the SI must push the settings into the underlying solver. In
	the context of switching solvers in cbc, this means that command line
	settings will get lost. Stash the command line somewhere and reread it
	here, maybe?

  TODO: More generally, how much state should be transferred from the old
	solver to the new solver? Best perhaps to see how usage develops.
	What's done here mimics the CbcModel(OsiSolverInterface) constructor.
*/
void CbcModel::assignSolver(OsiSolverInterface *&solver, bool deleteSolver)

{
  // resize stuff if exists
  if (solver && solver_) {
    int nOld = solver_->getNumCols();
    int nNew = solver->getNumCols();
    if (nNew > nOld) {
      originalColumns_ = resizeInt(originalColumns_, nOld, nNew);
      usedInSolution_ = resizeInt(usedInSolution_, nOld, nNew);
      continuousSolution_ = resizeDouble(continuousSolution_, nOld, nNew);
      hotstartSolution_ = resizeDouble(hotstartSolution_, nOld, nNew);
      bestSolution_ = resizeDouble(bestSolution_, nOld, nNew);
      currentSolution_ = resizeDouble(currentSolution_, nOld, nNew);
      if (savedSolutions_) {
        for (int i = 0; i < maximumSavedSolutions_; i++)
          savedSolutions_[i] = resizeDouble(savedSolutions_[i], nOld, nNew);
      }
    }
  }
  // Keep the current message level for solver (if solver exists)
  if (solver_)
    solver->messageHandler()->setLogLevel(solver_->messageHandler()->logLevel());

  if (modelOwnsSolver() && deleteSolver) {
    solverCharacteristics_ = NULL;
    delete solver_;
  }
  solver_ = solver;
  solver = NULL;
  setModelOwnsSolver(true);
  /*
      Basis information is solver-specific.
    */
  if (emptyWarmStart_) {
    delete emptyWarmStart_;
    emptyWarmStart_ = 0;
  }
  bestSolutionBasis_ = CoinWarmStartBasis();
  /*
      Initialize integer variable vector.
    */
  numberIntegers_ = 0;
  int numberColumns = solver_->getNumCols();
  int iColumn;
  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    if (solver_->isInteger(iColumn))
      numberIntegers_++;
  }
  delete[] integerVariable_;
  if (numberIntegers_) {
    integerVariable_ = new int[numberIntegers_];
    numberIntegers_ = 0;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (solver_->isInteger(iColumn))
        integerVariable_[numberIntegers_++] = iColumn;
    }
  } else {
    integerVariable_ = NULL;
  }

  return;
}

// Cloning method

CbcModel *CbcModel::clone(bool cloneHandler)
{
  return new CbcModel(*this, cloneHandler);
}

// Copy constructor.

CbcModel::CbcModel(const CbcModel &rhs, bool cloneHandler)
  : continuousSolver_(NULL)
  , referenceSolver_(NULL)
  , defaultHandler_(rhs.defaultHandler_)
  , emptyWarmStart_(NULL)
  , bestObjective_(rhs.bestObjective_)
  , bestPossibleObjective_(rhs.bestPossibleObjective_)
  , sumChangeObjective1_(rhs.sumChangeObjective1_)
  , sumChangeObjective2_(rhs.sumChangeObjective2_)
  , globalConflictCuts_(NULL)
  , minimumDrop_(rhs.minimumDrop_)
  , numberSolutions_(rhs.numberSolutions_)
  , numberSavedSolutions_(rhs.numberSavedSolutions_)
  , maximumSavedSolutions_(rhs.maximumSavedSolutions_)
  , stateOfSearch_(rhs.stateOfSearch_)
  , whenCuts_(rhs.whenCuts_)
  , numberHeuristicSolutions_(rhs.numberHeuristicSolutions_)
  , numberNodes_(rhs.numberNodes_)
  , numberNodes2_(rhs.numberNodes2_)
  , numberIterations_(rhs.numberIterations_)
  , numberSolves_(rhs.numberSolves_)
  , status_(rhs.status_)
  , secondaryStatus_(rhs.secondaryStatus_)
  , preProcess_(rhs.preProcess_)
  , specialOptions_(rhs.specialOptions_)
  , moreSpecialOptions_(rhs.moreSpecialOptions_)
  , moreSpecialOptions2_(rhs.moreSpecialOptions2_)
  , topOfTree_(NULL)
  , subTreeModel_(rhs.subTreeModel_)
  , heuristicModel_(NULL)
  , numberStoppedSubTrees_(rhs.numberStoppedSubTrees_)
  , presolve_(rhs.presolve_)
  , numberStrong_(rhs.numberStrong_)
  , numberBeforeTrust_(rhs.numberBeforeTrust_)
  , numberPenalties_(rhs.numberPenalties_)
  , stopNumberIterations_(rhs.stopNumberIterations_)
  , penaltyScaleFactor_(rhs.penaltyScaleFactor_)
  , numberAnalyzeIterations_(rhs.numberAnalyzeIterations_)
  , analyzeResults_(NULL)
  , numberInfeasibleNodes_(rhs.numberInfeasibleNodes_)
  , problemType_(rhs.problemType_)
  , printFrequency_(rhs.printFrequency_)
  , fastNodeDepth_(rhs.fastNodeDepth_)
  , howOftenGlobalScan_(rhs.howOftenGlobalScan_)
  , numberGlobalViolations_(rhs.numberGlobalViolations_)
  , numberExtraIterations_(rhs.numberExtraIterations_)
  , numberExtraNodes_(rhs.numberExtraNodes_)
  , numberFathoms_(rhs.numberFathoms_)
  , continuousObjective_(rhs.continuousObjective_)
  , originalContinuousObjective_(rhs.originalContinuousObjective_)
  , continuousInfeasibilities_(rhs.continuousInfeasibilities_)
  , maximumCutPassesAtRoot_(rhs.maximumCutPassesAtRoot_)
  , maximumCutPasses_(rhs.maximumCutPasses_)
  , preferredWay_(rhs.preferredWay_)
  , currentPassNumber_(rhs.currentPassNumber_)
  , maximumWhich_(rhs.maximumWhich_)
  , maximumRows_(0)
  , randomSeed_(rhs.randomSeed_)
  , multipleRootTries_(rhs.multipleRootTries_)
  , currentDepth_(0)
  , whichGenerator_(NULL)
  , maximumStatistics_(0)
  , statistics_(NULL)
  , maximumDepthActual_(0)
  , numberDJFixed_(0.0)
  , probingInfo_(NULL)
  , numberFixedAtRoot_(rhs.numberFixedAtRoot_)
  , numberFixedNow_(rhs.numberFixedNow_)
  , stoppedOnGap_(rhs.stoppedOnGap_)
  , eventHappened_(rhs.eventHappened_)
  , numberLongStrong_(rhs.numberLongStrong_)
  , numberOldActiveCuts_(rhs.numberOldActiveCuts_)
  , numberNewCuts_(rhs.numberNewCuts_)
  , searchStrategy_(rhs.searchStrategy_)
  , strongStrategy_(rhs.strongStrategy_)
  , numberStrongIterations_(rhs.numberStrongIterations_)
  , resolveAfterTakeOffCuts_(rhs.resolveAfterTakeOffCuts_)
  , maximumNumberIterations_(rhs.maximumNumberIterations_)
  , continuousPriority_(rhs.continuousPriority_)
  , numberUpdateItems_(rhs.numberUpdateItems_)
  , maximumNumberUpdateItems_(rhs.maximumNumberUpdateItems_)
  , updateItems_(NULL)
  , storedRowCuts_(NULL)
  , numberThreads_(rhs.numberThreads_)
  , threadMode_(rhs.threadMode_)
  , numberGlobalCutsIn_(rhs.numberGlobalCutsIn_)
  , master_(NULL)
  , masterThread_(NULL)
{
  memcpy(intParam_, rhs.intParam_, sizeof(intParam_));
  memcpy(dblParam_, rhs.dblParam_, sizeof(dblParam_));
  strongInfo_[0] = rhs.strongInfo_[0];
  strongInfo_[1] = rhs.strongInfo_[1];
  strongInfo_[2] = rhs.strongInfo_[2];
  strongInfo_[3] = rhs.strongInfo_[3];
  strongInfo_[4] = rhs.strongInfo_[4];
  strongInfo_[5] = rhs.strongInfo_[5];
  strongInfo_[6] = rhs.strongInfo_[6];
  keepNamesPreproc = rhs.keepNamesPreproc;
  solverCharacteristics_ = NULL;
  if (rhs.emptyWarmStart_)
    emptyWarmStart_ = rhs.emptyWarmStart_->clone();
  if (defaultHandler_ || cloneHandler) {
    handler_ = new CoinMessageHandler();
    handler_->setLogLevel(2);
  } else {
    handler_ = rhs.handler_;
  }
  messageHandler()->setLogLevel(rhs.messageHandler()->logLevel());
  numberCutGenerators_ = rhs.numberCutGenerators_;
  if (numberCutGenerators_) {
    generator_ = new CbcCutGenerator *[numberCutGenerators_];
    virginGenerator_ = new CbcCutGenerator *[numberCutGenerators_];
    int i;
    for (i = 0; i < numberCutGenerators_; i++) {
      generator_[i] = new CbcCutGenerator(*rhs.generator_[i]);
      virginGenerator_[i] = new CbcCutGenerator(*rhs.virginGenerator_[i]);
    }
  } else {
    generator_ = NULL;
    virginGenerator_ = NULL;
  }
  globalCuts_ = rhs.globalCuts_;
  numberHeuristics_ = rhs.numberHeuristics_;
  if (numberHeuristics_) {
    heuristic_ = new CbcHeuristic *[numberHeuristics_];
    int i;
    for (i = 0; i < numberHeuristics_; i++) {
      heuristic_[i] = rhs.heuristic_[i]->clone();
    }
  } else {
    heuristic_ = NULL;
  }
  lastHeuristic_ = NULL;
  if (rhs.eventHandler_) {
    eventHandler_ = rhs.eventHandler_->clone();
  } else {
    eventHandler_ = NULL;
  }
  ownObjects_ = rhs.ownObjects_;
  if (ownObjects_) {
    numberObjects_ = rhs.numberObjects_;
    if (numberObjects_) {
      object_ = new OsiObject *[numberObjects_];
      int i;
      for (i = 0; i < numberObjects_; i++) {
        object_[i] = (rhs.object_[i])->clone();
        CbcObject *obj = dynamic_cast< CbcObject * >(object_[i]);
        // Could be OsiObjects
        if (obj)
          obj->setModel(this);
      }
    } else {
      object_ = NULL;
    }
  } else {
    // assume will be redone
    numberObjects_ = 0;
    object_ = NULL;
  }
  if (rhs.continuousSolver_) {
    continuousSolver_ = rhs.continuousSolver_->clone();
  } else {
    continuousSolver_ = NULL;
  }
  if (rhs.referenceSolver_)
    referenceSolver_ = rhs.referenceSolver_->clone();
  else
    referenceSolver_ = NULL;
  solver_ = rhs.solver_->clone();
  if (rhs.originalColumns_) {
    int numberColumns = solver_->getNumCols();
    originalColumns_ = new int[numberColumns];
    memcpy(originalColumns_, rhs.originalColumns_, numberColumns * sizeof(int));
  } else {
    originalColumns_ = NULL;
  }
  if (maximumNumberUpdateItems_) {
    updateItems_ = new CbcObjectUpdateData[maximumNumberUpdateItems_];
    for (int i = 0; i < maximumNumberUpdateItems_; i++)
      updateItems_[i] = rhs.updateItems_[i];
  }
  if (maximumWhich_ && rhs.whichGenerator_)
    whichGenerator_ = CoinCopyOfArray(rhs.whichGenerator_, maximumWhich_);
  nodeCompare_ = rhs.nodeCompare_->clone();
  problemFeasibility_ = rhs.problemFeasibility_->clone();
  tree_ = rhs.tree_->clone();
  if (rhs.branchingMethod_)
    branchingMethod_ = rhs.branchingMethod_->clone();
  else
    branchingMethod_ = NULL;
  if (rhs.cutModifier_)
    cutModifier_ = rhs.cutModifier_->clone();
  else
    cutModifier_ = NULL;
  cbcColLower_ = NULL;
  cbcColUpper_ = NULL;
  cbcRowLower_ = NULL;
  cbcRowUpper_ = NULL;
  cbcColSolution_ = NULL;
  cbcRowPrice_ = NULL;
  cbcReducedCost_ = NULL;
  cbcRowActivity_ = NULL;
  if (rhs.strategy_)
    strategy_ = rhs.strategy_->clone();
  else
    strategy_ = NULL;
  parentModel_ = rhs.parentModel_;
  appData_ = rhs.appData_;
  messages_ = rhs.messages_;
  ownership_ = rhs.ownership_ | 0x80000000;
  messageHandler()->setLogLevel(rhs.messageHandler()->logLevel());
  numberIntegers_ = rhs.numberIntegers_;
  randomNumberGenerator_ = rhs.randomNumberGenerator_;
  if (numberIntegers_) {
    integerVariable_ = new int[numberIntegers_];
    memcpy(integerVariable_, rhs.integerVariable_, numberIntegers_ * sizeof(int));
    integerInfo_ = CoinCopyOfArray(rhs.integerInfo_, solver_->getNumCols());
  } else {
    integerVariable_ = NULL;
    integerInfo_ = NULL;
  }
  if (rhs.hotstartSolution_) {
    int numberColumns = solver_->getNumCols();
    hotstartSolution_ = CoinCopyOfArray(rhs.hotstartSolution_, numberColumns);
    hotstartPriorities_ = CoinCopyOfArray(rhs.hotstartPriorities_, numberColumns);
  } else {
    hotstartSolution_ = NULL;
    hotstartPriorities_ = NULL;
  }
  if (rhs.bestSolution_) {
    int numberColumns = solver_->getNumCols();
    bestSolution_ = new double[numberColumns];
    memcpy(bestSolution_, rhs.bestSolution_, numberColumns * sizeof(double));
  } else {
    bestSolution_ = NULL;
  }
  int numberColumns = solver_->getNumCols();
  if (maximumSavedSolutions_ && rhs.savedSolutions_) {
    savedSolutions_ = new double *[maximumSavedSolutions_];
    for (int i = 0; i < maximumSavedSolutions_; i++)
      savedSolutions_[i] = CoinCopyOfArray(rhs.savedSolutions_[i], numberColumns + 2);
  } else {
    savedSolutions_ = NULL;
  }
  // Space for current solution
  if (numberColumns) {
    currentSolution_ = new double[numberColumns];
    continuousSolution_ = CoinCopyOfArray(solver_->getColSolution(), numberColumns);
    usedInSolution_ = new int[numberColumns];
    CoinZeroN(usedInSolution_, numberColumns);
  } else {
    currentSolution_ = NULL;
    continuousSolution_ = NULL;
    usedInSolution_ = NULL;
  }
  testSolution_ = currentSolution_;
  numberRowsAtContinuous_ = rhs.numberRowsAtContinuous_;
  cutoffRowNumber_ = rhs.cutoffRowNumber_;
  maximumNumberCuts_ = rhs.maximumNumberCuts_;
  phase_ = rhs.phase_;
  currentNumberCuts_ = rhs.currentNumberCuts_;
  maximumDepth_ = rhs.maximumDepth_;
  // These are only used as temporary arrays so need not be filled
  if (maximumNumberCuts_) {
    addedCuts_ = new CbcCountRowCut *[maximumNumberCuts_];
  } else {
    addedCuts_ = NULL;
  }
  bestSolutionBasis_ = rhs.bestSolutionBasis_;
  nextRowCut_ = NULL;
  currentNode_ = NULL;
  if (maximumDepth_) {
    walkback_ = new CbcNodeInfo *[maximumDepth_];
    lastNodeInfo_ = new CbcNodeInfo *[maximumDepth_];
    lastNumberCuts_ = new int[maximumDepth_];
  } else {
    walkback_ = NULL;
    lastNodeInfo_ = NULL;
    lastNumberCuts_ = NULL;
  }
  maximumCuts_ = rhs.maximumCuts_;
  if (maximumCuts_) {
    lastCut_ = new const OsiRowCut *[maximumCuts_];
  } else {
    lastCut_ = NULL;
  }
#ifdef COIN_HAS_NTY
  if (rhs.symmetryInfo_)
    symmetryInfo_ = new CbcSymmetry(*rhs.symmetryInfo_);
  else
    symmetryInfo_ = NULL;
#endif
  synchronizeModel();
  if (cloneHandler && !defaultHandler_) {
    delete handler_;
    /* We have to clone handlers - otherwise will all be
	   writing to same buffer.  So if threads user will
	   have to sychronize */
    CoinMessageHandler *handler = rhs.handler_->clone();
    passInMessageHandler(handler);
    defaultHandler_ = true;
  }
}

// Assignment operator
CbcModel &
CbcModel::operator=(const CbcModel &rhs)
{
  if (this != &rhs) {
    if (modelOwnsSolver()) {
      solverCharacteristics_ = NULL;
      delete solver_;
      solver_ = NULL;
    }
    gutsOfDestructor();
    if (defaultHandler_) {
      delete handler_;
      handler_ = NULL;
    }
    defaultHandler_ = rhs.defaultHandler_;
    if (defaultHandler_) {
      handler_ = new CoinMessageHandler();
      handler_->setLogLevel(2);
    } else {
      handler_ = rhs.handler_;
    }
    messages_ = rhs.messages_;
    messageHandler()->setLogLevel(rhs.messageHandler()->logLevel());
    if (rhs.solver_) {
      solver_ = rhs.solver_->clone();
    } else {
      solver_ = 0;
    }
    ownership_ = 0x80000000;
    delete continuousSolver_;
    if (rhs.continuousSolver_) {
      continuousSolver_ = rhs.continuousSolver_->clone();
    } else {
      continuousSolver_ = 0;
    }
    delete referenceSolver_;
    if (rhs.referenceSolver_) {
      referenceSolver_ = rhs.referenceSolver_->clone();
    } else {
      referenceSolver_ = NULL;
    }

    delete emptyWarmStart_;
    if (rhs.emptyWarmStart_) {
      emptyWarmStart_ = rhs.emptyWarmStart_->clone();
    } else {
      emptyWarmStart_ = 0;
    }

    bestObjective_ = rhs.bestObjective_;
    bestPossibleObjective_ = rhs.bestPossibleObjective_;
    sumChangeObjective1_ = rhs.sumChangeObjective1_;
    sumChangeObjective2_ = rhs.sumChangeObjective2_;
    delete[] bestSolution_;
    if (rhs.bestSolution_) {
      int numberColumns = rhs.getNumCols();
      bestSolution_ = new double[numberColumns];
      memcpy(bestSolution_, rhs.bestSolution_, numberColumns * sizeof(double));
    } else {
      bestSolution_ = NULL;
    }
    for (int i = 0; i < maximumSavedSolutions_; i++)
      delete[] savedSolutions_[i];
    delete[] savedSolutions_;
    savedSolutions_ = NULL;
    int numberColumns = rhs.getNumCols();
    if (numberColumns) {
      // Space for current solution
      currentSolution_ = new double[numberColumns];
      continuousSolution_ = CoinCopyOfArray(solver_->getColSolution(), numberColumns);
      usedInSolution_ = new int[numberColumns];
      CoinZeroN(usedInSolution_, numberColumns);
    } else {
      currentSolution_ = NULL;
      continuousSolution_ = NULL;
      usedInSolution_ = NULL;
    }
    if (maximumSavedSolutions_) {
      savedSolutions_ = new double *[maximumSavedSolutions_];
      for (int i = 0; i < maximumSavedSolutions_; i++)
        savedSolutions_[i] = CoinCopyOfArray(rhs.savedSolutions_[i], numberColumns + 2);
    } else {
      savedSolutions_ = NULL;
    }
    testSolution_ = currentSolution_;
    minimumDrop_ = rhs.minimumDrop_;
    numberSolutions_ = rhs.numberSolutions_;
    numberSavedSolutions_ = rhs.numberSavedSolutions_;
    maximumSavedSolutions_ = rhs.maximumSavedSolutions_;
    stateOfSearch_ = rhs.stateOfSearch_;
    whenCuts_ = rhs.whenCuts_;
    numberHeuristicSolutions_ = rhs.numberHeuristicSolutions_;
    numberNodes_ = rhs.numberNodes_;
    numberNodes2_ = rhs.numberNodes2_;
    numberIterations_ = rhs.numberIterations_;
    numberSolves_ = rhs.numberSolves_;
    status_ = rhs.status_;
    secondaryStatus_ = rhs.secondaryStatus_;
    specialOptions_ = rhs.specialOptions_;
    moreSpecialOptions_ = rhs.moreSpecialOptions_;
    moreSpecialOptions2_ = rhs.moreSpecialOptions2_;
    subTreeModel_ = rhs.subTreeModel_;
    heuristicModel_ = NULL;
    numberStoppedSubTrees_ = rhs.numberStoppedSubTrees_;
    presolve_ = rhs.presolve_;
    numberStrong_ = rhs.numberStrong_;
    numberBeforeTrust_ = rhs.numberBeforeTrust_;
    numberPenalties_ = rhs.numberPenalties_;
    stopNumberIterations_ = rhs.stopNumberIterations_;
    penaltyScaleFactor_ = rhs.penaltyScaleFactor_;
    numberAnalyzeIterations_ = rhs.numberAnalyzeIterations_;
    delete[] analyzeResults_;
    analyzeResults_ = NULL;
    numberInfeasibleNodes_ = rhs.numberInfeasibleNodes_;
    problemType_ = rhs.problemType_;
    printFrequency_ = rhs.printFrequency_;
    howOftenGlobalScan_ = rhs.howOftenGlobalScan_;
    numberGlobalViolations_ = rhs.numberGlobalViolations_;
    numberExtraIterations_ = rhs.numberExtraIterations_;
    numberExtraNodes_ = rhs.numberExtraNodes_;
    preProcess_ = rhs.preProcess_;
    numberFathoms_ = rhs.numberFathoms_;
    continuousObjective_ = rhs.continuousObjective_;
    originalContinuousObjective_ = rhs.originalContinuousObjective_;
    continuousInfeasibilities_ = rhs.continuousInfeasibilities_;
    maximumCutPassesAtRoot_ = rhs.maximumCutPassesAtRoot_;
    maximumCutPasses_ = rhs.maximumCutPasses_;
    randomSeed_ = rhs.randomSeed_;
    multipleRootTries_ = rhs.multipleRootTries_;
    preferredWay_ = rhs.preferredWay_;
    currentPassNumber_ = rhs.currentPassNumber_;
    memcpy(intParam_, rhs.intParam_, sizeof(intParam_));
    memcpy(dblParam_, rhs.dblParam_, sizeof(dblParam_));
    globalCuts_ = rhs.globalCuts_;
    delete globalConflictCuts_;
    globalConflictCuts_ = NULL;
    int i;
    for (i = 0; i < numberCutGenerators_; i++) {
      delete generator_[i];
      delete virginGenerator_[i];
    }
    delete[] generator_;
    delete[] virginGenerator_;
    delete[] heuristic_;
    maximumWhich_ = rhs.maximumWhich_;
    delete[] whichGenerator_;
    whichGenerator_ = NULL;
    if (maximumWhich_ && rhs.whichGenerator_)
      whichGenerator_ = CoinCopyOfArray(rhs.whichGenerator_, maximumWhich_);
    maximumRows_ = 0;
    currentDepth_ = 0;
    randomNumberGenerator_ = rhs.randomNumberGenerator_;
    workingBasis_ = CoinWarmStartBasis();
    for (i = 0; i < maximumStatistics_; i++)
      delete statistics_[i];
    delete[] statistics_;
    maximumStatistics_ = 0;
    statistics_ = NULL;
    delete probingInfo_;
    probingInfo_ = NULL;
    numberFixedAtRoot_ = rhs.numberFixedAtRoot_;
    numberFixedNow_ = rhs.numberFixedNow_;
    stoppedOnGap_ = rhs.stoppedOnGap_;
    eventHappened_ = rhs.eventHappened_;
    numberLongStrong_ = rhs.numberLongStrong_;
    numberOldActiveCuts_ = rhs.numberOldActiveCuts_;
    numberNewCuts_ = rhs.numberNewCuts_;
    resolveAfterTakeOffCuts_ = rhs.resolveAfterTakeOffCuts_;
    maximumNumberIterations_ = rhs.maximumNumberIterations_;
    continuousPriority_ = rhs.continuousPriority_;
    numberUpdateItems_ = rhs.numberUpdateItems_;
    maximumNumberUpdateItems_ = rhs.maximumNumberUpdateItems_;
    delete[] updateItems_;
    if (maximumNumberUpdateItems_) {
      updateItems_ = new CbcObjectUpdateData[maximumNumberUpdateItems_];
      for (i = 0; i < maximumNumberUpdateItems_; i++)
        updateItems_[i] = rhs.updateItems_[i];
    } else {
      updateItems_ = NULL;
    }
    numberThreads_ = rhs.numberThreads_;
    threadMode_ = rhs.threadMode_;
    numberGlobalCutsIn_ = rhs.numberGlobalCutsIn_;
    delete master_;
    master_ = NULL;
    masterThread_ = NULL;
    searchStrategy_ = rhs.searchStrategy_;
    strongStrategy_ = rhs.strongStrategy_;
    numberStrongIterations_ = rhs.numberStrongIterations_;
    strongInfo_[0] = rhs.strongInfo_[0];
    strongInfo_[1] = rhs.strongInfo_[1];
    strongInfo_[2] = rhs.strongInfo_[2];
    strongInfo_[3] = rhs.strongInfo_[3];
    strongInfo_[4] = rhs.strongInfo_[4];
    strongInfo_[5] = rhs.strongInfo_[5];
    strongInfo_[6] = rhs.strongInfo_[6];
    solverCharacteristics_ = NULL;
    lastHeuristic_ = NULL;
    numberCutGenerators_ = rhs.numberCutGenerators_;
    if (numberCutGenerators_) {
      generator_ = new CbcCutGenerator *[numberCutGenerators_];
      virginGenerator_ = new CbcCutGenerator *[numberCutGenerators_];
      int i;
      for (i = 0; i < numberCutGenerators_; i++) {
        generator_[i] = new CbcCutGenerator(*rhs.generator_[i]);
        virginGenerator_[i] = new CbcCutGenerator(*rhs.virginGenerator_[i]);
      }
    } else {
      generator_ = NULL;
      virginGenerator_ = NULL;
    }
    numberHeuristics_ = rhs.numberHeuristics_;
    if (numberHeuristics_) {
      heuristic_ = new CbcHeuristic *[numberHeuristics_];
      memcpy(heuristic_, rhs.heuristic_,
        numberHeuristics_ * sizeof(CbcHeuristic *));
    } else {
      heuristic_ = NULL;
    }
    lastHeuristic_ = NULL;
    if (eventHandler_)
      delete eventHandler_;
    if (rhs.eventHandler_) {
      eventHandler_ = rhs.eventHandler_->clone();
    } else {
      eventHandler_ = NULL;
    }
    fastNodeDepth_ = rhs.fastNodeDepth_;
    if (ownObjects_) {
      for (i = 0; i < numberObjects_; i++)
        delete object_[i];
      delete[] object_;
      numberObjects_ = rhs.numberObjects_;
      if (numberObjects_) {
        object_ = new OsiObject *[numberObjects_];
        int i;
        for (i = 0; i < numberObjects_; i++)
          object_[i] = (rhs.object_[i])->clone();
      } else {
        object_ = NULL;
      }
    } else {
      // assume will be redone
      numberObjects_ = 0;
      object_ = NULL;
    }
    delete[] originalColumns_;
    if (rhs.originalColumns_) {
      int numberColumns = rhs.getNumCols();
      originalColumns_ = new int[numberColumns];
      memcpy(originalColumns_, rhs.originalColumns_, numberColumns * sizeof(int));
    } else {
      originalColumns_ = NULL;
    }
    nodeCompare_ = rhs.nodeCompare_->clone();
    problemFeasibility_ = rhs.problemFeasibility_->clone();
    delete tree_;
    tree_ = rhs.tree_->clone();
    if (rhs.branchingMethod_)
      branchingMethod_ = rhs.branchingMethod_->clone();
    else
      branchingMethod_ = NULL;
    if (rhs.cutModifier_)
      cutModifier_ = rhs.cutModifier_->clone();
    else
      cutModifier_ = NULL;

    if (strategy_)
        delete strategy_;
    if (rhs.strategy_)
      strategy_ = rhs.strategy_->clone();
    else
      strategy_ = NULL;
    parentModel_ = rhs.parentModel_;
    appData_ = rhs.appData_;

    delete[] integerVariable_;
    numberIntegers_ = rhs.numberIntegers_;
    if (numberIntegers_) {
      integerVariable_ = new int[numberIntegers_];
      memcpy(integerVariable_, rhs.integerVariable_,
        numberIntegers_ * sizeof(int));
      integerInfo_ = CoinCopyOfArray(rhs.integerInfo_, rhs.getNumCols());
    } else {
      integerVariable_ = NULL;
      integerInfo_ = NULL;
    }
    if (rhs.hotstartSolution_) {
      int numberColumns = solver_->getNumCols();
      hotstartSolution_ = CoinCopyOfArray(rhs.hotstartSolution_, numberColumns);
      hotstartPriorities_ = CoinCopyOfArray(rhs.hotstartPriorities_, numberColumns);
    } else {
      hotstartSolution_ = NULL;
      hotstartPriorities_ = NULL;
    }
    numberRowsAtContinuous_ = rhs.numberRowsAtContinuous_;
    cutoffRowNumber_ = rhs.cutoffRowNumber_;
    maximumNumberCuts_ = rhs.maximumNumberCuts_;
    phase_ = rhs.phase_;
    currentNumberCuts_ = rhs.currentNumberCuts_;
    maximumDepth_ = rhs.maximumDepth_;
    keepNamesPreproc = rhs.keepNamesPreproc;
    mipStart_ = rhs.mipStart_;
    delete[] addedCuts_;
    delete[] walkback_;
    // These are only used as temporary arrays so need not be filled
    if (maximumNumberCuts_) {
      addedCuts_ = new CbcCountRowCut *[maximumNumberCuts_];
    } else {
      addedCuts_ = NULL;
    }
    delete[] lastNodeInfo_;
    delete[] lastNumberCuts_;
    delete[] lastCut_;
    bestSolutionBasis_ = rhs.bestSolutionBasis_;
    nextRowCut_ = NULL;
    currentNode_ = NULL;
    if (maximumDepth_) {
      walkback_ = new CbcNodeInfo *[maximumDepth_];
      lastNodeInfo_ = new CbcNodeInfo *[maximumDepth_];
      lastNumberCuts_ = new int[maximumDepth_];
    } else {
      walkback_ = NULL;
      lastNodeInfo_ = NULL;
      lastNumberCuts_ = NULL;
    }
    maximumCuts_ = rhs.maximumCuts_;
    if (maximumCuts_) {
      lastCut_ = new const OsiRowCut *[maximumCuts_];
    } else {
      lastCut_ = NULL;
    }
#ifdef COIN_HAS_NTY
    if (rhs.symmetryInfo_)
      symmetryInfo_ = new CbcSymmetry(*rhs.symmetryInfo_);
    else
      symmetryInfo_ = NULL;
#endif
    synchronizeModel();
    cbcColLower_ = NULL;
    cbcColUpper_ = NULL;
    cbcRowLower_ = NULL;
    cbcRowUpper_ = NULL;
    cbcColSolution_ = NULL;
    cbcRowPrice_ = NULL;
    cbcReducedCost_ = NULL;
    cbcRowActivity_ = NULL;
  }
  return *this;
}
// Destructor
CbcModel::~CbcModel()
{
  if (defaultHandler_) {
    delete handler_;
    handler_ = NULL;
  }
  delete tree_;
  tree_ = NULL;
  if (modelOwnsSolver()) {
    delete solver_;
    solver_ = NULL;
  }
  gutsOfDestructor();
  delete eventHandler_;
  eventHandler_ = NULL;
#ifdef CBC_THREAD
  // Get rid of all threaded stuff
  delete master_;
#endif
}
// Clears out as much as possible (except solver)
void CbcModel::gutsOfDestructor()
{
  delete referenceSolver_;
  referenceSolver_ = NULL;
  int i;
  for (i = 0; i < numberCutGenerators_; i++) {
    delete generator_[i];
    delete virginGenerator_[i];
  }
  delete[] generator_;
  delete[] virginGenerator_;
  generator_ = NULL;
  virginGenerator_ = NULL;
  for (i = 0; i < numberHeuristics_; i++)
    delete heuristic_[i];
  delete[] heuristic_;
  heuristic_ = NULL;
  delete nodeCompare_;
  nodeCompare_ = NULL;
  delete problemFeasibility_;
  problemFeasibility_ = NULL;
  delete[] originalColumns_;
  originalColumns_ = NULL;
  delete strategy_;
  delete[] updateItems_;
  updateItems_ = NULL;
  numberUpdateItems_ = 0;
  maximumNumberUpdateItems_ = 0;
  gutsOfDestructor2();
}
// Clears out enough to reset CbcModel
void CbcModel::gutsOfDestructor2()
{
  delete[] integerInfo_;
  integerInfo_ = NULL;
  delete[] integerVariable_;
  integerVariable_ = NULL;
  int i;
  if (ownObjects_) {
    for (i = 0; i < numberObjects_; i++)
      delete object_[i];
    delete[] object_;
  }
  ownObjects_ = true;
  object_ = NULL;
  numberIntegers_ = 0;
  numberObjects_ = 0;
  // Below here is whatever consensus is
  ownership_ = 0x80000000;
  delete branchingMethod_;
  branchingMethod_ = NULL;
  delete cutModifier_;
  cutModifier_ = NULL;
  topOfTree_ = NULL;
  resetModel();
#ifdef COIN_HAS_NTY
  delete symmetryInfo_;
  symmetryInfo_ = NULL;
#endif
}
// Clears out enough to reset CbcModel
void CbcModel::resetModel()
{
  delete emptyWarmStart_;
  emptyWarmStart_ = NULL;
  delete continuousSolver_;
  continuousSolver_ = NULL;
  numberSavedSolutions_ = 0;
  delete[] bestSolution_;
  bestSolution_ = NULL;
  if (savedSolutions_) {
    for (int i = 0; i < maximumSavedSolutions_; i++)
      delete[] savedSolutions_[i];
    delete[] savedSolutions_;
    savedSolutions_ = NULL;
  }
  delete[] currentSolution_;
  currentSolution_ = NULL;
  delete[] continuousSolution_;
  continuousSolution_ = NULL;
  solverCharacteristics_ = NULL;
  delete[] usedInSolution_;
  usedInSolution_ = NULL;
  testSolution_ = NULL;
  lastHeuristic_ = NULL;
  delete[] addedCuts_;
  addedCuts_ = NULL;
  nextRowCut_ = NULL;
  currentNode_ = NULL;
  delete[] walkback_;
  walkback_ = NULL;
  delete[] lastNodeInfo_;
  lastNodeInfo_ = NULL;
  delete[] lastNumberCuts_;
  lastNumberCuts_ = NULL;
  delete[] lastCut_;
  lastCut_ = NULL;
  delete[] whichGenerator_;
  whichGenerator_ = NULL;
  for (int i = 0; i < maximumStatistics_; i++)
    delete statistics_[i];
  delete[] statistics_;
  statistics_ = NULL;
  maximumDepthActual_ = 0;
  numberDJFixed_ = 0.0;
  if (probingInfo_) {
    delete probingInfo_;
    probingInfo_ = NULL;
    if (!generator_)
      numberCutGenerators_ = 0;
    // also get rid of cut generator
    int n = 0;
    for (int i = 0; i < numberCutGenerators_; i++) {
      CglImplication *cutGen;
      cutGen = dynamic_cast< CglImplication * >(generator_[i]->generator());
      if (!cutGen) {
        generator_[n] = generator_[i];
        virginGenerator_[n] = virginGenerator_[i];
        n++;
      } else {
        cutGen->setProbingInfo(NULL);
        delete generator_[i];
        cutGen = dynamic_cast< CglImplication * >(virginGenerator_[i]->generator());
        assert(cutGen);
        cutGen->setProbingInfo(NULL);
        delete virginGenerator_[i];
      }
    }
    numberCutGenerators_ = n;
  }
  maximumStatistics_ = 0;
  delete[] analyzeResults_;
  analyzeResults_ = NULL;
  bestObjective_ = COIN_DBL_MAX;
  bestPossibleObjective_ = COIN_DBL_MAX;
  sumChangeObjective1_ = 0.0;
  sumChangeObjective2_ = 0.0;
  numberSolutions_ = 0;
  stateOfSearch_ = 0;
  delete[] hotstartSolution_;
  hotstartSolution_ = NULL;
  delete[] hotstartPriorities_;
  hotstartPriorities_ = NULL;
  numberHeuristicSolutions_ = 0;
  numberNodes_ = 0;
  numberNodes2_ = 0;
  numberIterations_ = 0;
  numberSolves_ = 0;
  status_ = -1;
  secondaryStatus_ = -1;
  maximumNumberCuts_ = 0;
  phase_ = 0;
  currentNumberCuts_ = 0;
  maximumDepth_ = 0;
  nextRowCut_ = NULL;
  currentNode_ = NULL;
  // clear out tree
  if (tree_ && tree_->size())
    tree_->cleanTree(this, -1.0e100, bestPossibleObjective_);
  subTreeModel_ = NULL;
  heuristicModel_ = NULL;
  numberStoppedSubTrees_ = 0;
  numberInfeasibleNodes_ = 0;
  numberGlobalViolations_ = 0;
  numberExtraIterations_ = 0;
  numberExtraNodes_ = 0;
  numberFathoms_ = 0;
  continuousObjective_ = 0.0;
  originalContinuousObjective_ = 0.0;
  continuousInfeasibilities_ = 0;
  numberFixedAtRoot_ = 0;
  numberFixedNow_ = 0;
  stoppedOnGap_ = false;
  eventHappened_ = false;
  numberLongStrong_ = 0;
  numberOldActiveCuts_ = 0;
  numberNewCuts_ = 0;
  searchStrategy_ = -1;
  strongStrategy_ = 0;
  numberStrongIterations_ = 0;
  // Parameters which need to be reset
  setCutoff(COIN_DBL_MAX);
  dblParam_[CbcCutoffIncrement] = 1e-5;
  dblParam_[CbcCurrentCutoff] = 1.0e100;
  dblParam_[CbcCurrentObjectiveValue] = 1.0e100;
  dblParam_[CbcCurrentMinimizationObjectiveValue] = 1.0e100;
  delete globalConflictCuts_;
  globalConflictCuts_ = NULL;
}
/* Most of copy constructor
      mode - 0 copy but don't delete before
             1 copy and delete before
	     2 copy and delete before (but use virgin generators)
*/
void CbcModel::gutsOfCopy(const CbcModel &rhs, int mode)
{
  minimumDrop_ = rhs.minimumDrop_;
  specialOptions_ = rhs.specialOptions_;
  moreSpecialOptions_ = rhs.moreSpecialOptions_;
  moreSpecialOptions2_ = rhs.moreSpecialOptions2_;
  numberStrong_ = rhs.numberStrong_;
  numberBeforeTrust_ = rhs.numberBeforeTrust_;
  numberPenalties_ = rhs.numberPenalties_;
  printFrequency_ = rhs.printFrequency_;
  fastNodeDepth_ = rhs.fastNodeDepth_;
  howOftenGlobalScan_ = rhs.howOftenGlobalScan_;
  maximumCutPassesAtRoot_ = rhs.maximumCutPassesAtRoot_;
  maximumCutPasses_ = rhs.maximumCutPasses_;
  randomSeed_ = rhs.randomSeed_;
  multipleRootTries_ = rhs.multipleRootTries_;
  preferredWay_ = rhs.preferredWay_;
  resolveAfterTakeOffCuts_ = rhs.resolveAfterTakeOffCuts_;
  maximumNumberIterations_ = rhs.maximumNumberIterations_;
  numberSavedSolutions_ = rhs.numberSavedSolutions_;
  maximumSavedSolutions_ = rhs.maximumSavedSolutions_;
  if (maximumSavedSolutions_) {
    int n = solver_->getNumCols();
    savedSolutions_ = new double *[maximumSavedSolutions_];
    for (int i = 0; i < maximumSavedSolutions_; i++)
      savedSolutions_[i] = CoinCopyOfArray(rhs.savedSolutions_[i], n + 2);
  }
  continuousPriority_ = rhs.continuousPriority_;
  numberThreads_ = rhs.numberThreads_;
  threadMode_ = rhs.threadMode_;
  numberGlobalCutsIn_ = rhs.numberGlobalCutsIn_;
  delete master_;
  master_ = NULL;
  masterThread_ = NULL;
  memcpy(intParam_, rhs.intParam_, sizeof(intParam_));
  memcpy(dblParam_, rhs.dblParam_, sizeof(dblParam_));
  int i;
  if (mode) {
    for (i = 0; i < numberCutGenerators_; i++) {
      delete generator_[i];
      delete virginGenerator_[i];
    }
    delete[] generator_;
    delete[] virginGenerator_;
    for (i = 0; i < numberHeuristics_; i++) {
      delete heuristic_[i];
    }
    delete[] heuristic_;
    delete eventHandler_;
    delete branchingMethod_;
  }
  numberCutGenerators_ = rhs.numberCutGenerators_;
  if (numberCutGenerators_) {
    generator_ = new CbcCutGenerator *[numberCutGenerators_];
    virginGenerator_ = new CbcCutGenerator *[numberCutGenerators_];
    int i;
    for (i = 0; i < numberCutGenerators_; i++) {
      if (mode < 2) {
        generator_[i] = new CbcCutGenerator(*rhs.generator_[i]);
      } else {
        generator_[i] = new CbcCutGenerator(*rhs.virginGenerator_[i]);
        // But copy across maximumTries and switches
        generator_[i]->setMaximumTries(rhs.generator_[i]->maximumTries());
        generator_[i]->setSwitches(rhs.generator_[i]->switches());
      }
      virginGenerator_[i] = new CbcCutGenerator(*rhs.virginGenerator_[i]);
    }
  } else {
    generator_ = NULL;
    virginGenerator_ = NULL;
  }
  numberHeuristics_ = rhs.numberHeuristics_;
  if (numberHeuristics_) {
    heuristic_ = new CbcHeuristic *[numberHeuristics_];
    int i;
    for (i = 0; i < numberHeuristics_; i++) {
      heuristic_[i] = rhs.heuristic_[i]->clone();
    }
  } else {
    heuristic_ = NULL;
  }
  if (rhs.eventHandler_)
    eventHandler_ = rhs.eventHandler_->clone();
  else
    eventHandler_ = NULL;
  if (rhs.branchingMethod_)
    branchingMethod_ = rhs.branchingMethod_->clone();
  else
    branchingMethod_ = NULL;
  messageHandler()->setLogLevel(rhs.messageHandler()->logLevel());
  whenCuts_ = rhs.whenCuts_;
#ifdef COIN_HAS_NTY
  if (rhs.symmetryInfo_)
    symmetryInfo_ = new CbcSymmetry(*rhs.symmetryInfo_);
  else
    symmetryInfo_ = NULL;
#endif
  synchronizeModel();
}
// Move status, nodes etc etc across
void CbcModel::moveInfo(const CbcModel &rhs)
{
  bestObjective_ = rhs.bestObjective_;
  bestPossibleObjective_ = rhs.bestPossibleObjective_;
  numberSolutions_ = rhs.numberSolutions_;
  numberHeuristicSolutions_ = rhs.numberHeuristicSolutions_;
  numberNodes_ = rhs.numberNodes_;
  numberNodes2_ = rhs.numberNodes2_;
  numberIterations_ = rhs.numberIterations_;
  numberSolves_ = rhs.numberSolves_;
  status_ = rhs.status_;
  secondaryStatus_ = rhs.secondaryStatus_;
  numberStoppedSubTrees_ = rhs.numberStoppedSubTrees_;
  numberInfeasibleNodes_ = rhs.numberInfeasibleNodes_;
  continuousObjective_ = rhs.continuousObjective_;
  originalContinuousObjective_ = rhs.originalContinuousObjective_;
  continuousInfeasibilities_ = rhs.continuousInfeasibilities_;
  numberFixedAtRoot_ = rhs.numberFixedAtRoot_;
  numberFixedNow_ = rhs.numberFixedNow_;
  stoppedOnGap_ = rhs.stoppedOnGap_;
  eventHappened_ = rhs.eventHappened_;
  numberLongStrong_ = rhs.numberLongStrong_;
  numberStrongIterations_ = rhs.numberStrongIterations_;
  strongInfo_[0] = rhs.strongInfo_[0];
  strongInfo_[1] = rhs.strongInfo_[1];
  strongInfo_[2] = rhs.strongInfo_[2];
  strongInfo_[3] = rhs.strongInfo_[3];
  strongInfo_[4] = rhs.strongInfo_[4];
  strongInfo_[5] = rhs.strongInfo_[5];
  strongInfo_[6] = rhs.strongInfo_[6];
  numberRowsAtContinuous_ = rhs.numberRowsAtContinuous_;
  cutoffRowNumber_ = rhs.cutoffRowNumber_;
  maximumDepth_ = rhs.maximumDepth_;
}
// Save a copy of the current solver so can be reset to
void CbcModel::saveReferenceSolver()
{
  delete referenceSolver_;
  referenceSolver_ = solver_->clone();
}

// Uses a copy of reference solver to be current solver
void CbcModel::resetToReferenceSolver()
{
  delete solver_;
  solver_ = referenceSolver_->clone();
  // clear many things
  gutsOfDestructor2();
  // Reset cutoff
  // Solvers know about direction
  double direction = solver_->getObjSense();
  double value;
  solver_->getDblParam(OsiDualObjectiveLimit, value);
  setCutoff(value * direction);
}

// Are there a numerical difficulties?
bool CbcModel::isAbandoned() const
{
  return status_ == 2;
}
// Is optimality proven?
bool CbcModel::isProvenOptimal() const
{
  if (!status_ && bestObjective_ < 1.0e30)
    return true;
  else
    return false;
}
// Is  infeasiblity proven (or none better than cutoff)?
bool CbcModel::isProvenInfeasible() const
{
  if (!status_ && (bestObjective_ >= 1.0e30 && (secondaryStatus_ == 0 || secondaryStatus_ == 1)))
    return true;
  else
    return false;
}
// Was continuous solution unbounded
bool CbcModel::isContinuousUnbounded() const
{
  if (!status_ && secondaryStatus_ == 7)
    return true;
  else
    return false;
}
// Was continuous solution unbounded
bool CbcModel::isProvenDualInfeasible() const
{
  if (!status_ && secondaryStatus_ == 7)
    return true;
  else
    return false;
}
// Node limit reached?
bool CbcModel::isNodeLimitReached() const
{
  return numberNodes_ >= intParam_[CbcMaxNumNode];
}
// Time limit reached?
bool CbcModel::isSecondsLimitReached() const
{
  if (status_ == 1 && secondaryStatus_ == 4)
    return true;
  else
    return false;
}
// Solution limit reached?
bool CbcModel::isSolutionLimitReached() const
{
  return numberSolutions_ >= intParam_[CbcMaxNumSol];
}
// Set language
void CbcModel::newLanguage(CoinMessages::Language language)
{
  messages_ = CbcMessage(language);
}
void CbcModel::setNumberStrong(int number)
{
  if (number < 0)
    numberStrong_ = 0;
  else
    numberStrong_ = number;
}
void CbcModel::setNumberBeforeTrust(int number)
{
  if (number < -3) {
    numberBeforeTrust_ = 0;
  } else {
    numberBeforeTrust_ = number;
    //numberStrong_ = CoinMax(numberStrong_,1);
  }
}
void CbcModel::setNumberPenalties(int number)
{
  if (number <= 0) {
    numberPenalties_ = 0;
  } else {
    numberPenalties_ = number;
  }
}
void CbcModel::setPenaltyScaleFactor(double value)
{
  if (value <= 0) {
    penaltyScaleFactor_ = 3.0;
  } else {
    penaltyScaleFactor_ = value;
  }
}
void CbcModel::setHowOftenGlobalScan(int number)
{
  if (number < -1)
    howOftenGlobalScan_ = 0;
  else
    howOftenGlobalScan_ = number;
}

// Add one generator
void CbcModel::addCutGenerator(CglCutGenerator *generator,
  int howOften, const char *name,
  bool normal, bool atSolution,
  bool whenInfeasible, int howOftenInSub,
  int whatDepth, int whatDepthInSub)
{
  CbcCutGenerator **temp = generator_;
  generator_ = new CbcCutGenerator *[numberCutGenerators_ + 1];
  memcpy(generator_, temp, numberCutGenerators_ * sizeof(CbcCutGenerator *));
  delete[] temp;
  generator_[numberCutGenerators_] = new CbcCutGenerator(this, generator, howOften, name,
    normal, atSolution, whenInfeasible, howOftenInSub,
    whatDepth, whatDepthInSub);
  // and before any changes
  temp = virginGenerator_;
  virginGenerator_ = new CbcCutGenerator *[numberCutGenerators_ + 1];
  memcpy(virginGenerator_, temp, numberCutGenerators_ * sizeof(CbcCutGenerator *));
  delete[] temp;
  virginGenerator_[numberCutGenerators_++] = new CbcCutGenerator(this, generator, howOften, name,
    normal, atSolution, whenInfeasible, howOftenInSub,
    whatDepth, whatDepthInSub);
}
// Add one heuristic
void CbcModel::addHeuristic(CbcHeuristic *generator, const char *name,
  int before)
{
  CbcHeuristic **temp = heuristic_;
  heuristic_ = new CbcHeuristic *[numberHeuristics_ + 1];
  memcpy(heuristic_, temp, numberHeuristics_ * sizeof(CbcHeuristic *));
  delete[] temp;
  int where;
  if (before < 0 || before >= numberHeuristics_) {
    where = numberHeuristics_;
  } else {
    // move up
    for (int i = numberHeuristics_; i > before; i--)
      heuristic_[i] = heuristic_[i - 1];
    where = before;
  }
  heuristic_[where] = generator->clone();
  if (name)
    heuristic_[where]->setHeuristicName(name);
#ifndef SAME_HEURISTIC_SEED
  heuristic_[where]->setSeed(987654321 + where);
#else
  heuristic_[where]->setSeed(987654321);
#endif
  numberHeuristics_++;
}

/*
  The last subproblem handled by the solver is not necessarily related to the
  one being recreated, so the first action is to remove all cuts from the
  constraint system.  Next, traverse the tree from node to the root to
  determine the basis size required for this subproblem and create an empty
  basis with the right capacity.  Finally, traverse the tree from root to
  node, adjusting bounds in the constraint system, adjusting the basis, and
  collecting the cuts that must be added to the constraint system.
  applyToModel does the heavy lifting.

  addCuts1 is used in contexts where all that's desired is the list of cuts:
  the node is already fathomed, and we're collecting cuts so that we can
  adjust reference counts as we prune nodes. Arguably the two functions
  should be separated. The culprit is applyToModel, which performs cut
  collection and model adjustment.

  Certainly in the contexts where all we need is a list of cuts, there's no
  point in passing in a valid basis --- an empty basis will do just fine.
*/
bool CbcModel::addCuts1(CbcNode *node, CoinWarmStartBasis *&lastws)
{
  int nNode = 0;
  CbcNodeInfo *nodeInfo = node->nodeInfo();
  int numberColumns = getNumCols();

  /*
      Accumulate the path from node to the root in walkback_, and accumulate a
      cut count in currentNumberCuts.

      original comment: when working then just unwind until where new node joins
      old node (for cuts?)
    */
  int currentNumberCuts = 0;
  while (nodeInfo) {
    //printf("nNode = %d, nodeInfo = %x\n",nNode,nodeInfo);
    walkback_[nNode++] = nodeInfo;
    currentNumberCuts += nodeInfo->numberCuts();
    nodeInfo = nodeInfo->parent();
    if (nNode == maximumDepth_) {
      redoWalkBack();
    }
  }
  resizeWhichGenerator(currentNumberCuts_, currentNumberCuts);
  currentNumberCuts_ = currentNumberCuts;
  if (currentNumberCuts > maximumNumberCuts_) {
    maximumNumberCuts_ = currentNumberCuts;
    delete[] addedCuts_;
    addedCuts_ = new CbcCountRowCut *[maximumNumberCuts_];
  }
  /*
      This last bit of code traverses the path collected in walkback_ from the
      root back to node. At the end of the loop,
       * lastws will be an appropriate basis for node;
       * variable bounds in the constraint system will be set to be correct for
         node; and
       * addedCuts_ will be set to a list of cuts that need to be added to the
         constraint system at node.
      applyToModel does all the heavy lifting.
    */
  bool sameProblem = false;
  if ((specialOptions_ & 4096) == 0) {
#if 0
        {
            int n1 = numberRowsAtContinuous_;
            for (int i = 0; i < lastDepth_; i++)
                n1 += lastNumberCuts_[i];
            int n2 = numberRowsAtContinuous_;
            for (int i = 0; i < nNode; i++)
                n2 += walkback_[i]->numberCuts();
            //printf("ROWS a %d - old thinks %d new %d\n",solver_->getNumRows(),n1,n2);
        }
#endif
    int nDel = 0;
    int nAdd = 0;
    int n = CoinMin(lastDepth_, nNode);
    int i;
    int difference = lastDepth_ - nNode;
    int iZ = lastDepth_;
    int iN = 0;
    // Last is reversed to minimize copying
    if (difference > 0) {
      for (i = 0; i < difference; i++) {
        // delete rows
        nDel += lastNumberCuts_[--iZ];
      }
    } else if (difference < 0) {
      for (i = 0; i < -difference; i++) {
        // add rows
        nAdd += walkback_[i]->numberCuts();
      }
      iN = -difference;
    }
    for (i = 0; i < n; i++) {
      iZ--;
      if (lastNodeInfo_[iZ] == walkback_[iN]) {
        break;
      } else {
        // delete rows
        nDel += lastNumberCuts_[iZ];
        // add rows
        nAdd += walkback_[iN++]->numberCuts();
      }
    }
    assert(i < n || lastDepth_ == 0);
    //printf("lastDepth %d thisDepth %d match at %d, rows+-= %d %d\n",
    //   lastDepth_,nNode,n-i,nAdd,nDel);
    sameProblem = (!nAdd) && (!nDel);
    if (lastDepth_) {
      while (iN >= 0) {
        lastNumberCuts_[iZ] = walkback_[iN]->numberCuts();
        lastNodeInfo_[iZ++] = walkback_[iN--];
      }
    } else {
      lastNumberCuts_[0] = walkback_[0]->numberCuts();
      lastNodeInfo_[0] = walkback_[0];
    }
    lastDepth_ = nNode;
  }
  currentDepth_ = nNode;
  /*
      Remove all cuts from the constraint system.
      (original comment includes ``see note below for later efficiency'', but
      the reference isn't clear to me).
    */
  /*
      Create an empty basis with sufficient capacity for the constraint system
      we'll construct: original system plus cuts. Make sure we have capacity to
      record those cuts in addedCuts_.

      The method of adjusting the basis at a FullNodeInfo object (the root, for
      example) is to use a copy constructor to duplicate the basis held in the
      nodeInfo, then resize it and return the new basis object. Guaranteed,
      lastws will point to a different basis when it returns. We pass in a basis
      because we need the parameter to return the allocated basis, and it's an
      easy way to pass in the size. But we take a hit for memory allocation.
    */
  if (lastws)
    lastws->setSize(numberColumns, numberRowsAtContinuous_ + currentNumberCuts);
  currentNumberCuts = 0;
  while (nNode) {
    --nNode;
    walkback_[nNode]->applyToModel(this, lastws,
      addedCuts_, currentNumberCuts);
  }
#ifndef NDEBUG
  if (lastws && !lastws->fullBasis()) {
#ifdef COIN_DEVELOP
    printf("******* bad basis\n");
#endif
    int numberRows = lastws->getNumArtificial();
    int i;
    for (i = 0; i < numberRows; i++)
      lastws->setArtifStatus(i, CoinWarmStartBasis::basic);
    int numberColumns = lastws->getNumStructural();
    for (i = 0; i < numberColumns; i++) {
      if (lastws->getStructStatus(i) == CoinWarmStartBasis::basic)
        lastws->setStructStatus(i, CoinWarmStartBasis::atLowerBound);
    }
  }
#endif
  return sameProblem;
}

/*
  adjustCuts might be a better name: If the node is feasible, we sift through
  the cuts collected by addCuts1, add the ones that are tight and omit the
  ones that are loose. If the node is infeasible, we just adjust the
  reference counts to reflect that we're about to prune this node and its
  descendants.
*/
int CbcModel::addCuts(CbcNode *node, CoinWarmStartBasis *&lastws)
{
  /*
      addCuts1 performs step 1 of restoring the subproblem at this node; see the
      comments there.
    */
  bool sameProblem = addCuts1(node, lastws);
  int i;
  int numberColumns = getNumCols();
  if (solver_->getNumRows() > maximumRows_) {
    maximumRows_ = solver_->getNumRows();
    workingBasis_.resize(maximumRows_, numberColumns);
  }
  CbcNodeInfo *nodeInfo = node->nodeInfo();
  double cutoff = getCutoff();
  int currentNumberCuts = currentNumberCuts_;
  /*
      If the node can't be fathomed by bound, reinstall tight cuts in the
      constraint system. Even if there are no cuts, we'll want to set the
      reconstructed basis in the solver.
    */
  if (node->objectiveValue() < cutoff || numberThreads_) {
    //#   define CBC_CHECK_BASIS
#ifdef CBC_CHECK_BASIS
    printf("addCuts: expanded basis; rows %d+%d\n",
      numberRowsAtContinuous_, currentNumberCuts);
    lastws->print();
#endif
    /*
          Adjust the basis and constraint system so that we retain only active cuts.
          There are three steps:
            1) Scan the basis. Sort the cuts into effective cuts to be kept and
               loose cuts to be dropped.
            2) Drop the loose cuts and resize the basis to fit.
            3) Install the tight cuts in the constraint system (applyRowCuts) and
               and install the basis (setWarmStart).
          Use of compressRows conveys we're compressing the basis and not just
          tweaking the artificialStatus_ array.
        */
    if (currentNumberCuts > 0) {
      int numberToAdd = 0;
      const OsiRowCut **addCuts;
      int numberToDrop = 0;
      int *cutsToDrop;
      addCuts = new const OsiRowCut *[currentNumberCuts];
      cutsToDrop = new int[currentNumberCuts];
      assert(currentNumberCuts + numberRowsAtContinuous_ <= lastws->getNumArtificial());
      assert(currentNumberCuts <= maximumWhich_); // we will read from whichGenerator_[0..currentNumberCuts-1] below, so should have all these entries
      for (i = 0; i < currentNumberCuts; i++) {
        CoinWarmStartBasis::Status status = lastws->getArtifStatus(i + numberRowsAtContinuous_);
        if (addedCuts_[i] && (status != CoinWarmStartBasis::basic || (addedCuts_[i]->effectiveness() > 1.0e10 && !addedCuts_[i]->canDropCut(solver_, i + numberRowsAtContinuous_)))) {
#ifdef CHECK_CUT_COUNTS
          printf("Using cut %d %x as row %d\n", i, addedCuts_[i],
            numberRowsAtContinuous_ + numberToAdd);
#endif
          assert(i < maximumWhich_);
          whichGenerator_[numberToAdd] = whichGenerator_[i];
          addCuts[numberToAdd++] = addedCuts_[i];
#if 1
          if ((specialOptions_ & 1) != 0) {
            const OsiRowCutDebugger *debugger = solver_->getRowCutDebugger();
            if (debugger)
              CoinAssert(!debugger->invalidCut(*addedCuts_[i]));
          }
#endif
        } else {
#ifdef CHECK_CUT_COUNTS
          printf("Dropping cut %d %x\n", i, addedCuts_[i]);
#endif
          addedCuts_[i] = NULL;
          cutsToDrop[numberToDrop++] = numberRowsAtContinuous_ + i;
        }
      }
      assert(lastws->fullBasis());
      int numberRowsNow = numberRowsAtContinuous_ + numberToAdd;
      lastws->compressRows(numberToDrop, cutsToDrop);
      lastws->resize(numberRowsNow, numberColumns);
      // Take out as local search can give bad basisassert (lastws->fullBasis());
      bool canMissStuff = false;
      if ((specialOptions_ & 4096) == 0) {
        bool redoCuts = true;
        if (CoinAbs(lastNumberCuts2_ - numberToAdd) < 5) {
          int numberToCheck = CoinMin(lastNumberCuts2_, numberToAdd);
          int i1 = 0;
          int i2 = 0;
          int nDiff = 0;
          int nSame = 0;
          if (lastNumberCuts2_ == numberToAdd) {
            for (int i = 0; i < numberToCheck; i++) {
              if (lastCut_[i1++] != addCuts[i2++]) {
                nDiff++;
              } else {
                nSame++;
              }
            }
          } else if (lastNumberCuts2_ > numberToAdd) {
            int nDiff2 = lastNumberCuts2_ - numberToAdd;
            for (int i = 0; i < numberToCheck; i++) {
              if (lastCut_[i1] != addCuts[i2]) {
                nDiff++;
                while (nDiff2) {
                  i1++;
                  nDiff2--;
                  if (lastCut_[i1] == addCuts[i2]) {
                    nSame++;
                    break;
                  } else {
                    nDiff++;
                  }
                }
              } else {
                nSame++;
              }
            }
            nDiff += nDiff2;
          } else {
            int nDiff2 = numberToAdd - lastNumberCuts2_;
            for (int i = 0; i < numberToCheck; i++) {
              if (lastCut_[i1] != addCuts[i2]) {
                nDiff++;
                while (nDiff2) {
                  i2++;
                  nDiff2--;
                  if (lastCut_[i1] == addCuts[i2]) {
                    nSame++;
                    break;
                  } else {
                    nDiff++;
                  }
                }
              } else {
                nSame++;
              }
            }
            nDiff += nDiff2;
          }
          canMissStuff = !nDiff && sameProblem;
          // But only if number of rows looks OK
          if (numberRowsAtContinuous_ + numberToAdd != solver_->getNumRows())
            canMissStuff = false;
        } else {
          //printf("add now %d add last %d NO2\n",numberToAdd,lastNumberCuts2_);
        }
        assert(lastws->fullBasis() && numberRowsAtContinuous_ + numberToAdd == numberRowsNow);
        if (redoCuts) {
          if (numberToAdd > maximumCuts_) {
            delete[] lastCut_;
            maximumCuts_ = 2 * numberToAdd + 10;
            lastCut_ = new const OsiRowCut *[maximumCuts_];
          }
          lastNumberCuts2_ = numberToAdd;
          for (int i = 0; i < numberToAdd; i++)
            lastCut_[i] = addCuts[i];
        }
      }
      if (!canMissStuff) {
        //if (canMissStuff)
        //solver_->writeMps("before");
        //printf("Not Skipped\n");
        //int n1=solver_->getNumRows();
        if ((specialOptions_ & 4096) == 0) {
          solver_->restoreBaseModel(numberRowsAtContinuous_);
        } else {
          // *** Fix later
          int numberCuts = solver_->getNumRows() - numberRowsAtContinuous_;
          int *which = new int[numberCuts];
          for (i = 0; i < numberCuts; i++)
            which[i] = i + numberRowsAtContinuous_;
          solver_->deleteRows(numberCuts, which);
          delete[] which;
        }
        //#define CHECK_DEBUGGER
#ifdef CHECK_DEBUGGER
        if ((specialOptions_ & 1) != 0) {
          const OsiRowCutDebugger *debugger = solver_->getRowCutDebugger();
          if (debugger) {
            for (int j = 0; j < numberToAdd; j++)
              CoinAssert(!debugger->invalidCut(*addCuts[j]));
            //addCuts[j]->print();
          }
        }
#endif
        solver_->applyRowCuts(numberToAdd, addCuts);
      }
#ifdef CBC_CHECK_BASIS
      printf("addCuts: stripped basis; rows %d + %d\n",
        numberRowsAtContinuous_, numberToAdd);
      lastws->print();
#endif
      delete[] addCuts;
      delete[] cutsToDrop;
    }
    /*
          Set the basis in the solver.
        */
    solver_->setWarmStart(lastws);
    /*
          Clean up and we're out of here.
        */
    numberNodes_++;
    return 0;
  }
  /*
      This node has been fathomed by bound as we try to revive it out of the live
      set. Adjust the cut reference counts to reflect that we no longer need to
      explore the remaining branch arms, hence they will no longer reference any
      cuts. Cuts whose reference count falls to zero are deleted.
    */
  else {
    int i;
    if (currentNumberCuts) {
      lockThread();
      int numberLeft = nodeInfo->numberBranchesLeft();
      for (i = 0; i < currentNumberCuts; i++) {
        if (addedCuts_[i]) {
          if (!addedCuts_[i]->decrement(numberLeft)) {
            delete addedCuts_[i];
            addedCuts_[i] = NULL;
          }
        }
      }
      unlockThread();
    }
    return 1;
  }
}
/* Makes all handlers same.  If makeDefault 1 then makes top level
   default and rest point to that.  If 2 then each is copy
*/
void CbcModel::synchronizeHandlers(int /*makeDefault*/)
{
  bool defaultHandler = defaultHandler_;
  if (!defaultHandler_) {
    // Must have clone
    handler_ = handler_->clone();  // Not sure - worst is small memory leak
    defaultHandler_ = true;
  }
#ifdef COIN_HAS_CLP
  if (!defaultHandler) {
    OsiClpSolverInterface *solver;
    solver = dynamic_cast< OsiClpSolverInterface * >(solver_);
    if (solver) {
      solver->passInMessageHandler(handler_);
      solver->getModelPtr()->passInMessageHandler(handler_);
    }
    solver = dynamic_cast< OsiClpSolverInterface * >(continuousSolver_);
    if (solver) {
      solver->passInMessageHandler(handler_);
      solver->getModelPtr()->passInMessageHandler(handler_);
    }
  }
#endif
}

/*
  Perform reduced cost fixing on integer variables.

  The variables in question are already nonbasic at bound. We're just nailing
  down the current situation.
*/
int CbcModel::reducedCostFix()

{
  if (!solverCharacteristics_->reducedCostsAccurate())
    return 0; //NLP
  double cutoff = getCutoff();
  double direction = solver_->getObjSense();
  double gap = cutoff - solver_->getObjValue() * direction;
  double tolerance;
  solver_->getDblParam(OsiDualTolerance, tolerance);
  if (gap <= 0.0)
    gap = tolerance; //return 0;
  gap += 100.0 * tolerance;
  double integerTolerance = getDblParam(CbcIntegerTolerance);

  const double *lower = solver_->getColLower();
  const double *upper = solver_->getColUpper();
  const double *solution = solver_->getColSolution();
  const double *reducedCost = solver_->getReducedCost();

  int numberFixed = 0;
  int numberTightened = 0;

#ifdef COIN_HAS_CLP
  OsiClpSolverInterface *clpSolver
    = dynamic_cast< OsiClpSolverInterface * >(solver_);
  ClpSimplex *clpSimplex = NULL;
  if (clpSolver)
    clpSimplex = clpSolver->getModelPtr();
#endif
  for (int i = 0; i < numberIntegers_; i++) {
    int iColumn = integerVariable_[i];
    double djValue = direction * reducedCost[iColumn];
    double boundGap = upper[iColumn] - lower[iColumn];
    if (boundGap > integerTolerance) {
      if (solution[iColumn] < lower[iColumn] + integerTolerance
        && djValue * boundGap > gap) {
#ifdef COIN_HAS_CLP
        // may just have been fixed before
        if (clpSimplex) {
          if (clpSimplex->getColumnStatus(iColumn) == ClpSimplex::basic) {
#ifdef COIN_DEVELOP
            printf("DJfix %d has status of %d, dj of %g gap %g, bounds %g %g\n",
              iColumn, clpSimplex->getColumnStatus(iColumn),
              djValue, gap, lower[iColumn], upper[iColumn]);
#endif
          } else {
            assert(clpSimplex->getColumnStatus(iColumn) == ClpSimplex::atLowerBound || clpSimplex->getColumnStatus(iColumn) == ClpSimplex::isFixed);
          }
        }
#endif
        double newBound = lower[iColumn];
        if (boundGap > 1.99) {
          boundGap = gap / djValue + 1.0e-4 * boundGap;
          newBound = lower[iColumn] + floor(boundGap);
          numberTightened++;
          //if (newBound)
          //printf("tighter - gap %g dj %g newBound %g\n",
          //   gap,djValue,newBound);
        }
        solver_->setColUpper(iColumn, newBound);
        numberFixed++;
      } else if (solution[iColumn] > upper[iColumn] - integerTolerance && -djValue > boundGap * gap) {
#ifdef COIN_HAS_CLP
        // may just have been fixed before
        if (clpSimplex) {
          if (clpSimplex->getColumnStatus(iColumn) == ClpSimplex::basic) {
#ifdef COIN_DEVELOP
            printf("DJfix %d has status of %d, dj of %g gap %g, bounds %g %g\n",
              iColumn, clpSimplex->getColumnStatus(iColumn),
              djValue, gap, lower[iColumn], upper[iColumn]);
#endif
          } else {
            assert(clpSimplex->getColumnStatus(iColumn) == ClpSimplex::atUpperBound || clpSimplex->getColumnStatus(iColumn) == ClpSimplex::isFixed);
          }
        }
#endif
        double newBound = upper[iColumn];
        if (boundGap > 1.99) {
          boundGap = -gap / djValue + 1.0e-4 * boundGap;
          newBound = upper[iColumn] - floor(boundGap);
          //if (newBound)
          //printf("tighter - gap %g dj %g newBound %g\n",
          //   gap,djValue,newBound);
          numberTightened++;
        }
        solver_->setColLower(iColumn, newBound);
        numberFixed++;
      }
    }
  }
  numberDJFixed_ += numberFixed - numberTightened;
#ifdef SWITCH_VARIABLES
  if (numberFixed)
    fixAssociated(NULL, 0);
#endif
  return numberFixed;
}
// Collect coding to replace whichGenerator
void CbcModel::resizeWhichGenerator(int numberNow, int numberAfter)
{
  if (numberAfter > maximumWhich_) {
#define MAXIMUM_WHICH_INCREMENT 100
#define MAXIMUM_WHICH_MULTIPLIER 2
    //printf("maximumWhich from %d to %d (%d needed)\n",maximumWhich_,
    //   CoinMax(maximumWhich_ * MAXIMUM_WHICH_MULTIPLIER + MAXIMUM_WHICH_INCREMENT, numberAfter),
    //   numberAfter);
    maximumWhich_ = CoinMax(maximumWhich_ * MAXIMUM_WHICH_MULTIPLIER + MAXIMUM_WHICH_INCREMENT, numberAfter);
    //maximumWhich_ = numberAfter ;
    int *temp = new int[2 * maximumWhich_];
    memcpy(temp, whichGenerator_, numberNow * sizeof(int));
    delete[] whichGenerator_;
    whichGenerator_ = temp;
    memset(whichGenerator_ + numberNow, 0, (maximumWhich_ - numberNow) * sizeof(int));
  }
}

/** Solve the model using cuts

  This version takes off redundant cuts from node.
  Returns true if feasible.

  \todo
  Why do I need to resolve the problem? What has been done between the last
  relaxation and calling solveWithCuts?

  If numberTries == 0 then user did not want any cuts.
*/

bool CbcModel::solveWithCuts(OsiCuts &cuts, int numberTries, CbcNode *node)
/*
  Parameters:
    numberTries: (i) the maximum number of iterations for this round of cut
		     generation; if negative then we don't mind if drop is tiny.

    cuts:	(o) all cuts generated in this round of cut generation

    node: (i)     So we can update dynamic pseudo costs
*/

{
#ifdef JJF_ZERO
  if (node && numberTries > 1) {
    if (currentDepth_ < 5)
      numberTries *= 4; // boost
    else if (currentDepth_ < 10)
      numberTries *= 2; // boost
  }
#endif
#define CUT_HISTORY 7
  double cut_obj[CUT_HISTORY];
  for (int j = 0; j < CUT_HISTORY; j++)
    cut_obj[j] = -COIN_DBL_MAX;
#ifdef COIN_HAS_CLP
  OsiClpSolverInterface *clpSolver
    = dynamic_cast< OsiClpSolverInterface * >(solver_);
  int saveClpOptions = 0;
  if (clpSolver)
    saveClpOptions = clpSolver->specialOptions();
#endif
    //solver_->writeMps("saved");
#ifdef CBC_THREAD
  /*
      Thread mode makes a difference here only when it specifies using separate
      threads to generate cuts at the root (bit 2^1 set in threadMode_). In which
      case we'll create an array of empty CbcModels (!). Solvers will be cloned
      later.

      Don't start up threads here if we're already threaded.
    */
  CbcBaseModel *master = NULL;
  if (numberThreads_ && (threadMode_ & 2) != 0 && !numberNodes_) {
    master = new CbcBaseModel(*this, -1);
  }
#endif
  bool feasible = true;
  int violated = 0;
  int numberRowsAtStart = solver_->getNumRows();
  //printf("solver had %d rows\n",numberRowsAtStart);
  int numberColumns = solver_->getNumCols();
  CoinBigIndex numberElementsAtStart = solver_->getNumElements();

  numberOldActiveCuts_ = numberRowsAtStart - numberRowsAtContinuous_;
#ifndef NDEBUG
  {
    int n = 0;
    for (int i = 0; i < currentNumberCuts_; i++) {
      if (addedCuts_[i])
        n++;
    }
    assert(n == numberOldActiveCuts_);
  }
#endif
  numberNewCuts_ = cuts.sizeRowCuts();
  int lastNumberCuts = numberNewCuts_;
  if (numberNewCuts_) {
    // from multiple root passes
    const OsiRowCut **addCuts = new const OsiRowCut *[numberNewCuts_];
    for (int i = 0; i < numberNewCuts_; i++) {
      addCuts[i] = cuts.rowCutPtr(i);
    }
    solver_->applyRowCuts(numberNewCuts_, addCuts);
    delete[] addCuts;
  }

  bool onOptimalPath = false;
  const OsiRowCutDebugger *debugger = NULL;
  if ((specialOptions_ & 1) != 0) {
    /*
          See OsiRowCutDebugger for details. In a nutshell, make sure that current
          variable values do not conflict with a known optimal solution. (Obviously
          this can be fooled when there are multiple solutions.)
        */
    debugger = solver_->getRowCutDebugger();
    if (debugger)
      onOptimalPath = (debugger->onOptimalPath(*solver_));
  }
  /*
      As the final action in each round of cut generation (the numberTries loop),
      we'll call takeOffCuts to remove slack cuts. These are saved into slackCuts
      and rechecked immediately after the cut generation phase of the loop.
    */
  OsiCuts slackCuts;
  /*
    lh:
      Resolve the problem

    The resolve will also refresh cached copies of the solver solution
    (cbcColLower_, ...) held by CbcModel.
    This resolve looks like the best point to capture a warm start for use in
    the case where cut generation proves ineffective and we need to back out
    a few tight cuts.
    I've always maintained that this resolve is unnecessary. Let's put in a hook
    to report if it's every nontrivial. -lh

    Resolve the problem. If we've lost feasibility, might as well bail out right
      after the debug stuff. The resolve will also refresh cached copies of the
      solver solution (cbcColLower_, ...) held by CbcModel.
    */
  double objectiveValue = solver_->getObjValue() * solver_->getObjSense();
  if (node) {
    objectiveValue = node->objectiveValue();
  }
  int save = moreSpecialOptions_;
  if ((moreSpecialOptions_ & 4194304) != 0)
    moreSpecialOptions_ |= 8388608;
  int returnCode = resolve(node ? node->nodeInfo() : NULL, 1);
  moreSpecialOptions_ = save;
#ifdef CONFLICT_CUTS
#ifdef COIN_HAS_CLP
  // if infeasible conflict analysis
  if (solver_->isProvenPrimalInfeasible() && !parentModel_ && (moreSpecialOptions_ & 4194304) != 0 && clpSolver) {
    if (!topOfTree_ && masterThread_)
      topOfTree_ = masterThread_->master_->baseModel_->topOfTree_;
    assert(topOfTree_);
    int iType = 0;
    OsiRowCut *cut = clpSolver->modelCut(topOfTree_->lower(),
      topOfTree_->upper(),
      numberRowsAtContinuous_, whichGenerator_, iType);
    if (cut) {
      //cut->print();
      if (!iType) {
        int badCut = makeGlobalCut(cut);
        if (!badCut) {
#if PRINT_CONFLICT == 1
          numberConflictCuts++;
          lengthConflictCuts += cut->row().getNumElements();
#endif
#if PRINT_CONFLICT < 2
          if (handler_->logLevel() > 1) {
#endif
            printf("Conflict cut at depth %d (%d elements)\n",
              currentDepth_, cut->row().getNumElements());
            if (cut->row().getNumElements() < 3)
              cut->print();
#if PRINT_CONFLICT < 2
          }
#endif
        }
        if ((specialOptions_ & 1) != 0) {
          debugger = continuousSolver_->getRowCutDebugger();
          if (debugger) {
            if (debugger->invalidCut(*cut)) {
              continuousSolver_->applyRowCuts(1, cut);
              continuousSolver_->writeMps("bad");
            }
            CoinAssert(!debugger->invalidCut(*cut));
          }
        }
      } else {
        makePartialCut(cut);
      }
      delete cut;
    }
  }
  if ((moreSpecialOptions_ & 4194304) != 0 && solver_->isProvenPrimalInfeasible()
    && clpSolver && clpSolver->lastAlgorithm() == 2 && clpSolver->getModelPtr()->infeasibilityRay() && !parentModel_) {
    printf("ray exists\n");
  }
#endif
#endif
  double lastObjective = solver_->getObjValue() * solver_->getObjSense();
  cut_obj[CUT_HISTORY - 1] = lastObjective;
  //double firstObjective = lastObjective+1.0e-8+1.0e-12*fabs(lastObjective);
  /*
      Contemplate the result of the resolve.
        - CbcModel::resolve() has a hook that calls CbcStrategy::status to look
          over the solution. The net result is that resolve can return
          0 (infeasible), 1 (feasible), or -1 (feasible, but do no further work).
        - CbcFeasbililityBase::feasible() can return 0 (no comment),
          1 (pretend this is an integer solution), or -1 (pretend this is
          infeasible). As of 080104, this seems to be a stub to allow overrides,
          with a default implementation that always returns 0.

      Setting numberTries = 0 for `do no more work' is problematic. The main cut
      generation loop will still execute once, so we do not observe the `no
      further work' semantics.

      As best I can see, allBranchesGone is a null function as of 071220.
    */
  if (node && node->nodeInfo() && !node->nodeInfo()->numberBranchesLeft())
    node->nodeInfo()->allBranchesGone(); // can clean up
  feasible = returnCode != 0;
  if (returnCode < 0)
    numberTries = 0;
  if (problemFeasibility_->feasible(this, 0) < 0) {
    feasible = false; // pretend infeasible
  }
  //#define CHECK_KNOWN_SOLUTION
#ifdef CHECK_KNOWN_SOLUTION
  if (onOptimalPath && (solver_->isDualObjectiveLimitReached() || !feasible)) {
    printf("help 1\n");
  }
#endif
  /*
      NEW_UPDATE_OBJECT is defined to 0 when unthreaded (CBC_THREAD undefined), 2
      when threaded. No sign of 1 as of 071220.

      At present, there are two sets of hierarchies for branching classes. Call
      them CbcHier and OsiHier. For example, we have OsiBranchingObject, with
      children CbcBranchingObject and OsiTwoWayBranchingObject. All
      specialisations descend from one of these two children. Similarly, there is
      OsiObject, with children CbcObject and OsiObject2.

      In the original setup, there's a single CbcBranchDecision object attached
      to CbcModel (branchingMethod_). It has a field to hold the current CbcHier
      branching object, and the updateInformation routine reaches through the
      branching object to update the underlying CbcHier object.

      NEW_UPDATE_OBJECT = 0 would seem to assume the original setup. But,
      if we're using the OSI hierarchy for objects and branching, a call to a
      nontrivial branchingMethod_->updateInformation would have no effect (it
      would expect a CbcObject to work on) or perhaps crash.  For the
      default CbcBranchDefaultDecision, updateInformation is a noop (actually
      defined in the base CbcBranchDecision class).

      NEW_UPDATE_OBJECT = 2 looks like it's prepared to cope with either CbcHier or
      OsiHier, but it'll be executed only when threads are activated. See the
      comments below. The setup is scary.

      But ... if the OsiHier update actually reaches right through to the object
      list in the solver, it should work just fine in unthreaded mode. It would
      seem that the appropriate thing to do in unthreaded mode would be to choose
      between the existing code for NEW_UPDATE_OBJECT = 0 and the OsiHier code for
      NEW_UPDATE_OBJECT = 2. But I'm going to let John hash that out. The worst
      that can happen is inefficiency because I'm not properly updating an object.
    */

  // Update branching information if wanted
  if (node && branchingMethod_) {
    OsiBranchingObject *bobj = node->modifiableBranchingObject();
    CbcBranchingObject *cbcobj = dynamic_cast< CbcBranchingObject * >(bobj);
    if (cbcobj && cbcobj->object()) {
      CbcObject *object = cbcobj->object();
      CbcObjectUpdateData update = object->createUpdateInformation(solver_, node, cbcobj);
      // have to compute object number as not saved
      CbcSimpleInteger *simpleObject = static_cast< CbcSimpleInteger * >(object);
      int iObject = simpleObject->position();
#ifndef NDEBUG
      int iColumn = simpleObject->columnNumber();
      int jObject;
      for (jObject = 0; jObject < numberObjects_; jObject++) {
        simpleObject = static_cast< CbcSimpleInteger * >(object_[jObject]);
        if (simpleObject->columnNumber() == iColumn)
          break;
      }
      assert(jObject < numberObjects_ && iObject == jObject);
#else
#ifdef CBCMODEL_TIGHTEN_BOUNDS
      int iColumn = simpleObject->columnNumber();
#endif
#endif
      update.objectNumber_ = iObject;
      // Care! We must be careful not to update the same variable in parallel threads.
      addUpdateInformation(update);
      // update here
      {
        CbcObject *object = dynamic_cast< CbcObject * >(update.object_);
        if (object)
          object->updateInformation(update);
      }
      //#define CBCMODEL_TIGHTEN_BOUNDS
#ifdef CBCMODEL_TIGHTEN_BOUNDS
      double cutoff = getCutoff();
      if (feasible && cutoff < 1.0e20) {
        int way = cbcobj->way();
        // way is what will be taken next
        way = -way;
        double value = cbcobj->value();
        //const double * lower = solver_->getColLower();
        //const double * upper = solver_->getColUpper();
        double objectiveChange = lastObjective - objectiveValue;
        if (objectiveChange > 1.0e-5) {
          CbcIntegerBranchingObject *branch = dynamic_cast< CbcIntegerBranchingObject * >(cbcobj);
          assert(branch);
          if (way < 0) {
            double down = value - floor(value);
            double changePer = objectiveChange / (down + 1.0e-7);
            double distance = (cutoff - objectiveValue) / changePer;
            distance += 1.0e-3;
            if (distance < 5.0) {
              double newLower = ceil(value - distance);
              const double *downBounds = branch->downBounds();
              if (newLower > downBounds[0]) {
                //printf("%d way %d bounds %g %g value %g\n",
                //     iColumn,way,lower[iColumn],upper[iColumn],value);
                //printf("B Could increase lower bound on %d from %g to %g\n",
                //     iColumn,downBounds[0],newLower);
                solver_->setColLower(iColumn, newLower);
              }
            }
          } else {
            double up = ceil(value) - value;
            double changePer = objectiveChange / (up + 1.0e-7);
            double distance = (cutoff - objectiveValue) / changePer;
            distance += 1.0e-3;
            if (distance < 5.0) {
              double newUpper = floor(value + distance);
              const double *upBounds = branch->upBounds();
              if (newUpper < upBounds[1]) {
                //printf("%d way %d bounds %g %g value %g\n",
                //     iColumn,way,lower[iColumn],upper[iColumn],value);
                //printf("B Could decrease upper bound on %d from %g to %g\n",
                //     iColumn,upBounds[1],newUpper);
                solver_->setColUpper(iColumn, newUpper);
              }
            }
          }
        }
      }
#endif
    } else {
      OsiIntegerBranchingObject *obj = dynamic_cast< OsiIntegerBranchingObject * >(bobj);
      if (obj) {
        const OsiObject *object = obj->originalObject();
        // have to compute object number as not saved
        int iObject;
        int iColumn = object->columnNumber();
        for (iObject = 0; iObject < numberObjects_; iObject++) {
          if (object_[iObject]->columnNumber() == iColumn)
            break;
        }
        assert(iObject < numberObjects_);
        int branch = obj->firstBranch();
        if (obj->branchIndex() == 2)
          branch = 1 - branch;
        assert(branch == 0 || branch == 1);
        double originalValue = node->objectiveValue();
        double objectiveValue = solver_->getObjValue() * solver_->getObjSense();
        double changeInObjective = CoinMax(0.0, objectiveValue - originalValue);
        double value = obj->value();
        double movement;
        if (branch)
          movement = ceil(value) - value;
        else
          movement = value - floor(value);
        branchingMethod_->chooseMethod()->updateInformation(iObject, branch, changeInObjective,
          movement, 0 /*(feasible) ? 0 : 1; */);
      }
    }
  }

#ifdef CBC_DEBUG
  if (feasible) {
    printf("Obj value %g (%s) %d rows\n", solver_->getObjValue(),
      (solver_->isProvenOptimal()) ? "proven" : "unproven",
      solver_->getNumRows());
  }

  else {
    printf("Infeasible %d rows\n", solver_->getNumRows());
  }
#endif
  if ((specialOptions_ & 1) != 0) {
    /*
          If the RowCutDebugger said we were compatible with the optimal solution,
          and now we're suddenly infeasible, we might be confused. Then again, we
          may have fathomed by bound, heading for a rediscovery of an optimal solution.
        */
    if (onOptimalPath && !solver_->isDualObjectiveLimitReached()) {
      if (!feasible) {
        solver_->writeMpsNative("infeas.mps", NULL, NULL, 2);
        solver_->getRowCutDebuggerAlways()->printOptimalSolution(*solver_);
        CoinWarmStartBasis *slack = dynamic_cast< CoinWarmStartBasis * >(solver_->getEmptyWarmStart());
        solver_->setWarmStart(slack);
        delete slack;
        solver_->setHintParam(OsiDoReducePrint, false, OsiHintDo, 0);
        solver_->initialSolve();
      }
      assert(feasible);
    }
  }

  if (!feasible) {
    numberInfeasibleNodes_++;
#ifdef COIN_HAS_CLP
    if (clpSolver)
      clpSolver->setSpecialOptions(saveClpOptions);
#endif
    return (false);
  }
  double change = lastObjective - objectiveValue;
  if (change > 1.0e-10) {
    dblParam_[CbcSmallestChange] = CoinMin(dblParam_[CbcSmallestChange], change);
    dblParam_[CbcSumChange] += change;
    dblParam_[CbcLargestChange] = CoinMax(dblParam_[CbcLargestChange], change);
    intParam_[CbcNumberBranches]++;
  }
  sumChangeObjective1_ += solver_->getObjValue() * solver_->getObjSense()
    - objectiveValue;
  if (maximumSecondsReached())
    numberTries = 0; // exit
  if ((moreSpecialOptions2_ & (2048 | 4096)) != 0 && currentDepth_ > 5) {
    // howOftenGlobalScan_ = 10;
    int type = (moreSpecialOptions2_ & (2048 | 4096)) >> 11;
    if (type == 1) {
      int n = 0;
      int k = currentDepth_;
      while (k) {
        if ((k & 1) != 0)
          n++;
        k = k >> 1;
      }
      if (n > 1)
        numberTries = 0;
    } else if (type == 2) {
      if ((currentDepth_ % 4) != 0)
        numberTries = 0;
    } else {
      if ((currentDepth_ % 8) != 0)
        numberTries = 0;
    }
  }
  //if ((numberNodes_%100)==0)
  //printf("XXa sum obj changed by %g\n",sumChangeObjective1_);
  objectiveValue = solver_->getObjValue() * solver_->getObjSense();
  // Return at once if numberTries zero
  if (!numberTries) {
    cuts = OsiCuts();
    numberNewCuts_ = 0;
#ifdef COIN_HAS_CLP
    if (clpSolver)
      clpSolver->setSpecialOptions(saveClpOptions);
#endif
    setPointers(solver_);
    return true;
  }
  /*
      Do reduced cost fixing.
    */
  int xxxxxx = 0;
  if (xxxxxx)
    solver_->resolve();
  reducedCostFix();
  /*
      Set up for at most numberTries rounds of cut generation. If numberTries is
      negative, we'll ignore the minimumDrop_ cutoff and keep generating cuts for
      the specified number of rounds.
    */
  double minimumDrop = minimumDrop_;
  bool allowZeroIterations = false;
  int maximumBadPasses = 0;
  if (numberTries < 0) {
    numberTries = -numberTries;
    // minimumDrop *= 1.0e-5 ;
    // if (numberTries >= -1000000) {
    //numberTries=100;
    minimumDrop = -1.0;
    // }
    //numberTries=CoinMax(numberTries,100);
    allowZeroIterations = true;
  }
  int saveNumberTries = numberTries;
  /*
      Is it time to scan the cuts in order to remove redundant cuts? If so, set
      up to do it.
    */
  int fullScan = 0;
  if ((numberNodes_ % SCANCUTS) == 0 || (specialOptions_ & 256) != 0) {
    fullScan = 1;
    if (!numberNodes_ || (specialOptions_ & 256) != 0)
      fullScan = 2;
    specialOptions_ &= ~256; // mark as full scan done
  }

  double direction = solver_->getObjSense();
  double startObjective = solver_->getObjValue() * direction;

  currentPassNumber_ = 0;
  // Really primalIntegerTolerance; relates to an illposed problem with various
  // integer solutions depending on integer tolerance.
  //double primalTolerance = 1.0e-7 ;
  // We may need to keep going on
  bool keepGoing = false;
  // Say we have not tried one last time
  int numberLastAttempts = 0;
  /* Get experimental option as to when to stop adding cuts
       0 - old style
       1 - new style
       2 - new style plus don't break if zero cuts first time
       3 - as 2 but last drop has to be >0.1*min to say OK
    */
  int experimentBreak = (moreSpecialOptions_ >> 11) & 3;
  // Whether to increase minimum drop
  bool increaseDrop = (moreSpecialOptions_ & 8192) != 0;
  for (int i = 0; i < numberCutGenerators_; i++)
    generator_[i]->setWhetherInMustCallAgainMode(false);
  /*
      Begin cut generation loop. Cuts generated during each iteration are
      collected in theseCuts. The loop can be divided into four phases:
       1) Prep: Fix variables using reduced cost. In the first iteration only,
          consider scanning globalCuts_ and activating any applicable cuts.
       2) Cut Generation: Call each generator and heuristic registered in the
          generator_ and heuristic_ arrays. Newly generated global cuts are
          copied to globalCuts_ at this time.
       3) Cut Installation and Reoptimisation: Install column and row cuts in
          the solver. Copy row cuts to cuts (parameter). Reoptimise.
       4) Cut Purging: takeOffCuts() removes inactive cuts from the solver, and
          does the necessary bookkeeping in the model.
    */
  do {
    currentPassNumber_++;
    numberTries--;
    if (numberTries < 0 && keepGoing) {
      // switch off all normal generators (by the generator's opinion of normal)
      // Intended for situations where the primal problem really isn't complete,
      // and there are `not normal' cut generators that will augment.
      for (int i = 0; i < numberCutGenerators_; i++) {
        if (!generator_[i]->mustCallAgain())
          generator_[i]->setSwitchedOff(true);
        else
          generator_[i]->setWhetherInMustCallAgainMode(true);
      }
    }
    keepGoing = false;
    OsiCuts theseCuts;
    /*
          Scan previously generated global column and row cuts to see if any are
          useful.
        */
    int numberViolated = 0;
    if ((currentPassNumber_ == 1 || !numberNodes_) && howOftenGlobalScan_ > 0 && (numberNodes_ % howOftenGlobalScan_) == 0 && (doCutsNow(1) || true)) {
      // global column cuts now done in node at top of tree
      int numberCuts = numberCutGenerators_ ? globalCuts_.sizeRowCuts() : 0;
      if (numberCuts) {
        // possibly extend whichGenerator
        resizeWhichGenerator(numberViolated, numberViolated + numberCuts);
        // only add new cuts up to 10% of current elements
        CoinBigIndex numberElements = solver_->getNumElements();
        int numberColumns = solver_->getNumCols();
        CoinBigIndex maximumAdd = CoinMax(numberElements / 10,
                                    static_cast< CoinBigIndex >(2 * numberColumns))
          + 100;
        double *violations = new double[numberCuts];
        int *which = new int[numberCuts];
        int numberPossible = 0;
        for (int i = 0; i < numberCuts; i++) {
          OsiRowCut *thisCut = globalCuts_.rowCutPtr(i);
          double violation = thisCut->violated(cbcColSolution_);
          if (thisCut->effectiveness() == COIN_DBL_MAX) {
            // see if already there
            int j;
            for (j = 0; j < currentNumberCuts_; j++) {
              if (addedCuts_[j] == thisCut)
                break;
            }
            if (j == currentNumberCuts_)
              violation = COIN_DBL_MAX;
            //else
            //printf("already done??\n");
          }
          if (violation > 0.005) {
            violations[numberPossible] = -violation;
            which[numberPossible++] = i;
          }
        }
        CoinSort_2(violations, violations + numberPossible, which);
        for (int i = 0; i < numberPossible; i++) {
          int k = which[i];
          OsiRowCut *thisCut = globalCuts_.rowCutPtr(k);
          assert(thisCut->violated(cbcColSolution_) > 0.005 /*primalTolerance*/ || thisCut->effectiveness() == COIN_DBL_MAX);
#define CHECK_DEBUGGER
#ifdef CHECK_DEBUGGER
          if ((specialOptions_ & 1) != 0 && !parentModel_) {
            CoinAssert(!solver_->getRowCutDebuggerAlways()->invalidCut(*thisCut));
          }
#endif
#if 0 //ndef NDEBUG
		printf("Global cut added - violation %g\n",
		       thisCut->violated(cbcColSolution_)) ;
#endif
          whichGenerator_[numberViolated++] = 20099;
#ifndef GLOBAL_CUTS_JUST_POINTERS
          theseCuts.insert(*thisCut);
#else
          theseCuts.insert(thisCut);
#endif
          if (violations[i] != -COIN_DBL_MAX)
            maximumAdd -= thisCut->row().getNumElements();
          if (maximumAdd < 0)
            break;
        }
        delete[] which;
        delete[] violations;
        numberGlobalViolations_ += numberViolated;
      }
    }
    /*
          Generate new cuts (global and/or local) and/or apply heuristics.  If
          CglProbing is used, then it should be first as it can fix continuous
          variables.

          At present, CglProbing is the only case where generateCuts will return
          true. generateCuts actually modifies variable bounds in the solver when
          CglProbing indicates that it can fix a variable. Reoptimisation is required
          to take full advantage.

          The need to resolve here should only happen after a heuristic solution.
        optimalBasisIsAvailable resolves to basisIsAvailable, which seems to be part
        of the old OsiSimplex API. Doc'n says `Returns true if a basis is available
        and the problem is optimal. Should be used to see if the BinvARow type
        operations are possible and meaningful.' Which means any solver other the clp
        is probably doing a lot of unnecessary resolves right here.
          (Note default OSI implementation of optimalBasisIsAvailable always returns
          false.)
        */
    if (solverCharacteristics_->warmStart() && !solver_->optimalBasisIsAvailable()) {
      //printf("XXXXYY no opt basis\n");
#ifdef JJF_ZERO //def COIN_HAS_CLP
      //OsiClpSolverInterface * clpSolver
      //= dynamic_cast<OsiClpSolverInterface *> (solver_);
      int save = 0;
      if (clpSolver) {
        save = clpSolver->specialOptions();
        clpSolver->setSpecialOptions(save | 2048 /*4096*/); // Bonmin <something>
      }
#endif
      resolve(node ? node->nodeInfo() : NULL, 3);
#ifdef JJF_ZERO //def COIN_HAS_CLP
      if (clpSolver)
        clpSolver->setSpecialOptions(save);
#if CBC_USEFUL_PRINTING > 1
      if (clpSolver->getModelPtr()->numberIterations())
        printf("ITS %d pass %d\n",
          clpSolver->getModelPtr()->numberIterations(),
          currentPassNumber_);
#endif
#endif
    }
    if (nextRowCut_) {
      // branch was a cut - add it
      theseCuts.insert(*nextRowCut_);
      if (handler_->logLevel() > 1)
        nextRowCut_->print();
      const OsiRowCut *cut = nextRowCut_;
      double lb = cut->lb();
      double ub = cut->ub();
      int n = cut->row().getNumElements();
      const int *column = cut->row().getIndices();
      const double *element = cut->row().getElements();
      double sum = 0.0;
      for (int i = 0; i < n; i++) {
        int iColumn = column[i];
        double value = element[i];
        //if (cbcColSolution_[iColumn]>1.0e-7)
        //printf("value of %d is %g\n",iColumn,cbcColSolution_[iColumn]);
        sum += value * cbcColSolution_[iColumn];
      }
      delete nextRowCut_;
      nextRowCut_ = NULL;
      if (handler_->logLevel() > 1)
        printf("applying branch cut, sum is %g, bounds %g %g\n", sum, lb, ub);
      // possibly extend whichGenerator
      resizeWhichGenerator(numberViolated, numberViolated + 1);
      // set whichgenerator (also serves as marker to say don't delete0
      whichGenerator_[numberViolated++] = 20098;
    }

    // reset probing info
    //if (probingInfo_)
    //probingInfo_->initializeFixing();
    int i;
    // If necessary make cut generators work harder
    bool strongCuts = (!node && cut_obj[CUT_HISTORY - 1] != -COIN_DBL_MAX && fabs(cut_obj[CUT_HISTORY - 1] - cut_obj[CUT_HISTORY - 2]) < 1.0e-7 + 1.0e-6 * fabs(cut_obj[CUT_HISTORY - 1]));
    for (i = 0; i < numberCutGenerators_; i++)
      generator_[i]->setIneffectualCuts(strongCuts);
    // Print details
    if (!node) {
      handler_->message(CBC_ROOT_DETAIL, messages_)
        << currentPassNumber_
        << solver_->getNumRows()
        << solver_->getNumRows() - numberRowsAtContinuous_
        << solver_->getObjValue()
        << CoinMessageEol;
    }
    //Is Necessary for Bonmin? Always keepGoing if cuts have been generated in last iteration (taken from similar code in Cbc-2.4)
    if (solverCharacteristics_->solutionAddsCuts() && numberViolated) {
      for (i = 0; i < numberCutGenerators_; i++) {
        if (generator_[i]->mustCallAgain()) {
          keepGoing = true; // say must go round
          break;
        }
      }
    }
    if (!keepGoing) {
      // Status for single pass of cut generation
      int status = 0;
      /*
          threadMode with bit 2^1 set indicates we should use threads for root cut
          generation.
        */
      if ((threadMode_ & 2) == 0 || numberNodes_) {
        status = serialCuts(theseCuts, node, slackCuts, lastNumberCuts);
      } else {
        // do cuts independently
#ifdef CBC_THREAD
        status = parallelCuts(master, theseCuts, node, slackCuts, lastNumberCuts);
#endif
      }
      // Do we need feasible and violated?
      feasible = (status >= 0);
      if (status == 1)
        keepGoing = true;
      else if (status == 2)
        numberTries = 0;
      if (!feasible)
        violated = -2;
    }
    //if (!feasible)
    //break;
    /*
          End of the loop to exercise each generator - try heuristics
          - unless at root node and first pass
        */
    if ((numberNodes_ || currentPassNumber_ != 1) && (!this->maximumSecondsReached())) {
      double *newSolution = new double[numberColumns];
      double heuristicValue = getCutoff();
      int found = -1; // no solution found
      int whereFrom = numberNodes_ ? 4 : 1;
      for (i = 0; i < numberHeuristics_; i++) {
        // skip if can't run here
        if (!heuristic_[i]->shouldHeurRun(whereFrom))
          continue;
        // see if heuristic will do anything
        double saveValue = heuristicValue;
        int ifSol = heuristic_[i]->solution(heuristicValue,
          newSolution);
        //theseCuts) ;
        if (ifSol > 0) {
          // better solution found
          heuristic_[i]->incrementNumberSolutionsFound();
          found = i;
          incrementUsed(newSolution);
          lastHeuristic_ = heuristic_[found];
#ifdef HEURISTIC_INFORM
          printf("HEUR %s where %d A\n",
            lastHeuristic_->heuristicName(), whereFrom);
#endif
          // CBC_ROUNDING is symbolic; just says found by heuristic
          setBestSolution(CBC_ROUNDING, heuristicValue, newSolution);
          whereFrom |= 8; // say solution found
        } else if (ifSol < 0) {
          heuristicValue = saveValue;
        }
      }
      /*
              Did any of the heuristics turn up a new solution? Record it before we free
              the vector.
            */
      if (found >= 0) {
        phase_ = 4;
        CbcTreeLocal *tree
          = dynamic_cast< CbcTreeLocal * >(tree_);
        if (tree)
          tree->passInSolution(bestSolution_, heuristicValue);
      }
      delete[] newSolution;
    }
    CbcEventHandler *eventHandler = getEventHandler();
    if (eventHandler) {
      // Massage cuts??
      // save appData
      void *saveAppData = getApplicationData();
      // point to cuts
      setApplicationData(&theseCuts);
      eventHandler->event(CbcEventHandler::generatedCuts);
      setApplicationData(saveAppData);
    }
#ifdef JJF_ZERO
    // switch on to get all cuts printed
    theseCuts.printCuts();
#endif
    int numberColumnCuts = theseCuts.sizeColCuts();
    int numberRowCuts = theseCuts.sizeRowCuts();
    if (violated >= 0)
      violated = numberRowCuts + numberColumnCuts;
    /*
          Apply column cuts (aka bound tightening). This may be partially redundant
          for column cuts returned by CglProbing, as generateCuts installs bounds
          from CglProbing when it determines it can fix a variable.

          TODO: Looks like the use of violated has evolved. The value set above is
        	completely ignored. All that's left is violated == -1 indicates some
        	cut is violated, violated == -2 indicates infeasibility. Only
        	infeasibility warrants exceptional action.

          TODO: Strikes me that this code will fail to detect infeasibility, because
        	the breaks escape the inner loops but the outer loop keeps going.
        	Infeasibility in an early cut will be overwritten if a later cut is
        	merely violated.
        */
    if (numberColumnCuts) {

#ifdef CBC_DEBUG
      double *oldLower = new double[numberColumns];
      double *oldUpper = new double[numberColumns];
      memcpy(oldLower, cbcColLower_, numberColumns * sizeof(double));
      memcpy(oldUpper, cbcColUpper_, numberColumns * sizeof(double));
#endif

      double integerTolerance = getDblParam(CbcIntegerTolerance);
      for (int i = 0; i < numberColumnCuts; i++) {
        const OsiColCut *thisCut = theseCuts.colCutPtr(i);
        const CoinPackedVector &lbs = thisCut->lbs();
        const CoinPackedVector &ubs = thisCut->ubs();
        int j;
        int n;
        const int *which;
        const double *values;
        n = lbs.getNumElements();
        which = lbs.getIndices();
        values = lbs.getElements();
        for (j = 0; j < n; j++) {
          int iColumn = which[j];
          double value = cbcColSolution_[iColumn];
#if CBC_DEBUG > 1
          printf("%d %g %g %g %g\n", iColumn, oldLower[iColumn],
            cbcColSolution_[iColumn], oldUpper[iColumn], values[j]);
#endif
          solver_->setColLower(iColumn, values[j]);
          if (value < values[j] - integerTolerance)
            violated = -1;
          if (values[j] > cbcColUpper_[iColumn] + integerTolerance) {
            // infeasible
            violated = -2;
            break;
          }
        }
        n = ubs.getNumElements();
        which = ubs.getIndices();
        values = ubs.getElements();
        for (j = 0; j < n; j++) {
          int iColumn = which[j];
          double value = cbcColSolution_[iColumn];
#if CBC_DEBUG > 1
          printf("%d %g %g %g %g\n", iColumn, oldLower[iColumn],
            cbcColSolution_[iColumn], oldUpper[iColumn], values[j]);
#endif
          solver_->setColUpper(iColumn, values[j]);
          if (value > values[j] + integerTolerance)
            violated = -1;
          if (values[j] < cbcColLower_[iColumn] - integerTolerance) {
            // infeasible
            violated = -2;
            break;
          }
        }
      }
#ifdef CBC_DEBUG
      delete[] oldLower;
      delete[] oldUpper;
#endif
    }
    /*
          End installation of column cuts. The break here escapes the numberTries
          loop.
        */
    if (violated == -2 || !feasible) {
      // infeasible
      feasible = false;
      violated = -2;
      if (!numberNodes_)
        messageHandler()->message(CBC_INFEAS,
          messages())
          << CoinMessageEol;
      break;
    }
    /*
          Now apply the row (constraint) cuts. This is a bit more work because we need
          to obtain and augment the current basis.

          TODO: Why do this work, if there are no row cuts? The current basis will do
        	just fine.
        */
    int numberRowsNow = solver_->getNumRows();
#ifndef NDEBUG
    assert(numberRowsNow == numberRowsAtStart + lastNumberCuts);
#else
    // ? maybe clue to threaded problems
    if (numberRowsNow != numberRowsAtStart + lastNumberCuts) {
      fprintf(stderr, "*** threaded error - numberRowsNow(%d) != numberRowsAtStart(%d)+lastNumberCuts(%d)\n",
        numberRowsNow, numberRowsAtStart, lastNumberCuts);
      fprintf(stdout, "*** threaded error - numberRowsNow(%d) != numberRowsAtStart(%d)+lastNumberCuts(%d)\n",
        numberRowsNow, numberRowsAtStart, lastNumberCuts);
      abort();
    }
#endif
    int numberToAdd = theseCuts.sizeRowCuts();
    numberNewCuts_ = lastNumberCuts + numberToAdd;
    // resize whichGenerator
    resizeWhichGenerator(lastNumberCuts, numberNewCuts_);
    /*
          Now actually add the row cuts and reoptimise.

          Install the cuts in the solver using applyRowCuts and
          augment the basis with the corresponding slack. We also add each row cut to
          the set of row cuts (cuts.insert()) supplied as a parameter. The new basis
          must be set with setWarmStart().

          TODO: Seems to me the original code could allocate addCuts with size 0, if
        	numberRowCuts was 0 and numberColumnCuts was nonzero. That might
        	explain the memory fault noted in the comment by AJK.  Unfortunately,
        	just commenting out the delete[] results in massive memory leaks. Try
        	a revision to separate the row cut case. Why do we need addCuts at
        	all? A typing issue, apparently: OsiCut vs. OsiRowCut.

          TODO: It looks to me as if numberToAdd and numberRowCuts are identical at
        	this point. Confirm & get rid of one of them.

          TODO: Any reason why the three loops can't be consolidated?
        */
    const OsiRowCut **addCuts = NULL;
    if (numberRowCuts > 0 || numberColumnCuts > 0) {
      if (numberToAdd > 0) {
        int i;
        int *whichGenerator = whichGenerator_ + lastNumberCuts;
        // Faster to add all at once
        addCuts = new const OsiRowCut *[numberToAdd];
        for (i = 0; i < numberToAdd; i++) {
          addCuts[i] = &theseCuts.rowCut(i);
          whichGenerator[i] = 90;
        }
        if ((specialOptions_ & 262144) != 0 && !parentModel_) {
          //save
          for (i = 0; i < numberToAdd; i++)
            storedRowCuts_->addCut(*addCuts[i]);
        }
        solver_->applyRowCuts(numberToAdd, addCuts);
        CoinWarmStartBasis *basis = dynamic_cast< CoinWarmStartBasis * >(solver_->getWarmStart());
        assert(basis != NULL); // make sure not volume
        /* dylp bug

                  Consistent size used by OsiDylp as sanity check. Implicit resize seen
                  as an error. Hence this call to resize is necessary.
                */
        basis->resize(numberRowsAtStart + numberNewCuts_, numberColumns);
        for (i = 0; i < numberToAdd; i++) {
          basis->setArtifStatus(numberRowsNow + i,
            CoinWarmStartBasis::basic);
        }
        if (solver_->setWarmStart(basis) == false) {
          throw CoinError("Fail setWarmStart() after cut installation.",
            "solveWithCuts", "CbcModel");
        }
        delete basis;
      }
      //solver_->setHintParam(OsiDoDualInResolve,false,OsiHintTry);
      feasible = (resolve(node ? node->nodeInfo() : NULL, 2) != 0);
      //solver_->setHintParam(OsiDoDualInResolve,true,OsiHintTry);
      if (maximumSecondsReached()) {
        numberTries = -1000; // exit
        feasible = false;
        delete[] addCuts;
        break;
      }
#ifdef CBC_DEBUG
      printf("Obj value after cuts %g %d rows\n", solver_->getObjValue(),
        solver_->getNumRows());
      if (onOptimalPath && !solver_->isDualObjectiveLimitReached())
        assert(feasible);
#endif
    }
    /*
          No cuts. Cut short the cut generation (numberTries) loop.
        */
    else if (numberLastAttempts > 2 || experimentBreak < 2) {
      numberTries = 0;
    }
    /*
          If the problem is still feasible, first, call takeOffCuts() to remove cuts
          that are now slack. takeOffCuts() will call the solver to reoptimise if
          that's needed to restore a valid solution.

          Next, see if we should quit due to diminishing returns:
            * we've tried three rounds of cut generation and we're getting
              insufficient improvement in the objective; or
            * we generated no cuts; or
            * the solver declared optimality with 0 iterations after we added the
              cuts generated in this round.
          If we decide to keep going, prep for the next iteration.

          It just seems more safe to tell takeOffCuts() to call resolve(), even if
          we're not continuing cut generation. Otherwise code executed between here
          and final disposition of the node will need to be careful not to access the
          lp solution. It can happen that we lose feasibility in takeOffCuts ---
          numerical jitters when the cutoff bound is epsilon less than the current
          best, and we're evaluating an alternative optimum.

          TODO: After successive rounds of code motion, there seems no need to
        	distinguish between the two checks for aborting the cut generation
        	loop. Confirm and clean up.
        */
    if (feasible) {
      int cutIterations = solver_->getIterationCount();
      if (numberOldActiveCuts_ + numberNewCuts_
        && (numberNewCuts_ || doCutsNow(1))) {
        OsiCuts *saveCuts = node ? NULL : &slackCuts;
        int nDel = takeOffCuts(cuts, resolveAfterTakeOffCuts_, saveCuts, numberToAdd, addCuts);
        if (nDel)
          lastNumberCuts2_ = 0;
        if (solver_->isDualObjectiveLimitReached() && resolveAfterTakeOffCuts_) {
          feasible = false;
#ifdef CBC_DEBUG
          double z = solver_->getObjValue();
          double cut = getCutoff();
          printf("Lost feasibility by %g in takeOffCuts; z = %g, cutoff = %g\n",
            z - cut, z, cut);
#endif
        }
      }
      delete[] addCuts;
      if (feasible) {
        numberRowsAtStart = numberOldActiveCuts_ + numberRowsAtContinuous_;
        lastNumberCuts = numberNewCuts_;
        double thisObj = direction * solver_->getObjValue();
        bool badObj = (allowZeroIterations) ? thisObj < cut_obj[0] + minimumDrop
                                            : thisObj < cut_obj[CUT_HISTORY - 1] + minimumDrop;
#ifdef JJF_ZERO // probably not a good idea
        if (!badObj)
          numberLastAttempts = CoinMax(0, numberLastAttempts - 1);
#endif
        // Compute maximum number of bad passes
        if (minimumDrop > 0.0) {
          if (increaseDrop) {
            // slowly increase minimumDrop; breakpoints are rule-of-thumb
            if (currentPassNumber_ == 13)
              minimumDrop = CoinMax(1.5 * minimumDrop, 1.0e-5 * fabs(thisObj));
            else if (currentPassNumber_ > 20 && (currentPassNumber_ % 5) == 0)
              minimumDrop = CoinMax(1.1 * minimumDrop, 1.0e-5 * fabs(thisObj));
            else if (currentPassNumber_ > 50)
              minimumDrop = CoinMax(1.1 * minimumDrop, 1.0e-5 * fabs(thisObj));
          }
          int nBadPasses = 0;
          // The standard way of determining escape
          if (!experimentBreak) {
            double test = 0.01 * minimumDrop;
            double goodDrop = COIN_DBL_MAX;
            for (int j = CUT_HISTORY - 1; j >= 0; j--) {
              if (thisObj - cut_obj[j] < test) {
                nBadPasses++;
              } else {
                goodDrop = (thisObj - cut_obj[j]) / static_cast< double >(nBadPasses + 1);
                break;
              }
            }
            maximumBadPasses = CoinMax(maximumBadPasses, nBadPasses);
            if (nBadPasses < maximumBadPasses && goodDrop > minimumDrop)
              badObj = false; // carry on
          } else {
            // Experimental escape calculations
            //if (currentPassNumber_==13||currentPassNumber_>50)
            //minimumDrop = CoinMax(1.5*minimumDrop,1.0e-5*fabs(thisObj));
            double test = 0.1 * minimumDrop;
            double goodDrop = (thisObj - cut_obj[0]) / static_cast< double >(CUT_HISTORY);
            double objValue = thisObj;
            for (int j = CUT_HISTORY - 1; j >= 0; j--) {
              if (objValue - cut_obj[j] < test) {
                nBadPasses++;
                objValue = cut_obj[j];
              } else {
                break;
              }
            }
#if CBC_USEFUL_PRINTING > 12
            if (!parentModel_ && !numberNodes_)
              printf("badObj %s nBad %d maxBad %d goodDrop %g minDrop %g thisDrop %g obj %g\n",
                badObj ? "true" : "false",
                nBadPasses, maximumBadPasses, goodDrop, minimumDrop,
                thisObj - cut_obj[CUT_HISTORY - 1],
                solver_->getObjValue());
#endif
            maximumBadPasses = CoinMax(maximumBadPasses, nBadPasses);
            if (nBadPasses < 2 || goodDrop > 2.0 * minimumDrop) {
              if (experimentBreak <= 2 || goodDrop > 0.1 * minimumDrop)
                badObj = false; // carry on
            }
            if (experimentBreak > 1 && goodDrop < minimumDrop)
              numberLastAttempts++;
          }
        }
        // magic numbers, they seemed reasonable; there's a possibility here of going more than
        // nominal number of passes if we're doing really well.
        if (numberTries == 1 && currentDepth_ < 12 && currentPassNumber_ < 10) {
          double drop[12] = { 1.0, 2.0, 3.0, 10.0, 10.0, 10.0, 10.0, 20.0, 100.0, 100.0, 1000.0, 1000.0 };
          if (thisObj - lastObjective > drop[currentDepth_] * minimumDrop) {
            numberTries++;
#if CBC_USEFUL_PRINTING > 1
            //printf("drop %g %g %d\n",thisObj,lastObjective,currentPassNumber_);
#endif
          }
        }
        for (int j = 0; j < CUT_HISTORY - 1; j++)
          cut_obj[j] = cut_obj[j + 1];
        cut_obj[CUT_HISTORY - 1] = thisObj;
        bool allowEarlyTermination = currentPassNumber_ >= 10;
        if (currentDepth_ > 10 || (currentDepth_ > 5 && numberColumns > 200))
          allowEarlyTermination = true;
        //if (badObj && (currentPassNumber_ >= 10 || (currentDepth_>10))
        if (badObj && allowEarlyTermination
          //&&(currentPassNumber_>=10||lastObjective>firstObjective)
          && !keepGoing) {
          numberTries = 0;
        }
        if (numberRowCuts + numberColumnCuts == 0 || (cutIterations == 0 && !allowZeroIterations)) {
          // maybe give it one more try
          if (numberLastAttempts > 2 || currentDepth_ || experimentBreak < 2)
            numberTries = 0;
          else
            numberLastAttempts++;
        }
        if (numberTries > 0) {
          reducedCostFix();
          lastObjective = direction * solver_->getObjValue();
        }
      }
    } else {
      // not feasible
      delete[] addCuts;
    }
    /*
          We've lost feasibility --- this node won't be referencing the cuts we've
          been collecting, so decrement the reference counts.
        */
    if (!feasible) {
      int i;
      if (currentNumberCuts_) {
        lockThread();
        for (i = 0; i < currentNumberCuts_; i++) {
          // take off node
          if (addedCuts_[i]) {
            if (!addedCuts_[i]->decrement())
              delete addedCuts_[i];
            addedCuts_[i] = NULL;
          }
        }
        unlockThread();
      }
      numberTries = 0;
      keepGoing = false;
    }
    if (numberTries == 0 && feasible && !keepGoing && !parentModel_ && !numberNodes_) {
      for (int i = 0; i < numberCutGenerators_; i++) {
        if (generator_[i]->whetherCallAtEnd()
          && !generator_[i]->whetherInMustCallAgainMode()) {
          // give it some goes and switch off
          numberTries = (saveNumberTries + 4) / 5;
          generator_[i]->setWhetherCallAtEnd(false);
        }
      }
    }
  } while ( (numberTries > 0 || keepGoing) && (!this->maximumSecondsReached()) );
  /*
      End cut generation loop.
    */
  {
    // switch on
    for (int i = 0; i < numberCutGenerators_; i++)
      generator_[i]->setSwitchedOff(false);
  }
  //check feasibility.
  //If solution seems to be integer feasible calling setBestSolution
  //will eventually add extra global cuts which we need to install at
  //the nodes

  if (feasible && solverCharacteristics_->solutionAddsCuts()) { //check integer feasibility
    bool integerFeasible = true;
    const double *save = testSolution_;
    testSolution_ = solver_->getColSolution();
    // point to useful information
    OsiBranchingInformation usefulInfo = usefulInformation();
    for (int i = 0; i < numberObjects_ && integerFeasible; i++) {
      double infeasibility = object_[i]->checkInfeasibility(&usefulInfo);
      if (infeasibility)
        integerFeasible = false;
    }
    testSolution_ = save;
    // Consider the possibility that some alternatives here only make sense in context
    // of bonmin.
    if (integerFeasible) { //update
      double objValue = solver_->getObjValue();
      int numberGlobalBefore = globalCuts_.sizeRowCuts();
      // SOLUTION2 so won't up cutoff or print message
      setBestSolution(CBC_SOLUTION2, objValue,
        solver_->getColSolution(), 0);
      int numberGlobalAfter = globalCuts_.sizeRowCuts();
      int numberToAdd = numberGlobalAfter - numberGlobalBefore;
      if (numberToAdd > 0)
      //We have added some cuts say they are tight at that node
      //Basis and lp should already have been updated
      {
        feasible = (solver_->isProvenOptimal() && !solver_->isDualObjectiveLimitReached());
        if (feasible) {
          int numberCuts = numberNewCuts_ = cuts.sizeRowCuts();
          // possibly extend whichGenerator
          resizeWhichGenerator(numberCuts, numberToAdd + numberCuts);

          for (int i = numberGlobalBefore; i < numberGlobalAfter; i++) {
            whichGenerator_[numberNewCuts_++] = 20099;
#ifndef GLOBAL_CUTS_JUST_POINTERS
            cuts.insert(*globalCuts_.rowCutPtr(i));
#else
            OsiRowCut *rowCutPointer = globalCuts_.rowCutPtr(i);
            cuts.insert(rowCutPointer);
#endif
          }
          numberNewCuts_ = lastNumberCuts + numberToAdd;
          //now take off the cuts which are not tight anymore
          takeOffCuts(cuts, resolveAfterTakeOffCuts_, NULL);
          if (solver_->isDualObjectiveLimitReached() && resolveAfterTakeOffCuts_) {
            feasible = false;
          }
        }
        if (!feasible) { //node will be fathomed
          lockThread();
          for (int i = 0; i < currentNumberCuts_; i++) {
            // take off node
            if (addedCuts_[i]) {
              if (!addedCuts_[i]->decrement())
                delete addedCuts_[i];
              addedCuts_[i] = NULL;
            }
          }
          unlockThread();
        }
      }
    }
  }
  /*
    End of code block to check for a solution, when cuts may be added as a result
    of a feasible solution.

      Reduced cost fix at end. Must also check feasible, in case we've popped out
      because a generator indicated we're infeasible.
    */
  if (feasible && solver_->isProvenOptimal())
    reducedCostFix();
  // If at root node do heuristics
  if (!numberNodes_ && !maximumSecondsReached()) {
    // First see if any cuts are slack
    int numberRows = solver_->getNumRows();
    int numberAdded = numberRows - numberRowsAtContinuous_;
    if (numberAdded) {
      CoinWarmStartBasis *basis = dynamic_cast< CoinWarmStartBasis * >(solver_->getWarmStart());
      assert(basis != NULL);
      int *added = new int[numberAdded];
      int nDelete = 0;
      for (int j = numberRowsAtContinuous_; j < numberRows; j++) {
        if (basis->getArtifStatus(j) == CoinWarmStartBasis::basic) {
          //printf("%d slack!\n",j);
          added[nDelete++] = j;
        }
      }
      if (nDelete) {
        solver_->deleteRows(nDelete, added);
      }
      delete[] added;
      delete basis;
    }
    // mark so heuristics can tell
    int savePass = currentPassNumber_;
    currentPassNumber_ = 999999;
    double *newSolution = new double[numberColumns];
    double heuristicValue = getCutoff();
    int found = -1; // no solution found
    if (feasible) {
      int whereFrom = node ? 3 : 2;
      for (int i = 0; i < numberHeuristics_; i++) {
        // skip if can't run here
        if (!heuristic_[i]->shouldHeurRun(whereFrom))
          continue;
        // see if heuristic will do anything
        double saveValue = heuristicValue;
        int ifSol = heuristic_[i]->solution(heuristicValue,
          newSolution);
        if (ifSol > 0) {
          // better solution found
          heuristic_[i]->incrementNumberSolutionsFound();
          found = i;
          incrementUsed(newSolution);
          lastHeuristic_ = heuristic_[found];
#ifdef HEURISTIC_INFORM
          printf("HEUR %s where %d B\n",
            lastHeuristic_->heuristicName(), whereFrom);
#endif
          setBestSolution(CBC_ROUNDING, heuristicValue, newSolution);
          whereFrom |= 8; // say solution found
        } else {
          heuristicValue = saveValue;
        }
      }
    }
    currentPassNumber_ = savePass;
    if (found >= 0) {
      phase_ = 4;
    }
    delete[] newSolution;
  }
  // Up change due to cuts
  if (feasible)
    sumChangeObjective2_ += solver_->getObjValue() * solver_->getObjSense()
      - objectiveValue;
    //if ((numberNodes_%100)==0)
    //printf("XXb sum obj changed by %g\n",sumChangeObjective2_);
    /*
      End of cut generation loop.

      Now, consider if we want to disable or adjust the frequency of use for any
      of the cut generators. If the client specified a positive number for
      howOften, it will never change. If the original value was negative, it'll
      be converted to 1000000+|howOften|, and this value will be adjusted each
      time fullScan is true. Actual cut generation is performed every
      howOften%1000000 nodes; the 1000000 offset is just a convenient way to
      specify that the frequency is adjustable.

      During cut generation, we recorded the number of cuts produced by each
      generator for this node. For all cuts, whichGenerator records the generator
      that produced a cut.

      TODO: All this should probably be hidden in a method of the CbcCutGenerator
      class.
    lh:
    TODO: Can the loop that scans over whichGenerator to accumulate per
    generator counts be replaced by values in countRowCuts and
    countColumnCuts?

    << I think the answer is yes, but not the other way 'round. Row and
       column cuts are block interleaved in whichGenerator. >>

    The root is automatically a full scan interval. At the root, decide if
    we're going to do cuts in the tree, and whether we should keep the cuts we
    have.

    Codes for willBeCutsInTree:
    -1: no cuts in tree and currently active cuts seem ineffective; delete
    them
     0: no cuts in tree but currently active cuts seem effective; make them
    into architecturals (faster than treating them as cuts)
     1: cuts will be generated in the tree; currently active cuts remain as
    cuts
    -lh
    */
#ifdef NODE_LOG
  int fatherNum = (node == NULL) ? -1 : node->nodeNumber();
  double value = (node == NULL) ? -1 : node->branchingObject()->value();
  string bigOne = (solver_->getIterationCount() > 30) ? "*******" : "";
  string way = (node == NULL) ? "" : (node->branchingObject()->way()) == 1 ? "Down" : "Up";
  std::cout << "Node " << numberNodes_ << ", father " << fatherNum << ", #iterations " << solver_->getIterationCount() << ", sol value : " << solver_->getObjValue() << std::endl;
#endif
  if (fullScan && numberCutGenerators_) {
    /* If cuts just at root node then it will probably be faster to
           update matrix and leave all in */
    int willBeCutsInTree = 0;
    double thisObjective = solver_->getObjValue() * direction;
    // get sizes
    int numberRowsAdded = solver_->getNumRows() - numberRowsAtStart;
    CoinBigIndex numberElementsAdded = solver_->getNumElements() - numberElementsAtStart;
    double densityOld = static_cast< double >(numberElementsAtStart) / static_cast< double >(numberRowsAtStart);
    double densityNew = numberRowsAdded ? (static_cast< double >(numberElementsAdded)) / static_cast< double >(numberRowsAdded)
                                        : 0.0;
    /*
          If we're at the root, and we added cuts, and the cuts haven't changed the
          objective, and the cuts resulted in a significant increase (> 20%) in nonzero
          coefficients, do no cuts in the tree and ditch the current cuts. They're not
          cost-effective.
        */
    if (!numberNodes_) {
      if (!parentModel_) {
        //printf("%d global cuts\n",globalCuts_.sizeRowCuts()) ;
        if ((specialOptions_ & 1) != 0) {
          //specialOptions_ &= ~1;
          int numberCuts = globalCuts_.sizeRowCuts();
          const OsiRowCutDebugger *debugger = continuousSolver_->getRowCutDebugger();
          if (debugger) {
            for (int i = 0; i < numberCuts; i++) {
              OsiRowCut *cut = globalCuts_.rowCutPtr(i);
              if (debugger->invalidCut(*cut)) {
                continuousSolver_->applyRowCuts(1, cut);
                continuousSolver_->writeMps("bad");
                printf("BAD cut\n");
              }
              //CoinAssert (!debugger->invalidCut(*cut));
            }
          }
        }
      }
      //solver_->writeMps("second");
      if (numberRowsAdded)
        handler_->message(CBC_CUTS_STATS, messages_)
          << numberRowsAdded
          << densityNew
          << CoinMessageEol;
      if (thisObjective - startObjective < 1.0e-5 && numberElementsAdded > 0.2 * numberElementsAtStart)
        willBeCutsInTree = -1;
      int whenC = whenCuts_;
      if (whenC == 999999 || whenC == 999998) {
        int size = continuousSolver_->getNumRows() + continuousSolver_->getNumCols();
        bool smallProblem = size <= 550;
        smallProblem = false;
#if CBC_USEFUL_PRINTING > 1
        int maxPass = maximumCutPasses_;
#endif
        if (thisObjective - startObjective < 1.0e-5) {
          // No change in objective function
          if (numberElementsAdded > 0.2 * numberElementsAtStart) {
            if (whenCuts_ == 999999) {
              whenCuts_ = 5000010;
              if (!smallProblem)
                maximumCutPasses_ = CoinMax(maximumCutPasses_ >> 1, 1);
            } else if (whenCuts_ == 999998) {
              whenCuts_ = 5000010;
              if (!smallProblem)
                maximumCutPasses_ = CoinMax(maximumCutPasses_ >> 1, 1);
            }
#ifdef JJF_ZERO
          } else if (currentPassNumber_ < CoinMin(CoinAbs(maximumCutPassesAtRoot_), 8)) {
            if (whenCuts_ == 999999) {
              whenCuts_ = 8000008;
              maximumCutPasses_ = 1;
            } else if (whenCuts_ == 999998) {
              whenCuts_ = 10000008;
              maximumCutPasses_ = 1;
            }
          } else if (currentPassNumber_ < CoinMin(CoinAbs(maximumCutPassesAtRoot_), 50)) {
            if (whenCuts_ == 999999) {
              whenCuts_ = 8000008;
              maximumCutPasses_ = 1;
            } else if (whenCuts_ == 999998) {
              whenCuts_ = 10000006;
              maximumCutPasses_ = 1;
            }
          } else if (currentPassNumber_ < CoinAbs(maximumCutPassesAtRoot_)) {
            if (whenCuts_ == 999999) {
              whenCuts_ = 8000008;
              maximumCutPasses_ = 1;
            } else if (whenCuts_ == 999998) {
              whenCuts_ = 10000004;
              maximumCutPasses_ = 1;
            }
#endif
          } else {
            if (whenCuts_ == 999999) {
              whenCuts_ = 8000008;
              if (!smallProblem)
                maximumCutPasses_ = CoinMax(maximumCutPasses_ >> 1, 1);
            } else if (whenCuts_ == 999998) {
              whenCuts_ = 10000004;
              if (!smallProblem)
                maximumCutPasses_ = CoinMax(maximumCutPasses_ >> 1, 1);
            }
          }
        } else {
          // Objective changed
#ifdef JJF_ZERO
          if (currentPassNumber_ < CoinMin(CoinAbs(maximumCutPassesAtRoot_), 8)) {
            if (whenCuts_ == 999999) {
              whenCuts_ = 8000008;
              maximumCutPasses_ = 1;
            } else if (whenCuts_ == 999998) {
              whenCuts_ = 10000008;
              maximumCutPasses_ = 1;
            }
          } else if (currentPassNumber_ < CoinMin(CoinAbs(maximumCutPassesAtRoot_), 50)) {
            if (whenCuts_ == 999999) {
              whenCuts_ = 8000008;
              maximumCutPasses_ = 1;
            } else if (whenCuts_ == 999998) {
              whenCuts_ = 10000004;
              maximumCutPasses_ = 1;
            }
          } else
#endif
            if (currentPassNumber_ < CoinAbs(maximumCutPassesAtRoot_)) {
            if (whenCuts_ == 999999) {
              whenCuts_ = 8000008;
              if (!smallProblem)
                maximumCutPasses_ = CoinMax(maximumCutPasses_ >> 1, 1);
            } else if (whenCuts_ == 999998) {
              whenCuts_ = 10000004;
              if (!smallProblem)
                maximumCutPasses_ = CoinMax(maximumCutPasses_ >> 1, 1);
            }
          } else {
            if (whenCuts_ == 999999) {
              whenCuts_ = 10000004;
              maximumCutPasses_ = CoinMax(maximumCutPasses_, 2);
            } else if (whenCuts_ == 999998) {
              whenCuts_ = 11000002;
              maximumCutPasses_ = CoinMax(maximumCutPasses_, 2);
            }
          }
        }
        // Set bit to say don't try too hard if seems reasonable
        if (maximumCutPasses_ <= 5)
          whenCuts_ += 100000;
          //// end
#if CBC_USEFUL_PRINTING > 1
        printf("changing whenCuts from %d to %d and cutPasses from %d to %d objchange %g\n",
          whenC, whenCuts_, maxPass, maximumCutPasses_, thisObjective - startObjective);
#endif
      }
    }
    /*
          Noop block 071219.
        */
    if ((numberRowsAdded > 100 + 0.5 * numberRowsAtStart
          || numberElementsAdded > 0.5 * numberElementsAtStart)
      && (densityNew > 200.0 && numberRowsAdded > 100 && densityNew > 2.0 * densityOld)) {
      // much bigger
      //if (thisObjective-startObjective<0.1*fabs(startObjective)+1.0e-5)
      //willBeCutsInTree=-1;
      //printf("Cuts will be taken off , %d rows added with density %g\n",
      //     numberRowsAdded,densityNew);
    }
    /*
          Noop block 071219.
        */
    if (densityNew > 100.0 && numberRowsAdded > 2 && densityNew > 2.0 * densityOld) {
      //if (thisObjective-startObjective<0.1*fabs(startObjective)+1.0e-5)
      //willBeCutsInTree=-2;
      //printf("Density says no cuts ? , %d rows added with density %g\n",
      //     numberRowsAdded,densityNew);
    }
    // Root node or every so often - see what to turn off
    /*
          Hmmm ... > -90 for any generator will overrule previous decision to do no
          cuts in tree and delete existing cuts.
        */
    int i;
    for (i = 0; i < numberCutGenerators_; i++) {
      int howOften = generator_[i]->howOften();
      if (howOften > -90)
        willBeCutsInTree = 0;
    }
    if (!numberNodes_) {
      handler_->message(CBC_ROOT, messages_)
        << numberNewCuts_
        << startObjective << thisObjective
        << currentPassNumber_
        << CoinMessageEol;
      // do heuristics again! if feasibility pump still exists
      if ((specialOptions_ & 33554432) != 0 && !parentModel_) {
        specialOptions_ &= ~33554432;
        doHeuristicsAtRoot();
      }
    }
    /*
          Count the number of cuts produced by each cut generator on this call. Not
          clear to me that the accounting is equivalent here. whichGenerator_ records
          the generator for column and row cuts. So unless numberNewCuts is row cuts
          only, we're double counting for JUST_ACTIVE. Note too the multiplier applied
          to column cuts.
        */
    if (!numberNodes_) {
      double value = CoinMax(minimumDrop_, 0.005 * (thisObjective - startObjective) / static_cast< double >(currentPassNumber_));
      if (numberColumns < 200)
        value = CoinMax(minimumDrop_, 0.1 * value);
#if CBC_USEFUL_PRINTING > 1
      printf("Minimum drop for cuts was %g, now is %g\n", minimumDrop_, value);
#endif
      minimumDrop_ = value;
    }
    int *count = new int[numberCutGenerators_];
    memset(count, 0, numberCutGenerators_ * sizeof(int));
    int numberActiveGenerators = 0;
    for (i = 0; i < numberNewCuts_; i++) {
      int iGenerator = whichGenerator_[i];
      //assert (iGenerator>=0);
      if (iGenerator >= 0)
        iGenerator = iGenerator % 10000;
      if (iGenerator >= 0 && iGenerator < numberCutGenerators_)
        count[iGenerator]++;
    }
    // add in any active cuts if at root node (for multiple solvers)
#ifdef CHECK_KNOWN_SOLUTION
    if (onOptimalPath && (solver_->isDualObjectiveLimitReached() || !feasible)) {
      printf("help 2\n");
    }
#endif
    if (!numberNodes_) {
      for (i = 0; i < numberCutGenerators_; i++)
        count[i] += generator_[i]->numberCutsActive();
    }
    double totalCuts = 0.0;
    //#define JUST_ACTIVE
    for (i = 0; i < numberCutGenerators_; i++) {
      if (generator_[i]->numberCutsInTotal() || generator_[i]->numberColumnCuts())
        numberActiveGenerators++;
#ifdef JUST_ACTIVE
      double value = count[i];
#else
      double value = generator_[i]->numberCutsInTotal();
#endif
      totalCuts += value;
    }
    /*
          Open up a loop to step through the cut generators and decide what (if any)
          adjustment should be made for calling frequency.
        */
    int iProbing = -1;
    double smallProblem = (0.2 * totalCuts) / static_cast< double >(numberActiveGenerators + 1.0e-100);
    for (i = 0; i < numberCutGenerators_; i++) {
      int howOften = generator_[i]->howOften();
      /*  Probing can be set to just do column cuts in treee.
            But if doing good then leave as on
            Ok, let me try to explain this. rowCuts = 3 says do disaggregation (1<<0) and
            coefficient (1<<1) cuts. But if the value is negative, there's code at the
            entry to generateCuts, and generateCutsAndModify, that temporarily changes
            the value to 4 (1<<2) if we're in a search tree.

            Which does nothing to explain this next bit. We set a boolean, convert
            howOften to the code for `generate while objective is improving', and change
            over to `do everywhere'. Hmmm ... now I write it out, this makes sense in the
            context of the original comment. If we're doing well (objective improving)
            we'll keep probing fully active.

            */
      bool probingWasOnBut = false;
      CglProbing *probing = dynamic_cast< CglProbing * >(generator_[i]->generator());
      if (probing && !numberNodes_) {
        if (generator_[i]->numberCutsInTotal()) {
          // If large number of probing - can be biased
          smallProblem = (0.2 * (totalCuts - generator_[i]->numberCutsInTotal())) / static_cast< double >(numberActiveGenerators - 1 + 1.0e-100);
        }
        iProbing = i;
        if (probing->rowCuts() == -3) {
          probingWasOnBut = true;
          howOften = -98;
          probing->setRowCuts(3);
        }
      }
      /*
              Convert `as long as objective is improving' into `only at root' if we've
              decided cuts just aren't worth it.
            */
      if (willBeCutsInTree < 0 && howOften == -98)
        howOften = -99;
      /*
              And check to see if the objective is improving. But don't do the check if
              the user has specified some minimum number of cuts.

              This exclusion seems bogus, or at least counterintuitive. Why would a user
              suspect that setting a minimum cut limit would invalidate the objective
              check? Nor do I see the point in comparing the number of rows and columns
              in the second test.
            */
      if (!probing && howOften == -98 && !generator_[i]->numberShortCutsAtRoot() && generator_[i]->numberCutsInTotal()) {
        // switch off as no short cuts generated
        //printf("Switch off %s?\n",generator_[i]->cutGeneratorName());
        howOften = -99;
      }
      if (howOften == -98 && generator_[i]->switchOffIfLessThan() > 0) {
        if (thisObjective - startObjective < 0.005 * fabs(startObjective) + 1.0e-5)
          howOften = -99; // switch off
        if (thisObjective - startObjective < 0.1 * fabs(startObjective) + 1.0e-5
          && 5 * solver_->getNumRows() < solver_->getNumCols())
          howOften = -99; // switch off
      }
      if (generator_[i]->maximumTries() != -1)
        howOften = CoinMin(howOften, -99); // switch off
      /*
              Below -99, this generator is switched off. There's no need to consider
              further. Then again, there was no point in persisting this far!
            */
      if (howOften < -99) {
        // may have been switched off - report
        if (!numberNodes_) {
          int n = generator_[i]->numberCutsInTotal();
          if (n) {
            double average = 0.0;
            average = generator_[i]->numberElementsInTotal();
            average /= n;
            handler_->message(CBC_GENERATOR, messages_)
              << i
              << generator_[i]->cutGeneratorName()
              << n
              << average
              << generator_[i]->numberColumnCuts()
              << generator_[i]->numberCutsActive()
                + generator_[i]->numberColumnCuts();
            handler_->printing(generator_[i]->timing())
              << generator_[i]->timeInCutGenerator();
            handler_->message()
              << -100
              << CoinMessageEol;
          }
        }
        continue;
      }
      /*
              Adjust, if howOften is adjustable.
            */
      if (howOften < 0 || howOften >= 1000000) {
        if (!numberNodes_) {
          /*
                      If root only, or objective improvement but no cuts generated, switch off. If
                      it's just that the generator found no cuts at the root, give it one more
                      chance.
                    */
          // If small number switch mostly off
#ifdef JUST_ACTIVE
          double thisCuts = count[i] + 5.0 * generator_[i]->numberColumnCuts();
#else
          double thisCuts = generator_[i]->numberCutsInTotal() + 5.0 * generator_[i]->numberColumnCuts();
#endif
          // Allow on smaller number if <-1
          if (generator_[i]->switchOffIfLessThan() < 0) {
            double multiplier[] = { 2.0, 5.0 };
            int iSwitch = -generator_[i]->switchOffIfLessThan() - 1;
            assert(iSwitch >= 0 && iSwitch < 2);
            thisCuts *= multiplier[iSwitch];
          }
          if (!thisCuts || howOften == -99) {
            if (howOften == -99 || howOften == -98) {
              howOften = -100;
            } else {
              howOften = 1000000 + SCANCUTS; // wait until next time
              if (probing) {
                // not quite so drastic
                howOften = 1000000 + 1;
                probing->setMaxLook(1);
                probing->setMaxProbe(123);
              }
            }
            /*
                          Not productive, but not zero either.
                        */
          } else if ((thisCuts + generator_[i]->numberColumnCuts() < smallProblem)
            && !generator_[i]->whetherToUse()) {
            /*
                          Not unadjustable every node, and not strong probing.
                        */
            if (howOften != 1 && !probingWasOnBut) {
              /*
                              No depth spec, or not adjustable every node.
                            */
              if (generator_[i]->whatDepth() < 0 || howOften != -1) {
                int k = static_cast< int >(sqrt(smallProblem / thisCuts));
                /*
                                  Not objective improvement, set to new frequency, otherwise turn off.
                                */
                if (howOften != -98)
                  howOften = k + 1000000;
                else
                  howOften = -100;
                /*
                                  Depth spec, or adjustable every node. Force to unadjustable every node.
                                */
              } else {
                howOften = 1;
              }
              /*
                              Unadjustable every node, or strong probing. Force unadjustable every node and
                              force not strong probing? I don't understand.
                            */
            } else {
              howOften = 1;
              // allow cuts
              probingWasOnBut = false;
            }
            /*
                          Productive cut generator. Say we'll do it every node, adjustable. But if the
                          objective isn't improving, restrict that to every fifth depth level
                          (whatDepth overrides howOften in generateCuts).
                        */
          } else {
            if (thisObjective - startObjective < 0.1 * fabs(startObjective) + 1.0e-5 && generator_[i]->whatDepth() < 0)
              generator_[i]->setWhatDepth(5);
            howOften = 1 + 1000000;
          }
        }
        /*
                  End root actions.

                  sumChangeObjective2_ is the objective change due to cuts. If we're getting
                  much better results from branching over a large number of nodes, switch off
                  cuts.

                  Except it doesn't, really --- it just puts off the decision 'til the
                  next full scan, when it'll put it off again unless cuts look better.
                */
        // If cuts useless switch off
        if (numberNodes_ >= 100000 && sumChangeObjective1_ > 2.0e2 * (sumChangeObjective2_ + 1.0e-12)) {
          howOften = 1000000 + SCANCUTS; // wait until next time
          //printf("switch off cut %d due to lack of use\n",i);
        }
      }
      /*
              Ok, that's the frequency adjustment bit.

              Now, if we're at the root, force probing back on at every node, for column
              cuts at least, even if it looks useless for row cuts. Notice that if it
              looked useful, the values set above mean we'll be doing strong probing in
              the tree subject to objective improvement.
            */
      if (!numberNodes_) {
        if (probingWasOnBut && howOften == -100) {
          probing->setRowCuts(-3);
          howOften = 1;
        }
        if (howOften == 1)
          generator_[i]->setWhatDepth(1);

        if (howOften >= 0 && generator_[i]->generator()->mayGenerateRowCutsInTree())
          willBeCutsInTree = 1;
      }
      /*
              Set the new frequency in the generator. If this is an adjustable frequency,
              use the value to set whatDepth.

              Hey! Seems like this could override the user's depth setting.
            */
      generator_[i]->setHowOften(howOften);
      if (howOften >= 1000000 && howOften < 2000000 && 0) {
        // Go to depth
        int bias = 1;
        if (howOften == 1 + 1000000)
          generator_[i]->setWhatDepth(bias + 1);
        else if (howOften <= 10 + 1000000)
          generator_[i]->setWhatDepth(bias + 2);
        else
          generator_[i]->setWhatDepth(bias + 1000);
      }
      int newFrequency = generator_[i]->howOften() % 1000000;
      // increment cut counts
      generator_[i]->incrementNumberCutsActive(count[i]);
      CglStored *stored = dynamic_cast< CglStored * >(generator_[i]->generator());
      if (stored && !generator_[i]->numberCutsInTotal())
        continue;
      double average = 0.0;
      int n = generator_[i]->numberCutsInTotal();
      if (n) {
        average = generator_[i]->numberElementsInTotal();
        average /= n;
      }
      if (handler_->logLevel() > 1 || !numberNodes_) {
        handler_->message(CBC_GENERATOR, messages_)
          << i
          << generator_[i]->cutGeneratorName()
          //<<generator_[i]->numberCutsInTotal()<<count[i]
          << n
          << average
          << generator_[i]->numberColumnCuts()
          << generator_[i]->numberCutsActive()
            + generator_[i]->numberColumnCuts();
        handler_->printing(!numberNodes_ && generator_[i]->timing())
          << generator_[i]->timeInCutGenerator();
        handler_->message()
          << newFrequency
          << CoinMessageEol;
      }
    }
    /*
          End loop to adjust cut generator frequency of use.
        */
    delete[] count;
    if (!numberNodes_) {
      // save statistics
      for (i = 0; i < numberCutGenerators_; i++) {
        generator_[i]->setNumberCutsAtRoot(generator_[i]->numberCutsInTotal());
        generator_[i]->setNumberActiveCutsAtRoot(generator_[i]->numberCutsActive());
      }
      /*
              Garbage code 071219
            */
      // decide on pseudo cost strategy
      int howOften = iProbing >= 0 ? generator_[iProbing]->howOften() : 0;
      if ((howOften % 1000000) != 1)
        howOften = 0;
      //if (howOften) {
      //CglProbing * probing = dynamic_cast<CglProbing*>(generator_[iProbing]->generator());
      //}
      howOften = 0;
      if (howOften) {
        COIN_DETAIL_PRINT(printf("** method 1\n"));
        //CglProbing * probing = dynamic_cast<CglProbing*>(generator_[iProbing]->generator());
        generator_[iProbing]->setWhatDepth(1);
        // could set no row cuts
        //if (thisObjective-startObjective<0.001*fabs(startObjective)+1.0e-5)
        // probing->setRowCuts(0);
        for (int i = 0; i < numberObjects_; i++) {
          CbcSimpleIntegerDynamicPseudoCost *obj = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object_[i]);
          if (obj)
            obj->setMethod(1);
        }
      }
      if (willBeCutsInTree == -2)
        willBeCutsInTree = 0;
      /*
              End garbage code.

              Now I've reached the problem area. This is a problem only at the root node,
              so that should simplify the issue of finding a workable basis? Or maybe not.
            */
      if (willBeCutsInTree <= 0) {
        // Take off cuts
        cuts = OsiCuts();
        numberNewCuts_ = 0;
        if (!willBeCutsInTree) {
          // update size of problem
          numberRowsAtContinuous_ = solver_->getNumRows();
        } else {
          // take off cuts
          int numberRows = solver_->getNumRows();
          int numberAdded = numberRows - numberRowsAtContinuous_;
          if (numberAdded) {
            int *added = new int[numberAdded];
            for (int i = 0; i < numberAdded; i++)
              added[i] = i + numberRowsAtContinuous_;
            solver_->deleteRows(numberAdded, added);
            delete[] added;
            // resolve so optimal
            resolve(solver_);
          }
        }
#ifdef COIN_HAS_CLP
        OsiClpSolverInterface *clpSolver
          = dynamic_cast< OsiClpSolverInterface * >(solver_);
        if (clpSolver) {
          // Maybe solver might like to know only column bounds will change
          //int options = clpSolver->specialOptions();
          //clpSolver->setSpecialOptions(options|128);
          clpSolver->synchronizeModel();
        }
#endif
      } else {
#ifdef COIN_HAS_CLP
        OsiClpSolverInterface *clpSolver
          = dynamic_cast< OsiClpSolverInterface * >(solver_);
        if (clpSolver) {
          // make sure factorization can't carry over
          int options = clpSolver->specialOptions();
          clpSolver->setSpecialOptions(options & (~8));
        }
#endif
      }
    }
  } else {
#ifdef COIN_HAS_CLP
    OsiClpSolverInterface *clpSolver
      = dynamic_cast< OsiClpSolverInterface * >(solver_);
    if (clpSolver) {
      // Maybe solver might like to know only column bounds will change
      //int options = clpSolver->specialOptions();
      //clpSolver->setSpecialOptions(options|128);
      clpSolver->synchronizeModel();
    }
#endif
    if (numberCutGenerators_) {
      int i;
      // What if not feasible as cuts may have helped
      if (feasible) {
        for (i = 0; i < numberNewCuts_; i++) {
          int iGenerator = whichGenerator_[i];
#ifdef CONFLICT_CUTS
          assert(iGenerator >= 0);
#endif
          if (iGenerator >= 0)
            iGenerator = iGenerator % 10000;
          if (iGenerator >= 0 && iGenerator < numberCutGenerators_)
            generator_[iGenerator]->incrementNumberCutsActive();
        }
      }
    }
  }

#ifdef CHECK_CUT_COUNTS
  if (feasible) {
    CoinWarmStartBasis *basis = dynamic_cast< CoinWarmStartBasis * >(solver_->getWarmStart());
    printf("solveWithCuts: Number of rows at end (only active cuts) %d\n",
      numberRowsAtContinuous_ + numberNewCuts_ + numberOldActiveCuts_);
    basis->print();
    delete basis;
  }
#endif
#ifdef CHECK_KNOWN_SOLUTION
  if (onOptimalPath && (solver_->isDualObjectiveLimitReached() || !feasible)) {
    printf("help\n");
  }
#endif
#ifdef CBC_DEBUG
  if (onOptimalPath && !solver_->isDualObjectiveLimitReached())
    assert(feasible);
#endif
#ifdef COIN_HAS_CLP
  if (clpSolver)
    clpSolver->setSpecialOptions(saveClpOptions);
#endif
#ifdef CBC_THREAD
  // Get rid of all threaded stuff
  if (master) {
    master->stopThreads(0);
    delete master;
  }
#endif
  // make sure pointers are up to date
  setPointers(solver_);

  return feasible;
}

// Generate one round of cuts - serial mode
int CbcModel::serialCuts(OsiCuts &theseCuts, CbcNode *node, OsiCuts &slackCuts, int lastNumberCuts)
{
  /*
      Is it time to scan the cuts in order to remove redundant cuts? If so, set
      up to do it.
    */
  int fullScan = 0;
  if ((numberNodes_ % SCANCUTS) == 0 || (specialOptions_ & 256) != 0) {
    fullScan = 1;
    if (!numberNodes_ || (specialOptions_ & 256) != 0)
      fullScan = 2;
    specialOptions_ &= ~256; // mark as full scan done
  }
#if 0 //def COIN_HAS_CLP
    // check basis
    OsiClpSolverInterface * clpSolver
      = dynamic_cast<OsiClpSolverInterface *> (solver_);
    if (clpSolver) {
      ClpSimplex * simplex = clpSolver->getModelPtr();
      int numberTotal=simplex->numberRows()+simplex->numberColumns();
      int superbasic=0;
      for (int i=0;i<numberTotal;i++) {
	if (simplex->getStatus(i)==ClpSimplex::superBasic)
	  superbasic++;
      }
      if (superbasic) {
	printf("%d superbasic!\n",superbasic);
	clpSolver->resolve();
	superbasic=0;
	for (int i=0;i<numberTotal;i++) {
	  if (simplex->getStatus(i)==ClpSimplex::superBasic)
	    superbasic++;
	}
	assert (!superbasic);
      }
    }
#endif
  int switchOff = (!doCutsNow(1) && !fullScan) ? 1 : 0;
  int status = 0;
  int i;
  for (i = 0; i < numberCutGenerators_ && (!this->maximumSecondsReached()) ; i++) {
    int numberRowCutsBefore = theseCuts.sizeRowCuts();
    int numberColumnCutsBefore = theseCuts.sizeColCuts();
    int numberRowCutsAfter = numberRowCutsBefore;
    int numberColumnCutsAfter = numberColumnCutsBefore;
    /*printf("GEN %d %s switches %d\n",
	       i,generator_[i]->cutGeneratorName(),
	       generator_[i]->switches());*/
    bool generate = generator_[i]->normal();
    // skip if not optimal and should be (maybe a cut generator has fixed variables)
    if (generator_[i]->howOften() == -100 || (generator_[i]->needsOptimalBasis() && !solver_->basisIsAvailable())
      || generator_[i]->switchedOff())
      generate = false;
    if (switchOff && !generator_[i]->mustCallAgain()) {
      // switch off if default
      if (generator_[i]->howOften() == 1 && generator_[i]->whatDepth() < 0) {
        generate = false;
      } else if (currentDepth_ > -10 && switchOff == 2) {
        generate = false;
      }
    }
    if (generator_[i]->whetherCallAtEnd())
      generate = false;
    const OsiRowCutDebugger *debugger = NULL;
    bool onOptimalPath = false;
    if (generate) {
      bool mustResolve = generator_[i]->generateCuts(theseCuts, fullScan, solver_, node);
      numberRowCutsAfter = theseCuts.sizeRowCuts();
      if (fullScan && generator_[i]->howOften() == 1000000 + SCANCUTS_PROBING) {
        CglProbing *probing = dynamic_cast< CglProbing * >(generator_[i]->generator());
        if (probing && (numberRowCutsBefore < numberRowCutsAfter || numberColumnCutsBefore < theseCuts.sizeColCuts())) {
          // switch on
          generator_[i]->setHowOften(1);
        }
      }
      if (numberRowCutsBefore < numberRowCutsAfter && generator_[i]->mustCallAgain() && status >= 0)
        /*printf("%s before %d after %d must %c atend %c off %c endmode %c\n",
		   generator_[i]->cutGeneratorName(),
		   numberRowCutsBefore,numberRowCutsAfter,
		   generator_[i]->mustCallAgain() ? 'Y': 'N',
		   generator_[i]->whetherCallAtEnd() ? 'Y': 'N',
		   generator_[i]->switchedOff() ? 'Y': 'N',
		   generator_[i]->whetherInMustCallAgainMode() ? 'Y': 'N');*/
        if (numberRowCutsBefore < numberRowCutsAfter && generator_[i]->mustCallAgain() && status >= 0)
          status = 1; // say must go round
      // Check last cut to see if infeasible
      /*
            The convention is that if the generator proves infeasibility, it should
            return as its last cut something with lb > ub.
            */
      if (numberRowCutsBefore < numberRowCutsAfter) {
        const OsiRowCut *thisCut = theseCuts.rowCutPtr(numberRowCutsAfter - 1);
        if (thisCut->lb() > thisCut->ub()) {
          status = -1; // sub-problem is infeasible
          break;
        }
      }
#ifdef CBC_DEBUG
      {
        int k;
        for (k = numberRowCutsBefore; k < numberRowCutsAfter; k++) {
          OsiRowCut thisCut = theseCuts.rowCut(k);
          /* check size of elements.
                       We can allow smaller but this helps debug generators as it
                       is unsafe to have small elements */
          int n = thisCut.row().getNumElements();
          const int *column = thisCut.row().getIndices();
          const double *element = thisCut.row().getElements();
          //assert (n);
          for (int i = 0; i < n; i++) {
            double value = element[i];
            assert(fabs(value) > 1.0e-12 && fabs(value) < 1.0e20);
          }
        }
      }
#endif
      if (mustResolve /*|| (specialOptions_&1) != 0*/) {
        int returnCode = resolve(node ? node->nodeInfo() : NULL, 2);
        if (returnCode == 0)
          status = -1;
        if (returnCode < 0 && !status)
          status = 2;
        if ((specialOptions_ & 1) != 0) {
          debugger = solver_->getRowCutDebugger();
          if (debugger)
            onOptimalPath = (debugger->onOptimalPath(*solver_));
          else
            onOptimalPath = false;
          if (onOptimalPath && !solver_->isDualObjectiveLimitReached())
            assert(status >= 0);
        }
        if (status < 0)
          break;
      }
    }
    numberRowCutsAfter = theseCuts.sizeRowCuts();
    numberColumnCutsAfter = theseCuts.sizeColCuts();
    if ((specialOptions_ & 1) != 0) {
      if (onOptimalPath) {
        int k;
        for (k = numberRowCutsBefore; k < numberRowCutsAfter; k++) {
          OsiRowCut thisCut = theseCuts.rowCut(k);
          if (debugger->invalidCut(thisCut)) {
            solver_->getRowCutDebuggerAlways()->printOptimalSolution(*solver_);
            solver_->writeMpsNative("badCut.mps", NULL, NULL, 2);
            printf("Cut generator %d (%s) produced invalid cut (%dth in this go)\n",
              i, generator_[i]->cutGeneratorName(), k - numberRowCutsBefore);
            const double *lower = getColLower();
            const double *upper = getColUpper();
            int numberColumns = solver_->getNumCols();
            if (numberColumns < 200) {
              for (int i = 0; i < numberColumns; i++)
                printf("%d bounds %g,%g\n", i, lower[i], upper[i]);
            }
            abort();
          }
          assert(!debugger->invalidCut(thisCut));
        }
      }
    }
    /*
          The cut generator has done its thing, and maybe it generated some
          cuts.  Do a bit of bookkeeping: load
          whichGenerator[i] with the index of the generator responsible for a cut,
          and place cuts flagged as global in the global cut pool for the model.

          lastNumberCuts is the sum of cuts added in previous iterations; it's the
          offset to the proper starting position in whichGenerator.
        */
    int numberBefore = numberRowCutsBefore + lastNumberCuts;
    int numberAfter = numberRowCutsAfter + lastNumberCuts;
    // possibly extend whichGenerator
    resizeWhichGenerator(numberBefore, numberAfter);
    int j;

    /*
          Look for numerically unacceptable cuts.
        */
    bool dodgyCuts = false;
    for (j = numberRowCutsBefore; j < numberRowCutsAfter; j++) {
      const OsiRowCut *thisCut = theseCuts.rowCutPtr(j);
      if (thisCut->lb() > 1.0e10 || thisCut->ub() < -1.0e10) {
        dodgyCuts = true;
        break;
      }
      whichGenerator_[numberBefore++] = i + 20000;
      if (!numberNodes_ || generator_[i]->globalCuts())
        whichGenerator_[numberBefore - 1] = i + 10000;
      if (thisCut->lb() > thisCut->ub())
        status = -1; // sub-problem is infeasible
      if (thisCut->globallyValid() || !numberNodes_) {
        // add to global list
        OsiRowCut newCut(*thisCut);
        newCut.setGloballyValid(true);
        newCut.mutableRow().setTestForDuplicateIndex(false);
        globalCuts_.addCutIfNotDuplicate(newCut);
        whichGenerator_[numberBefore - 1] = i + 10000;
      }
    }
    if (dodgyCuts) {
      for (int k = numberRowCutsAfter - 1; k >= j; k--) {
        const OsiRowCut *thisCut = theseCuts.rowCutPtr(k);
        if (thisCut->lb() > thisCut->ub())
          status = -1; // sub-problem is infeasible
        if (thisCut->lb() > 1.0e10 || thisCut->ub() < -1.0e10)
          theseCuts.eraseRowCut(k);
      }
      numberRowCutsAfter = theseCuts.sizeRowCuts();
      for (; j < numberRowCutsAfter; j++) {
        const OsiRowCut *thisCut = theseCuts.rowCutPtr(j);
        whichGenerator_[numberBefore++] = i + 20000;
        if (!numberNodes_ || generator_[i]->globalCuts())
          whichGenerator_[numberBefore - 1] = i + 10000;
        if (thisCut->globallyValid()) {
          // add to global list
          OsiRowCut newCut(*thisCut);
          newCut.setGloballyValid(true);
          newCut.mutableRow().setTestForDuplicateIndex(false);
          globalCuts_.addCutIfNotDuplicate(newCut);
          whichGenerator_[numberBefore - 1] = i + 10000;
        }
      }
    }
    for (j = numberColumnCutsBefore; j < numberColumnCutsAfter; j++) {
      //whichGenerator_[numberBefore++] = i ;
      const OsiColCut *thisCut = theseCuts.colCutPtr(j);
      if (thisCut->globallyValid()) {
        // fix
        makeGlobalCut(thisCut);
      }
    }
  }
  /*
      End of loop to run each cut generator.
    */
  if (status >= 0) {
    // delete null cuts
    int nCuts = theseCuts.sizeRowCuts();
    int k;
    for (k = nCuts - 1; k >= 0; k--) {
      const OsiRowCut *thisCut = theseCuts.rowCutPtr(k);
      int n = thisCut->row().getNumElements();
      if (!n)
        theseCuts.eraseRowCut(k);
    }
  }
  // Add in any violated saved cuts
  if (!theseCuts.sizeRowCuts() && !theseCuts.sizeColCuts()) {
    int numberOld = theseCuts.sizeRowCuts() + lastNumberCuts;
    int numberCuts = slackCuts.sizeRowCuts();
    int i;
    // possibly extend whichGenerator
    resizeWhichGenerator(numberOld, numberOld + numberCuts);
    double primalTolerance;
    solver_->getDblParam(OsiPrimalTolerance, primalTolerance);
    for (i = 0; i < numberCuts; i++) {
      const OsiRowCut *thisCut = slackCuts.rowCutPtr(i);
      if (thisCut->violated(cbcColSolution_) > 100.0 * primalTolerance) {
        if (messageHandler()->logLevel() > 2)
          printf("Old cut added - violation %g\n",
            thisCut->violated(cbcColSolution_));
        whichGenerator_[numberOld++] = 20097;
        theseCuts.insert(*thisCut);
      }
    }
  }
  return status;
}

/*
  Remove slack cuts. We obtain a basis and scan it. Cuts with basic slacks
  are purged. If any cuts are purged, resolve() is called to restore the
  solution held in the solver.	If resolve() pivots, there's the possibility
  that a slack may be pivoted in (trust me :-), so the process iterates.
  Setting allowResolve to false will suppress reoptimisation (but see note
  below).

  At the level of the solver's constraint system, loose cuts are really
  deleted.  There's an implicit assumption that deleteRows will also update
  the active basis in the solver.

  At the level of nodes and models, it's more complicated.

  New cuts exist only in the collection of cuts passed as a parameter. They
  are deleted from the collection and that's the end of them.

  Older cuts have made it into addedCuts_. Two separate actions are needed.
  The reference count for the CbcCountRowCut object is decremented. If this
  count falls to 0, the node which owns the cut is located, the reference to
  the cut is removed, and then the cut object is destroyed (courtesy of the
  CbcCountRowCut destructor). We also need to set the addedCuts_ entry to
  NULL. This is important so that when it comes time to generate basis edits
  we can tell this cut was dropped from the basis during processing of the
  node.

  NOTE: In general, it's necessary to call resolve() after purging slack
	cuts.  Deleting constraints constitutes a change in the problem, and
	an OSI is not required to maintain a valid solution when the problem
	is changed. But ... it's really useful to maintain the active basis,
	and the OSI is supposed to do that. (Yes, it's splitting hairs.) In
	some places, it's possible to know that the solution will never be
	consulted after this call, only the basis.  (E.g., this routine is
	called as a last act before generating info to place the node in the
	live set.) For such use, set allowResolve to false.

  TODO: No real harm would be done if we just ignored the rare occasion when
	the call to resolve() pivoted a slack back into the basis. It's a
	minor inefficiency, at worst. But it does break assertions which
	check that there are no loose cuts in the basis. It might be better
	to remove the assertions.
*/

int CbcModel::takeOffCuts(OsiCuts &newCuts,
  bool allowResolve, OsiCuts *saveCuts,
  int numberNewCuts, const OsiRowCut **addedCuts)

{ // int resolveIterations = 0 ;
  int numberDropped = 0;
  int firstOldCut = numberRowsAtContinuous_;
  int totalNumberCuts = numberNewCuts_ + numberOldActiveCuts_;
  assert(numberRowsAtContinuous_ + totalNumberCuts == solver_->getNumRows());
  int *solverCutIndices = new int[totalNumberCuts];
  int *newCutIndices = new int[numberNewCuts_];
  const CoinWarmStartBasis *ws;
  CoinWarmStartBasis::Status status;
  //#define COIN_HAS_CLP_KEEP_STATUS
#ifdef COIN_HAS_CLP_KEEP_STATUS
  int problemStatus = -1;
  OsiClpSolverInterface *clpSolver
    = dynamic_cast< OsiClpSolverInterface * >(solver_);
  if (clpSolver)
    problemStatus = clpSolver->getModelPtr()->status();
#endif
  bool needPurge = true;
  /*
      The outer loop allows repetition of purge in the event that reoptimisation
      changes the basis. To start an iteration, clear the deletion counts and grab
      the current basis.
    */
  while (needPurge) {
    int numberNewToDelete = 0;
    int numberOldToDelete = 0;
    int i;
    int kCut = 0;
    ws = dynamic_cast< const CoinWarmStartBasis * >(solver_->getWarmStart());
    /*
          Scan the basis entries of the old cuts generated prior to this round of cut
          generation.  Loose cuts are `removed' by decrementing their reference count
          and setting the addedCuts_ entry to NULL. (If the reference count falls to
          0, they're really deleted.  See CbcModel and CbcCountRowCut doc'n for
          principles of cut handling.)
        */
    int oldCutIndex = 0;
    if (numberOldActiveCuts_) {
      lockThread();
      for (i = 0; i < numberOldActiveCuts_; i++) {
        status = ws->getArtifStatus(i + firstOldCut);
        while (!addedCuts_[oldCutIndex])
          oldCutIndex++;
        assert(oldCutIndex < currentNumberCuts_);
        // always leave if from nextRowCut_
        if (status == CoinWarmStartBasis::basic && (addedCuts_[oldCutIndex]->effectiveness() <= 1.0e10 || addedCuts_[oldCutIndex]->canDropCut(solver_, i + firstOldCut))) {
          solverCutIndices[numberOldToDelete++] = i + firstOldCut;
          if (saveCuts) {
            // send to cut pool
            OsiRowCut *slackCut = addedCuts_[oldCutIndex];
            if (slackCut->effectiveness() != -1.234) {
              slackCut->setEffectiveness(-1.234);
              slackCut->setGloballyValid();
              saveCuts->insert(*slackCut);
            }
          }
          if (addedCuts_[oldCutIndex]->decrement() == 0)
            delete addedCuts_[oldCutIndex];
          addedCuts_[oldCutIndex] = NULL;
          oldCutIndex++;
        } else {
          int iGenerator = addedCuts_[oldCutIndex]->whichCutGenerator();
          if (iGenerator == -1)
            iGenerator = 100;
          whichGenerator_[kCut++] = iGenerator;
          oldCutIndex++;
        }
      }
      unlockThread();
    }
    /*
          Scan the basis entries of the new cuts generated with this round of cut
          generation.  At this point, newCuts is the only record of the new cuts, so
          when we delete loose cuts from newCuts, they're really gone. newCuts is a
          vector, so it's most efficient to compress it (eraseRowCut) from back to
          front.
        */
    int firstNewCut = firstOldCut + numberOldActiveCuts_;
    int nCuts = newCuts.sizeRowCuts();
    for (i = 0; i < nCuts; i++) {
      status = ws->getArtifStatus(i + firstNewCut);
      if (status == CoinWarmStartBasis::basic &&
        /*whichGenerator_[i]!=-2*/ newCuts.rowCutPtr(i)->effectiveness() < 1.0e20) {
        solverCutIndices[numberNewToDelete + numberOldToDelete] = i + firstNewCut;
        newCutIndices[numberNewToDelete++] = i;
      } else { // save which generator did it
        // 20098 means branch cut! assert (whichGenerator_[i]!=20098); // ?? what if it is - memory leak?
        whichGenerator_[kCut++] = whichGenerator_[i];
      }
    }
    int baseRow = firstNewCut + nCuts;
    //OsiRowCut ** mutableAdded = const_cast<OsiRowCut **>(addedCuts);
    int numberTotalToDelete = numberNewToDelete + numberOldToDelete;
    for (i = 0; i < numberNewCuts; i++) {
      status = ws->getArtifStatus(i + baseRow);
      if (status != CoinWarmStartBasis::basic ||
        /*whichGenerator_[i+nCuts]==-2*/ addedCuts[i]->effectiveness() >= 1.0e20) {
        newCuts.insert(*addedCuts[i]);
        //newCuts.insert(mutableAdded[i]) ;
        //mutableAdded[i]=NULL;
        //if (status == CoinWarmStartBasis::basic&&whichGenerator_[i]!=-2) {
        // save which generator did it
        //whichGenerator_[k++] = whichGenerator_[i+nCuts] ;
        //}
      } else {
        solverCutIndices[numberTotalToDelete++] = i + baseRow;
      }
    }
    numberNewCuts = 0;
    numberNewCuts_ = newCuts.sizeRowCuts();
    delete ws;
    for (i = numberNewToDelete - 1; i >= 0; i--) {
      int iCut = newCutIndices[i];
      if (saveCuts) {
        // send to cut pool
        OsiRowCut *slackCut = newCuts.rowCutPtrAndZap(iCut);
        if (slackCut->effectiveness() != -1.234) {
          slackCut->setEffectiveness(-1.234);
          slackCut->setGloballyValid();
          saveCuts->insert(slackCut);
        } else {
          delete slackCut;
        }
      } else {
        newCuts.eraseRowCut(iCut);
      }
    }
    /*
          Did we delete anything? If so, delete the cuts from the constraint system
          held in the solver and reoptimise unless we're forbidden to do so. If the
          call to resolve() results in pivots, there's the possibility we again have
          basic slacks. Repeat the purging loop.
        */
    if (numberTotalToDelete > 0) {
      solver_->deleteRows(numberTotalToDelete,
        solverCutIndices);
      numberDropped += numberTotalToDelete;
      numberNewCuts_ -= numberNewToDelete;
      assert(numberNewCuts_ == newCuts.sizeRowCuts());
      numberOldActiveCuts_ -= numberOldToDelete;
#ifdef CBC_DEBUG
      printf("takeOffCuts: purged %d+%d cuts\n", numberOldToDelete,
        numberNewToDelete);
#endif
      if (allowResolve) {
        phase_ = 3;
        // can do quick optimality check
        int easy = 2;
        solver_->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, &easy);
        resolve(solver_);
        setPointers(solver_);
        solver_->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, NULL);
        if (solver_->getIterationCount() == 0) {
          needPurge = false;
        }
#ifdef CBC_DEBUG
        else {
          printf("Repeating purging loop. %d iters.\n",
            solver_->getIterationCount());
        }
#endif
      } else {
        needPurge = false;
      }
    } else {
      needPurge = false;
    }
  }

#ifdef COIN_HAS_CLP_KEEP_STATUS
  // need to check further that only zero duals dropped
  if (clpSolver) // status may have got to -1
    clpSolver->getModelPtr()->setProblemStatus(problemStatus);
#endif
  /*
      Clean up and return.
    */
  delete[] solverCutIndices;
  delete[] newCutIndices;
  return numberDropped;
}
/*
  Return values:
    1:	feasible
    0:	infeasible
   -1:	feasible and finished (do no more work on this subproblem)
*/
int CbcModel::resolve(CbcNodeInfo *parent, int whereFrom,
  double *saveSolution,
  double *saveLower,
  double *saveUpper)
{
#ifdef CBC_STATISTICS
  void cbc_resolve_check(const OsiSolverInterface *solver);
  cbc_resolve_check(solver_);
#endif
  bool onOptimalPath = false;
  if ((specialOptions_ & 1) != 0) {
    const OsiRowCutDebugger *debugger = solver_->getRowCutDebugger();
    if (debugger) {
      onOptimalPath = true;
      printf("On optimal path d\n");
    }
  }
  // We may have deliberately added in violated cuts - check to avoid message
  int iRow;
  int numberRows = solver_->getNumRows();
  const double *rowLower = solver_->getRowLower();
  const double *rowUpper = solver_->getRowUpper();
  bool feasible = true;
  for (iRow = numberRowsAtContinuous_; iRow < numberRows; iRow++) {
    if (rowLower[iRow] > rowUpper[iRow] + 1.0e-8)
      feasible = false;
  }
  // Can't happen if strong branching as would have been found before
  if ((!numberStrong_ || (moreSpecialOptions_ & 1073741824) != 0)
    && numberObjects_ > numberIntegers_) {
    int iColumn;
    int numberColumns = solver_->getNumCols();
    const double *columnLower = solver_->getColLower();
    const double *columnUpper = solver_->getColUpper();
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (columnLower[iColumn] > columnUpper[iColumn] + 1.0e-5)
        feasible = false;
    }
  }
#ifdef COIN_HAS_CLP
  OsiClpSolverInterface *clpSolver
    = dynamic_cast< OsiClpSolverInterface * >(solver_);
#endif
  /*
      Reoptimize. Consider the possibility that we should fathom on bounds. But be
      careful --- where the objective takes on integral values, we may want to keep
      a solution where the objective is right on the cutoff.
    */
  if (feasible) {
    int nTightened = 0;
#ifdef COIN_HAS_CLP
    // Pierre pointed out that this is not valid for all solvers
    // so just do if Clp
    if ((specialOptions_ & 1) != 0 && onOptimalPath) {
      solver_->writeMpsNative("before-tighten.mps", NULL, NULL, 2);
    }
    if (clpSolver && (!currentNode_ || (currentNode_->depth() & 2) != 0) && !solverCharacteristics_->solutionAddsCuts() && (moreSpecialOptions_ & 1073741824) == 0)
      nTightened = clpSolver->tightenBounds();
    if (nTightened) {
      //printf("%d bounds tightened\n",nTightened);
      if ((specialOptions_ & 1) != 0 && onOptimalPath) {
        const OsiRowCutDebugger *debugger = solver_->getRowCutDebugger();
        if (!debugger) {
          // tighten did something???
          solver_->getRowCutDebuggerAlways()->printOptimalSolution(*solver_);
          solver_->writeMpsNative("infeas4.mps", NULL, NULL, 2);
          printf("Not on optimalpath aaaa\n");
          //abort();
          onOptimalPath = false;
        }
      }
    }
#endif
    if (nTightened >= 0) {
      resolve(solver_);
      numberIterations_ += solver_->getIterationCount();
      feasible = (solver_->isProvenOptimal() && !solver_->isDualObjectiveLimitReached());
      if (feasible) {
        // double check
        double testValue = solver_->getObjSense() * solver_->getObjValue();
        //double cutoff = getCutoff();
        if (bestObjective_ - getCutoffIncrement() < testValue) {
#if CBC_USEFUL_PRINTING > 1
          double value;
          solver_->getDblParam(OsiDualObjectiveLimit, value);
          printf("Should cutoff as obj %.18g, best %.18g, inc %.18g - solver cutoff %.18g model cutoff %.18g\n",
            testValue, bestObjective_, getCutoffIncrement(),
            value, getCutoff());
#endif
          feasible = false;
        }
      } else if (solver_->isAbandoned()) {
        setMaximumSeconds(-COIN_DBL_MAX);
      }
#ifdef COIN_HAS_CLP
      if (clpSolver && feasible && !numberNodes_ && false) {
        double direction = solver_->getObjSense();
        double tolerance;
        solver_->getDblParam(OsiDualTolerance, tolerance);
        double primalTolerance;
        solver_->getDblParam(OsiPrimalTolerance, primalTolerance);

        const double *lower = solver_->getColLower();
        const double *upper = solver_->getColUpper();
        const double *solution = solver_->getColSolution();
        const double *reducedCost = solver_->getReducedCost();
        ClpSimplex *clpSimplex = clpSolver->getModelPtr();
        double *rowLower = clpSimplex->rowLower();
        double *rowUpper = clpSimplex->rowUpper();
        int numberRows = clpSimplex->numberRows();
        double *saveRowLower = CoinCopyOfArray(rowLower, numberRows);
        double *saveRowUpper = CoinCopyOfArray(rowUpper, numberRows);
        {
          const double *dual = clpSimplex->dualRowSolution();
          const double *rowActivity = clpSimplex->primalRowSolution();
          for (int iRow = 0; iRow < numberRows; iRow++) {
            double djValue = direction * dual[iRow];
            double lowerValue = rowLower[iRow];
            double upperValue = rowUpper[iRow];
            if (rowActivity[iRow] < lowerValue + primalTolerance && djValue > tolerance) {
              rowUpper[iRow] = lowerValue;
              assert(clpSimplex->getRowStatus(iRow) != ClpSimplex::basic);
            } else if (rowActivity[iRow] > upperValue - primalTolerance && djValue < -tolerance) {
              rowLower[iRow] = upperValue;
              assert(clpSimplex->getRowStatus(iRow) != ClpSimplex::basic);
            }
          }
        }
        int numberColumns = solver_->getNumCols();
        double *objective = clpSimplex->objective();
        double *saveObj = CoinCopyOfArray(objective, numberColumns);
        double objValue = 0.01;
        bool someFree = false;
        for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
          double djValue = direction * reducedCost[iColumn];
          double lowerValue = lower[iColumn];
          double upperValue = upper[iColumn];
          if (solution[iColumn] < lowerValue + primalTolerance && djValue > tolerance) {
            objective[iColumn] = 1.0e8 * direction;
            assert(clpSimplex->getColumnStatus(iColumn) != ClpSimplex::basic);
          } else if (solution[iColumn] > upperValue - primalTolerance && djValue < -tolerance) {
            objective[iColumn] = -1.0e8 * direction;
            assert(clpSimplex->getColumnStatus(iColumn) != ClpSimplex::basic);
          } else if (lowerValue > -1.0e20 || upperValue < 1.0e20) {
            assert(fabs(djValue) <= tolerance);
            if (fabs(lowerValue) < fabs(upperValue))
              objective[iColumn] = objValue * direction;
            else
              objective[iColumn] = -objValue * direction;
            objValue += 0.01;
          } else {
            objective[iColumn] = 0.0;
            someFree = true;
          }
        }
        if (!someFree)
          clpSimplex->primal(1);
        memcpy(objective, saveObj, numberColumns * sizeof(double));
        delete[] saveObj;
        memcpy(rowLower, saveRowLower, numberRows * sizeof(double));
        delete[] saveRowLower;
        memcpy(rowUpper, saveRowUpper, numberRows * sizeof(double));
        delete[] saveRowUpper;
        if (!someFree) {
          clpSimplex->primal(1);
          //assert (clpSimplex->numberIterations()<10);
        }
        //clpSimplex->writeMps("xx");
        //clpSimplex->primal(1);
        clpSolver->setWarmStart(NULL);
      }
#endif
      if ((specialOptions_ & 1) != 0 && onOptimalPath) {
        if (!solver_->getRowCutDebugger()) {
          // tighten did something???
          solver_->getRowCutDebuggerAlways()->printOptimalSolution(*solver_);
          solver_->writeMpsNative("infeas4.mps", NULL, NULL, 2);
          //assert (solver_->getRowCutDebugger()) ;
          printf("Not on optimalpath e\n");
          //abort();
        }
      }
    } else {
      feasible = false;
    }
  }
  if (0 && feasible) {
    const double *lb = solver_->getColLower();
    const double *ub = solver_->getColUpper();
    const double *x = solver_->getColSolution();
    const double *dj = solver_->getReducedCost();
    int numberColumns = solver_->getNumCols();
    for (int i = 0; i < numberColumns; i++) {
      if (dj[i] > 1.0e-4 && ub[i] - lb[i] > 1.0e-4 && x[i] > lb[i] + 1.0e-4)
        printf("error %d %g %g %g %g\n", i, dj[i], lb[i], x[i], ub[i]);
      if (dj[i] < -1.0e-4 && ub[i] - lb[i] > 1.0e-4 && x[i] < ub[i] - 1.0e-4)
        printf("error %d %g %g %g %g\n", i, dj[i], lb[i], x[i], ub[i]);
    }
  }
  if (false && !feasible && continuousObjective_ < -1.0e30) {
    // at root node - double double check
    bool saveTakeHint;
    OsiHintStrength saveStrength;
    solver_->getHintParam(OsiDoDualInResolve, saveTakeHint, saveStrength);
    if (saveTakeHint || saveStrength == OsiHintIgnore) {
      solver_->setHintParam(OsiDoDualInResolve, false, OsiHintDo);
      resolve(solver_);
      solver_->setHintParam(OsiDoDualInResolve, saveTakeHint, saveStrength);
      numberIterations_ += solver_->getIterationCount();
      feasible = solver_->isProvenOptimal();
      //      solver_->writeMps("infeas");
    }
  }
#ifdef JJF_ZERO
  if (cutModifier_ && feasible && !solverCharacteristics_->solutionAddsCuts()) {
    //double increment = getDblParam(CbcModel::CbcCutoffIncrement) ;
    double cutoff;
    solver_->getDblParam(OsiDualObjectiveLimit, cutoff);
    double distance = fabs(cutoff - solver_->getObjValue());
    if (distance < 10.0 * trueIncrement) {
      double offset;
      solver_->getDblParam(OsiObjOffset, offset);
      double objFixedValue = -offset;
      double objValue = 0.0;
      double direction = solver_->getObjSense();
      const double *solution = solver_->getColSolution();
      const double *objective = solver_->getObjCoefficients();
      const double *columnLower = solver_->getColLower();
      const double *columnUpper = solver_->getColUpper();
      int numberColumns = solver_->getNumCols();
      int increment = 0;
      double multiplier = 1.0 / trueIncrement;
      int bigIntegers = 0; // Count of large costs which are integer
      for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
        double value = solution[iColumn];
        // make sure clean
        value = CoinMin(value, columnUpper[iColumn]);
        value = CoinMax(value, columnLower[iColumn]);
        double cost = direction * objective[iColumn];
        if (cost) {
          if (columnLower[iColumn] < columnUpper[iColumn]) {
            objValue += value * cost;
            value = fabs(cost) * multiplier;
            if (value < 2.1e9) {
              int nearest = static_cast< int >(floor(value + 0.5));
              assert(fabs(value - floor(value + 0.5)) < 1.0e-8);
              if (!increment)
                increment = nearest;
              else
                increment = gcd(increment, nearest);
            } else {
              // large value - may still be multiple of 1.0
              value = fabs(objective[iColumn]);
              assert(fabs(value - floor(value + 0.5)) < 1.0e-8);
              bigIntegers++;
            }
          } else {
            // fixed
            objFixedValue += value * cost;
          }
        }
      }
      if (increment) {
        double value = increment;
        value /= multiplier;
        if (value > trueIncrement) {
          double x = objValue / value;
          x = ceil(x - 1.0e-5);
          x *= value;
          //printf("fixed %g, variable %g -> %g, sum %g - cutoff %g\n",
          // objFixedValue,objValue,x,x+objFixedValue,cutoff);
          x += objFixedValue;
          if (x > cutoff + 1.0e-5 * fabs(cutoff) + 1.0e-5) {
            //printf("Node cutoff\n");
            feasible = false;
          }
        } else {
          value = trueIncrement;
          double x = objValue / value;
          x = ceil(x - 1.0e-5);
          x *= value;
          x += objFixedValue;
          if (x > cutoff + 1.0e-5 * fabs(cutoff) + 1.0e-5) {
            //printf("Node cutoff\n");
            feasible = false;
          }
        }
      }
    }
  }
#endif

  setPointers(solver_);
  if (feasible && saveSolution) {
    // called from CbcNode
    assert(saveLower);
    assert(saveUpper);
    int numberColumns = solver_->getNumCols();
    memcpy(saveSolution, solver_->getColSolution(), numberColumns * sizeof(double));
    reserveCurrentSolution(saveSolution);
    memcpy(saveLower, solver_->getColLower(), numberColumns * sizeof(double));
    memcpy(saveUpper, solver_->getColUpper(), numberColumns * sizeof(double));
  }
#ifdef COIN_HAS_CLP
  if (clpSolver && !feasible) {
    // make sure marked infeasible
    if (!clpSolver->isProvenDualInfeasible())
      clpSolver->getModelPtr()->setProblemStatus(1);
  }
#endif
  int returnStatus = feasible ? 1 : 0;
  if (strategy_) {
    /*
          Possible returns from status:
            -1:	no recommendation
             0: treat as optimal
             1: treat as optimal and finished (no more resolves, cuts, etc.)
             2: treat as infeasible.
        */
    // user can play clever tricks here
    int status = strategy_->status(this, parent, whereFrom);
    if (status >= 0) {
      if (status == 0)
        returnStatus = 1;
      else if (status == 1)
        returnStatus = -1;
      else
        returnStatus = 0;
    }
  }
#if 0
    if ((specialOptions_&1) != 0 && onOptimalPath) {
      const OsiRowCutDebugger *debugger = solver_->getRowCutDebugger() ;
      if (!debugger) {
	// tighten did something???
	solver_->getRowCutDebuggerAlways()->printOptimalSolution(*solver_);
	solver_->writeMpsNative("infeas4.mps", NULL, NULL, 2);
	printf("Not on optimalpath aaaa\n");
	//abort();
      } else {
	printf("Still on optimal path\n");
      }
    }
#endif
  return returnStatus;
}

/* Set up objects.  Only do ones whose length is in range.
   If makeEquality true then a new model may be returned if
   modifications had to be made, otherwise "this" is returned.

   Could use Probing at continuous to extend objects
*/
CbcModel *
CbcModel::findCliques(bool makeEquality,
  int atLeastThisMany, int lessThanThis,
  int /*defaultValue*/)
{
  // No objects are allowed to exist
  assert(numberObjects_ == numberIntegers_ || !numberObjects_);
  CoinPackedMatrix matrixByRow(*solver_->getMatrixByRow());
  int numberRows = solver_->getNumRows();
  int numberColumns = solver_->getNumCols();

  // We may want to add columns
  int numberSlacks = 0;
  int *rows = new int[numberRows];
  double *element = new double[numberRows];

  int iRow;

  findIntegers(true);
  numberObjects_ = numberIntegers_;

  int numberCliques = 0;
  OsiObject **object = new OsiObject *[numberRows];
  int *which = new int[numberIntegers_];
  char *type = new char[numberIntegers_];
  int *lookup = new int[numberColumns];
  int i;
  for (i = 0; i < numberColumns; i++)
    lookup[i] = -1;
  for (i = 0; i < numberIntegers_; i++)
    lookup[integerVariable_[i]] = i;

  // Statistics
  int totalP1 = 0, totalM1 = 0;
  int numberBig = 0, totalBig = 0;
  int numberFixed = 0;

  // Row copy
  const double *elementByRow = matrixByRow.getElements();
  const int *column = matrixByRow.getIndices();
  const CoinBigIndex *rowStart = matrixByRow.getVectorStarts();
  const int *rowLength = matrixByRow.getVectorLengths();

  // Column lengths for slacks
  const int *columnLength = solver_->getMatrixByCol()->getVectorLengths();

  const double *lower = getColLower();
  const double *upper = getColUpper();
  const double *rowLower = getRowLower();
  const double *rowUpper = getRowUpper();
  /*
      Scan the rows, looking for individual rows that are clique constraints.
    */

  for (iRow = 0; iRow < numberRows; iRow++) {
    int numberP1 = 0, numberM1 = 0;
    CoinBigIndex j;
    double upperValue = rowUpper[iRow];
    double lowerValue = rowLower[iRow];
    bool good = true;
    int slack = -1;
    /*
          Does this row qualify? All variables must be binary and all coefficients
          +/- 1.0. Variables with positive coefficients are recorded at the low end of
          which, variables with negative coefficients the high end.
        */
    for (j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
      int iColumn = column[j];
      int iInteger = lookup[iColumn];
      if (upper[iColumn] - lower[iColumn] < 1.0e-8) {
        // fixed
        upperValue -= lower[iColumn] * elementByRow[j];
        lowerValue -= lower[iColumn] * elementByRow[j];
        continue;
      } else if (upper[iColumn] != 1.0 || lower[iColumn] != 0.0) {
        good = false;
        break;
      } else if (iInteger < 0) {
        good = false;
        break;
      } else {
        if (columnLength[iColumn] == 1)
          slack = iInteger;
      }
      if (fabs(elementByRow[j]) != 1.0) {
        good = false;
        break;
      } else if (elementByRow[j] > 0.0) {
        which[numberP1++] = iInteger;
      } else {
        numberM1++;
        which[numberIntegers_ - numberM1] = iInteger;
      }
    }
    int iUpper = static_cast< int >(floor(upperValue + 1.0e-5));
    int iLower = static_cast< int >(ceil(lowerValue - 1.0e-5));
    /*
          What do we have? If the row upper bound is greater than 1-numberM1, this
          isn't a clique. If the row upper bound is 1-numberM1, we have the classic
          clique (an SOS1 on binary variables, if numberM1 = 0). If the upper bound
          equals numberM1, we can fix all variables. If the upper bound is less than
          numberM1, we're infeasible.

          A similar analysis applies using numberP1 against the lower bound.
        */
    int state = 0;
    if (upperValue < 1.0e6) {
      if (iUpper == 1 - numberM1)
        state = 1;
      else if (iUpper == -numberM1)
        state = 2;
      else if (iUpper < -numberM1)
        state = 3;
    }
    if (!state && lowerValue > -1.0e6) {
      if (-iLower == 1 - numberP1)
        state = -1;
      else if (-iLower == -numberP1)
        state = -2;
      else if (-iLower < -numberP1)
        state = -3;
    }
    /*
        What to do? If we learned nothing, move on to the next iteration. If we're
        infeasible, we're outta here. If we decided we can fix variables, do it.
        */
    if (good && state) {
      if (abs(state) == 3) {
        // infeasible
        numberObjects_ = -1;
        break;
      } else if (abs(state) == 2) {
        // we can fix all
        numberFixed += numberP1 + numberM1;
        if (state > 0) {
          // fix all +1 at 0, -1 at 1
          for (i = 0; i < numberP1; i++)
            solver_->setColUpper(integerVariable_[which[i]], 0.0);
          for (i = 0; i < numberM1; i++)
            solver_->setColLower(integerVariable_[which[numberIntegers_ - i - 1]],
              1.0);
        } else {
          // fix all +1 at 1, -1 at 0
          for (i = 0; i < numberP1; i++)
            solver_->setColLower(integerVariable_[which[i]], 1.0);
          for (i = 0; i < numberM1; i++)
            solver_->setColUpper(integerVariable_[which[numberIntegers_ - i - 1]],
              0.0);
        }
      } else {
        /*
                  And the final case: we have a clique constraint. If it's within the allowed
                  size range, make a clique object.
                */
        int length = numberP1 + numberM1;
        if (length >= atLeastThisMany && length < lessThanThis) {
          // create object
          bool addOne = false;
          int objectType;
          /*
                      Choose equality (type 1) or inequality (type 0). If we're forcing equalities,
                      add a slack.
                    */
          if (iLower == iUpper) {
            objectType = 1;
          } else {
            if (makeEquality) {
              objectType = 1;
              element[numberSlacks] = state;
              rows[numberSlacks++] = iRow;
              addOne = true;
            } else {
              objectType = 0;
            }
          }
          /*
                      Record the strong values for the variables. Variables with positive
                      coefficients force all others when set to 1; variables with negative
                      coefficients force when set to 0. If the clique is formed against the row
                      lower bound, convert to the canonical form of a clique against the row
                      upper bound.
                    */
          if (state > 0) {
            totalP1 += numberP1;
            totalM1 += numberM1;
            for (i = 0; i < numberP1; i++)
              type[i] = 1;
            for (i = 0; i < numberM1; i++) {
              which[numberP1] = which[numberIntegers_ - i - 1];
              type[numberP1++] = 0;
            }
          } else {
            totalP1 += numberM1;
            totalM1 += numberP1;
            for (i = 0; i < numberP1; i++)
              type[i] = 0;
            for (i = 0; i < numberM1; i++) {
              which[numberP1] = which[numberIntegers_ - i - 1];
              type[numberP1++] = 1;
            }
          }
          if (addOne) {
            // add in slack
            which[numberP1] = numberIntegers_ + numberSlacks - 1;
            slack = numberP1;
            type[numberP1++] = 1;
          } else if (slack >= 0) {
            for (i = 0; i < numberP1; i++) {
              if (which[i] == slack) {
                slack = i;
              }
            }
          }
          object[numberCliques] = new CbcClique(this, objectType, numberP1,
            which, type,
            1000000 + numberCliques, slack);
          numberCliques++;
        } else if (numberP1 + numberM1 >= lessThanThis) {
          // too big
          numberBig++;
          totalBig += numberP1 + numberM1;
        }
      }
    }
  }
  delete[] which;
  delete[] type;
  delete[] lookup;
#if COIN_DEVELOP > 1
  if (numberCliques < 0) {
    printf("*** Problem infeasible\n");
  } else {
    if (numberCliques)
      printf("%d cliques of average size %g found, %d P1, %d M1\n",
        numberCliques,
        (static_cast< double >(totalP1 + totalM1)) / (static_cast< double >(numberCliques)),
        totalP1, totalM1);
    else
      printf("No cliques found\n");
    if (numberBig)
      printf("%d large cliques ( >= %d) found, total %d\n",
        numberBig, lessThanThis, totalBig);
    if (numberFixed)
      printf("%d variables fixed\n", numberFixed);
  }
#endif
  /*
      If required, augment the constraint matrix with clique slacks. Seems like we
      should be able to add the necessary integer objects without a complete
      rebuild of existing integer objects, but I'd need to look further to confirm
      that (lh, 071219). Finally, add the clique objects.
    */
  if (numberCliques > 0 && numberSlacks && makeEquality) {
    COIN_DETAIL_PRINT(printf("adding %d integer slacks\n", numberSlacks));
    // add variables to make equality rows
    int *temp = new int[numberIntegers_ + numberSlacks];
    memcpy(temp, integerVariable_, numberIntegers_ * sizeof(int));
    // Get new model
    CbcModel *newModel = new CbcModel(*this);
    OsiSolverInterface *newSolver = newModel->solver();
    for (i = 0; i < numberSlacks; i++) {
      temp[i + numberIntegers_] = i + numberColumns;
      int iRow = rows[i];
      double value = element[i];
      double lowerValue = 0.0;
      double upperValue = 1.0;
      double objValue = 0.0;
      CoinPackedVector column(1, &iRow, &value);
      newSolver->addCol(column, lowerValue, upperValue, objValue);
      // set integer
      newSolver->setInteger(numberColumns + i);
      if (value > 0)
        newSolver->setRowLower(iRow, rowUpper[iRow]);
      else
        newSolver->setRowUpper(iRow, rowLower[iRow]);
    }
    // replace list of integers
    for (i = 0; i < newModel->numberObjects_; i++)
      delete newModel->object_[i];
    newModel->numberObjects_ = 0;
    delete[] newModel->object_;
    newModel->object_ = NULL;
    newModel->findIntegers(true); //Set up all integer objects
    for (i = 0; i < numberIntegers_; i++) {
      newModel->modifiableObject(i)->setPriority(object_[i]->priority());
    }
    if (originalColumns_) {
      // old model had originalColumns
      delete[] newModel->originalColumns_;
      newModel->originalColumns_ = new int[numberColumns + numberSlacks];
      memcpy(newModel->originalColumns_, originalColumns_, numberColumns * sizeof(int));
      // mark as not in previous model
      for (i = numberColumns; i < numberColumns + numberSlacks; i++)
        newModel->originalColumns_[i] = -1;
    }
    delete[] rows;
    delete[] element;
    newModel->addObjects(numberCliques, object);
    assert(ownObjects_);
    for (; i < numberCliques; i++)
      delete object[i];
    delete[] object;
    newModel->synchronizeModel();
    return newModel;
  } else {
    assert(ownObjects_);
    if (numberCliques > 0) {
      addObjects(numberCliques, object);
      for (; i < numberCliques; i++)
        delete object[i];
      synchronizeModel();
    }
    delete[] object;
    delete[] rows;
    delete[] element;
    return this;
  }
}
// Fill in useful estimates
void CbcModel::pseudoShadow(int iActive)
{
  assert(iActive < 2 * 8 * 32 && iActive > -3);
  if (iActive == -1) {
    if (numberNodes_) {
      // zero out
      for (int i = 0; i < numberObjects_; i++) {
        CbcSimpleIntegerDynamicPseudoCost *obj1 = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object_[i]);
        if (obj1) {
          //assert (obj1->downShadowPrice()>0.0);
#define P_FACTOR 1.0
#ifndef JJF_ONE
          obj1->setDownShadowPrice(-P_FACTOR * obj1->downShadowPrice());
          obj1->setUpShadowPrice(-P_FACTOR * obj1->upShadowPrice());
#else
          double pCost;
          double sCost;
          pCost = obj1->downDynamicPseudoCost();
          sCost = P_FACTOR * obj1->downShadowPrice();
          if (!obj1->numberTimesDown() || sCost > pCost)
            obj1->updateDownDynamicPseudoCost(sCost);
          obj1->setDownShadowPrice(0.0);
          pCost = obj1->upDynamicPseudoCost();
          sCost = P_FACTOR * obj1->upShadowPrice();
          if (!obj1->numberTimesUp() || sCost > pCost)
            obj1->updateUpDynamicPseudoCost(sCost);
          obj1->setUpShadowPrice(0.0);
#endif
        }
      }
    }
    return;
  }
  bool doShadow = false;
  if (!iActive || iActive >= 32) {
    doShadow = true;
    if (iActive >= 32)
      iActive -= 32;
  }
  double *rowWeight = NULL;
  double *columnWeight = NULL;
  int numberColumns = solver_->getNumCols();
  int numberRows = solver_->getNumRows();
  // Column copy of matrix
  const double *element = solver_->getMatrixByCol()->getElements();
  const int *row = solver_->getMatrixByCol()->getIndices();
  const CoinBigIndex *columnStart = solver_->getMatrixByCol()->getVectorStarts();
  const int *columnLength = solver_->getMatrixByCol()->getVectorLengths();
  const double *dual = solver_->getRowPrice();
  const double *solution = solver_->getColSolution();
  const double *dj = solver_->getReducedCost();
  bool useMax = false;
  bool useAlpha = false;
  if (iActive) {
    // Use Patel and Chinneck ideas
    rowWeight = new double[numberRows];
    columnWeight = new double[numberColumns];
    // add in active constraints
    double tolerance = 1.0e-5;
    const double *rowLower = getRowLower();
    const double *rowUpper = getRowUpper();
    const double *rowActivity = solver_->getRowActivity();
    const double *lower = getColLower();
    const double *upper = getColUpper();
    CoinZeroN(rowWeight, numberRows);
    /* 1 A weight 1
           2 B weight 1/sum alpha
           3 L weight 1/number integer
           4 M weight 1/number active integer
           7 O weight 1/number integer and use alpha
           8 P weight 1/number active integer and use alpha
           9 up subtract 8 and use maximum
        */
    if (iActive > 8) {
      iActive -= 8;
      useMax = true;
    }
    if (iActive > 4) {
      iActive -= 4;
      useAlpha = true;
    }
    switch (iActive) {
      // A
    case 1:
      for (int iRow = 0; iRow < numberRows; iRow++) {
        if (rowActivity[iRow] > rowUpper[iRow] - tolerance || rowActivity[iRow] < rowLower[iRow] + tolerance) {
          rowWeight[iRow] = 1.0;
        }
      }
      break;
      // B
    case 2:
      for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
        if (upper[iColumn] > lower[iColumn]) {
          CoinBigIndex start = columnStart[iColumn];
          CoinBigIndex end = start + columnLength[iColumn];
          for (CoinBigIndex j = start; j < end; j++) {
            int iRow = row[j];
            rowWeight[iRow] += fabs(element[j]);
          }
        }
      }
      for (int iRow = 0; iRow < numberRows; iRow++) {
        if (rowWeight[iRow])
          rowWeight[iRow] = 1.0 / rowWeight[iRow];
      }
      break;
      // L
    case 3:
      for (int jColumn = 0; jColumn < numberIntegers_; jColumn++) {
        int iColumn = integerVariable_[jColumn];
        if (upper[iColumn] > lower[iColumn]) {
          CoinBigIndex start = columnStart[iColumn];
          CoinBigIndex end = start + columnLength[iColumn];
          for (CoinBigIndex j = start; j < end; j++) {
            int iRow = row[j];
            rowWeight[iRow]++;
          }
        }
      }
      for (int iRow = 0; iRow < numberRows; iRow++) {
        if (rowWeight[iRow])
          rowWeight[iRow] = 1.0 / rowWeight[iRow];
      }
      break;
      // M
    case 4:
      for (int jColumn = 0; jColumn < numberIntegers_; jColumn++) {
        int iColumn = integerVariable_[jColumn];
        double value = solution[iColumn];
        if (fabs(value - floor(value + 0.5)) > 1.0e-5) {
          CoinBigIndex start = columnStart[iColumn];
          CoinBigIndex end = start + columnLength[iColumn];
          for (CoinBigIndex j = start; j < end; j++) {
            int iRow = row[j];
            rowWeight[iRow]++;
          }
        }
      }
      for (int iRow = 0; iRow < numberRows; iRow++) {
        if (rowWeight[iRow])
          rowWeight[iRow] = 1.0 / rowWeight[iRow];
      }
      break;
    }
    if (doShadow) {
      for (int iRow = 0; iRow < numberRows; iRow++) {
        rowWeight[iRow] *= dual[iRow];
      }
    }
    dual = rowWeight;
  }
  const double *objective = solver_->getObjCoefficients();
  double direction = solver_->getObjSense();
  double *down = new double[numberColumns];
  double *up = new double[numberColumns];
  double upSum = 1.0e-20;
  double downSum = 1.0e-20;
  int numberIntegers = 0;
  if (doShadow) {
    // shadow prices
    if (!useMax) {
      for (int jColumn = 0; jColumn < numberIntegers_; jColumn++) {
        int iColumn = integerVariable_[jColumn];
        CoinBigIndex start = columnStart[iColumn];
        CoinBigIndex end = start + columnLength[iColumn];
        double upValue = 0.0;
        double downValue = 0.0;
        double value = direction * objective[iColumn];
        if (value) {
          if (value > 0.0)
            upValue += value;
          else
            downValue -= value;
        }
        for (CoinBigIndex j = start; j < end; j++) {
          int iRow = row[j];
          value = -dual[iRow];
          assert(fabs(dual[iRow]) < 1.0e50);
          if (value) {
            value *= element[j];
            if (value > 0.0)
              upValue += value;
            else
              downValue -= value;
          }
        }
        up[iColumn] = upValue;
        down[iColumn] = downValue;
        if (solver_->isInteger(iColumn)) {
          if (!numberNodes_ && handler_->logLevel() > 1)
            printf("%d - up %g down %g cost %g\n",
              iColumn, upValue, downValue, objective[iColumn]);
          upSum += upValue;
          downSum += downValue;
          numberIntegers++;
        }
      }
    } else {
      for (int jColumn = 0; jColumn < numberIntegers_; jColumn++) {
        int iColumn = integerVariable_[jColumn];
        CoinBigIndex start = columnStart[iColumn];
        CoinBigIndex end = start + columnLength[iColumn];
        double upValue = 0.0;
        double downValue = 0.0;
        double value = direction * objective[iColumn];
        if (value) {
          if (value > 0.0)
            upValue += value;
          else
            downValue -= value;
        }
        for (CoinBigIndex j = start; j < end; j++) {
          int iRow = row[j];
          value = -dual[iRow];
          if (value) {
            value *= element[j];
            if (value > 0.0)
              upValue = CoinMax(upValue, value);
            else
              downValue = CoinMax(downValue, -value);
          }
        }
        up[iColumn] = upValue;
        down[iColumn] = downValue;
        if (solver_->isInteger(iColumn)) {
          if (!numberNodes_ && handler_->logLevel() > 1)
            printf("%d - up %g down %g cost %g\n",
              iColumn, upValue, downValue, objective[iColumn]);
          upSum += upValue;
          downSum += downValue;
          numberIntegers++;
        }
      }
    }
  } else {
    for (int jColumn = 0; jColumn < numberIntegers_; jColumn++) {
      int iColumn = integerVariable_[jColumn];
      CoinBigIndex start = columnStart[iColumn];
      CoinBigIndex end = start + columnLength[iColumn];
      double upValue = 0.0;
      double downValue = 0.0;
      double value = direction * objective[iColumn];
      if (value) {
        if (value > 0.0)
          upValue += value;
        else
          downValue -= value;
      }
      double weight = 0.0;
      for (CoinBigIndex j = start; j < end; j++) {
        int iRow = row[j];
        value = -dual[iRow];
        double thisWeight = rowWeight[iRow];
        if (useAlpha)
          thisWeight *= fabs(element[j]);
        if (!useMax)
          weight += thisWeight;
        else
          weight = CoinMax(weight, thisWeight);
        if (value) {
          value *= element[j];
          if (value > 0.0)
            upValue += value;
          else
            downValue -= value;
        }
      }
      columnWeight[iColumn] = weight;
      // use dj if bigger
      double djValue = dj[iColumn];
      upValue = CoinMax(upValue, djValue);
      downValue = CoinMax(downValue, -djValue);
      up[iColumn] = upValue;
      down[iColumn] = downValue;
      if (solver_->isInteger(iColumn)) {
        if (!numberNodes_ && handler_->logLevel() > 1)
          printf("%d - dj %g up %g down %g cost %g\n",
            iColumn, djValue, upValue, downValue, objective[iColumn]);
        upSum += upValue;
        downSum += downValue;
        numberIntegers++;
      }
    }
    if (numberIntegers) {
      double averagePrice = (0.5 * (upSum + downSum)) / static_cast< double >(numberIntegers);
      //averagePrice *= 0.1;
      averagePrice *= 100.0;
      for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
        double weight = columnWeight[iColumn];
        up[iColumn] += averagePrice * weight;
        down[iColumn] += averagePrice * weight;
      }
    }
  }
  delete[] rowWeight;
  delete[] columnWeight;
  if (numberIntegers) {
    double smallDown = 0.0001 * (downSum / static_cast< double >(numberIntegers));
    double smallUp = 0.0001 * (upSum / static_cast< double >(numberIntegers));
#define PSEUDO_FACTOR 5.0e-1
    double pseudoFactor = PSEUDO_FACTOR;
    //if (!numberNodes_)
    //pseudoFactor=0.0;
    for (int i = 0; i < numberObjects_; i++) {
      CbcSimpleIntegerDynamicPseudoCost *obj1 = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object_[i]);
      if (obj1 && obj1->upShadowPrice() >= 0.0) {
        int iColumn = obj1->columnNumber();
        double upPseudoCost = obj1->upDynamicPseudoCost();
        double saveUp = upPseudoCost;
        upPseudoCost = CoinMax(pseudoFactor * upPseudoCost, smallUp);
        upPseudoCost = CoinMax(upPseudoCost, up[iColumn]);
        upPseudoCost = CoinMax(upPseudoCost, 0.001 * down[iColumn]);
        obj1->setUpShadowPrice(upPseudoCost);
        if (upPseudoCost > saveUp && !numberNodes_ && handler_->logLevel() > 1)
          printf("For %d up went from %g to %g\n",
            iColumn, saveUp, upPseudoCost);
        double downPseudoCost = obj1->downDynamicPseudoCost();
        double saveDown = downPseudoCost;
        downPseudoCost = CoinMax(pseudoFactor * downPseudoCost, smallDown);
        downPseudoCost = CoinMax(downPseudoCost, down[iColumn]);
        downPseudoCost = CoinMax(downPseudoCost, 0.001 * up[iColumn]);
        obj1->setDownShadowPrice(downPseudoCost);
        if (downPseudoCost > saveDown && !numberNodes_ && handler_->logLevel() > 1)
          printf("For %d down went from %g to %g\n",
            iColumn, saveDown, downPseudoCost);
      }
    }
  }
  delete[] down;
  delete[] up;
}

/*
  Set branching priorities.

  Setting integer priorities looks pretty robust; the call to findIntegers
  makes sure that SimpleInteger objects are in place. Setting priorities for
  other objects is entirely dependent on their existence, and the routine may
  quietly fail in several directions.
*/

void CbcModel::passInPriorities(const int *priorities,
  bool ifObject)
{
  findIntegers(false);
  int i;
  if (priorities) {
    int i0 = 0;
    int i1 = numberObjects_ - 1;
    if (ifObject) {
      for (i = numberIntegers_; i < numberObjects_; i++) {
        object_[i]->setPriority(priorities[i - numberIntegers_]);
      }
      i0 = numberIntegers_;
    } else {
      for (i = 0; i < numberIntegers_; i++) {
        object_[i]->setPriority(priorities[i]);
      }
      i1 = numberIntegers_ - 1;
    }
    messageHandler()->message(CBC_PRIORITY,
      messages())
      << i0 << i1 << numberObjects_ << CoinMessageEol;
  }
}

// Delete all object information
void CbcModel::deleteObjects(bool getIntegers)
{
  if (ownObjects_) {
    int i;
    for (i = 0; i < numberObjects_; i++)
      delete object_[i];
    delete[] object_;
  }
  object_ = NULL;
  numberObjects_ = 0;
  if (getIntegers && ownObjects_)
    findIntegers(true);
}

/*!
  Ensure all attached objects (OsiObjects, heuristics, and cut
  generators) point to this model.
*/
void CbcModel::synchronizeModel()
{
  if (!numberObjects_)
    return;
  int i;
  for (i = 0; i < numberHeuristics_; i++)
    heuristic_[i]->setModel(this);
  for (i = 0; i < numberObjects_; i++) {
    CbcObject *obj = dynamic_cast< CbcObject * >(object_[i]);
    if (obj) {
      obj->setModel(this);
      obj->setPosition(i);
    }
  }
  for (i = 0; i < numberCutGenerators_; i++)
    generator_[i]->refreshModel(this);

  if (!solverCharacteristics_) {
    OsiBabSolver *solverCharacteristics = dynamic_cast< OsiBabSolver * >(solver_->getAuxiliaryInfo());
    if (solverCharacteristics) {
      solverCharacteristics_ = solverCharacteristics;
    } else {
      // replace in solver
      OsiBabSolver defaultC;
      solver_->setAuxiliaryInfo(&defaultC);
      solverCharacteristics_ = dynamic_cast< OsiBabSolver * >(solver_->getAuxiliaryInfo());
    }
  }

  solverCharacteristics_->setSolver(solver_);
}

// Fill in integers and create objects

/**
  The routine first does a scan to count the number of integer variables.
  It then creates an array, integerVariable_, to store the indices of the
  integer variables, and an array of `objects', one for each variable.

  The scan is repeated, this time recording the index of each integer
  variable in integerVariable_, and creating an CbcSimpleInteger object that
  contains information about the integer variable. Initially, this is just
  the index and upper & lower bounds.

  \todo
  Note the assumption in cbc that the first numberIntegers_ objects are
  CbcSimpleInteger. In particular, the code which handles the startAgain
  case assumes that if the object_ array exists it can simply replace the first
  numberInteger_ objects. This is arguably unsafe.

  I am going to re-order if necessary
*/

void CbcModel::findIntegers(bool startAgain, int type)
{
  assert(solver_);
  /*
      No need to do this if we have previous information, unless forced to start
      over.
    */
  if (numberIntegers_ && !startAgain && object_)
    return;
  /*
      Clear out the old integer variable list, then count the number of integer
      variables.
    */
  delete[] integerVariable_;
  integerVariable_ = NULL;
  numberIntegers_ = 0;
  int numberColumns = getNumCols();
  int iColumn;
  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    if (isInteger(iColumn))
      numberIntegers_++;
  }
  // Find out how many old non-integer objects there are
  int nObjects = 0;
  OsiObject **oldObject = object_;
  int iObject;
  // also see where old ones were
  char *mark = new char[numberColumns];
  CoinZeroN(mark, numberColumns);
  int iPriority = -100000;
  for (iObject = 0; iObject < numberObjects_; iObject++) {
    iPriority = CoinMax(iPriority, object_[iObject]->priority());
    CbcSimpleInteger *obj = dynamic_cast< CbcSimpleInteger * >(oldObject[iObject]);
    if (obj) {
      int iColumn = obj->columnNumber();
      if (iColumn >= 0 && iColumn < numberColumns)
        mark[iColumn] = 1;
      delete oldObject[iObject];
    } else {
      oldObject[nObjects++] = oldObject[iObject];
    }
  }
  // See if there any SOS
#ifdef COIN_HAS_CLP
  if (!nObjects) {
    OsiClpSolverInterface *clpSolver
      = dynamic_cast< OsiClpSolverInterface * >(solver_);
    if (clpSolver && (clpSolver->numberSOS() || clpSolver->numberObjects())) {
      // deal with sos
      const CoinSet *setInfo = clpSolver->setInfo();
      int numberSOS = clpSolver->numberSOS();
      if (numberSOS) {
        nObjects = 0;
        delete[] oldObject;
        oldObject = new OsiObject *[numberSOS];
        for (int i = 0; i < numberSOS; i++) {
          int type = setInfo[i].setType();
          int n = setInfo[i].numberEntries();
          const int *which = setInfo[i].which();
          const double *weights = setInfo[i].weights();
          oldObject[nObjects++] = new CbcSOS(this, n, which, weights, i, type);
        }
      } else {
        // objects - only works with SOS at present
        int numberObjects = clpSolver->numberObjects();
        nObjects = 0;
        delete[] oldObject;
        oldObject = new OsiObject *[numberObjects];
        OsiObject **osiObjects = clpSolver->objects();
        for (int i = 0; i < numberObjects; i++) {
          OsiSOS *obj = dynamic_cast< OsiSOS * >(osiObjects[i]);
          if (obj) {
            int type = obj->setType();
            int n = obj->numberMembers();
            const int *which = obj->members();
            const double *weights = obj->weights();
            oldObject[nObjects++] = new CbcSOS(this, n, which, weights, i, type);
          }
        }
      }
    }
  }
#endif

  /*
      Found any? Allocate an array to hold the indices of the integer variables.
      Make a large enough array for all objects
    */
  delete[] integerVariable_;
  object_ = new OsiObject *[numberIntegers_ + nObjects];
  numberObjects_ = numberIntegers_ + nObjects;
  integerVariable_ = new int[numberIntegers_];
  /*
      Walk the variables again, filling in the indices and creating objects for
      the integer variables. Initially, the objects hold the index and upper &
      lower bounds.
    */
  numberIntegers_ = 0;
  if (type == 2)
    continuousPriority_ = iPriority + 1;
  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    if (isInteger(iColumn)) {
      if (!type) {
        object_[numberIntegers_] = new CbcSimpleInteger(this, iColumn);
      } else if (type == 1) {
        object_[numberIntegers_] = new CbcSimpleIntegerPseudoCost(this, iColumn, 0.3);
      } else if (type == 2) {
        object_[numberIntegers_] = new CbcSimpleInteger(this, iColumn);
        if (mark[iColumn]) {
          // could up priority on costs if all costs same??
        } else {
          object_[numberIntegers_]->setPriority(iPriority + 1);
        }
      }
      integerVariable_[numberIntegers_++] = iColumn;
    }
  }
  delete[] mark;
  // Now append other objects
  memcpy(object_ + numberIntegers_, oldObject, nObjects * sizeof(OsiObject *));
  // Delete old array (just array)
  delete[] oldObject;

  if (!numberObjects_)
    handler_->message(CBC_NOINT, messages_) << CoinMessageEol;
}
/* If numberBeforeTrust >0 then we are going to use CbcBranchDynamic.
   Scan and convert CbcSimpleInteger objects
*/
void CbcModel::convertToDynamic()
{
  int iObject;
  const double *cost = solver_->getObjCoefficients();
  bool allDynamic = true;
  for (iObject = 0; iObject < numberObjects_; iObject++) {
    CbcSimpleInteger *obj1 = dynamic_cast< CbcSimpleInteger * >(object_[iObject]);
    CbcSimpleIntegerPseudoCost *obj1a = dynamic_cast< CbcSimpleIntegerPseudoCost * >(object_[iObject]);
    CbcSimpleIntegerDynamicPseudoCost *obj2 = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object_[iObject]);
    if (obj1 && !obj2) {
      // replace
      int iColumn = obj1->columnNumber();
      int priority = obj1->priority();
      int preferredWay = obj1->preferredWay();
      double costValue = CoinMax(1.0e-5, fabs(cost[iColumn]));
      // treat as if will cost what it says up
      double upCost = costValue;
#ifndef BRANCH_BREAKEVEN
#define BRANCH_BREAKEVEN 0.3
#else
      preferredWay = 1;
#endif
      // and balance at breakeven of 0.3
      double downCost = ((1.0 - BRANCH_BREAKEVEN) * upCost) / BRANCH_BREAKEVEN;
      if (obj1a) {
        upCost = obj1a->upPseudoCost();
        downCost = obj1a->downPseudoCost();
      }
      delete object_[iObject];
      CbcSimpleIntegerDynamicPseudoCost *newObject = new CbcSimpleIntegerDynamicPseudoCost(this, iColumn, 1.0e0 * downCost, 1.0e0 * upCost);
      //newObject->setNumberBeforeTrust(numberBeforeTrust_);
      newObject->setPriority(priority);
      newObject->setPosition(iObject);
      newObject->setPreferredWay(preferredWay);
      object_[iObject] = newObject;
    } else if (!obj2) {
      CbcObject *obj3 = dynamic_cast< CbcObject * >(object_[iObject]);
      if (!obj3 || !obj3->optionalObject())
        allDynamic = false;
    } else {
      // synchronize trust
      //obj2->setNumberBeforeTrust(numberBeforeTrust_);
    }
  }
  if (branchingMethod_) {
    if ((branchingMethod_->whichMethod() & 1) == 0 && !branchingMethod_->chooseMethod()) {
      // Need a method which can do better
      delete branchingMethod_;
      branchingMethod_ = NULL;
    }
  }
  if (allDynamic)
    ownership_ |= 0x40000000;
  if (!branchingMethod_ && allDynamic) {
    // create one
    branchingMethod_ = new CbcBranchDynamicDecision();
  }
#ifdef SWITCH_VARIABLES
  // see if any switching variables
  if (numberIntegers_ < solver_->getNumCols())
    findSwitching();
#endif
  synchronizeNumberBeforeTrust();
}
#ifdef SWITCH_VARIABLES
// Convert Dynamic to Switching
int CbcModel::findSwitching()
{
  if ((moreSpecialOptions2_ & 1) == 0)
    return 0;
  const CoinPackedMatrix *rowCopy = solver_->getMatrixByRow();
  const int *column = rowCopy->getIndices();
  const int *rowLength = rowCopy->getVectorLengths();
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
  const double *rowLower = solver_->getRowLower();
  const double *rowUpper = solver_->getRowUpper();
  const double *columnLower = solver_->getColLower();
  const double *columnUpper = solver_->getColUpper();
  const double *element = rowCopy->getElements();
  //const double * element = solver_->getMatrixByCol()->getElements();
  const int *row = solver_->getMatrixByCol()->getIndices();
  const CoinBigIndex *columnStart = solver_->getMatrixByCol()->getVectorStarts();
  const int *columnLength = solver_->getMatrixByCol()->getVectorLengths();
  int numberRows = solver_->getNumRows();
  int numberColumns = solver_->getNumCols();
  int *sort = new int[2 * numberRows + 2 + numberColumns];
  int *whichRow = sort + numberRows + 1;
  int *marked = whichRow + numberRows + 1;
  memset(marked, 0, numberColumns * sizeof(int));
  int nnSwitch = 0;
  int nnSwitchTotal = 0;
  int n2Switch = 0;
  double largeRatio1 = 1000.0;
  double largeRatio2 = 100.0;
  for (int i = 0; i < numberIntegers_; i++) {
    int iColumn = integerVariable_[i];
    if (columnLower[iColumn] || columnUpper[iColumn] != 1.0)
      continue;
    if (!dynamic_cast< CbcSimpleInteger * >(object_[i]))
      continue;
    int nAdd = 0;
    bool takeThis = false;
    CoinBigIndex start = columnStart[iColumn];
    CoinBigIndex end = start + columnLength[iColumn];
    for (CoinBigIndex j = start; j < end; j++) {
      int iRow = row[j];
      if (rowLength[iRow] != 2) {
        continue;
      }
      // for now just 0.0 in rhs
      if (!rowLower[iRow]) {
        if (rowUpper[iRow] != COIN_DBL_MAX)
          continue;
      } else if (rowLower[iRow] != -COIN_DBL_MAX) {
        continue;
      } else if (rowUpper[iRow]) {
        continue;
      }
      CoinBigIndex k = rowStart[iRow];
      double bValue, cValue;
      int cColumn;
      if (column[k] == iColumn) {
        bValue = element[k];
        cValue = element[k + 1];
        cColumn = column[k + 1];
      } else {
        bValue = element[k + 1];
        cValue = element[k];
        cColumn = column[k];
      }
      if (solver_->isInteger(cColumn))
        continue;
      if (columnLower[cColumn] < 0.0)
        continue;
      if (bValue * cValue > 0.0)
        continue;
      if (fabs(bValue) > largeRatio1 * fabs(cValue))
        takeThis = true;
      // add to list
      whichRow[nAdd] = iRow;
      sort[nAdd++] = cColumn;
    }
    if (nAdd) {
      n2Switch++;
      CoinSort_2(sort, sort + nAdd, whichRow);
      int last = sort[0];
      for (int k = 1; k < nAdd; k++) {
        if (sort[k] == last)
          takeThis = true;
        else
          last = sort[k];
      }
      if (takeThis) {
        int last = sort[0];
        marked[last]++;
        for (int k = 1; k < nAdd; k++) {
          if (sort[k] != last) {
            last = sort[k];
            marked[last]++;
          }
        }
        //printf("Column %d has %d other columns\n",iColumn,nAdd);
        sort[nAdd] = COIN_INT_MAX;
        whichRow[nAdd] = COIN_INT_MAX;
        CbcSimpleIntegerDynamicPseudoCost *thisOne = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object_[i]);
        if (thisOne) {
          assert(iColumn == thisOne->columnNumber());
          object_[i] = new CbcSwitchingBinary(thisOne, nAdd, sort, whichRow);
          delete thisOne;
        } else {
          CbcSimpleInteger *thisOne = dynamic_cast< CbcSimpleInteger * >(object_[i]);
          assert(thisOne);
          assert(iColumn == thisOne->columnNumber());
          CbcSimpleIntegerDynamicPseudoCost tempObj(this, iColumn, 0.1);
          object_[i] = new CbcSwitchingBinary(&tempObj, nAdd, sort, whichRow);
          delete thisOne;
        }
      }
    }
    // see if there is an interesting row
    for (CoinBigIndex j = start; j < end; j++) {
      int iRow = row[j];
      // for now just 0.0 in rhs
      if (!rowLower[iRow]) {
        if (rowUpper[iRow] != COIN_DBL_MAX)
          continue;
      } else if (rowLower[iRow] != -COIN_DBL_MAX) {
        continue;
      } else if (rowUpper[iRow]) {
        continue;
      }
      int nOther = 0;
      double bEl = 0.0;
      double cMax = -COIN_DBL_MAX;
      double cMin = COIN_DBL_MAX;
      for (CoinBigIndex k = rowStart[iRow];
           k < rowStart[iRow] + rowLength[iRow]; k++) {
        int jColumn = column[k];
        if (jColumn == iColumn) {
          bEl = element[k];
        } else {
          sort[nOther++] = jColumn;
          if (solver_->isInteger(jColumn)) {
            cMin = -1.0;
            cMax = 1.0;
            break;
          } else {
            cMax = CoinMax(cMax, element[k]);
            cMin = CoinMin(cMin, element[k]);
            if (columnLower[jColumn] < 0.0) {
              cMin = -1.0;
              cMax = 1.0;
              break;
            }
          }
        }
      }
      double largestC = CoinMax(fabs(cMin), fabs(cMax));
      if (((cMin > 0.0 && bEl < 0.0 && !rowUpper[iRow]) || (cMin < 0.0 && bEl > 0.0 && !rowLower[iRow])) && cMin * cMax > 0.0 && fabs(bEl) > largeRatio2 * largestC) {
        // forces to zero
        CbcSwitchingBinary *object = dynamic_cast< CbcSwitchingBinary * >(object_[i]);
        if (!object) {
          // create empty one
          CbcSimpleIntegerDynamicPseudoCost *thisOne = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object_[i]);
          if (thisOne) {
            assert(iColumn == thisOne->columnNumber());
            object = new CbcSwitchingBinary(thisOne, 0, sort, whichRow);
            delete thisOne;
          } else {
            CbcSimpleInteger *thisOne = dynamic_cast< CbcSimpleInteger * >(object_[i]);
            assert(thisOne);
            assert(iColumn == thisOne->columnNumber());
            CbcSimpleIntegerDynamicPseudoCost tempObj(this, iColumn, 0.1);
            object = new CbcSwitchingBinary(&tempObj, 0, sort, whichRow);
            delete thisOne;
          }
          object_[i] = object;
        }
        object->addZeroSwitches(nOther, sort);
        nnSwitch++;
        nnSwitchTotal += nOther;
      }
    }
  }
  if (n2Switch + nnSwitch) {
    if (handler_->logLevel() > 2)
      printf("%d two switch variables - %d multi (total multi %d)\n",
        n2Switch, nnSwitch, nnSwitchTotal);
    memset(whichRow, 0, (numberRows + 1) * sizeof(int));
    for (int i = 0; i < numberColumns; i++) {
      whichRow[marked[i]]++;
    }
    if (handler_->logLevel() > 2) {
      for (int i = 0; i < numberRows + 1; i++) {
        if (whichRow[i])
          printf("%d variables have %d switches\n", whichRow[i], i);
      }
    }
  }
  delete[] sort;
  // say switches exist
  if (n2Switch + nnSwitch)
    moreSpecialOptions2_ |= 4;
  return n2Switch + nnSwitch;
}
// Fix associated variables
int CbcModel::fixAssociated(OsiSolverInterface *solver, int cleanBasis)
{
  int nChanged = 0;
  if ((moreSpecialOptions2_ & 4) != 0) {
    int n = -1;
    while (n) {
      n = 0;
      for (int i = 0; i < numberObjects_; i++) {
        CbcSwitchingBinary *object = dynamic_cast< CbcSwitchingBinary * >(object_[i]);
        if (object) {
          n += object->setAssociatedBounds(solver, cleanBasis);
        }
      }
      nChanged += n;
    }
  }
  return nChanged;
}
/* Debug associated variables
   printLevel - 1 summary if bad on fixed
                2 summary if bad on satisfied
                3 for individuals
 */
int CbcModel::checkAssociated(const OsiSolverInterface *solver,
  const double *solution, int printLevel)
{
  int nBad = 0;
  int nBadFixed = 0;
  if ((moreSpecialOptions2_ & 4) != 0) {
    int nAt0 = 0;
    int nAt1 = 0;
    int nBetween = 0;
    for (int i = 0; i < numberObjects_; i++) {
      CbcSwitchingBinary *object = dynamic_cast< CbcSwitchingBinary * >(object_[i]);
      if (object) {
        int state[3];
        nBad += object->checkAssociatedBounds(solver, solution, printLevel, state,
          nBadFixed);
        if (state[0] == 0)
          nBetween++;
        else if (state[0] == -1)
          nAt0++;
        else
          nAt1++;
      }
    }
    if (handler_->logLevel() > 2) {
      if (printLevel > 1 || (printLevel == 1 && nBadFixed)) {
        printf("%d switches, %d at 0, %d at 1, %d between - %d bad values (%d when fixed)\n",
          nBetween + nAt0 + nAt1, nAt0, nAt1, nBetween, nBad, nBadFixed);
        if (nBadFixed && printLevel != 3)
          checkAssociated(solver, solution, 3);
      }
    }
  }
  return nBad;
}
#endif
// Set numberBeforeTrust in all objects
void CbcModel::synchronizeNumberBeforeTrust(int type)
{
  int iObject;
  for (iObject = 0; iObject < numberObjects_; iObject++) {
    CbcSimpleIntegerDynamicPseudoCost *obj2 = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object_[iObject]);
    if (obj2) {
      // synchronize trust
      if (!type) {
        obj2->setNumberBeforeTrust(numberBeforeTrust_);
      } else if (type == 1) {
        int value = obj2->numberBeforeTrust();
        value = (value * 11) / 10 + 1;
        value = CoinMax(numberBeforeTrust_, value);
        obj2->setNumberBeforeTrust(value);
      } else {
        assert(type == 2);
        int value = obj2->numberBeforeTrust();
        int n = CoinMax(obj2->numberTimesDown(),
          obj2->numberTimesUp());
        if (n >= value) {
          value = CoinMin(CoinMin(n + 1, 3 * (value + 1) / 2), 5 * numberBeforeTrust_);
          obj2->setNumberBeforeTrust(value);
        }
      }
    }
  }
}

/* Add in any object information (objects are cloned - owner can delete
   originals */
void CbcModel::addObjects(int numberObjects, CbcObject **objects)
{
  // If integers but not enough objects fudge
  if (numberIntegers_ > numberObjects_ || !numberObjects_)
    findIntegers(true);
  /* But if incoming objects inherit from simple integer we just want
       to replace */
  int numberColumns = solver_->getNumCols();
  /** mark is -1 if not integer, >=0 if using existing simple integer and
        >=numberColumns if using new integer */
  int *mark = new int[numberColumns];
  int i;
  for (i = 0; i < numberColumns; i++)
    mark[i] = -1;
  int newNumberObjects = numberObjects;
  int newIntegers = 0;
  for (i = 0; i < numberObjects; i++) {
    CbcSimpleInteger *obj = dynamic_cast< CbcSimpleInteger * >(objects[i]);
    if (obj) {
      int iColumn = obj->columnNumber();
      assert(iColumn >= 0);
      mark[iColumn] = i + numberColumns;
      newIntegers++;
    }
  }
  // and existing
  for (i = 0; i < numberObjects_; i++) {
    CbcSimpleInteger *obj = dynamic_cast< CbcSimpleInteger * >(object_[i]);
    if (obj) {
      int iColumn = obj->columnNumber();
      if (mark[iColumn] < 0) {
        newIntegers++;
        newNumberObjects++;
        mark[iColumn] = i;
      }
    } else {
      // some other object - keep
      newNumberObjects++;
    }
  }
  delete[] integerVariable_;
  integerVariable_ = NULL;
#if COIN_DEVELOP > 1
  if (newIntegers != numberIntegers_)
    printf("changing number of integers from %d to %d\n",
      numberIntegers_, newIntegers);
#endif
  numberIntegers_ = newIntegers;
  integerVariable_ = new int[numberIntegers_];
  OsiObject **temp = new OsiObject *[newNumberObjects];
  // Put integers first
  newIntegers = 0;
  numberIntegers_ = 0;
  for (i = 0; i < numberColumns; i++) {
    int which = mark[i];
    if (which >= 0) {
      if (!isInteger(i)) {
        newIntegers++;
        solver_->setInteger(i);
      }
      if (which < numberColumns) {
        temp[numberIntegers_] = object_[which];
        object_[which] = NULL;
      } else {
        temp[numberIntegers_] = objects[which - numberColumns]->clone();
      }
      integerVariable_[numberIntegers_++] = i;
    }
  }
#if COIN_DEVELOP > 1
  if (newIntegers)
    printf("%d variables were declared integer\n", newIntegers);
#endif
  int n = numberIntegers_;
  // Now rest of old
  for (i = 0; i < numberObjects_; i++) {
    if (object_[i]) {
      CbcSimpleInteger *obj = dynamic_cast< CbcSimpleInteger * >(object_[i]);
      if (obj) {
        delete object_[i];
      } else {
        temp[n++] = object_[i];
      }
    }
  }
  // and rest of new
  for (i = 0; i < numberObjects; i++) {
    CbcSimpleInteger *obj = dynamic_cast< CbcSimpleInteger * >(objects[i]);
    if (!obj) {
      temp[n] = objects[i]->clone();
      CbcObject *obj = dynamic_cast< CbcObject * >(temp[n]);
      if (obj)
        obj->setModel(this);
      n++;
    }
  }
  delete[] mark;
  assert(ownObjects_);
  delete[] object_;
  object_ = temp;
  assert(n == newNumberObjects);
  numberObjects_ = newNumberObjects;
}
/* Add in any object information (objects are cloned - owner can delete
   originals */
void CbcModel::addObjects(int numberObjects, OsiObject **objects)
{
  // If integers but not enough objects fudge
  if (numberIntegers_ > numberObjects_)
    findIntegers(true);
  /* But if incoming objects inherit from simple integer we just want
       to replace */
  int numberColumns = solver_->getNumCols();
  /** mark is -1 if not integer, >=0 if using existing simple integer and
        >=numberColumns if using new integer */
  int *mark = new int[numberColumns];
  int i;
  for (i = 0; i < numberColumns; i++)
    mark[i] = -1;
  int newNumberObjects = numberObjects;
  int newIntegers = 0;
  for (i = 0; i < numberObjects; i++) {
    CbcSimpleInteger *obj = dynamic_cast< CbcSimpleInteger * >(objects[i]);
    if (obj) {
      int iColumn = obj->columnNumber();
      mark[iColumn] = i + numberColumns;
      newIntegers++;
    } else {
      OsiSimpleInteger *obj2 = dynamic_cast< OsiSimpleInteger * >(objects[i]);
      if (obj2) {
        // Osi takes precedence
        int iColumn = obj2->columnNumber();
        mark[iColumn] = i + numberColumns;
        newIntegers++;
      }
    }
  }
  // and existing
  for (i = 0; i < numberObjects_; i++) {
    CbcSimpleInteger *obj = dynamic_cast< CbcSimpleInteger * >(object_[i]);
    if (obj) {
      int iColumn = obj->columnNumber();
      if (mark[iColumn] < 0) {
        newIntegers++;
        newNumberObjects++;
        mark[iColumn] = i;
      }
    } else {
      newNumberObjects++;
    }
  }
  delete[] integerVariable_;
  integerVariable_ = NULL;
#if COIN_DEVELOP > 1
  if (newIntegers != numberIntegers_)
    printf("changing number of integers from %d to %d\n",
      numberIntegers_, newIntegers);
#endif
  numberIntegers_ = newIntegers;
  integerVariable_ = new int[numberIntegers_];
  OsiObject **temp = new OsiObject *[newNumberObjects];
  // Put integers first
  newIntegers = 0;
  numberIntegers_ = 0;
  for (i = 0; i < numberColumns; i++) {
    int which = mark[i];
    if (which >= 0) {
      if (!isInteger(i)) {
        newIntegers++;
        solver_->setInteger(i);
      }
      if (which < numberColumns) {
        temp[numberIntegers_] = object_[which];
        object_[which] = NULL;
      } else {
        temp[numberIntegers_] = objects[which - numberColumns]->clone();
      }
      integerVariable_[numberIntegers_++] = i;
    }
  }
#if COIN_DEVELOP > 1
  if (newIntegers)
    printf("%d variables were declared integer\n", newIntegers);
#endif
  int n = numberIntegers_;
  // Now rest of old
  for (i = 0; i < numberObjects_; i++) {
    if (object_[i]) {
      CbcSimpleInteger *obj = dynamic_cast< CbcSimpleInteger * >(object_[i]);
      if (obj) {
        delete object_[i];
      } else {
        temp[n++] = object_[i];
      }
    }
  }
  // and rest of new
  for (i = 0; i < numberObjects; i++) {
    CbcSimpleInteger *obj = dynamic_cast< CbcSimpleInteger * >(objects[i]);
    OsiSimpleInteger *obj2 = dynamic_cast< OsiSimpleInteger * >(objects[i]);
    if (!obj && !obj2) {
      temp[n] = objects[i]->clone();
      CbcObject *obj = dynamic_cast< CbcObject * >(temp[n]);
      if (obj)
        obj->setModel(this);
      n++;
    }
  }
  delete[] mark;
  assert(ownObjects_);
  delete[] object_;
  object_ = temp;
  assert(n == newNumberObjects);
  numberObjects_ = newNumberObjects;
}

/**
  This routine sets the objective cutoff value used for fathoming and
  determining monotonic variables.

  If the fathoming discipline is strict, a small tolerance is added to the
  new cutoff. This avoids problems due to roundoff when the target value
  is exact. The common example would be an IP with only integer variables in
  the objective. If the target is set to the exact value z of the optimum,
  it's possible to end up fathoming an ancestor of the solution because the
  solver returns z+epsilon.

  Determining if strict fathoming is needed is best done by analysis.
  In cbc, that's analyseObjective. The default is false.

  In cbc we always minimize so add epsilon
*/

void CbcModel::setCutoff(double value)

{
#ifdef JJF_ZERO
  double tol = 0;
  int fathomStrict = getIntParam(CbcFathomDiscipline);
  if (fathomStrict == 1) {
    solver_->getDblParam(OsiDualTolerance, tol);
    tol = tol * (1 + fabs(value));

    value += tol;
  }
#endif
  dblParam_[CbcCurrentCutoff] = value;
  if (solver_) {
    // Solvers know about direction
    // but Clp tries to be too clever and flips twice!
#ifndef COIN_HAS_CLP
    double direction = solver_->getObjSense();
#else
    double direction = 1.0;
    if (!dynamic_cast< OsiClpSolverInterface * >(solver_))
      direction = solver_->getObjSense();
#endif
    solver_->setDblParam(OsiDualObjectiveLimit, value * direction);
  }
}

/*
  Call this to really test if a valid solution can be feasible. The cutoff is
  passed in as a parameter so that we don't need to worry here after swapping
  solvers.  The solution is assumed to be numberColumns in size.  If
  fixVariables is true then the bounds of the continuous solver are updated.
  The routine returns the objective value determined by reoptimizing from
  scratch. If the solution is rejected, this will be worse than the cutoff.

  TODO: There's an issue with getting the correct cutoff value: We update the
	cutoff in the regular solver, but not in continuousSolver_. But our only
	use for continuousSolver_ is verifying candidate solutions. Would it
	make sense to update the cutoff? Then we wouldn't need to step around
	isDualObjectiveLimitReached().
*/
double
CbcModel::checkSolution(double cutoff, double *solution,
  int fixVariables, double objectiveValue)

{
  int numberContinuousColumns = continuousSolver_->getNumCols();
  if (!solverCharacteristics_->solutionAddsCuts()) {
    // Can trust solution
    int numberColumns = solver_->getNumCols();
#ifdef COIN_HAS_CLP
    OsiClpSolverInterface *clpContinuousSolver
      = dynamic_cast< OsiClpSolverInterface * >(continuousSolver_);
    int modifiedTolerances = 0;
#ifndef CBC_LEAVE_PERTURBATION_ON_CHECK_SOLUTION
    int savePerturbation = -1;
#endif
#ifndef CBC_LEAVE_TOLERANCE_ON_CHECK_SOLUTION
    double savePrimalTolerance = 0.0;
#endif
#ifndef CBC_LEAVE_SCALING_ON_CHECK_SOLUTION
    int saveScaling = -1;
#endif
#define CBC_LEAVE_CRUNCH_ON_CHECK_SOLUTION // for now
#ifndef CBC_LEAVE_CRUNCH_ON_CHECK_SOLUTION
    int saveSpecialOptions=0;
#endif
    if (clpContinuousSolver) {
      // be more accurate if possible
      ClpSimplex *clp = clpContinuousSolver->getModelPtr();
#ifndef CBC_LEAVE_PERTURBATION_ON_CHECK_SOLUTION
      savePerturbation = clp->perturbation();
#endif
#ifndef CBC_LEAVE_TOLERANCE_ON_CHECK_SOLUTION
      savePrimalTolerance = clp->primalTolerance();
#endif
#ifndef CBC_LEAVE_SCALING_ON_CHECK_SOLUTION
      saveScaling = clp->scalingFlag();
#endif
#ifndef CBC_LEAVE_TOLERANCE_ON_CHECK_SOLUTION
      if (savePrimalTolerance > 0.9999999e-7) {
        modifiedTolerances |= 1;
        clp->setPrimalTolerance(1.0e-8);
      }
#endif
#ifndef CBC_LEAVE_PERTURBATION_ON_CHECK_SOLUTION
      if (savePerturbation < 100) {
        modifiedTolerances |= 2;
        clp->setPerturbation(100);
      }
#endif
#ifndef CBC_LEAVE_SCALING_ON_CHECK_SOLUTION
      if (saveScaling) {
        modifiedTolerances |= 4;
        clp->scaling(0);
        clpContinuousSolver->setHintParam(OsiDoScale, false, OsiHintTry);
      }
#endif
#ifndef CBC_LEAVE_CRUNCH_ON_CHECK_SOLUTION
      modifiedTolerances |= 8;
      saveSpecialOptions = clpContinuousSolver->specialOptions();
      clpContinuousSolver->setSpecialOptions(saveSpecialOptions&(~1)); // switch off crunch 
#endif
    }
#endif

    /*
          Grab the continuous solver (the pristine copy of the problem, made before
          starting to work on the root node). Save the bounds on the variables.
          Install the solution passed as a parameter, and copy it to the model's
          currentSolution_.

          TODO: This is a belt-and-suspenders approach. Once the code has settled
          a bit, we can cast a critical eye here.
        */
    OsiSolverInterface *saveSolver = solver_;
    if (continuousSolver_)
      solver_ = continuousSolver_;
    // save basis and solution
    CoinWarmStartBasis *basis = dynamic_cast< CoinWarmStartBasis * >(solver_->getWarmStart());
    assert(basis != NULL);
    double *saveSolution = CoinCopyOfArray(solver_->getColSolution(),
      solver_->getNumCols());
    // move solution to continuous copy
    solver_->setColSolution(solution);
    // Put current solution in safe place
    // Point to current solution
    const double *save = testSolution_;
    // Safe as will be const inside infeasibility()
    testSolution_ = solver_->getColSolution();
    //memcpy(currentSolution_,solver_->getColSolution(),
    // numberColumns*sizeof(double));
    //solver_->messageHandler()->setLogLevel(4);

    // save original bounds
    double *saveUpper = new double[numberColumns];
    double *saveLower = new double[numberColumns];
    memcpy(saveUpper, getColUpper(), numberColumns * sizeof(double));
    memcpy(saveLower, getColLower(), numberColumns * sizeof(double));
    //#define CLP_INVESTIGATE4
#if CBC_USEFUL_PRINTING > 14
    {
      int nBad = checkAssociated(solver_, solver_->getColSolution(), 1);
      if (nBad)
        checkAssociated(solver_, solver_->getColSolution(), 3);
      double largestInfeasibility = 0.0;
      double primalTolerance;
      double offset;
      solver_->getDblParam(OsiObjOffset, offset);
      solver_->getDblParam(OsiPrimalTolerance, primalTolerance);
      const double *objective = getObjCoefficients();
      const double *rowLower = solver_->getRowLower();
      const double *rowUpper = solver_->getRowUpper();
      const double *columnLower = solver_->getColLower();
      const double *columnUpper = solver_->getColUpper();
      int numberRows = solver_->getNumRows();
      double *rowActivity = new double[numberRows];
      memset(rowActivity, 0, numberRows * sizeof(double));
      double *rowSum = new double[numberRows];
      memset(rowSum, 0, numberRows * sizeof(double));
      int *marked = new int[numberColumns];
      for (int i = 0; i < numberColumns; i++)
        marked[i] = -1;
      for (int i = 0; i < numberIntegers_; i++)
        marked[integerVariable_[i]] = -2;
      if ((moreSpecialOptions2_ & 4) != 0) {
        for (int i = 0; i < numberObjects_; i++) {
          CbcSwitchingBinary *object = dynamic_cast< CbcSwitchingBinary * >(object_[i]);
          if (object) {
            int iColumn = object->columnNumber();
            const int *other = object->otherVariable();
            marked[iColumn] = -3 - other[0];
            int n = object->numberOther();
            for (int k = 0; k < n; k++)
              marked[other[k]] = iColumn;
          }
        }
      }
      const double *element = solver_->getMatrixByCol()->getElements();
      const int *row = solver_->getMatrixByCol()->getIndices();
      const CoinBigIndex *columnStart = solver_->getMatrixByCol()->getVectorStarts();
      const int *columnLength = solver_->getMatrixByCol()->getVectorLengths();
      const CoinPackedMatrix *rowCopy = solver_->getMatrixByRow();
      const int *column = rowCopy->getIndices();
      const int *rowLength = rowCopy->getVectorLengths();
      const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
      const double *elementByRow = rowCopy->getElements();
      double objValue = -offset;
      for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
        double value = solution[iColumn];
        objValue += value * objective[iColumn];
        if (value > columnUpper[iColumn]) {
          if (value - columnUpper[iColumn] > 1.0e-8)
            printf("column %d has value %.12g above %.12g\n", iColumn, value, columnUpper[iColumn]);
          value = columnUpper[iColumn];
        } else if (value < columnLower[iColumn]) {
          if (value - columnLower[iColumn] < -1.0e-8)
            printf("column %d has value %.12g below %.12g\n", iColumn, value, columnLower[iColumn]);
          value = columnLower[iColumn];
        }
        if (value) {
          CoinBigIndex start = columnStart[iColumn];
          CoinBigIndex end = start + columnLength[iColumn];
          for (CoinBigIndex j = start; j < end; j++) {
            int iRow = row[j];
            rowActivity[iRow] += value * element[j];
            rowSum[iRow] += fabs(value * element[j]);
          }
        }
      }
      for (int i = 0; i < numberRows; i++) {
#if 0 //def CLP_INVESTIGATE
	    double inf;
	    inf = rowLower[i] - rowActivity[i];
	    if (inf > primalTolerance)
	      printf("Row %d inf %g sum %g %g <= %g <= %g\n",
		     i, inf, rowSum[i], rowLower[i], rowActivity[i], rowUpper[i]);
	    inf = rowActivity[i] - rowUpper[i];
	    if (inf > primalTolerance)
	      printf("Row %d inf %g sum %g %g <= %g <= %g\n",
		     i, inf, rowSum[i], rowLower[i], rowActivity[i], rowUpper[i]);
#endif
        double infeasibility = CoinMax(rowActivity[i] - rowUpper[i],
          rowLower[i] - rowActivity[i]);
        // but allow for errors
        double factor = CoinMax(1.0, rowSum[i] * 1.0e-3);
        if (infeasibility > largestInfeasibility * factor) {
          largestInfeasibility = infeasibility / factor;
          printf("Ainf of %g on row %d sum %g scaled %g\n",
            infeasibility, i, rowSum[i], largestInfeasibility);
          if (infeasibility > 1.0e10) {
            for (CoinBigIndex j = rowStart[i];
                 j < rowStart[i] + rowLength[i]; j++) {
              printf("col %d element %g marked %d\n",
                column[j], elementByRow[j], marked[column[j]]);
            }
          }
        }
      }
      delete[] rowActivity;
      delete[] rowSum;
      delete[] marked;
      if (largestInfeasibility > 10.0 * primalTolerance)
        printf("Alargest infeasibility is %g - obj %g\n", largestInfeasibility, objValue);
      else
        printf("Afeasible (%g) - obj %g\n", largestInfeasibility, objValue);
    }
#endif
    // point to useful information
    OsiBranchingInformation usefulInfo = usefulInformation();

    /*
          Run through the objects and use feasibleRegion() to set variable bounds
          so as to fix the variables specified in the objects at their value in this
          solution. Since the object list contains (at least) one object for every
          integer variable, this has the effect of fixing all integer variables.
        */
    int i;
    for (i = 0; i < numberObjects_; i++)
      object_[i]->feasibleRegion(solver_, &usefulInfo);
    // If SOS then might have been declared infeasible (bad heuristic)
    {
      int numberColumns = solver_->getNumCols();
      const double *columnLower = solver_->getColLower();
      const double *columnUpper = solver_->getColUpper();
      bool looksGood = true;
      for (int i = 0; i < numberColumns; i++) {
        if (columnUpper[i] < columnLower[i])
          looksGood = false;
      }
      if (!looksGood) {
        // not good
        messageHandler()->message(CBC_FPUMP2, messages())
          << "On closer inspection - solution discarded"
          << CoinMessageEol;
        objectiveValue = 1.0e50;
        for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
          solver_->setColLower(iColumn, saveLower[iColumn]);
          solver_->setColUpper(iColumn, saveUpper[iColumn]);
        }
        delete[] saveLower;
        delete[] saveUpper;

        solver_->setColSolution(saveSolution);
        delete[] saveSolution;
        solver_->setWarmStart(basis);
        delete basis;
        /*
	      Restore the usual solver.
	    */
        solver_ = saveSolver;
        testSolution_ = save;
#ifdef COIN_HAS_CLP
        if (modifiedTolerances) {
          // Restore
          ClpSimplex *clp = clpContinuousSolver->getModelPtr();
#ifndef CBC_LEAVE_TOLERANCE_ON_CHECK_SOLUTION
          clp->setPrimalTolerance(savePrimalTolerance);
#endif
#ifndef CBC_LEAVE_PERTURBATION_ON_CHECK_SOLUTION
          clp->setPerturbation(savePerturbation);
#endif
#ifndef CBC_LEAVE_SCALING_ON_CHECK_SOLUTION
          if (saveScaling) {
            clp->scaling(saveScaling);
            clpContinuousSolver->setHintParam(OsiDoScale, true, OsiHintTry);
          }
#endif
#ifndef CBC_LEAVE_CRUNCH_ON_CHECK_SOLUTION
	  // Restore
	  clpContinuousSolver->setSpecialOptions(saveSpecialOptions);
#endif
        }
#endif
        return 1.0e50;
      }
    }
#if CBC_USEFUL_PRINTING > 14
    {
      int nBad = checkAssociated(solver_, solver_->getColSolution(), 1);
      if (nBad)
        checkAssociated(solver_, solver_->getColSolution(), 3);
      double largestInfeasibility = 0.0;
      double primalTolerance;
      double offset;
      solver_->getDblParam(OsiObjOffset, offset);
      solver_->getDblParam(OsiPrimalTolerance, primalTolerance);
      const double *objective = getObjCoefficients();
      const double *rowLower = solver_->getRowLower();
      const double *rowUpper = solver_->getRowUpper();
      const double *columnLower = solver_->getColLower();
      const double *columnUpper = solver_->getColUpper();
      int numberRows = solver_->getNumRows();
      double *rowActivity = new double[numberRows];
      memset(rowActivity, 0, numberRows * sizeof(double));
      double *rowSum = new double[numberRows];
      memset(rowSum, 0, numberRows * sizeof(double));
      int *marked = new int[numberColumns];
      for (int i = 0; i < numberColumns; i++)
        marked[i] = -1;
      for (int i = 0; i < numberIntegers_; i++)
        marked[integerVariable_[i]] = -2;
      if ((moreSpecialOptions2_ & 4) != 0) {
        for (int i = 0; i < numberObjects_; i++) {
          CbcSwitchingBinary *object = dynamic_cast< CbcSwitchingBinary * >(object_[i]);
          if (object) {
            int iColumn = object->columnNumber();
            const int *other = object->otherVariable();
            marked[iColumn] = -3 - other[0];
            int n = object->numberOther();
            for (int k = 0; k < n; k++)
              marked[other[k]] = iColumn;
          }
        }
      }
      const double *element = solver_->getMatrixByCol()->getElements();
      const int *row = solver_->getMatrixByCol()->getIndices();
      const CoinBigIndex *columnStart = solver_->getMatrixByCol()->getVectorStarts();
      const int *columnLength = solver_->getMatrixByCol()->getVectorLengths();
      const CoinPackedMatrix *rowCopy = solver_->getMatrixByRow();
      const int *column = rowCopy->getIndices();
      const int *rowLength = rowCopy->getVectorLengths();
      const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
      const double *elementByRow = rowCopy->getElements();
      double objValue = -offset;
      for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
        double value = solution[iColumn];
        objValue += value * objective[iColumn];
        if (value > columnUpper[iColumn]) {
          if (value - columnUpper[iColumn] > 1.0e-8)
            printf("column %d has value %.12g above %.12g\n", iColumn, value, columnUpper[iColumn]);
          value = columnUpper[iColumn];
        } else if (value < columnLower[iColumn]) {
          if (value - columnLower[iColumn] < -1.0e-8)
            printf("column %d has value %.12g below %.12g\n", iColumn, value, columnLower[iColumn]);
          value = columnLower[iColumn];
        }
        if (value) {
          CoinBigIndex start = columnStart[iColumn];
          CoinBigIndex end = start + columnLength[iColumn];
          for (CoinBigIndex j = start; j < end; j++) {
            int iRow = row[j];
            rowActivity[iRow] += value * element[j];
            rowSum[iRow] += fabs(value * element[j]);
          }
        }
      }
      for (int i = 0; i < numberRows; i++) {
#if 0 //def CLP_INVESTIGATE
	    double inf;
	    inf = rowLower[i] - rowActivity[i];
	    if (inf > primalTolerance)
	      printf("Row %d inf %g sum %g %g <= %g <= %g\n",
		     i, inf, rowSum[i], rowLower[i], rowActivity[i], rowUpper[i]);
	    inf = rowActivity[i] - rowUpper[i];
	    if (inf > primalTolerance)
	      printf("Row %d inf %g sum %g %g <= %g <= %g\n",
		     i, inf, rowSum[i], rowLower[i], rowActivity[i], rowUpper[i]);
#endif
        double infeasibility = CoinMax(rowActivity[i] - rowUpper[i],
          rowLower[i] - rowActivity[i]);
        // but allow for errors
        double factor = CoinMax(1.0, rowSum[i] * 1.0e-3);
        if (infeasibility > largestInfeasibility * factor) {
          largestInfeasibility = infeasibility / factor;
          printf("inf of %g on row %d sum %g scaled %g\n",
            infeasibility, i, rowSum[i], largestInfeasibility);
          if (infeasibility > 1.0e10) {
            for (CoinBigIndex j = rowStart[i];
                 j < rowStart[i] + rowLength[i]; j++) {
              printf("col %d element %g marked %d\n",
                column[j], elementByRow[j], marked[column[j]]);
            }
          }
        }
      }
      delete[] rowActivity;
      delete[] rowSum;
      delete[] marked;
      if (largestInfeasibility > 10.0 * primalTolerance)
        printf("Largest infeasibility is %g - obj %g\n", largestInfeasibility, objValue);
      else
        printf("Feasible (%g) - obj %g\n", largestInfeasibility, objValue);
    }
#endif
    // If relaxed then leave bounds on basic variables
    if (fixVariables == -1 && (specialOptions_ & 16) == 0) {
      CoinWarmStartBasis *basis = dynamic_cast< CoinWarmStartBasis * >(saveSolver->getWarmStart());
      assert(basis != NULL);
#ifdef JJF_ZERO //ndef CBC_OTHER_SOLVER
      for (i = 0; i < numberObjects_; i++) {
        CbcSimpleInteger *obj = dynamic_cast< CbcSimpleInteger * >(object_[i]);
        if (obj) {
          int iColumn = obj->columnNumber();
          if (basis->getStructStatus(iColumn) == CoinWarmStartBasis::basic) {
            solver_->setColLower(iColumn, saveLower[iColumn]);
            solver_->setColUpper(iColumn, saveUpper[iColumn]);
          }
        }
      }
#endif
      delete basis;
    }
    // We can switch off check
    if ((specialOptions_ & 4) == 0 && (moreSpecialOptions2_ & 10) != 8) {
      // Be on safe side - unless very few integers and large
      bool allSlack = (specialOptions_ & 2) == 0 && solverCharacteristics_->warmStart();
      if (numberIntegers_ * 4 > solver_->getNumCols() || solver_->getNumCols() < 10000)
        allSlack = true;
      if (allSlack) {
        /*
                  Remove any existing warm start information to be sure there is no
                  residual influence on initialSolve().
                */
        CoinWarmStartBasis *slack = dynamic_cast< CoinWarmStartBasis * >(solver_->getEmptyWarmStart());
        solver_->setWarmStart(slack);
        delete slack;
      } else {
        if (bestSolutionBasis_.getNumStructural() == solver_->getNumCols() && bestSolutionBasis_.getNumArtificial() == solver_->getNumRows())
          solver_->setWarmStart(&bestSolutionBasis_);
      }
      // Give a hint to do dual
      bool saveTakeHint;
      OsiHintStrength saveStrength;
#ifndef NDEBUG
      bool gotHint = (solver_->getHintParam(OsiDoDualInInitial, saveTakeHint, saveStrength));
      assert(gotHint);
#else
      (solver_->getHintParam(OsiDoDualInInitial, saveTakeHint, saveStrength));
#endif
      solver_->setHintParam(OsiDoDualInInitial, true, OsiHintTry);
      solver_->initialSolve();
#ifdef SWITCH_VARIABLES
      if (solver_->isProvenOptimal()) {
        int nBad = checkAssociated(solver_, solver_->getColSolution(), 1);
        if (nBad)
          checkAssociated(solver_, solver_->getColSolution(), 3);
      }
#endif
#ifdef JJF_ZERO
      if (solver_->isProvenOptimal()) {
        solver_->writeMpsNative("feasible.mps", NULL, NULL, 2);
#ifdef COIN_HAS_CLP
        OsiClpSolverInterface *clpSolver
          = dynamic_cast< OsiClpSolverInterface * >(solver_);
        if (clpSolver) {
          clpSolver->getModelPtr()->writeBasis("feasible.bas", true);
        }
#endif
        printf("XXXXXXXXXXXX - saving feasible\n");
      }
#endif
      if (!solver_->isProvenOptimal()) {
#if CBC_FEASIBILITY_INVESTIGATE
        printf("checkSolution infeas! Retrying with primal.\n");
#endif
        //bool saveTakeHint;
        //OsiHintStrength saveStrength;
        //bool savePrintHint;
        //solver_->writeMpsNative("infeas.mps", NULL, NULL, 2);
        //bool gotHint = (solver_->getHintParam(OsiDoReducePrint,savePrintHint,saveStrength));
        //gotHint = (solver_->getHintParam(OsiDoScale,saveTakeHint,saveStrength));
        //solver_->setHintParam(OsiDoScale,false,OsiHintTry);
        //solver_->setHintParam(OsiDoReducePrint,false,OsiHintTry) ;
        solver_->setHintParam(OsiDoDualInInitial, false, OsiHintTry);
        solver_->initialSolve();
        //solver_->setHintParam(OsiDoScale,saveTakeHint,saveStrength);
        //solver_->setHintParam(OsiDoReducePrint,savePrintHint,OsiHintTry) ;
        // go from all slack now
        specialOptions_ &= ~2;
        if (!solver_->isProvenOptimal()) {
          CoinWarmStartBasis *slack = dynamic_cast< CoinWarmStartBasis * >(solver_->getEmptyWarmStart());
          solver_->setWarmStart(slack);
          delete slack;
#if CBC_FEASIBILITY_INVESTIGATE
          printf("checkSolution infeas! Retrying wihout basis and with primal.\n");
#endif
          solver_->initialSolve();
          //solver_->writeMps("bad");
#ifdef COIN_HAS_CLP
          if (!solver_->isProvenOptimal() && modifiedTolerances) {
            // Restore
            ClpSimplex *clp = clpContinuousSolver->getModelPtr();
#ifndef CBC_LEAVE_TOLERANCE_ON_CHECK_SOLUTION
            clp->setPrimalTolerance(savePrimalTolerance);
#endif
#ifndef CBC_LEAVE_PERTURBATION_ON_CHECK_SOLUTION
            clp->setPerturbation(savePerturbation);
#endif
#ifndef CBC_LEAVE_SCALING_ON_CHECK_SOLUTION
            if (saveScaling) {
              clp->scaling(saveScaling);
              clpContinuousSolver->setHintParam(OsiDoScale, true, OsiHintTry);
            }
#endif
            solver_->resolve();
          }
#endif
#if CBC_FEASIBILITY_INVESTIGATE
          if (!solver_->isProvenOptimal()) {
            printf("checkSolution still infeas!\n");
          }
#endif
        }
      }
      //assert(solver_->isProvenOptimal());
      solver_->setHintParam(OsiDoDualInInitial, saveTakeHint, saveStrength);
      objectiveValue = solver_->isProvenOptimal() ? solver_->getObjValue() * solver_->getObjSense() : 1.0e50;
    }
    bestSolutionBasis_ = CoinWarmStartBasis();

    /*
          Check that the solution still beats the objective cutoff.

          If it passes, make a copy of the primal variable values and do some
          cleanup and checks:
          + Values of all variables are are within original bounds and values of
          all integer variables are within tolerance of integral.
          + There are no constraint violations.
          There really should be no need for the check against original bounds.
          Perhaps an opportunity for a sanity check?
        */
    if (objectiveValue > cutoff && objectiveValue < cutoff + 1.0e-8 + 1.0e-8 * fabs(cutoff))
      cutoff = objectiveValue; // relax
    if ((solver_->isProvenOptimal() || (specialOptions_ & 4) != 0) && objectiveValue <= cutoff) {
      memcpy(solution, solver_->getColSolution(), numberColumns * sizeof(double));
      int iColumn;
#ifndef NDEBUG
      double integerTolerance = getIntegerTolerance();
#endif
#if CBC_FEASIBILITY_INVESTIGATE
      const double *dj = solver_->getReducedCost();
      const double *colLower = saveSolver->getColLower();
      const double *colUpper = saveSolver->getColUpper();
      int nAtLbNatural = 0;
      int nAtUbNatural = 0;
      int nAtLbNaturalZero = 0;
      int nAtUbNaturalZero = 0;
      int nAtLbFixed = 0;
      int nAtUbFixed = 0;
      int nAtOther = 0;
      int nAtOtherNatural = 0;
      int nNotNeeded = 0;
#endif
      for (iColumn = 0; iColumn < numberContinuousColumns; iColumn++) {
        double value = solution[iColumn];
        value = CoinMax(value, saveLower[iColumn]);
        value = CoinMin(value, saveUpper[iColumn]);
        if (solver_->isInteger(iColumn)) {
          assert(fabs(value - solution[iColumn]) <= 100.0 * integerTolerance);
#if CBC_FEASIBILITY_INVESTIGATE
          double value2 = floor(value + 0.5);
          if (dj[iColumn] < -1.0e-6) {
            // negative dj
            //assert (value2==colUpper[iColumn]);
            if (saveUpper[iColumn] == colUpper[iColumn]) {
              nAtUbNatural++;
              if (saveLower[iColumn] != colLower[iColumn])
                nNotNeeded++;
            } else if (saveLower[iColumn] == colUpper[iColumn]) {
              nAtLbFixed++;
            } else {
              nAtOther++;
              if (saveLower[iColumn] != colLower[iColumn] && saveUpper[iColumn] != colUpper[iColumn])
                nNotNeeded++;
            }
          } else if (dj[iColumn] > 1.0e-6) {
            // positive dj
            //assert (value2==colLower[iColumn]);
            if (saveLower[iColumn] == colLower[iColumn]) {
              nAtLbNatural++;
              if (saveUpper[iColumn] != colUpper[iColumn])
                nNotNeeded++;
            } else if (saveUpper[iColumn] == colLower[iColumn]) {
              nAtUbFixed++;
            } else {
              nAtOther++;
              if (saveLower[iColumn] != colLower[iColumn] && saveUpper[iColumn] != colUpper[iColumn])
                nNotNeeded++;
            }
          } else {
            // zero dj
            if (value2 == saveUpper[iColumn]) {
              nAtUbNaturalZero++;
              if (saveLower[iColumn] != colLower[iColumn])
                nNotNeeded++;
            } else if (value2 == saveLower[iColumn]) {
              nAtLbNaturalZero++;
            } else {
              nAtOtherNatural++;
              if (saveLower[iColumn] != colLower[iColumn] && saveUpper[iColumn] != colUpper[iColumn])
                nNotNeeded++;
            }
          }
#endif
        }
        solution[iColumn] = value;
      }
#if CBC_FEASIBILITY_INVESTIGATE
      printf("nAtLbNat %d,nAtUbNat %d,nAtLbNatZero %d,nAtUbNatZero %d,nAtLbFixed %d,nAtUbFixed %d,nAtOther %d,nAtOtherNat %d, useless %d\n",
        nAtLbNatural,
        nAtUbNatural,
        nAtLbNaturalZero,
        nAtUbNaturalZero,
        nAtLbFixed,
        nAtUbFixed,
        nAtOther,
        nAtOtherNatural, nNotNeeded);
      //if (currentNode_)
      //printf(" SOL at depth %d\n",currentNode_->depth());
      //else
      //printf(" SOL at unknown depth\n");
#endif
      if ((specialOptions_ & 16) == 0) {
#ifdef JJF_ZERO
        // check without scaling
        bool saveTakeHint;
        OsiHintStrength saveStrength;
        solver_->getHintParam(OsiDoScale, saveTakeHint, saveStrength);
        solver_->setHintParam(OsiDoScale, false, OsiHintTry);
        solver_->resolve();
        solver_->setHintParam(OsiDoScale, saveTakeHint, saveStrength);
#endif
        double largestInfeasibility = 0.0;
#ifdef COIN_HAS_CLP
        if (clpContinuousSolver) {
          ClpSimplex *clp = clpContinuousSolver->getModelPtr();
          if ((modifiedTolerances & 1) != 0)
            clp->setPrimalTolerance(savePrimalTolerance);
          assert(savePrimalTolerance);
        }
#endif
        double primalTolerance;
        solver_->getDblParam(OsiPrimalTolerance, primalTolerance);
        const double *rowLower = solver_->getRowLower();
        const double *rowUpper = solver_->getRowUpper();
        int numberRows = solver_->getNumRows();
        double *rowActivity = new double[numberRows];
        memset(rowActivity, 0, numberRows * sizeof(double));
        double *rowSum = new double[numberRows];
        memset(rowSum, 0, numberRows * sizeof(double));
        const double *element = solver_->getMatrixByCol()->getElements();
        const int *row = solver_->getMatrixByCol()->getIndices();
        const CoinBigIndex *columnStart = solver_->getMatrixByCol()->getVectorStarts();
        const int *columnLength = solver_->getMatrixByCol()->getVectorLengths();
        double offset;
        solver_->getDblParam(OsiObjOffset, offset);
        double objValue = -offset;
        const double *objective = getObjCoefficients();
        for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
          double value = solution[iColumn];
          objValue += value * objective[iColumn];
          if (value) {
            CoinBigIndex start = columnStart[iColumn];
            CoinBigIndex end = start + columnLength[iColumn];
            for (CoinBigIndex j = start; j < end; j++) {
              int iRow = row[j];
              rowActivity[iRow] += value * element[j];
              rowSum[iRow] += fabs(value * element[j]);
            }
          }
        }
        for (i = 0; i < numberRows; i++) {
#if CBC_FEASIBILITY_INVESTIGATE > 1
          double inf;
          inf = rowLower[i] - rowActivity[i];
          if (inf > primalTolerance)
            printf("Row %d inf %g sum %g %g <= %g <= %g\n",
              i, inf, rowSum[i], rowLower[i], rowActivity[i], rowUpper[i]);
          inf = rowActivity[i] - rowUpper[i];
          if (inf > primalTolerance)
            printf("Row %d inf %g sum %g %g <= %g <= %g\n",
              i, inf, rowSum[i], rowLower[i], rowActivity[i], rowUpper[i]);
#endif
          double infeasibility = CoinMax(rowActivity[i] - rowUpper[i],
            rowLower[i] - rowActivity[i]);
          // but allow for errors
          double factor = CoinMax(1.0, rowSum[i] * 1.0e-3);
          if (infeasibility > largestInfeasibility * factor) {
            largestInfeasibility = infeasibility / factor;
            //printf("inf of %g on row %d sum %g scaled %g\n",
            //     infeasibility,i,rowSum[i],largestInfeasibility);
          }
        }
        delete[] rowActivity;
        delete[] rowSum;
#if CBC_FEASIBILITY_INVESTIGATE == 0
        if (handler_->logLevel() > 2) {
#endif
          if (largestInfeasibility > 10.0 * primalTolerance)
            printf("BLargest infeasibility is %g - obj %g (%g)\n", largestInfeasibility, objValue, objectiveValue);
          else
            printf("BFeasible (%g) - obj %g %g\n", largestInfeasibility, objValue, objectiveValue);
#if CBC_FEASIBILITY_INVESTIGATE == 0
        }
#else
        solver_->writeMpsNative("BFeasible.mps", NULL, NULL, 2);
#endif
        //if (fabs(objValue-objectiveValue)>1.0e-7*fabs(objectiveValue)) {
        //printf("Bad obj values\n");
        objectiveValue = objValue;
        //}
#if CBC_FEASIBILITY_INVESTIGATE
        if (largestInfeasibility > 10.0 * primalTolerance)
          printf("XX largest infeasibility is %g\n", largestInfeasibility);
#endif
        if (largestInfeasibility > 200.0 * primalTolerance) {
          handler_->message(CBC_NOTFEAS3, messages_)
            << largestInfeasibility << CoinMessageEol;
          objectiveValue = 1.0e50;
        }
      }
    } else {
      objectiveValue = 1.0e50;
    }
    /*
          Regardless of what we think of the solution, we may need to restore the
          original bounds of the continuous solver. Unfortunately, const'ness
          prevents us from simply reversing the memcpy used to make these snapshots.
        */
    if (fixVariables <= 0) {
      for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
        solver_->setColLower(iColumn, saveLower[iColumn]);
        solver_->setColUpper(iColumn, saveUpper[iColumn]);
      }
    }
    delete[] saveLower;
    delete[] saveUpper;

    solver_->setColSolution(saveSolution);
    delete[] saveSolution;
    solver_->setWarmStart(basis);
    delete basis;
    /*
          Restore the usual solver.
        */
    solver_ = saveSolver;
    testSolution_ = save;
#ifdef COIN_HAS_CLP
    if (modifiedTolerances) {
      // Restore
      ClpSimplex *clp = clpContinuousSolver->getModelPtr();
#ifndef CBC_LEAVE_TOLERANCE_ON_CHECK_SOLUTION
      clp->setPrimalTolerance(savePrimalTolerance);
#endif
#ifndef CBC_LEAVE_PERTURBATION_ON_CHECK_SOLUTION
      clp->setPerturbation(savePerturbation);
#endif
#ifndef CBC_LEAVE_SCALING_ON_CHECK_SOLUTION
      if (saveScaling) {
        clp->scaling(saveScaling);
        clpContinuousSolver->setHintParam(OsiDoScale, true, OsiHintTry);
      }
#endif
#ifndef CBC_LEAVE_CRUNCH_ON_CHECK_SOLUTION
      clpContinuousSolver->setSpecialOptions(saveSpecialOptions);
#endif
    }
#endif
    return objectiveValue;
  } else {
    // Outer approximation or similar
    //If this is true then the solution comes from the nlp we don't need to resolve the same nlp with ipopt
    //solverCharacteristics_->setSolver(solver_);
    bool solutionComesFromNlp = solverCharacteristics_->bestObjectiveValue() < cutoff;
    double objectiveValue;
    int numberColumns = solver_->getNumCols();
    double *saveLower = NULL;
    double *saveUpper = NULL;

    if (!solutionComesFromNlp) { //Otherwise solution already comes from ipopt and cuts are known
      if (fixVariables > 0) { //Will temporarily fix all integer valued var
        // save original bounds
        saveUpper = new double[numberColumns];
        saveLower = new double[numberColumns];
        memcpy(saveUpper, solver_->getColUpper(), numberColumns * sizeof(double));
        memcpy(saveLower, solver_->getColLower(), numberColumns * sizeof(double));
        //in any case solution should be already loaded into solver_
        /*
                  Run through the objects and use feasibleRegion() to set variable bounds
                  so as to fix the variables specified in the objects at their value in this
                  solution. Since the object list contains (at least) one object for every
                  integer variable, this has the effect of fixing all integer variables.
                */
        const double *save = testSolution_;
        testSolution_ = solution;
        // point to useful information
        OsiBranchingInformation usefulInfo = usefulInformation();
        for (int i = 0; i < numberObjects_; i++)
          object_[i]->feasibleRegion(solver_, &usefulInfo);
        testSolution_ = save;
        resolve(solver_);
      }

      /*
              Now step through the cut generators and see if any of them are flagged to
              run when a new solution is discovered. Only global cuts are useful.
              (The solution being evaluated may not correspond to the current location in the
              search tree --- discovered by heuristic, for example.)
            */
      OsiCuts theseCuts;
      int i;
      int lastNumberCuts = 0;
      // reset probing info
      //if (probingInfo_)
      //probingInfo_->initializeFixing();
      for (i = 0; i < numberCutGenerators_; i++) {
        if (generator_[i]->atSolution()) {
          generator_[i]->generateCuts(theseCuts, 1, solver_, NULL);
          int numberCuts = theseCuts.sizeRowCuts();
          for (int j = lastNumberCuts; j < numberCuts; j++) {
            const OsiRowCut *thisCut = theseCuts.rowCutPtr(j);
            if (thisCut->globallyValid()) {
              //           if ((specialOptions_&1)!=0)
              //           {
              //             /* As these are global cuts -
              //             a) Always get debugger object
              // b) Not fatal error to cutoff optimal (if we have just got optimal)
              // */
              //             const OsiRowCutDebugger *debugger = solver_->getRowCutDebuggerAlways() ;
              //             if (debugger)
              //             {
              //               if(debugger->invalidCut(*thisCut))
              //                 printf("ZZZZ Global cut - cuts off optimal solution!\n");
              //             }
              //           }
              // add to global list
              OsiRowCut newCut(*thisCut);
              newCut.setGloballyValid(true);
              newCut.mutableRow().setTestForDuplicateIndex(false);
              globalCuts_.addCutIfNotDuplicate(newCut);
            } else {
              // obviously wrong
              if (handler_->logLevel() > 1)
                printf("Cut generator %s set to run on new solution but NOT globally valid!!\n",
                  generator_[i]->cutGeneratorName());
            }
          }
        }
      }
      //   int numberCuts = theseCuts.sizeColCuts();
      //   for (i=0;i<numberCuts;i++) {
      //     const OsiColCut * thisCut = theseCuts.colCutPtr(i);
      //     if (thisCut->globallyValid()) {
      //       // add to global list
      //       globalCuts_.insert(*thisCut);
      //     }
      //   }
      //have to retrieve the solution and its value from the nlp
    }
    double newObjectiveValue = cutoff;
    if (solverCharacteristics_->solution(newObjectiveValue,
          const_cast< double * >(solution),
          numberColumns)) {
      objectiveValue = newObjectiveValue;
    } else {
      objectiveValue = 2e50;
    }
    if (!solutionComesFromNlp && fixVariables > 0) {
      for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
        solver_->setColLower(iColumn, saveLower[iColumn]);
        solver_->setColUpper(iColumn, saveUpper[iColumn]);
      }
      delete[] saveLower;
      delete[] saveUpper;
      solver_->resolve();
    }
    //If the variables were fixed the cutting plane procedure may have believed that the node could be fathomed
    //re-establish truth.- should do no harm for non nlp
    if (!solutionComesFromNlp && fixVariables > 0)
      solverCharacteristics_->setMipBound(-COIN_DBL_MAX);
    return objectiveValue;
  }
}

/*
  Call this routine from anywhere when a solution is found. The solution
  vector is assumed to contain one value for each structural variable.

  The first action is to run checkSolution() to confirm the objective and
  feasibility. If this check causes the solution to be rejected, we're done.
  If fixVariables = true, the variable bounds held by the continuous solver
  will be left fixed to the values in the solution; otherwise they are
  restored to the original values.

  If the solution is accepted, install it as the best solution.

  The routine also contains a hook to run any cut generators that are flagged
  to run when a new solution is discovered. There's a potential hazard because
  the cut generators see the continuous solver >after< possible restoration of
  original bounds (which may well invalidate the solution).
*/

void CbcModel::setBestSolution(CBC_Message how,
  double &objectiveValue, const double *solutionIn,
  int fixVariables)

{

  double *solution = CoinCopyOfArray(solutionIn, solver_->getNumCols());
#ifdef JJF_ZERO
  {
    double saveOffset;
    solver_->getDblParam(OsiObjOffset, saveOffset);
    const double *obj = solver_->getObjCoefficients();
    double newTrueSolutionValue = -saveOffset;
    double newSumInfeas = 0.0;
    int numberColumns = solver_->getNumCols();
    for (int i = 0; i < numberColumns; i++) {
      if (solver_->isInteger(i)) {
        double value = solution[i];
        double nearest = floor(value + 0.5);
        newSumInfeas += fabs(value - nearest);
      }
      if (solution[i])
        printf("%d obj %g val %g - total %g true\n", i, obj[i], solution[i],
          newTrueSolutionValue);
      newTrueSolutionValue += obj[i] * solution[i];
    }
    printf("obj %g\n", newTrueSolutionValue);
  }
#endif
  if (!solverCharacteristics_->solutionAddsCuts()) {
    // Can trust solution
    double cutoff = getCutoff();
    if (cutoff < 1.0e30)
      cutoff = CoinMin(cutoff, bestObjective_);

    /*
          Double check the solution to catch pretenders.
        */
    double saveObjectiveValue = objectiveValue;
    // save basis
    CoinWarmStartBasis *basis = dynamic_cast< CoinWarmStartBasis * >(solver_->getWarmStart());
    assert(basis != NULL);
    objectiveValue = checkSolution(cutoff, solution, fixVariables, objectiveValue);
    if (cutoff > 1.0e40 && objectiveValue < 1.0e10)
      saveObjectiveValue = objectiveValue; // take anyway
    if (saveObjectiveValue + 1.0e-3 + 1.0e-7 * fabs(saveObjectiveValue)
      < objectiveValue) {
#if CBC_FEASIBILITY_INVESTIGATE
      printf("First try at solution had objective %.16g, rechecked as %.16g\n",
        saveObjectiveValue, objectiveValue);
#endif
      // try again with basic variables with original bounds
      // save basis
      CoinWarmStartBasis *basis2 = dynamic_cast< CoinWarmStartBasis * >(solver_->getWarmStart());
      assert(basis2 != NULL);
      solver_->setWarmStart(basis);
      int numberColumns = solver_->getNumCols();
      double *solution2 = CoinCopyOfArray(solutionIn, numberColumns);
      double objectiveValue2 = saveObjectiveValue;
      objectiveValue2 = checkSolution(cutoff, solution2, -1, objectiveValue2);
#if CBC_FEASIBILITY_INVESTIGATE
      printf("Relaxed second try had objective of %.16g\n",
        objectiveValue2);
#endif
      if (objectiveValue2 + 1.0e-7 < objectiveValue) {
        // Now check tolerances
        double integerTolerance = dblParam_[CbcIntegerTolerance];
        double tolerance;
        solver_->getDblParam(OsiPrimalTolerance, tolerance);
        double largestAway = 0.0;
        int iAway = -1;
        double largestInfeasibility = tolerance;
#if CBC_FEASIBILITY_INVESTIGATE
        int iInfeas = -1;
#endif
        const double *columnLower = continuousSolver_->getColLower();
        const double *columnUpper = continuousSolver_->getColUpper();
        int i;
        for (i = 0; i < numberColumns; i++) {
          double value = solution2[i];
          if (value > columnUpper[i] + largestInfeasibility) {
#if CBC_FEASIBILITY_INVESTIGATE
            iInfeas = i;
#endif
            largestInfeasibility = value - columnUpper[i];
          } else if (value < columnLower[i] - largestInfeasibility) {
#if CBC_FEASIBILITY_INVESTIGATE
            iInfeas = i;
#endif
            largestInfeasibility = columnLower[i] - value;
          }
        }
        for (i = 0; i < numberObjects_; i++) {
          CbcSimpleInteger *obj = dynamic_cast< CbcSimpleInteger * >(object_[i]);
          if (obj) {
            int iColumn = obj->columnNumber();
            double value = solution2[iColumn];
            value = fabs(floor(value + 0.5) - value);
            if (value > largestAway) {
              iAway = iColumn;
              largestAway = value;
            }
          }
        }
#if CBC_FEASIBILITY_INVESTIGATE
        if (iInfeas >= 0)
          printf("Largest infeasibility of %g on column %d - tolerance %g\n",
            largestInfeasibility, iInfeas, tolerance);
#endif
        if (largestAway > integerTolerance) {
          handler_->message(CBC_RELAXED1, messages_)
            << objectiveValue2
            << iAway
            << largestAway
            << integerTolerance
            << CoinMessageEol;
        } else {
          handler_->message(CBC_RELAXED2, messages_)
            << objectiveValue2
            << integerTolerance
            << CoinMessageEol;
          // take
          CoinCopyN(solution2, numberColumns, solution);
          objectiveValue = objectiveValue2;
        }
      } else if (!parentModel_) {
        // not good
        messageHandler()->message(CBC_FPUMP2, messages())
          << "On closer inspection - solution discarded"
          << CoinMessageEol;
      }
      delete[] solution2;
      solver_->setWarmStart(basis2);
      delete basis2;
    }
    delete basis;
    if (objectiveValue > cutoff && objectiveValue < cutoff + 1.0e-8 + 1.0e-8 * fabs(cutoff))
      cutoff = objectiveValue; // relax
    CbcEventHandler::CbcAction action = dealWithEventHandler(CbcEventHandler::beforeSolution2,
      objectiveValue, solution);
    if (action == CbcEventHandler::killSolution) {
      // Pretend solution never happened
      objectiveValue = cutoff + 1.0e30;
    }
    if (objectiveValue > cutoff || objectiveValue > 1.0e30) {
      if (objectiveValue > 1.0e30)
        handler_->message(CBC_NOTFEAS1, messages_) << CoinMessageEol;
      else
        handler_->message(CBC_NOTFEAS2, messages_)
          << objectiveValue << cutoff << CoinMessageEol;
    } else if (objectiveValue < bestObjective_) {
      /*
              We have a winner. Install it as the new incumbent.
              Bump the objective cutoff value and solution counts. Give the user the
              good news.
            */
      specialOptions_ |= 256; // mark as full cut scan should be done
      saveBestSolution(solution, objectiveValue);
      //bestObjective_ = objectiveValue;
      //int numberColumns = solver_->getNumCols();
      //if (!bestSolution_)
      //bestSolution_ = new double[numberColumns];
      //CoinCopyN(solution,numberColumns,bestSolution_);

      cutoff = bestObjective_ - dblParam_[CbcCutoffIncrement];
      // But allow for rounding errors
      if (dblParam_[CbcCutoffIncrement] == 1e-5) {
#if CBC_FEASIBILITY_INVESTIGATE
        if (saveObjectiveValue + 1.0e-7 < bestObjective_)
          printf("First try at solution had objective %.16g, rechecked as %.16g\n",
            saveObjectiveValue, bestObjective_);
#endif
        saveObjectiveValue = CoinMax(saveObjectiveValue, bestObjective_ - 0.0000001 * fabs(bestObjective_));
        cutoff = CoinMin(bestObjective_, saveObjectiveValue) - 1.0e-5;
        if (fabs(cutoff + 1.0e-5 - floor(cutoff + 0.5)) < 1.0e-8)
          cutoff -= 2.0e-5;
      }
      if (!parentModel_ && (moreSpecialOptions2_ & 2) != 0) {
        // put back objective
        solver_->setObjective(continuousSolver_->getObjCoefficients());
        double offset;
        continuousSolver_->getDblParam(OsiObjOffset, offset);
        solver_->setDblParam(OsiObjOffset, offset);
        moreSpecialOptions2_ &= ~2;
      }
      // This is not correct - that way cutoff can go up if maximization
      //double direction = solver_->getObjSense();
      //setCutoff(cutoff*direction);
      setCutoff(cutoff);
      // change cutoff as constraint if wanted
      if (cutoffRowNumber_ >= 0) {
        if (solver_->getNumRows() > cutoffRowNumber_) {
          double offset;
          solver_->getDblParam(OsiObjOffset, offset);
          solver_->setRowUpper(cutoffRowNumber_, cutoff + offset);
        }
      }

      if (how == CBC_ROUNDING)
        numberHeuristicSolutions_++;
      numberSolutions_++;

      if (how != CBC_ROUNDING) {
        handler_->message(how, messages_)
          << bestObjective_ << numberIterations_
          << numberNodes_ << getCurrentSeconds()
          << CoinMessageEol;
        dealWithEventHandler(CbcEventHandler::solution,
          objectiveValue, solution);
      } else {
        const char *name;
        if (lastHeuristic_)
          name = lastHeuristic_->heuristicName();
        else
          name = "Reduced search";
        handler_->message(CBC_ROUNDING, messages_)
          << bestObjective_
          << name
          << numberIterations_
          << numberNodes_ << getCurrentSeconds()
          << CoinMessageEol;
        dealWithEventHandler(CbcEventHandler::heuristicSolution,
          objectiveValue, solution);
      }
      /*
              Now step through the cut generators and see if any of them are flagged to
              run when a new solution is discovered. Only global cuts are useful. (The
              solution being evaluated may not correspond to the current location in the
              search tree --- discovered by heuristic, for example.)
            */
      OsiCuts theseCuts;
      int i;
      int lastNumberCuts = 0;
      // reset probing info
      //if (probingInfo_)
      //probingInfo_->initializeFixing();
      for (i = 0; i < numberCutGenerators_; i++) {
        bool generate = generator_[i]->atSolution();
        // skip if not optimal and should be (maybe a cut generator has fixed variables)
        if (generator_[i]->needsOptimalBasis() && !solver_->basisIsAvailable())
          generate = false;
        if (generate) {
          generator_[i]->generateCuts(theseCuts, 1, solver_, NULL);
          int numberCuts = theseCuts.sizeRowCuts();
          for (int j = lastNumberCuts; j < numberCuts; j++) {
            const OsiRowCut *thisCut = theseCuts.rowCutPtr(j);
            if (thisCut->globallyValid()) {
              if ((specialOptions_ & 1) != 0) {
                /* As these are global cuts -
                                   a) Always get debugger object
                                   b) Not fatal error to cutoff optimal (if we have just got optimal)
                                */
                const OsiRowCutDebugger *debugger = solver_->getRowCutDebuggerAlways();
                if (debugger) {
                  if (debugger->invalidCut(*thisCut))
                    printf("ZZZZ Global cut - cuts off optimal solution!\n");
                }
              }
              // add to global list
              OsiRowCut newCut(*thisCut);
              newCut.setGloballyValid(true);
              newCut.mutableRow().setTestForDuplicateIndex(false);
              globalCuts_.addCutIfNotDuplicate(newCut);
              generator_[i]->incrementNumberCutsInTotal();
            }
          }
        }
      }
      int numberCuts = theseCuts.sizeColCuts();
      for (i = 0; i < numberCuts; i++) {
        const OsiColCut *thisCut = theseCuts.colCutPtr(i);
        if (thisCut->globallyValid()) {
          // fix
          makeGlobalCut(thisCut);
        }
      }
    }
  } else {
    // Outer approximation or similar
    double cutoff = getCutoff();

    /*
          Double check the solution to catch pretenders.
        */

    int numberRowBefore = solver_->getNumRows();
    int numberColBefore = solver_->getNumCols();
    double *saveColSol = NULL;

    CoinWarmStart *saveWs = NULL;
    // if(how!=CBC_SOLUTION) return;
    if (how == CBC_ROUNDING) //We don't want to make any change to solver_
    //take a snapshot of current state
    {
      //save solution
      saveColSol = new double[numberColBefore];
      CoinCopyN(solver_->getColSolution(), numberColBefore, saveColSol);
      //save warm start
      saveWs = solver_->getWarmStart();
    }

    //run check solution this will eventually generate cuts
    //if in strongBranching or heuristic will do only one cut generation iteration
    // by fixing variables.
    if (!fixVariables && ((how == CBC_ROUNDING) || (how == CBC_STRONGSOL)))
      fixVariables = 1;
    double *candidate = new double[numberColBefore];
    CoinCopyN(solution, numberColBefore, candidate);
    objectiveValue = checkSolution(cutoff, candidate, fixVariables, objectiveValue);

    //If it was an heuristic solution we have to clean up the solver
    if (how == CBC_ROUNDING) {
      //delete the cuts
      int currentNumberRowCuts = solver_->getNumRows() - numberRowBefore;
      int currentNumberColCuts = solver_->getNumCols() - numberColBefore;
      if (CoinMax(currentNumberColCuts, currentNumberRowCuts) > 0) {
        int *which = new int[CoinMax(currentNumberColCuts, currentNumberRowCuts)];
        if (currentNumberRowCuts) {
          for (int i = 0; i < currentNumberRowCuts; i++)
            which[i] = i + numberRowBefore;

          solver_->deleteRows(currentNumberRowCuts, which);
        }
        if (currentNumberColCuts) {
          for (int i = 0; i < currentNumberColCuts; i++)
            which[i] = i + numberColBefore;
          solver_->deleteCols(currentNumberColCuts, which);
        }
        delete[] which;
      }
      // Reset solution and warm start info
      solver_->setColSolution(saveColSol);
      solver_->setWarmStart(saveWs);
      delete[] saveColSol;
      delete saveWs;
    }

    if (objectiveValue > cutoff) {
      // message only for solution
      if (how == CBC_SOLUTION) {
        if (!solverCharacteristics_->solutionAddsCuts()) {
          if (objectiveValue > 1.0e30)
            handler_->message(CBC_NOTFEAS1, messages_) << CoinMessageEol;
          else
            handler_->message(CBC_NOTFEAS2, messages_)
              << objectiveValue << cutoff << CoinMessageEol;
        }
      }
    } else {
      /*
              We have a winner. Install it as the new incumbent.
              Bump the objective cutoff value and solution counts. Give the user the
              good news.
              NB - Not all of this if from solve with cuts
            */
      saveBestSolution(candidate, objectiveValue);
      //bestObjective_ = objectiveValue;
      //int numberColumns = solver_->getNumCols();
      //if (!bestSolution_)
      //bestSolution_ = new double[numberColumns];
      //CoinCopyN(candidate,numberColumns,bestSolution_);

      // don't update if from solveWithCuts
      if (how != CBC_SOLUTION2) {
        if (how == CBC_ROUNDING)
          numberHeuristicSolutions_++;
        cutoff = bestObjective_ - dblParam_[CbcCutoffIncrement];
        // This is not correct - that way cutoff can go up if maximization
        //double direction = solver_->getObjSense();
        //setCutoff(cutoff*direction);
        setCutoff(cutoff);
        // change cutoff as constraint if wanted
        if (cutoffRowNumber_ >= 0) {
          if (solver_->getNumRows() > cutoffRowNumber_) {
            double offset;
            solver_->getDblParam(OsiObjOffset, offset);
            solver_->setRowUpper(cutoffRowNumber_, cutoff + offset);
          }
        }

        numberSolutions_++;

        if (how != CBC_ROUNDING) {
          handler_->message(how, messages_)
            << bestObjective_ << numberIterations_
            << numberNodes_ << getCurrentSeconds()
            << CoinMessageEol;
        } else {
          assert(lastHeuristic_);
          const char *name = lastHeuristic_->heuristicName();
          handler_->message(CBC_ROUNDING, messages_)
            << bestObjective_
            << name
            << numberIterations_
            << numberNodes_ << getCurrentSeconds()
            << CoinMessageEol;
        }
      }
    }
    delete[] candidate;
  }
  delete[] solution;
  return;
}
// Deals with event handler and solution
CbcEventHandler::CbcAction
CbcModel::dealWithEventHandler(CbcEventHandler::CbcEvent event,
  double objValue,
  const double *solution)
{
  CbcEventHandler *eventHandler = getEventHandler();
  if (eventHandler) {
    // Temporarily put in best
    double saveObj = bestObjective_;
    int numberColumns = solver_->getNumCols();
    double *saveSol = CoinCopyOfArray(bestSolution_, numberColumns);
    if (!saveSol)
      bestSolution_ = new double[numberColumns];
    bestObjective_ = objValue;
    memcpy(bestSolution_, solution, numberColumns * sizeof(double));
    CbcEventHandler::CbcAction action = eventHandler->event(event);
    bestObjective_ = saveObj;
    if (saveSol) {
      memcpy(bestSolution_, saveSol, numberColumns * sizeof(double));
      delete[] saveSol;
    } else {
      delete[] bestSolution_;
      bestSolution_ = NULL;
    }
    return action;
  } else {
    return CbcEventHandler::noAction;
  }
}

/* Test the current solution for feasibility.

   Calculate the number of standard integer infeasibilities, then scan the
   remaining objects to see if any of them report infeasibilities.

   Currently (2003.08) the only object besides SimpleInteger is Clique, hence
   the comments about `odd ones' infeasibilities.
*/
bool CbcModel::feasibleSolution(int &numberIntegerInfeasibilities,
  int &numberObjectInfeasibilities) const
{
  int numberUnsatisfied = 0;
  //double sumUnsatisfied=0.0;
  int j;
  // Point to current solution
  const double *save = testSolution_;
  // Safe as will be const inside infeasibility()
  testSolution_ = solver_->getColSolution();
  // Put current solution in safe place
  //memcpy(currentSolution_,solver_->getColSolution(),
  // solver_->getNumCols()*sizeof(double));
  // point to useful information
  OsiBranchingInformation usefulInfo = usefulInformation();
#define SIMPLE_INTEGER
#ifdef SIMPLE_INTEGER
  const double *solution = usefulInfo.solution_;
  const double *lower = usefulInfo.lower_;
  const double *upper = usefulInfo.upper_;
  double tolerance = usefulInfo.integerTolerance_;
#endif
  for (j = 0; j < numberIntegers_; j++) {
#ifndef SIMPLE_INTEGER
    const OsiObject *object = object_[j];
    double infeasibility = object->checkInfeasibility(&usefulInfo);
    if (infeasibility) {
      assert(infeasibility > 0);
      numberUnsatisfied++;
      //sumUnsatisfied += infeasibility;
    }
#else
    int iColumn = integerVariable_[j];
    double value = solution[iColumn];
    value = CoinMax(value, lower[iColumn]);
    value = CoinMin(value, upper[iColumn]);
    double nearest = floor(value + 0.5);
    if (fabs(value - nearest) > tolerance) {
      numberUnsatisfied++;
    }
#endif
  }
  numberIntegerInfeasibilities = numberUnsatisfied;
  for (; j < numberObjects_; j++) {
    const OsiObject *object = object_[j];
    double infeasibility = object->checkInfeasibility(&usefulInfo);
    if (infeasibility) {
      assert(infeasibility > 0);
      numberUnsatisfied++;
      //sumUnsatisfied += infeasibility;
    }
  }
  // and restore
  testSolution_ = save;
  numberObjectInfeasibilities = numberUnsatisfied - numberIntegerInfeasibilities;
  return (!numberUnsatisfied);
}

/* For all vubs see if we can tighten bounds by solving Lp's
   type - 0 just vubs
   1 all (could be very slow)
   -1 just vubs where variable away from bound
   Returns false if not feasible
*/
bool CbcModel::tightenVubs(int type, bool allowMultipleBinary, double useCutoff)
{

  CoinPackedMatrix matrixByRow(*solver_->getMatrixByRow());
  int numberRows = solver_->getNumRows();
  int numberColumns = solver_->getNumCols();

  int iRow, iColumn;

  // Row copy
  //const double * elementByRow = matrixByRow.getElements();
  const int *column = matrixByRow.getIndices();
  const CoinBigIndex *rowStart = matrixByRow.getVectorStarts();
  const int *rowLength = matrixByRow.getVectorLengths();

  const double *colUpper = solver_->getColUpper();
  const double *colLower = solver_->getColLower();
  //const double * rowUpper = solver_->getRowUpper();
  //const double * rowLower = solver_->getRowLower();

  const double *objective = solver_->getObjCoefficients();
  //double direction = solver_->getObjSense();
  const double *colsol = solver_->getColSolution();

  int numberVub = 0;
  int *continuous = new int[numberColumns];
  if (type >= 0) {
    double *sort = new double[numberColumns];
    for (iRow = 0; iRow < numberRows; iRow++) {
      CoinBigIndex j;
      int numberBinary = 0;
      int numberUnsatisfiedBinary = 0;
      int numberContinuous = 0;
      int iCont = -1;
      double weight = 1.0e30;
      for (j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
        int iColumn = column[j];
        if (colUpper[iColumn] - colLower[iColumn] > 1.0e-8) {
          if (solver_->isFreeBinary(iColumn)) {
            numberBinary++;
            /* For sort I make naive assumption:
                           x - a * delta <=0 or
                           -x + a * delta >= 0
                        */
            if (colsol[iColumn] > colLower[iColumn] + 1.0e-6 && colsol[iColumn] < colUpper[iColumn] - 1.0e-6) {
              numberUnsatisfiedBinary++;
              weight = CoinMin(weight, fabs(objective[iColumn]));
            }
          } else {
            numberContinuous++;
            iCont = iColumn;
          }
        }
      }
      if (numberContinuous == 1 && numberBinary) {
        if (numberBinary == 1 || allowMultipleBinary) {
          // treat as vub
          if (!numberUnsatisfiedBinary)
            weight = -1.0; // at end
          sort[numberVub] = -weight;
          continuous[numberVub++] = iCont;
        }
      }
    }
    if (type > 0) {
      // take so many
      CoinSort_2(sort, sort + numberVub, continuous);
      numberVub = CoinMin(numberVub, type);
    }
    delete[] sort;
  } else {
    for (iColumn = 0; iColumn < numberColumns; iColumn++)
      continuous[iColumn] = iColumn;
    numberVub = numberColumns;
  }
  bool feasible = tightenVubs(numberVub, continuous, useCutoff);
  delete[] continuous;

  return feasible;
}
// This version is just handed a list of variables
bool CbcModel::tightenVubs(int numberSolves, const int *which,
  double useCutoff)
{

  int numberColumns = solver_->getNumCols();

  int iColumn;

  OsiSolverInterface *solver = solver_;
  double saveCutoff = getCutoff();

  double *objective = new double[numberColumns];
  memcpy(objective, solver_->getObjCoefficients(), numberColumns * sizeof(double));
  double direction = solver_->getObjSense();

  // add in objective if there is a cutoff
  if (useCutoff < 1.0e30) {
    // get new version of model
    solver = solver_->clone();
    CoinPackedVector newRow;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      solver->setObjCoeff(iColumn, 0.0); // zero out in new model
      if (objective[iColumn])
        newRow.insert(iColumn, direction * objective[iColumn]);
    }
    solver->addRow(newRow, -COIN_DBL_MAX, useCutoff);
    // signal no objective
    delete[] objective;
    objective = NULL;
  }
  setCutoff(COIN_DBL_MAX);

  bool *vub = new bool[numberColumns];
  int iVub;

  // mark vub columns
  for (iColumn = 0; iColumn < numberColumns; iColumn++)
    vub[iColumn] = false;
  for (iVub = 0; iVub < numberSolves; iVub++)
    vub[which[iVub]] = true;
  OsiCuts cuts;
  // First tighten bounds anyway if CglProbing there
  CglProbing *generator = NULL;
  int iGen;
  // reset probing info
  //if (probingInfo_)
  //probingInfo_->initializeFixing();
  for (iGen = 0; iGen < numberCutGenerators_; iGen++) {
    generator = dynamic_cast< CglProbing * >(generator_[iGen]->generator());
    if (generator)
      break;
  }
  int numberFixed = 0;
  int numberTightened = 0;
  int numberFixedByProbing = 0;
  int numberTightenedByProbing = 0;
  int printFrequency = (numberSolves + 19) / 20; // up to 20 messages
  int save[4] = { 0, 0, 0, 0 };
  if (generator) {
    // set to cheaper and then restore at end
    save[0] = generator->getMaxPass();
    save[1] = generator->getMaxProbe();
    save[2] = generator->getMaxLook();
    save[3] = generator->rowCuts();
    generator->setMaxPass(1);
    generator->setMaxProbe(10);
    generator->setMaxLook(50);
    generator->setRowCuts(0);

    // Probing - return tight column bounds
    CglTreeInfo info;
    generator->generateCutsAndModify(*solver, cuts, &info);
    const double *tightLower = generator->tightLower();
    const double *lower = solver->getColLower();
    const double *tightUpper = generator->tightUpper();
    const double *upper = solver->getColUpper();
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      double newUpper = tightUpper[iColumn];
      double newLower = tightLower[iColumn];
      if (newUpper < upper[iColumn] - 1.0e-8 * (fabs(upper[iColumn]) + 1) || newLower > lower[iColumn] + 1.0e-8 * (fabs(lower[iColumn]) + 1)) {
        if (newUpper < newLower) {
          fprintf(stderr, "Problem is infeasible\n");
          return false;
        }
        if (newUpper == newLower) {
          numberFixed++;
          numberFixedByProbing++;
          solver->setColLower(iColumn, newLower);
          solver->setColUpper(iColumn, newUpper);
          COIN_DETAIL_PRINT(printf("Column %d, new bounds %g %g\n", iColumn,
            newLower, newUpper));
        } else if (vub[iColumn]) {
          numberTightened++;
          numberTightenedByProbing++;
          if (!solver->isInteger(iColumn)) {
            // relax
            newLower = CoinMax(lower[iColumn],
              newLower
                - 1.0e-5 * (fabs(lower[iColumn]) + 1));
            newUpper = CoinMin(upper[iColumn],
              newUpper
                + 1.0e-5 * (fabs(upper[iColumn]) + 1));
          }
          solver->setColLower(iColumn, newLower);
          solver->setColUpper(iColumn, newUpper);
        }
      }
    }
  }
  CoinWarmStart *ws = solver->getWarmStart();
  double *solution = new double[numberColumns];
  memcpy(solution, solver->getColSolution(), numberColumns * sizeof(double));
  for (iColumn = 0; iColumn < numberColumns; iColumn++)
    solver->setObjCoeff(iColumn, 0.0);
  //solver->messageHandler()->setLogLevel(2);
  for (iVub = 0; iVub < numberSolves; iVub++) {
    iColumn = which[iVub];
    int iTry;
    for (iTry = 0; iTry < 2; iTry++) {
      double saveUpper = solver->getColUpper()[iColumn];
      double saveLower = solver->getColLower()[iColumn];
      double value;
      if (iTry == 1) {
        // try all way up
        solver->setObjCoeff(iColumn, -1.0);
      } else {
        // try all way down
        solver->setObjCoeff(iColumn, 1.0);
      }
      solver->initialSolve();
      setPointers(continuousSolver_);
      value = solver->getColSolution()[iColumn];
      bool change = false;
      if (iTry == 1) {
        if (value < saveUpper - 1.0e-4) {
          if (solver->isInteger(iColumn)) {
            value = floor(value + 0.00001);
          } else {
            // relax a bit
            value = CoinMin(saveUpper, value + 1.0e-8 * (fabs(saveUpper) + 1));
          }
          if (value - saveLower < 1.0e-7)
            value = saveLower; // make sure exactly same
          solver->setColUpper(iColumn, value);
          saveUpper = value;
          change = true;
        }
      } else {
        if (value > saveLower + 1.0e-4) {
          if (solver->isInteger(iColumn)) {
            value = ceil(value - 0.00001);
          } else {
            // relax a bit
            value = CoinMax(saveLower, value - 1.0e-8 * (fabs(saveLower) + 1));
          }
          if (saveUpper - value < 1.0e-7)
            value = saveUpper; // make sure exactly same
          solver->setColLower(iColumn, value);
          saveLower = value;
          change = true;
        }
      }
      solver->setObjCoeff(iColumn, 0.0);
      if (change) {
        if (saveUpper == saveLower)
          numberFixed++;
        else
          numberTightened++;
        int saveFixed = numberFixed;

        int jColumn;
        if (generator) {
          // Probing - return tight column bounds
          cuts = OsiCuts();
          CglTreeInfo info;
          generator->generateCutsAndModify(*solver, cuts, &info);
          const double *tightLower = generator->tightLower();
          const double *lower = solver->getColLower();
          const double *tightUpper = generator->tightUpper();
          const double *upper = solver->getColUpper();
          for (jColumn = 0; jColumn < numberColumns; jColumn++) {
            double newUpper = tightUpper[jColumn];
            double newLower = tightLower[jColumn];
            if (newUpper < upper[jColumn] - 1.0e-8 * (fabs(upper[jColumn]) + 1) || newLower > lower[jColumn] + 1.0e-8 * (fabs(lower[jColumn]) + 1)) {
              if (newUpper < newLower) {
                fprintf(stderr, "Problem is infeasible\n");
                delete[] solution;
                return false;
              }
              if (newUpper == newLower) {
                numberFixed++;
                numberFixedByProbing++;
                solver->setColLower(jColumn, newLower);
                solver->setColUpper(jColumn, newUpper);
              } else if (vub[jColumn]) {
                numberTightened++;
                numberTightenedByProbing++;
                if (!solver->isInteger(jColumn)) {
                  // relax
                  newLower = CoinMax(lower[jColumn],
                    newLower
                      - 1.0e-8 * (fabs(lower[jColumn]) + 1));
                  newUpper = CoinMin(upper[jColumn],
                    newUpper
                      + 1.0e-8 * (fabs(upper[jColumn]) + 1));
                }
                solver->setColLower(jColumn, newLower);
                solver->setColUpper(jColumn, newUpper);
              }
            }
          }
        }
        if (numberFixed > saveFixed) {
          // original solution may not be feasible
          // go back to true costs to solve if exists
          if (objective) {
            for (jColumn = 0; jColumn < numberColumns; jColumn++)
              solver->setObjCoeff(jColumn, objective[jColumn]);
          }
          solver->setColSolution(solution);
          solver->setWarmStart(ws);
          solver->resolve();
          if (!solver->isProvenOptimal()) {
            fprintf(stderr, "Problem is infeasible\n");
            if (vub)
                delete[] vub;
            return false;
          }
          delete ws;
          ws = solver->getWarmStart();
          memcpy(solution, solver->getColSolution(),
            numberColumns * sizeof(double));
          for (jColumn = 0; jColumn < numberColumns; jColumn++)
            solver->setObjCoeff(jColumn, 0.0);
        }
      }
      solver->setColSolution(solution);
      solver->setWarmStart(ws);
    }
    if (iVub % printFrequency == 0)
      handler_->message(CBC_VUB_PASS, messages_)
        << iVub + 1 << numberFixed << numberTightened
        << CoinMessageEol;
  }
  handler_->message(CBC_VUB_END, messages_)
    << numberFixed << numberTightened
    << CoinMessageEol;
  delete ws;
  delete[] solution;
  // go back to true costs to solve if exists
  if (objective) {
    for (iColumn = 0; iColumn < numberColumns; iColumn++)
      solver_->setObjCoeff(iColumn, objective[iColumn]);
    delete[] objective;
  }
  delete[] vub;
  if (generator) {
    /*printf("Probing fixed %d and tightened %d\n",
           numberFixedByProbing,
           numberTightenedByProbing);*/
    if (generator_[iGen]->howOften() == -1 && (numberFixedByProbing + numberTightenedByProbing) * 5 > (numberFixed + numberTightened))
      generator_[iGen]->setHowOften(1000000 + 1);
    generator->setMaxPass(save[0]);
    generator->setMaxProbe(save[1]);
    generator->setMaxLook(save[2]);
    generator->setRowCuts(save[3]);
  }

  if (solver != solver_) {
    // move bounds across
    const double *lower = solver->getColLower();
    const double *upper = solver->getColUpper();
    const double *lowerOrig = solver_->getColLower();
    const double *upperOrig = solver_->getColUpper();
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      solver_->setColLower(iColumn, CoinMax(lower[iColumn], lowerOrig[iColumn]));
      solver_->setColUpper(iColumn, CoinMin(upper[iColumn], upperOrig[iColumn]));
    }
    delete solver;
  }
  setCutoff(saveCutoff);
  return true;
}
// Pass in Message handler (not deleted at end)
void CbcModel::passInMessageHandler(CoinMessageHandler *handler)
{
  if (defaultHandler_) {
    delete handler_;
    handler_ = NULL;
  }
  defaultHandler_ = false;
  handler_ = handler;
  if (solver_)
    solver_->passInMessageHandler(handler);
  if (continuousSolver_)
    continuousSolver_->passInMessageHandler(handler);
  if (referenceSolver_)
    referenceSolver_->passInMessageHandler(handler);
}
void CbcModel::passInTreeHandler(CbcTree &tree)
{
  delete tree_;
  tree_ = tree.clone();
}
// Make sure region there
void CbcModel::reserveCurrentSolution(const double *solution)
{
  int numberColumns = getNumCols();
  if (!currentSolution_)
    currentSolution_ = new double[numberColumns];
  testSolution_ = currentSolution_;
  if (solution)
    memcpy(currentSolution_, solution, numberColumns * sizeof(double));
}
/* For passing in an CbcModel to do a sub Tree (with derived tree handlers).
   Passed in model must exist for duration of branch and bound
*/
void CbcModel::passInSubTreeModel(CbcModel &model)
{
  subTreeModel_ = &model;
}
// For retrieving a copy of subtree model with given OsiSolver or NULL
CbcModel *
CbcModel::subTreeModel(OsiSolverInterface *solver) const
{
  const CbcModel *subModel = subTreeModel_;
  if (!subModel)
    subModel = this;
  // Get new copy
  CbcModel *newModel = new CbcModel(*subModel);
  if (solver)
    newModel->assignSolver(solver);
  return newModel;
}
//#############################################################################
// Set/Get Application Data
// This is a pointer that the application can store into and retrieve
// from the solverInterface.
// This field is the application to optionally define and use.
//#############################################################################

void CbcModel::setApplicationData(void *appData)
{
  appData_ = appData;
}
//-----------------------------------------------------------------------------
void *CbcModel::getApplicationData() const
{
  return appData_;
}
// Set a pointer to a row cut which will be added instead of normal branching.
void CbcModel::setNextRowCut(const OsiRowCut &cut)
{
  nextRowCut_ = new OsiRowCut(cut);
  nextRowCut_->setEffectiveness(COIN_DBL_MAX); // mark so will always stay
}
// Just update objectiveValue
void CbcModel::setBestObjectiveValue(double objectiveValue)
{
  bestObjective_ = objectiveValue;
}
double
CbcModel::getBestPossibleObjValue() const
{
  return CoinMin(bestPossibleObjective_, bestObjective_) * solver_->getObjSense();
}
// Make given rows (L or G) into global cuts and remove from lp
void CbcModel::makeGlobalCuts(int number, const int *which)
{
  const double *rowLower = solver_->getRowLower();
  const double *rowUpper = solver_->getRowUpper();

  int numberRows = solver_->getNumRows();

  // Row copy
  const double *elementByRow = solver_->getMatrixByRow()->getElements();
  const int *column = solver_->getMatrixByRow()->getIndices();
  const CoinBigIndex *rowStart = solver_->getMatrixByRow()->getVectorStarts();
  const int *rowLength = solver_->getMatrixByRow()->getVectorLengths();

  // Not all rows may be good so we need new array
  int *whichDelete = new int[numberRows];
  int nDelete = 0;
  for (int i = 0; i < number; i++) {
    int iRow = which[i];
    if (iRow >= 0 && iRow < numberRows) {
      if (rowLower[iRow] < -1.0e20 || rowUpper[iRow] > 1.0e20) {
        whichDelete[nDelete++] = iRow;
        OsiRowCut thisCut;
        thisCut.setLb(rowLower[iRow]);
        thisCut.setUb(rowUpper[iRow]);
        CoinBigIndex start = rowStart[iRow];
        thisCut.setRow(rowLength[iRow], column + start, elementByRow + start, false);
        thisCut.setGloballyValid(true);
        globalCuts_.addCutIfNotDuplicate(thisCut);
      }
    }
  }
  if (nDelete)
    solver_->deleteRows(nDelete, whichDelete);
  delete[] whichDelete;
}
// Make given cut into a global cut
int CbcModel::makeGlobalCut(const OsiRowCut *cut)
{
  if (cut->row().getNumElements() > 1 - 1) {
    OsiRowCut newCut(*cut);
    newCut.setGloballyValidAsInteger(2);
    newCut.mutableRow().setTestForDuplicateIndex(false);
    return globalCuts_.addCutIfNotDuplicate(newCut, 1);
  } else {
    assert(cut->row().getNumElements() == 1);
    int iColumn = cut->row().getIndices()[0];
    double value = cut->row().getElements()[0];
    double lb = cut->lb();
    double ub = cut->ub();
    if (value > 0) {
      if (lb > -COIN_DBL_MAX)
        lb /= value;
      if (ub < COIN_DBL_MAX)
        ub /= value;
    } else {
      double saveUb = ub;
      if (lb > -COIN_DBL_MAX)
        ub = lb / value;
      else
        ub = COIN_DBL_MAX;
      if (saveUb < COIN_DBL_MAX)
        lb = saveUb / value;
      else
        lb = -COIN_DBL_MAX;
    }
#if PRINT_CONFLICT == 0
    if (handler_->logLevel() > 1) {
#endif
      printf("Conflict cut at depth %d (%d elements)\n",
        currentDepth_, cut->row().getNumElements());
      cut->print();
#if PRINT_CONFLICT == 0
    }
#endif
    const double *lower;
    const double *upper;
    if (topOfTree_) {
      lower = topOfTree_->lower();
      upper = topOfTree_->upper();
      lb = CoinMax(lb, lower[iColumn]);
      topOfTree_->setColLower(iColumn, lb);
      ub = CoinMin(ub, upper[iColumn]);
      topOfTree_->setColUpper(iColumn, ub);
    } else {
      lower = solver_->getColLower();
      upper = solver_->getColUpper();
      lb = CoinMax(lb, lower[iColumn]);
      solver_->setColLower(iColumn, lb);
      ub = CoinMin(ub, upper[iColumn]);
      solver_->setColUpper(iColumn, ub);
    }
    return 1;
  }
}
// Make given cut into a global cut
int CbcModel::makeGlobalCut(const OsiRowCut &cut)
{
  OsiRowCut newCut(cut);
  newCut.setGloballyValid(true);
  newCut.mutableRow().setTestForDuplicateIndex(false);
  return globalCuts_.addCutIfNotDuplicate(newCut);
}
// Make given column cut into a global cut
void CbcModel::makeGlobalCut(const OsiColCut *cut)
{
  const double *lower;
  const double *upper;
  if (topOfTree_) {
    lower = topOfTree_->lower();
    upper = topOfTree_->upper();
  } else {
    lower = solver_->getColLower();
    upper = solver_->getColUpper();
  }
  int nLower = cut->lbs().getNumElements();
  const int *indexLower = cut->lbs().getIndices();
  const double *boundLower = cut->lbs().getElements();
  for (int i = 0; i < nLower; i++) {
    int iColumn = indexLower[i];
    double newValue = CoinMax(lower[iColumn], boundLower[iColumn]);
    if (topOfTree_)
      topOfTree_->setColLower(iColumn, newValue);
    else
      solver_->setColLower(iColumn, newValue);
  }
  int nUpper = cut->ubs().getNumElements();
  const int *indexUpper = cut->ubs().getIndices();
  const double *boundUpper = cut->ubs().getElements();
  for (int i = 0; i < nUpper; i++) {
    int iColumn = indexUpper[i];
    double newValue = CoinMin(upper[iColumn], boundUpper[iColumn]);
    if (topOfTree_)
      topOfTree_->setColUpper(iColumn, newValue);
    else
      solver_->setColUpper(iColumn, newValue);
  }
}
// Make given column cut into a global cut
void CbcModel::makeGlobalCut(const OsiColCut &cut)
{
  const double *lower;
  const double *upper;
  if (topOfTree_) {
    lower = topOfTree_->lower();
    upper = topOfTree_->upper();
  } else {
    lower = solver_->getColLower();
    upper = solver_->getColUpper();
  }
  int nLower = cut.lbs().getNumElements();
  const int *indexLower = cut.lbs().getIndices();
  const double *boundLower = cut.lbs().getElements();
  for (int i = 0; i < nLower; i++) {
    int iColumn = indexLower[i];
    double newValue = CoinMax(lower[iColumn], boundLower[iColumn]);
    if (topOfTree_)
      topOfTree_->setColLower(iColumn, newValue);
    else
      solver_->setColLower(iColumn, newValue);
  }
  int nUpper = cut.ubs().getNumElements();
  const int *indexUpper = cut.ubs().getIndices();
  const double *boundUpper = cut.ubs().getElements();
  for (int i = 0; i < nUpper; i++) {
    int iColumn = indexUpper[i];
    double newValue = CoinMin(upper[iColumn], boundUpper[iColumn]);
    if (topOfTree_)
      topOfTree_->setColUpper(iColumn, newValue);
    else
      solver_->setColUpper(iColumn, newValue);
  }
}
// Make partial cut into a global cut and save
void CbcModel::makePartialCut(const OsiRowCut *partialCut,
  const OsiSolverInterface *solver)
{
  // get greedy cut
  double bSum = partialCut->lb();
  assert(bSum < 0.0);
  if (!solver)
    solver = solver_;
  int nConflict = partialCut->row().getNumElements();
  const int *column = partialCut->row().getIndices();
  const double *element = partialCut->row().getElements();
  double *originalLower = topOfTree_->mutableLower();
  const double *columnLower = solver->getColLower();
  double *originalUpper = topOfTree_->mutableUpper();
  const double *columnUpper = solver->getColUpper();
  int nC = nConflict;
  while (nConflict) {
    int iColumn = column[nConflict - 1];
    double farkasValue = element[nConflict - 1];
    double change;
    if (farkasValue > 0.0) {
      change = farkasValue * (originalUpper[iColumn] - columnUpper[iColumn]);
    } else {
      change = farkasValue * (originalLower[iColumn] - columnLower[iColumn]);
    }
    if (bSum + change > -1.0e-4)
      break;
    nConflict--;
    bSum += change;
  }
  OsiRowCut newCut;
  newCut.setUb(COIN_DBL_MAX);
  double lo = 1.0;
  double *values = new double[nConflict];
  for (int i = 0; i < nConflict; i++) {
    int iColumn = column[i];
    if (originalLower[iColumn] == columnLower[iColumn]) {
      // must be at least one higher
      values[i] = 1.0;
      lo += originalLower[iColumn];
    } else {
      // must be at least one lower
      values[i] = -1.0;
      lo -= originalUpper[iColumn];
    }
  }
  newCut.setLb(lo);
  newCut.setRow(nConflict, column, values);
  printf("CUTa has %d (started at %d) - final bSum %g - depth %d\n", nConflict, nC, bSum, currentDepth_);
  if (nConflict > 1) {
    if ((specialOptions_ & 1) != 0) {
      const OsiRowCutDebugger *debugger = continuousSolver_->getRowCutDebugger();
      if (debugger) {
        if (debugger->invalidCut(newCut)) {
          continuousSolver_->applyRowCuts(1, &newCut);
          continuousSolver_->writeMps("bad");
        }
        CoinAssert(!debugger->invalidCut(newCut));
      }
    }
    newCut.setGloballyValidAsInteger(2);
    newCut.mutableRow().setTestForDuplicateIndex(false);
    globalCuts_.addCutIfNotDuplicate(newCut);
  } else {
    // change bounds
    int iColumn = column[0];
    if (values[0] < 0.0) {
      // change upper bound
      double newUpper = -lo;
      assert(newUpper < originalUpper[iColumn]);
      printf("Changing upper bound on %d from %g to %g\n",
        iColumn, originalUpper[iColumn], newUpper);
      originalUpper[iColumn] = newUpper;
    } else {
      // change lower bound
      double newLower = lo;
      assert(newLower > originalLower[iColumn]);
      printf("Changing lower bound on %d from %g to %g\n",
        iColumn, originalLower[iColumn], newLower);
      originalLower[iColumn] = newLower;
    }
  }
  // add to partial cuts
  if (globalConflictCuts_) {
    globalConflictCuts_->addCutIfNotDuplicateWhenGreedy(*partialCut, 2);
  }
  delete[] values;
}
// Make partial cuts into global cuts
void CbcModel::makeGlobalCuts()
{
}
void CbcModel::setNodeComparison(CbcCompareBase *compare)
{
  delete nodeCompare_;
  nodeCompare_ = compare->clone();
}
void CbcModel::setNodeComparison(CbcCompareBase &compare)
{
  delete nodeCompare_;
  nodeCompare_ = compare.clone();
}
void CbcModel::setProblemFeasibility(CbcFeasibilityBase *feasibility)
{
  delete problemFeasibility_;
  problemFeasibility_ = feasibility->clone();
}
void CbcModel::setProblemFeasibility(CbcFeasibilityBase &feasibility)
{
  delete problemFeasibility_;
  problemFeasibility_ = feasibility.clone();
}
// Set the strategy. Clones
void CbcModel::setStrategy(CbcStrategy &strategy)
{
  delete strategy_;
  strategy_ = strategy.clone();
}
// Increases usedInSolution for nonzeros
void CbcModel::incrementUsed(const double *solution)
{
  if (usedInSolution_) {
    // might as well mark all including continuous
    int numberColumns = solver_->getNumCols();
    for (int i = 0; i < numberColumns; i++) {
      if (solution[i])
        usedInSolution_[i]++;
    }
  }
}
// Are there numerical difficulties (for initialSolve) ?
bool CbcModel::isInitialSolveAbandoned() const
{
  if (status_ != -1) {
    return false;
  } else {
    return solver_->isAbandoned();
  }
}
// Is optimality proven (for initialSolve) ?
bool CbcModel::isInitialSolveProvenOptimal() const
{
  if (status_ != -1) {
    return fabs(originalContinuousObjective_) < 1.0e50;
  } else {
    return solver_->isProvenOptimal();
  }
}
// Is primal infeasiblity proven (for initialSolve) ?
bool CbcModel::isInitialSolveProvenPrimalInfeasible() const
{
  if (status_ != -1) {
    if (status_ == 0 && secondaryStatus_ == 7)
      return false;
    else
      return originalContinuousObjective_ >= 1.0e50;
  } else {
    return solver_->isProvenPrimalInfeasible();
  }
}
// Is dual infeasiblity proven (for initialSolve) ?
bool CbcModel::isInitialSolveProvenDualInfeasible() const
{
  if (status_ != -1) {
    if (status_ == 0 && secondaryStatus_ == 7)
      return true;
    else
      return false;
  } else {
    return solver_->isProvenDualInfeasible();
  }
}
// Set pointers for speed
void CbcModel::setPointers(const OsiSolverInterface *solver)
{
  /// Pointer to array[getNumCols()] (for speed) of column lower bounds
  cbcColLower_ = solver_->getColLower();
  /// Pointer to array[getNumCols()] (for speed) of column upper bounds
  cbcColUpper_ = solver_->getColUpper();
  /// Pointer to array[getNumRows()] (for speed) of row lower bounds
  cbcRowLower_ = solver_->getRowLower();
  /// Pointer to array[getNumRows()] (for speed) of row upper bounds
  cbcRowUpper_ = solver_->getRowUpper();
  /// Pointer to array[getNumCols()] (for speed) of primal solution vector
  cbcColSolution_ = solver_->getColSolution();
  /// Pointer to array[getNumRows()] (for speed) of dual prices
  cbcRowPrice_ = solver_->getRowPrice();
  /// Get a pointer to array[getNumCols()] (for speed) of reduced costs
  if (solverCharacteristics_ && solverCharacteristics_->reducedCostsAccurate())
    cbcReducedCost_ = solver_->getReducedCost();
  else
    cbcReducedCost_ = NULL;
  /// Pointer to array[getNumRows()] (for speed) of row activity levels.
  cbcRowActivity_ = solver_->getRowActivity();
  dblParam_[CbcCurrentObjectiveValue] = solver->getObjValue();
  dblParam_[CbcCurrentMinimizationObjectiveValue] = dblParam_[CbcCurrentObjectiveValue] * dblParam_[CbcOptimizationDirection];
}

/*
  Delete any existing handler and create a clone of the one supplied.
*/
void CbcModel::passInEventHandler(const CbcEventHandler *eventHandler)
{
  delete eventHandler_;
  eventHandler_ = NULL;
  if (eventHandler) {
    eventHandler_ = eventHandler->clone();
    eventHandler_->setModel(this);
  }
}

/*
  CbcEventHandler* CbcModel::eventHandler is inlined in CbcModel.hpp.
*/

// Encapsulates solver resolve
int CbcModel::resolve(OsiSolverInterface *solver)
{
  numberSolves_++;
#ifdef COIN_HAS_CLP
  OsiClpSolverInterface *clpSolver
    = dynamic_cast< OsiClpSolverInterface * >(solver);
#endif
#ifdef CLIQUE_ANALYSIS
  if (probingInfo_ && currentDepth_ > 0) {
    int nFix = probingInfo_->fixColumns(*solver);
#ifdef SWITCH_VARIABLES
    if (nFix > 0)
      fixAssociated(solver_, 0);
#endif
    if (nFix < 0) {
#ifdef COIN_HAS_CLP
      if (clpSolver)
        clpSolver->getModelPtr()->setProblemStatus(1);
#endif
      return 0;
    }
  }
#endif
#ifdef COIN_HAS_CLP
  if (clpSolver) {
    /*bool takeHint;
        OsiHintStrength strength;
        bool gotHint = (clpSolver->getHintParam(OsiDoDualInResolve,takeHint,strength));
        assert (gotHint);
        int algorithm=-1;
        if (strength!=OsiHintIgnore)
          algorithm = takeHint ? -1 : 1;
          assert (algorithm==-1);*/
    //clpSolver->setHintParam(OsiDoDualInResolve,true,OsiHintTry);
    ClpSimplex *clpSimplex = clpSolver->getModelPtr();
    int save = clpSimplex->specialOptions();
    if ((moreSpecialOptions_ & 8388608) == 0)
      clpSimplex->setSpecialOptions(save | 0x11000000); // say is Cbc (and in branch and bound)
    else
      clpSimplex->setSpecialOptions(save | 0x11200000); // say is Cbc (and in branch and bound - but save ray)
    int save2 = clpSolver->specialOptions();
    if (false && (save2 & 2048) == 0) {
      // see if worthwhile crunching
      int nFixed = 0;
      const double *columnLower = clpSimplex->columnLower();
      const double *columnUpper = clpSimplex->columnUpper();
      for (int i = 0; i < numberIntegers_; i++) {
        int iColumn = integerVariable_[i];
        if (columnLower[iColumn] == columnUpper[iColumn])
          nFixed++;
      }
      if (nFixed * 20 < clpSimplex->numberColumns()) {
        double d = nFixed;
        printf("%d fixed out of %d - ratio %g\n",
          nFixed,
          clpSimplex->numberColumns(),
          d / clpSimplex->numberColumns());
        clpSolver->setSpecialOptions(save2 | 2048);
      }
    }
#ifdef CHECK_KNOWN_SOLUTION
    bool onOptimalPath = false;
    if ((specialOptions_ & 1) != 0) {
      const OsiRowCutDebugger *debugger = solver_->getRowCutDebugger();
      if (debugger) {
        onOptimalPath = true;
        printf("On optimal path before resolve\n");
      }
    }
#endif
    clpSolver->resolve();
#ifdef CHECK_RAY
    static int nSolves = 0;
    static int nInfSolves = 0;
    static int nRays = 0;
    nSolves++;
    if (!parentModel_ && clpSolver->getModelPtr()->problemStatus() == 1
      && (clpSolver->getModelPtr()->specialOptions() & 32) != 0) {
      nInfSolves++;
      if (clpSolver->getModelPtr()->infeasibilityRay())
        nRays++;
    }
    if ((nSolves % 1000) == 0)
      printf("ZZ %d solves, %d infeasible %d rays\n",
        nSolves, nInfSolves, nRays);
#endif
#ifdef CHECK_KNOWN_SOLUTION
    if ((specialOptions_ & 1) != 0 && onOptimalPath) {
      const OsiRowCutDebugger *debugger = solver_->getRowCutDebugger();
      if (debugger) {
        printf("On optimal path after resolve\n");
      } else {
        solver_->writeMpsNative("badSolve.mps", NULL, NULL, 2);
        printf("NOT on optimal path after resolve\n");
      }
    }
#endif
    if (!numberNodes_) {
#if 0
	  // make sure really feasible when not scaled
	  if (clpSimplex->secondaryStatus()==2) {
	    messageHandler()->message(CBC_FPUMP1, messages())
	      << "Continuous solution cleaned for scaling" 
	      << CoinMessageEol ;
	    // overkill but stop odd things happeneing
	    bool takeHint1,takeHint2;
	    OsiHintStrength strength;
	    clpSolver->getHintParam(OsiDoScale,takeHint1,strength);
	    clpSolver->getHintParam(OsiDoPresolveInInitial,takeHint2,strength);
	    //assert (takeHint);
	    clpSolver->setHintParam(OsiDoScale,false,OsiHintTry);
	    clpSolver->setHintParam(OsiDoPresolveInInitial,false,OsiHintTry);
	    int saveFlag=clpSimplex->scalingFlag();
	    clpSimplex->scaling(0);
	    clpSolver->initialSolve();
	    clpSolver->setHintParam(OsiDoScale,takeHint1,OsiHintTry);
	    clpSolver->setHintParam(OsiDoPresolveInInitial,takeHint2,OsiHintTry);
	    clpSimplex->scaling(saveFlag);
	  }
#endif
      double error = CoinMax(clpSimplex->largestDualError(),
        clpSimplex->largestPrimalError());
      if (error > 1.0e-2 || !clpSolver->isProvenOptimal()) {
#if CBC_USEFUL_PRINTING > 1
        printf("Problem was %s largest dual error %g largest primal %g - safer cuts\n",
          clpSolver->isProvenOptimal() ? "optimal" : "!infeasible",
          clpSimplex->largestDualError(),
          clpSimplex->largestPrimalError());
#endif
        if (!clpSolver->isProvenOptimal()) {
          // check if proven infeasible i.e. bad bounds
          int numberColumns = clpSolver->getNumCols();
          const double *columnLower = clpSolver->getColLower();
          const double *columnUpper = clpSolver->getColUpper();
          bool provenInfeasible = false;
          for (int i = 0; i < numberColumns; i++) {
            if (columnLower[i] > columnUpper[i]) {
              provenInfeasible = true;
            }
          }
          if (!provenInfeasible) {
            clpSolver->setSpecialOptions(save2 | 2048);
            clpSimplex->allSlackBasis(true);
            clpSolver->resolve();
            if (!clpSolver->isProvenOptimal()) {
              bool takeHint;
              OsiHintStrength strength;
              clpSolver->getHintParam(OsiDoDualInResolve, takeHint, strength);
              clpSolver->setHintParam(OsiDoDualInResolve, false, OsiHintDo);
              clpSolver->resolve();
              clpSolver->setHintParam(OsiDoDualInResolve, takeHint, strength);
            }
          }
        }
        // make cuts safer
        for (int iCutGenerator = 0; iCutGenerator < numberCutGenerators_; iCutGenerator++) {
          CglCutGenerator *generator = generator_[iCutGenerator]->generator();
          CglGomory *cgl1 = dynamic_cast< CglGomory * >(generator);
          if (cgl1) {
            cgl1->setLimitAtRoot(cgl1->getLimit());
          }
          CglTwomir *cgl2 = dynamic_cast< CglTwomir * >(generator);
          if (cgl2) {
            generator_[iCutGenerator]->setHowOften(-100);
          }
        }
      }
    }
    clpSolver->setSpecialOptions(save2);
#if CBC_USEFUL_PRINTING > 1
    if (clpSimplex->numberIterations() > 1000)
      printf("node %d took %d iterations\n", numberNodes_, clpSimplex->numberIterations());
#endif
    clpSimplex->setSpecialOptions(save);
    if (clpSimplex->status() == 4)
      clpSimplex->setProblemStatus(1);
  } else {
    solver->resolve();
  }
#else
  solver->resolve();
#endif
#ifdef SWITCH_VARIABLES
  if (solver_->isProvenOptimal()) {
    int nBad = checkAssociated(solver_, solver_->getColSolution(), 0);
    if (nBad)
      checkAssociated(solver_, solver_->getColSolution(), 1);
  }
#endif
  return solver->isProvenOptimal() ? 1 : 0;
}
#ifdef CLP_RESOLVE
// Special purpose resolve
int CbcModel::resolveClp(OsiClpSolverInterface *clpSolver, int type)
{
  numberSolves_++;
  ClpSimplex *clpSimplex = clpSolver->getModelPtr();
  int save = clpSimplex->specialOptions();
  clpSimplex->setSpecialOptions(save | 0x11000000); // say is Cbc (and in branch and bound)
  int save2 = clpSolver->specialOptions();
  clpSolver->resolve();
  if (!numberNodes_) {
    double error = CoinMax(clpSimplex->largestDualError(),
      clpSimplex->largestPrimalError());
    if (error > 1.0e-2 || !clpSolver->isProvenOptimal()) {
#if CBC_USEFUL_PRINTING > 1
      printf("Problem was %s largest dual error %g largest primal %g - safer cuts\n",
        clpSolver->isProvenOptimal() ? "optimal" : "!infeasible",
        clpSimplex->largestDualError(),
        clpSimplex->largestPrimalError());
#endif
      if (!clpSolver->isProvenOptimal()) {
        clpSolver->setSpecialOptions(save2 | 2048);
        clpSimplex->allSlackBasis(true);
        clpSolver->resolve();
        if (!clpSolver->isProvenOptimal()) {
          bool takeHint;
          OsiHintStrength strength;
          clpSolver->getHintParam(OsiDoDualInResolve, takeHint, strength);
          clpSolver->setHintParam(OsiDoDualInResolve, false, OsiHintDo);
          clpSolver->resolve();
          clpSolver->setHintParam(OsiDoDualInResolve, takeHint, strength);
        }
      }
      // make cuts safer
      for (int iCutGenerator = 0; iCutGenerator < numberCutGenerators_; iCutGenerator++) {
        CglCutGenerator *generator = generator_[iCutGenerator]->generator();
        CglGomory *cgl1 = dynamic_cast< CglGomory * >(generator);
        if (cgl1) {
          cgl1->setLimitAtRoot(cgl1->getLimit());
        }
        CglTwomir *cgl2 = dynamic_cast< CglTwomir * >(generator);
        if (cgl2) {
          generator_[iCutGenerator]->setHowOften(-100);
        }
      }
    }
  }
  clpSolver->setSpecialOptions(save2);
#if CBC_USEFUL_PRINTING > 1
  if (clpSimplex->numberIterations() > 1000)
    printf("node %d took %d iterations\n", numberNodes_, clpSimplex->numberIterations());
#endif
  if (type == 0 && clpSolver->isProvenOptimal()) {
    ClpSimplex newModel(*clpSimplex);
    newModel.primal();
    int numberColumns = newModel.numberColumns();
    int numberRows = newModel.numberRows();
    double *obj = new double[numberColumns];
    int *which = new int[numberColumns];
    const double *solution = clpSimplex->primalColumnSolution();
    double rhs = 1.0e-8;
    int numberObj = 0;
    double integerTolerance = getDblParam(CbcIntegerTolerance);
    double *objective = newModel.objective();
    for (int i = 0; i < numberColumns; i++) {
      if (objective[i]) {
        rhs += objective[i] * solution[i];
        obj[numberObj] = objective[i];
        which[numberObj++] = i;
        objective[i] = 0.0;
      }
    }
    if (numberObj) {
      newModel.addRow(numberObj, which, obj, -COIN_DBL_MAX, rhs);
    }
    delete[] obj;
    delete[] which;
    double *lower = newModel.columnLower();
    double *upper = newModel.columnUpper();
    int numberInf = 0;
    int numberLb = 0;
    int numberUb = 0;
    int numberInt = 0;
    double sumInf = 0.0;
    for (int i = 0; i < numberIntegers_; i++) {
      int iSequence = integerVariable_[i];
      double value = solution[iSequence];
      value = CoinMax(value, lower[iSequence]);
      value = CoinMin(value, upper[iSequence]);
      double nearest = floor(value + 0.5);
      if (value < lower[iSequence] + integerTolerance) {
        objective[iSequence] = 1.0;
        numberLb++;
      } else if (value > upper[iSequence] - integerTolerance) {
        objective[iSequence] = -1.0;
        numberUb++;
      } else if (fabs(value - nearest) <= integerTolerance) {
        // fix??
        lower[iSequence] = nearest;
        upper[iSequence] = nearest;
        numberInt++;
      } else {
        lower[iSequence] = floor(value);
        upper[iSequence] = ceil(value);
        if (value > nearest) {
          objective[iSequence] = 1.0;
          sumInf += value - nearest;
        } else {
          objective[iSequence] = -1.0;
          sumInf -= value - nearest;
        }
        numberInf++;
      }
    }
    printf("XX %d inf (sum %g), %d at lb %d at ub %d other integer\n",
      numberInf, sumInf, numberLb, numberUb, numberInt);
    if (numberInf) {
      newModel.primal(1);
      if (!newModel.isProvenOptimal()) {
        printf("not optimal - scaling issue - switch off\n");
        clpSimplex->setSpecialOptions(save);
        if (clpSimplex->status() == 4)
          clpSimplex->setProblemStatus(1);
        return clpSolver->isProvenOptimal() ? 1 : 0;
      }
      //newModel.writeMps("bad.mps");
      //assert (newModel.isProvenOptimal());
      printf("%d iterations\n", newModel.numberIterations());
      int numberInf2 = 0;
      int numberLb2 = 0;
      int numberUb2 = 0;
      int numberInt2 = 0;
      double sumInf2 = 0.0;
      const double *solution = newModel.primalColumnSolution();
      const double *lower = clpSimplex->columnLower();
      const double *upper = clpSimplex->columnUpper();
      for (int i = 0; i < numberIntegers_; i++) {
        int iSequence = integerVariable_[i];
        double value = solution[iSequence];
        value = CoinMax(value, lower[iSequence]);
        value = CoinMin(value, upper[iSequence]);
        double nearest = floor(value + 0.5);
        if (value < lower[iSequence] + integerTolerance) {
          numberLb2++;
        } else if (value > upper[iSequence] - integerTolerance) {
          numberUb2++;
        } else if (fabs(value - nearest) <= integerTolerance) {
          numberInt2++;
        } else {
          if (value > nearest) {
            sumInf2 += value - nearest;
          } else {
            sumInf2 -= value - nearest;
          }
          numberInf2++;
        }
      }
      printf("XXX %d inf (sum %g), %d at lb %d at ub %d other integer\n",
        numberInf2, sumInf2, numberLb2, numberUb2, numberInt2);
      if (sumInf2 < sumInf * 0.95) {
        printf("XXXX suminf reduced from %g (%d) to %g (%d)\n",
          sumInf, numberInf, sumInf2, numberInf2);
        if (numberObj) {
          newModel.deleteRows(1, &numberRows);
        }
        memcpy(newModel.objective(),
          clpSimplex->objective(),
          numberColumns * sizeof(double));
        memcpy(newModel.columnLower(),
          clpSimplex->columnLower(),
          numberColumns * sizeof(double));
        memcpy(newModel.columnUpper(),
          clpSimplex->columnUpper(),
          numberColumns * sizeof(double));
        newModel.setClpScaledMatrix(NULL);
        newModel.primal(1);
        printf("%d iterations\n", newModel.numberIterations());
        int numberInf3 = 0;
        int numberLb3 = 0;
        int numberUb3 = 0;
        int numberInt3 = 0;
        double sumInf3 = 0.0;
        const double *solution = newModel.primalColumnSolution();
        const double *lower = clpSimplex->columnLower();
        const double *upper = clpSimplex->columnUpper();
        for (int i = 0; i < numberIntegers_; i++) {
          int iSequence = integerVariable_[i];
          double value = solution[iSequence];
          value = CoinMax(value, lower[iSequence]);
          value = CoinMin(value, upper[iSequence]);
          double nearest = floor(value + 0.5);
          if (value < lower[iSequence] + integerTolerance) {
            numberLb3++;
          } else if (value > upper[iSequence] - integerTolerance) {
            numberUb3++;
          } else if (fabs(value - nearest) <= integerTolerance) {
            numberInt3++;
          } else {
            if (value > nearest) {
              sumInf3 += value - nearest;
            } else {
              sumInf3 -= value - nearest;
            }
            numberInf3++;
          }
        }
        printf("XXXXX %d inf (sum %g), %d at lb %d at ub %d other integer\n",
          numberInf3, sumInf3, numberLb3, numberUb3, numberInt3);
        if (sumInf3 < sumInf * 0.95) {
          memcpy(clpSimplex->primalColumnSolution(),
            newModel.primalColumnSolution(),
            numberColumns * sizeof(double));
          memcpy(clpSimplex->dualColumnSolution(),
            newModel.dualColumnSolution(),
            numberColumns * sizeof(double));
          memcpy(clpSimplex->primalRowSolution(),
            newModel.primalRowSolution(),
            numberRows * sizeof(double));
          memcpy(clpSimplex->dualRowSolution(),
            newModel.dualRowSolution(),
            numberRows * sizeof(double));
          memcpy(clpSimplex->statusArray(),
            newModel.statusArray(),
            (numberColumns + numberRows) * sizeof(unsigned char));
          clpSolver->setWarmStart(NULL);
        }
      }
    }
  }
  clpSimplex->setSpecialOptions(save);
  if (clpSimplex->status() == 4)
    clpSimplex->setProblemStatus(1);
  return clpSolver->isProvenOptimal() ? 1 : 0;
}
#endif
/*!
    \todo It'd be really nice if there were an overload for this method that
	  allowed a separate value for the underlying solver's log level. The
	  overload could be coded to allow an increase in the log level of the
	  underlying solver.

	  It's worth contemplating whether OSI should have a setLogLevel method
	  that's more specific than the hint mechanism.
*/

// Set log level
void CbcModel::setLogLevel(int value)
{
  handler_->setLogLevel(value);
  // Reduce print out in Osi
  if (solver_) {
    int oldLevel = solver_->messageHandler()->logLevel();
    if (value < oldLevel)
      solver_->messageHandler()->setLogLevel(value);
#ifdef COIN_HAS_CLP
    OsiClpSolverInterface *clpSolver
      = dynamic_cast< OsiClpSolverInterface * >(solver_);
    if (clpSolver) {
      ClpSimplex *clpSimplex = clpSolver->getModelPtr();
      oldLevel = clpSimplex->logLevel();
      if (value < oldLevel)
        clpSimplex->setLogLevel(value);
    }
#else // COIN_HAS_CLP
    /*
          For generic OSI solvers, if the new log level is 0, try the
          DoReducePrint hint for emphasis.
        */
    if (value == 0) {
      solver_->setHintParam(OsiDoReducePrint, true, OsiHintDo);
    }
#endif // COIN_HAS_CLP
  }
}

/* Pass in target solution and optional priorities.
   If priorities then >0 means only branch if incorrect
   while <0 means branch even if correct. +1 or -1 are
   highest priority */
void CbcModel::setHotstartSolution(const double *solution, const int *priorities)
{
  if (solution == NULL) {
    delete[] hotstartSolution_;
    hotstartSolution_ = NULL;
    delete[] hotstartPriorities_;
    hotstartPriorities_ = NULL;
  } else {
    int numberColumns = solver_->getNumCols();
    hotstartSolution_ = CoinCopyOfArray(solution, numberColumns);
    hotstartPriorities_ = CoinCopyOfArray(priorities, numberColumns);
    for (int i = 0; i < numberColumns; i++) {
      if (hotstartSolution_[i] == -COIN_DBL_MAX) {
        hotstartSolution_[i] = 0.0;
        hotstartPriorities_[i] += 10000;
      }
      if (solver_->isInteger(i))
        hotstartSolution_[i] = floor(hotstartSolution_[i] + 0.5);
    }
  }
}
// Increment strong info
void CbcModel::incrementStrongInfo(int numberTimes, int numberIterations,
  int numberFixed, bool ifInfeasible)
{
  strongInfo_[0] += numberTimes;
  numberStrongIterations_ += numberIterations;
  strongInfo_[1] += numberFixed;
  if (ifInfeasible)
    strongInfo_[2]++;
}
/* Set objective value in a node.  This is separated out so that
   odd solvers can use.  It may look at extra information in
   solverCharacteriscs_ and will also use bound from parent node
*/
void CbcModel::setObjectiveValue(CbcNode *thisNode, const CbcNode *parentNode) const
{
  double newObjValue = solver_->getObjSense() * solver_->getObjValue();
  // If odd solver take its bound
  if (solverCharacteristics_) {
    newObjValue = CoinMax(newObjValue, solverCharacteristics_->mipBound());
    // Reset bound anyway (no harm if not odd)
    solverCharacteristics_->setMipBound(-COIN_DBL_MAX);
  }
  // If not root then use max of this and parent
  if (parentNode)
    newObjValue = CoinMax(newObjValue, parentNode->objectiveValue());
  thisNode->setObjectiveValue(newObjValue);
}
// Current time since start of branchAndbound
double
CbcModel::getCurrentSeconds() const
{
  if (!useElapsedTime())
    return CoinCpuTime() - getDblParam(CbcStartSeconds);
  else
    return CoinGetTimeOfDay() - getDblParam(CbcStartSeconds);
}
/* Encapsulates choosing a variable -
   anyAction: -2 infeasible
	      -1 round again
	       0 done

   At the point where chooseBranch is called, we've decided that this problem
   will need to be placed in the live set and we need to choose a branching
   variable.

   Parameters:
     newNode:	the node just created for the active subproblem.
     oldNode:	newNode's parent.
     lastws:	oldNode's basis
     lowerBefore, upperBefore: column bound arrays for oldNode
     cuts:	list of cuts added to newNode.

     resolved:	(o)  set to true if newNode is resolved during processing
     branches:	(o) will be filled in with ... ? Null on entry
*/
int CbcModel::chooseBranch(CbcNode *&newNode, int numberPassesLeft,
  CbcNode *oldNode, OsiCuts &cuts,
  bool &resolved, CoinWarmStartBasis *lastws,
  const double *lowerBefore, const double *upperBefore,
  OsiSolverBranch *&branches)
{
  // Set state of search
  /*
      0 - outside CbcNode
      1 - no solutions
      2 - all heuristic solutions
      3 - a solution reached by branching (could be strong)
      4 - no solution but many nodes
         add 10 if depth >= K
    K is currently hardcoded to 8, a few lines below.

    CBCMODEL_DEBUG: Seems like stateOfSearch_ should be 2 if
           numberHeuristicSolutions_ == numberSolutions_.

    */
  stateOfSearch_ = 1;
  if (numberSolutions_ > 0) {
    if (numberHeuristicSolutions_ == numberSolutions_)
      stateOfSearch_ = 3;
    else
      stateOfSearch_ = 3;
  }
  if (numberNodes_ > 2 * numberObjects_ + 1000) {
    stateOfSearch_ = 4;
  }
  //stateOfSearch_=3;
  if (currentNode_ && currentNode_->depth() >= 8)
    stateOfSearch_ += 10;
  //printf("STate %d, %d nodes - parent %c - sol %d %d\n",
  // stateOfSearch_,numberNodes_,parentModel_ ? 'Y' :'N',
  // numberSolutions_,numberHeuristicSolutions_);
  int anyAction = -1;
  resolved = false;
  if (newNode->objectiveValue() >= getCutoff())
    anyAction = -2;
  branches = NULL;
  bool feasible = true;
  int branchingState = -1;
  // Compute "small" change in branch
  int nBranches = intParam_[CbcNumberBranches];
  if (nBranches) {
    double average = dblParam_[CbcSumChange] / static_cast< double >(nBranches);
    dblParam_[CbcSmallChange] = CoinMax(average * 1.0e-5, dblParam_[CbcSmallestChange]);
    dblParam_[CbcSmallChange] = CoinMax(dblParam_[CbcSmallChange], 1.0e-8);
  } else {
    dblParam_[CbcSmallChange] = 1.0e-8;
  }
#ifdef JJF_ZERO
  // Say not on optimal path
  bool onOptimalPath = false;
  if ((specialOptions_ & 1) != 0) {
    /*
          This doesn't work as intended --- getRowCutDebugger will return null
          unless the current feasible solution region includes the optimal solution
          that RowCutDebugger knows. There's no way to tell inactive from off the
          optimal path.
        */
    const OsiRowCutDebugger *debugger = solver_->getRowCutDebugger();
    if (debugger) {
      onOptimalPath = true;
      printf("On optimal path - choose\n");
    }
  }
#endif
  currentNode_ = newNode; // so can be used elsewhere
  // Remember number of rows to restore at the end of the loop
  int saveNumberRows = solver_->getNumRows();
  /*
      Enough preparation. Get down to the business of choosing a branching
      variable.
    */
  while (anyAction == -1) {
    // Set objective value (not so obvious if NLP etc)
    setObjectiveValue(newNode, oldNode);
    //if (numberPassesLeft<=0)
    //branchingState=1;
    /*
          Is there a CbcBranchDecision object installed? Does it specify a
          chooseVariable method? If not, we're using the old (Cbc) side of the branch
          decision hierarchy.  In quick summary, CbcNode::chooseBranch uses strong
          branching on any objects, while CbcNode::chooseDynamicBranch uses dynamic
          branching, but only on simple integers (-3 is the code for punt due to
          complex objects). Serious bugs remain on the Cbc side, particularly in
          chooseDynamicBranch.
        */
    if (!branchingMethod_ || !branchingMethod_->chooseMethod()) {
#ifdef COIN_HAS_CLP
      bool doClp = oldNode && (oldNode->depth() % 2) == 1;
      if (!doCutsNow(1))
        doClp = true;
      //doClp = true;
      int testDepth = 5;
      // Don't do if many iterations per node
      int totalNodes = numberNodes_ + numberExtraNodes_;
      int totalIterations = numberIterations_ + numberExtraIterations_;
      bool diving = false;
      if ((moreSpecialOptions_ & 33554432) != 0) {
        testDepth = COIN_INT_MAX;
        if (oldNode && (oldNode->depth() == -2 || oldNode->depth() == 4))
          diving = true;
      }
      if (totalNodes * 40 < totalIterations || numberNodes_ < 1000) {
        doClp = false;
        //} else if (oldNode&&fastNodeDepth_>=0&&oldNode->depth()>=testDepth&&(specialOptions_&2048)==0) {
        //printf("size %d %d - cuts %d - nodes %d its %d %c\n",solver_->getNumRows(),
        //     solver_->getNumCols(),cuts.sizeRowCuts(),
        //     totalNodes,totalIterations,doClp ? 'Y' : 'N');
      }
      if (oldNode && ((fastNodeDepth_ >= 0 && oldNode->depth() >= testDepth && doClp) || diving) && /*!parentModel_*/ (specialOptions_ & 2048) == 0
        && !cuts.sizeRowCuts()) {
        OsiClpSolverInterface *clpSolver
          = dynamic_cast< OsiClpSolverInterface * >(solver_);
        if (clpSolver) {
          anyAction = newNode->chooseClpBranch(this, oldNode);
	  currentNode_ = NULL;
          if (anyAction != -1)
            break;
        }
      }
#endif
#ifdef COIN_HAS_CLP
      // Deal with funny variables
      if ((moreSpecialOptions2_ & 32768) != 0)
        cleanBounds(solver_, NULL);
      int save = 0;
      OsiClpSolverInterface *clpSolver
        = dynamic_cast< OsiClpSolverInterface * >(solver_);
      if (clpSolver && (moreSpecialOptions_ & 4194304) != 0) {
        ClpSimplex *clpSimplex = clpSolver->getModelPtr();
        save = clpSimplex->specialOptions();
        clpSimplex->setSpecialOptions(save | 0x11200000); // say is Cbc (and in branch and bound - but save ray)
      }
#endif
#ifdef COIN_HAS_NTY
      if (symmetryInfo_) {
        CbcNodeInfo *infoX = oldNode ? oldNode->nodeInfo() : NULL;
        bool worthTrying = false;
        if (infoX) {
          CbcNodeInfo *info = infoX;
          for (int i = 0; i < NTY_BAD_DEPTH; i++) {
            if (!info->parent()) {
              worthTrying = true;
              break;
            }
            info = info->parent();
            if (info->symmetryWorked()) {
              worthTrying = true;
              break;
            }
          }
        } else {
          worthTrying = true;
        }
        if ((moreSpecialOptions2_ & (128 | 256)) == (128 | 256) && currentDepth_ > 5)
          worthTrying = false;
        if (worthTrying) {
          int n = symmetryInfo_->orbitalFixing(solver_);
          if (n) {
#if PRINT_MORE == 0
            if (logLevel() > 1)
              printf("%d orbital fixes\n", n);
#endif
            solver_->resolve();
            if (!isProvenOptimal()) {
              if (logLevel() > 1)
                printf("infeasible after orbital fixing\n");
            }
          }
        }
      }
#endif
      if (numberBeforeTrust_ == 0) {
        anyAction = newNode->chooseBranch(this, oldNode, numberPassesLeft);
      } else {
        anyAction = newNode->chooseDynamicBranch(this, oldNode, branches, numberPassesLeft);
        if (anyAction == -3)
          anyAction = newNode->chooseBranch(this, oldNode, numberPassesLeft); // dynamic did nothing
      }
      currentNode_ = NULL;
#ifdef COIN_HAS_CLP
      if (clpSolver && (moreSpecialOptions_ & 4194304) != 0) {
        ClpSimplex *clpSimplex = clpSolver->getModelPtr();
        clpSimplex->setSpecialOptions(save);
      }
#endif
      /*
              We're on the new (Osi) side of the branching hierarchy.
            */
    } else {
      OsiBranchingInformation usefulInfo = usefulInformation();
      anyAction = newNode->chooseOsiBranch(this, oldNode, &usefulInfo, branchingState);
      currentNode_ = NULL;
      //branchingState=0;
    }
    if (!oldNode) {
      if (numberUpdateItems_) {
        for (int i = 0; i < numberUpdateItems_; i++) {
          CbcObjectUpdateData *update = updateItems_ + i;
          CbcObject *object = dynamic_cast< CbcObject * >(update->object_);
#ifndef NDEBUG
          bool found = false;
          for (int j = 0; j < numberObjects_; j++) {
            if (update->object_ == object_[j]) {
              found = true;
              break;
            }
          }
          assert(found);
#endif
          if (object)
            object->updateInformation(*update);
        }
        numberUpdateItems_ = 0;
      }
    }
    if (solverCharacteristics_ && solverCharacteristics_->solutionAddsCuts() && // we are in some OA based bab
      feasible && (newNode->numberUnsatisfied() == 0) //solution has become integer feasible during strong branching
    ) {
      //in the present case we need to check here integer infeasibility if the node is not fathomed we will have to do the loop
      // again
      //std::cout<<solver_<<std::endl;

      OsiCuts feasCuts;

      for (int i = 0; i < numberCutGenerators_ && (feasCuts.sizeRowCuts() == 0); i++) {
        if (generator_[i]->normal() && (!generator_[i]->needsOptimalBasis() || solver_->basisIsAvailable()))
          generator_[i]->generateCuts(feasCuts, 1 /* = fullscan */, solver_, NULL);
      }
      solver_->applyCuts(feasCuts);

      resolve(solver_);
      double objval = solver_->getObjValue();
      lastHeuristic_ = NULL;
      setBestSolution(CBC_SOLUTION, objval,
        solver_->getColSolution());
      int easy = 2;
      if (!solverCharacteristics_->mipFeasible()) //did we prove that the node could be pruned?
        feasible = false;
      // Reset the bound now
      solverCharacteristics_->setMipBound(-COIN_DBL_MAX);

      solver_->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, &easy);
      feasible &= resolve(oldNode ? oldNode->nodeInfo() : NULL, 11) != 0;
      solver_->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, NULL);
      resolved = true;
      if (problemFeasibility_->feasible(this, 0) < 0) {
        feasible = false; // pretend infeasible
      }
      if (feasible)
        anyAction = -1;
      else
        anyAction = -2;
    }
    /*
          Yep, false positives for sure. And no easy way to distinguish honest
          infeasibility from `found a solution and tightened objective target.'

          if (onOptimalPath)
          assert (anyAction!=-2); // can be useful but gives false positives on strong
          */
    numberPassesLeft--;
    if (numberPassesLeft <= -1) {
      if (!numberLongStrong_ && !numberThreads_)
        messageHandler()->message(CBC_WARNING_STRONG,
          messages())
          << CoinMessageEol;
      numberLongStrong_++;
    }
    if (anyAction == -1) {
      // can do quick optimality check
      int easy = 2;
      solver_->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, &easy);
      feasible = resolve(oldNode ? oldNode->nodeInfo() : NULL, 11) != 0;
      solver_->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, NULL);
      resolved = true;
      if (problemFeasibility_->feasible(this, 0) < 0) {
        feasible = false; // pretend infeasible
      }
      if (feasible) {
        // Set objective value (not so obvious if NLP etc)
        setObjectiveValue(newNode, oldNode);
        reducedCostFix();
        if (newNode->objectiveValue() >= getCutoff())
          anyAction = -2;
      } else {
        anyAction = -2;
      }
    }
  }
  //A candidate has been found; restore the subproblem.
  if (saveNumberRows < solver_->getNumRows()) {
    // delete rows - but leave solution
    int n = solver_->getNumRows();
    int *del = new int[n - saveNumberRows];
    for (int i = saveNumberRows; i < n; i++)
      del[i - saveNumberRows] = i;
    solver_->deleteRows(n - saveNumberRows, del);
    delete[] del;
  }
  /*
      End main loop to choose a branching variable.
    */
  if (anyAction >= 0) {
    if (resolved) {
      /*
              Used to be that when the node was not fathomed (branching object present)
              the solution was not needed. But that's no longer the case --- heuristics
              are applied, and they may want the solution.
            */
      // bool needValidSolution = (newNode->branchingObject() == NULL) ;
      bool needValidSolution = true;
      takeOffCuts(cuts, needValidSolution, NULL);
#ifdef CHECK_CUT_COUNTS
      {
        printf("Number of rows after chooseBranch fix (node)"
               "(active only) %d\n",
          numberRowsAtContinuous_ + numberNewCuts_ + numberOldActiveCuts_);
        const CoinWarmStartBasis *debugws = dynamic_cast< const CoinWarmStartBasis * >(solver_->getWarmStart());
        debugws->print();
        delete debugws;
      }
#endif
    }
    {
      OsiBranchingObject *branchingObject = newNode->modifiableBranchingObject();
      CbcGeneralBranchingObject *generalBranch = dynamic_cast< CbcGeneralBranchingObject * >(branchingObject);
      if (generalBranch && false) {
        int numberProblems = generalBranch->numberSubProblems();
        for (int i = 0; i < numberProblems; i++) {
          double objectiveValue;
          double sumInfeasibilities;
          int numberUnsatisfied;
          generalBranch->state(objectiveValue, sumInfeasibilities,
            numberUnsatisfied, i);
          printf("node %d obj %g sumI %g numI %i rel depth %d\n",
            i, objectiveValue, sumInfeasibilities, numberUnsatisfied,
            generalBranch->subProblem(i)->depth_);
        }
      }
      if (generalBranch) {
        int numberProblems = generalBranch->numberSubProblems();
        newNode->setBranchingObject(NULL);
        CbcNode *newNode2 = NULL;
        assert(numberProblems);
        int nProbMinus1 = numberProblems - 1;
        lockThread();
        for (int i = 0; i < currentNumberCuts_; i++) {
          if (addedCuts_[i])
            addedCuts_[i]->increment(nProbMinus1);
        }
        unlockThread();
        for (int i = 0; i < numberProblems; i++) {
          double objectiveValue;
          double sumInfeasibilities;
          int numberUnsatisfied;
          generalBranch->state(objectiveValue, sumInfeasibilities,
            numberUnsatisfied, i);
          //printf("node %d obj %g sumI %g numI %i rel depth %d\n",
          // i,objectiveValue,sumInfeasibilities,numberUnsatisfied,
          // generalBranch->subProblem(i)->depth_);
          newNode2 = new CbcNode();
          newNode2->setDepth(generalBranch->subProblem(i)->depth_ + currentDepth_);
          generalBranch->subProblem(i)->apply(solver_, 8); // basis
          newNode2->setNumberUnsatisfied(numberUnsatisfied);
          newNode2->setSumInfeasibilities(sumInfeasibilities);
          newNode2->setGuessedObjectiveValue(objectiveValue);
          newNode2->setObjectiveValue(objectiveValue);
          CbcOneGeneralBranchingObject *object = new CbcOneGeneralBranchingObject(this, generalBranch, i);
          newNode2->setBranchingObject(object);
          assert(lastws->fullBasis());
          newNode2->createInfo(this, oldNode, lastws,
            lowerBefore, upperBefore,
            numberOldActiveCuts_, numberNewCuts_);
          newNode2->nodeInfo()->setNumberBranchesLeft(1);
          //newNode2->nodeInfo()->unsetParentBasedData();
          if (i < nProbMinus1) {
            //OsiBranchingObject * object = oldNode->modifiableBranchingObject();
            CbcNodeInfo *nodeInfo = oldNode->nodeInfo();
            //object->incrementNumberBranchesLeft();
            nodeInfo->incrementNumberPointingToThis();
            newNode2->nodeInfo()->setNodeNumber(numberNodes2_);
            //newNode2->nodeInfo()->setNumberBranchesLeft(1);
            newNode2->initializeInfo();
            numberNodes2_++;
            tree_->push(newNode2);
          }
        }
        delete newNode;
        newNode = newNode2;
      } else {
        if (lastws) {
          if (parallelMode() < -1) {
            lastws->fixFullBasis();
          } else {
            if ((specialOptions_ & 8192) == 0)
              assert(lastws->fullBasis());
            else
              lastws->fixFullBasis();
          }
        }
        newNode->createInfo(this, oldNode, lastws, lowerBefore, upperBefore,
          numberOldActiveCuts_, numberNewCuts_);
      }
    }
    if (newNode->numberUnsatisfied()) {
      maximumDepthActual_ = CoinMax(maximumDepthActual_, newNode->depth());
      // Number of branches is in oldNode!
      newNode->initializeInfo();
      if (cuts.sizeRowCuts()) {
        int initialNumber = ((threadMode_ & 1) == 0) ? 0 : 1000000000;
        lockThread();
        newNode->nodeInfo()->addCuts(cuts, newNode->numberBranches(),
          //whichGenerator_,
          initialNumber);
        unlockThread();
      }
    }
  } else {
    anyAction = -2;
    // Reset bound anyway (no harm if not odd)
    solverCharacteristics_->setMipBound(-COIN_DBL_MAX);
  }
  // May have slipped through i.e. anyAction == 0 and objective above cutoff
  // I think this will screw up cut reference counts if executed.
  // We executed addCuts just above. (lh)
  if (anyAction >= 0) {
    assert(newNode);
    if (newNode->objectiveValue() >= getCutoff()) {
      anyAction = -2; // say bad after all
      // zap parent nodeInfo
#ifdef COIN_DEVELOP
      printf("zapping3 CbcNodeInfo %x\n", newNode->nodeInfo()->parent());
#endif
      if (newNode->nodeInfo())
        newNode->nodeInfo()->nullParent();
    }
  }
  stateOfSearch_ = 0; // outside chooseBranch
#ifdef JJF_ZERO
  if (onOptimalPath) {
    const OsiRowCutDebugger *debugger = solver_->getRowCutDebugger();
    if (!debugger) {
      printf("NOT On optimal path - choose\n");
      abort();
    } else {
      printf("Still On optimal path - choose\n");
      if (anyAction == -2) {
        printf("anyAction 2!!\n");
        abort();
      }
    }
  }
#endif
  return anyAction;
}

/*
   For advanced applications you may wish to modify the behavior of Cbc
   e.g. if the solver is a NLP solver then you may not have an exact
   optimum solution at each step.  Information could be built into
   OsiSolverInterface but this is an alternative so that that interface
   does not have to be changed.  If something similar is useful to
   enough solvers then it could be migrated.
   You can also pass in by using solver->setAuxiliaryInfo.
   You should do that if solver is odd - if solver is normal simplex
   then use this
*/
void CbcModel::passInSolverCharacteristics(OsiBabSolver *solverCharacteristics)
{
  solverCharacteristics_ = solverCharacteristics;
}
// Generate an OsiBranchingInformation object
OsiBranchingInformation
CbcModel::usefulInformation() const
{
  OsiBranchingInformation usefulInfo(solver_, normalSolver(), false);
  // and modify
  usefulInfo.solution_ = testSolution_;
  usefulInfo.integerTolerance_ = dblParam_[CbcIntegerTolerance];
  usefulInfo.hotstartSolution_ = hotstartSolution_;
  usefulInfo.numberSolutions_ = numberSolutions_;
  usefulInfo.numberBranchingSolutions_ = numberSolutions_ - numberHeuristicSolutions_;
  usefulInfo.depth_ = -1;
  return usefulInfo;
}
void CbcModel::setBestSolution(const double *solution, int numberColumns,
  double objectiveValue, bool checkSolution)
{
  // May be odd discontinuities - so only check if asked
  if (checkSolution) {
    assert(numberColumns == solver_->getNumCols());
    double *saveLower = CoinCopyOfArray(solver_->getColLower(), numberColumns);
    double *saveUpper = CoinCopyOfArray(solver_->getColUpper(), numberColumns);
    // Fix integers
    int numberAway = 0;
    for (int i = 0; i < numberColumns; i++) {
      if (solver_->isInteger(i)) {
        double value = solution[i];
        double intValue = floor(value + 0.5);
        if (fabs(value - intValue) > 1.0e-4)
          numberAway++;
        solver_->setColLower(i, intValue);
        solver_->setColUpper(i, intValue);
      }
    }
    // Save basis
    CoinWarmStart *saveBasis = solver_->getWarmStart();
    // Solve
    solver_->initialSolve();
    char printBuffer[200];
    if (numberAway) {
      sprintf(printBuffer, "Warning %d integer variables were more than 1.0e-4 away from integer", numberAway);
      messageHandler()->message(CBC_GENERAL, messages())
        << printBuffer << CoinMessageEol;
    }
    bool looksGood = solver_->isProvenOptimal();
    if (looksGood) {
      double direction = solver_->getObjSense();
      double objValue = direction * solver_->getObjValue();
      if (objValue > objectiveValue + 1.0e-8 * (1.0 + fabs(objectiveValue))) {
        sprintf(printBuffer, "Given objective value %g, computed %g",
          objectiveValue, objValue);
        messageHandler()->message(CBC_GENERAL, messages())
          << printBuffer << CoinMessageEol;
      }
      // Use this as objective value and solution
      objectiveValue = objValue;
      solution = solver_->getColSolution();
      // Save current basis
      CoinWarmStartBasis *ws = dynamic_cast< CoinWarmStartBasis * >(solver_->getWarmStart());
      assert(ws);
      setBestSolutionBasis(*ws);
      delete ws;
    }
    // Restore basis
    solver_->setWarmStart(saveBasis);
    delete saveBasis;
    // Restore bounds
    solver_->setColLower(saveLower);
    delete[] saveLower;
    solver_->setColUpper(saveUpper);
    delete[] saveUpper;
    // Return if no good
    if (!looksGood) {
      messageHandler()->message(CBC_GENERAL, messages())
        << "Error solution not saved as not feasible" << CoinMessageEol;
      return;
    } else {
      // message
      sprintf(printBuffer, "Solution with objective value %g saved",
        objectiveValue);
      messageHandler()->message(CBC_GENERAL, messages())
        << printBuffer << CoinMessageEol;
    }
  }
  if (bestSolution_)
    saveExtraSolution(bestSolution_, bestObjective_);
  bestObjective_ = objectiveValue;
  // may be able to change cutoff now
  double cutoff = getCutoff();
  double increment = getDblParam(CbcModel::CbcCutoffIncrement);
  if (cutoff > objectiveValue - increment) {
    cutoff = objectiveValue - increment;
    setCutoff(cutoff);
    // change cutoff as constraint if wanted
    if (cutoffRowNumber_ >= 0) {
      if (solver_->getNumRows() > cutoffRowNumber_) {
        double offset;
        solver_->getDblParam(OsiObjOffset, offset);
        solver_->setRowUpper(cutoffRowNumber_, cutoff + offset);
        if (continuousSolver_ && solver_->getNumCols() > continuousSolver_->getNumCols()) {
          solver_->setRowUpper(cutoffRowNumber_, floor(cutoff) + offset);
          solver_->setRowLower(cutoffRowNumber_, floor(cutoff) + offset);
        }
      }
    }
  }
  int n = CoinMax(numberColumns, solver_->getNumCols());
  delete[] bestSolution_;
  bestSolution_ = new double[n];
  memset(bestSolution_, 0, n * sizeof(double));
  memcpy(bestSolution_, solution, numberColumns * sizeof(double));
}
/* Do heuristics at root.
   0 - don't delete
   1 - delete
      2 - just delete - don't even use
  Parameter of 2 means what it says --- the routine will do nothing except
  delete the existing heuristics. A feasibility pump is always deleted,
  independent of the parameter value, as it's only useful at the root.

  The routine is called from branchAndBound to process the root node. But it
  will also be called when we've recursed into branchAndBound via smallBaB.
*/
void CbcModel::doHeuristicsAtRoot(int deleteHeuristicsAfterwards)
{

  int numberColumns = getNumCols();
  double *newSolution = new double[numberColumns];
  int i;
  if (deleteHeuristicsAfterwards != 2) {
    /*
          If mode == 1, we delete and recreate here, then delete at the bottom. The
          create/delete part makes sense, but why delete the existing array? Seems like
          it should be preserved and restored.
        */
    if (deleteHeuristicsAfterwards) {
      delete[] usedInSolution_;
      usedInSolution_ = new int[numberColumns];
      CoinZeroN(usedInSolution_, numberColumns);
    }
    double heuristicValue = getCutoff();
    int found = -1; // no solution found
    CbcEventHandler *eventHandler = getEventHandler();
    if (eventHandler)
      eventHandler->setModel(this);
    /*
          currentPassNumber_ is described as `cut pass number'. Again, seems a bit
          cavalier to just change it.

          Whether this has any effect is determined by individual heuristics. Typically
          there will be a check at the front of the solution() routine that determines
          whether it will run or simply return. Root heuristics are characterised by
          node count of 0. In addition, currentPassNumber_ can be checked to to limit
          execution in terms of passes through cut generation / heuristic execution in
          solveWithCuts.
        */

    currentPassNumber_ = 1; // so root heuristics will run
    /*
          A loop to run the heuristics. incrementUsed will mark entries in
          usedInSolution corresponding to variables that are nonzero in the solution.
          CBC_ROUNDING just identifies a message template, not the heuristic.
        */
    // Modify based on size etc
    adjustHeuristics();
    // See if already within allowable gap
    bool exitNow = false;
    for (i = 0; i < numberHeuristics_; i++) {
      if (heuristic_[i]->exitNow(bestObjective_))
        exitNow = true;
    }
    if (!exitNow) {
      /** -1 first time otherwise number of solutions last time */
      int lastSolutionCount = -1;
      while (lastSolutionCount) {
        int thisSolutionCount = 0;
#ifdef CBC_THREAD
        if ((threadMode_ & 4) != 0) {
          typedef struct {
            double solutionValue;
            CbcModel *model;
            double *solution;
            int foundSol;
          } argBundle;
          int chunk;
          if (!numberThreads_)
            chunk = numberHeuristics_;
          else
            chunk = numberThreads_;
          for (int iChunk = 0; iChunk < numberHeuristics_; iChunk += chunk) {
            argBundle *parameters = new argBundle[chunk];
            for (int i = 0; i < chunk; i++)
              parameters[i].model = NULL;
            int nThisTime = CoinMin(numberHeuristics_ - iChunk, chunk);
            for (int i = iChunk; i < iChunk + nThisTime; i++) {
              // skip if can't run here
              if (!heuristic_[i]->shouldHeurRun(0))
                continue;
              if (lastSolutionCount > 0 && (heuristic_[i]->switches() & 16) == 0)
                continue; // no point
              parameters[i - iChunk].solutionValue = heuristicValue;
              // Don't want a strategy object
              CbcStrategy *saveStrategy = strategy_;
              strategy_ = NULL;
              CbcModel *newModel = new CbcModel(*this);
              strategy_ = saveStrategy;
              assert(!newModel->continuousSolver_);
              if (continuousSolver_)
                newModel->continuousSolver_ = continuousSolver_->clone();
              else
                newModel->continuousSolver_ = solver_->clone();
              parameters[i - iChunk].model = newModel;
              parameters[i - iChunk].solution = new double[numberColumns];
              ;
              parameters[i - iChunk].foundSol = 0;
              //newModel->gutsOfCopy(*this,-1);
              for (int j = 0; j < numberHeuristics_; j++)
                delete newModel->heuristic_[j];
              //newModel->heuristic_ = new CbcHeuristic * [1];
              newModel->heuristic_[0] = heuristic_[i]->clone();
              newModel->heuristic_[0]->setModel(newModel);
              newModel->heuristic_[0]->resetModel(newModel);
              newModel->numberHeuristics_ = 1;
            }
            void
            parallelHeuristics(int numberThreads,
              int sizeOfData,
              void *argBundle);
            parallelHeuristics(nThisTime,
              static_cast< int >(sizeof(argBundle)),
              parameters);
            double cutoff = heuristicValue;
            for (int i = 0; i < chunk; i++) {
              if (parameters[i].model) {
                if (parameters[i].foundSol > 0 && parameters[i].solutionValue < heuristicValue) {
                  memcpy(newSolution, parameters[i].solution,
                    numberColumns * sizeof(double));
                  lastHeuristic_ = heuristic_[i + iChunk];
                  double value = parameters[i].solutionValue;
                  setBestSolution(CBC_ROUNDING, value, newSolution);
                  // Double check valid
                  if (getCutoff() < cutoff) {
                    cutoff = getCutoff();
                    heuristicValue = value;
                    heuristic_[i + iChunk]->incrementNumberSolutionsFound();
                    incrementUsed(newSolution);
                    // increment number of solutions so other heuristics can test
                    thisSolutionCount++;
                    numberHeuristicSolutions_++;
                    found = i + iChunk;
                  }
                }
                if (heuristic_[i + iChunk]->exitNow(bestObjective_) || (parameters[i].model->heuristic(0)->switches() & (1024 + 2048)) == (1024 + 2048))
                  exitNow = true;
                delete[] parameters[i].solution;
                delete parameters[i].model;
              }
            }
            delete[] parameters;
            if (exitNow)
              break;
          }
        } else {
#endif
          int whereFrom = 0;
          for (i = 0; i < numberHeuristics_; i++) {
            // skip if can't run here
            if (!heuristic_[i]->shouldHeurRun(whereFrom))
              continue;
            if (lastSolutionCount > 0 && (heuristic_[i]->switches() & 16) == 0)
              continue; // no point
            if (maximumSecondsReached()) {
              thisSolutionCount = -1000000;
              break;
            }
            // see if heuristic will do anything
            double saveValue = heuristicValue;
            double before = getCurrentSeconds();
            int ifSol = heuristic_[i]->solution(heuristicValue,
              newSolution);
            if (handler_->logLevel() > 1) {
              char line[100];
              sprintf(line, "Heuristic %s took %g seconds (%s)",
                heuristic_[i]->heuristicName(),
                getCurrentSeconds() - before,
                ifSol ? "good" : "no good");
              handler_->message(CBC_GENERAL, messages_) << line << CoinMessageEol;
            }
            //#define DEBUG_BEST
#ifdef DEBUG_BEST
            FILE *fp = fopen("solution.data", "rb");
            if (!fp && ifSol > 0) {
              int numberColumns = getNumCols();
              fp = fopen("solution.data", "wb");
              printf("Solution data on file solution.data\n");
              size_t numberWritten;
              numberWritten = fwrite(&numberColumns, sizeof(int), 1, fp);
              assert(numberWritten == 1);
              numberWritten = fwrite(&heuristicValue, sizeof(double), 1, fp);
              assert(numberWritten == 1);
              numberWritten = fwrite(newSolution, sizeof(double), numberColumns, fp);
              assert(numberWritten == numberColumns);
              fclose(fp);
            } else if (fp) {
              int numberColumns = getNumCols();
              int numberColumnsX;
              size_t numberRead;
              numberRead = fread(&numberColumnsX, sizeof(int), 1, fp);
              assert(numberRead == 1);
              if (numberColumns == numberColumnsX) {
                numberRead = fread(&heuristicValue, sizeof(double), 1, fp);
                assert(numberRead == 1);
                numberRead = fread(newSolution, sizeof(double), numberColumns, fp);
                assert(numberRead == numberColumns);
                ifSol = 1;
              }
              fclose(fp);
            }
#endif
            if (ifSol > 0) {
              // better solution found
              double currentObjective = bestObjective_;
              CbcHeuristic *saveHeuristic = lastHeuristic_;
              lastHeuristic_ = heuristic_[i];
              setBestSolution(CBC_ROUNDING, heuristicValue, newSolution);
              if (bestObjective_ < currentObjective) {
                thisSolutionCount++;
                heuristic_[i]->incrementNumberSolutionsFound();
                found = i;
                incrementUsed(newSolution);
                // increment number of solutions so other heuristics can test
                //                            numberSolutions_++;
                numberHeuristicSolutions_++;
#ifdef HEURISTIC_INFORM
                printf("HEUR %s where %d C\n",
                  lastHeuristic_->heuristicName(), whereFrom);
#endif
                whereFrom |= 8; // say solution found
                if (heuristic_[i]->exitNow(bestObjective_)
                  || numberSolutions_ >= getMaximumSolutions()) {
                  thisSolutionCount = -1000000;
                  break;
                }
                if (eventHandler) {
                  if (!eventHandler->event(CbcEventHandler::heuristicSolution)) {
                    eventHappened_ = true; // exit
                    thisSolutionCount = -1000000;
                    break;
                  }
                }
                double testGap = CoinMax(dblParam_[CbcAllowableGap],
                  CoinMax(fabs(bestObjective_), fabs(bestPossibleObjective_))
                    * dblParam_[CbcAllowableFractionGap]);
                if (bestObjective_ - bestPossibleObjective_ < testGap && getCutoffIncrement() >= 0.0 && bestPossibleObjective_ < 1.0e30) {
                  if (bestPossibleObjective_ < getCutoff())
                    stoppedOnGap_ = true;
                  //eventHappened_=true; // stop as fast as possible
                  thisSolutionCount = -1000000;
                  break;
                }
                reducedCostFix();
              } else {
                // NOT better solution
#if CBC_USEFUL_PRINTING > 1
                printf("HEUR %s where %d REJECTED i==%d\n",
                  heuristic_[i]->heuristicName(), whereFrom, i);
#endif
                lastHeuristic_ = saveHeuristic;
                heuristicValue = saveValue;
              }
            } else {
              heuristicValue = saveValue;
            }
            if (eventHandler) {
              if (!eventHandler->event(CbcEventHandler::afterHeuristic)) {
                eventHappened_ = true; // exit
                thisSolutionCount = -1000000;
                break;
              }
            }
          }
#ifdef CBC_THREAD
        }
#endif
        if (thisSolutionCount <= 0)
          break;
        lastSolutionCount = thisSolutionCount;
      }
    }
    currentPassNumber_ = 0;
    /*
          Did any of the heuristics turn up a new solution? Record it before we free
        the vector. tree_ will not necessarily be a CbcTreeLocal; the main model gets
        a CbcTree by default. CbcTreeLocal actually implements a k-neighbourhood
        search heuristic. This initialises it with a solution and creates the
        k-neighbourhood cut.
        */
    if (found >= 0) {
      CbcTreeLocal *tree
        = dynamic_cast< CbcTreeLocal * >(tree_);
      if (tree)
        tree->passInSolution(bestSolution_, heuristicValue);
      if (eventHandler) {
        if (!eventHandler->event(CbcEventHandler::solution)) {
          eventHappened_ = true; // exit
        }
      }
    }
  }
  /*
      Cleanup. The feasibility pump heuristic is a root heuristic to look for an
      initial feasible solution. It's had its chance; remove it.

      For modes 1 and 2, all the heuristics are deleted.
    */
  if (!deleteHeuristicsAfterwards) {
    for (i = 0; i < numberHeuristics_; i++) {
      // delete FPump
      CbcHeuristicFPump *pump
        = dynamic_cast< CbcHeuristicFPump * >(heuristic_[i]);
      if (pump && pump->feasibilityPumpOptions() < 1000000
        && (specialOptions_ & 33554432) == 0) {
        delete pump;
        numberHeuristics_--;
        for (int j = i; j < numberHeuristics_; j++)
          heuristic_[j] = heuristic_[j + 1];
      }
    }
  } else {
    // delete all
    for (i = 0; i < numberHeuristics_; i++)
      delete heuristic_[i];
    numberHeuristics_ = 0;
    delete[] heuristic_;
    heuristic_ = NULL;
    delete[] usedInSolution_;
    usedInSolution_ = NULL;
  }
  delete[] newSolution;
}
// Zap integer information in problem (may leave object info)
void CbcModel::zapIntegerInformation(bool leaveObjects)
{
  numberIntegers_ = 0;
  delete[] integerVariable_;
  integerVariable_ = NULL;
  if (!leaveObjects && ownObjects_) {
    int i;
    for (i = 0; i < numberObjects_; i++)
      delete object_[i];
    delete[] object_;
    numberObjects_ = 0;
    object_ = NULL;
  }
}
// Create C++ lines to get to current state
void CbcModel::generateCpp(FILE *fp, int /*options*/)
{
  // Do cut generators
  int i;
  for (i = 0; i < numberCutGenerators_; i++) {
    CglCutGenerator *generator = generator_[i]->generator();
    std::string name = generator->generateCpp(fp);
    int howOften = generator_[i]->howOften();
    int howOftenInSub = generator_[i]->howOftenInSub();
    int whatDepth = generator_[i]->whatDepth();
    int whatDepthInSub = generator_[i]->whatDepthInSub();
    bool normal = generator_[i]->normal();
    bool atSolution = generator_[i]->atSolution();
    bool whenInfeasible = generator_[i]->whenInfeasible();
    bool timing = generator_[i]->timing();
    fprintf(fp, "3  cbcModel->addCutGenerator(&%s,%d,",
      name.c_str(), howOften);
    // change name
    name[0] = static_cast< char >(toupper(name[0]));
    fprintf(fp, "\"%s\",%s,%s,%s,%d,%d,%d);\n",
      name.c_str(), normal ? "true" : "false",
      atSolution ? "true" : "false",
      whenInfeasible ? "true" : "false",
      howOftenInSub, whatDepth, whatDepthInSub);
    fprintf(fp, "3  cbcModel->cutGenerator(%d)->setTiming(%s);\n",
      i, timing ? "true" : "false");
    fprintf(fp, "3  \n");
  }
  for (i = 0; i < numberHeuristics_; i++) {
    CbcHeuristic *heuristic = heuristic_[i];
    heuristic->generateCpp(fp);
    fprintf(fp, "3  \n");
  }
  if (nodeCompare_)
    nodeCompare_->generateCpp(fp);
  tree_->generateCpp(fp);
  CbcModel defaultModel;
  CbcModel *other = &defaultModel;
  int iValue1, iValue2;
  double dValue1, dValue2;
  iValue1 = this->getMaximumNodes();
  iValue2 = other->getMaximumNodes();
  fprintf(fp, "%d  int save_getMaximumNodes = cbcModel->getMaximumNodes();\n", iValue1 == iValue2 ? 2 : 1);
  fprintf(fp, "%d  cbcModel->setMaximumNodes(%d);\n", iValue1 == iValue2 ? 4 : 3, iValue1);
  fprintf(fp, "%d  cbcModel->setMaximumNodes(save_getMaximumNodes);\n", iValue1 == iValue2 ? 7 : 6);
  iValue1 = this->getMaximumSolutions();
  iValue2 = other->getMaximumSolutions();
  fprintf(fp, "%d  int save_getMaximumSolutions = cbcModel->getMaximumSolutions();\n", iValue1 == iValue2 ? 2 : 1);
  fprintf(fp, "%d  cbcModel->setMaximumSolutions(%d);\n", iValue1 == iValue2 ? 4 : 3, iValue1);
  fprintf(fp, "%d  cbcModel->setMaximumSolutions(save_getMaximumSolutions);\n", iValue1 == iValue2 ? 7 : 6);
  iValue1 = this->numberStrong();
  iValue2 = other->numberStrong();
  fprintf(fp, "%d  int save_numberStrong = cbcModel->numberStrong();\n", iValue1 == iValue2 ? 2 : 1);
  fprintf(fp, "%d  cbcModel->setNumberStrong(%d);\n", iValue1 == iValue2 ? 4 : 3, iValue1);
  fprintf(fp, "%d  cbcModel->setNumberStrong(save_numberStrong);\n", iValue1 == iValue2 ? 7 : 6);
  iValue1 = this->numberBeforeTrust();
  iValue2 = other->numberBeforeTrust();
  fprintf(fp, "%d  int save_numberBeforeTrust = cbcModel->numberBeforeTrust();\n", iValue1 == iValue2 ? 2 : 1);
  fprintf(fp, "%d  cbcModel->setNumberBeforeTrust(%d);\n", iValue1 == iValue2 ? 4 : 3, iValue1);
  fprintf(fp, "%d  cbcModel->setNumberBeforeTrust(save_numberBeforeTrust);\n", iValue1 == iValue2 ? 7 : 6);
  iValue1 = this->numberPenalties();
  iValue2 = other->numberPenalties();
  fprintf(fp, "%d  int save_numberPenalties = cbcModel->numberPenalties();\n", iValue1 == iValue2 ? 2 : 1);
  fprintf(fp, "%d  cbcModel->setNumberPenalties(%d);\n", iValue1 == iValue2 ? 4 : 3, iValue1);
  fprintf(fp, "%d  cbcModel->setNumberPenalties(save_numberPenalties);\n", iValue1 == iValue2 ? 7 : 6);
  iValue1 = this->howOftenGlobalScan();
  iValue2 = other->howOftenGlobalScan();
  fprintf(fp, "%d  int save_howOftenGlobalScan = cbcModel->howOftenGlobalScan();\n", iValue1 == iValue2 ? 2 : 1);
  fprintf(fp, "%d  cbcModel->setHowOftenGlobalScan(%d);\n", iValue1 == iValue2 ? 4 : 3, iValue1);
  fprintf(fp, "%d  cbcModel->setHowOftenGlobalScan(save_howOftenGlobalScan);\n", iValue1 == iValue2 ? 7 : 6);
  iValue1 = this->printFrequency();
  iValue2 = other->printFrequency();
  fprintf(fp, "%d  int save_printFrequency = cbcModel->printFrequency();\n", iValue1 == iValue2 ? 2 : 1);
  fprintf(fp, "%d  cbcModel->setPrintFrequency(%d);\n", iValue1 == iValue2 ? 4 : 3, iValue1);
  fprintf(fp, "%d  cbcModel->setPrintFrequency(save_printFrequency);\n", iValue1 == iValue2 ? 7 : 6);
  iValue1 = this->getPrintingMode();
  iValue2 = other->getPrintingMode();
  fprintf(fp, "%d  int save_printingMode = cbcModel->getPrintingMode();\n", iValue1 == iValue2 ? 2 : 1);
  fprintf(fp, "%d  cbcModel->setPrintingMode(%d);\n", iValue1 == iValue2 ? 4 : 3, iValue1);
  fprintf(fp, "%d  cbcModel->setPrintingMode(save_printingMode);\n", iValue1 == iValue2 ? 7 : 6);
  iValue1 = this->searchStrategy();
  iValue2 = other->searchStrategy();
  fprintf(fp, "%d  int save_searchStrategy = cbcModel->searchStrategy();\n", iValue1 == iValue2 ? 2 : 1);
  fprintf(fp, "%d  cbcModel->setSearchStrategy(%d);\n", iValue1 == iValue2 ? 4 : 3, iValue1);
  fprintf(fp, "%d  cbcModel->setSearchStrategy(save_searchStrategy);\n", iValue1 == iValue2 ? 7 : 6);
  iValue1 = this->specialOptions();
  iValue2 = other->specialOptions();
  fprintf(fp, "%d  int save_cbcSpecialOptions = cbcModel->specialOptions();\n", iValue1 == iValue2 ? 2 : 1);
  fprintf(fp, "%d  cbcModel->setSpecialOptions(%d);\n", iValue1 == iValue2 ? 4 : 3, iValue1);
  fprintf(fp, "%d  cbcModel->setSpecialOptions(save_cbcSpecialOptions);\n", iValue1 == iValue2 ? 7 : 6);
  iValue1 = this->messageHandler()->logLevel();
  iValue2 = other->messageHandler()->logLevel();
  fprintf(fp, "%d  int save_cbcMessageLevel = cbcModel->messageHandler()->logLevel();\n", iValue1 == iValue2 ? 2 : 1);
  fprintf(fp, "%d  cbcModel->messageHandler()->setLogLevel(%d);\n", iValue1 == iValue2 ? 4 : 3, iValue1);
  fprintf(fp, "%d  cbcModel->messageHandler()->setLogLevel(save_cbcMessageLevel);\n", iValue1 == iValue2 ? 7 : 6);
  iValue1 = this->getMaximumCutPassesAtRoot();
  iValue2 = other->getMaximumCutPassesAtRoot();
  fprintf(fp, "%d  int save_getMaximumCutPassesAtRoot = cbcModel->getMaximumCutPassesAtRoot();\n", iValue1 == iValue2 ? 2 : 1);
  fprintf(fp, "%d  cbcModel->setMaximumCutPassesAtRoot(%d);\n", iValue1 == iValue2 ? 4 : 3, iValue1);
  fprintf(fp, "%d  cbcModel->setMaximumCutPassesAtRoot(save_getMaximumCutPassesAtRoot);\n", iValue1 == iValue2 ? 7 : 6);
  iValue1 = this->getMaximumCutPasses();
  iValue2 = other->getMaximumCutPasses();
  fprintf(fp, "%d  int save_getMaximumCutPasses = cbcModel->getMaximumCutPasses();\n", iValue1 == iValue2 ? 2 : 1);
  fprintf(fp, "%d  cbcModel->setMaximumCutPasses(%d);\n", iValue1 == iValue2 ? 4 : 3, iValue1);
  fprintf(fp, "%d  cbcModel->setMaximumCutPasses(save_getMaximumCutPasses);\n", iValue1 == iValue2 ? 7 : 6);
  iValue1 = this->getPreferredWay();
  iValue2 = other->getPreferredWay();
  fprintf(fp, "%d  int save_getPreferredWay = cbcModel->getPreferredWay();\n", iValue1 == iValue2 ? 2 : 1);
  fprintf(fp, "%d  cbcModel->setPreferredWay(%d);\n", iValue1 == iValue2 ? 4 : 3, iValue1);
  fprintf(fp, "%d  cbcModel->setPreferredWay(save_getPreferredWay);\n", iValue1 == iValue2 ? 7 : 6);
  dValue1 = this->getMinimumDrop();
  dValue2 = other->getMinimumDrop();
  fprintf(fp, "%d  double save_getMinimumDrop = cbcModel->getMinimumDrop();\n", dValue1 == dValue2 ? 2 : 1);
  fprintf(fp, "%d  cbcModel->setMinimumDrop(%g);\n", dValue1 == dValue2 ? 4 : 3, dValue1);
  fprintf(fp, "%d  cbcModel->setMinimumDrop(save_getMinimumDrop);\n", dValue1 == dValue2 ? 7 : 6);
  dValue1 = this->getIntegerTolerance();
  dValue2 = other->getIntegerTolerance();
  fprintf(fp, "%d  double save_getIntegerTolerance = cbcModel->getIntegerTolerance();\n", dValue1 == dValue2 ? 2 : 1);
  fprintf(fp, "%d  cbcModel->setIntegerTolerance(%g);\n", dValue1 == dValue2 ? 4 : 3, dValue1);
  fprintf(fp, "%d  cbcModel->setIntegerTolerance(save_getIntegerTolerance);\n", dValue1 == dValue2 ? 7 : 6);
  dValue1 = this->getInfeasibilityWeight();
  dValue2 = other->getInfeasibilityWeight();
  fprintf(fp, "%d  double save_getInfeasibilityWeight = cbcModel->getInfeasibilityWeight();\n", dValue1 == dValue2 ? 2 : 1);
  fprintf(fp, "%d  cbcModel->setInfeasibilityWeight(%g);\n", dValue1 == dValue2 ? 4 : 3, dValue1);
  fprintf(fp, "%d  cbcModel->setInfeasibilityWeight(save_getInfeasibilityWeight);\n", dValue1 == dValue2 ? 7 : 6);
  dValue1 = this->getCutoffIncrement();
  dValue2 = other->getCutoffIncrement();
  fprintf(fp, "%d  double save_getCutoffIncrement = cbcModel->getCutoffIncrement();\n", dValue1 == dValue2 ? 2 : 1);
  fprintf(fp, "%d  cbcModel->setCutoffIncrement(%g);\n", dValue1 == dValue2 ? 4 : 3, dValue1);
  fprintf(fp, "%d  cbcModel->setCutoffIncrement(save_getCutoffIncrement);\n", dValue1 == dValue2 ? 7 : 6);
  dValue1 = this->getAllowableGap();
  dValue2 = other->getAllowableGap();
  fprintf(fp, "%d  double save_getAllowableGap = cbcModel->getAllowableGap();\n", dValue1 == dValue2 ? 2 : 1);
  fprintf(fp, "%d  cbcModel->setAllowableGap(%g);\n", dValue1 == dValue2 ? 4 : 3, dValue1);
  fprintf(fp, "%d  cbcModel->setAllowableGap(save_getAllowableGap);\n", dValue1 == dValue2 ? 7 : 6);
  dValue1 = this->getAllowableFractionGap();
  dValue2 = other->getAllowableFractionGap();
  fprintf(fp, "%d  double save_getAllowableFractionGap = cbcModel->getAllowableFractionGap();\n", dValue1 == dValue2 ? 2 : 1);
  fprintf(fp, "%d  cbcModel->setAllowableFractionGap(%g);\n", dValue1 == dValue2 ? 4 : 3, dValue1);
  fprintf(fp, "%d  cbcModel->setAllowableFractionGap(save_getAllowableFractionGap);\n", dValue1 == dValue2 ? 7 : 6);
  dValue1 = this->getMaximumSeconds();
  dValue2 = other->getMaximumSeconds();
  fprintf(fp, "%d  double save_cbcMaximumSeconds = cbcModel->getMaximumSeconds();\n", dValue1 == dValue2 ? 2 : 1);
  fprintf(fp, "%d  cbcModel->setMaximumSeconds(%g);\n", dValue1 == dValue2 ? 4 : 3, dValue1);
  fprintf(fp, "%d  cbcModel->setMaximumSeconds(save_cbcMaximumSeconds);\n", dValue1 == dValue2 ? 7 : 6);
}
// So we can use osiObject or CbcObject during transition
void getIntegerInformation(const OsiObject *object, double &originalLower,
  double &originalUpper)
{
  const CbcSimpleInteger *integerObject = dynamic_cast< const CbcSimpleInteger * >(object);
  if (integerObject) {
    // get original bounds
    originalLower = integerObject->originalLowerBound();
    originalUpper = integerObject->originalUpperBound();
  } else {
    const OsiSimpleInteger *integerObject = dynamic_cast< const OsiSimpleInteger * >(object);
    assert(integerObject);
    // get original bounds
    originalLower = integerObject->originalLowerBound();
    originalUpper = integerObject->originalUpperBound();
  }
}
// Set original columns as created by preprocessing
void CbcModel::setOriginalColumns(const int *originalColumns, int numberGood)
{
  int numberColumns = getNumCols();
  delete[] originalColumns_;
  originalColumns_ = new int[numberColumns];
  int numberCopy = CoinMin(numberColumns, numberGood);
  memcpy(originalColumns_, originalColumns, numberCopy * sizeof(int));
  for (int i = numberCopy; i < numberColumns; i++)
    originalColumns_[i] = -1;
}
// Set the cut modifier method
void CbcModel::setCutModifier(CbcCutModifier *modifier)
{
  delete cutModifier_;
  cutModifier_ = modifier->clone();
}
/* Set the cut modifier method
 */
void CbcModel::setCutModifier(CbcCutModifier &modifier)
{
  delete cutModifier_;
  cutModifier_ = modifier.clone();
}
/* Do one node - broken out for clarity?
   also for parallel (when baseModel!=this)
   Returns 1 if solution found
   node NULL on return if no branches left
   newNode NULL if no new node created
*/
int CbcModel::doOneNode(CbcModel *baseModel, CbcNode *&node, CbcNode *&newNode)
{
  int foundSolution = 0;
  int saveNumberCutGenerators = numberCutGenerators_;
  if ((moreSpecialOptions_ & 33554432) != 0 && (specialOptions_ & 2048) == 0) {
    if (node && (node->depth() == -2 || node->depth() == 4))
      numberCutGenerators_ = 0; // so can dive and branch
  }
  currentNode_ = node; // so can be accessed elsewhere
  double bestObjective = bestObjective_;
  numberUpdateItems_ = 0;
  // Say not on optimal path
  bool onOptimalPath = false;
#ifdef CHECK_NODE
  printf("Node %x popped from tree - %d left, %d count\n", node,
    node->nodeInfo()->numberBranchesLeft(),
    node->nodeInfo()->numberPointingToThis());
  printf("\tdepth = %d, z =  %g, unsat = %d\n", //var = %d.\n",
    node->depth(), node->objectiveValue(),
    node->numberUnsatisfied());
  //node->columnNumber()) ;
#endif

  /*
      Rebuild the subproblem for this node:	 Call addCuts() to adjust the model
      to recreate the subproblem for this node (set proper variable bounds, add
      cuts, create a basis).  This may result in the problem being fathomed by
      bound or infeasibility. Returns 1 if node is fathomed.
      Execute the current arm of the branch: If the problem survives, save the
      resulting variable bounds and call branch() to modify variable bounds
      according to the current arm of the branching object. If we're processing
      the final arm of the branching object, flag the node for removal from the
      live set.
    */
  /*
      Used to generate bound edits for CbcPartialNodeInfo.
    */
  int numberColumns = getNumCols();
  double *lowerBefore = new double[numberColumns];
  double *upperBefore = new double[numberColumns];
  if (parallelMode() >= 0)
    newNode = NULL;
  else
    newNode = new CbcNode();
  bool feasible = true;
  CoinWarmStartBasis *lastws = new CoinWarmStartBasis();
  lockThread();
  // point to genuine ones
  //int save1 = maximumNumberCuts_;
  //maximumNumberCuts_ = baseModel->maximumNumberCuts_;
  //addedCuts_ = baseModel->addedCuts_;
  if (parallelMode() >= 0) {
    maximumDepth_ = baseModel->maximumDepth_;
    walkback_ = baseModel->walkback_;
    lastNodeInfo_ = baseModel->lastNodeInfo_;
    lastNumberCuts_ = baseModel->lastNumberCuts_;
    lastCut_ = baseModel->lastCut_;
    lastNumberCuts2_ = baseModel->lastNumberCuts2_;
  }
  int save2 = maximumDepth_;
  int retCode = addCuts(node, lastws);
#ifdef SWITCH_VARIABLES
  fixAssociated(solver_, 0);
#endif
  //if (save1<maximumNumberCuts_) {
  // increased
  //baseModel->maximumNumberCuts_ = maximumNumberCuts_;
  //baseModel->addedCuts_ = addedCuts_;
  //}
  if (parallelMode() >= 0 && save2 < maximumDepth_) {
    // increased
    baseModel->maximumDepth_ = maximumDepth_;
    baseModel->walkback_ = walkback_;
    baseModel->lastNodeInfo_ = lastNodeInfo_;
    baseModel->lastNumberCuts_ = lastNumberCuts_;
    baseModel->lastCut_ = lastCut_;
    baseModel->lastNumberCuts2_ = lastNumberCuts2_;
  }
  int branchesLeft = 0;
  if (!retCode) {
    unlockThread();
    int i;
    const double *lower = getColLower();
    const double *upper = getColUpper();
    for (i = 0; i < numberColumns; i++) {
      lowerBefore[i] = lower[i];
      upperBefore[i] = upper[i];
    }
    if ((solverCharacteristics_->extraCharacteristics() & 2) != 0) {
      solverCharacteristics_->setBeforeLower(lowerBefore);
      solverCharacteristics_->setBeforeUpper(upperBefore);
    }
    lockThread();
    assert(node->objectiveValue() < 1.0e200);
    if (messageHandler()->logLevel() > 2)
      node->modifiableBranchingObject()->print();
    if (branchingMethod_ && branchingMethod_->chooseMethod()) {
      branchesLeft = node->branch(solver_); // new way
    } else {
      // old way so need to cheat
      OsiBranchingObject *branch2 = node->modifiableBranchingObject();
#ifndef NDEBUG
      CbcBranchingObject *branch = dynamic_cast< CbcBranchingObject * >(branch2);
      assert(branch);
#else
      CbcBranchingObject *branch = static_cast< CbcBranchingObject * >(branch2);
#endif
#if 1
      branch->setModel(this);
      branchesLeft = node->branch(NULL); // old way
#else
      branchesLeft = node->branch(solver_);
#endif
      if (parallelMode() >= 0)
        branch->setModel(baseModel);
    }
    assert(branchesLeft == node->nodeInfo()->numberBranchesLeft());
    if (parallelMode() > 0) {
      assert(masterThread_);
      assert(node->nodeInfo());
      node->nodeInfo()->increment();
      unlockThread();
    }
    if ((specialOptions_ & 1) != 0) {
      /*
            This doesn't work as intended --- getRowCutDebugger will return null
            unless the current feasible solution region includes the optimal solution
            that RowCutDebugger knows. There's no way to tell inactive from off the
            optimal path.
            */
      const OsiRowCutDebugger *debugger = solver_->getRowCutDebugger();
      if (debugger) {
        onOptimalPath = true;
        printf("On optimal path\n");
      }
    }

    /*
          Reoptimize, possibly generating cuts and/or using heuristics to find
          solutions.  Cut reference counts are unaffected unless we lose feasibility,
          in which case solveWithCuts() will make the adjustment.
        */
    phase_ = 2;
    OsiCuts cuts;
    int saveNumber = numberIterations_;
    if (solverCharacteristics_->solutionAddsCuts()) {
      int returnCode = resolve(node ? node->nodeInfo() : NULL, 1);
      feasible = returnCode != 0;
      if (feasible) {
        int iObject;
        int numberUnsatisfied = 0;
        memcpy(currentSolution_, solver_->getColSolution(),
          numberColumns * sizeof(double));
        // point to useful information
        OsiBranchingInformation usefulInfo = usefulInformation();

        for (iObject = 0; iObject < numberObjects_; iObject++) {
          double infeasibility = object_[iObject]->checkInfeasibility(&usefulInfo);
          if (infeasibility)
            numberUnsatisfied++;
        }
        if (returnCode > 0) {
          if (numberUnsatisfied) {
            feasible = solveWithCuts(cuts, maximumCutPasses_, node);
          } else {
            // may generate cuts and turn the solution
            //to an infeasible one
            feasible = solveWithCuts(cuts, 1,
              node);
          }
        }
        // check extra info on feasibility
        if (!solverCharacteristics_->mipFeasible()) {
          feasible = false;
          solverCharacteristics_->setMipBound(-COIN_DBL_MAX);
        }
      }
    } else {
      // normal
      if (false) {
        const double *lower = solver_->getColLower();
        const double *upper = solver_->getColUpper();
        printf("STATE before solve\n");
        for (int i = 0; i < 100; i++)
          if (lower[i] || !upper[i])
            printf("%d fixed to %g\n", i, lower[i]);
      }
#ifdef COIN_HAS_CLP
      OsiClpSolverInterface *clpSolver
        = dynamic_cast< OsiClpSolverInterface * >(solver_);
      if ((clpSolver || (specialOptions_ & 16384) != 0) && fastNodeDepth_ < -1
        && (specialOptions_ & 2048) == 0) {
#define FATHOM_BIAS -2
        if (numberNodes_ == 1) {
          int numberNodesBeforeFathom = 500;
          if (fastNodeDepth_ < -1000001) {
            numberNodesBeforeFathom = (-fastNodeDepth_) / 1000000;
            numberNodesBeforeFathom = 250 * numberNodesBeforeFathom;
          }
#ifdef COIN_DEVELOP
          int fastNodeDepth1 = -fastNodeDepth_ % 1000000;
          printf("initial depth %d after %d nodes\n",
            FATHOM_BIAS + fastNodeDepth1, numberNodesBeforeFathom);
#endif
        }
        //#endif
        ClpNodeStuff stuff;
        ClpNodeStuff *info = &stuff;
        /*
                  Used to generate bound edits for CbcPartialNodeInfo.
                */
        //double * lowerBefore = NULL;
        //double * upperBefore = NULL;
        int fastNodeDepth1 = -fastNodeDepth_ % 1000000;
        int numberNodesBeforeFathom = 500;
        if (fastNodeDepth_ < -1000001) {
          numberNodesBeforeFathom = (-fastNodeDepth_) / 1000000;
          numberNodesBeforeFathom = 100 * numberNodesBeforeFathom;
        }
        int go_fathom = FATHOM_BIAS + fastNodeDepth1;
        if ((specialOptions_ & 16384) != 0)
          numberNodesBeforeFathom = 0;
        if (node->depth() >= go_fathom && (specialOptions_ & 2048) == 0
          //if (node->depth()>=FATHOM_BIAS-fastNodeDepth_&&!parentModel_
          && numberNodes_ >= numberNodesBeforeFathom && !hotstartSolution_) {
#ifndef COIN_HAS_CPX
          specialOptions_ &= ~16384;
#endif
          if ((specialOptions_ & 16384) == 0) {
            info->integerTolerance_ = getIntegerTolerance();
            info->integerIncrement_ = getCutoffIncrement();
            info->numberBeforeTrust_ = numberBeforeTrust_;
            info->stateOfSearch_ = 1;
            if (numberSolutions_ > 0) {
              info->stateOfSearch_ = 3;
            }
            if (numberNodes_ > 2 * numberObjects_ + 1000) {
              info->stateOfSearch_ = 4;
            }
            // Compute "small" change in branch
            int nBranches = intParam_[CbcNumberBranches];
            if (nBranches) {
              double average = dblParam_[CbcSumChange] / static_cast< double >(nBranches);
              info->smallChange_ = CoinMax(average * 1.0e-5, dblParam_[CbcSmallestChange]);
              info->smallChange_ = CoinMax(info->smallChange_, 1.0e-8);
            } else {
              info->smallChange_ = 1.0e-8;
            }
            double *down = new double[numberIntegers_];
            double *up = new double[numberIntegers_];
            int *priority = new int[numberIntegers_];
            int *numberDown = new int[numberIntegers_];
            int *numberUp = new int[numberIntegers_];
            int *numberDownInfeasible = new int[numberIntegers_];
            int *numberUpInfeasible = new int[numberIntegers_];
            fillPseudoCosts(down, up, priority, numberDown, numberUp,
              numberDownInfeasible, numberUpInfeasible);
            // See if all priorities same
            bool allSame = true;
            int kPriority = priority[0];
            for (int i = 1; i < numberIntegers_; i++) {
              if (kPriority != priority[i]) {
                allSame = false;
                break;
              }
            }
            ClpSimplex *simplex = clpSolver->getModelPtr();
            double *saveLower = CoinCopyOfArray(solver_->getColLower(), numberColumns);
            double *saveUpper = CoinCopyOfArray(solver_->getColUpper(), numberColumns);
            if (allSame && false) {
              // change priorities on general
              const double *lower = simplex->columnLower();
              const double *upper = simplex->columnUpper();
              for (int i = 0; i < numberIntegers_; i++) {
                int iColumn = integerVariable_[i];
                if (upper[iColumn] > lower[iColumn] + 1.1)
                  priority[i] = kPriority + 1;
              }
            }
            info->fillPseudoCosts(down, up, priority, numberDown, numberUp,
              numberDownInfeasible,
              numberUpInfeasible, numberIntegers_);
            info->presolveType_ = 1;
            // for reduced costs and duals
            info->solverOptions_ |= 7;
            delete[] down;
            delete[] up;
            delete[] numberDown;
            delete[] priority;
            delete[] numberUp;
            delete[] numberDownInfeasible;
            delete[] numberUpInfeasible;
            bool takeHint;
            OsiHintStrength strength;
            solver_->getHintParam(OsiDoReducePrint, takeHint, strength);
            //printf("mod cutoff %g solver %g offset %g\n",
            //   getCutoff(),simplex->dualObjectiveLimit(),simplex->objectiveOffset());
            int saveLevel = simplex->logLevel();
            if (strength != OsiHintIgnore && takeHint && saveLevel == 1)
              simplex->setLogLevel(0);
            clpSolver->setBasis();
            int perturbation = simplex->perturbation();
            if ((specialOptions_ & 131072) != 0) {
              //assert (perturbation == 100);
              simplex->setPerturbation(50);
            }
            int saveMoreOptions = simplex->moreSpecialOptions();
            int flags = (moreSpecialOptions_ >> 18) & 3;
            simplex->setMoreSpecialOptions(saveMoreOptions | flags << 11);
#ifndef NO_FATHOM_PRINT
            info->startingDepth_ = node->depth();
            info->nodeCalled_ = numberNodes_;
            info->handler_ = handler_;
#endif
            feasible = simplex->fathom(info) != 0;
            simplex->setMoreSpecialOptions(saveMoreOptions);
            simplex->setPerturbation(perturbation);
            incrementExtra(info->numberNodesExplored_,
              info->numberIterations_);
            if (feasible && false) { // can mess up cuts
              double objValue = simplex->objectiveValue();
              feasible = solveWithCuts(cuts, 1, node);
              if (!feasible) {
                if ((specialOptions_ & 1) != 0)
                  printf("small was feasible %g now infeasible! - depths %d %d\n",
                    objValue, node->depth(), fastNodeDepth_);
                // switch off
                info->nNodes_ = -99;
                solver_->setColLower(saveLower);
                solver_->setColUpper(saveUpper);
              }
            }
            char general[200];
            int fathomStatus = info->nNodes_;
            if (feasible)
              fathomStatus = 1;
            sprintf(general, "fathom took %d nodes, %d iterations - status %d",
              info->numberNodesExplored_,
              info->numberIterations_, fathomStatus);
            messageHandler()->message(CBC_FPUMP2,
              messages())
              << general << CoinMessageEol;
            if (info->numberNodesExplored_ > 10000 /* && !feasible */
              && (moreSpecialOptions_ & 524288) == 0 && info->nNodes_ >= 0) {
              fastNodeDepth_--;
#ifndef NO_FATHOM_PRINT
              if ((moreSpecialOptions_ & 262144) != 0)
                handler_->message(CBC_FATHOM_CHANGE, messages_) << FATHOM_BIAS - fastNodeDepth_ << CoinMessageEol;
#endif
#if CBC_USEFUL_PRINTING > 0
              printf(">10000 - depth now %d so at depth >= %d\n",
                fastNodeDepth_, FATHOM_BIAS - fastNodeDepth_);
#endif
            }
            if (info->nNodes_ < 0) {
              // we gave up
              //abort();
              fastNodeDepth_ -= (info->nNodes_ == -10) ? 5 : 2;
              if (info->nNodes_ == -99)
                fastNodeDepth_ = -1; // switch off
#ifndef NO_FATHOM_PRINT
              if ((moreSpecialOptions_ & 262144) != 0)
                handler_->message(CBC_FATHOM_CHANGE, messages_) << FATHOM_BIAS - fastNodeDepth_ << CoinMessageEol;
#endif
#if CBC_USEFUL_PRINTING > 0
              printf("gave up fastNodeDepth now %d - so at depth >= %d\n",
                fastNodeDepth_, FATHOM_BIAS - fastNodeDepth_);
#endif
              if (feasible) {
                // Save bounds round bestSolution
                //double * saveLower = CoinCopyOfArray(solver_->getColLower(),
                //			     numberColumns);
                //double * saveUpper = CoinCopyOfArray(solver_->getColUpper(),
                //			     numberColumns);
                clpSolver->setWarmStart(NULL);
                // try and do solution
                double value = simplex->objectiveValue();
                double *newSolution = CoinCopyOfArray(simplex->primalColumnSolution(),
                  numberColumns);
                double saveBest = bestObjective_;
                setBestSolution(CBC_STRONGSOL, value, newSolution);
                delete[] newSolution;
                if (bestObjective_ == saveBest) {
                  if ((specialOptions_ & 1) != 0)
                    printf("small was feasible now just infeasible! - depths %d %d\n",
                      node->depth(), fastNodeDepth_);
                  fastNodeDepth_ = -1; // switch off
                  solver_->setColLower(saveLower);
                  solver_->setColUpper(saveUpper);
                }
              }
              // say feasible so will redo node
              feasible = true;
            } else {
              if (feasible) {
                clpSolver->setWarmStart(NULL);
                // try and do solution
                double value = simplex->objectiveValue();
                double *newSolution = CoinCopyOfArray(simplex->primalColumnSolution(),
                  numberColumns);
                double saveBest = bestObjective_;
                setBestSolution(CBC_STRONGSOL, value, newSolution);
                // in case of inaccuracy
                simplex->setObjectiveValue(CoinMax(bestObjective_,
                  simplex->objectiveValue()));
                delete[] newSolution;
                if (bestObjective_ == saveBest) {
                  if ((specialOptions_ & 1) != 0)
                    printf("small was feasible now just infeasible! - depths %d %d\n",
                      node->depth(), fastNodeDepth_);
                  fastNodeDepth_ = -1; // switch off
                  solver_->setColLower(saveLower);
                  solver_->setColUpper(saveUpper);
                }
              }
              // update pseudo costs
              double smallest = 1.0e50;
              double largest = -1.0;
              for (int i = 0; i < numberIntegers_; i++) {
                CbcSimpleIntegerDynamicPseudoCost *obj = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object_[i]);
                if (!obj)
                  continue;
                assert(obj->columnNumber() == integerVariable_[i]);
                if (info->numberUp_[i] > 0) {
                  if (info->downPseudo_[i] > largest)
                    largest = info->downPseudo_[i];
                  if (info->downPseudo_[i] < smallest)
                    smallest = info->downPseudo_[i];
                  if (info->upPseudo_[i] > largest)
                    largest = info->upPseudo_[i];
                  if (info->upPseudo_[i] < smallest)
                    smallest = info->upPseudo_[i];
                  obj->updateAfterMini(info->numberDown_[i],
                    info->numberDownInfeasible_[i],
                    info->downPseudo_[i],
                    info->numberUp_[i],
                    info->numberUpInfeasible_[i],
                    info->upPseudo_[i]);
                }
              }
              //printf("range of costs %g to %g\n",smallest,largest);
            }
            delete[] saveLower;
            delete[] saveUpper;
            simplex->setLogLevel(saveLevel);
#ifdef COIN_HAS_CPX
          } else {
            // try cplex
            OsiCpxSolverInterface cpxSolver;
            double direction = clpSolver->getObjSense();
            cpxSolver.setObjSense(direction);
            // load up cplex
            const CoinPackedMatrix *matrix = clpSolver->getMatrixByCol();
            const double *rowLower = clpSolver->getRowLower();
            const double *rowUpper = clpSolver->getRowUpper();
            const double *columnLower = clpSolver->getColLower();
            const double *columnUpper = clpSolver->getColUpper();
            const double *objective = clpSolver->getObjCoefficients();
            cpxSolver.loadProblem(*matrix, columnLower, columnUpper,
              objective, rowLower, rowUpper);
            double *setSol = new double[numberIntegers_];
            int *setVar = new int[numberIntegers_];
            // cplex doesn't know about objective offset
            double offset = clpSolver->getModelPtr()->objectiveOffset();
            for (int i = 0; i < numberIntegers_; i++) {
              int iColumn = integerVariable_[i];
              cpxSolver.setInteger(iColumn);
              if (bestSolution_) {
                setSol[i] = bestSolution_[iColumn];
                setVar[i] = iColumn;
              }
            }
            CPXENVptr env = cpxSolver.getEnvironmentPtr();
            CPXLPptr lpPtr = cpxSolver.getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL);
            cpxSolver.switchToMIP();
            if (bestSolution_) {
#if 0
                            CPXcopymipstart(env, lpPtr, numberIntegers_, setVar, setSol);
#else
              int zero = 0;
              CPXaddmipstarts(env, lpPtr, 1, numberIntegers_, &zero, setVar, setSol, NULL, NULL);
#endif
            }
            if (getCutoff() < 1.0e50) {
              double useCutoff = getCutoff() + offset;
              if (bestObjective_ < 1.0e50)
                useCutoff = bestObjective_ + offset + 1.0e-7;
              cpxSolver.setDblParam(OsiDualObjectiveLimit, useCutoff * direction);
              if (direction > 0.0)
                CPXsetdblparam(env, CPX_PARAM_CUTUP, useCutoff); // min
              else
                CPXsetdblparam(env, CPX_PARAM_CUTLO, useCutoff); // max
            }
            CPXsetdblparam(env, CPX_PARAM_EPGAP, dblParam_[CbcAllowableFractionGap]);
            delete[] setSol;
            delete[] setVar;
            if (offset) {
              char printBuffer[200];
              sprintf(printBuffer, "Add %g to all Cplex messages for true objective",
                -offset);
              messageHandler()->message(CBC_GENERAL, messages())
                << printBuffer << CoinMessageEol;
              cpxSolver.setDblParam(OsiObjOffset, offset);
            }
            cpxSolver.branchAndBound();
            numberExtraNodes_ += CPXgetnodecnt(env, lpPtr);
            numberExtraIterations_ += CPXgetmipitcnt(env, lpPtr);
            double value = cpxSolver.getObjValue() * direction;
            if (cpxSolver.isProvenOptimal() && value <= getCutoff()) {
              feasible = true;
              clpSolver->setWarmStart(NULL);
              // try and do solution
              double *newSolution = CoinCopyOfArray(cpxSolver.getColSolution(),
                getNumCols());
              setBestSolution(CBC_STRONGSOL, value, newSolution);
              delete[] newSolution;
            }
#endif
          }
        }
      }
      if (feasible) {
        //int numberPasses = doCutsNow(1) ? maximumCutPasses_ : 0;
        int numberPasses = /*doCutsNow(1) ?*/ maximumCutPasses_ /*: 0*/;
        feasible = solveWithCuts(cuts, numberPasses, node);
      }
#else
      feasible = solveWithCuts(cuts, maximumCutPasses_, node);
#endif
    }
    if ((specialOptions_ & 1) != 0 && onOptimalPath) {
      if (solver_->getRowCutDebuggerAlways()->optimalValue() < getCutoff()) {
        if (!solver_->getRowCutDebugger() || !feasible) {
          // dj fix did something???
          solver_->writeMpsNative("infeas2.mps", NULL, NULL, 2);
          solver_->getRowCutDebuggerAlways()->printOptimalSolution(*solver_);
#ifndef NDEBUG
          const OsiRowCutDebugger *debugger = solver_->getRowCutDebugger();
#endif
          assert(debugger);
          int numberRows0 = continuousSolver_->getNumRows();
          int numberRows = solver_->getNumRows();
          const CoinPackedMatrix *rowCopy = solver_->getMatrixByRow();
          const int *rowLength = rowCopy->getVectorLengths();
          const double *elements = rowCopy->getElements();
          const int *column = rowCopy->getIndices();
          const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
          const double *rowLower = solver_->getRowLower();
          const double *rowUpper = solver_->getRowUpper();
          for (int iRow = numberRows0; iRow < numberRows; iRow++) {
            OsiRowCut rc;
            rc.setLb(rowLower[iRow]);
            rc.setUb(rowUpper[iRow]);
            CoinBigIndex start = rowStart[iRow];
            rc.setRow(rowLength[iRow], column + start, elements + start, false);
            CoinAssert(!debugger->invalidCut(rc));
          }
          assert(feasible);
        }
      }
    }
    if (statistics_) {
      assert(numberNodes2_);
      assert(statistics_[numberNodes2_ - 1]);
      assert(statistics_[numberNodes2_ - 1]->node() == numberNodes2_ - 1);
      statistics_[numberNodes2_ - 1]->endOfBranch(numberIterations_ - saveNumber,
        feasible ? solver_->getObjValue()
                 : COIN_DBL_MAX);
    }
    /*
          Are we still feasible? If so, create a node and do the work to attach a
          branching object, reoptimising as needed if chooseBranch() identifies
          monotone objects.

          Finally, attach a partial nodeInfo object and store away any cuts that we
          created back in solveWithCuts. addCuts() will initialise the reference
          counts for these new cuts.

          This next test can be problematic if we've discovered an
          alternate equivalent answer and subsequently fathom the solution
          known to the row cut debugger due to bounds.
        */
    if (onOptimalPath) {
      bool objLim = solver_->isDualObjectiveLimitReached();
      if (!feasible && !objLim) {
        if (solver_->getRowCutDebuggerAlways()->optimalValue() < getCutoff()) {
          printf("infeas2\n");
          solver_->getRowCutDebuggerAlways()->printOptimalSolution(*solver_);
          solver_->writeMpsNative("infeas.mps", NULL, NULL, 2);
          CoinWarmStartBasis *slack = dynamic_cast< CoinWarmStartBasis * >(solver_->getEmptyWarmStart());
          solver_->setWarmStart(slack);
          delete slack;
          solver_->setHintParam(OsiDoReducePrint, false, OsiHintDo, 0);
          solver_->initialSolve();
          assert(!solver_->isProvenOptimal());
          assert(feasible || objLim);
        }
      }
    }
    bool checkingNode = false;
    if (feasible) {
#ifdef FUNNY_BRANCHING2
      // Far too clever
      if ((numberThreads_ == -10 || true) && node->numberBranches() == 2) {
        // see if any parent branches redundant
        // Look at state of "node"
        CbcNodeInfo *nodeInfo = node->nodeInfo();
        if (nodeInfo) {
          // See if any branched variables off bounds
          const double *dj = solver_->getReducedCost();
          const double *lower = solver_->getColLower();
          const double *upper = solver_->getColUpper();
          const double *solution = solver_->getColSolution();
          int numberColumns = solver_->getNumCols();
          double *currentLower = CoinCopyOfArray(lower, numberColumns);
          double *currentUpper = CoinCopyOfArray(upper, numberColumns);
          char *touched = new char[numberColumns];
          memset(touched, 0, numberColumns);
          double direction = solver_->getObjSense();
          bool canDelete = nodeInfo->numberBranchesLeft() > 0;
          //int numberBounds = nodeInfo->numberChangedBounds();
          //const int * which = nodeInfo->variables();
          //const double * bounds = nodeInfo->newBounds();
          const OsiBranchingObject *obj = node->branchingObject();
          const CbcIntegerBranchingObject *objectI = dynamic_cast< const CbcIntegerBranchingObject * >(obj);
          if (objectI) {
            const CbcSimpleInteger *object1 = dynamic_cast< const CbcSimpleInteger * >(objectI->object());
            int iColumn1 = -1;
            int way1 = 0;
            const double *bounds1 = NULL;
            bool zeroOne1 = false;
            if (object1) {
              iColumn1 = object1->columnNumber();
              double originalLower1 = object1->originalLowerBound();
              double originalUpper1 = object1->originalUpperBound();
              // Unset all bounds from parents
              CbcPartialNodeInfo *partial = dynamic_cast< CbcPartialNodeInfo * >(nodeInfo);
              touched[iColumn1] = 1;
              if (partial) {
                /* maybe don't do if obj hasn't changed
                                as then you might get loop
                                at present just 0-1
                                as need to know original bound
                                */
                int n = partial->numberChangedBounds();
                const int *which = partial->variables();
                const double *values = partial->newBounds();
                for (int i = 0; i < n; i++) {
                  int variable = which[i];
                  int k = variable & 0x3fffffff;
                  assert(k != iColumn1);
                  if (!touched[k]) {
                    if ((variable & 0x80000000) == 0) {
                      // lower bound changing
                      assert(currentLower[k] == 1.0);
                      currentLower[k] = 0.0;
                    } else {
                      // upper bound changing
                      assert(currentUpper[k] == 0.0);
                      currentUpper[k] = 1.0;
                    }
                  }
                }
              }
              zeroOne1 = originalLower1 == 0.0 && originalUpper1 == 1.0;
              way1 = objectI->way();
              assert(way1 == -1 || way1 == 1);
              int kWay = way1;
              //way1 = -way1; // what last branch did
              // work out using bounds
              if (objectI->downBounds()[1] >= upper[iColumn1] && objectI->downBounds()[0] <= lower[iColumn1])
                way1 = -1;
              else
                way1 = 1;
              assert(kWay == -way1);
              if (way1 < 0) {
                // must have been down branch
                bounds1 = objectI->downBounds();
              } else {
                // must have been up branch
                bounds1 = objectI->upBounds();
              }
              // double check bounds
              assert(bounds1[0] <= lower[iColumn1] && bounds1[1] >= upper[iColumn1]);
            }
            bool inBetween = false;
#ifdef CBC_PRINT2
            printf("%d (way %d) with down bounds %g, %g and up bounds %g, %g current bounds %g, %g solution %g dj %g (bleft %d)\n",
              iColumn1, way1, objectI->downBounds()[0], objectI->downBounds()[1],
              objectI->upBounds()[0], objectI->upBounds()[1],
              lower[iColumn1], upper[iColumn1], solution[iColumn1],
              dj[iColumn1], nodeInfo->numberBranchesLeft());
#endif
            while (nodeInfo->parent()) {
              nodeInfo = nodeInfo->parent();
              CbcNode *nodeLook = nodeInfo->mutableOwner();
              if (!nodeLook || nodeLook->objectiveValue() == 0.5 * COIN_DBL_MAX)
                continue;
              OsiBranchingObject *obj = nodeLook->modifiableBranchingObject();
              CbcIntegerBranchingObject *objectI = dynamic_cast< CbcIntegerBranchingObject * >(obj);
              //const OsiObject * object2a = obj->originalObject();
              //assert (object2a);
              const CbcSimpleInteger *object2 = dynamic_cast< const CbcSimpleInteger * >(objectI->object());
              if (nodeInfo->numberBranchesLeft() && object2) {
                int iColumn2 = object2->columnNumber();
                double originalLower = object2->originalLowerBound();
                double originalUpper = object2->originalUpperBound();
                bool zeroOne2 = originalLower == 0.0 && originalUpper == 1.0;
                zeroOne1 = true; // temp
                double newUpper = originalUpper;
                double newLower = originalLower;
                //double value = solution[iColumn2];
                double djValue = dj[iColumn2] * direction;
                int way = objectI->way();
                assert(way == -1 || way == 1);
                way = -way; // what last branch did
#ifdef CBC_PRINT2
                printf("%d (way %d) with down bounds %g, %g and up bounds %g, %g current bounds %g, %g solution %g dj %g (bleft %d)\n",
                  iColumn2, way, objectI->downBounds()[0], objectI->downBounds()[1],
                  objectI->upBounds()[0], objectI->upBounds()[1],
                  lower[iColumn2], upper[iColumn2], solution[iColumn2],
                  djValue, nodeInfo->numberBranchesLeft());
#endif
                /*if (objectI->downBounds()[0]==0&&objectI->downBounds()[1]==1&&
                                    objectI->upBounds()[0]==0&&objectI->upBounds()[1]==1)
                                    assert(lower[iColumn2]<upper[iColumn2]);*/
                if (way < 0) {
                  // must have been down branch
                  const double *bounds = objectI->downBounds();
                  if (djValue > 1.0e-3 || solution[iColumn2] < upper[iColumn2] - 1.0e-5) {
                    if (canDelete) {
                      //nRedundantDown++;
#ifndef JJF_ONE
                      COIN_DETAIL_PRINT(printf("%d redundant branch down with bounds %g, %g current upper %g solution %g dj %g\n",
                        iColumn2, bounds[0], bounds[1], upper[iColumn2], solution[iColumn2], djValue));
#endif
                      if (bounds[0] == bounds[1] || zeroOne2 || (bounds[0] == lower[iColumn2] && false)) {
                        {
                          // get rid of node as far as branching
                          nodeLook->setObjectiveValue(0.5 * COIN_DBL_MAX);
                          objectI->deactivate();
                        }
                        previousBounds(node, nodeInfo, iColumn2, newLower, newUpper, 2);
                        solver_->setColUpper(iColumn2, newUpper);
                        assert(newLower == lower[iColumn2]);
                      } else {
                        COIN_DETAIL_PRINT(printf("SKipping\n"));
                      }
                    } else if (iColumn1 >= 0 && iColumn1 != iColumn2 && (!inBetween || true) && zeroOne1 && zeroOne2 && false) {
#ifndef JJF_ONE
                      if (true) {
                        // add in bounds
                        newLower = bounds1[0];
                        newUpper = bounds1[1];
                        COIN_DETAIL_PRINT(printf("setting bounds of %g and %g (column %d) on other branch for column %d\n",
                          newLower, newUpper, iColumn1, iColumn2));
                        int infeasible = objectI->applyExtraBounds(iColumn1, newLower, newUpper, objectI->way());
                        if (infeasible) {
                          COIN_DETAIL_PRINT(printf("infeasa!\n"));
                          // get rid of node as far as branching
                          nodeLook->setObjectiveValue(0.5 * COIN_DBL_MAX);
                        }
                      }
#endif
                    }
                    //break;
                  } else {
                    inBetween = true;
                  }
                } else {
                  // must have been up branch
                  const double *bounds = objectI->upBounds();
                  if (djValue < -1.0e-3 || solution[iColumn2] > lower[iColumn2] + 1.0e-5) {
                    if (canDelete) {
                      //nRedundantUp++;
#ifndef JJF_ONE
                      COIN_DETAIL_PRINT(printf("%d redundant branch up with bounds %g, %g current lower %g solution %g dj %g\n",
                        iColumn2, bounds[0], bounds[1], lower[iColumn2], solution[iColumn2], djValue));
#endif
                      if (bounds[0] == bounds[1] || zeroOne2 || (bounds[1] == upper[iColumn2] && false)) {
                        {
                          // get rid of node as far as branching
                          nodeLook->setObjectiveValue(0.5 * COIN_DBL_MAX);
                          objectI->deactivate();
                        }
                        previousBounds(node, nodeInfo, iColumn2, newLower, newUpper, 1);
                        solver_->setColLower(iColumn2, newLower);
                        assert(newUpper == upper[iColumn2]);
                      } else {
                        COIN_DETAIL_PRINT(printf("SKipping\n"));
                      }
                    } else if (iColumn1 >= 0 && iColumn1 != iColumn2 && (!inBetween || true) && zeroOne1 && zeroOne2 && false) {
#ifndef JJF_ONE
                      // add in bounds
                      newLower = bounds1[0];
                      newUpper = bounds1[1];
                      COIN_DETAIL_PRINT(printf("setting bounds of %g and %g (column %d) on other branch for column %d\n",
                        newLower, newUpper, iColumn1, iColumn2));
                      int infeasible = objectI->applyExtraBounds(iColumn1, newLower, newUpper, objectI->way());
                      if (infeasible) {
                        COIN_DETAIL_PRINT(printf("infeasb!\n"));
                        // get rid of node as far as branching
                        nodeLook->setObjectiveValue(0.5 * COIN_DBL_MAX);
                      }
#endif
                    }
                    // break;
                  } else {
                    inBetween = true;
                  }
                }
              } else {
                // odd
                break;
              }
            }
          }
          delete[] currentLower;
          delete[] currentUpper;
        }
      }
#endif
      if (parallelMode() >= 0)
        newNode = new CbcNode();
#if 0
	    // Try diving
            if (parallelMode() >= 0 && (specialOptions_&2048) == 0) {
	      // See if any diving heuristics set to do dive+save
	      CbcHeuristicDive * dive=NULL;
	      for (int i = 0; i < numberHeuristics_; i++) {
		CbcHeuristicDive * possible = dynamic_cast<CbcHeuristicDive *>(heuristic_[i]);
		if (possible&&possible->maxSimplexIterations()==COIN_INT_MAX) {
		  // if more than one then rotate later?
		  //if (possible->canHeuristicRun()) {
		  if (node->depth()==0||node->depth()==5) {
		    dive=possible;
		    break;
		  }
		}
	      }
	      if (dive) {
		int numberNodes;
		CbcSubProblem ** nodes=NULL;
		int branchState=dive->fathom(this,numberNodes,nodes);
		if (branchState) {
		  printf("new solution\n");
		}
		if (0) {
		  for (int iNode=0;iNode<numberNodes;iNode++) {
		    //tree_->push(nodes[iNode]) ;
		  }
		  assert (node->nodeInfo());
		  if (node->nodeInfo()->numberBranchesLeft()) {
		    tree_->push(node) ;
		  } else {
		    node->setActive(false);
		  }
		}
		delete [] nodes;
	      }
	    }
	    // end try diving
#endif
      // Set objective value (not so obvious if NLP etc)
      setObjectiveValue(newNode, node);
      int anyAction = -1;
      bool resolved = false;
      if (newNode->objectiveValue() >= getCutoff()) {
        anyAction = -2;
      } else { // only allow at most a few passes
        int numberPassesLeft = 5;
        checkingNode = true;
        OsiSolverBranch *branches = NULL;
        // point to useful information
        anyAction = chooseBranch(newNode, numberPassesLeft, node, cuts, resolved,
          lastws, lowerBefore, upperBefore, branches);
      }
      /*
            If we end up infeasible, we can delete the new node immediately. Since this
            node won't be needing the cuts we collected, decrement the reference counts.
            If we are feasible, then we'll be placing this node into the live set, so
            increment the reference count in the current (parent) nodeInfo.
            */
      lockThread();
      if (anyAction == -2) {
        if (parallelMode() > 0) {
          assert(masterThread_);
          assert(node->nodeInfo());
          node->nodeInfo()->decrement();
          delete newNode;
          assert(node->nodeInfo());
          node->nodeInfo()->increment();
          newNode = NULL;
        } else if (parallelMode() == 0) {
          delete newNode;
          newNode = NULL;
        } else {
          //assert (newNode->active());
          newNode->setActive(false);
        }
        // say strong doing well
        if (checkingNode)
          setSpecialOptions(specialOptions_ | 8);
        for (i = 0; i < currentNumberCuts_; i++) {
          if (addedCuts_[i]) {
            if (!addedCuts_[i]->decrement(1)) {
              delete addedCuts_[i];
            }
            addedCuts_[i] = NULL;
            //}
          }
        }
      } else {
        assert(node->nodeInfo());
        if (parallelMode() >= 0)
          node->nodeInfo()->increment();
        if ((numberNodes_ % 20) == 0) {
          // say strong not doing as well
          setSpecialOptions(specialOptions_ & ~8);
        }
      }
      unlockThread();
    }
    /*
          At this point, there are three possibilities:
          * newNode is live and will require further branching to resolve
          (variable() >= 0). Increment the cut reference counts by
          numberBranches() to allow for use by children of this node, and
          decrement by 1 because we've executed one arm of the branch of our
          parent (consuming one reference). Before we push newNode onto the
          search tree, try for a heuristic solution.
          * We have a solution, in which case newNode is non-null but we have no
          branching variable. Decrement the cut counts and save the solution.
          * The node was found to be infeasible, in which case it's already been
          deleted, and newNode is null.
        */
    if (eventHandler_ && !eventHandler_->event(CbcEventHandler::node)) {
      eventHappened_ = true; // exit
    }
    if (parallelMode() >= 0)
      assert(!newNode || newNode->objectiveValue() <= getCutoff());
    else
      assert(!newNode->active() || newNode->objectiveValue() <= getCutoff());
    if (statistics_) {
      assert(numberNodes2_);
      assert(statistics_[numberNodes2_ - 1]);
      assert(statistics_[numberNodes2_ - 1]->node() == numberNodes2_ - 1);
      if (newNode && newNode->active())
        statistics_[numberNodes2_ - 1]->updateInfeasibility(newNode->numberUnsatisfied());
      else
        statistics_[numberNodes2_ - 1]->sayInfeasible();
    }
    lockThread();
    bool locked = true;
    if (parallelMode() <= 0) {
      if (numberUpdateItems_) {
        for (i = 0; i < numberUpdateItems_; i++) {
          CbcObjectUpdateData *update = updateItems_ + i;
          CbcObject *object = dynamic_cast< CbcObject * >(update->object_);
#ifndef NDEBUG
          bool found = false;
          for (int j = 0; j < numberObjects_; j++) {
            if (update->object_ == object_[j]) {
              found = true;
              break;
            }
          }
          assert(found);
#endif
          if (object)
            object->updateInformation(*update);
        }
        numberUpdateItems_ = 0;
      }
    }
    if (newNode)
      if (newNode && newNode->active()) {
        if (newNode->branchingObject() == NULL) {
          const double *solution = solver_->getColSolution();
          CbcEventHandler::CbcAction action = dealWithEventHandler(CbcEventHandler::beforeSolution1,
            getSolverObjValue(), solution);
          if (action == CbcEventHandler::addCuts || solverCharacteristics_->solverType() == 4) {
            // need to check if any cuts would do anything
            OsiCuts theseCuts;
            // reset probing info
            //if (probingInfo_)
            //probingInfo_->initializeFixing(solver_);
            for (int i = 0; i < numberCutGenerators_; i++) {
              bool generate = generator_[i]->normal();
              // skip if not optimal and should be (maybe a cut generator has fixed variables)
              if (generator_[i]->needsOptimalBasis() && !solver_->basisIsAvailable())
                generate = false;
              if (!generator_[i]->mustCallAgain())
                generate = false; // only special cuts
              if (generate) {
                generator_[i]->generateCuts(theseCuts, -1, solver_, NULL);
                int numberRowCutsAfter = theseCuts.sizeRowCuts();
                if (numberRowCutsAfter)
                  break;
              }
            }
            int numberRowCutsAfter = theseCuts.sizeRowCuts();
            if (numberRowCutsAfter || action == CbcEventHandler::addCuts) {
              // need dummy branch
              newNode->setBranchingObject(new CbcDummyBranchingObject(this));
              newNode->nodeInfo()->initializeInfo(1);
            }
          }
        }
        if (newNode->branchingObject()) {
          handler_->message(CBC_BRANCH, messages_)
            << numberNodes_ << newNode->objectiveValue()
            << newNode->numberUnsatisfied() << newNode->depth()
            << CoinMessageEol;
          // Increment cut counts (taking off current)
          int numberLeft = newNode->numberBranches();
          for (i = 0; i < currentNumberCuts_; i++) {
            if (addedCuts_[i]) {
#ifdef CHECK_CUT_COUNTS
              printf("Count on cut %x increased by %d\n", addedCuts_[i],
                numberLeft - 1);
#endif
              addedCuts_[i]->increment(numberLeft - 1);
            }
          }
          unlockThread();
          locked = false;
          double estValue = newNode->guessedObjectiveValue();
          int found = -1;
          double *newSolution = new double[numberColumns];
          double heurValue = getCutoff();
          int iHeur;
          int whereFrom = 3;
          // allow more heuristics
          currentPassNumber_ = 0;
          for (iHeur = 0; iHeur < numberHeuristics_; iHeur++) {
            // skip if can't run here
            if (!heuristic_[iHeur]->shouldHeurRun(whereFrom))
              continue;
            double saveValue = heurValue;
            int ifSol = heuristic_[iHeur]->solution(heurValue, newSolution);
            if (ifSol > 0) {
              // new solution found
              heuristic_[iHeur]->incrementNumberSolutionsFound();
              found = iHeur;
              if (parallelMode() > 0) {
                lockThread();
                baseModel->incrementUsed(newSolution);
                unlockThread();
              } else {
                lastHeuristic_ = heuristic_[found];
#ifdef HEURISTIC_INFORM
                printf("HEUR %s where %d D\n",
                  lastHeuristic_->heuristicName(), whereFrom);
#endif
                setBestSolution(CBC_ROUNDING, heurValue, newSolution);
                foundSolution = 1;
                whereFrom |= 8; // say solution found
              }
            } else if (ifSol < 0) { // just returning an estimate
              estValue = heurValue; //CoinMin(heurValue, estValue) ;
              heurValue = saveValue;
            }
          }
          if (found >= 0 && parallelMode() > 0) {
            lastHeuristic_ = heuristic_[found];
#if CBC_USEFUL_PRINTING > 1
            printf("HEUR %s where %d D\n",
              lastHeuristic_->heuristicName(), whereFrom);
#endif
            setBestSolution(CBC_ROUNDING, heurValue, newSolution);
            foundSolution = 1;
          }
          delete[] newSolution;
          newNode->setGuessedObjectiveValue(estValue);
          if (parallelMode() >= 0) {
            if (!masterThread_) // only if serial
              tree_->push(newNode);
          }
          if (statistics_) {
            if (numberNodes2_ == maximumStatistics_) {
              maximumStatistics_ = 2 * maximumStatistics_;
              CbcStatistics **temp = new CbcStatistics *[maximumStatistics_];
              memset(temp, 0, maximumStatistics_ * sizeof(CbcStatistics *));
              memcpy(temp, statistics_, numberNodes2_ * sizeof(CbcStatistics *));
              delete[] statistics_;
              statistics_ = temp;
            }
            assert(!statistics_[numberNodes2_]);
            statistics_[numberNodes2_] = new CbcStatistics(newNode, this);
          }
          numberNodes2_++;
#ifdef CHECK_NODE
          printf("Node %x pushed on tree c\n", newNode);
#endif
        } else {
          if (solverCharacteristics_ && //we may be in a non standard bab
            solverCharacteristics_->solutionAddsCuts() // we are in some kind of OA based bab.
          ) {

            std::cerr << "You should never get here" << std::endl;
            throw CoinError("Nodes should not be fathomed on integer infeasibility in this setting",
              "branchAndBound", "CbcModel");
          }
          for (i = 0; i < currentNumberCuts_; i++) {
            if (addedCuts_[i]) {
              if (!addedCuts_[i]->decrement(1)) {
                delete addedCuts_[i];
                addedCuts_[i] = NULL;
              }
            }
          }
          double objectiveValue = newNode->objectiveValue();
          lastHeuristic_ = NULL;
          // Just possible solver did not know about a solution from another thread!
          if (objectiveValue < getCutoff()) {
            incrementUsed(solver_->getColSolution());
            setBestSolution(CBC_SOLUTION, objectiveValue,
              solver_->getColSolution());
            // Check if was found
            if (bestObjective_ < getCutoff())
              foundSolution = 1;
          }
          //assert(nodeInfo->numberPointingToThis() <= 2) ;
          if (parallelMode() >= 0) {
            // avoid accidental pruning, if newNode was final branch arm
            node->nodeInfo()->increment();
            delete newNode;
            newNode = NULL;
            node->nodeInfo()->decrement();
          } else {
            newNode->setActive(false);
          }
        }
      }
    if (branchesLeft) {
      // set nodenumber correctly
      if (node->nodeInfo())
        node->nodeInfo()->setNodeNumber(numberNodes2_);
      if (parallelMode() >= 0) {
        if (!masterThread_) // only if serial
          tree_->push(node);
      }
      if (statistics_) {
        if (numberNodes2_ == maximumStatistics_) {
          maximumStatistics_ = 2 * maximumStatistics_;
          CbcStatistics **temp = new CbcStatistics *[maximumStatistics_];
          memset(temp, 0, maximumStatistics_ * sizeof(CbcStatistics *));
          memcpy(temp, statistics_, numberNodes2_ * sizeof(CbcStatistics *));
          delete[] statistics_;
          statistics_ = temp;
        }
        assert(!statistics_[numberNodes2_]);
        statistics_[numberNodes2_] = new CbcStatistics(node, this);
      }
      numberNodes2_++;
      //nodeOnTree=true; // back on tree
      //deleteNode = false ;
#ifdef CHECK_NODE
      printf("Node %x pushed back on tree - %d left, %d count\n", node,
        node->nodeInfo()->numberBranchesLeft(),
        node->nodeInfo()->numberPointingToThis());
#endif
      if (parallelMode() > 0) {
        assert(node->nodeInfo());
        node->nodeInfo()->decrement();
      }
    } else {
      /*
              This node has been completely expanded and can be removed from the live
              set.
            */
      if (parallelMode() > 0) {
        assert(masterThread_);
        assert(node->nodeInfo());
        node->nodeInfo()->decrement();
      }
      assert(node->nodeInfo());
      if (parallelMode() >= 0) {
        if (!node->nodeInfo()->numberBranchesLeft())
          node->nodeInfo()->allBranchesGone(); // can clean up
        deleteNode(node);
        node = NULL;
      } else {
        node->setActive(false);
      }
    }
    if (locked)
      unlockThread();
  } else {
    // add cuts found to be infeasible (on bound)!
    COIN_DETAIL_PRINT(printf("found to be infeas! - branches left %d - cutoff %g\n", node->nodeInfo()->numberBranchesLeft(),
      getCutoff()));
#ifdef COIN_DETAIL
    node->print();
#endif
    //abort();
    assert(node->nodeInfo());
    if (parallelMode() >= 0) {
      if (!node->nodeInfo()->numberBranchesLeft())
        node->nodeInfo()->allBranchesGone(); // can clean up
      delete node;
      node = NULL;
    } else {
      node->setActive(false);
    }
  }
  /*
      Delete cuts to get back to the original system.

      I'm thinking this is redundant --- the call to addCuts that conditions entry
      to this code block also performs this action.
    */
#ifndef JJF_ONE
  //if (numberThreads_)
  {
    int numberToDelete = getNumRows() - numberRowsAtContinuous_;
    if (numberToDelete) {
      int *delRows = new int[numberToDelete];
      int i;
      for (i = 0; i < numberToDelete; i++)
        delRows[i] = i + numberRowsAtContinuous_;
      solver_->deleteRows(numberToDelete, delRows);
      delete[] delRows;
    }
    numberNewCuts_ = 0;
  }
#endif
  delete lastws;
  delete[] lowerBefore;
  delete[] upperBefore;
  if (bestObjective > bestObjective_)
    foundSolution = 2;
  if (parallelMode() > 0 && foundSolution) {
    lockThread();
    // might as well mark all including continuous
    int numberColumns = solver_->getNumCols();
    for (int i = 0; i < numberColumns; i++) {
      baseModel->usedInSolution_[i] += usedInSolution_[i];
      usedInSolution_[i] = 0;
    }
    if (bestObjective_ < baseModel->bestObjective_ && bestObjective_ < baseModel->getCutoff()) {
      baseModel->bestObjective_ = bestObjective_;
      int numberColumns = solver_->getNumCols();
      if (!baseModel->bestSolution_)
        baseModel->bestSolution_ = new double[numberColumns];
      CoinCopyN(bestSolution_, numberColumns, baseModel->bestSolution_);
      baseModel->setCutoff(getCutoff());
      baseModel->handler_->message(CBC_ROUNDING, messages_)
        << bestObjective_
        << "heuristic"
        << baseModel->numberIterations_
        << baseModel->numberNodes_ << getCurrentSeconds()
        << CoinMessageEol;
    }
    baseModel->numberSolutions_++;
    unlockThread();
  }
  numberCutGenerators_ = saveNumberCutGenerators;
  return foundSolution;
}
// Adds an update information object
void CbcModel::addUpdateInformation(const CbcObjectUpdateData &data)
{
  if (numberUpdateItems_ == maximumNumberUpdateItems_) {
    maximumNumberUpdateItems_ += 10;
    CbcObjectUpdateData *temp = new CbcObjectUpdateData[maximumNumberUpdateItems_];
    for (int i = 0; i < maximumNumberUpdateItems_ - 10; i++)
      temp[i] = updateItems_[i];
    delete[] updateItems_;
    updateItems_ = temp;
  }
  updateItems_[numberUpdateItems_++] = data;
}
// Returns bounds just before where - initially original bounds - also sets bounds
void CbcModel::previousBounds(CbcNode *node, CbcNodeInfo *where, int iColumn,
  double &lower, double &upper, int force)
{
  int i;
  int nNode = 0;
  CbcNodeInfo *nodeInfo = node->nodeInfo();
  int nWhere = -1;

  /*
      Accumulate the path from node to the root in walkback_
    */
  while (nodeInfo) {
    //printf("nNode = %d, nodeInfo = %x\n",nNode,nodeInfo);
    walkback_[nNode++] = nodeInfo;
    nodeInfo = nodeInfo->parent();
    if (nNode == maximumDepth_) {
      redoWalkBack();
    }
    if (nodeInfo == where)
      nWhere = nNode;
  }
  assert(nWhere >= 0);
  nWhere = nNode - nWhere;
  for (i = 0; i < nWhere; i++) {
    --nNode;
    walkback_[nNode]->applyBounds(iColumn, lower, upper, 0);
  }
  // correct bounds
  walkback_[nNode]->applyBounds(iColumn, lower, upper, 3);
  CbcNode *nodeLook = walkback_[nNode]->mutableOwner();
  if (nodeLook) {
    OsiBranchingObject *obj = nodeLook->modifiableBranchingObject();
    CbcIntegerBranchingObject *objectI = dynamic_cast< CbcIntegerBranchingObject * >(obj);
    //const OsiObject * object2 = obj->orig
#ifndef NDEBUG
    const CbcSimpleInteger *object2 = dynamic_cast< const CbcSimpleInteger * >(objectI->object());
    assert(object2);
    assert(iColumn == object2->columnNumber());
#endif
    double bounds[2];
    bounds[0] = lower;
    bounds[1] = upper;
    objectI->setDownBounds(bounds);
    objectI->setUpBounds(bounds);
  }
  while (nNode) {
    --nNode;
    walkback_[nNode]->applyBounds(iColumn, lower, upper, force);
#ifdef JJF_ZERO
    CbcNode *nodeLook = walkback_[nNode]->mutableOwner();
    if (nodeLook) {
      const OsiBranchingObject *obj = nodeLook->branchingObject();
      const CbcIntegerBranchingObject *objectI = dynamic_cast< const CbcIntegerBranchingObject * >(obj);
      //const OsiObject * object2 = obj->orig
      const CbcSimpleInteger *object2 = dynamic_cast< const CbcSimpleInteger * >(objectI->object());
      assert(object2);
      int iColumn2 = object2->columnNumber();
      assert(iColumn != iColumn2);
    }
#endif
  }
}
/* Return pseudo costs
   If not all integers or not pseudo costs - returns all zero
   Length of arrays are numberIntegers() and entries
      correspond to integerVariable()[i]
      User must allocate arrays before call
*/
void CbcModel::fillPseudoCosts(double *downCosts, double *upCosts,
  int *priority,
  int *numberDown, int *numberUp,
  int *numberDownInfeasible,
  int *numberUpInfeasible) const
{
  CoinFillN(downCosts, numberIntegers_, 1.0);
  CoinFillN(upCosts, numberIntegers_, 1.0);
  if (priority) {
    CoinFillN(priority, numberIntegers_, 1000000);
  }
  if (numberDown) {
    CoinFillN(numberDown, numberIntegers_, 1);
    CoinFillN(numberUp, numberIntegers_, 1);
  }
  if (numberDownInfeasible) {
    CoinZeroN(numberDownInfeasible, numberIntegers_);
    CoinZeroN(numberUpInfeasible, numberIntegers_);
  }
  int numberColumns = getNumCols();
  int *back = new int[numberColumns];
  int i;
  for (i = 0; i < numberColumns; i++)
    back[i] = -1;
  for (i = 0; i < numberIntegers_; i++)
    back[integerVariable_[i]] = i;
#if CBC_USEFUL_PRINTING > 1
  int numberNot = 0;
#endif
  for (i = 0; i < numberObjects_; i++) {
    CbcSimpleIntegerDynamicPseudoCost *obj = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(object_[i]);
    if (!obj)
      continue;
#if CBC_USEFUL_PRINTING > 1
    if (obj->numberTimesDown() < numberBeforeTrust_ || obj->numberTimesUp() < numberBeforeTrust_)
      numberNot++;
#endif
    int iColumn = obj->columnNumber();
    iColumn = back[iColumn];
    assert(iColumn >= 0);
    if (priority)
      priority[iColumn] = obj->priority();
    downCosts[iColumn] = obj->downDynamicPseudoCost();
    upCosts[iColumn] = obj->upDynamicPseudoCost();
    if (numberDown) {
      numberDown[iColumn] = obj->numberTimesDown();
      numberUp[iColumn] = obj->numberTimesUp();
    }
    if (numberDownInfeasible) {
      numberDownInfeasible[iColumn] = obj->numberTimesDownInfeasible();
      numberUpInfeasible[iColumn] = obj->numberTimesUpInfeasible();
    }
  }
#if CBC_USEFUL_PRINTING > 5
  if (priority)
    printf("Before fathom %d not trusted out of %d\n",
      numberNot, numberIntegers_);
#endif
  delete[] back;
}
// Redo walkback arrays
void CbcModel::redoWalkBack()
{
  int nNode = maximumDepth_;
  maximumDepth_ *= 2;
  CbcNodeInfo **temp = new CbcNodeInfo *[maximumDepth_];
  CbcNodeInfo **temp2 = new CbcNodeInfo *[maximumDepth_];
  int *temp3 = new int[maximumDepth_];
  for (int i = 0; i < nNode; i++) {
    temp[i] = walkback_[i];
    temp2[i] = lastNodeInfo_[i];
    temp3[i] = lastNumberCuts_[i];
  }
  delete[] walkback_;
  walkback_ = temp;
  delete[] lastNodeInfo_;
  lastNodeInfo_ = temp2;
  delete[] lastNumberCuts_;
  lastNumberCuts_ = temp3;
}
/* Return true if we want to do cuts
   If allowForTopOfTree zero then just does on multiples of depth
   if 1 then allows for doing at top of tree
   if 2 then says if cuts allowed anywhere apart from root
   if 3 then gives smallest valid depth >shallow
*/
bool CbcModel::doCutsNow(int allowForTopOfTree) const
{
  int whenCutsUse = whenCuts_;
  int alwaysReturnAt10 = whenCutsUse % 100000;
  if (whenCutsUse > 0 && alwaysReturnAt10) {
    whenCutsUse -= alwaysReturnAt10;
    if (currentDepth_ > 10)
      return false;
  }
  //if (currentDepth_>10)
  //return false;
#define TRY_IDEA1 2
  int size = continuousSolver_->getNumRows() + continuousSolver_->getNumCols();

  if (true && (whenCutsUse < 0 || (size <= 500 - 500 * TRY_IDEA1 && allowForTopOfTree != 3))) {
    int whenCuts = (size <= 500) ? -1 : 1;
    //whenCuts = (size<=500) ? 1 :1;
    if (parentModel_)
      whenCuts = 1;
    //int nodeDepth = currentDepth_-1;
    bool doCuts2 = !(currentDepth_ > 11 && (currentDepth_ & 1) == whenCuts);
    if (fastNodeDepth_ > 0 && currentDepth_ > 10)
      doCuts2 = false;
    //printf("when %d node %d depth %d size %d doing cuts %s\n",whenCutsUse,
    //   numberNodes_,currentDepth_,size,doCuts2 ? "yes" : "no");
    return doCuts2;
  }
  //if (!parentModel_&&currentDepth_==7)
  //printf("q\n");
  int top = whenCutsUse / 1000000;
  int shallow = top ? (top - 1) : 9;
  int when = whenCutsUse - 1000000 * top;
#if TRY_IDEA1
  if (when < 15 && when > 1 && size <= 500)
    when /= 2;
#endif
  if ((when > 15 || (top && top < 5)) && currentDepth_ > when)
    when = 100000; // off
  bool doCuts = when ? ((currentDepth_ % when) == 0) || (when == 1) : false;
  if (allowForTopOfTree == 1 && currentDepth_ <= shallow) {
    doCuts = true;
  } else if (allowForTopOfTree == 2 && shallow >= 1) {
    doCuts = true;
#if TRY_IDEA1 < 2
  } else if (allowForTopOfTree == 3 && doCuts) {
    // only if first
    if (currentDepth_ <= shallow || currentDepth_ - when > shallow)
      doCuts = false;
#else
  } else if (allowForTopOfTree == 3) {
    // only exactly at 10
    doCuts = (currentDepth_ == 10);
#endif
  }
  //if (!doCuts&&currentDepth_&&!parentModel_)
  //printf("zzz\n");
  return doCuts;
}
// See if can stop on gap
bool CbcModel::canStopOnGap() const
{
  bool returnCode = false;
  if (bestObjective_ < 1.0e50) {
    double testGap = CoinMax(dblParam_[CbcAllowableGap],
      CoinMax(fabs(bestObjective_), fabs(bestPossibleObjective_))
        * dblParam_[CbcAllowableFractionGap]);
    returnCode = (bestObjective_ - bestPossibleObjective_ < testGap && getCutoffIncrement() >= 0.0);
  }
#if 0
  if (returnCode) {
    if (fabs(bestObjective_+1469650.0)<1.0) {
      fprintf(stderr,"BAD - cr to continue\n");
      fflush(stdout);
      char xx;
      xx=getc(stdin);
    }
  }
#endif
  return returnCode;
}
// Adjust heuristics based on model
void CbcModel::adjustHeuristics()
{
  int numberRows = solver_->getNumRows();
  int numberColumns = solver_->getNumCols();
  int nTree = CoinMax(10000, 2 * numberRows + numberColumns);
  int nRoot = CoinMax(40000, 8 * numberRows + 4 * numberColumns);
  for (int i = 0; i < numberHeuristics_; i++) {
    CbcHeuristicDive *heuristic = dynamic_cast< CbcHeuristicDive * >(heuristic_[i]);
    if (heuristic && heuristic->maxSimplexIterations() != COIN_INT_MAX) {
      heuristic->setMaxSimplexIterations(nTree);
      heuristic->setMaxSimplexIterationsAtRoot(nRoot);
    }
  }
}
// Number of saved solutions (including best)
int CbcModel::numberSavedSolutions() const
{
  if (!bestSolution_)
    return 0;
  else
    return numberSavedSolutions_ + 1;
}
// Set maximum number of extra saved solutions
void CbcModel::setMaximumSavedSolutions(int value)
{
  if (value < maximumSavedSolutions_) {
    for (int i = value; i < maximumSavedSolutions_; i++)
      delete[] savedSolutions_[i];
    maximumSavedSolutions_ = value;
    numberSavedSolutions_ = CoinMin(numberSavedSolutions_,
      maximumSavedSolutions_);
    if (!maximumSavedSolutions_)
      delete[] savedSolutions_;
  } else if (value > maximumSavedSolutions_) {
    double **temp = new double *[value];
    int i;
    for (i = 0; i < maximumSavedSolutions_; i++)
      temp[i] = savedSolutions_[i];
    for (; i < value; i++)
      temp[i] = NULL;
    delete[] savedSolutions_;
    maximumSavedSolutions_ = value;
    savedSolutions_ = temp;
  }
}
// Return a saved solution objective (0==best) - COIN_DBL_MAX if off end
double
CbcModel::savedSolutionObjective(int which) const
{
  if (which == 0) {
    return bestObjective_;
  } else if (which <= numberSavedSolutions_) {
    double *sol = savedSolutions_[which - 1];
    assert(static_cast< int >(sol[0]) == solver_->getNumCols());
    return sol[1];
  } else {
    return COIN_DBL_MAX;
  }
}
// Return a saved solution (0==best) - NULL if off end
const double *
CbcModel::savedSolution(int which) const
{
  if (which == 0) {
    return bestSolution_;
  } else if (which <= numberSavedSolutions_) {
    double *sol = savedSolutions_[which - 1];
    assert(static_cast< int >(sol[0]) == solver_->getNumCols());
    return sol + 2;
  } else {
    return NULL;
  }
}
// Save a solution
void CbcModel::saveExtraSolution(const double *solution, double objectiveValue)
{
  double *save = NULL;
  if (maximumSavedSolutions_) {
    if (!savedSolutions_) {
      savedSolutions_ = new double *[maximumSavedSolutions_];
      for (int i = 0; i < maximumSavedSolutions_; i++)
        savedSolutions_[i] = NULL;
    }
    int n = solver_->getNumCols();
    int k;
    for (k = numberSavedSolutions_ - 1; k >= 0; k--) {
      double *sol = savedSolutions_[k];
      assert(static_cast< int >(sol[0]) == n);
      if (objectiveValue > sol[1])
        break;
    }
    k++; // where to put
    if (k < maximumSavedSolutions_) {
      if (numberSavedSolutions_ == maximumSavedSolutions_) {
        save = savedSolutions_[numberSavedSolutions_ - 1];
      } else {
        save = new double[n + 2];
        numberSavedSolutions_++;
      }
      // move up
      for (int j = maximumSavedSolutions_ - 1; j > k; j--)
        savedSolutions_[j] = savedSolutions_[j - 1];
      savedSolutions_[k] = save;
      save[0] = n;
      save[1] = objectiveValue;
      memcpy(save + 2, solution, n * sizeof(double));
    }
  }
}
// Save a solution to best and move current to saved
void CbcModel::saveBestSolution(const double *solution, double objectiveValue)
{
  int n = solver_->getNumCols();
  if (bestSolution_)
    saveExtraSolution(bestSolution_, bestObjective_);
  else
    bestSolution_ = new double[n];
  bestObjective_ = objectiveValue;
  memcpy(bestSolution_, solution, n * sizeof(double));
}
// Delete best and saved solutions
void CbcModel::deleteSolutions()
{
  delete[] bestSolution_;
  bestSolution_ = NULL;
  for (int i = 0; i < maximumSavedSolutions_; i++) {
    delete[] savedSolutions_[i];
    savedSolutions_[i] = NULL;
  }
  numberSavedSolutions_ = 0;
}
// Delete a saved solution and move others up
void CbcModel::deleteSavedSolution(int which)
{
  if (which > 0 && which <= numberSavedSolutions_) {
    delete[] savedSolutions_[which - 1];
    // move up
    numberSavedSolutions_--;
    for (int j = which - 1; j < numberSavedSolutions_; j++) {
      savedSolutions_[j] = savedSolutions_[j + 1];
    }
    savedSolutions_[numberSavedSolutions_] = NULL;
  }
}
#ifdef COIN_HAS_CLP
void CbcModel::goToDantzig(int numberNodes, ClpDualRowPivot *&savePivotMethod)
{
  // Possible change of pivot method
  if (!savePivotMethod && !parentModel_) {
    OsiClpSolverInterface *clpSolver
      = dynamic_cast< OsiClpSolverInterface * >(solver_);
    if (clpSolver && numberNodes_ >= numberNodes && numberNodes_ < 2 * numberNodes && clpSolver->getNumRows() < 10000) {
      if (numberIterations_ < (numberSolves_ + numberNodes_) * 10) {
        //if (numberIterations_<numberNodes_*20) {
        ClpSimplex *simplex = clpSolver->getModelPtr();
        ClpDualRowPivot *pivotMethod = simplex->dualRowPivot();
        ClpDualRowDantzig *pivot = dynamic_cast< ClpDualRowDantzig * >(pivotMethod);
        if (!pivot) {
          savePivotMethod = pivotMethod->clone(true);
          ClpDualRowDantzig dantzig;
          simplex->setDualRowPivotAlgorithm(dantzig);
#ifdef COIN_DEVELOP
          printf("%d node, %d iterations ->Dantzig\n", numberNodes_,
            numberIterations_);
#endif
#ifdef CBC_THREAD
          if (master_)
            master_->setDantzigState();
#endif
        }
      }
    }
  }
}
#else
CbcModel::goToDantzig(int numberNodes, ClpDualRowPivot *&savePivotMethod)
{
  printf("Need Clp to go to Dantzig\n");
  abort();
}
#endif
// Below this is deprecated or at least fairly deprecated
/*
   Do Integer Presolve. Returns new model.
   I have to work out cleanest way of getting solution to
   original problem at end.  So this is very preliminary.
*/
CbcModel *
CbcModel::integerPresolve(bool weak)
{
  status_ = 0;
  // solve LP
  //solver_->writeMps("bad");
  bool feasible = (resolve(NULL, 3) != 0);

  CbcModel *newModel = NULL;
  if (feasible) {

    // get a new model
    newModel = new CbcModel(*this);
    newModel->messageHandler()->setLogLevel(messageHandler()->logLevel());

    feasible = newModel->integerPresolveThisModel(solver_, weak);
  }
  if (!feasible) {
    handler_->message(CBC_INFEAS, messages_)
      << CoinMessageEol;
    status_ = 0;
    secondaryStatus_ = 1;
    delete newModel;
    return NULL;
  } else {
    newModel->synchronizeModel(); // make sure everything that needs solver has it
    return newModel;
  }
}
/*
   Do Integer Presolve - destroying current model
*/
bool CbcModel::integerPresolveThisModel(OsiSolverInterface *originalSolver,
  bool weak)
{
  printf("DEPRECATED\n");
  status_ = 0;
  // solve LP
  bool feasible = (resolve(NULL, 3) != 0);

  bestObjective_ = 1.0e50;
  numberSolutions_ = 0;
  numberHeuristicSolutions_ = 0;
  double cutoff = getCutoff();
  double direction = solver_->getObjSense();
  if (cutoff < 1.0e20 && direction < 0.0)
    messageHandler()->message(CBC_CUTOFF_WARNING1,
      messages())
      << cutoff << -cutoff << CoinMessageEol;
  if (cutoff > bestObjective_)
    cutoff = bestObjective_;
  setCutoff(cutoff);
  int iColumn;
  int numberColumns = getNumCols();
  int originalNumberColumns = numberColumns;
  currentPassNumber_ = 0;
  synchronizeModel(); // make sure everything that needs solver has it
  if (!solverCharacteristics_) {
    OsiBabSolver *solverCharacteristics = dynamic_cast< OsiBabSolver * >(solver_->getAuxiliaryInfo());
    if (solverCharacteristics) {
      solverCharacteristics_ = solverCharacteristics;
    } else {
      // replace in solver
      OsiBabSolver defaultC;
      solver_->setAuxiliaryInfo(&defaultC);
      solverCharacteristics_ = dynamic_cast< OsiBabSolver * >(solver_->getAuxiliaryInfo());
    }
  }
  solverCharacteristics_->setSolver(solver_);
  // just point to solver_
  delete continuousSolver_;
  continuousSolver_ = solver_;
  // get a copy of original so we can fix bounds
  OsiSolverInterface *cleanModel = originalSolver->clone();
#ifdef CBC_DEBUG
  std::string problemName;
  cleanModel->getStrParam(OsiProbName, problemName);
  printf("Problem name - %s\n", problemName.c_str());
  cleanModel->activateRowCutDebugger(problemName.c_str());
  const OsiRowCutDebugger *debugger = cleanModel->getRowCutDebugger();
#endif

  // array which points from original columns to presolved
  int *original = new int[numberColumns];
  // arrays giving bounds - only ones found by probing
  // rest will be found by presolve
  double *originalLower = new double[numberColumns];
  double *originalUpper = new double[numberColumns];
  {
    const double *lower = getColLower();
    const double *upper = getColUpper();
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      original[iColumn] = iColumn;
      originalLower[iColumn] = lower[iColumn];
      originalUpper[iColumn] = upper[iColumn];
    }
  }
  findIntegers(true);
  // save original integers
  int *originalIntegers = new int[numberIntegers_];
  int originalNumberIntegers = numberIntegers_;
  memcpy(originalIntegers, integerVariable_, numberIntegers_ * sizeof(int));

  int todo = 20;
  if (weak)
    todo = 1;
  while (currentPassNumber_ < todo) {

    currentPassNumber_++;
    numberSolutions_ = 0;
    // this will be set false to break out of loop with presolved problem
    bool doIntegerPresolve = (currentPassNumber_ != 20);

    // Current number of free integer variables
    // Get increment in solutions
    {
      const double *objective = cleanModel->getObjCoefficients();
      const double *lower = cleanModel->getColLower();
      const double *upper = cleanModel->getColUpper();
      double maximumCost = 0.0;
      bool possibleMultiple = true;
      int numberChanged = 0;
      for (iColumn = 0; iColumn < originalNumberColumns; iColumn++) {
        if (originalUpper[iColumn] > originalLower[iColumn]) {
          if (cleanModel->isInteger(iColumn)) {
            maximumCost = CoinMax(maximumCost, fabs(objective[iColumn]));
          } else if (objective[iColumn]) {
            possibleMultiple = false;
          }
        }
        if (originalUpper[iColumn] < upper[iColumn]) {
#ifdef CBC_DEBUG
          printf("Changing upper bound on %d from %g to %g\n",
            iColumn, upper[iColumn], originalUpper[iColumn]);
#endif
          cleanModel->setColUpper(iColumn, originalUpper[iColumn]);
          numberChanged++;
        }
        if (originalLower[iColumn] > lower[iColumn]) {
#ifdef CBC_DEBUG
          printf("Changing lower bound on %d from %g to %g\n",
            iColumn, lower[iColumn], originalLower[iColumn]);
#endif
          cleanModel->setColLower(iColumn, originalLower[iColumn]);
          numberChanged++;
        }
      }
      // if first pass - always try
      if (currentPassNumber_ == 1)
        numberChanged += 1;
      if (possibleMultiple && maximumCost) {
        int increment = 0;
        double multiplier = 2520.0;
        while (10.0 * multiplier * maximumCost < 1.0e8)
          multiplier *= 10.0;
        for (int j = 0; j < originalNumberIntegers; j++) {
          iColumn = originalIntegers[j];
          if (originalUpper[iColumn] > originalLower[iColumn]) {
            if (objective[iColumn]) {
              double value = fabs(objective[iColumn]) * multiplier;
              int nearest = static_cast< int >(floor(value + 0.5));
              if (fabs(value - floor(value + 0.5)) > 1.0e-8 || value > 2.1e9) {
                increment = 0;
                break; // no good
              } else if (!increment) {
                // first
                increment = nearest;
              } else {
                increment = gcd(increment, nearest);
              }
            }
          }
        }
        if (increment) {
          double value = increment;
          value /= multiplier;
          if (value * 0.999 > dblParam_[CbcCutoffIncrement]) {
            messageHandler()->message(CBC_INTEGERINCREMENT, messages())
              << value
              << CoinMessageEol;
            dblParam_[CbcCutoffIncrement] = value * 0.999;
          }
        }
      }
      if (!numberChanged) {
        doIntegerPresolve = false; // not doing any better
      }
    }
#ifdef CBC_DEBUG
    if (debugger)
      assert(debugger->onOptimalPath(*cleanModel));
#endif
#ifdef COIN_HAS_CLP
    // do presolve - for now just clp but easy to get osi interface
    OsiClpSolverInterface *clpSolver
      = dynamic_cast< OsiClpSolverInterface * >(cleanModel);
    if (clpSolver) {
      ClpSimplex *clp = clpSolver->getModelPtr();
      clp->messageHandler()->setLogLevel(cleanModel->messageHandler()->logLevel());
      ClpPresolve pinfo;
      //printf("integerPresolve - temp switch off doubletons\n");
      //pinfo.setPresolveActions(4);
      ClpSimplex *model2 = pinfo.presolvedModel(*clp, 1.0e-8);
      if (!model2) {
        // presolve found to be infeasible
        feasible = false;
      } else {
        // update original array
        const int *originalColumns = pinfo.originalColumns();
        // just slot in new solver
        OsiClpSolverInterface *temp = new OsiClpSolverInterface(model2, true);
        numberColumns = temp->getNumCols();
        for (iColumn = 0; iColumn < originalNumberColumns; iColumn++)
          original[iColumn] = -1;
        for (iColumn = 0; iColumn < numberColumns; iColumn++)
          original[originalColumns[iColumn]] = iColumn;
        // copy parameters
        temp->copyParameters(*solver_);
        // and specialized ones
        temp->setSpecialOptions(clpSolver->specialOptions());
        delete solver_;
        solver_ = temp;
        setCutoff(cutoff);
        deleteObjects();
        if (!numberObjects_) {
          // Nothing left
          doIntegerPresolve = false;
          weak = true;
          break;
        }
        synchronizeModel(); // make sure everything that needs solver has it
        // just point to solver_
        continuousSolver_ = solver_;
        feasible = (resolve(NULL, 3) != 0);
        if (!feasible || !doIntegerPresolve || weak)
          break;
        // see if we can get solution by heuristics
        int found = -1;
        int iHeuristic;
        double *newSolution = new double[numberColumns];
        double heuristicValue = getCutoff();
        int whereFrom = 0;
        for (iHeuristic = 0; iHeuristic < numberHeuristics_; iHeuristic++) {
          // skip if can't run here
          if (!heuristic_[iHeuristic]->shouldHeurRun(whereFrom))
            continue;
          double saveValue = heuristicValue;
          int ifSol = heuristic_[iHeuristic]->solution(heuristicValue,
            newSolution);
          if (ifSol > 0) {
            // better solution found
            heuristic_[iHeuristic]->incrementNumberSolutionsFound();
            found = iHeuristic;
            incrementUsed(newSolution);
            whereFrom |= 8; // say solution found
          } else if (ifSol < 0) {
            heuristicValue = saveValue;
          }
        }
        if (found >= 0) {
          // We probably already have a current solution, but just in case ...
          int numberColumns = getNumCols();
          if (!currentSolution_)
            currentSolution_ = new double[numberColumns];
          testSolution_ = currentSolution_;
          // better solution save
          lastHeuristic_ = heuristic_[found];
#ifdef HEURISTIC_INFORM
          printf("HEUR %s where %d oddE\n",
            lastHeuristic_->heuristicName(), whereFrom);
#endif
          setBestSolution(CBC_ROUNDING, heuristicValue,
            newSolution);
          // update cutoff
          cutoff = getCutoff();
        }
        delete[] newSolution;
        // Space for type of cuts
        maximumWhich_ = INITIAL_MAXIMUM_WHICH;
        delete[] whichGenerator_;
        whichGenerator_ = new int[maximumWhich_];
        // save number of rows
        numberRowsAtContinuous_ = getNumRows();
        maximumNumberCuts_ = 0;
        currentNumberCuts_ = 0;
        delete[] addedCuts_;
        addedCuts_ = NULL;

        // maximum depth for tree walkback
        maximumDepth_ = 10;
        delete[] walkback_;
        walkback_ = new CbcNodeInfo *[maximumDepth_];
        lastDepth_ = 0;
        delete[] lastNodeInfo_;
        lastNodeInfo_ = new CbcNodeInfo *[maximumDepth_];
        delete[] lastNumberCuts_;
        lastNumberCuts_ = new int[maximumDepth_];
        maximumCuts_ = 100;
        delete[] lastCut_;
        lastCut_ = new const OsiRowCut *[maximumCuts_];

        OsiCuts cuts;
        numberOldActiveCuts_ = 0;
        numberNewCuts_ = 0;
        feasible = solveWithCuts(cuts, maximumCutPassesAtRoot_, NULL);
        currentNumberCuts_ = numberNewCuts_;
        delete[] whichGenerator_;
        whichGenerator_ = NULL;
        delete[] walkback_;
        walkback_ = NULL;
        delete[] addedCuts_;
        addedCuts_ = NULL;
        if (feasible) {
          // fix anything in original which integer presolve fixed
          // for now just integers
          const double *lower = solver_->getColLower();
          const double *upper = solver_->getColUpper();
          int i;
          for (i = 0; i < originalNumberIntegers; i++) {
            iColumn = originalIntegers[i];
            int jColumn = original[iColumn];
            if (jColumn >= 0) {
              if (upper[jColumn] < originalUpper[iColumn])
                originalUpper[iColumn] = upper[jColumn];
              if (lower[jColumn] > originalLower[iColumn])
                originalLower[iColumn] = lower[jColumn];
            }
          }
        }
      }
    }
#endif
    if (!feasible || !doIntegerPresolve) {
      break;
    }
  }
  //solver_->writeMps("xx");
  delete cleanModel;
  delete[] originalIntegers;
  numberColumns = getNumCols();
  delete[] originalColumns_;
  originalColumns_ = new int[numberColumns];
  numberColumns = 0;
  for (iColumn = 0; iColumn < originalNumberColumns; iColumn++) {
    int jColumn = original[iColumn];
    if (jColumn >= 0)
      originalColumns_[numberColumns++] = iColumn;
  }
  delete[] original;
  delete[] originalLower;
  delete[] originalUpper;

  deleteObjects();
  synchronizeModel(); // make sure everything that needs solver has it
  continuousSolver_ = NULL;
  currentNumberCuts_ = 0;
  return feasible;
}
// Put back information into original model - after integerpresolve
void CbcModel::originalModel(CbcModel *presolvedModel, bool weak)
{
  solver_->copyParameters(*(presolvedModel->solver_));
  bestObjective_ = presolvedModel->bestObjective_;
  delete[] bestSolution_;
  findIntegers(true);
  if (presolvedModel->bestSolution_) {
    int numberColumns = getNumCols();
    int numberOtherColumns = presolvedModel->getNumCols();
    //bestSolution_ = new double[numberColumns];
    // set up map
    int *back = new int[numberColumns];
    int i;
    for (i = 0; i < numberColumns; i++)
      back[i] = -1;
    for (i = 0; i < numberOtherColumns; i++)
      back[presolvedModel->originalColumns_[i]] = i;
    int iColumn;
    // set ones in presolved model to values
    double *otherSolution = presolvedModel->bestSolution_;
    //const double * lower = getColLower();
    for (i = 0; i < numberIntegers_; i++) {
      iColumn = integerVariable_[i];
      int jColumn = back[iColumn];
      //bestSolution_[iColumn]=lower[iColumn];
      if (jColumn >= 0) {
        double value = floor(otherSolution[jColumn] + 0.5);
        solver_->setColLower(iColumn, value);
        solver_->setColUpper(iColumn, value);
        //bestSolution_[iColumn]=value;
      }
    }
    delete[] back;
#ifdef JJF_ZERO
    // ** looks as if presolve needs more intelligence
    // do presolve - for now just clp but easy to get osi interface
    OsiClpSolverInterface *clpSolver
      = dynamic_cast< OsiClpSolverInterface * >(solver_);
    assert(clpSolver);
    ClpSimplex *clp = clpSolver->getModelPtr();
    Presolve pinfo;
    ClpSimplex *model2 = pinfo.presolvedModel(*clp, 1.0e-8);
    model2->primal(1);
    pinfo.postsolve(true);
    const double *solution = solver_->getColSolution();
    for (i = 0; i < numberIntegers_; i++) {
      iColumn = integerVariable_[i];
      double value = floor(solution[iColumn] + 0.5);
      solver_->setColLower(iColumn, value);
      solver_->setColUpper(iColumn, value);
    }
#else
    if (!weak) {
      // for now give up
      int save = numberCutGenerators_;
      numberCutGenerators_ = 0;
      bestObjective_ = 1.0e100;
      branchAndBound();
      numberCutGenerators_ = save;
    }
#endif
    if (bestSolution_) {
      // solve problem
      resolve(NULL, 3);
      // should be feasible
      if (!currentSolution_)
        currentSolution_ = new double[numberColumns];
      testSolution_ = currentSolution_;
#ifndef NDEBUG
      int numberIntegerInfeasibilities;
      int numberObjectInfeasibilities;
      assert(feasibleSolution(numberIntegerInfeasibilities,
        numberObjectInfeasibilities));
#endif
    }
  } else {
    bestSolution_ = NULL;
  }
  numberSolutions_ = presolvedModel->numberSolutions_;
  numberHeuristicSolutions_ = presolvedModel->numberHeuristicSolutions_;
  numberNodes_ = presolvedModel->numberNodes_;
  numberIterations_ = presolvedModel->numberIterations_;
  status_ = presolvedModel->status_;
  secondaryStatus_ = presolvedModel->secondaryStatus_;
  synchronizeModel();
}
void CbcModel::setOptionalInteger(int index)
{
#ifdef COIN_HAS_CLP
  OsiClpSolverInterface *clpSolver
    = dynamic_cast< OsiClpSolverInterface * >(solver_);
  if (clpSolver)
    clpSolver->setOptionalInteger(index);
  else
#endif
    solver_->setInteger(index);
}
// Return true if maximum time reached
bool CbcModel::maximumSecondsReached() const
{
  double totalTime = getCurrentSeconds();
  double maxSeconds = getMaximumSeconds();
  bool hitMaxTime = (totalTime >= maxSeconds);
  if (parentModel_ && !hitMaxTime) {
    // In a sub tree
    assert(parentModel_);
    maxSeconds = parentModel_->getMaximumSeconds();
    hitMaxTime = (totalTime >= maxSeconds);
  }
  if (hitMaxTime) {
    // Set eventHappened_ so will by-pass as much stuff as possible
    eventHappened_ = true;
  }
  return hitMaxTime;
}
// Check original model before it gets messed up
void CbcModel::checkModel()
{
  int iColumn;
  int numberColumns = getNumCols();
  const double *lower = getColLower();
  const double *upper = getColUpper();
  int setFlag = 65536;
  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    if (upper[iColumn] > lower[iColumn] + 1.0e-8) {
      double value;
      value = fabs(lower[iColumn]);
      if (floor(value + 0.5) != value) {
        setFlag = 0;
        break;
      }
      value = fabs(upper[iColumn]);
      if (floor(value + 0.5) != value) {
        setFlag = 0;
        break;
      }
    }
  }
  specialOptions_ |= setFlag;
}
static void flipSolver(OsiSolverInterface *solver, double newCutoff)
{
  if (solver) {
    double objValue = solver->getObjValue();
    double objectiveOffset;
    solver->setObjSense(-solver->getObjSense());
    solver->getDblParam(OsiObjOffset, objectiveOffset);
    solver->setDblParam(OsiObjOffset, -objectiveOffset);
    int numberColumns = solver->getNumCols();
    double *array = CoinCopyOfArray(solver->getObjCoefficients(), numberColumns);
    for (int i = 0; i < numberColumns; i++)
      array[i] = -array[i];
    solver->setObjective(array);
    delete[] array;
    solver->setDblParam(OsiDualObjectiveLimit, newCutoff);
#ifdef COIN_HAS_CLP
    OsiClpSolverInterface *clpSolver
      = dynamic_cast< OsiClpSolverInterface * >(solver);
    if (clpSolver) {
      double *dj = clpSolver->getModelPtr()->dualColumnSolution();
      for (int i = 0; i < numberColumns; i++)
        dj[i] = -dj[i];
      int numberRows = clpSolver->getNumRows();
      double *pi = clpSolver->getModelPtr()->dualRowSolution();
      for (int i = 0; i < numberRows; i++)
        pi[i] = -pi[i];
      clpSolver->getModelPtr()->setObjectiveValue(-objValue);
    } else {
#endif
      // update values
      solver->resolve();
#ifdef COIN_HAS_CLP
    }
#endif
  }
}
/*
  Flip direction of optimization on all models
*/
void CbcModel::flipModel()
{
  if (parentModel_)
    return;
  // I think cutoff is always minimization
  double cutoff = getCutoff();
  flipSolver(referenceSolver_, cutoff);
  flipSolver(continuousSolver_, cutoff);
  flipSolver(solver_, cutoff);
}
#ifdef CBC_KEEP_DEPRECATED
/* preProcess problem - replacing solver
   If makeEquality true then <= cliques converted to ==.
   Presolve will be done numberPasses times.

   Returns NULL if infeasible

   If makeEquality is 1 add slacks to get cliques,
   if 2 add slacks to get sos (but only if looks plausible) and keep sos info
*/
CglPreProcess *
CbcModel::preProcess(int makeEquality, int numberPasses, int tuning)
{
  CglPreProcess *process = new CglPreProcess();
  // Default set of cut generators
  CglProbing generator1;
  generator1.setUsingObjective(true);
  generator1.setMaxPass(3);
  generator1.setMaxProbeRoot(solver_->getNumCols());
  generator1.setMaxElements(100);
  generator1.setMaxLookRoot(50);
  generator1.setRowCuts(3);
  // Add in generators
  process->addCutGenerator(&generator1);
  process->messageHandler()->setLogLevel(this->logLevel());
  /* model may not have created objects
       If none then create
    */
  if (!numberIntegers_ || !numberObjects_) {
    this->findIntegers(true, 1);
  }
  // Do SOS
  int i;
  int numberSOS2 = 0;
  for (i = 0; i < numberObjects_; i++) {
    CbcSOS *objSOS = dynamic_cast< CbcSOS * >(object_[i]);
    if (objSOS) {
      int type = objSOS->sosType();
      if (type == 2)
        numberSOS2++;
    }
  }
  if (numberSOS2) {
    // SOS
    int numberColumns = solver_->getNumCols();
    char *prohibited = new char[numberColumns];
    memset(prohibited, 0, numberColumns);
    for (i = 0; i < numberObjects_; i++) {
      CbcSOS *objSOS = dynamic_cast< CbcSOS * >(object_[i]);
      if (objSOS) {
        int type = objSOS->sosType();
        if (type == 2) {
          int n = objSOS->numberMembers();
          const int *which = objSOS->members();
          for (int j = 0; j < n; j++) {
            int iColumn = which[j];
            prohibited[iColumn] = 1;
          }
        }
      }
    }
    process->passInProhibited(prohibited, numberColumns);
    delete[] prohibited;
  }
  // Tell solver we are not in Branch and Cut
  solver_->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo);
  OsiSolverInterface *newSolver = process->preProcessNonDefault(*solver_, makeEquality,
    numberPasses, tuning);
  // Tell solver we are not in Branch and Cut
  solver_->setHintParam(OsiDoInBranchAndCut, false, OsiHintDo);
  if (newSolver) {
    int numberOriginalObjects = numberObjects_;
    OsiSolverInterface *originalSolver = solver_;
    solver_ = newSolver->clone(); // clone as process owns solver
    // redo sequence
    numberIntegers_ = 0;
    int numberColumns = solver_->getNumCols();
    int nOrig = originalSolver->getNumCols();
    const int *originalColumns = process->originalColumns();
    // allow for cliques etc
    nOrig = CoinMax(nOrig, originalColumns[numberColumns - 1] + 1);
    OsiObject **originalObject = object_;
    // object number or -1
    int *temp = new int[nOrig];
    int iColumn;
    for (iColumn = 0; iColumn < nOrig; iColumn++)
      temp[iColumn] = -1;
    int iObject;
    numberObjects_ = 0;
    int nNonInt = 0;
    for (iObject = 0; iObject < numberOriginalObjects; iObject++) {
      iColumn = originalObject[iObject]->columnNumber();
      if (iColumn < 0) {
        nNonInt++;
      } else {
        temp[iColumn] = iObject;
      }
    }
    int numberNewIntegers = 0;
    int numberOldIntegers = 0;
    int numberOldOther = 0;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      int jColumn = originalColumns[iColumn];
      if (temp[jColumn] >= 0) {
        int iObject = temp[jColumn];
        CbcSimpleInteger *obj = dynamic_cast< CbcSimpleInteger * >(originalObject[iObject]);
        if (obj)
          numberOldIntegers++;
        else
          numberOldOther++;
      } else if (isInteger(iColumn)) {
        numberNewIntegers++;
      }
    }
    /*
          Allocate an array to hold the indices of the integer variables.
          Make a large enough array for all objects
        */
    numberObjects_ = numberNewIntegers + numberOldIntegers + numberOldOther + nNonInt;
    object_ = new OsiObject *[numberObjects_];
    delete[] integerVariable_;
    integerVariable_ = new int[numberNewIntegers + numberOldIntegers];
    /*
          Walk the variables again, filling in the indices and creating objects for
          the integer variables. Initially, the objects hold the index and upper &
          lower bounds.
        */
    numberIntegers_ = 0;
    int n = originalColumns[numberColumns - 1] + 1;
    int *backward = new int[n];
    int i;
    for (i = 0; i < n; i++)
      backward[i] = -1;
    for (i = 0; i < numberColumns; i++)
      backward[originalColumns[i]] = i;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      int jColumn = originalColumns[iColumn];
      if (temp[jColumn] >= 0) {
        int iObject = temp[jColumn];
        CbcSimpleInteger *obj = dynamic_cast< CbcSimpleInteger * >(originalObject[iObject]);
        if (obj) {
          object_[numberIntegers_] = originalObject[iObject]->clone();
          // redo ids etc
          //object_[numberIntegers_]->resetSequenceEtc(numberColumns,originalColumns);
          object_[numberIntegers_]->resetSequenceEtc(numberColumns, backward);
          integerVariable_[numberIntegers_++] = iColumn;
        }
      } else if (isInteger(iColumn)) {
        object_[numberIntegers_] = new CbcSimpleInteger(this, iColumn);
        integerVariable_[numberIntegers_++] = iColumn;
      }
    }
    delete[] backward;
    numberObjects_ = numberIntegers_;
    // Now append other column stuff
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      int jColumn = originalColumns[iColumn];
      if (temp[jColumn] >= 0) {
        int iObject = temp[jColumn];
        CbcSimpleInteger *obj = dynamic_cast< CbcSimpleInteger * >(originalObject[iObject]);
        if (!obj) {
          object_[numberObjects_] = originalObject[iObject]->clone();
          // redo ids etc
          CbcObject *obj = dynamic_cast< CbcObject * >(object_[numberObjects_]);
          assert(obj);
          obj->redoSequenceEtc(this, numberColumns, originalColumns);
          numberObjects_++;
        }
      }
    }
    // now append non column stuff
    for (iObject = 0; iObject < numberOriginalObjects; iObject++) {
      iColumn = originalObject[iObject]->columnNumber();
      if (iColumn < 0) {
        object_[numberObjects_] = originalObject[iObject]->clone();
        // redo ids etc
        CbcObject *obj = static_cast< CbcObject * >(object_[numberObjects_]);
        assert(obj);
        obj->redoSequenceEtc(this, numberColumns, originalColumns);
        numberObjects_++;
      }
      delete originalObject[iObject];
    }
    delete[] originalObject;
    delete[] temp;
    if (!numberObjects_)
      handler_->message(CBC_NOINT, messages_) << CoinMessageEol;
    return process;
  } else {
    // infeasible
    delete process;
    return NULL;
  }
}
/* Does postprocessing - original solver back.
   User has to delete process */
void CbcModel::postProcess(CglPreProcess *process)
{
  process->postProcess(*solver_);
  delete solver_;
  solver_ = process->originalModel();
}

/* Process root node and return a strengthened model

The method assumes that initialSolve() has been called to solve the
LP relaxation. It processes the root node and then returns a pointer
to the strengthened model (or NULL if infeasible)
*/
OsiSolverInterface *
CbcModel::strengthenedModel()
{
  /*
      Switch off heuristics
    */
  int saveNumberHeuristics = numberHeuristics_;
  numberHeuristics_ = 0;
  /*
      Scan the variables, noting the integer variables. Create an
      CbcSimpleInteger object for each integer variable.
    */
  findIntegers(false);
  /*
      Ensure that objects on the lists of OsiObjects, heuristics, and cut
      generators attached to this model all refer to this model.
    */
  synchronizeModel();

  // Set so we can tell we are in initial phase in resolve
  continuousObjective_ = -COIN_DBL_MAX;
  /*
      Solve the relaxation.

      Apparently there are circumstances where this will be non-trivial --- i.e.,
      we've done something since initialSolve that's trashed the solution to the
      continuous relaxation.
    */
  bool feasible = resolve(NULL, 0) != 0;
  /*
      If the linear relaxation of the root is infeasible, bail out now. Otherwise,
      continue with processing the root node.
    */
  if (!feasible) {
    handler_->message(CBC_INFEAS, messages_) << CoinMessageEol;
    return NULL;
  }
  // Save objective (just so user can access it)
  originalContinuousObjective_ = solver_->getObjValue();

  /*
      Begin setup to process a feasible root node.
    */
  bestObjective_ = CoinMin(bestObjective_, 1.0e50);
  numberSolutions_ = 0;
  numberHeuristicSolutions_ = 0;
  // Everything is minimization
  double cutoff = getCutoff();
  double direction = solver_->getObjSense();
  if (cutoff < 1.0e20 && direction < 0.0)
    messageHandler()->message(CBC_CUTOFF_WARNING1,
      messages())
      << cutoff << -cutoff << CoinMessageEol;
  if (cutoff > bestObjective_)
    cutoff = bestObjective_;
  setCutoff(cutoff);
  /*
      We probably already have a current solution, but just in case ...
    */
  int numberColumns = getNumCols();
  if (!currentSolution_)
    currentSolution_ = new double[numberColumns];
  testSolution_ = currentSolution_;
  /*
      Create a copy of the solver, thus capturing the original (root node)
      constraint system (aka the continuous system).
    */
  continuousSolver_ = solver_->clone();
  numberRowsAtContinuous_ = getNumRows();
  /*
      Check the objective to see if we can deduce a nontrivial increment. If
      it's better than the current value for CbcCutoffIncrement, it'll be
      installed.
    */
  analyzeObjective();
  /*
      Set up for cut generation. addedCuts_ holds the cuts which are relevant for
      the active subproblem. whichGenerator will be used to record the generator
      that produced a given cut.
    */
  maximumWhich_ = INITIAL_MAXIMUM_WHICH;
  delete[] whichGenerator_;
  whichGenerator_ = new int[maximumWhich_];
  maximumNumberCuts_ = 0;
  currentNumberCuts_ = 0;
  delete[] addedCuts_;
  addedCuts_ = NULL;
  /*
    Generate cuts at the root node and reoptimise. solveWithCuts does the heavy
    lifting. It will iterate a generate/reoptimise loop (including reduced cost
    fixing) until no cuts are generated, the change in objective falls off,  or
    the limit on the number of rounds of cut generation is exceeded.

    At the end of all this, any cuts will be recorded in cuts and also
    installed in the solver's constraint system. We'll have reoptimised, and
    removed any slack cuts (numberOldActiveCuts_ and numberNewCuts_ have been
    adjusted accordingly).

    Tell cut generators they can be a bit more aggressive at root node

    */
  int iCutGenerator;
  for (iCutGenerator = 0; iCutGenerator < numberCutGenerators_; iCutGenerator++) {
    CglCutGenerator *generator = generator_[iCutGenerator]->generator();
    generator->setAggressiveness(generator->getAggressiveness() + 100);
  }
  OsiCuts cuts;
  numberOldActiveCuts_ = 0;
  numberNewCuts_ = 0;
  {
    int iObject;
    int numberUnsatisfied = 0;
    memcpy(currentSolution_, solver_->getColSolution(),
      numberColumns * sizeof(double));

    // point to useful information
    OsiBranchingInformation usefulInfo = usefulInformation();
    for (iObject = 0; iObject < numberObjects_; iObject++) {
      double infeasibility = object_[iObject]->checkInfeasibility(&usefulInfo);
      if (infeasibility)
        numberUnsatisfied++;
    }
    if (numberUnsatisfied) {
      feasible = solveWithCuts(cuts, maximumCutPassesAtRoot_,
        NULL);
    }
  }
  /*
      We've taken the continuous relaxation as far as we can.
    */

  OsiSolverInterface *newSolver = NULL;
  if (feasible) {
    // make copy of current solver
    newSolver = solver_->clone();
  }
  /*
      Clean up dangling objects. continuousSolver_ may already be toast.
    */
  delete[] whichGenerator_;
  whichGenerator_ = NULL;
  delete[] walkback_;
  walkback_ = NULL;
  delete[] lastNodeInfo_;
  lastNodeInfo_ = NULL;
  delete[] lastNumberCuts_;
  lastNumberCuts_ = NULL;
  delete[] lastCut_;
  lastCut_ = NULL;
  delete[] addedCuts_;
  addedCuts_ = NULL;
  if (continuousSolver_) {
    delete continuousSolver_;
    continuousSolver_ = NULL;
  }
  /*
      Destroy global cuts by replacing with an empty OsiCuts object.
    */
  globalCuts_ = OsiCuts();
  delete globalConflictCuts_;
  globalConflictCuts_ = NULL;
  numberHeuristics_ = saveNumberHeuristics;

  return newSolver;
}
/*  create a submodel from partially fixed problem

The method creates a new clean model with given bounds.
*/
CbcModel *
CbcModel::cleanModel(const double *lower, const double *upper)
{
  OsiSolverInterface *solver = continuousSolver_->clone();

  int numberIntegers = numberIntegers_;
  const int *integerVariable = integerVariable_;

  int i;
  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    const OsiObject *object = object_[i];
#ifndef NDEBUG
    const CbcSimpleInteger *integerObject = dynamic_cast< const CbcSimpleInteger * >(object);
    assert(integerObject);
#else
      const CbcSimpleInteger *integerObject = static_cast< const CbcSimpleInteger * >(object);
#endif
    // get original bounds
    double originalLower = integerObject->originalLowerBound();
    double originalUpper = integerObject->originalUpperBound();
    solver->setColLower(iColumn, CoinMax(lower[iColumn], originalLower));
    solver->setColUpper(iColumn, CoinMin(upper[iColumn], originalUpper));
  }
  CbcModel *model = new CbcModel(*solver);
  // off some messages
  if (handler_->logLevel() <= 1) {
    model->messagesPointer()->setDetailMessage(3, 9);
    model->messagesPointer()->setDetailMessage(3, 6);
    model->messagesPointer()->setDetailMessage(3, 4);
    model->messagesPointer()->setDetailMessage(3, 1);
    model->messagesPointer()->setDetailMessage(3, 13);
    model->messagesPointer()->setDetailMessage(3, 14);
    model->messagesPointer()->setDetailMessage(3, 3007);
  }
  // Cuts
  for (i = 0; i < numberCutGenerators_; i++) {
    int howOften = generator_[i]->howOftenInSub();
    if (howOften > -100) {
      CbcCutGenerator *generator = virginGenerator_[i];
      CglCutGenerator *cglGenerator = generator->generator();
      model->addCutGenerator(cglGenerator, howOften,
        generator->cutGeneratorName(),
        generator->normal(),
        generator->atSolution(),
        generator->whenInfeasible(),
        -100, generator->whatDepthInSub(), -1);
    }
  }
  double cutoff = getCutoff();
  model->setCutoff(cutoff);
  return model;
}
/* Invoke the branch & cut algorithm on partially fixed problem

   The method uses a subModel created by cleanModel. The search
   ends when the tree is exhausted or maximum nodes is reached.

   If better solution found then it is saved.

   Returns 0 if search completed and solution, 1 if not completed and solution,
   2 if completed and no solution, 3 if not completed and no solution.

   Normally okay to do subModel immediately followed by subBranchandBound
   (== other form of subBranchAndBound)
   but may need to get at model for advanced features.

   Deletes model

*/

int CbcModel::subBranchAndBound(CbcModel *model,
  CbcModel *presolvedModel,
  int maximumNodes)
{
  int i;
  double cutoff = model->getCutoff();
  CbcModel *model2;
  if (presolvedModel)
    model2 = presolvedModel;
  else
    model2 = model;
  // Do complete search

  for (i = 0; i < numberHeuristics_; i++) {
    model2->addHeuristic(heuristic_[i]);
    model2->heuristic(i)->resetModel(model2);
  }
  // Definition of node choice
  model2->setNodeComparison(nodeCompare_->clone());
  //model2->solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
  model2->messageHandler()->setLogLevel(CoinMax(0, handler_->logLevel() - 1));
  //model2->solver()->messageHandler()->setLogLevel(2);
  model2->setMaximumCutPassesAtRoot(maximumCutPassesAtRoot_);
  model2->setPrintFrequency(50);
  model2->setIntParam(CbcModel::CbcMaxNumNode, maximumNodes);
  model2->branchAndBound();
  delete model2->nodeComparison();
  if (model2->getMinimizationObjValue() > cutoff) {
    // no good
    if (model != model2)
      delete model2;
    delete model;
    return 2;
  }
  if (model != model2) {
    // get back solution
    model->originalModel(model2, false);
    delete model2;
  }
  int status;
  if (model->getMinimizationObjValue() < cutoff && model->bestSolution()) {
    double objValue = model->getObjValue();
    const double *solution = model->bestSolution();
    setBestSolution(CBC_TREE_SOL, objValue, solution);
    status = 0;
  } else {
    status = 2;
  }
  if (model->status())
    status++; // not finished search
  delete model;
  return status;
}
/* Invoke the branch & cut algorithm on partially fixed problem

The method creates a new model with given bounds, presolves it
then proceeds to explore the branch & cut search tree. The search
ends when the tree is exhausted or maximum nodes is reached.
Returns 0 if search completed and solution, 1 if not completed and solution,
2 if completed and no solution, 3 if not completed and no solution.
*/
int CbcModel::subBranchAndBound(const double *lower, const double *upper,
  int maximumNodes)
{
  OsiSolverInterface *solver = continuousSolver_->clone();

  int numberIntegers = numberIntegers_;
  const int *integerVariable = integerVariable_;

  int i;
  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    const OsiObject *object = object_[i];
#ifndef NDEBUG
    const CbcSimpleInteger *integerObject = dynamic_cast< const CbcSimpleInteger * >(object);
    assert(integerObject);
#else
      const CbcSimpleInteger *integerObject = static_cast< const CbcSimpleInteger * >(object);
#endif
    // get original bounds
    double originalLower = integerObject->originalLowerBound();
    double originalUpper = integerObject->originalUpperBound();
    solver->setColLower(iColumn, CoinMax(lower[iColumn], originalLower));
    solver->setColUpper(iColumn, CoinMin(upper[iColumn], originalUpper));
  }
  CbcModel model(*solver);
  // off some messages
  if (handler_->logLevel() <= 1) {
    model.messagesPointer()->setDetailMessage(3, 9);
    model.messagesPointer()->setDetailMessage(3, 6);
    model.messagesPointer()->setDetailMessage(3, 4);
    model.messagesPointer()->setDetailMessage(3, 1);
    model.messagesPointer()->setDetailMessage(3, 3007);
  }
  double cutoff = getCutoff();
  model.setCutoff(cutoff);
  // integer presolve
  CbcModel *model2 = model.integerPresolve(false);
  if (!model2 || !model2->getNumRows()) {
    delete model2;
    delete solver;
    return 2;
  }
  if (handler_->logLevel() > 1)
    printf("Reduced model has %d rows and %d columns\n",
      model2->getNumRows(), model2->getNumCols());
  // Do complete search

  // Cuts
  for (i = 0; i < numberCutGenerators_; i++) {
    int howOften = generator_[i]->howOftenInSub();
    if (howOften > -100) {
      CbcCutGenerator *generator = virginGenerator_[i];
      CglCutGenerator *cglGenerator = generator->generator();
      model2->addCutGenerator(cglGenerator, howOften,
        generator->cutGeneratorName(),
        generator->normal(),
        generator->atSolution(),
        generator->whenInfeasible(),
        -100, generator->whatDepthInSub(), -1);
    }
  }
  for (i = 0; i < numberHeuristics_; i++) {
    model2->addHeuristic(heuristic_[i]);
    model2->heuristic(i)->resetModel(model2);
  }
  // Definition of node choice
  model2->setNodeComparison(nodeCompare_->clone());
  //model2->solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
  model2->messageHandler()->setLogLevel(CoinMax(0, handler_->logLevel() - 1));
  //model2->solver()->messageHandler()->setLogLevel(2);
  model2->setMaximumCutPassesAtRoot(maximumCutPassesAtRoot_);
  model2->setPrintFrequency(50);
  model2->setIntParam(CbcModel::CbcMaxNumNode, maximumNodes);
  model2->branchAndBound();
  delete model2->nodeComparison();
  if (model2->getMinimizationObjValue() > cutoff) {
    // no good
    delete model2;
    delete solver;
    return 2;
  }
  // get back solution
  model.originalModel(model2, false);
  delete model2;
  int status;
  if (model.getMinimizationObjValue() < cutoff && model.bestSolution()) {
    double objValue = model.getObjValue();
    const double *solution = model.bestSolution();
    setBestSolution(CBC_TREE_SOL, objValue, solution);
    status = 0;
  } else {
    status = 2;
  }
  if (model.status())
    status++; // not finished search
  delete solver;
  return status;
}
#endif

static void *doRootCbcThread(void *voidInfo)
{
  CbcModel *model = reinterpret_cast< CbcModel * >(voidInfo);
#ifdef COIN_HAS_CLP
  OsiClpSolverInterface *clpSolver
    = dynamic_cast< OsiClpSolverInterface * >(model->solver());
  char general[200];
  if (clpSolver) {
    sprintf(general, "Starting multiple root solver");
    model->messageHandler()->message(CBC_GENERAL,
      model->messages())
      << general << CoinMessageEol;
    clpSolver->setHintParam(OsiDoReducePrint, true, OsiHintTry);
    ClpSimplex *simplex = clpSolver->getModelPtr();
    int logLevel = simplex->logLevel();
    if (logLevel <= 1)
      simplex->setLogLevel(0);
    simplex->dual();
    simplex->setLogLevel(logLevel);
    clpSolver->setWarmStart(NULL);
  } else {
    model->initialSolve();
    sprintf(general, "Solver did %d iterations in initialSolve\n",
      model->solver()->getIterationCount());
    model->messageHandler()->message(CBC_GENERAL,
      model->messages())
      << general << CoinMessageEol;
  }
#endif
  model->branchAndBound();
  sprintf(general, "Ending multiple root solver");
  model->messageHandler()->message(CBC_GENERAL,
    model->messages())
    << general << CoinMessageEol;
  return NULL;
}
OsiRowCut *
CbcModel::conflictCut(const OsiSolverInterface *solver, bool &localCuts)
{
  OsiRowCut *cut = NULL;
  localCuts = false;
#ifdef COIN_HAS_CLP
  const OsiClpSolverInterface *clpSolver
    = dynamic_cast< const OsiClpSolverInterface * >(solver);
  if (clpSolver && topOfTree_) {
    int debugMode = 0;
    const double *originalLower = topOfTree_->lower();
    const double *originalUpper = topOfTree_->upper();
    int typeCut = 1;
    ClpSimplex *simplex = clpSolver->getModelPtr();
    assert(simplex->status() == 1);
    if (simplex->ray()) {
      {
        int numberRows = simplex->numberRows();
        double *saveRay = CoinCopyOfArray(simplex->ray(), numberRows);
#define SAFE_RAY
#ifdef SAFE_RAY
        ClpSimplex &tempSimplex = *simplex;
#else
          ClpSimplex tempSimplex = *simplex;
#endif
        int logLevel = simplex->logLevel();
        tempSimplex.setLogLevel(63);
        tempSimplex.scaling(0);
        tempSimplex.dual();
        tempSimplex.setLogLevel(logLevel);
        if (!tempSimplex.numberIterations()) {
          double *ray = tempSimplex.ray();
          int nBad = 0;
          for (int i = 0; i < numberRows; i++) {
            if (fabs(ray[i] - saveRay[i]) > 1.0e-3) {
              if (debugMode)
                printf("row %d true %g bad %g - diff %g\n",
                  i, ray[i], saveRay[i], ray[i] - saveRay[i]);
              nBad++;
            }
          }
          if (nBad)
            printf("%d mismatch crunch ray values\n", nBad);
        }
        delete[] saveRay;
      }
      // make sure we use non-scaled versions
      ClpPackedMatrix *saveMatrix = simplex->swapScaledMatrix(NULL);
      double *saveScale = simplex->swapRowScale(NULL);
      //printf("Could do normal cut\n");
      // could use existing arrays
      int numberRows = simplex->numberRows();
      int numberColumns = simplex->numberColumns();
      double *farkas = new double[2 * numberColumns + numberRows];
      double *bound = farkas + numberColumns;
      double *effectiveRhs = bound + numberColumns;
      // sign as internally for dual - so swap if primal
      /*const*/ double *ray = simplex->ray();
      // have to get rid of local cut rows
      if (whichGenerator_) {
        const int *whichGenerator = whichGenerator_ - numberRowsAtContinuous_;
        int badRows = 0;
        for (int iRow = numberRowsAtContinuous_; iRow < numberRows; iRow++) {
          int iType = whichGenerator[iRow];
          if ((iType >= 0 && iType < 20000)) {
            if (fabs(ray[iRow]) > 1.0e-10) {
              badRows++;
            } else {
              ray[iRow] = 0.0;
            }
          }
        }
        if (badRows) {
          if ((debugMode & 1) != 0)
            printf("%d rows from local cuts\n", badRows);
          localCuts = true;
        }
      }
      // get farkas row
      memset(farkas, 0, (2 * numberColumns + numberRows) * sizeof(double));
      simplex->transposeTimes(-1.0, ray, farkas);
      //const char * integerInformation = simplex->integerType_;
      //assert (integerInformation);

      int sequenceOut = simplex->sequenceOut();
      // Put nonzero bounds in bound
      const double *columnLower = simplex->columnLower();
      const double *columnUpper = simplex->columnUpper();
      int numberBad = 0;
      for (int i = 0; i < numberColumns; i++) {
        double value = farkas[i];
        double boundValue = 0.0;
        if (simplex->getStatus(i) == ClpSimplex::basic) {
          // treat as zero if small
          if (fabs(value) < 1.0e-8) {
            value = 0.0;
            farkas[i] = 0.0;
          }
          if (value) {
            //printf("basic %d direction %d farkas %g\n",
            //	   i,simplex->directionOut(),value);
            if (value < 0.0)
              boundValue = columnLower[i];
            else
              boundValue = columnUpper[i];
          }
        } else if (fabs(value) > 1.0e-10) {
          if (value < 0.0)
            boundValue = columnLower[i];
          else
            boundValue = columnUpper[i];
        }
        bound[i] = boundValue;
        if (fabs(boundValue) > 1.0e10)
          numberBad++;
      }
      const double *rowLower = simplex->rowLower();
      const double *rowUpper = simplex->rowUpper();
      //int pivotRow = simplex->spareIntArray_[3];
      //bool badPivot=pivotRow<0;
      for (int i = 0; i < numberRows; i++) {
        double value = ray[i];
        double rhsValue = 0.0;
        if (simplex->getRowStatus(i) == ClpSimplex::basic) {
          // treat as zero if small
          if (fabs(value) < 1.0e-8) {
            value = 0.0;
            ray[i] = 0.0;
          }
          if (value) {
            //printf("row basic %d direction %d ray %g\n",
            //	   i,simplex->directionOut(),value);
            if (value < 0.0)
              rhsValue = rowLower[i];
            else
              rhsValue = rowUpper[i];
          }
        } else if (fabs(value) > 1.0e-10) {
          if (value < 0.0)
            rhsValue = rowLower[i];
          else
            rhsValue = rowUpper[i];
        }
        effectiveRhs[i] = rhsValue;
      }
      simplex->times(-1.0, bound, effectiveRhs);
      simplex->swapRowScale(saveScale);
      simplex->swapScaledMatrix(saveMatrix);
      double bSum = 0.0;
      for (int i = 0; i < numberRows; i++) {
        bSum += effectiveRhs[i] * ray[i];
      }
      if (numberBad || bSum > -1.0e-4) {
#ifndef NDEBUG
        printf("bad BOUND bSum %g  - %d bad\n",
          bSum, numberBad);
#endif
      } else {
        const char *integerInformation = simplex->integerInformation();
        assert(integerInformation);
        int *conflict = new int[numberColumns];
        double *sort = new double[numberColumns];
        double relax = 0.0;
        int nConflict = 0;
        int nOriginal = 0;
        int nFixed = 0;
        for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
          if (integerInformation[iColumn]) {
            if ((debugMode & 1) != 0)
              printf("%d status %d %g <= %g <=%g (orig %g, %g) farkas %g\n",
                iColumn, simplex->getStatus(iColumn), columnLower[iColumn],
                simplex->primalColumnSolution()[iColumn], columnUpper[iColumn],
                originalLower[iColumn], originalUpper[iColumn],
                farkas[iColumn]);
            double gap = originalUpper[iColumn] - originalLower[iColumn];
            if (!gap)
              continue;
            if (gap == columnUpper[iColumn] - columnLower[iColumn])
              nOriginal++;
            if (columnUpper[iColumn] == columnLower[iColumn])
              nFixed++;
            if (fabs(farkas[iColumn]) < 1.0e-15) {
              farkas[iColumn] = 0.0;
              continue;
            }
            // temp
            if (gap >= 20000.0 && false) {
              // can't use
              if (farkas[iColumn] < 0.0) {
                assert(originalLower[iColumn] - columnLower[iColumn] <= 0.0);
                // farkas is negative - relax lower bound all way
                relax += farkas[iColumn] * (originalLower[iColumn] - columnLower[iColumn]);
              } else {
                assert(originalUpper[iColumn] - columnUpper[iColumn] >= 0.0);
                // farkas is positive - relax upper bound all way
                relax += farkas[iColumn] * (originalUpper[iColumn] - columnUpper[iColumn]);
              }
              continue;
            }
            if (originalLower[iColumn] == columnLower[iColumn]) {
              if (farkas[iColumn] > 0.0 && (simplex->getStatus(iColumn) == ClpSimplex::atUpperBound || simplex->getStatus(iColumn) == ClpSimplex::isFixed || iColumn == sequenceOut)) {
                // farkas is positive - add to list
                gap = originalUpper[iColumn] - columnUpper[iColumn];
                if (gap) {
                  sort[nConflict] = -farkas[iColumn] * gap;
                  conflict[nConflict++] = iColumn;
                }
                //assert (gap>columnUpper[iColumn]-columnLower[iColumn]);
              }
            } else if (originalUpper[iColumn] == columnUpper[iColumn]) {
              if (farkas[iColumn] < 0.0 && (simplex->getStatus(iColumn) == ClpSimplex::atLowerBound || simplex->getStatus(iColumn) == ClpSimplex::isFixed || iColumn == sequenceOut)) {
                // farkas is negative - add to list
                gap = columnLower[iColumn] - originalLower[iColumn];
                if (gap) {
                  sort[nConflict] = farkas[iColumn] * gap;
                  conflict[nConflict++] = iColumn;
                }
                //assert (gap>columnUpper[iColumn]-columnLower[iColumn]);
              }
            } else {
              // can't use
              if (farkas[iColumn] < 0.0) {
                assert(originalLower[iColumn] - columnLower[iColumn] <= 0.0);
                // farkas is negative - relax lower bound all way
                relax += farkas[iColumn] * (originalLower[iColumn] - columnLower[iColumn]);
              } else {
                assert(originalUpper[iColumn] - columnUpper[iColumn] >= 0.0);
                // farkas is positive - relax upper bound all way
                relax += farkas[iColumn] * (originalUpper[iColumn] - columnUpper[iColumn]);
              }
            }
            assert(relax >= 0.0);
          } else {
            // not integer - but may have been got at
            double gap = originalUpper[iColumn] - originalLower[iColumn];
            if (gap > columnUpper[iColumn] - columnLower[iColumn]) {
              // can't use
              if (farkas[iColumn] < 0.0) {
                assert(originalLower[iColumn] - columnLower[iColumn] <= 0.0);
                // farkas is negative - relax lower bound all way
                relax += farkas[iColumn] * (originalLower[iColumn] - columnLower[iColumn]);
              } else {
                assert(originalUpper[iColumn] - columnUpper[iColumn] >= 0.0);
                // farkas is positive - relax upper bound all way
                relax += farkas[iColumn] * (originalUpper[iColumn] - columnUpper[iColumn]);
              }
            }
          }
        }
        if (relax + bSum > -1.0e-4 || !nConflict) {
          if (relax + bSum > -1.0e-4) {
#ifndef NDEBUG
            printf("General integers relax bSum to %g\n", relax + bSum);
#endif
          } else {
            printf("All variables relaxed and still infeasible - what does this mean?\n");
            int nR = 0;
            for (int i = 0; i < numberRows; i++) {
              if (fabs(ray[i]) > 1.0e-10)
                nR++;
              else
                ray[i] = 0.0;
            }
            int nC = 0;
            for (int i = 0; i < numberColumns; i++) {
              if (fabs(farkas[i]) > 1.0e-10)
                nC++;
              else
                farkas[i] = 0.0;
            }
            if (nR < 3 && nC < 5) {
              printf("BAD %d nonzero rows, %d nonzero columns\n", nR, nC);
            }
          }
        } else {
          printf("BOUNDS violation bSum %g (relaxed %g) - %d at original bounds, %d fixed - %d in conflict\n", bSum,
            relax + bSum, nOriginal, nFixed, nConflict);
          CoinSort_2(sort, sort + nConflict, conflict);
          int nC = nConflict;
          bSum += relax;
          double saveBsum = bSum;
          while (nConflict) {
            //int iColumn=conflict[nConflict-1];
            double change = -sort[nConflict - 1];
            if (bSum + change > -1.0e-4)
              break;
            nConflict--;
            bSum += change;
          }
          if (!nConflict) {
            int nR = 0;
            for (int i = 0; i < numberRows; i++) {
              if (fabs(ray[i]) > 1.0e-10)
                nR++;
              else
                ray[i] = 0.0;
            }
            int nC = 0;
            for (int i = 0; i < numberColumns; i++) {
              if (fabs(farkas[i]) > 1.0e-10)
                nC++;
              else
                farkas[i] = 0.0;
            }
            if (nR < 3 && nC < 5) {
              printf("BAD2 %d nonzero rows, %d nonzero columns\n", nR, nC);
            }
          }
          // no point doing if no reduction (or big?) ?
          if (nConflict < nC + 1 && nConflict < 500) {
            cut = new OsiRowCut();
            cut->setUb(COIN_DBL_MAX);
            if (!typeCut) {
              double lo = 1.0;
              for (int i = 0; i < nConflict; i++) {
                int iColumn = conflict[i];
                if (originalLower[iColumn] == columnLower[iColumn]) {
                  // must be at least one higher
                  sort[i] = 1.0;
                  lo += originalLower[iColumn];
                } else {
                  // must be at least one lower
                  sort[i] = -1.0;
                  lo -= originalUpper[iColumn];
                }
              }
              cut->setLb(lo);
              cut->setRow(nConflict, conflict, sort);
              printf("CUT has %d (started at %d) - final bSum %g\n", nConflict, nC, bSum);
            } else {
              // just save for use later
              // first take off small
              int nC2 = nC;
              while (nC2) {
                //int iColumn=conflict[nConflict-1];
                double change = -sort[nC2 - 1];
                if (saveBsum + change > -1.0e-4 || change > 1.0e-4)
                  break;
                nC2--;
                saveBsum += change;
              }
              cut->setLb(saveBsum);
              for (int i = 0; i < nC2; i++) {
                int iColumn = conflict[i];
                sort[i] = farkas[iColumn];
              }
              cut->setRow(nC2, conflict, sort);
              printf("Stem CUT has %d (greedy %d - with small %d) - saved bSum %g final greedy bSum %g\n",
                nC2, nConflict, nC, saveBsum, bSum);
            }
          }
        }
        delete[] conflict;
        delete[] sort;
      }
      delete[] farkas;
    } else {
      printf("No dual ray\n");
    }
  }
#endif
  return cut;
}

void CbcModel::setMIPStart(int count, const char **colNames, const double colValues[])
{
  mipStart_.clear();
  for (int i = 0; (i < count); ++i)
    mipStart_.push_back(std::pair< std::string, double >(std::string(colNames[i]), colValues[i]));
}
#ifdef COIN_HAS_NTY
// get rid of all
void CbcModel::zapSymmetry()
{
  delete symmetryInfo_;
  symmetryInfo_ = NULL;
}
#endif
/* Add SOS info to solver -
   Overwrites SOS information in solver with information
   in CbcModel.  Has no effect with some solvers. 
   Also updates integer info. */
void CbcModel::addSOSEtcToSolver()
{
  // at present just for OsiClp
#ifdef COIN_HAS_CLP
  OsiClpSolverInterface *clpSolver
    = dynamic_cast< OsiClpSolverInterface * >(solver_);
  if (clpSolver) {
    int numberColumns = clpSolver->getNumCols();
    for (int i = 0; i < numberColumns; i++)
      clpSolver->setContinuous(i);
    int nOdd = 0;
    int numberSOS = 0;
    for (int i = 0; i < numberObjects_; i++) {
      CbcObject *obj = dynamic_cast< CbcObject * >(object_[i]);
      CbcSimpleInteger *thisInt = dynamic_cast< CbcSimpleInteger * >(obj);
      OsiSOS *objSOS1 = dynamic_cast< OsiSOS * >(obj);
      CbcSOS *objSOS2 = dynamic_cast< CbcSOS * >(obj);
      if (thisInt) {
        clpSolver->setInteger(thisInt->columnNumber());
      } else if (objSOS1) {
        numberSOS++;
      } else if (objSOS2) {
        numberSOS++;
      } else {
        nOdd++;
      }
    }
    if (nOdd) {
      char general[200];
      sprintf(general, "%d objects not SOS or Integer - can't move to Osi",
        nOdd);
      messageHandler()->message(CBC_GENERAL,
        messages())
        << general << CoinMessageEol;
    }
    if (numberSOS) {
      CoinSet *setInfo = new CoinSet[numberSOS];
      numberSOS = 0;
      for (int i = 0; i < numberObjects_; i++) {
        CbcObject *obj = dynamic_cast< CbcObject * >(object_[i]);
        OsiSOS *objSOS1 = dynamic_cast< OsiSOS * >(obj);
        CbcSOS *objSOS2 = dynamic_cast< CbcSOS * >(obj);
        if (objSOS1 || objSOS2) {
          int numberMembers;
          const int *members;
          int type;
          const double *weights;
          if (objSOS1) {
            numberMembers = objSOS1->numberMembers();
            members = objSOS1->members();
            type = objSOS1->sosType();
            weights = objSOS1->weights();
          } else {
            numberMembers = objSOS2->numberMembers();
            members = objSOS2->members();
            type = objSOS2->sosType();
            weights = objSOS2->weights();
          }
          CoinSosSet info(numberMembers, members,
            weights, type);
          //info.setSetType(type);
          //memcpy(info.modifiableWeights(),weights,
          //	 numberMembers*sizeof(double));
          setInfo[numberSOS++] = info;
        }
      }
      clpSolver->replaceSetInfo(numberSOS, setInfo);
    }
  }
#endif
}
/*
  Clean model i.e. make SOS/integer variables exactly at bound if needed
  Only if moreSpecialOptions2_ 32768
  Returns number of variables forced out
  cleanVariables array will be used if exists
*/
int CbcModel::cleanBounds(OsiSolverInterface *solver, char *cleanIn)
{
  int numberBad = 0;
#ifdef COIN_HAS_CLP
#ifndef ZERO_ODD_TOLERANCE
#define ZERO_ODD_TOLERANCE 1.0e-14
#endif
  OsiClpSolverInterface *osiclp = dynamic_cast< OsiClpSolverInterface * >(solver);
  if (osiclp && osiclp->isProvenOptimal()) {
    int numberColumns = osiclp->getNumCols();
    char *cleanVariables;
    if (!cleanIn) {
      cleanVariables = setupCleanVariables();
    } else {
      cleanVariables = cleanIn;
    }
    // If any small values re-do
    ClpSimplex *clp = osiclp->getModelPtr();
    double *solution = clp->primalColumnSolution();
    const double *columnLower = clp->columnLower();
    const double *columnUpper = clp->columnUpper();
    //#define SHOW_BAD
#ifdef SHOW_BAD
    double sumBadLow = 0.0;
    double sumBadHigh = 0.0;
    double maxBadLow = 0.0;
    double maxBadHigh = 0.0;
#endif
    for (int i = 0; i < numberColumns; i++) {
      if (cleanVariables[i]) {
        if (solution[i] > columnUpper[i] + ZERO_ODD_TOLERANCE) {
          numberBad++;
#ifdef SHOW_BAD
          sumBadHigh += solution[i] - columnUpper[i];
          maxBadHigh = CoinMax(maxBadHigh, solution[i] - columnUpper[i]);
#endif
        } else if (solution[i] < columnLower[i] - ZERO_ODD_TOLERANCE) {
          numberBad++;
#ifdef SHOW_BAD
          sumBadLow += columnLower[i] - solution[i];
          maxBadLow = CoinMax(maxBadLow, columnLower[i] - solution[i]);
#endif
        }
      }
    }
    if (numberBad) {
#ifdef SHOW_BAD
      printf("%d small variables low (%g max %g), high (%g max %g)\n", numberBad, sumBadLow, maxBadLow, sumBadHigh, maxBadHigh);
#endif
      for (int i = 0; i < numberColumns; i++) {
        if (cleanVariables[i]) {
          if (solution[i] > columnUpper[i] + ZERO_ODD_TOLERANCE) {
            solution[i] = columnUpper[i];
            clp->setColumnStatus(i, ClpSimplex::atUpperBound);
          } else if (solution[i] < columnLower[i] - ZERO_ODD_TOLERANCE) {
            solution[i] = columnLower[i];
            clp->setColumnStatus(i, ClpSimplex::atLowerBound);
          }
        }
      }
      int saveLevel = clp->logLevel();
      clp->setLogLevel(0);
      clp->dual();
      clp->setLogLevel(saveLevel);
    }
    if (!cleanIn)
      delete[] cleanVariables;
  }
#endif
  return numberBad;
}
// Sets up cleanVariables array (i.e. ones to be careful about)
char *
CbcModel::setupCleanVariables()
{
  OsiClpSolverInterface *osiclp = dynamic_cast< OsiClpSolverInterface * >(solver_);
  int numberColumns = osiclp->getNumCols();
  char *cleanVariables = NULL;
  if (osiclp) {
    cleanVariables = new char[numberColumns];
    memset(cleanVariables, 0, numberColumns);
    for (int i = 0; i < numberObjects_; i++) {
      const CbcSimpleInteger *intvar = dynamic_cast< const CbcSimpleInteger * >(object_[i]);
      const CbcSOS *sos = dynamic_cast< const CbcSOS * >(object_[i]);
      if (intvar) {
#ifdef CLEAN_INTEGER_VARIABLES
        cleanVariables[intvar->columnNumber()] = 1;
#endif
      } else if (sos) {
        int n = sos->numberMembers();
        const int *members = sos->members();
        for (int j = 0; j < n; j++)
          cleanVariables[members[j]] = 2;
      }
    }
  }
  return cleanVariables;
}
/* Returns postProcessed solution in solver(called from event handler)
   Normally used for integer solution (not really tested otherwise)
   solutionType 1 is best integer so far, 0 is current solution 
   (may not be integer) */
const OsiSolverInterface *
CbcModel::postProcessedSolver(int solutionType)
{
  CbcModel *model = this;
  CglPreProcess *processPointer = model->preProcess();
  OsiSolverInterface *originalModel = NULL;
  const double *solution = bestSolution();
  while (processPointer) {
    int numberSolvers = processPointer->numberSolvers();
    OsiSolverInterface *solver = processPointer->presolve(numberSolvers - 1)->presolvedModel();
    if (solutionType) {
      // massage solution
      int numberColumns = solver->getNumCols();
      double *saveLower = CoinCopyOfArray(model->solver()->getColLower(), numberColumns);
      double *saveUpper = CoinCopyOfArray(model->solver()->getColUpper(), numberColumns);
      const double *saveSolution = testSolution_;
      setTestSolution(solution);
      model->solver()->setColLower(solution);
      model->solver()->setColUpper(solution);
      OsiBranchingInformation usefulInfo = model->usefulInformation();
      /*
	Run through the objects and use feasibleRegion() to set variable bounds
	so as to fix the variables specified in the objects at their value in this
	solution. Since the object list contains (at least) one object for every
	integer variable, this has the effect of fixing all integer variables.
      */
      for (int i = 0; i < model->numberObjects(); i++)
        model->object(i)->feasibleRegion(solver, &usefulInfo);
      setTestSolution(saveSolution);
      model->solver()->setColLower(saveLower);
      model->solver()->setColUpper(saveUpper);
      delete[] saveLower;
      delete[] saveUpper;
    }
    solver->resolve();
    processPointer->postProcess(*solver, false);
    originalModel = processPointer->originalModel();
    solution = originalModel->getColSolution();
    processPointer = NULL;
    model = model->parentModel();
    processPointer = model ? model->preProcess() : NULL;
  }
  return originalModel;
}
// Delete a node and possibly null out currentNode_
void
CbcModel::deleteNode(CbcNode * node)
{
  delete node;
  if (node==currentNode_)
    currentNode_ = NULL;
}

