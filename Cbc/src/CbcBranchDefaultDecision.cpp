// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/10/2009-- carved out of CbcBranchActual

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>
//#define CBC_DEBUG

#include "CoinTypes.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiSolverBranch.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcBranchDefaultDecision.hpp"
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

//##############################################################################

// Default Constructor
CbcBranchDefaultDecision::CbcBranchDefaultDecision()
        : CbcBranchDecision()
{
    bestCriterion_ = 0.0;
    bestChangeUp_ = 0.0;
    bestNumberUp_ = 0;
    bestChangeDown_ = 0.0;
    bestObject_ = NULL;
    bestNumberDown_ = 0;
}

// Copy constructor
CbcBranchDefaultDecision::CbcBranchDefaultDecision (
    const CbcBranchDefaultDecision & rhs)
        : CbcBranchDecision(rhs)
{
    bestCriterion_ = rhs.bestCriterion_;
    bestChangeUp_ = rhs.bestChangeUp_;
    bestNumberUp_ = rhs.bestNumberUp_;
    bestChangeDown_ = rhs.bestChangeDown_;
    bestNumberDown_ = rhs.bestNumberDown_;
    bestObject_ = rhs.bestObject_;
    model_ = rhs.model_;
}

CbcBranchDefaultDecision::~CbcBranchDefaultDecision()
{
}

// Clone
CbcBranchDecision *
CbcBranchDefaultDecision::clone() const
{
    return new CbcBranchDefaultDecision(*this);
}

// Initialize i.e. before start of choosing at a node
void
CbcBranchDefaultDecision::initialize(CbcModel * model)
{
    bestCriterion_ = 0.0;
    bestChangeUp_ = 0.0;
    bestNumberUp_ = 0;
    bestChangeDown_ = 0.0;
    bestNumberDown_ = 0;
    bestObject_ = NULL;
    model_ = model;
}


/*
  Simple default decision algorithm. Compare based on infeasibility (numInfUp,
  numInfDn) until a solution is found by search, then switch to change in
  objective (changeUp, changeDn). Note that bestSoFar is remembered in
  bestObject_, so the parameter bestSoFar is unused.
*/

int
CbcBranchDefaultDecision::betterBranch(CbcBranchingObject * thisOne,
                                       CbcBranchingObject * /*bestSoFar*/,
                                       double changeUp, int numInfUp,
                                       double changeDn, int numInfDn)
{
    bool beforeSolution = cbcModel()->getSolutionCount() ==
                          cbcModel()->getNumberHeuristicSolutions();;
    int betterWay = 0;
    if (beforeSolution) {
        if (!bestObject_) {
            bestNumberUp_ = COIN_INT_MAX;
            bestNumberDown_ = COIN_INT_MAX;
        }
        // before solution - choose smallest number
        // could add in depth as well
        int bestNumber = CoinMin(bestNumberUp_, bestNumberDown_);
        if (numInfUp < numInfDn) {
            if (numInfUp < bestNumber) {
                betterWay = 1;
            } else if (numInfUp == bestNumber) {
                if (changeUp < bestCriterion_)
                    betterWay = 1;
            }
        } else if (numInfUp > numInfDn) {
            if (numInfDn < bestNumber) {
                betterWay = -1;
            } else if (numInfDn == bestNumber) {
                if (changeDn < bestCriterion_)
                    betterWay = -1;
            }
        } else {
            // up and down have same number
            bool better = false;
            if (numInfUp < bestNumber) {
                better = true;
            } else if (numInfUp == bestNumber) {
                if (CoinMin(changeUp, changeDn) < bestCriterion_)
                    better = true;;
            }
            if (better) {
                // see which way
                if (changeUp <= changeDn)
                    betterWay = 1;
                else
                    betterWay = -1;
            }
        }
    } else {
        if (!bestObject_) {
            bestCriterion_ = -1.0;
        }
        // got a solution
        if (changeUp <= changeDn) {
            if (changeUp > bestCriterion_)
                betterWay = 1;
        } else {
            if (changeDn > bestCriterion_)
                betterWay = -1;
        }
    }
    if (betterWay) {
        bestCriterion_ = CoinMin(changeUp, changeDn);
        bestChangeUp_ = changeUp;
        bestNumberUp_ = numInfUp;
        bestChangeDown_ = changeDn;
        bestNumberDown_ = numInfDn;
        bestObject_ = thisOne;
        // See if user is overriding way
        if (thisOne->object() && thisOne->object()->preferredWay())
            betterWay = thisOne->object()->preferredWay();
    }
    return betterWay;
}
/* Sets or gets best criterion so far */
void
CbcBranchDefaultDecision::setBestCriterion(double value)
{
    bestCriterion_ = value;
}
double
CbcBranchDefaultDecision::getBestCriterion() const
{
    return bestCriterion_;
}

/* Compare N branching objects. Return index of best
   and sets way of branching in chosen object.

   This routine is used only after strong branching.
*/

int
CbcBranchDefaultDecision::bestBranch (CbcBranchingObject ** objects, int numberObjects,
                                      int numberUnsatisfied,
                                      double * changeUp, int * numberInfeasibilitiesUp,
                                      double * changeDown, int * numberInfeasibilitiesDown,
                                      double objectiveValue)
{

    int bestWay = 0;
    int whichObject = -1;
    if (numberObjects) {
        CbcModel * model = cbcModel();
        // at continuous
        //double continuousObjective = model->getContinuousObjective();
        //int continuousInfeasibilities = model->getContinuousInfeasibilities();

        // average cost to get rid of infeasibility
        //double averageCostPerInfeasibility =
        //(objectiveValue-continuousObjective)/
        //(double) (abs(continuousInfeasibilities-numberUnsatisfied)+1);
        /* beforeSolution is :
           0 - before any solution
           n - n heuristic solutions but no branched one
           -1 - branched solution found
        */
        int numberSolutions = model->getSolutionCount();
        double cutoff = model->getCutoff();
        int method = 0;
        int i;
        if (numberSolutions) {
            int numberHeuristic = model->getNumberHeuristicSolutions();
            if (numberHeuristic < numberSolutions) {
                method = 1;
            } else {
                method = 2;
                // look further
                for ( i = 0 ; i < numberObjects ; i++) {
                    int numberNext = numberInfeasibilitiesUp[i];

                    if (numberNext < numberUnsatisfied) {
                        int numberUp = numberUnsatisfied - numberInfeasibilitiesUp[i];
                        double perUnsatisfied = changeUp[i] / static_cast<double> (numberUp);
                        double estimatedObjective = objectiveValue + numberUnsatisfied * perUnsatisfied;
                        if (estimatedObjective < cutoff)
                            method = 3;
                    }
                    numberNext = numberInfeasibilitiesDown[i];
                    if (numberNext < numberUnsatisfied) {
                        int numberDown = numberUnsatisfied - numberInfeasibilitiesDown[i];
                        double perUnsatisfied = changeDown[i] / static_cast<double> (numberDown);
                        double estimatedObjective = objectiveValue + numberUnsatisfied * perUnsatisfied;
                        if (estimatedObjective < cutoff)
                            method = 3;
                    }
                }
            }
            method = 2;
        } else {
            method = 0;
        }
        // Uncomment next to force method 4
        //method=4;

	// FIXME This should be an enum.  It will be easier to
	// understand in the code than numbers.
        /* Methods :
           0 - fewest infeasibilities
           1 - largest min change in objective
           2 - as 1 but use sum of changes if min close
           3 - predicted best solution
           4 - take cheapest up branch if infeasibilities same
        */
        int bestNumber = COIN_INT_MAX;
        double bestCriterion = -1.0e50;
        double alternativeCriterion = -1.0;
        double bestEstimate = 1.0e100;
        switch (method) {
        case 0:
            // could add in depth as well
            for ( i = 0 ; i < numberObjects ; i++) {
                int thisNumber = CoinMin(numberInfeasibilitiesUp[i], numberInfeasibilitiesDown[i]);
                if (thisNumber <= bestNumber) {
                    int betterWay = 0;
                    if (numberInfeasibilitiesUp[i] < numberInfeasibilitiesDown[i]) {
                        if (numberInfeasibilitiesUp[i] < bestNumber) {
                            betterWay = 1;
                        } else {
                            if (changeUp[i] < bestCriterion)
                                betterWay = 1;
                        }
                    } else if (numberInfeasibilitiesUp[i] > numberInfeasibilitiesDown[i]) {
                        if (numberInfeasibilitiesDown[i] < bestNumber) {
                            betterWay = -1;
                        } else {
                            if (changeDown[i] < bestCriterion)
                                betterWay = -1;
                        }
                    } else {
                        // up and down have same number
                        bool better = false;
                        if (numberInfeasibilitiesUp[i] < bestNumber) {
                            better = true;
                        } else if (numberInfeasibilitiesUp[i] == bestNumber) {
                            if (CoinMin(changeUp[i], changeDown[i]) < bestCriterion)
                                better = true;;
                        }
                        if (better) {
                            // see which way
                            if (changeUp[i] <= changeDown[i])
                                betterWay = 1;
                            else
                                betterWay = -1;
                        }
                    }
                    if (betterWay) {
                        bestCriterion = CoinMin(changeUp[i], changeDown[i]);
                        bestNumber = thisNumber;
                        whichObject = i;
                        bestWay = betterWay;
                    }
                }
            }
            break;
        case 1:
            for ( i = 0 ; i < numberObjects ; i++) {
                int betterWay = 0;
                if (changeUp[i] <= changeDown[i]) {
                    if (changeUp[i] > bestCriterion)
                        betterWay = 1;
                } else {
                    if (changeDown[i] > bestCriterion)
                        betterWay = -1;
                }
                if (betterWay) {
                    bestCriterion = CoinMin(changeUp[i], changeDown[i]);
                    whichObject = i;
                    bestWay = betterWay;
                }
            }
            break;
        case 2:
            for ( i = 0 ; i < numberObjects ; i++) {
                double change = CoinMin(changeUp[i], changeDown[i]);
                double sum = changeUp[i] + changeDown[i];
                bool take = false;
                if (change > 1.1*bestCriterion)
                    take = true;
                else if (change > 0.9*bestCriterion && sum + change > bestCriterion + alternativeCriterion)
                    take = true;
                if (take) {
                    if (changeUp[i] <= changeDown[i]) {
                        if (changeUp[i] > bestCriterion)
                            bestWay = 1;
                    } else {
                        if (changeDown[i] > bestCriterion)
                            bestWay = -1;
                    }
                    bestCriterion = change;
                    alternativeCriterion = sum;
                    whichObject = i;
                }
            }
            break;
        case 3:
            for ( i = 0 ; i < numberObjects ; i++) {
                int numberNext = numberInfeasibilitiesUp[i];

                if (numberNext < numberUnsatisfied) {
                    int numberUp = numberUnsatisfied - numberInfeasibilitiesUp[i];
                    double perUnsatisfied = changeUp[i] / static_cast<double> (numberUp);
                    double estimatedObjective = objectiveValue + numberUnsatisfied * perUnsatisfied;
                    if (estimatedObjective < bestEstimate) {
                        bestEstimate = estimatedObjective;
                        bestWay = 1;
                        whichObject = i;
                    }
                }
                numberNext = numberInfeasibilitiesDown[i];
                if (numberNext < numberUnsatisfied) {
                    int numberDown = numberUnsatisfied - numberInfeasibilitiesDown[i];
                    double perUnsatisfied = changeDown[i] / static_cast<double> (numberDown);
                    double estimatedObjective = objectiveValue + numberUnsatisfied * perUnsatisfied;
                    if (estimatedObjective < bestEstimate) {
                        bestEstimate = estimatedObjective;
                        bestWay = -1;
                        whichObject = i;
                    }
                }
            }
            break;
        case 4:
            // if number infeas same then cheapest up
            // first get best number or when going down
            // now choose smallest change up amongst equal number infeas
            for ( i = 0 ; i < numberObjects ; i++) {
                int thisNumber = CoinMin(numberInfeasibilitiesUp[i], numberInfeasibilitiesDown[i]);
                if (thisNumber <= bestNumber) {
                    int betterWay = 0;
                    if (numberInfeasibilitiesUp[i] < numberInfeasibilitiesDown[i]) {
                        if (numberInfeasibilitiesUp[i] < bestNumber) {
                            betterWay = 1;
                        } else {
                            if (changeUp[i] < bestCriterion)
                                betterWay = 1;
                        }
                    } else if (numberInfeasibilitiesUp[i] > numberInfeasibilitiesDown[i]) {
                        if (numberInfeasibilitiesDown[i] < bestNumber) {
                            betterWay = -1;
                        } else {
                            if (changeDown[i] < bestCriterion)
                                betterWay = -1;
                        }
                    } else {
                        // up and down have same number
                        bool better = false;
                        if (numberInfeasibilitiesUp[i] < bestNumber) {
                            better = true;
                        } else if (numberInfeasibilitiesUp[i] == bestNumber) {
                            if (CoinMin(changeUp[i], changeDown[i]) < bestCriterion)
                                better = true;;
                        }
                        if (better) {
                            // see which way
                            if (changeUp[i] <= changeDown[i])
                                betterWay = 1;
                            else
                                betterWay = -1;
                        }
                    }
                    if (betterWay) {
                        bestCriterion = CoinMin(changeUp[i], changeDown[i]);
                        bestNumber = thisNumber;
                        whichObject = i;
                        bestWay = betterWay;
                    }
                }
            }
            bestCriterion = 1.0e50;
            for ( i = 0 ; i < numberObjects ; i++) {
                int thisNumber = numberInfeasibilitiesUp[i];
                if (thisNumber == bestNumber && changeUp) {
                    if (changeUp[i] < bestCriterion) {
                        bestCriterion = changeUp[i];
                        whichObject = i;
                        bestWay = 1;
                    }
                }
            }
            break;
        }
        // set way in best
        if (whichObject >= 0) {
            CbcBranchingObject * bestObject = objects[whichObject];
            if (bestObject->object() && bestObject->object()->preferredWay())
                bestWay = bestObject->object()->preferredWay();
            bestObject->way(bestWay);
        } else {
	  COIN_DETAIL_PRINT(printf("debug\n"));
        }
    }
    return whichObject;
}

