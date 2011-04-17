/* $Id: CbcHeuristicPivotAndFix.cpp 1200 2009-07-25 08:44:13Z forrest $ */
// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <vector>

#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcHeuristicPivotAndFix.hpp"
#include "OsiClpSolverInterface.hpp"
#include  "CoinTime.hpp"

//#define FORNOW

// Default Constructor
CbcHeuristicPivotAndFix::CbcHeuristicPivotAndFix()
        : CbcHeuristic()
{
}

// Constructor with model - assumed before cuts

CbcHeuristicPivotAndFix::CbcHeuristicPivotAndFix(CbcModel & model)
        : CbcHeuristic(model)
{
}

// Destructor
CbcHeuristicPivotAndFix::~CbcHeuristicPivotAndFix ()
{
}

// Clone
CbcHeuristic *
CbcHeuristicPivotAndFix::clone() const
{
    return new CbcHeuristicPivotAndFix(*this);
}
// Create C++ lines to get to current state
void
CbcHeuristicPivotAndFix::generateCpp( FILE * fp)
{
    CbcHeuristicPivotAndFix other;
    fprintf(fp, "0#include \"CbcHeuristicPivotAndFix.hpp\"\n");
    fprintf(fp, "3  CbcHeuristicPivotAndFix heuristicPFX(*cbcModel);\n");
    CbcHeuristic::generateCpp(fp, "heuristicPFX");
    fprintf(fp, "3  cbcModel->addHeuristic(&heuristicPFX);\n");
}

// Copy constructor
CbcHeuristicPivotAndFix::CbcHeuristicPivotAndFix(const CbcHeuristicPivotAndFix & rhs)
        :
        CbcHeuristic(rhs)
{
}

// Assignment operator
CbcHeuristicPivotAndFix &
CbcHeuristicPivotAndFix::operator=( const CbcHeuristicPivotAndFix & rhs)
{
    if (this != &rhs) {
        CbcHeuristic::operator=(rhs);
    }
    return *this;
}

// Resets stuff if model changes
void
CbcHeuristicPivotAndFix::resetModel(CbcModel * /*model*/)
{
    //CbcHeuristic::resetModel(model);
}
/*
  Comments needed
  Returns 1 if solution, 0 if not */
int
CbcHeuristicPivotAndFix::solution(double & /*solutionValue*/,
                                  double * /*betterSolution*/)
{

    numCouldRun_++; // Todo: Ask JJHF what this for.
    std::cout << "Entering Pivot-and-Fix Heuristic" << std::endl;

#ifdef FORNOW
    std::cout << "Lucky you! You're in the Pivot-and-Fix Heuristic" << std::endl;
    // The struct should be moved to member data



    typedef struct {
        int numberSolutions;
        int maximumSolutions;
        int numberColumns;
        double ** solution;
        int * numberUnsatisfied;
    } clpSolution;

    double start = CoinCpuTime();

    OsiClpSolverInterface * clpSolverOriginal
    = dynamic_cast<OsiClpSolverInterface *> (model_->solver());
    assert (clpSolverOriginal);
    OsiClpSolverInterface *clpSolver(clpSolverOriginal);

    ClpSimplex * simplex = clpSolver->getModelPtr();

    // Initialize the structure holding the solutions
    clpSolution solutions;
    // Set typeStruct field of ClpTrustedData struct to one.
    // This tells Clp it's "Mahdi!"
    ClpTrustedData trustedSolutions;
    trustedSolutions.typeStruct = 1;
    trustedSolutions.data = &solutions;
    solutions.numberSolutions = 0;
    solutions.maximumSolutions = 0;
    solutions.numberColumns = simplex->numberColumns();
    solutions.solution = NULL;
    solutions.numberUnsatisfied = NULL;
    simplex->setTrustedUserPointer(&trustedSolutions);

    // Solve from all slack to get some points
    simplex->allSlackBasis();
    simplex->primal();

    // -------------------------------------------------
    // Get the problem information
    // - get the number of cols and rows
    int numCols = clpSolver->getNumCols();
    int numRows = clpSolver->getNumRows();

    // - get the right hand side of the rows
    const double * rhs = clpSolver->getRightHandSide();

    // - find the integer variables
    bool * varClassInt = new bool[numCols];
    int numInt = 0;
    for (int i = 0; i < numCols; i++) {
        if (clpSolver->isContinuous(i))
            varClassInt[i] = 0;
        else {
            varClassInt[i] = 1;
            numInt++;
        }
    }

    // -Get the rows sense
    const char * rowSense;
    rowSense = clpSolver->getRowSense();

    // -Get the objective coefficients
    const double *objCoefficients = clpSolver->getObjCoefficients();
    double *originalObjCoeff = new double [numCols];
    for (int i = 0; i < numCols; i++)
        originalObjCoeff[i] = objCoefficients[i];

    // -Get the matrix of the problem
    double ** matrix = new double * [numRows];
    for (int i = 0; i < numRows; i++) {
        matrix[i] = new double[numCols];
        for (int j = 0; j < numCols; j++)
            matrix[i][j] = 0;
    }
    const CoinPackedMatrix* matrixByRow = clpSolver->getMatrixByRow();
    const double * matrixElements = matrixByRow->getElements();
    const int * matrixIndices = matrixByRow->getIndices();
    const int * matrixStarts = matrixByRow->getVectorStarts();
    for (int j = 0; j < numRows; j++) {
        for (int i = matrixStarts[j]; i < matrixStarts[j+1]; i++) {
            matrix[j][matrixIndices[i]] = matrixElements[i];
        }
    }

    // The newObj is the randomly perturbed constraint used to find new
    // corner points
    double * newObj = new double [numCols];

    // Set the random seed
    srand ( time(NULL) + 1);
    int randNum;

    // We're going to add a new row to the LP formulation
    // after finding each new solution.
    // Adding a new row requires the new elements and the new indices.
    // The elements are original objective function coefficients.
    // The indicies are the (dense) columns indices stored in addRowIndex.
    // The rhs is the value of the new solution stored in solutionValue.
    int * addRowIndex = new int[numCols];
    for (int i = 0; i < numCols; i++)
        addRowIndex[i] = i;

    // The number of feasible solutions found by the PF heuristic.
    // This controls the return code of the solution() method.
    int numFeasibles = 0;

    // Shuffle the rows
    int * index = new int [numRows];
    for (int i = 0; i < numRows; i++)
        index[i] = i;
    for (int i = 0; i < numRows; i++) {
        int temp = index[i];
        int randNumTemp = i + (rand() % (numRows - i));
        index[i] = index[randNumTemp];
        index[randNumTemp] = temp;
    }

    // In the clpSolution struct, we store a lot of column solutions.
    // For each perturb objective, we store the solution from each
    // iteration of the LP solve.
    // For each perturb objective, we look at the collection of
    // solutions to do something extremly intelligent :-)
    // We could (and should..and will :-) wipe out the block of
    // solutions when we're done with them. But for now, we just move on
    // and store the next block of solutions for the next (perturbed)
    // objective.
    // The variable startIndex tells us where the new block begins.
    int startIndex = 0;

    // At most "fixThreshold" number of integer variables can be unsatisfied
    // for calling smallBranchAndBound().
    // The PF Heuristic only fixes fixThreshold number of variables to
    // their integer values. Not more. Not less. The reason is to give
    // the smallBB some opportunity to find better solutions. If we fix
    // everything it might be too many (leading the heuristic to come up
    // with infeasibility rather than a useful result).
    // (This is an important paramater. And it is dynamically set.)
    double fixThreshold;
    /*
        if(numInt > 400)
        fixThreshold = 17*sqrt(numInt);
        if(numInt<=400 && numInt>100)
        fixThreshold = 5*sqrt(numInt);
        if(numInt<=100)
        fixThreshold = 4*sqrt(numInt);
    */
    // Initialize fixThreshold based on the number of integer
    // variables
    if (numInt <= 100)
        fixThreshold = .35 * numInt;
    if (numInt > 100 && numInt < 1000)
        fixThreshold = .85 * numInt;
    if (numInt >= 1000)
        fixThreshold = .1 * numInt;

    // Whenever the dynamic system for changing fixThreshold
    // kicks in, it changes the parameter by the
    // fixThresholdChange amount.
    // (The 25% should be member data and tuned. Another paper!)
    double fixThresholdChange = 0.25 * fixThreshold;

    // maxNode is the maximum number of nodes we allow smallBB to
    // search. It's initialized to 400 and changed dynamically.
    // The 400 should be member data, if we become virtuous.
    int maxNode = 400;

    // We control the decision to change maxNode through the boolean
    // variable  changeMaxNode. The boolean variable is initialized to
    // true and gets set to false under a condition (and is never true
    // again.)
    // It's flipped off and stays off (in the current incarnation of PF)
    bool changeMaxNode = 1;

    // The sumReturnCode is used for the dynamic system that sets
    // fixThreshold and changeMaxNode.
    //
    // We track what's happening in sumReturnCode. There are 8 switches.
    // The first 5 switches corresponds to a return code for smallBB.
    //
    // We want to know how many times we consecutively get the same
    // return code.
    //
    // If "good" return codes are happening often enough, we're happy.
    //
    // If a "bad" returncodes happen consecutively, we want to
    // change something.
    //
    // The switch 5 is the number of times PF didn't call smallBB
    // becuase the number of integer variables that took integer values
    // was less than fixThreshold.
    //
    // The swicth 6 was added for a brilliant idea...to be announced
    // later (another paper!)
    //
    // The switch 7 is the one that changes the max node. Read the
    // code. (Todo: Verbalize the brilliant idea for the masses.)
    //
    int sumReturnCode[8];
    /*
      sumReturnCode[0] ~ -1 --> problem too big for smallBB
      sumReturnCode[1] ~ 0  --> smallBB not finshed and no soln
      sumReturnCode[2] ~ 1  --> smallBB not finshed and there is a soln
      sumReturnCode[3] ~ 2  --> smallBB finished and no soln
      sumReturnCode[4] ~ 3  --> smallBB finished and there is a soln
      sumReturnCode[5] ~ didn't call smallBranchAndBound too few to fix
      sumReturnCode[6] ~ didn't call smallBranchAndBound too many unsatisfied
      sumReturnCode[7] ~ the same as sumReturnCode[1] but becomes zero just if the returnCode is not 0
    */

    for (int i = 0; i < 8; i++)
        sumReturnCode[i] = 0;
    int * colIndex = new int[numCols];
    for (int i = 0; i < numCols; i++)
        colIndex[i] = i;
    double cutoff = COIN_DBL_MAX;
    bool didMiniBB;

    // Main loop
    for (int i = 0; i < numRows; i++) {
        // track the number of mini-bb for the dynamic threshold setting
        didMiniBB = 0;

        for (int k = startIndex; k < solutions.numberSolutions; k++)
            //if the point has 0 unsatisfied variables; make sure it is
            //feasible. Check integer feasiblity and constraints.
            if (solutions.numberUnsatisfied[k] == 0) {
                double feasibility = 1;
                //check integer feasibility
                for (int icol = 0; icol < numCols; icol++) {
                    double closest = floor(solutions.solution[k][icol] + 0.5);
                    if (varClassInt[icol] && (fabs(solutions.solution[k][icol] - closest) > 1e-6)) {
                        feasibility = 0;
                        break;
                    }
                }
                //check if the solution satisfies the constraints
                for (int irow = 0; irow < numRows; irow++) {
                    double lhs = 0;
                    for (int j = 0; j < numCols; j++)
                        lhs += matrix[irow][j] * solutions.solution[k][j];
                    if (rowSense[irow] == 'L' && lhs > rhs[irow] + 1e-6) {
                        feasibility = 0;
                        break;
                    }
                    if (rowSense[irow] == 'G' && lhs < rhs[irow] - 1e-6) {
                        feasibility = 0;
                        break;
                    }
                    if (rowSense[irow] == 'E' && (lhs - rhs[irow] > 1e-6 || lhs - rhs[irow] < -1e-6)) {
                        feasibility = 0;
                        break;
                    }
                }

                //if feasible, find the objective value and set the cutoff
                // for the smallBB and add a new constraint to the LP
                // (and update the best solution found so far for the
                // return arguments)
                if (feasibility) {
                    double objectiveValue = 0;
                    for (int j = 0; j < numCols; j++)
                        objectiveValue += solutions.solution[k][j] * originalObjCoeff[j];
                    cutoff = objectiveValue;
                    clpSolver->addRow(numCols, addRowIndex, originalObjCoeff, -COIN_DBL_MAX, cutoff);

                    // Todo: pick up the best solution in the block (not
                    // the last).
                    solutionValue = objectiveValue;
                    for (int m = 0; m < numCols; m++)
                        betterSolution[m] = solutions.solution[k][m];
                    numFeasibles++;
                }
            }

        // Go through the block of solution and decide if to call smallBB
        for (int k = startIndex; k < solutions.numberSolutions; k++) {
            if (solutions.numberUnsatisfied[k] <= fixThreshold) {
                // get new copy
                OsiSolverInterface * newSolver;
                newSolver = new OsiClpSolverInterface(*clpSolver);
                newSolver->setObjSense(1);
                newSolver->setObjective(originalObjCoeff);
                int numberColumns = newSolver->getNumCols();
                int numFixed = 0;

                // Fix the first fixThreshold number of integer vars
                // that are satisfied
                for (int iColumn = 0 ; iColumn < numberColumns ; iColumn++) {
                    if (newSolver->isInteger(iColumn)) {
                        double value = solutions.solution[k][iColumn];
                        double intValue = floor(value + 0.5);
                        if (fabs(value - intValue) < 1.0e-5) {
                            newSolver->setColLower(iColumn, intValue);
                            newSolver->setColUpper(iColumn, intValue);
                            numFixed++;
                            if (numFixed > numInt - fixThreshold)
                                break;
                        }
                    }
                }
                COIN_DETAIL_PRINT(printf("numFixed: %d\n", numFixed));
                COIN_DETAIL_PRINT(printf("fixThreshold: %f\n", fixThreshold));
		COIN_DETAIL_PRINT(printf("numInt: %d\n", numInt));
                double *newSolution = new double[numCols];
                double newSolutionValue;

                // Call smallBB on the modified problem
                int returnCode = smallBranchAndBound(newSolver, maxNode, newSolution,
                                                     newSolutionValue, cutoff, "mini");

                // If smallBB found a solution, update the better
                // solution and solutionValue (we gave smallBB our
                // cutoff, so it only finds improving solutions)
                if (returnCode == 1 || returnCode == 3) {
                    numFeasibles ++;
                    solutionValue = newSolutionValue;
                    for (int m = 0; m < numCols; m++)
                        betterSolution[m] = newSolution[m];
                    COIN_DETAIL_PRINT(printf("cutoff: %f\n", newSolutionValue));
                    COIN_DETAIL_PRINT(printf("time: %.2lf\n", CoinCpuTime() - start));
                }
                didMiniBB = 1;
                COIN_DETAIL_PRINT(printf("returnCode: %d\n", returnCode));

                //Update sumReturnCode array
                for (int iRC = 0; iRC < 6; iRC++) {
                    if (iRC == returnCode + 1)
                        sumReturnCode[iRC]++;
                    else
                        sumReturnCode[iRC] = 0;
                }
                if (returnCode != 0)
                    sumReturnCode[7] = 0;
                else
                    sumReturnCode[7]++;
                if (returnCode == 1 || returnCode == 3) {
                    cutoff = newSolutionValue;
                    clpSolver->addRow(numCols, addRowIndex, originalObjCoeff, -COIN_DBL_MAX, cutoff);
                    COIN_DETAIL_PRINT(printf("******************\n\n*****************\n"));
                }
                break;
            }
        }

        if (!didMiniBB && solutions.numberSolutions - startIndex > 0) {
            sumReturnCode[5]++;
            for (int iRC = 0; iRC < 5; iRC++)
                sumReturnCode[iRC] = 0;
        }

        //Change "fixThreshold" if needed
        // using the data we've recorded in sumReturnCode
        if (sumReturnCode[1] >= 3)
            fixThreshold -= fixThresholdChange;
        if (sumReturnCode[7] >= 3 && changeMaxNode) {
            maxNode *= 5;
            changeMaxNode = 0;
        }
        if (sumReturnCode[3] >= 3 && fixThreshold < 0.95 * numInt)
            fixThreshold += fixThresholdChange;
        if (sumReturnCode[5] >= 4)
            fixThreshold += fixThresholdChange;
        if (sumReturnCode[0] > 3)
            fixThreshold -= fixThresholdChange;

        startIndex = solutions.numberSolutions;

        //Check if the maximum iterations limit is reached
        // rlh: Ask John how this is working with the change to trustedUserPtr.
        if (solutions.numberSolutions > 20000)
            break;

        // The first time in this loop PF solves orig LP.

        //Generate the random objective function
        randNum = rand() % 10 + 1;
        randNum = fmod(randNum, 2);
        for (int j = 0; j < numCols; j++) {
            if (randNum == 1)
                if (fabs(matrix[index[i]][j]) < 1e-6)
                    newObj[j] = 0.1;
                else
                    newObj[j] = matrix[index[i]][j] * 1.1;
            else if (fabs(matrix[index[i]][j]) < 1e-6)
                newObj[j] = -0.1;
            else
                newObj[j] = matrix[index[i]][j] * 0.9;
        }
        clpSolver->setObjective(newObj);
        if (rowSense[i] == 'L')
            clpSolver->setObjSense(-1);
        else
            // Todo #1: We don't need to solve the LPs to optimality.
            // We just need corner points.
            // There's a problem in stopping Clp that needs to be looked
            // into. So for now, we solve optimality.
            clpSolver->setObjSense(1);
        //	  simplex->setMaximumIterations(100);
        clpSolver->getModelPtr()->primal(1);
        //	  simplex->setMaximumIterations(100000);
#ifdef COIN_DETAIL
        printf("cutoff: %f\n", cutoff);
        printf("time: %.2f\n", CoinCpuTime() - start);
        for (int iRC = 0; iRC < 8; iRC++)
            printf("%d ", sumReturnCode[iRC]);
        printf("\nfixThreshold: %f\n", fixThreshold);
        printf("numInt: %d\n", numInt);
        printf("\n---------------------------------------------------------------- %d\n", i);
#endif

        //temp:
        if (i > 3) break;

    }

    COIN_DETAIL_PRINT(printf("Best Feasible Found: %f\n", cutoff));
    COIN_DETAIL_PRINT(printf("Total time: %.2f\n", CoinCpuTime() - start));

    if (numFeasibles == 0) {
        return 0;
    }



    // We found something better
    std::cout << "See you soon! You're leaving the Pivot-and-Fix Heuristic" << std::endl;
    std::cout << std::endl;

    return 1;
#endif

    return 0;

}




// update model
void CbcHeuristicPivotAndFix::setModel(CbcModel * )
{
    // probably same as resetModel
}


