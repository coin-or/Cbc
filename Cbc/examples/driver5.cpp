// $Id: driver5.cpp 2101 2014-12-03 17:43:20Z forrest $
// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>
#include <iomanip>

#include "CoinPragma.hpp"
#include "CbcModel.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CbcSolver.hpp"
#include "CbcSOS.hpp"

#include "CoinTime.hpp"

//#############################################################################

/************************************************************************

This main program shows how to take advantage of the standalone cbc in your program,
while still making major modifications.
First it reads in an integer model from an mps file
Then it initializes the integer model with cbc defaults
Then it calls CbcMain1 passing all parameters apart from first
Then it modifies all SOS objects on first CbcEvent
Then it solves
Finally it prints solution

************************************************************************/
/*
  This is obviously stupid - just to show what can be done -
  just branch halfway rather than use weights and solution values
 */
class CbcUserSOS : public CbcSOS {

public:
  // Default Constructor
  CbcUserSOS();

  /* Constructor from CbcSOS
     Just adds dummy integer
  */

  CbcUserSOS(CbcSOS object, int dummy);

  // Copy constructor
  CbcUserSOS(const CbcUserSOS &);

  /// Clone
  virtual CbcObject *clone() const;

  // Assignment operator
  CbcUserSOS &operator=(const CbcUserSOS &rhs);

  // Destructor
  virtual ~CbcUserSOS();

  /// Infeasibility - large is 0.5
  virtual double infeasibility(const OsiBranchingInformation *info,
    int &preferredWay) const;

  /// Creates a branching object
  virtual CbcBranchingObject *createCbcBranch(OsiSolverInterface *solver, const OsiBranchingInformation *info, int way);

private:
  /// data

  /// Dummy integer
  int dummyInteger_;
};
//##############################################################################

// Default Constructor
CbcUserSOS::CbcUserSOS()
  : CbcSOS()
  , dummyInteger_(0)
{
}

// Constructor from CbcSOS
CbcUserSOS::CbcUserSOS(CbcSOS object, int dummyInteger)
  : CbcSOS(object)
  , dummyInteger_(dummyInteger)
{
}

// Copy constructor
CbcUserSOS::CbcUserSOS(const CbcUserSOS &rhs)
  : CbcSOS(rhs)
{
  dummyInteger_ = rhs.dummyInteger_;
}

// Clone
CbcObject *
CbcUserSOS::clone() const
{
  return new CbcUserSOS(*this);
}

// Assignment operator
CbcUserSOS &
CbcUserSOS::operator=(const CbcUserSOS &rhs)
{
  if (this != &rhs) {
    CbcSOS::operator=(rhs);
    dummyInteger_ = rhs.dummyInteger_;
  }
  return *this;
}

// Destructor
CbcUserSOS::~CbcUserSOS()
{
}
/*
  Routine to calculate standard infeasibility of an SOS set and return a
  preferred branching direction.
  This is just a copy of CbcSOS 
*/
double
CbcUserSOS::infeasibility(const OsiBranchingInformation *info,
  int &preferredWay) const
{
  int j;
  int firstNonZero = -1;
  int lastNonZero = -1;
  OsiSolverInterface *solver = model_->solver();
  const double *solution = model_->testSolution();
  const double *lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  //double largestValue=0.0;
#ifndef ZERO_ODD_TOLERANCE
#define ZERO_SOS_TOLERANCE 1.0e-14
#else
#define ZERO_SOS_TOLERANCE ZERO_ODD_TOLERANCE
#endif
  //double integerTolerance =
  //  model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double integerTolerance = ZERO_SOS_TOLERANCE;
  double weight = 0.0;
  double sum = 0.0;

  // check bounds etc
  double lastWeight = -1.0e100;
  for (j = 0; j < numberMembers_; j++) {
    int iColumn = members_[j];
    /*
          The value used here (1.0e-7) is larger than the value enforced in the
          constructor.
        */

    if (lastWeight >= weights_[j] - 1.0e-7)
      throw CoinError("Weights too close together in SOS", "infeasibility", "CbcSOS");
    double value = CoinMax(lower[iColumn], solution[iColumn]);
    value = CoinMin(upper[iColumn], value);
    sum += value;
    if (fabs(value) > integerTolerance && (upper[iColumn] > 0.0 || oddValues_)) {
      weight += weights_[j] * value;
      if (firstNonZero < 0)
        firstNonZero = j;
      lastNonZero = j;
    }
  }
  /* ?? */
  preferredWay = 1;
  /*
  SOS1 allows one nonzero; SOS2 allows two consecutive nonzeros. Infeasibility
  is calculated as (.5)(range of nonzero values)/(number of members). So if
  the first and last elements of the set are nonzero, we have maximum
  infeasibility.
*/
  if (lastNonZero - firstNonZero >= sosType_) {
    // find where to branch
    if (!oddValues_)
      weight /= sum;
    else
      weight = 0.5 * (weights_[firstNonZero] + weights_[lastNonZero]);
    if (info->defaultDual_ >= 0.0 && info->usefulRegion_ && info->columnStart_) {
      assert(sosType_ == 1);
      int iWhere;
      for (iWhere = firstNonZero; iWhere < lastNonZero - 1; iWhere++) {
        if (weight < weights_[iWhere + 1]) {
          break;
        }
      }
      int jColumnDown = members_[iWhere];
      int jColumnUp = members_[iWhere + 1];
      int n = 0;
      CoinBigIndex j;
      double objMove = info->objective_[jColumnDown];
      for (j = info->columnStart_[jColumnDown];
           j < info->columnStart_[jColumnDown] + info->columnLength_[jColumnDown]; j++) {
        double value = info->elementByColumn_[j];
        int iRow = info->row_[j];
        info->indexRegion_[n++] = iRow;
        info->usefulRegion_[iRow] = value;
      }
      for (iWhere = firstNonZero; iWhere < lastNonZero; iWhere++) {
        int jColumn = members_[iWhere];
        double solValue = info->solution_[jColumn];
        if (!solValue)
          continue;
        objMove -= info->objective_[jColumn] * solValue;
        for (j = info->columnStart_[jColumn];
             j < info->columnStart_[jColumn] + info->columnLength_[jColumn]; j++) {
          double value = -info->elementByColumn_[j] * solValue;
          int iRow = info->row_[j];
          double oldValue = info->usefulRegion_[iRow];
          if (!oldValue) {
            info->indexRegion_[n++] = iRow;
          } else {
            value += oldValue;
            if (!value)
              value = 1.0e-100;
          }
          info->usefulRegion_[iRow] = value;
        }
      }
      const double *pi = info->pi_;
      const double *activity = info->rowActivity_;
      const double *lower = info->rowLower_;
      const double *upper = info->rowUpper_;
      double tolerance = info->primalTolerance_;
      double direction = info->direction_;
      shadowEstimateDown_ = objMove * direction;
      bool infeasible = false;
      for (int k = 0; k < n; k++) {
        int iRow = info->indexRegion_[k];
        double movement = info->usefulRegion_[iRow];
        // not this time info->usefulRegion_[iRow]=0.0;
        double valueP = pi[iRow] * direction;
        // if move makes infeasible then make at least default
        double newValue = activity[iRow] + movement;
        if (newValue > upper[iRow] + tolerance || newValue < lower[iRow] - tolerance) {
          shadowEstimateDown_ += fabs(movement) * CoinMax(fabs(valueP), info->defaultDual_);
          infeasible = true;
        }
      }
      if (shadowEstimateDown_ < info->integerTolerance_) {
        if (!infeasible) {
          shadowEstimateDown_ = 1.0e-10;
        } else
          shadowEstimateDown_ = info->integerTolerance_;
      }
      // And other way
      // take off
      objMove -= info->objective_[jColumnDown];
      for (j = info->columnStart_[jColumnDown];
           j < info->columnStart_[jColumnDown] + info->columnLength_[jColumnDown]; j++) {
        double value = -info->elementByColumn_[j];
        int iRow = info->row_[j];
        double oldValue = info->usefulRegion_[iRow];
        if (!oldValue) {
          info->indexRegion_[n++] = iRow;
        } else {
          value += oldValue;
          if (!value)
            value = 1.0e-100;
        }
        info->usefulRegion_[iRow] = value;
      }
      // add on
      objMove += info->objective_[jColumnUp];
      for (j = info->columnStart_[jColumnUp];
           j < info->columnStart_[jColumnUp] + info->columnLength_[jColumnUp]; j++) {
        double value = info->elementByColumn_[j];
        int iRow = info->row_[j];
        double oldValue = info->usefulRegion_[iRow];
        if (!oldValue) {
          info->indexRegion_[n++] = iRow;
        } else {
          value += oldValue;
          if (!value)
            value = 1.0e-100;
        }
        info->usefulRegion_[iRow] = value;
      }
      shadowEstimateUp_ = objMove * direction;
      infeasible = false;
      for (int k = 0; k < n; k++) {
        int iRow = info->indexRegion_[k];
        double movement = info->usefulRegion_[iRow];
        info->usefulRegion_[iRow] = 0.0;
        double valueP = pi[iRow] * direction;
        // if move makes infeasible then make at least default
        double newValue = activity[iRow] + movement;
        if (newValue > upper[iRow] + tolerance || newValue < lower[iRow] - tolerance) {
          shadowEstimateUp_ += fabs(movement) * CoinMax(fabs(valueP), info->defaultDual_);
          infeasible = true;
        }
      }
      if (shadowEstimateUp_ < info->integerTolerance_) {
        if (!infeasible) {
          shadowEstimateUp_ = 1.0e-10;
        } else
          shadowEstimateUp_ = info->integerTolerance_;
      }
      // adjust
      double downCost = shadowEstimateDown_;
      double upCost = shadowEstimateUp_;
      if (numberTimesDown_)
        downCost *= downDynamicPseudoRatio_ / static_cast< double >(numberTimesDown_);
      if (numberTimesUp_)
        upCost *= upDynamicPseudoRatio_ / static_cast< double >(numberTimesUp_);
#define WEIGHT_AFTER 0.7
#define WEIGHT_BEFORE 0.1
      int stateOfSearch = model_->stateOfSearch() % 10;
      double returnValue = 0.0;
      double minValue = CoinMin(downCost, upCost);
      double maxValue = CoinMax(downCost, upCost);
      if (stateOfSearch <= 2) {
        // no branching solution
        returnValue = WEIGHT_BEFORE * minValue + (1.0 - WEIGHT_BEFORE) * maxValue;
      } else {
        returnValue = WEIGHT_AFTER * minValue + (1.0 - WEIGHT_AFTER) * maxValue;
      }
      return returnValue;
    } else {
      double value = lastNonZero - firstNonZero + 1;
      value *= 0.5 / static_cast< double >(numberMembers_);
      return value;
    }
  } else {
    return 0.0; // satisfied
  }
}

CbcBranchingObject *
CbcUserSOS::createCbcBranch(OsiSolverInterface *solver, const OsiBranchingInformation * /*info*/, int way)
{
  int j;
  const double *solution = model_->testSolution();
#ifndef ZERO_SOS_TOLERANCE
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
#else
  double integerTolerance = ZERO_SOS_TOLERANCE;
#endif
  const double *lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  int firstNonZero = -1;
  int lastNonZero = -1;
  for (j = 0; j < numberMembers_; j++) {
    int iColumn = members_[j];
    double value = CoinMax(lower[iColumn], solution[iColumn]);
    value = CoinMin(upper[iColumn], value);
    if (fabs(value) > integerTolerance) {
      if (firstNonZero < 0)
        firstNonZero = j;
      lastNonZero = j;
    }
  }
  assert(lastNonZero - firstNonZero >= sosType_);
  // find where to branch
  // stupid but just to show
  double separator;
  int iWhere = (firstNonZero + lastNonZero) / 2;
  if (sosType_ == 1) {
    // SOS 1
    separator = 0.5 * (weights_[iWhere] + weights_[iWhere + 1]);
  } else {
    // SOS 2
    if (iWhere + 1 == lastNonZero)
      iWhere--;
    separator = weights_[iWhere + 1];
  }
  // create object
  CbcBranchingObject *branch;
  branch = new CbcSOSBranchingObject(model_, this, way, separator);
  branch->setOriginalObject(this);
  return branch;
}
#include "CbcEventHandler.hpp"
/** This is so user can trap events and do useful stuff.  

    CbcModel model_ is available as well as anything else you care 
    to pass in
*/

class MyEventHandler3 : public CbcEventHandler {

public:
  /**@name Overrides */
  //@{
  virtual CbcAction event(CbcEvent whichEvent);
  //@}

  /**@name Constructors, destructor etc*/
  //@{
  /** Default constructor. */
  MyEventHandler3();
  /// Constructor with pointer to model (redundant as setEventHandler does)
  MyEventHandler3(CbcModel *model);
  /** Destructor */
  virtual ~MyEventHandler3();
  /** The copy constructor. */
  MyEventHandler3(const MyEventHandler3 &rhs);
  /// Assignment
  MyEventHandler3 &operator=(const MyEventHandler3 &rhs);
  /// Clone
  virtual CbcEventHandler *clone() const;
  //@}

protected:
  // data goes here
};
//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
MyEventHandler3::MyEventHandler3()
  : CbcEventHandler()
{
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
MyEventHandler3::MyEventHandler3(const MyEventHandler3 &rhs)
  : CbcEventHandler(rhs)
{
}

// Constructor with pointer to model
MyEventHandler3::MyEventHandler3(CbcModel *model)
  : CbcEventHandler(model)
{
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
MyEventHandler3::~MyEventHandler3()
{
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
MyEventHandler3 &
MyEventHandler3::operator=(const MyEventHandler3 &rhs)
{
  if (this != &rhs) {
    CbcEventHandler::operator=(rhs);
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CbcEventHandler *MyEventHandler3::clone() const
{
  return new MyEventHandler3(*this);
}
static bool firstTime = true;
CbcEventHandler::CbcAction
MyEventHandler3::event(CbcEvent whichEvent)
{
  // If in sub tree carry on
  if (!model_->parentModel()) {
    if (firstTime) {
      // Replace all SOS with home grown version
      // Could of course do same for integers
      firstTime = false;
      int numberObjects = model_->numberObjects();
      OsiObject **objects = model_->objects();
      for (int i = 0; i < numberObjects; i++) {
        CbcSOS *sosObj = dynamic_cast< CbcSOS * >(objects[i]);
        if (sosObj) {
          objects[i] = new CbcUserSOS(*sosObj, i);
          delete sosObj;
        }
      }
    }
  }
  return noAction; // carry on
}

int main(int argc, const char *argv[])
{

  OsiClpSolverInterface solver1;
  // Read in model using argv[1]
  // and assert that it is a clean model
  if (argc < 2) {
    fprintf(stderr, "Need input file name.\n");
    exit(1);
  }
  int numReadErrors;
  if (!strstr(argv[1], ".lp"))
    numReadErrors = solver1.readMps(argv[1], "");
  else
    numReadErrors = solver1.readLp(argv[1], 1.0e-10);
  if (numReadErrors != 0) {
    printf("%d errors reading file\n", numReadErrors);
    return numReadErrors;
  }

  // Messy code below copied from CbcSolver.cpp
  // Pass to Cbc initialize defaults
  CbcModel modelA(solver1);
  CbcModel *model = &modelA;
  CbcMain0(modelA);
  // Event handler
  MyEventHandler3 eventHandler;
  model->passInEventHandler(&eventHandler);
  /* Now go into code for standalone solver
     Could copy arguments and add -quit at end to be safe
     but this will do
  */
  if (argc > 2) {
    CbcMain1(argc - 1, argv + 1, modelA);
  } else {
    const char *argv2[] = { "driver5", "-solve", "-quit" };
    CbcMain1(3, argv2, modelA);
  }
  // Solver was cloned so get current copy
  OsiSolverInterface *solver = model->solver();
  // Print solution if finished (could get from model->bestSolution() as well

  if (model->bestSolution()) {

    const double *solution = solver->getColSolution();

    int iColumn;
    int numberColumns = solver->getNumCols();
    std::cout << std::setiosflags(std::ios::fixed | std::ios::showpoint) << std::setw(14);

    std::cout << "--------------------------------------" << std::endl;
    // names may not be in current solver - use original

    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value = solution[iColumn];
      if (fabs(value) > 1.0e-7)
        std::cout << std::setw(6) << iColumn << " " << std::setw(8) << setiosflags(std::ios::left) << solver1.getModelPtr()->columnName(iColumn)
                  << resetiosflags(std::ios::adjustfield) << std::setw(14) << " " << value << std::endl;
    }
    std::cout << "--------------------------------------" << std::endl;

    std::cout << std::resetiosflags(std::ios::fixed | std::ios::showpoint | std::ios::scientific);
  } else {
    std::cout << " No solution!" << std::endl;
  }
  return 0;
}
