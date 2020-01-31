/* $Id$ */
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif

#include <cassert>
#include <cstdio>

#include "CbcStatistics.hpp"
CbcStatistics &
CbcStatistics::operator=(const CbcStatistics &rhs)
{
  if (this != &rhs) {
    value_ = rhs.value_;
    startingObjective_ = rhs.startingObjective_;
    endingObjective_ = rhs.endingObjective_;
    id_ = rhs.id_;
    parentId_ = rhs.parentId_;
    way_ = rhs.way_;
    sequence_ = rhs.sequence_;
    depth_ = rhs.depth_;
    startingInfeasibility_ = rhs.startingInfeasibility_;
    endingInfeasibility_ = rhs.endingInfeasibility_;
    numberIterations_ = rhs.numberIterations_;
  }
  return *this;
}

CbcStatistics::CbcStatistics()
  : value_(0.0)
  , startingObjective_(0.0)
  , endingObjective_(COIN_DBL_MAX)
  , id_(-1)
  , parentId_(-1)
  , way_(0)
  , sequence_(-1)
  , depth_(0)
  , startingInfeasibility_(-1)
  , endingInfeasibility_(0)
  , numberIterations_(0)
{
}
// First or second branch
CbcStatistics::CbcStatistics(CbcNode *node, CbcModel *model)
  : endingObjective_(COIN_DBL_MAX)
  , endingInfeasibility_(0)
  , numberIterations_(0)
{
  CbcNodeInfo *nodeInfo = node->nodeInfo();
  CbcNodeInfo *parent = nodeInfo->parent();
  int numberBranches = nodeInfo->numberBranchesLeft();
  const CbcBranchingObject *branch = dynamic_cast< const CbcBranchingObject * >(node->branchingObject());
  const OsiTwoWayBranchingObject *branch2 = dynamic_cast< const OsiTwoWayBranchingObject * >(node->branchingObject());
  startingObjective_ = node->objectiveValue();
  way_ = node->way();
  depth_ = node->depth();
  startingInfeasibility_ = node->numberUnsatisfied();
  if (branch) {
    sequence_ = branch->variable();
    value_ = branch->value();
  } else {
    const OsiSimpleInteger *obj = dynamic_cast< const OsiSimpleInteger * >(branch2->originalObject());
    assert(obj);
    sequence_ = obj->columnNumber();
    value_ = branch2->value();
  }
  if (parent)
    parentId_ = parent->nodeNumber();
  else
    parentId_ = -1;
  if (numberBranches == 2) {
    id_ = nodeInfo->nodeNumber();
  } else {
    way_ *= 10;
    id_ = model->getNodeCount2();
  }
}

CbcStatistics::CbcStatistics(const CbcStatistics &rhs)
  : value_(rhs.value_)
  , startingObjective_(rhs.startingObjective_)
  , endingObjective_(rhs.endingObjective_)
  , id_(rhs.id_)
  , parentId_(rhs.parentId_)
  , way_(rhs.way_)
  , sequence_(rhs.sequence_)
  , depth_(rhs.depth_)
  , startingInfeasibility_(rhs.startingInfeasibility_)
  , endingInfeasibility_(rhs.endingInfeasibility_)
  , numberIterations_(rhs.numberIterations_)
{
}

CbcStatistics::~CbcStatistics()
{
}
// Update at end of branch
void CbcStatistics::endOfBranch(int numberIterations, double objectiveValue)
{
  numberIterations_ = numberIterations;
  endingObjective_ = objectiveValue;
}
// Update number of infeasibilities
void CbcStatistics::updateInfeasibility(int numberInfeasibilities)
{
  endingInfeasibility_ = numberInfeasibilities;
}
// Branch found to be infeasible by chooseBranch
void CbcStatistics::sayInfeasible()
{
  endingObjective_ = COIN_DBL_MAX;
}
// Just prints
void CbcStatistics::print(const int *sequenceLookup) const
{
  int sequence = -1;
  if (sequence_ >= 0)
    sequence = sequenceLookup ? sequenceLookup[sequence_] : sequence_;
  printf("%6d %6d %5d %6d %7.3f %s %s %13.7g (%5d) -> ",
    id_, parentId_, depth_, sequence, value_, abs(way_) == 1 ? " left" : "right",
    way_ < 0 ? "down" : " up ", startingObjective_, startingInfeasibility_);
  if (endingObjective_ != COIN_DBL_MAX)
    if (endingInfeasibility_)
      printf("%13.7g (%5d)\n", endingObjective_, endingInfeasibility_);
    else
      printf("%13.7g ** Solution\n", endingObjective_);
  else
    printf("cutoff\n");
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
