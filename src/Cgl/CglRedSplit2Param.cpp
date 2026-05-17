// Name:     CglRedSplit2Param.hpp
// Author:   Giacomo Nannicini
//           Singapore University of Technology and Design
//           Singapore
//           email: nannicini@sutd.edu.sg
// Date:     03/09/09
//-----------------------------------------------------------------------------
// Copyright (C) 2010, Giacomo Nannicini and others.  All Rights Reserved.

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <iostream>

#include "CoinPragma.hpp"
#include "CglRedSplit2Param.hpp"

/***********************************************************************/
void CglRedSplit2Param::setAway(double value)
{
  if (value > 0.0 && value <= 0.5)
    away_ = value;
}

/***********************************************************************/
void CglRedSplit2Param::setEPS_ELIM(double eps_el)
{
  if(eps_el >= 0)
    EPS_ELIM = eps_el;
} /* setEPS_ELIM */

/***********************************************************************/
void CglRedSplit2Param::setEPS_RELAX_ABS(double eps_ra)
{
  if(eps_ra >= 0)
    EPS_RELAX_ABS = eps_ra;
} /* setEPS_RELAX_ABS */

/***********************************************************************/
void CglRedSplit2Param::setEPS_RELAX_REL(double eps_rr)
{
  if(eps_rr >= 0)
    EPS_RELAX_REL = eps_rr;
} /* setEPS_RELAX_REL */

/***********************************************************************/
void CglRedSplit2Param::setMAXDYN(double value)
{
    if (value > 1.0) {
    MAXDYN = value;
  }
  else {
    printf("### WARNING: CglRedSplit2::setMAXDYN(): value: %f ignored\n", 
	   value);
  }
} /* setMAXDYN */

/***********************************************************************/
void CglRedSplit2Param::setMINVIOL(double value)
{
  if (value > 0.0 && value <= 0.1) {
    MINVIOL = value;
  }
  else {
    printf("### WARNING: CglRedSplit2Param::setMINVIOL(): value: %f ignored\n", 
	   value);
  }
} /* setMINVIOL */

/***********************************************************************/
void CglRedSplit2Param::setMAX_SUPP_REL(double value)
{
  if (value > 0.0 && value <= 0.1) {
    MINVIOL = value;
  }
  else {
    printf("### WARNING: CglRedSplit2Param::setMINVIOL(): value: %f ignored\n", 
	   value);
  }
} /* setMAX_SUPP_REL */

/***********************************************************************/
void CglRedSplit2Param::setUSE_INTSLACKS(int value)
{
  USE_INTSLACKS = value;
} /* setUSE_INTSLACKS */

/***********************************************************************/
void CglRedSplit2Param::setNormIsZero(double value)
{
  if (value > 0.0 && value <= 1) {
    normIsZero_ = value;
  }
  else {
    printf("### WARNING: CglRedSplit2Param::setNormIsZero(): value: %f ignored\n",
	   value);
  }
} /* setNormIsZero */

/***********************************************************************/
void CglRedSplit2Param::setMinNormReduction(double value)
{
  if (value > 0.0 && value <= 1) {
    minNormReduction_ = value;
  }
  else {
    printf("### WARNING: CglRedSplit2Param::setMinNormReduction(): value: %f ignored\n",
	   value);
  }
} /* setMinNormReduction */

/***********************************************************************/
void CglRedSplit2Param::setMaxSumMultipliers(int value)
{
  if (value > 1) {
    maxSumMultipliers_ = value;
  }
  else {
    printf("### WARNING: CglRedSplit2Param::setMaxSumMultipliers(): value: %d ignored\n",
	   value);
  }
} /* setMaxSumMultipliers */

/***********************************************************************/
void CglRedSplit2Param::setNormalization(double value)
{
  if (value >= 0) {
    normalization_ = value;
  }
  else {
    printf("### WARNING: CglRedSplit2Param::setNormalization(): value: %f ignored\n",
	   value);
  }
} /* setNormalization */

/***********************************************************************/
void CglRedSplit2Param::addNumRowsReduction(int value)
{
  if (value >= 0) {
    numRowsReduction_.push_back(value);
  }
  else {
    printf("### WARNING: CglRedSplit2Param::addNumRowsReduction(): value: %d ignored\n",
	   value);
  }
} /* addNumRowsReduction */

/***********************************************************************/
void CglRedSplit2Param::addColumnSelectionStrategy(ColumnSelectionStrategy value)
{
  if (value != CS_ALL && value != CS_BEST && value != CS_LAP_NONBASICS) {
    columnSelectionStrategy_.push_back(value);
  }
  else if (value == CS_ALL) {
    // CS_ALL means all strategies.
    columnSelectionStrategy_.push_back(CS1);
    columnSelectionStrategy_.push_back(CS2);
    columnSelectionStrategy_.push_back(CS3);
    columnSelectionStrategy_.push_back(CS4);
    columnSelectionStrategy_.push_back(CS5);
    columnSelectionStrategy_.push_back(CS6);
    columnSelectionStrategy_.push_back(CS7);
    columnSelectionStrategy_.push_back(CS8);
    columnSelectionStrategy_.push_back(CS9);
    columnSelectionStrategy_.push_back(CS10);
    columnSelectionStrategy_.push_back(CS11);
    columnSelectionStrategy_.push_back(CS12);
    columnSelectionStrategy_.push_back(CS13);
    columnSelectionStrategy_.push_back(CS14);
    columnSelectionStrategy_.push_back(CS15);
    columnSelectionStrategy_.push_back(CS16);
    columnSelectionStrategy_.push_back(CS17);
    columnSelectionStrategy_.push_back(CS18);
    columnSelectionStrategy_.push_back(CS19);
    columnSelectionStrategy_.push_back(CS20);
    columnSelectionStrategy_.push_back(CS21);
  }
  else if (value == CS_BEST){
    // Select C-5P, I-2P-2/3, I-2P-4/5, I-3P
    columnSelectionStrategy_.push_back(CS4);
    columnSelectionStrategy_.push_back(CS5);
    columnSelectionStrategy_.push_back(CS6);
    columnSelectionStrategy_.push_back(CS7);
    columnSelectionStrategy_.push_back(CS8);

    columnSelectionStrategy_.push_back(CS9);
    columnSelectionStrategy_.push_back(CS10);
    columnSelectionStrategy_.push_back(CS11);
    columnSelectionStrategy_.push_back(CS12);

    columnSelectionStrategy_.push_back(CS18);
    columnSelectionStrategy_.push_back(CS19);
    columnSelectionStrategy_.push_back(CS20);
    columnSelectionStrategy_.push_back(CS21);
  }

} /* addColumnSelectionStrategy */

/***********************************************************************/
void CglRedSplit2Param::addRowSelectionStrategy(RowSelectionStrategy value)
{
  if (value != RS_ALL && value != RS_BEST){
    rowSelectionStrategy_.push_back(value);
  }
  else if (value == RS_ALL){
    rowSelectionStrategy_.push_back(RS1);
    rowSelectionStrategy_.push_back(RS2);
    rowSelectionStrategy_.push_back(RS3);
    rowSelectionStrategy_.push_back(RS4);
    rowSelectionStrategy_.push_back(RS5);
    rowSelectionStrategy_.push_back(RS6);
    rowSelectionStrategy_.push_back(RS7);
    rowSelectionStrategy_.push_back(RS8);
  }
  else if (value == RS_BEST){
    rowSelectionStrategy_.push_back(RS7);
    rowSelectionStrategy_.push_back(RS8);
  }
} /* addRowSelectionStrategy */

/***********************************************************************/

void CglRedSplit2Param::addNumRowsReductionLAP(int value)
{
  if (value >= 0) {
    numRowsReductionLAP_.push_back(value);
  }
  else {
    printf("### WARNING: CglRedSplit2Param::addNumRowsReductionLAP(): value: %d ignored\n", value);
  }
} /* addNumRowsReductionLAP */

/***********************************************************************/
void CglRedSplit2Param::addColumnSelectionStrategyLAP(ColumnSelectionStrategy value)
{
  if (value != CS_ALL && value != CS_BEST) {
    columnSelectionStrategyLAP_.push_back(value);
  }
  else if (value == CS_BEST){
    columnSelectionStrategyLAP_.push_back(CS1);
  }
  else{
    printf("### WARNING: CglRedSplit2Param::addColumnSelectionStrategyLAP(): value: %d ignored\n", value);
  }
} /* addColumnSelectionStrategyLAP */

/***********************************************************************/
void CglRedSplit2Param::addRowSelectionStrategyLAP(RowSelectionStrategy value)
{
  if (value != RS_ALL && value != RS_BEST){
    rowSelectionStrategyLAP_.push_back(value);
  }
  else if (value == RS_BEST){
    rowSelectionStrategyLAP_.push_back(RS8);
  }
  else{
    printf("### WARNING: CglRedSplit2Param::addRowSelectionStrategyLAP(): value: %d ignored\n", value);
  }
} /* addRowSelectionStrategyLAP */

/***********************************************************************/

void CglRedSplit2Param::setTimeLimit(double value)
{
  if (value >= 0.0){
    timeLimit_ = value;
  }
  else {
    timeLimit_ = 0.0;
  }
} /* setTimeLimit */

/***********************************************************************/

void CglRedSplit2Param::setMaxNumCuts(int value)
{
  if (value >= 0){
    maxNumCuts_ = value;
  }
  else {
    printf("### WARNING: CglRedSplit2Param::maxNumCuts(): value: %d ignored\n",
	   value);
  }
} /* setMaxNumCuts */

/***********************************************************************/

void CglRedSplit2Param::setMaxNumComputedCuts(int value)
{
  if (value >= 0){
    maxNumComputedCuts_ = value;
  }
  else {
    printf("### WARNING: CglRedSplit2Param::maxNumComputedCuts(): value: %d ignored\n",
	   value);
  }
} /* setMaxNumComputedCuts */

/***********************************************************************/

void CglRedSplit2Param::setMaxNonzeroesTab(int value)
{
  if (value >= 0){
    maxNonzeroesTab_ = value;
  }
  else {
    printf("### WARNING: CglRedSplit2Param::maxNonzeroesTab(): value: %d ignored\n",
	   value);
  }
} /* setMaxNumComputedCuts */

/***********************************************************************/

void CglRedSplit2Param::setColumnScalingStrategyLAP(ColumnScalingStrategy value)
{
  columnScalingStrategyLAP_ = value;
} /* setColumnScalingStrategyLAP */

/***********************************************************************/

void CglRedSplit2Param::setColumnScalingBoundLAP(double value)
{
  if (value >= 0){
    columnScalingBoundLAP_ = value;
  }
  else {
    printf("### WARNING: CglRedSplit2Param::columnScalingBoundLAP(): value: %f ignored\n",
	   value);
  }
} /* setColumnScalingBoundLAP */

/***********************************************************************/

void CglRedSplit2Param::setSkipGomory(int value)
{
  if (value >= 0 && value <= 1){
    skipGomory_ = value;
  }
  else{
    printf("### WARNING: CglRedSplit2Param::skipGomory(): value: %d ignored\n",
	   value);
  }
}

/***********************************************************************/
CglRedSplit2Param::CglRedSplit2Param(bool use_default_strategies,
				     double eps,
				     double eps_coeff,
				     double eps_elim,
				     double eps_relax_abs,
				     double eps_relax_rel,
				     double max_dyn,
				     double min_viol,
				     int max_supp_abs,
				     double max_supp_rel,
				     int use_int_slacks,
				     double norm_zero,
				     double minNormReduction,
				     int maxSumMultipliers,
				     double normalization,
				     double away,
				     double timeLimit,
				     int maxNumCuts,
				     int maxNumComputedCuts,
				     int maxNonzeroesTab,
				     double columnScalingBoundLAP,
				     int skipGomory) :
  CglParam(COIN_DBL_MAX, eps, eps_coeff, max_supp_abs),
  EPS_ELIM(eps_elim),
  EPS_RELAX_ABS(eps_relax_abs),
  EPS_RELAX_REL(eps_relax_rel),
  MAXDYN(max_dyn),
  MINVIOL(min_viol),
  MAX_SUPP_REL(max_supp_rel),
  USE_INTSLACKS(use_int_slacks),
  normIsZero_(norm_zero),
  minNormReduction_(minNormReduction),
  maxSumMultipliers_(maxSumMultipliers),
  normalization_(normalization),
  away_(away),
  columnScalingBoundLAP_(columnScalingBoundLAP),
  timeLimit_(timeLimit),
  maxNumCuts_(maxNumCuts),
  maxNumComputedCuts_(maxNumComputedCuts),
  maxNonzeroesTab_(maxNonzeroesTab),
  skipGomory_(skipGomory)
{
  if (use_default_strategies) {
    addNumRowsReduction(5);
    addColumnSelectionStrategy(CglRedSplit2Param::CS_BEST);
    addRowSelectionStrategy(CglRedSplit2Param::RS_BEST);
    addNumRowsReductionLAP(3);
    addColumnSelectionStrategyLAP(CglRedSplit2Param::CS1);
    addRowSelectionStrategyLAP(CglRedSplit2Param::RS8);
    setColumnScalingStrategyLAP(CglRedSplit2Param::SC_UNIFORM_NZ);
  }
}

/***********************************************************************/
CglRedSplit2Param::CglRedSplit2Param(const CglParam &source,
				     bool use_default_strategies,
				     double eps_elim, 
				     double eps_relax_abs, 
				     double eps_relax_rel, 
				     double max_dyn,
				     double min_viol,
				     double max_supp_rel,
				     int use_int_slacks,
				     double norm_zero,
				     double minNormReduction,
				     int maxSumMultipliers,
				     double normalization,
				     double away,
				     double timeLimit,
				     int maxNumCuts,
				     int maxNumComputedCuts,
				     int maxNonzeroesTab,
				     double columnScalingBoundLAP,
				     int skipGomory) :
  CglParam(source),
  EPS_ELIM(eps_elim),
  EPS_RELAX_ABS(eps_relax_abs),
  EPS_RELAX_REL(eps_relax_rel),
  MAXDYN(max_dyn),
  MINVIOL(min_viol),
  MAX_SUPP_REL(max_supp_rel),
  USE_INTSLACKS(use_int_slacks),
  normIsZero_(norm_zero),
  minNormReduction_(minNormReduction),
  maxSumMultipliers_(maxSumMultipliers),
  normalization_(normalization),
  away_(away),
  columnScalingBoundLAP_(columnScalingBoundLAP),
  timeLimit_(timeLimit),
  maxNumCuts_(maxNumCuts),
  maxNumComputedCuts_(maxNumComputedCuts),
  maxNonzeroesTab_(maxNonzeroesTab),
  skipGomory_(skipGomory)
{
  if (use_default_strategies) {
    addNumRowsReduction(5);
    addColumnSelectionStrategy(CglRedSplit2Param::CS_BEST);
    addRowSelectionStrategy(CglRedSplit2Param::RS_BEST);
    addNumRowsReductionLAP(3);
    addColumnSelectionStrategyLAP(CglRedSplit2Param::CS1);
    addRowSelectionStrategyLAP(CglRedSplit2Param::RS8);
    setColumnScalingStrategyLAP(CglRedSplit2Param::SC_UNIFORM_NZ);
  }
}

/***********************************************************************/
CglRedSplit2Param::CglRedSplit2Param(const CglRedSplit2Param &source) :
  CglParam(source),
  EPS_ELIM(source.EPS_ELIM),
  EPS_RELAX_ABS(source.EPS_RELAX_ABS),
  EPS_RELAX_REL(source.EPS_RELAX_REL),
  MAXDYN(source.MAXDYN),
  MINVIOL(source.MINVIOL),
  MAX_SUPP_REL(source.MAX_SUPP_REL),
  USE_INTSLACKS(source.USE_INTSLACKS),
  normIsZero_(source.normIsZero_),
  minNormReduction_(source.minNormReduction_),
  maxSumMultipliers_(source.maxSumMultipliers_),
  normalization_(source.normalization_),
  away_(source.away_),
  numRowsReduction_(source.numRowsReduction_),
  columnSelectionStrategy_(source.columnSelectionStrategy_),
  rowSelectionStrategy_(source.rowSelectionStrategy_),
  numRowsReductionLAP_(source.numRowsReductionLAP_),
  columnSelectionStrategyLAP_(source.columnSelectionStrategyLAP_),
  rowSelectionStrategyLAP_(source.rowSelectionStrategyLAP_),
  columnScalingStrategyLAP_(source.columnScalingStrategyLAP_),
  columnScalingBoundLAP_(source.columnScalingBoundLAP_),
  timeLimit_(source.timeLimit_),
  maxNumCuts_(source.maxNumCuts_),
  maxNumComputedCuts_(source.maxNumComputedCuts_),
  maxNonzeroesTab_(source.maxNonzeroesTab_),
  skipGomory_(source.skipGomory_)
{}

/***********************************************************************/
CglRedSplit2Param* CglRedSplit2Param::clone() const
{
  return new CglRedSplit2Param(*this);
}

/***********************************************************************/
CglRedSplit2Param& CglRedSplit2Param::operator=(const CglRedSplit2Param &rhs)
{
  if(this != &rhs) {
    CglParam::operator=(rhs);

    EPS_ELIM = rhs.EPS_ELIM;
    EPS_RELAX_ABS = rhs.EPS_RELAX_ABS;
    EPS_RELAX_REL = rhs.EPS_RELAX_REL;
    MAXDYN = rhs.MAXDYN;
    MINVIOL = rhs.MINVIOL;
    MAX_SUPP_REL = rhs.MAX_SUPP_REL;
    USE_INTSLACKS = rhs.USE_INTSLACKS;
    normIsZero_ = rhs.normIsZero_;
    minNormReduction_ = rhs.minNormReduction_;
    maxSumMultipliers_ = rhs.maxSumMultipliers_;
    normalization_ = rhs.normalization_;
    away_ = rhs.away_;
    numRowsReduction_ = rhs.numRowsReduction_;
    columnSelectionStrategy_ = rhs.columnSelectionStrategy_;
    rowSelectionStrategy_ = rhs.rowSelectionStrategy_;
    numRowsReductionLAP_ = rhs.numRowsReductionLAP_;
    columnSelectionStrategyLAP_ = rhs.columnSelectionStrategyLAP_;
    rowSelectionStrategyLAP_ = rhs.rowSelectionStrategyLAP_;
    columnScalingStrategyLAP_ = rhs.columnScalingStrategyLAP_;
    columnScalingBoundLAP_ = rhs.columnScalingBoundLAP_;
    timeLimit_ = rhs.timeLimit_;
    maxNumCuts_ = rhs.maxNumCuts_;
    maxNumComputedCuts_ = rhs.maxNumComputedCuts_;
    maxNonzeroesTab_ = rhs.maxNonzeroesTab_;
    skipGomory_ = rhs.skipGomory_;
  }
  return *this;
}

/***********************************************************************/
CglRedSplit2Param::~CglRedSplit2Param()
{}
