// Name:     CglRedSplitParam.cpp
// Author:   Francois Margot                                                  
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     11/24/06
//
//---------------------------------------------------------------------------
// Copyright (C) 2006, Francois Margot and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <iostream>

#include "CoinPragma.hpp"
#include "CglRedSplitParam.hpp"

/***********************************************************************/
void CglRedSplitParam::setAway(const double value)
{
  if (value > 0.0 && value <= 0.5)
    away_ = value;
}

/***********************************************************************/
void CglRedSplitParam::setMaxTab(const double value)
{
  if (value > 10) {
    maxTab_ = value;
  }
  else {
    printf("### WARNING: CglRedSplitParam::setMaxTab(): value: %f ignored\n", 
	   value);
  }
}

/***********************************************************************/
void CglRedSplitParam::setLUB(const double value)
{
  if (value > 0.0) {
    LUB = value;
  }
  else {
    printf("### WARNING: CglRedSplitParam::setLUB(): value: %f ignored\n", value);
  }
} /* setLUB */

/***********************************************************************/
void CglRedSplitParam::setEPS_ELIM(const double eps_el)
{
  if(eps_el >= 0)
    EPS_ELIM = eps_el;
} /* setEPS_ELIM */

/***********************************************************************/
void CglRedSplitParam::setEPS_RELAX_ABS(const double eps_ra)
{
  if(eps_ra >= 0)
    EPS_RELAX_ABS = eps_ra;
} /* setEPS_RELAX_ABS */

/***********************************************************************/
void CglRedSplitParam::setEPS_RELAX_REL(const double eps_rr)
{
  if(eps_rr >= 0)
    EPS_RELAX_REL = eps_rr;
} /* setEPS_RELAX_REL */

/***********************************************************************/
void CglRedSplitParam::setMAXDYN(double value)
{
    if (value > 1.0) {
    MAXDYN = value;
  }
  else {
    printf("### WARNING: CglRedSplit::setMAXDYN(): value: %f ignored\n", 
	   value);
  }
} /* setMAXDYN */

/***********************************************************************/
void CglRedSplitParam::setMAXDYN_LUB(double value)
{
  if (value > 1.0) {
    MAXDYN_LUB = value;
  }
  else {
    printf("### WARNING: CglRedSplit::setMAXDYN_LUB(): value: %f ignored\n", 
	   value);
  }
} /* setMAXDYN_LUB */

/***********************************************************************/
void CglRedSplitParam::setEPS_COEFF_LUB(const double value)
{
  if (value > 0.0 && value <= 0.1) {
    EPS_COEFF_LUB = value;
  }
  else {
    printf("### WARNING: CglRedSplitParam::setEPS_COEFF_LUB(): value: %f ignored\n", 
	   value);
  }
} /* setEPS_COEFF_LUB */

/***********************************************************************/
void CglRedSplitParam::setMINVIOL(double value)
{
  if (value > 0.0 && value <= 0.1) {
    MINVIOL = value;
  }
  else {
    printf("### WARNING: CglRedSplitParam::setMINVIOL(): value: %f ignored\n", 
	   value);
  }
} /* setMINVIOL */

/***********************************************************************/
void CglRedSplitParam::setUSE_INTSLACKS(int value)
{
  USE_INTSLACKS = value;
} /* setUSE_INTSLACKS */

/***********************************************************************/
void CglRedSplitParam::setUSE_CG2(int value)
{
  USE_CG2 = value;
} /* setUSE_CG2 */

/***********************************************************************/
void CglRedSplitParam::setNormIsZero(const double value)
{
  if (value > 0.0 && value <= 1) {
    normIsZero = value;
  }
  else {
    printf("### WARNING: CglRedSplitParam::setNormIsZero(): value: %f ignored\n",
	   value);
  }
} /* setNormIsZero */

/***********************************************************************/
void CglRedSplitParam::setMinReduc(const double value)
{
  if (value > 0.0 && value <= 1) {
    minReduc = value;
  }
  else {
    printf("### WARNING: CglRedSplitParam::MinReduc(): value: %f ignored\n",
	   value);
  }
} /* setMinReduc */

/***********************************************************************/
CglRedSplitParam::CglRedSplitParam(const double lub,
				   const double eps_el,
				   const double eps_relax_abs,
				   const double eps_relax_rel,
				   const double max_dyn,
				   const double max_dyn_lub,
				   const double eps_coeff_lub,
				   const double min_viol,
				   const int use_int_slacks,
				   const int use_cg2,
				   const double norm_zero,
				   const double min_reduc,
				   const double away,
				   const double max_tab) :
  CglParam(),
  LUB(lub),
  EPS_ELIM(eps_el),
  EPS_RELAX_ABS(eps_relax_abs),
  EPS_RELAX_REL(eps_relax_rel),
  MAXDYN(max_dyn),
  MAXDYN_LUB(max_dyn_lub),
  EPS_COEFF_LUB(eps_coeff_lub),
  MINVIOL(min_viol),
  USE_INTSLACKS(use_int_slacks),
  USE_CG2(use_cg2),
  normIsZero(norm_zero),
  minReduc(min_reduc),
  away_(away),
  maxTab_(max_tab)
{}

/***********************************************************************/
CglRedSplitParam::CglRedSplitParam(const CglParam &source,
				   const double lub,
				   const double eps_el, 
				   const double eps_ra, 
				   const double eps_rr, 
				   const double max_dyn,
				   const double max_dyn_lub,
				   const double eps_coeff_lub,
				   const double min_viol,
				   const int use_int_slacks,
				   const int use_cg2,
				   const double norm_zero,
				   const double min_reduc,
				   const double away,
				   const double max_tab) :

  CglParam(source), 
  LUB(lub),
  EPS_ELIM(eps_el),
  EPS_RELAX_ABS(eps_ra),
  EPS_RELAX_REL(eps_rr),
  MAXDYN(max_dyn),
  MAXDYN_LUB(max_dyn_lub),
  EPS_COEFF_LUB(eps_coeff_lub),
  MINVIOL(min_viol),
  USE_INTSLACKS(use_int_slacks),
  USE_CG2(use_cg2),
  normIsZero(norm_zero),
  minReduc(min_reduc),
  away_(away),
  maxTab_(max_tab)
{}

/***********************************************************************/
CglRedSplitParam::CglRedSplitParam(const CglRedSplitParam &source) :
  CglParam(source),
  LUB(source.LUB),
  EPS_ELIM(source.EPS_ELIM),
  EPS_RELAX_ABS(source.EPS_RELAX_ABS),
  EPS_RELAX_REL(source.EPS_RELAX_REL),
  MAXDYN(source.MAXDYN),
  MAXDYN_LUB(source.MAXDYN_LUB),
  EPS_COEFF_LUB(source.EPS_COEFF_LUB),
  MINVIOL(source.MINVIOL),
  USE_INTSLACKS(source.USE_INTSLACKS),
  USE_CG2(source.USE_CG2),
  normIsZero(source.normIsZero),
  minReduc(source.minReduc),
  away_(source.away_),
  maxTab_(source.maxTab_)
{}

/***********************************************************************/
CglRedSplitParam* CglRedSplitParam::clone() const
{
  return new CglRedSplitParam(*this);
}

/***********************************************************************/
CglRedSplitParam& CglRedSplitParam::operator=(const CglRedSplitParam &rhs)
{
  if(this != &rhs) {
    CglParam::operator=(rhs);

    LUB = rhs.LUB;
    EPS_ELIM = rhs.EPS_ELIM;
    EPS_RELAX_ABS = rhs.EPS_RELAX_ABS;
    EPS_RELAX_REL = rhs.EPS_RELAX_REL;
    MAXDYN = rhs.MAXDYN;
    MAXDYN_LUB = rhs.MAXDYN_LUB;
    EPS_COEFF_LUB = rhs.EPS_COEFF_LUB;
    MINVIOL = rhs.MINVIOL;
    USE_INTSLACKS = rhs.USE_INTSLACKS;
    USE_CG2 = rhs.USE_CG2;
    normIsZero = rhs.normIsZero;
    minReduc = rhs.minReduc;
    away_ = rhs.away_;
    maxTab_ = rhs.maxTab_;
  }
  return *this;
}

/***********************************************************************/
CglRedSplitParam::~CglRedSplitParam()
{}
