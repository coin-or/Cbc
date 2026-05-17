// Name:     CglGMIParam.hpp
// Author:   Giacomo Nannicini
//           Singapore University of Technology and Design
//           email: nannicini@sutd.edu.sg
//           based on CglRedSplitParam.cpp by Francois Margot
// Date:     11/17/09
//-----------------------------------------------------------------------------
// Copyright (C) 2009, Giacomo Nannicini and others.  All Rights Reserved.

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <iostream>

#include "CoinPragma.hpp"
#include "CglGMIParam.hpp"

/***********************************************************************/
void CglGMIParam::setAway(double value) {
  if (value > 0.0 && value <= 0.5) {
    AWAY = value;
  }
  else {
    printf("### WARNING: CglGMIParam::setAWAY(): value: %f ignored\n", value);
  }
}

/***********************************************************************/
void CglGMIParam::setEPS_ELIM(double value) {
  if (value >= 0) {
    EPS_ELIM = value;
  }
  else {
    printf("### WARNING: CglGMIParam::setEPS_ELIM(): value: %f ignored\n", value);
  }
} /* setEPS_ELIM */

/***********************************************************************/
void CglGMIParam::setEPS_RELAX_ABS(double value) {
  if (value >= 0) {
    EPS_RELAX_ABS = value;
  }
  else {
    printf("### WARNING: CglGMIParam::setEPS_RELAX_ABS(): value: %f ignored\n", value);
  }
} /* setEPS_RELAX_ABS */

/***********************************************************************/
void CglGMIParam::setEPS_RELAX_REL(double value) {
  if (value >= 0) {
    EPS_RELAX_REL = value;
  }
  else {
    printf("### WARNING: CglGMIParam::setEPS_RELAX_REL(): value: %f ignored\n", value);
  }
} /* setEPS_RELAX_REL */

/***********************************************************************/
void CglGMIParam::setMAXDYN(double value) {
  if (value >= 1.0) {
    MAXDYN = value;
  }
  else {
    printf("### WARNING: CglGMI::setMAXDYN(): value: %f ignored\n", 
	   value);
  }
} /* setMAXDYN */

/***********************************************************************/
void CglGMIParam::setMINVIOL(double value) {
  if (value >= 0.0) {
    MINVIOL = value;
  }
  else {
    printf("### WARNING: CglGMIParam::setMINVIOL(): value: %f ignored\n", 
	   value);
  }
} /* setMINVIOL */

/***********************************************************************/
void CglGMIParam::setMAX_SUPPORT_REL(double value) {
  if (value >= 0.0 && value <= 1.0) {
    MAX_SUPPORT_REL = value;
  }
  else {
    printf("### WARNING: CglGMIParam::setMAX_SUPPORT_REL(): value: %f ignored\n", 
	   value);
  }
} /* setMAX_SUPPORT_REL */

/***********************************************************************/
void CglGMIParam::setUSE_INTSLACKS(bool value) {
  USE_INTSLACKS = value;
} /* setUSE_INTSLACKS */

/***********************************************************************/
void CglGMIParam::setCHECK_DUPLICATES(bool value) {
  CHECK_DUPLICATES = value;
} /* setCHECK_DUPLICATES */

/***********************************************************************/
void CglGMIParam::setINTEGRAL_SCALE_CONT(bool value) {
  INTEGRAL_SCALE_CONT = value;
} /* setINTEGRAL_SCALE_CONT */

/***********************************************************************/
void CglGMIParam::setENFORCE_SCALING(bool value) {
  ENFORCE_SCALING = value;
} /* setENFORCE_SCALING */

/***********************************************************************/
void CglGMIParam::setCLEAN_PROC(CleaningProcedure value) {
  CLEAN_PROC = value;
} /* setCLEAN_PROC */

/***********************************************************************/
CglGMIParam::CglGMIParam(double eps,
			 double away,
			 double eps_coeff,
			 double eps_el,
			 double eps_relax_abs,
			 double eps_relax_rel,
			 double max_dyn,
			 double min_viol,
			 int max_supp_abs,
			 double max_supp_rel,
			 CleaningProcedure clean_proc,
			 bool use_int_slacks,
			 bool check_duplicates,
			 bool integral_scale_cont,
			 bool enforce_scaling) :
  CglParam(COIN_DBL_MAX, eps, eps_coeff, max_supp_abs),
  AWAY(away),
  EPS_ELIM(eps_el),
  EPS_RELAX_ABS(eps_relax_abs),
  EPS_RELAX_REL(eps_relax_rel),
  MAXDYN(max_dyn),
  MINVIOL(min_viol),
  MAX_SUPPORT_REL(max_supp_rel),
  CLEAN_PROC(clean_proc),
  USE_INTSLACKS(use_int_slacks),
  CHECK_DUPLICATES(check_duplicates),
  INTEGRAL_SCALE_CONT(integral_scale_cont),
  ENFORCE_SCALING(enforce_scaling)
{}

/***********************************************************************/
CglGMIParam::CglGMIParam(CglParam &source,
			 double away,
			 double eps_el, 
			 double eps_ra, 
			 double eps_rr, 
			 double max_dyn,
			 double min_viol,
			 double max_supp_rel,
			 CleaningProcedure clean_proc,
			 bool use_int_slacks,
			 bool check_duplicates,
			 bool integral_scale_cont,
			 bool enforce_scaling) :
  CglParam(source), 
  AWAY(away),
  EPS_ELIM(eps_el),
  EPS_RELAX_ABS(eps_ra),
  EPS_RELAX_REL(eps_rr),
  MAXDYN(max_dyn),
  MINVIOL(min_viol),
  MAX_SUPPORT_REL(max_supp_rel),
  CLEAN_PROC(clean_proc),
  USE_INTSLACKS(use_int_slacks),
  CHECK_DUPLICATES(check_duplicates),
  INTEGRAL_SCALE_CONT(integral_scale_cont),
  ENFORCE_SCALING(enforce_scaling)
{}

/***********************************************************************/
CglGMIParam::CglGMIParam(const CglGMIParam &source) :
  CglParam(source),
  AWAY(source.AWAY),
  EPS_ELIM(source.EPS_ELIM),
  EPS_RELAX_ABS(source.EPS_RELAX_ABS),
  EPS_RELAX_REL(source.EPS_RELAX_REL),
  MAXDYN(source.MAXDYN),
  MINVIOL(source.MINVIOL),
  MAX_SUPPORT_REL(source.MAX_SUPPORT_REL),
  CLEAN_PROC(source.CLEAN_PROC),
  USE_INTSLACKS(source.USE_INTSLACKS),
  CHECK_DUPLICATES(source.CHECK_DUPLICATES),
  INTEGRAL_SCALE_CONT(source.INTEGRAL_SCALE_CONT),
  ENFORCE_SCALING(source.ENFORCE_SCALING)
{}

/***********************************************************************/
CglGMIParam* CglGMIParam::clone() const
{
  return new CglGMIParam(*this);
}

/***********************************************************************/
CglGMIParam& CglGMIParam::operator=(const CglGMIParam &rhs)
{
  if(this != &rhs) {
    CglParam::operator=(rhs);

    AWAY = rhs.AWAY;
    EPS_ELIM = rhs.EPS_ELIM;
    EPS_RELAX_ABS = rhs.EPS_RELAX_ABS;
    EPS_RELAX_REL = rhs.EPS_RELAX_REL;
    MAXDYN = rhs.MAXDYN;
    MINVIOL = rhs.MINVIOL;
    MAX_SUPPORT_REL = rhs.MAX_SUPPORT_REL;
    CLEAN_PROC = rhs.CLEAN_PROC;
    USE_INTSLACKS = rhs.USE_INTSLACKS;
    CHECK_DUPLICATES = rhs.CHECK_DUPLICATES;
    INTEGRAL_SCALE_CONT = rhs.INTEGRAL_SCALE_CONT;
    ENFORCE_SCALING = rhs.ENFORCE_SCALING;
  }
  return *this;
}

/***********************************************************************/
CglGMIParam::~CglGMIParam()
{}
