// Name:     CglGMIParam.hpp
// Author:   Giacomo Nannicini
//           Singapore University of Technology and Design
//           email: nannicini@sutd.edu.sg
//           based on CglRedSplitParam.hpp by Francois Margot
// Date:     11/17/09
//-----------------------------------------------------------------------------
// Copyright (C) 2009, Giacomo Nannicini and others.  All Rights Reserved.

#ifndef CglGMIParam_H
#define CglGMIParam_H

#include "CglParam.hpp"


  /**@name CglGMI Parameters */
  //@{

  /** Class collecting parameters for the GMI cut generator.

      Parameters of the generator are listed below. Modifying the default 
      values for parameters other than the last four might result in 
      invalid cuts.

      - MAXDYN: Maximum ratio between largest and smallest non zero 
                coefficients in a cut. See method setMAXDYN().
      - EPS_ELIM: Precision for deciding if a coefficient is zero when 
                  eliminating slack variables. See method setEPS_ELIM().
      - MINVIOL: Minimum violation for the current basic solution in 
                 a generated cut. See method setMINVIOL().
      - USE_INTSLACKS: Use integer slacks to generate cuts. 
                       (not implemented yet, will be in the future).
                       See method setUSE_INTSLACKS().
      - AWAY: Look only at basic integer variables whose current value is at
              least this value away from being integer. See method setAway().
      - CHECK_DUPLICATES: Should we check for duplicates when adding a cut
                          to the collection? Can be slow. 
                          Default 0 - do not check, add cuts anyway.
      - CLEAN_PROC: Cleaning procedure that should be used. Look below at the
                    enumeration CleaningProcedure for possible values.
      - INTEGRAL_SCALE_CONT: If we try to scale cut coefficients so that
                             they become integral, do we also scale on 
			     continuous variables?
			     Default 0 - do not scale continuous vars.
			     Used only if CLEAN_PROC does integral scaling.
      - ENFORCE_SCALING: Discard badly scaled cuts, or keep them (unscaled).
                         Default 1 - yes.

  */
  //@}

class CGLLIB_EXPORT CglGMIParam : public CglParam {

public:

  /**@name Enumerations */
  enum CleaningProcedure{
    /* CglLandP procedure I */
    CP_CGLLANDP1,
    /* CglLandP procedure II */
    CP_CGLLANDP2,
    /* CglRedSplit procedure I */
    CP_CGLREDSPLIT,
    /* Only integral cuts, i.e. cuts with integral coefficients */
    CP_INTEGRAL_CUTS,
    /* CglLandP procedure I with integral scaling */
    CP_CGLLANDP1_INT,
    /* CglLandP procedure I with scaling of the max element to 1 if possible */
    CP_CGLLANDP1_SCALEMAX,
    /* CglLandP procedure I with scaling of the rhs to 1 if possible */
    CP_CGLLANDP1_SCALERHS
  };

  /**@name Set/get methods */

  //@{
  /** Aliases for parameter get/set method in the base class CglParam */
  
  /** Value for Infinity. Default: DBL_MAX */
  inline void setInfinity(double value) {setINFINIT(value);}
  inline double getInfinity() const {return INFINIT;}

  /** Epsilon for comparing numbers. Default: 1.0e-6  */
  inline void setEps(double value) {setEPS(value);}
  inline double getEps() const {return EPS;}

  /** Epsilon for zeroing out coefficients. Default: 1.0e-5 */
  inline void setEpsCoeff(double value) {setEPS_COEFF(value);}
  inline double getEpsCoeff() const {return EPS_COEFF;}

  /** Maximum support of the cutting planes. Default: INT_MAX */
  inline void setMaxSupport(int value) {setMAX_SUPPORT(value);}
  inline int getMaxSupport() const {return MAX_SUPPORT;}
  /** Alias for consistency with our naming scheme */
  inline void setMaxSupportAbs(int value) {setMAX_SUPPORT(value);}
  inline int getMaxSupportAbs() const {return MAX_SUPPORT;}
  inline int getMAX_SUPPORT_ABS() const {return MAX_SUPPORT;}

  /** Set AWAY, the minimum distance from being integer used for selecting 
      rows for cut generation;  all rows whose pivot variable should be 
      integer but is more than away from integrality will be selected; 
      Default: 0.005 */
  virtual void setAway(double value);
  /** Get value of away */
  inline double getAway() const {return AWAY;}
  /// Aliases
  inline void setAWAY(double value) {setAway(value);}
  inline double getAWAY() const {return AWAY;}

  /** Set the value of EPS_ELIM, epsilon for values of coefficients when 
      eliminating slack variables;
      Default: 0 */
  virtual void setEPS_ELIM(double value);
  /** Get the value of EPS_ELIM */
  inline double getEPS_ELIM() const {return EPS_ELIM;}
  /// Aliases
  inline void setEpsElim(double value) {setEPS_ELIM(value);}
  inline double getEpsElim() const {return EPS_ELIM;}
  
  /** Set EPS_RELAX_ABS */
  virtual void setEPS_RELAX_ABS(double value);
  /** Get value of EPS_RELAX_ABS */
  inline double getEPS_RELAX_ABS() const {return EPS_RELAX_ABS;}
  /// Aliases
  inline void setEpsRelaxAbs(double value) {setEPS_RELAX_ABS(value);}
  inline double getEpsRelaxAbs() const {return EPS_RELAX_ABS;}

  /** Set EPS_RELAX_REL */
  virtual void setEPS_RELAX_REL(double value);
  /** Get value of EPS_RELAX_REL */
  inline double getEPS_RELAX_REL() const {return EPS_RELAX_REL;}
  /// Aliases
  inline void setEpsRelaxRel(double value) {setEPS_RELAX_REL(value);}
  inline double getEpsRelaxRel() const {return EPS_RELAX_REL;}

  // Set the maximum ratio between largest and smallest non zero 
  // coefficients in a cut. Default: 1e6.
  virtual void setMAXDYN(double value);
  /** Get the value of MAXDYN */
  inline double getMAXDYN() const {return MAXDYN;}
  /// Aliases
  inline void setMaxDyn(double value) {setMAXDYN(value);}
  inline double getMaxDyn() const {return MAXDYN;}

  /** Set the value of MINVIOL, the minimum violation for the current 
      basic solution in a generated cut. Default: 1e-7 */
  virtual void setMINVIOL(double value);
  /** Get the value of MINVIOL */
  inline double getMINVIOL() const {return MINVIOL;}
  /// Aliases
  inline void setMinViol(double value) {setMINVIOL(value);}
  inline double getMinViol() const {return MINVIOL;}

  /** Set the value of MAX_SUPPORT_REL, the factor contributing to the
      maximum support relative to the number of columns. Maximum
      allowed support is: MAX_SUPPORT_ABS +
      MAX_SUPPORT_REL*ncols. Default: 0.1 */
  virtual void setMAX_SUPPORT_REL(double value);
  /** Get the value of MINVIOL */
  inline double getMAX_SUPPORT_REL() const {return MAX_SUPPORT_REL;}
  /// Aliases
  inline void setMaxSupportRel(double value) {setMAX_SUPPORT_REL(value);}
  inline double getMaxSupportRel() const {return MAX_SUPPORT_REL;}

  /** Set the value of USE_INTSLACKS. Default: 0 */
  virtual void setUSE_INTSLACKS(bool value);
  /** Get the value of USE_INTSLACKS */
  inline bool getUSE_INTSLACKS() const {return USE_INTSLACKS;}
  /// Aliases
  inline void setUseIntSlacks(bool value) {setUSE_INTSLACKS(value);}
  inline int getUseIntSlacks() const {return USE_INTSLACKS;}

  /** Set the value of CHECK_DUPLICATES. Default: 0 */
  virtual void setCHECK_DUPLICATES(bool value);
  /** Get the value of CHECK_DUPLICATES */
  inline bool getCHECK_DUPLICATES() const {return CHECK_DUPLICATES;}
  /// Aliases
  inline void setCheckDuplicates(bool value) {setCHECK_DUPLICATES(value);}
  inline bool getCheckDuplicates() const {return CHECK_DUPLICATES;}

  /** Set the value of CLEAN_PROC. Default: CP_CGLLANDP1 */
  virtual void setCLEAN_PROC(CleaningProcedure value);
  /** Get the value of CLEAN_PROC. */
  inline CleaningProcedure getCLEAN_PROC() const {return CLEAN_PROC;}
  /// Aliases
  inline void setCleanProc(CleaningProcedure value) {setCLEAN_PROC(value);}
  inline CleaningProcedure getCleaningProcedure() const {return CLEAN_PROC;}

  /** Set the value of INTEGRAL_SCALE_CONT. Default: 0 */
  virtual void setINTEGRAL_SCALE_CONT(bool value);
  /** Get the value of INTEGRAL_SCALE_CONT. */
  inline bool getINTEGRAL_SCALE_CONT() const {return INTEGRAL_SCALE_CONT;}
  /// Aliases
  inline void setIntegralScaleCont(bool value) {setINTEGRAL_SCALE_CONT(value);}
  inline bool getIntegralScaleCont() const {return INTEGRAL_SCALE_CONT;}

  /** Set the value of ENFORCE_SCALING. Default: 1 */
  virtual void setENFORCE_SCALING(bool value);
  /** Get the value of ENFORCE_SCALING. */
  inline bool getENFORCE_SCALING() const {return ENFORCE_SCALING;}
  /// Aliases
  inline void setEnforceScaling(bool value) {setENFORCE_SCALING(value);}
  inline bool getEnforcescaling() const {return ENFORCE_SCALING;}

  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglGMIParam(double eps = 1e-12,
	      double away = 0.005,
	      double eps_coeff = 1e-11,
	      double eps_elim = 0,
	      double eps_relax_abs = 1e-11,
	      double eps_relax_rel = 1e-13,
	      double max_dyn = 1e6,
	      double min_viol = 1e-4,
	      int max_supp_abs = 1000,
	      double max_supp_rel = 0.1,
	      CleaningProcedure clean_proc = CP_CGLLANDP1,
	      bool use_int_slacks = false,
	      bool check_duplicates = false,
	      bool integral_scale_cont = false,
	      bool enforce_scaling = true);

   /// Constructor from CglParam
  CglGMIParam(CglParam &source,
	      double away = 0.005,
	      double eps_elim = 1e-12,
	      double eps_relax_abs = 1e-11,
	      double eps_relax_rel = 1e-13,
	      double max_dyn = 1e6,
	      double min_viol = 1e-4,
	      double max_supp_rel = 0.1,
	      CleaningProcedure clean_proc = CP_CGLLANDP1,
	      bool use_int_slacks = false,
	      bool check_duplicates = false,
	      bool integral_scale_cont = false,
	      bool enforce_scaling = true);
  
  /// Copy constructor 
  CglGMIParam(const CglGMIParam &source);

  /// Clone
  virtual CglGMIParam* clone() const;

  /// Assignment operator 
  virtual CglGMIParam& operator=(const CglGMIParam &rhs);

  /// Destructor 
  virtual ~CglGMIParam();
  //@}

protected:

  /**@name Parameters */
  //@{

  /** Use row only if pivot variable should be integer but is more 
      than AWAY from being integer. */
  double AWAY;

  /** Epsilon for value of coefficients when eliminating slack variables. 
      Default: 0. */
  double EPS_ELIM;

  /** Value added to the right hand side of each generated cut to relax it.
      Default: 1e-11 */
  double EPS_RELAX_ABS;

  /** For a generated cut with right hand side rhs_val, 
      EPS_RELAX_EPS * fabs(rhs_val) is used to relax the constraint.
      Default: 1.e-13 */
  double EPS_RELAX_REL;

  /** Maximum ratio between largest and smallest non zero 
      coefficients in a cut. Default: 1e6. */
  double MAXDYN;

  /** Minimum violation for the current basic solution in a generated cut.
      Default: 1e-4. */
  double MINVIOL;

  /** Maximum support relative to number of columns. Must be between 0
      and 1. Default: 0. */
  double MAX_SUPPORT_REL;

  /** Which cleaning procedure should be used? */
  CleaningProcedure CLEAN_PROC;

  /** Use integer slacks to generate cuts if USE_INTSLACKS = 1. Default: 0. */
  bool USE_INTSLACKS;

  /** Check for duplicates when adding the cut to the collection? */
  bool CHECK_DUPLICATES;

  /** Should we try to rescale cut coefficients on continuous variables
      so that they become integral, or do we only rescale coefficients
      on integral variables? Used only by cleaning procedure that try
      the integral scaling. */
  bool INTEGRAL_SCALE_CONT;

  /** Should we discard badly scaled cuts (according to the scaling
      procedure in use)? If false, CglGMI::scaleCut always returns
      true, even though it still scales cuts whenever possible, but
      not cut is rejected for scaling. Default true. Used only by
      cleaning procedure that try to scale. */
  bool ENFORCE_SCALING;

  //@}
};

#endif
