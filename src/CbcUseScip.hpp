// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcUseScip_H
#define CbcUseScip_H
#if CBC_TRY_SCIP
// I have put Scip "interface" here to encourage modifying it
// It would be good if someone sensible could clean the whole thing up
/* to get CBC_TRY_SCIP > 1 to work i.e. presolve
   I do not really understand Scip code so used modified code which
   writes out .lp file.
   On restoring I was unable to free ALL Scip memory so gave up
   - so fine for one run but horrific memory leak
   SCIP_EXPORT had to be added in some files
added misc.h
+   SCIP_EXPORT
  SCIP_RETCODE SCIPrealarrayExtend(
     SCIP_REALARRAY*       realarray, 
+   SCIP_EXPORT
  SCIP_RETCODE SCIPboolarrayExtend(
     SCIP_BOOLARRAY*       boolarray, 
added prob.h
+  SCIP_EXPORT
  void SCIPprobMarkNConss(
+ SCIP_EXPORT
  void SCIPprobAddObjoffset(
+   SCIP_EXPORT
  SCIP_RETCODE SCIPprobInitSolve(
added sol.h
+   SCIP_EXPORT
  SCIP_RETCODE SCIPsolCreate(
+   SCIP_EXPORT
  SCIP_RETCODE SCIPsolRetransform(
added var.h
+   SCIP_EXPORT
  SCIP_RETCODE SCIPvarGetProbvarSum(
 */
#include <scip/scip.h>
#include <blockmemshell/memory.h>
#include "scip/reader_nl.h"
#include <scip/type_cons.h>
#include <scip/primal.h>
#include <scip/type_primal.h>
#include "scip/cons_linear.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"
#include "scip/cons.h"
#include "scip/struct_scip.h"
#include "scip/struct_cons.h"
#include "scip/struct_var.h"
#include "scip/struct_mem.h"
#include "scip/struct_sol.h"
#include "scip/struct_misc.h"
#include "scip/scip_numerics.h"
#include <scip/scip_prob.h>
#include "scip/set.h"
#include "scip/dcmp.h"
#include "scip/cutpool.h"
#include "scip/struct_reader.h"
#include "scip/reader.h"
#include "scip/tree.h"
#include "scip/solve.h"
#include "scip/sol.h"
#include "scip/scip_reader.h"
#include <scip/prob.h>
#include <scip/type_prob.h>
#include <scip/type_heur.h>
#include <scip/type_cons.h>
#include <scip/pricestore.h>
#include <scip/sepastore.h>
#include <scip/type_var.h>
#include <scip/var.h>
#include <scip/misc.h>
#include <scip/lp.h>
#include <scip/heur_trysol.h>
#include <scip/scipdefplugins.h>
/** constraint data for AND-constraints */
struct SCIP_ConsData
{
   SCIP_VAR**            vars;               /**< variables in the AND-constraint */
   SCIP_VAR*             resvar;             /**< resultant variable */
   SCIP_ROW**            rows;               /**< rows for linear relaxation of AND-constraint */
   SCIP_ROW*             aggrrow;            /**< aggregated row for linear relaxation of AND-constraint */
   SCIP_NLROW*           nlrow;              /**< row for representation in nonlinear relaxation */
   int                   nvars;              /**< number of variables in AND-constraint */
   int                   varssize;           /**< size of vars array */
   int                   nrows;              /**< number of rows for linear relaxation of AND-constraint */
   int                   watchedvar1;        /**< position of first watched operator variable */
   int                   watchedvar2;        /**< position of second watched operator variable */
   int                   filterpos1;         /**< event filter position of first watched operator variable */
   int                   filterpos2;         /**< event filter position of second watched operator variable */
   unsigned int          propagated:1;       /**< is constraint already preprocessed/propagated? */
   unsigned int          nofixedzero:1;      /**< is none of the operator variables fixed to FALSE? */
   unsigned int          impladded:1;        /**< were the implications of the constraint already added? */
   unsigned int          opimpladded:1;      /**< was the implication for 2 operands with fixed resultant added? */
   unsigned int          sorted:1;           /**< are the constraint's variables sorted? */
   unsigned int          changed:1;          /**< was constraint changed since last pair preprocessing round? */
   unsigned int          merged:1;           /**< are the constraint's equal variables already merged? */
   unsigned int          checkwhenupgr:1;    /**< if AND-constraint is upgraded to an logicor constraint or the and-
                                              *   constraint is linearized, should the check flag be set to true, even
                                              *   if the AND-constraint has a check flag set to false? */
   unsigned int          notremovablewhenupgr:1;/**< if AND-constraint is upgraded to an logicor constraint or the and-
                                              *   constraint is linearized, should the removable flag be set to false,
                                              *   even if the AND-constraint has a removable flag set to true? */
};
#endif

#include "CbcModel.hpp"
/** 
    This is for playing with Scip in any way
*/

class CBCLIB_EXPORT CbcUseScip {

public:
  // Default Constructor
  CbcUseScip();

  /** Useful constructor
      Just needs to point to model.
      Later may do more - so that is why it inherits
    */
  CbcUseScip(CbcModel *model);

  // Copy constructor
  CbcUseScip(const CbcUseScip &);

  // Assignment operator
  CbcUseScip &operator=(const CbcUseScip &rhs);

  // Destructor
  ~CbcUseScip();

  /// Solver
  inline OsiSolverInterface * modifiedSolver()
  { return modSolver_;}

  /** type -
      0 - use continuous solver
      1 - use current solver
  */
  int tryScip(int type);

  /// Create presolved model
  int presolveModel(OsiSolverInterface * solver);

  /// postprocess model
  // saveOriginalSolver has been saved
  int afterSolve(const double * solution,double objValue);

protected:
  /// data
  /// Cbc model
  CbcModel * model_;
  /// Osi model
  OsiSolverInterface * originalSolver_;
  /// Osi model
  OsiSolverInterface * modSolver_;
  /// Scip stuff
  SCIP * scip_;
};

#endif
