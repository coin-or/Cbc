
/* include the COIN-OR-wide system specific configure header */
#include "configall_system.h"

/* this needs to come before the include of config_cbc_default.h */
#ifndef CBCLIB_EXPORT
#if defined(_WIN32) && defined(DLL_EXPORT)
#define CBCLIB_EXPORT __declspec(dllexport)
#else
#define CBCLIB_EXPORT
#endif
#endif

/* include the public project specific macros */
#include "config_cbc_default.h"

/***************************************************************************/
/*             HERE DEFINE THE PROJECT SPECIFIC MACROS                     */
/*    These are only in effect in a setting that doesn't use configure     */
/***************************************************************************/

/* Define to 1 if the CoinUtils package is used */
#define CBC_HAS_COINUTILS 1

/* Define to 1 if the Osi package is used */
#define CBC_HAS_OSI 1

/* Define to 1 if the Clp package is used */
#define CBC_HAS_CLP 1

/* Define to 1 if the Clp package is used */
#define CBC_HAS_OSICLP 1

/* Define to 1 if the Cgl package is used */
#define COIN_HAS_CGL 1

/* Define to 1 if the Vol package is used */
/* #define COIN_HAS_VOL 1 */

/* Define to 1 if the Cplex package is used */
/* #undef COIN_HAS_CPX */

/* Define to 1 if the Dylp package is used */
/* #undef COIN_HAS_DYLP */

/* Define to 1 if the Mosek package is used */
/* #undef COIN_HAS_MSK */

/* Define to 1 if the Soplex package is used */
/* #undef COIN_HAS_SPX */

/* Define to 1 if the Sym package is used */
/* #undef COIN_HAS_SYM */

/* Define to 1 if the Xpress package is used */
/* #undef COIN_HAS_XPR */

/*
  For additional information about how to set OSICBC_DFLT_SOLVER,
  OSICBC_DFLT_SOLVER_CLP, and OSICBC_DFLT_SOLVER_HPP, please see comments at
  the beginning of OsiCbcSolverInterface.cpp. Unless you know what you're
  doing, you should use clp with OsiCbc. Just uncomment the next three
  defines.
*/
/*
  Define to the name of the default solver interface class, e.g.,
  OsiClpSolverInterface.
*/
/* #define OSICBC_DFLT_SOLVER OsiClpSolverInterface */

/* Define this symbol if clp is the default solver. */
/* #define OSICBC_DFLT_SOLVER_CLP 1 */

/*
  Define to the name of the .hpp file for the default solver interface class,
  e.g., "OsiClpSolverInterface.hpp" (include quotes)
*/
/* #define OSICBC_DFLT_SOLVER_HPP "OsiClpSolverInterface.hpp" */

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
