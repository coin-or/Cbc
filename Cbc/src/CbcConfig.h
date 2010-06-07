/* $Id$ */
/*
 * Include file for the configuration of Cbc.
 *
 * On systems where the code is configured with the configure script
 * (i.e., compilation is always done with HAVE_CONFIG_H defined), this
 * header file includes the automatically generated header file, and
 * undefines macros that might configure with other Config.h files.
 *
 * On systems that are compiled in other ways (e.g., with the
 * Developer Studio), a header files is included to define those
 * macros that depend on the operating system and the compiler.  The
 * macros that define the configuration of the particular user setting
 * (e.g., presence of other COIN packages or third party code) are set
 * here.  The project maintainer needs to remember to update this file
 * and choose reasonable defines.  A user can modify the default
 * setting by editing this file here.
 *
 */

#ifndef __CBCCONFIG_H__
#define __CBCCONFIG_H__

#ifdef HAVE_CONFIG_H
#include "config_cbc.h"

/* undefine macros that could conflict with those in other config.h
   files */
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION

#else /* HAVE_CONFIG_H */

/* include the COIN-wide system specific configure header */
#include "configall_system.h"

/***************************************************************************/
/*             HERE DEFINE THE CONFIGURATION SPECIFIC MACROS               */
/***************************************************************************/

/* Define to the debug sanity check level (0 is no test) */
#define COIN_CBC_CHECKLEVEL 0

/* Define to the debug verbosity level (0 is no output) */
#define COIN_CBC_VERBOSITY 0

/* Define to 1 if the Cbc package is used */
#define COIN_HAS_CBC 1

/* Define to 1 if the Cgl package is used */
#define COIN_HAS_CGL 1

/* Define to 1 if the Clp package is used */
#define COIN_HAS_CLP 1

/* Define to 1 if the CoinUtils package is used */
#define COIN_HAS_COINUTILS 1

/* Define to 1 if the Osi package is used */
#define COIN_HAS_OSI 1

/* Define to 1 if the Vol package is used */
#define COIN_HAS_VOL 1

/* Define to 1 if the Cplex package is used */
/* #undef COIN_HAS_CPX */

/* Define to 1 if the Dylp package is used */
/* #undef COIN_HAS_DYLP */

/* Define to 1 if the FortMP package is used */
/* #undef COIN_HAS_FMP */

/* Define to 1 if the Glpk package is used */
/* #undef COIN_HAS_GLPK */

/* Define to 1 if the Mosek package is used */
/* #undef COIN_HAS_MSK */

/* Define to 1 if the Osl package is used */
/* #undef COIN_HAS_OSL */

/* Define to 1 if the Soplex package is used */
/* #undef COIN_HAS_SPX */

/* Define to 1 if the Sym package is used */
/* #undef COIN_HAS_SYM */

/* Define to 1 if the Xpress package is used */
/* #undef COIN_HAS_XPR */

#define CBC_VERSION "trunk"

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

#endif /* HAVE_CONFIG_H */

#endif /* __CBCCONFIG_H__ */
