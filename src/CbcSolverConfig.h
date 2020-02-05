/* Copyright (C) 2011
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * $Id: CbcConfig.h 2583 2019-06-02 23:23:54Z lou $
 *
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
 * (e.g., presence of other COIN-OR packages or third party code) are set
 * by the files config_*default.h. The project maintainer needs to remember
 * to update these file and choose reasonable defines.
 * A user can modify the default setting by editing the config_*default.h files.
 *
 */

#ifndef __CBCSOLVERCONFIG_H__
#define __CBCSOLVERCONFIG_H__

#ifdef HAVE_CONFIG_H
#ifdef CBCSOLVER_BUILD
#include "config.h"

/* overwrite CBCSOLVERLIB_EXPORT from config.h when building Cbc
 * we want it to be __declspec(dllexport) when building a DLL on Windows
 * we want it to be __attribute__((__visibility__("default"))) when building with GCC,
 *   so user can compile with -fvisibility=hidden
 */
#ifdef DLL_EXPORT
#undef CBCSOLVERLIB_EXPORT
#define CBCSOLVERLIB_EXPORT __declspec(dllexport)
#elif defined(__GNUC__) && __GNUC__ >= 4
#undef CBCSOLVERLIB_EXPORT
#define CBCSOLVERLIB_EXPORT __attribute__((__visibility__("default")))
#endif

#else
#include "config_cbcsolver.h"
#endif

#else /* HAVE_CONFIG_H */

#ifdef CBCSOLVER_BUILD
#include "config_default.h"
#else
#include "config_cbc_default.h"
#endif

#ifndef CBCSOLVERLIB_EXPORT
# if defined(_WIN32) && defined(DLL_EXPORT)
#  ifdef CBCSOLVER_BUILD
#   define CBCSOLVERLIB_EXPORT __declspec(dllexport)
#  else
#   define CBCSOLVERLIB_EXPORT __declspec(dllimport)
#  endif
# elif defined(__GNUC__) && __GNUC__ >= 4
#  define CBCSOLVERLIB_EXPORT __attribute__((__visibility__("default")))
# endif
#endif


#endif /* HAVE_CONFIG_H */

#endif /*__CBCSOLVERCONFIG_H__*/

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
