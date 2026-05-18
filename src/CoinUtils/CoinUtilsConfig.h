/* src/config_coinutils.h.in.  Generated from configure.ac by autoheader.  */
#ifndef COINUTILSCONFIG_H
#define COINUTILSCONFIG_H

#define COINUTILSLIB_EXPORT

#define COINUTILS_VERSION "devel"
#define COINUTILS_VERSION_MAJOR 9999
#define COINUTILS_VERSION_MINOR 9999
#define COINUTILS_VERSION_RELEASE 9999

#define COINUTILS_C_FINITE std::isfinite
#define COINUTILS_C_ISNAN std::isnan

#define COINUTILS_HAS_ZLIB 1
#define COINUTILS_HAS_LAPACK 1
#define COINUTILS_LAPACK_FUNC(name,NAME) name ## _
#define COINUTILS_LAPACK_FUNC_(name,NAME) name ## _

#define COINUTILS_INT64_T long long
#define COINUTILS_INTPTR_T intptr_t
#define COINUTILS_UINT64_T unsigned long long

#define COINUTILS_BIGINDEX_T int
#define COINUTILS_BIGINDEX_IS_INT 1

#define COINUTILS_HAS_STDINT_H 1
#define COINUTILS_HAS_CSTDINT 1
#define COINUTILS_CPLUSPLUS11 1
#define COINUTILS_CHECKLEVEL 0
#define COINUTILS_VERBOSITY 0
#define COINUTILS_MEMPOOL_MAXPOOLED -1

#endif
