/* This is the header file for the Microsoft compiler, defining all
 * system and compiler dependent configuration macros */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

#ifndef COIN_USE_F2C
/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
# define F77_FUNC(name,NAME) NAME

/* As F77_FUNC, but for C identifiers containing underscores. */
# define F77_FUNC_(name,NAME) NAME
#else
/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
# define F77_FUNC(name,NAME) name ## _

/* As F77_FUNC, but for C identifiers containing underscores. */
# define F77_FUNC_(name,NAME) name ## __
#endif

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to 1 if you have the <assert.h> header file. */
/* #undef HAVE_ASSERT_H */

/* Define to 1 if you have the <cassert> header file. */
#define HAVE_CASSERT 1

/* Define to 1 if you have the <cctype> header file. */
#define HAVE_CCTYPE 1

/* Define to 1 if you have the <cfloat> header file. */
#define HAVE_CFLOAT 1

/* Define to 1 if you have the <cieeefp> header file. */
/* #undef HAVE_CIEEEFP */

/* Define to 1 if you have the <cmath> header file. */
#define HAVE_CMATH 1

/* Define to 1 if you have the <cstdarg> header file. */
#define HAVE_CSTDARG 1

/* Define to 1 if you have the <cstdio> header file. */
#define HAVE_CSTDIO 1

/* Define to 1 if you have the <cstdlib> header file. */
#define HAVE_CSTDLIB 1

/* Define to 1 if you have the <cstring> header file. */
#define HAVE_CSTRING 1

/* Define to 1 if you have the <ctime> header file. */
#define HAVE_CTIME 1

/* Define to 1 if you have the <cstddef> header file. */
#define HAVE_CSTDDEF 1

/* Define to 1 if you have the <ctype.h> header file. */
/* #undef HAVE_CTYPE_H */

/* Define to 1 if you have the <dlfcn.h> header file. */
/* #undef HAVE_DLFCN_H */

/* Define to 1 if function drand48 is available */
/* #undef HAVE_DRAND48 */

/* Define to 1 if you have the <float.h> header file. */
/* #undef HAVE_FLOAT_H */

/* Define to 1 if you have the <ieeefp.h> header file. */
/* #undef HAVE_IEEEFP_H */

/* Define to 1 if you have the <inttypes.h> header file. */
/* #define HAVE_INTTYPES_H */

/* Define to 1 if you have the <math.h> header file. */
/* #undef HAVE_MATH_H */

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if function rand is available */
#define HAVE_RAND 1

/* Define to 1 if you have the <stdarg.h> header file. */
/* #undef HAVE_STDARG_H */

/* Define to 1 if you have the <stdint.h> header file. */
/* #undef HAVE_STDINT_H */

/* Define to 1 if you have the <stdio.h> header file. */
/* #undef HAVE_STDIO_H */

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if function std::rand is available */
#define HAVE_STD__RAND 1

/* Define to 1 if you have the <strings.h> header file. */
/* #define HAVE_STRINGS_H */

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <time.h> header file. */
/* #undef HAVE_TIME_H */

/* Define to 1 if you have the <unistd.h> header file. */
/* #define HAVE_UNISTD_H */

/* Define to 1 if va_copy is avaliable */
/* #undef HAVE_VA_COPY */

/* Define to 1 if you have the <windows.h> header file. */
/* #undef HAVE_WINDOWS_H */

/* Define to 1 if you have the `_snprintf' function. */
#define HAVE__SNPRINTF 1

/* The size of a `double', as computed by sizeof. */
#define SIZEOF_DOUBLE 8

/* The size of a `int', as computed by sizeof. */
#define SIZEOF_INT 4

/* The size of a `int *', as computed by sizeof. */
#define SIZEOF_INT_P 8

/* The size of a `long', as computed by sizeof. */
#define SIZEOF_LONG 4

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1
