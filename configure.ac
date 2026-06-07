# Copyright (C) 2006, 2007 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_INIT([Cbc],[0.2.3],[https://github.com/h-g-s/mipster/issues/new],[mipster],[https://github.com/h-g-s/mipster])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Eclipse Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there.
AC_CONFIG_SRCDIR(src/CbcTree.hpp)

# Do some initialization work (version numbers, change prefix default, ...)
AC_CONFIG_MACRO_DIRS([BuildTools])
AC_COIN_INITIALIZE

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Initialize libtool
AC_COIN_PROG_LIBTOOL

# set RPATH_FLAGS to the compiler link flags required to hardcode location
# of the shared objects (expanded_libdir is set somewhere in configure before)
# (use in examples Makefile)
AC_COIN_RPATH_FLAGS([$expanded_libdir])

#############################################################################
#                   Math headers                                            #
#############################################################################

# Check for cmath/math.h, cfloat/float.h, cieeefp/ieeefp.h
AC_COIN_CHECK_MATH_HDRS

#############################################################################
#                   Debug info                                              #
#############################################################################

AC_CHECK_HEADERS([execinfo.h])
AC_CHECK_HEADERS([signal.h])
AC_CHECK_HEADERS([immintrin.h])

#############################################################################
#                            COIN-OR components                             #
#############################################################################

AC_DEFINE([CBC_HAS_COINUTILS],[1],[Define to 1 because CoinUtils is merged into Cbc.])
AC_COIN_CHK_LAPACK(CbcLib)
AC_COIN_CHK_LIBM(CbcLib)
AC_COIN_CHK_ZLIB(CbcLib)
AC_COIN_CHK_BZLIB(CbcLib)

# AMD (Approximate Minimum Degree ordering) from SuiteSparse — BSD 3-clause license.
# Used by merged Clp sources for the barrier (interior-point) LP solver.
AC_COIN_CHK_LIBHDR(AMD,[CbcLib],[-lamd],[-I/usr/include/suitesparse],[],
  [[double Control[AMD_CONTROL];
   amd_defaults(Control);]],
  [#include "amd.h"])

# Optional package for command line edit support

#############################################################################
#                       Other third party software                          #
#############################################################################

# Ampl Solver library
AC_COIN_CHK_PKG(ASL,[CbcLib],[coinasl])

# Nauty library (for symmetry detection).

#############################################################################

# The rest are not required for linking against Cbc but are used for the unit
# test. For the data packages, we just want XXX_DATA for use in Makefile.am.
AC_COIN_CHK_PKG(OsiTests,[CbcUnitTest],[osi-unittests])
AC_COIN_CHK_PKG(Sample,,[coindatasample],,dataonly)
AC_COIN_CHK_PKG(Netlib,,[coindatanetlib],,dataonly)
AC_COIN_CHK_PKG(Miplib3,,[coindatamiplib3],,dataonly)

#############################################################################
#                    LP solvers other than CLP                              #
#############################################################################

# Check which other LP solvers are available, some of them become a dependency
# of CbcGeneric
#AC_COIN_CHK_PKG(OsiDyLP,[CbcGeneric OsiCbcLib],[osi-dylp])
#AC_COIN_CHK_PKG(OsiCpx,[CbcLib CbcGeneric OsiCbcLib],[osi-cplex])
#AC_COIN_CHK_PKG(OsiGlpk,[CbcGeneric OsiCbcLib],[osi-glpk])
##AC_COIN_CHK_PKG(OsiGrb,,[osi-gurobi])
#AC_COIN_CHK_PKG(OsiMsk,[CbcGeneric OsiCbcLib],[osi-mosek])
#AC_COIN_CHK_PKG(OsiSpx,[CbcGeneric OsiCbcLib],[osi-soplex])
#AC_COIN_CHK_PKG(OsiXpr,,[osi-xpress])

AC_COIN_FINALIZE_FLAGS([CbcLib CbcUnitTest])

#############################################################################
#                       CbcGeneric configuration                            #
#############################################################################

# Are we building cbc-generic? Default is no (false).
dnl AC_ARG_ENABLE([cbc-generic],
dnl   AS_HELP_STRING([--enable-cbc-generic],[enables build of cbc-generic]),
dnl   [cbc_generic=$enableval],
dnl   [cbc_generic=no])
dnl AM_CONDITIONAL(CBC_BUILD_CBC_GENERIC, test "$cbc_generic" = yes)

dnl if test "$cbc_generic" = yes ; then
dnl   # Set the default solver for cbc-generic. In theory, any OsiXXX should work.
dnl   # In practice, only thoselisted below have had any testing long long ago or never.
dnl   AC_MSG_CHECKING(for cbc-generic default solver)
dnl   AC_ARG_WITH([cbc-generic-solver],
dnl               AS_HELP_STRING([--with-cbc-generic-solver],[specify default solver for cbc-generic in lower case letters (if not given, clp is assumed)]),
dnl               [cbc_default_solver=$withval],
dnl               [cbc_default_solver=clp])
dnl   AC_MSG_RESULT($cbc_default_solver)

dnl   case $cbc_default_solver in
dnl       clp) ;;
dnl     cplex) ;;
dnl      dylp) ;;
dnl      glpk) ;;
dnl     highs) ;;
dnl     mosek) ;;
dnl    soplex) ;;
dnl         *) AC_MSG_WARN([Unrecognised solver $cbc_default_solver chosen as cbc-generic default solver.])
dnl           ;;
dnl   esac
dnl   AC_DEFINE_UNQUOTED([CBC_DEFAULT_SOLVER],"$cbc_default_solver",
dnl     [Define to the name of the default solver to be used in cbc-generic in small letters])

dnl   AC_MSG_NOTICE([building cbc-generic with default solver $cbc_default_solver])
dnl fi

#############################################################################
#                       Cbc parallel configuration	                         #
#############################################################################

# First we define a new option, --enable-cbc-parallel
AC_ARG_ENABLE([cbc-parallel],
  [AS_HELP_STRING([--disable-cbc-parallel],[disables compilation of the SMP version of Cbc])],
  [enable_cbc_parallel=$enableval],
  [enable_cbc_parallel=auto])

if test "$enable_cbc_parallel" != no; then
  # Define the preprocessor macro
  case "$CXX" in
  *cl) # Windows / cl: look for pthreadsw32 lib
    AC_COIN_CHK_LIB(pthreadsw32,[CbcLib CbcSolverLib],[pthreadVCE2.lib],[],[],[pthread_create])
    if test "$coin_has_pthreadsw32" = yes ; then
      enable_cbc_parallel=yes
    elif test "$enable_cbc_parallel" = yes ; then
      AC_MSG_ERROR([--enable-cbc-parallel selected, but pthreads-w32 unavailable])
    else
      enable_cbc_parallel=no
    fi
    ;;
  *) # not Windows / cl: look for pthread lib
    AC_CHECK_LIB([pthread],[pthread_create],
      [CBCLIB_LFLAGS="-lpthread $CBCLIB_LFLAGS"
       enable_cbc_parallel=yes],
      [if test "$enable_cbc_parallel" = yes ; then
         AC_MSG_ERROR([--enable-cbc-parallel selected, but -lpthreads unavailable])
       else
         enable_cbc_parallel=no
       fi])
    if test "$enable_cbc_parallel" = yes ; then
      AC_CHECK_LIB([rt],[clock_gettime],
        [CBCLIB_LFLAGS="-lrt $CBCLIB_LFLAGS"
         AC_DEFINE([HAVE_CLOCK_GETTIME],[1],[Define if clock_gettime and rt library is available])],
        [AC_MSG_WARN([SMP version of Cbc selected, but -lrt unavailable; using gettimeofday instead of clock_gettime])])
    fi
    AM_CONDITIONAL(COIN_HAS_PTHREADSW32, false)
  ;;
  esac
else
  AM_CONDITIONAL(COIN_HAS_PTHREADSW32, false)
fi
if test "$enable_cbc_parallel" = yes ; then
  AC_DEFINE([CBC_THREAD],[1],[Define to 1 if the SMP version of Cbc should be compiled])
  AC_MSG_NOTICE([Cbc multithreading enabled]);
fi

##############################################################################
#                   VPATH links for example input files                      #
##############################################################################

# In case this is a VPATH configuration we need to make sure that the
# input files for the examples are available in the VPATH directory.
AC_COIN_VPATH_LINK(examples/sudoku_sample.csv)

#############################################################################
#                  Check for doxygen                                        #
#############################################################################

AC_COIN_DOXYGEN([CoinUtils])

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 examples/Makefile
                 src/CoinUtils/Makefile
                 src/Clp/Makefile
                 src/Cgl/Makefile
                 src/Makefile
                 test/Makefile
                 mipster.pc])
AC_CONFIG_FILES([doxydoc/doxygen.conf])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADERS([src/config.h src/config_cbc.h])

AC_COIN_FINALIZE_FLAGS([CbcLib CbcSolverLib CbcGeneric CbcUnitTest])

# Capture git hash for development builds (records state at configure time)
AC_MSG_CHECKING([git revision])
CBC_GIT_HASH=`cd $srcdir && (git describe --always --dirty 2>/dev/null || echo "unknown")`
AC_MSG_RESULT([$CBC_GIT_HASH])
AC_DEFINE_UNQUOTED([CBC_GIT_HASH], ["$CBC_GIT_HASH"], [Git hash at configure time])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
