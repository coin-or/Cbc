# Copyright (C) 2006, 2007 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Cbc],[2.8],[cbc@lists.coin-or.org])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Eclipse Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/CbcTree.hpp)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

# If this project depends on external projects, the Externals file in
# the source root directory contains definition of where to find those
# externals.  The following macro ensures that those externals are
# retrieved by svn if they are not there yet.
AC_COIN_PROJECTDIR_INIT(Cbc)

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(Cbc)

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

#############################################################################
#                            COIN-OR components                             #
#############################################################################

AC_COIN_CHECK_PACKAGE(CoinDepend, [cgl osi coinutils], [CbcLib CbcGeneric])
if test $coin_has_coindepend != yes ; then
  AC_MSG_ERROR([One or more required packages CoinUtils, Osi, and Cgl are not available.])
fi

# Clp and OsiClp are inseparable here.

AC_COIN_CHECK_PACKAGE(Clp, [osi-clp], [CbcLib CbcGeneric])
if test $coin_has_clp != yes ; then
  AC_MSG_ERROR("Required package Clp unavailable.")
fi

# The rest are not required for linking against Cbc

AC_COIN_CHECK_PACKAGE(Sample,   [coindatasample])
AC_COIN_CHECK_PACKAGE(Miplib3,  [coindatamiplib3])

# Cbc really should make these tests; right now it assumes far too much about
# its environment. Needed for cbc-generic.

AC_COIN_CHECK_CXX_CHEADER(math)
AC_COIN_CHECK_CXX_CHEADER(float)
AC_COIN_CHECK_CXX_CHEADER(ieeefp)

AC_COIN_CHECK_GNU_READLINE(CbcLib CbcGeneric)

#############################################################################
#                    LP solvers other than CLP                              #
#############################################################################

# Check which other LP solvers are available, some of them become a dependency of CbcGeneric
AC_COIN_CHECK_PACKAGE(DyLP, [osi-dylp],   [CbcGeneric])
AC_COIN_CHECK_PACKAGE(Vol,  [osi-vol])
AC_COIN_CHECK_PACKAGE(Cpx,  [osi-cplex],  [CbcLib CbcGeneric])
AC_COIN_CHECK_PACKAGE(Glpk, [osi-glpk],   [CbcGeneric])
AC_COIN_CHECK_PACKAGE(Grb,  [osi-gurobi])
AC_COIN_CHECK_PACKAGE(Msk,  [osi-mosek],  [CbcGeneric])
AC_COIN_CHECK_PACKAGE(Spx,  [osi-soplex], [CbcGeneric])
AC_COIN_CHECK_PACKAGE(Xpr,  [osi-xpress])

#############################################################################
#                       Other third party software                          #
#############################################################################

# Ampl Solver library
AC_COIN_CHECK_PACKAGE(ASL, [coinasl], [CbcLib CbcGeneric])

#############################################################################
#                       CbcGeneric configuration			    #
#############################################################################

# Are we building cbc-generic? Default is no (false).

AC_ARG_WITH([cbc-generic],
  AS_HELP_STRING([--with-cbc-generic],
		 [specify whether to build cbc-generic (default: no)]),
  [case $withval in
    yes) cbc_with_cbc_generic=true ;;
      *) cbc_with_cbc_generic=false ;;
   esac],
  [cbc_with_cbc_generic=false])
if test x"$cbc_with_cbc_generic" = xyes ; then
  AC_MSG_NOTICE([building cbc-generic])
fi
AM_CONDITIONAL(CBC_BUILD_CBC_GENERIC,
	       test x"$cbc_with_cbc_generic" = xtrue)

# Set the default solver for cbc-generic. In theory, any OsiXXX should work.
# In practice, only the three listed below have had any testing.

AC_MSG_CHECKING(for cbc-generic default solver)
AC_ARG_WITH([cbc-generic-solver],
            AC_HELP_STRING([--with-cbc-generic-solver],
                           [specify default solver for cbc-generic in lower
			    case letters (if not given, clp is assumed)]),
            [cbc_default_solver=$withval],
            [cbc_default_solver=clp])
AC_MSG_RESULT($cbc_default_solver)

# FIXME what about cplex and mosek here?
case $cbc_default_solver in
    clp)
      ;;
   dylp)
      ;;
   glpk)
      ;;
 soplex)
      ;;
      *)
      AC_MSG_WARN([Unrecognised solver $cbc_default_solver chosen as cbc-generic default solver.])
      ;;
esac
AC_DEFINE_UNQUOTED([CBC_DEFAULT_SOLVER],"$cbc_default_solver",
  [Define to the name of the default solver to be used in cbc-generic in small letters])


#############################################################################
#                       Cbc parallel configuration			    #
#############################################################################

# First we define a new option, --enable-cbc-parallel

AC_ARG_ENABLE([cbc-parallel],
[AC_HELP_STRING([--enable-cbc-parallel],
                [enables compilation of the SMP version of Cbc])])

if test "$enable_cbc_parallel" = yes; then
  # Define the preprocessor macro
  AC_DEFINE([CBC_THREAD],[1],[Define to 1 if the SMP version of Cbc should be compiled])
  if test $coin_cxx_is_cl = true ;
  then
    # TODO we should check whether the library works and pthread.h is indeed there
    AC_ARG_WITH(pthreadsw32-lib, 
                AC_HELP_STRING([--with-pthreadsw32-lib],
                               [specify the name of the pthreads win32 library (for cl/icl builds only)]),
                [CBCLIB_LIBS="$withval $CBCLIB_LIBS"
                 CBCLIB_PCLIBS="$withval $CBCLIB_PCLIBS"],
                [AC_MSG_WARN([--enable-cbc-parallel selected, but --with-pthreadsw32-lib not given])
                 CBCLIB_LIBS="pthreadVCE2.lib $CBCLIB_LIBS"])
    
    AC_ARG_WITH(pthreadsw32-incdir,
                AC_HELP_STRING([--with-pthreadsw32-incdir],
                               [specify the path of the pthreads win32 header file (for cl/icl builds only)]),
                [CXXFLAGS="$CXXFLAGS -I$withval"],
                [AC_MSG_WARN([--enable-cbc-parallel selected, but --with-pthreadsw32-incdir not given])])
  else
    AC_CHECK_LIB([rt],[clock_gettime],
                 [CBCLIB_LIBS="-lrt $CBCLIB_LIBS"
                  CBCLIB_PCLIBS="-lrt $CBCLIB_PCLIBS"
                  AC_DEFINE([HAVE_CLOCK_GETTIME],[1],[Define if clock_gettime and rt library is available])
                 ],
                 [AC_MSG_WARN([--enable-cbc-parallel selected, but -lrt unavailable; using gettimeofday instead of clock_gettime])])
    AC_CHECK_LIB([pthread],[pthread_create],
                 [CBCLIB_LIBS="-lpthread $CBCLIB_LIBS"
                  CBCLIB_PCLIBS="-lpthread $CBCLIB_PCLIBS"
                 ],
                 [AC_MSG_ERROR([--enable-cbc-parallel selected, but -lpthreads unavailable])])
  fi
  AC_MSG_NOTICE([Cbc multithreading enabled]);
fi

##############################################################################
#                   VPATH links for example input files                      #
##############################################################################

# In case this is a VPATH configuration we need to make sure that the
# input files for the examples are available in the VPATH directory.

AC_COIN_VPATH_LINK(examples/sudoku_sample.csv)
AC_COIN_VPATH_LINK(examples/quad.mps)
AC_COIN_VPATH_LINK(examples/quad2.mps)

#############################################################################
#                  Check for doxygen                                        #
#############################################################################

AC_COIN_DOXYGEN(CoinUtils Osi Clp DyLP Cgl Vol)

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 examples/Makefile
                 src/Makefile
                 test/Makefile
		 cbc.pc cbc-uninstalled.pc])
AC_CONFIG_FILES([doxydoc/doxygen.conf])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([src/config.h src/config_cbc.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
