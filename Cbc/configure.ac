# Copyright (C) 2006, 2007 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Cbc],[2.3.2],[coin-cbc@list.coin-or.org])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Common Public License.])

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
AC_COIN_PROJECTDIR_INIT

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(Cbc)

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

#############################################################################
#                              COIN components                              #
#############################################################################

AC_COIN_HAS_PROJECT(Cbc)
AC_COIN_HAS_PROJECT(Cgl)
AC_COIN_HAS_PROJECT(Clp)
AC_COIN_HAS_PROJECT(CoinUtils)
AC_COIN_HAS_PROJECT(DyLP)
AC_COIN_HAS_PROJECT(Osi)
AC_COIN_HAS_PROJECT(Vol)
AC_COIN_HAS_PROJECT(SYMPHONY)
AC_COIN_HAS_PROJECT(Cpx)

#############################################################################
#                            Third party solvers                            #
#############################################################################

# Check which third party solvers are available ToDo!
AC_COIN_HAS_USER_LIBRARY([Cplex],[CPX],[cplex.h],[CPXgetstat])
AC_COIN_HAS_USER_LIBRARY([FortMP],[FMP])  #ToDo! is Fortran?
AC_COIN_HAS_USER_LIBRARY([Mosek],[MSK],[mosek.h],[MSK_openmosek])
AC_COIN_HAS_USER_LIBRARY([Osl],[OSL],[ekk_c_api.h],[ekk_simplex])
AC_COIN_HAS_USER_LIBRARY([Soplex],[SPX],[spxsolver.h]) # ToDo
AC_COIN_HAS_USER_LIBRARY([Xpress],[XPR],[xprs.h],[XPRSsetintcontrol])

AC_COIN_HAS_GLPK(Glpk)

#############################################################################
#                       Other third party software                          #
#############################################################################

# Ampl Solver library
AC_COIN_HAS_ASL

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

case $cbc_default_solver in
   clp)
    ;;
  dylp)
    ;;
  glpk)
    ;;
     *)
    AC_MSG_WARN([Unrecognised solver $cbc_default_solver chosen as cbc-generic default solver.])
    ;;
esac
AC_DEFINE_UNQUOTED([CBC_DEFAULT_SOLVER],"$cbc_default_solver",
  [Define to the name of the default solver to be used in cbc-generic in small letters])
AC_DEFINE_UNQUOTED(CBCVERSION, "$PACKAGE_VERSION", "Version of Cbc")


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
  AC_CHECK_LIB([rt],[clock_gettime],
               [ADDLIBS="-lrt $ADDLIBS"
                AC_DEFINE([HAVE_CLOCK_GETTIME],[1],[Define if clock_gettime and rt library is available])],
               [AC_MSG_WARN([--enable-cbc-parallel selected, but -lrt unavailable; using gettimeofday instead of clock_gettime])])
  AC_CHECK_LIB([pthread],[pthread_create],
               [ADDLIBS="-lpthread $ADDLIBS"],
	       [AC_MSG_ERROR([--enable-cbc-parallel selected, but -lpthreads unavailable])])
fi

##############################################################################
#                   VPATH links for example input files                      #
##############################################################################

# In case this is a VPATH configuration we need to make sure that the
# input files for the examples are available in the VPATH directory.

AC_COIN_VPATH_LINK(examples/sudoku_sample.csv)
AC_COIN_VPATH_LINK(examples/quad.mps)
AC_COIN_VPATH_LINK(examples/quad2.mps)

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 examples/Makefile
                 src/Makefile
                 test/Makefile
		 cbc_addlibs.txt])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([inc/config_cbc.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
