# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Cbc],[0.1],[coin-cbc@list.coin-or.org])

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
AC_COIN_HAS_PROJECT(Sym)
AC_COIN_HAS_PROJECT(Vol)

# ToDo: find out when CBC_HAS_CLP should be defined!  For now, I
# replaced all CBC_USE by COIN_HAS #ifdef's

# Decide which one is the default solver
AC_MSG_CHECKING(for CBC default solver)
AC_ARG_WITH([cbc-default-solver],
            AC_HELP_STRING([--with-cbc-default-solver],
                           [specify default solver for Cbc in lower case letters (if not given, clp is assumed)]),
            [cbc_default_solver=$withval],
            [cbc_default_solver=clp])
AC_MSG_RESULT($cbc_default_solver)

case $cbc_default_solver in
  clp)
    AC_DEFINE([CBC_DEFAULT_SOLVER],"clp",
              [Define to the name of the default solver to be used in CBC in small letters])
    ;;
  dylp)
    AC_DEFINE([CBC_DEFAULT_SOLVER],"dylp")
    ;;
  *)
    AC_MSG_ERROR([Unknown solver $cbc_default_solver chosen as CBC's default solver])
    ;;
esac

# Steal some files from Clp - they are going to be deleted by a
# distclean in src (see the Makefile.am there)
if test x"$coin_has_clp" != xunavailable; then
  AC_MSG_NOTICE([Stealing CbcOrClpParam's files from Clp])
  AC_CONFIG_LINKS([src/CbcOrClpParam.cpp:../Clp/src/CbcOrClpParam.cpp
                   src/CbcOrClpParam.hpp:../Clp/src/CbcOrClpParam.hpp])
fi

#############################################################################
#                            Third party solvers                            #
#############################################################################

# Check which third party solvers are available ToDo!
AC_COIN_HAS_USER_LIBRARY([Cplex],[CPX],[cplex.h],[CPXgetstat])
AC_COIN_HAS_USER_LIBRARY([Glpk],[GLPK],[glpk.h],[glp_lpx_simplex])
AC_COIN_HAS_USER_LIBRARY([FortMP],[FMP])  #ToDo! is Fortran?
AC_COIN_HAS_USER_LIBRARY([Mosek],[MSK],[mosekdl.h],[MSK_openmosek])
AC_COIN_HAS_USER_LIBRARY([Osl],[OSL],[ekk_c_api.h],[ekk_simplex])
AC_COIN_HAS_USER_LIBRARY([Soplex],[SPX],[spxsolver.h]) # ToDo
AC_COIN_HAS_USER_LIBRARY([Xpress],[XPR],[xprs.h],[XPRSsetintcontrol])

#############################################################################
#                       Other third party software                          #
#############################################################################

# Ampl Solver library
AC_COIN_HAS_ASL

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
