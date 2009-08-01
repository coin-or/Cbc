# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Cbc],[2.3stable],[coin-cbc@lists.coin-or.org])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Common Public License.])

AC_REVISION(0.9)

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(configure.ac)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#  Do the tests necessary to configure compilers and initialise autotools   #
#############################################################################

AC_COIN_CREATE_LIBTOOL

#############################################################################
#             Check which ThirdParty sources are available                  #
#############################################################################

#AC_COIN_THIRDPARTY_SUBDIRS([ASL],[ThirdParty/ASL],[solvers/asl.h])
#AC_COIN_THIRDPARTY_SUBDIRS([Blas],[ThirdParty/Blas],[daxpy.f])
#AC_COIN_THIRDPARTY_SUBDIRS([Lapack],[ThirdParty/Lapack],[LAPACK/SRC/dlarf.f])
#AC_COIN_THIRDPARTY_SUBDIRS([Glpk],[ThirdParty/Glpk],[glpk/src/glplpx01.c])

#############################################################################
#                  Check which subprojects are there                        #
#############################################################################

AC_COIN_MAIN_SUBDIR(CoinUtils)
AC_COIN_MAIN_SUBDIR(Netlib,Data)
AC_COIN_MAIN_SUBDIR(Sample,Data)
AC_COIN_MAIN_SUBDIR(miplib3,Data)
AC_COIN_MAIN_SUBDIR(ASL,ThirdParty,solvers/asl.h)
AC_COIN_MAIN_SUBDIR(Blas,ThirdParty,daxpy.f)
AC_COIN_MAIN_SUBDIR(Lapack,ThirdParty,LAPACK/SRC/dlarf.f)
AC_COIN_MAIN_SUBDIR(Glpk,ThirdParty,glpk/include/glpk.h)
AC_COIN_MAIN_SUBDIR(Clp)
AC_COIN_MAIN_SUBDIR(DyLP)
AC_COIN_MAIN_SUBDIR(Vol)
AC_COIN_MAIN_SUBDIR(Osi)
AC_COIN_MAIN_SUBDIR(Cgl)
AC_COIN_MAIN_SUBDIR(Cbc)

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
