# Copyright (C) 2006, 2007 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Cbc],[trunk],[cbc@lists.coin-or.org],[coin-or-cbc],
	[https://projects.coin-or.org/Cbc])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Eclipse Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there.

AC_CONFIG_SRCDIR(src/CbcTree.hpp)

# Do some initialization work (version numbers, change prefix default, ...)

AC_COIN_INITIALIZE(Cbc)

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
#                            COIN-OR components                             #
#############################################################################

m4_foreach_w([myvar],[CoinUtils Osi Cgl Clp],
  [ AC_COIN_CHK_PKG(myvar,[CbcLib CbcGeneric])
    if test $coin_has_[]m4_tolower(myvar) != yes ; then
      AC_MSG_ERROR([Required package myvar is not available.])
    fi
  ])

AC_COIN_CHK_PKG(OsiClp,[CbcLib CbcGeneric],[osi-clp])
if test $coin_has_osiclp != yes ; then
  AC_MSG_ERROR("Required package OsiClp unavailable.")
fi

# Optional package for command line edit support

AC_COIN_CHK_GNU_READLINE([CbcLib CbcGeneric])

#############################################################################
#                       Other third party software                          #
#############################################################################

# Ampl Solver library

AC_COIN_CHK_PKG(ASL,[CbcSolverLib],[coinasl])

# Nauty library (for symmetry detection).

AC_COIN_CHK_LIB([Nauty],[CbcLib CbcGeneric],[-lnauty],[],[],[nauty])

#############################################################################

# Finalize flags for CbcLib, add the result to OsiCbc, and finalize flags for
# it.

# The rest are not required for linking against Cbc but are used for the unit
# test. For the data packages, we just want XXX_DATA for use in Makefile.am.

AC_COIN_CHK_PKG(OsiTests,[CbcUnitTest],[osi-unittests])
AC_COIN_CHK_PKG(Sample,,[coindatasample],,dataonly)
AC_COIN_CHK_PKG(Netlib,,[coindatanetlib],,dataonly)
AC_COIN_CHK_PKG(Miplib3,,[coindatamiplib3],,dataonly)

AC_COIN_FINALIZE_FLAGS(CbcUnitTest)

#############################################################################
#                    LP solvers other than CLP                              #
#############################################################################

# Check which other LP solvers are available, some of them become a dependency
# of CbcGeneric

AC_COIN_CHK_PKG(OsiDyLP,[CbcGeneric],[osi-dylp])
AC_COIN_CHK_PKG(OsiCpx,[CbcLib CbcGeneric CbcSolverLib],[osi-cplex])
AC_COIN_CHK_PKG(OsiGlpk,[CbcGeneric],[osi-glpk])
#AC_COIN_CHK_PKG(OsiGrb,,[osi-gurobi])
AC_COIN_CHK_PKG(OsiMsk,[CbcGeneric],[osi-mosek])
AC_COIN_CHK_PKG(OsiSpx,[CbcGeneric],[osi-soplex])
#AC_COIN_CHK_PKG(OsiXpr,,[osi-xpress])

AC_COIN_CHK_LIB(Cplex,[CbcLib],[-lcplex -lpthread -lm -ldl],[],[],[CPXgetstat])
AC_COIN_CHK_PKG(Glpk,[CbcSolverLib],,[build])

#TODO organize this mess

AC_COIN_FINALIZE_FLAGS(CbcGeneric)

AC_COIN_FINALIZE_FLAGS([CbcLib])
AC_COIN_CHK_HERE([CbcLib],[OsiCbcLib])
AC_COIN_FINALIZE_FLAGS(OsiCbcLib)

AC_COIN_FINALIZE_FLAGS([CbcSolverLib])
AC_COIN_CHK_HERE([CbcLib],[CbcSolverLib],[cbc])

#############################################################################
#                       CbcGeneric configuration                            #
#############################################################################

# Are we building cbc-generic? Default is no (false).

AC_ARG_WITH([cbc-generic],
  AS_HELP_STRING([--with-cbc-generic],[specify whether to build cbc-generic (default: no)]),
  [case $withval in
    yes) cbc_with_cbc_generic=true ;;
      *) cbc_with_cbc_generic=false ;;
   esac],
  [cbc_with_cbc_generic=false])
if test x"$cbc_with_cbc_generic" = xyes ; then
  AC_MSG_NOTICE([building cbc-generic])
fi
AM_CONDITIONAL(CBC_BUILD_CBC_GENERIC, test x"$cbc_with_cbc_generic" = xtrue)

# Set the default solver for cbc-generic. In theory, any OsiXXX should work.
# In practice, only the three listed below have had any testing.

AC_MSG_CHECKING(for cbc-generic default solver)
AC_ARG_WITH([cbc-generic-solver],
            AC_HELP_STRING([--with-cbc-generic-solver],
                           [specify default solver for cbc-generic in lower case letters (if not given, clp is assumed)]),
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
  highs)
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
#                       Cbc parallel configuration	                         #
#############################################################################

# First we define a new option, --enable-cbc-parallel

AC_ARG_ENABLE([cbc-parallel],
[AC_HELP_STRING([--enable-cbc-parallel],
                [enables compilation of the SMP version of Cbc])])

if test "$enable_cbc_parallel" = yes; then
  # Define the preprocessor macro
  AC_DEFINE([CBC_THREAD],[1],[Define to 1 if the SMP version of Cbc should be compiled])
  case "$CXX" in
  *cl) # Windows / cl
    # TODO we should check whether the library works and pthread.h is indeed there
    AC_ARG_WITH(pthreadsw32-lib, 
      AC_HELP_STRING([--with-pthreadsw32-lib],
        [specify the name of the pthreads win32 library (for cl/icl builds only)]),
      [CBCLIB_LFLAGS="$withval $CBCLIB_LFLAGS"],
      [AC_MSG_WARN([--enable-cbc-parallel selected, but --with-pthreadsw32-lib not given])
       CBCLIB_LFLAGS="pthreadVCE2.lib $CBCLIB_LFLAGS"])
    
    AC_ARG_WITH(pthreadsw32-incdir,
      AC_HELP_STRING([--with-pthreadsw32-incdir],
        [specify the path of the pthreads win32 header file (for cl/icl builds only)]),
      [CBCLIB_CFLAGS="$CBCLIB_CFLAGS -I$withval"],
      [AC_MSG_WARN([--enable-cbc-parallel selected, but --with-pthreadsw32-incdir not given])])
    ;;
  *) # not Windows
    AC_CHECK_LIB([pthread],[pthread_create],
      [CBCLIB_LFLAGS="-lpthread $CBCLIB_LFLAGS"],
      [AC_MSG_ERROR([--enable-cbc-parallel selected, but -lpthreads unavailable])])
    AC_CHECK_LIB([rt],[clock_gettime],
      [CBCLIB_LFLAGS="-lrt $CBCLIB_LFLAGS"
       AC_DEFINE([HAVE_CLOCK_GETTIME],[1],[Define if clock_gettime and rt library is available])
                 ],
      [AC_MSG_WARN([--enable-cbc-parallel selected, but -lrt unavailable; using gettimeofday instead of clock_gettime])])
  ;;
  esac
  AC_MSG_NOTICE([Cbc multithreading enabled]);
fi

##############################################################################
#                            OsiCbc configuration                            #
##############################################################################

# Handles configuration of the underlying default solver in OsiCbc. The issue
# is that OsiCbc defines a default solver, used when the client does not     
# specify a solver in the constructor. The default solver must, therefore,   
# be present in the build. This macro checks that this is true, and sets the 
# compile-time symbols OSICBC_DFLT_SOLVER, OSICBC_CLP_DFLT_SOLVER, and
# OSICBC_DFLT_SOLVER_HPP that control the build. The parameter default_solver
# should normally be clp, unless you're working on some other solver.

# Ideally, this macro would not require enumeration of solvers, but the
# effort required to avoid it is just not justified at present. One enumeration
# is hidden in AC_OSI_CANONICAL. The other is visible, the AM_CONDITIONAL
# list at the end.                                                             

# In an ideal world this macro would do absolutely nothing if Cbc is not
# present, but autotools is not an ideal world and we have to satisfy its
# requirements. In particular, the AM_CONDITIONAL macros need to execute or
# automake will complain. Really the only thing we need to suppress is the 
# check that the default solver exists. All the rest is irrelevant when Cbc
# isn't present (hence OsiCbc will be configured but not actually compiled).

# Process the with-osicbc-default-solver option.

AC_ARG_WITH([osicbc-default-solver],
  AS_HELP_STRING([--with-osicbc-default-solver],
     [specify underlying solver for OsiCbc (default clp)]),
     [osicbc_with_solver=$withval],
     [osicbc_with_solver=clp])

case "${osicbc_with_solver}" in
  @<:@Cc@:>@@<:@Ll@:>@@<:@Pp@:>@*)
    osi_lc_solver=clp ;           
    osi_mc_solver=Clp ;           
    osi_uc_solver=CLP ;           
    osi_exists_solver=${coin_has_osiclp-"unavailable"}
    ;;                                             
  @<:@Dd@:>@@<:@Yy@:>@@<:@Ll@:>@@<:@Pp@:>@*)       
    osi_lc_solver=dylp ;                           
    osi_mc_solver=Dylp                             
    osi_uc_solver=DYLP                             
    osi_exists_solver=${coin_has_osidylp-"unavailable"}
    ;;                                              
  @<:@Cc@:>@@<:@Pp@:>@@<:@Xx@:>@*)                  
    osi_lc_solver=cpx ;                             
    osi_mc_solver=Cpx                               
    osi_uc_solver=CPX                               
    osi_exists_solver=${coin_has_osicpx-"unavailable"} 
    ;;
  @<:@Gg@:>@@<:@Ll@:>@@<:@Pp@:>@@<:@Kk@:>@*)        
    osi_lc_solver=glpk ;                            
    osi_mc_solver=Glpk                              
    osi_uc_solver=GLPK                              
    osi_exists_solver=${coin_has_osiglpk-"unavailable"}
    ;;
  @<:@Gg@:>@@<:@Rr@:>@@<:@Bb@:>@*)
    osi_lc_solver=grb ;
    osi_mc_solver=Grb ;           
    osi_uc_solver=GRB ;           
    osi_exists_solver=${coin_has_osigrb-"unavailable"}
    ;;                                             
  @<:@H@:>@@<:@i@:>@@<:@g@:>@@<:@h@:>@@<:@s@:>@*)        
    osi_lc_solver=highs ;                            
    osi_mc_solver=HiGHS ;                            
    osi_uc_solver=HIGHS ;                              
    osi_exists_solver=${coin_has_highs-"unavailable"}
    ;;
  @<:@Mm@:>@@<:@Ss@:>@@<:@Kk@:>@*)
    osi_lc_solver=msk ;
    osi_mc_solver=Msk
    osi_uc_solver=MSK
    osi_exists_solver=${coin_has_osimsk-"unavailable"}
    ;;
  @<:@Ss@:>@@<:@Pp@:>@@<:@Xx@:>@*)
    osi_lc_solver=spx ;
    osi_mc_solver=Spx
    osi_uc_solver=SPX
    osi_exists_solver=${coin_has_osispx-"unavailable"}
    ;;
  @<:@Ss@:>@@<:@Yy@:>@@<:@Mm@:>@*)
    osi_lc_solver=sym ;
    osi_mc_solver=Sym
    osi_uc_solver=SYM
    osi_exists_solver=${coin_has_osisym-"unavailable"}
    ;;
  @<:@Xx@:>@@<:@Pp@:>@@<:@Rr@:>@*)
    osi_lc_solver=xpr ;
    osi_mc_solver=Xpr
    osi_uc_solver=XPR
    osi_exists_solver=${coin_has_osixpr-"unavailable"}
    ;;
  *)
    osi_lc_solver=clp ;
    osi_mc_solver=Clp ;
    osi_uc_solver=CLP ;
    osi_exists_solver=${coin_has_osiclp-"unavailable"}
    AC_MSG_WARN([Unrecognised solver $1; defaulting to $osi_lc_solver.])
    ;;
esac

# Check that the requested solver is available.
if test $osi_exists_solver = no; then
  AC_MSG_ERROR([selected default solver $osicbc_with_solver is unavailable.
  Please select an available solver using the --with-osicbc-default-solver option.])
fi

# State the result.
AC_MSG_NOTICE([OsiCbc default solver is $osi_lc_solver])

# And set the configuration variables.
AC_DEFINE_UNQUOTED([OSICBC_DFLT_SOLVER],
  [Osi${osi_mc_solver}SolverInterface],
  [define to the name of the default solver interface class, e.g., OsiClpSolverInterface])
AC_DEFINE_UNQUOTED([OSICBC_DFLT_SOLVER_HPP],
  ["Osi${osi_mc_solver}SolverInterface.hpp"],
  [define to the name of the .hpp file for the default solver interface class, e.g., "OsiClpSolverInterface.hpp" (include quotes)])
if test $osi_mc_solver = "Clp"; then
  AC_DEFINE([OSICBC_DFLT_SOLVER_CLP],[1],[define this symbol if clp is the default solver])
fi

# Last but not least, we need automake conditionals.
AM_CONDITIONAL([OSICBC_DFLT_SOLVER_CLP],[test $osi_mc_solver = Clp])
AM_CONDITIONAL([OSICBC_DFLT_SOLVER_CPX],[test $osi_mc_solver = Cpx])
AM_CONDITIONAL([OSICBC_DFLT_SOLVER_DYLP],[test $osi_mc_solver = Dylp])
AM_CONDITIONAL([OSICBC_DFLT_SOLVER_GLPK],[test $osi_mc_solver = Glpk])
#AM_CONDITIONAL([OSICBC_DFLT_SOLVER_GRB],[test $osi_mc_solver = Grb])
AM_CONDITIONAL([OSICBC_DFLT_SOLVER_HIGHS],[test $osi_mc_solver = HiGHS])
AM_CONDITIONAL([OSICBC_DFLT_SOLVER_MSK],[test $osi_mc_solver = Msk])
AM_CONDITIONAL([OSICBC_DFLT_SOLVER_SPX],[test $osi_mc_solver = Spx])
AM_CONDITIONAL([OSICBC_DFLT_SOLVER_SYM],[test $osi_mc_solver = Sym])
#AM_CONDITIONAL([OSICBC_DFLT_SOLVER_XPR],[test $osi_mc_solver = Xpr])

##############################################################################
#                   VPATH links for example input files                      #
##############################################################################

# In case this is a VPATH configuration we need to make sure that the
# input files for the examples are available in the VPATH directory.

# AC_COIN_VPATH_LINK(examples/sudoku_sample.csv)
# AC_COIN_VPATH_LINK(examples/quad.mps)
# AC_COIN_VPATH_LINK(examples/quad2.mps)

#############################################################################
#                  Check for doxygen                                        #
#############################################################################

AC_COIN_DOXYGEN([CoinUtils Osi Clp DyLP Cgl])

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 examples/Makefile
                 src/Makefile
                 src/OsiCbc/Makefile
                 test/Makefile
                 cbc.pc
                 osi-cbc.pc])
AC_CONFIG_FILES([doxydoc/doxygen.conf])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([src/config.h src/config_cbc.h src/config_cbcsolver.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
