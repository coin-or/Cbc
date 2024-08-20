# Copyright (C) 2006, 2007 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Cbc],[2.10.12],[cbc@lists.coin-or.org])

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
AC_COIN_PROJECTDIR_INIT(Cbc,13:12:10)

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

AC_COIN_CHECK_PACKAGE(Clp, [clp osi-clp], [CbcLib CbcGeneric])
if test $coin_has_clp != yes ; then
  AC_MSG_ERROR("Required package Clp unavailable.")
fi

# The rest are not required for linking against Cbc

AC_COIN_CHECK_PACKAGE(OsiTests, [osi-unittests])
AC_COIN_CHECK_PACKAGE(Sample,   [coindatasample])
AC_COIN_CHECK_PACKAGE(Netlib,   [coindatanetlib])
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
AC_COIN_CHECK_PACKAGE(Highs,[osi-highs],  [CbcGeneric])
AC_COIN_CHECK_PACKAGE(Msk,  [osi-mosek],  [CbcGeneric])
AC_COIN_CHECK_PACKAGE(Spx,  [osi-soplex], [CbcGeneric])
AC_COIN_CHECK_PACKAGE(Xpr,  [osi-xpress])

#############################################################################
#                       Other third party software                          #
#############################################################################

# Ampl Solver library
AC_COIN_CHECK_PACKAGE(ASL, [coinasl], [CbcLib CbcGeneric])

# MUMPS
AC_COIN_CHECK_PACKAGE(Mumps, [coinmumps])

# Nauty library (for symmetry detection)
AC_LANG_PUSH(C)
AC_COIN_CHECK_USER_LIBRARY([Nauty],[NTY],[nauty.h],[nauty],,[CbcLib CbcGeneric])
AC_LANG_POP(C)

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
    osi_exists_solver=${coin_has_clp-"unavailable"}
    ;;                                             
  @<:@Dd@:>@@<:@Yy@:>@@<:@Ll@:>@@<:@Pp@:>@*)       
    osi_lc_solver=dylp ;                           
    osi_mc_solver=Dylp                             
    osi_uc_solver=DYLP                             
    osi_exists_solver=${coin_has_dylp-"unavailable"}
    ;;                                              
  @<:@Cc@:>@@<:@Pp@:>@@<:@Xx@:>@*)                  
    osi_lc_solver=cpx ;                             
    osi_mc_solver=Cpx                               
    osi_uc_solver=CPX                               
    osi_exists_solver=${coin_has_cpx-"unavailable"} 
    ;;
  @<:@Gg@:>@@<:@Ll@:>@@<:@Pp@:>@@<:@Kk@:>@*)        
    osi_lc_solver=glpk ;                            
    osi_mc_solver=Glpk                              
    osi_uc_solver=GLPK                              
    osi_exists_solver=${coin_has_glpk-"unavailable"}
    ;;
  @<:@Gg@:>@@<:@Rr@:>@@<:@Bb@:>@*)
    osi_lc_solver=grb ;
    osi_mc_solver=Grb ;           
    osi_uc_solver=GRB ;           
    osi_exists_solver=${coin_has_grb-"unavailable"}
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
    osi_exists_solver=${coin_has_msk-"unavailable"}
    ;;
  @<:@Ss@:>@@<:@Pp@:>@@<:@Xx@:>@*)
    osi_lc_solver=spx ;
    osi_mc_solver=Spx
    osi_uc_solver=SPX
    osi_exists_solver=${coin_has_spx-"unavailable"}
    ;;
  @<:@Ss@:>@@<:@Yy@:>@@<:@Mm@:>@*)
    osi_lc_solver=sym ;
    osi_mc_solver=Sym
    osi_uc_solver=SYM
    osi_exists_solver=${coin_has_sym-"unavailable"}
    ;;
  @<:@Vv@:>@@<:@Oo@:>@@<:@Ll@:>@*)
    osi_lc_solver=vol ;
    osi_mc_solver=Vol
    osi_uc_solver=VOL
    osi_exists_solver=${coin_has_vol-"unavailable"}
    ;;
  @<:@Xx@:>@@<:@Pp@:>@@<:@Rr@:>@*)
    osi_lc_solver=xpr ;
    osi_mc_solver=Xpr
    osi_uc_solver=XPR
    osi_exists_solver=${coin_has_xpr-"unavailable"}
    ;;
  *)
    osi_lc_solver=clp ;
    osi_mc_solver=Clp ;
    osi_uc_solver=CLP ;
    osi_exists_solver=${coin_has_clp-"unavailable"}
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
AM_CONDITIONAL([OSICBC_DFLT_SOLVER_GRB],[test $osi_mc_solver = Grb])
AM_CONDITIONAL([OSICBC_DFLT_SOLVER_HIGHS],[test $osi_mc_solver = HiGHS])
AM_CONDITIONAL([OSICBC_DFLT_SOLVER_MSK],[test $osi_mc_solver = Msk])
AM_CONDITIONAL([OSICBC_DFLT_SOLVER_SPX],[test $osi_mc_solver = Spx])
AM_CONDITIONAL([OSICBC_DFLT_SOLVER_SYM],[test $osi_mc_solver = Sym])
AM_CONDITIONAL([OSICBC_DFLT_SOLVER_VOL],[test $osi_mc_solver = Vol])
AM_CONDITIONAL([OSICBC_DFLT_SOLVER_XPR],[test $osi_mc_solver = Xpr])

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
                 src/OsiCbc/Makefile
                 test/Makefile
		 cbc.pc cbc-uninstalled.pc
		 osi-cbc.pc osi-cbc-uninstalled.pc])
AC_CONFIG_FILES([doxydoc/doxygen.conf])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([src/config.h src/config_cbc.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
