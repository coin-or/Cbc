# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                             libOsiCbc                                #
########################################################################

# Name of the library compiled in this directory.
lib_LTLIBRARIES = libOsiCbc.la

# List all source files for this library, including headers
libOsiCbc_la_SOURCES = OsiCbcSolverInterface.cpp OsiCbcSolverInterface.hpp

# List all additionally required libraries
if DEPENDENCY_LINKING
libOsiCbc_la_LIBADD = $(CBCLIB_LIBS) ../libCbc.la
endif

# This is for libtool (on Windows)
libOsiCbc_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
AM_CPPFLAGS = -I`$(CYGPATH_W) $(srcdir)/..` $(COINDEPEND_CFLAGS)

# Plus the appropriate include directories for the default solver.

if OSICBC_DFLT_SOLVER_CLP
  AM_CPPFLAGS += ${CLP_CFLAGS}
endif

if OSICBC_DFLT_SOLVER_CPX
  AM_CPPFLAGS += ${CPX_CFLAGS}
endif

if OSICBC_DFLT_SOLVER_DYLP
  AM_CPPFLAGS += ${DYLP_CFLAGS}
endif

if OSICBC_DFLT_SOLVER_GLPK
  AM_CPPFLAGS += ${GLPK_CFLAGS}
endif

if OSICBC_DFLT_SOLVER_GRB
  AM_CPPFLAGS += ${GRB_CFLAGS}
endif

if OSICBC_DFLT_SOLVER_MSK
  AM_CPPFLAGS += ${MSK_CFLAGS}
endif

if OSICBC_DFLT_SOLVER_HIGHS
  AM_CPPFLAGS += ${HIGHS_CFLAGS}
endif

if OSICBC_DFLT_SOLVER_SYM
  AM_CPPFLAGS += ${SYM_CFLAGS}
endif

if OSICBC_DFLT_SOLVER_XPR
  AM_CPPFLAGS += ${XPR_CFLAGS}
endif

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I..

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = OsiCbcSolverInterface.hpp
