# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

include $(top_srcdir)/BuildTools/Makemain.inc

########################################################################
#                             libOsiCbc                                #
########################################################################

# Name of the library compiled in this directory.

lib_LTLIBRARIES = libOsiCbc.la

# List all source files for this library, including headers
libOsiCbc_la_SOURCES = OsiCbcSolverInterface.cpp OsiCbcSolverInterface.hpp

# List all additionally required libraries
libOsiCbc_la_LIBADD = ../libCbc.la

# This is for libtool
libOsiCbc_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = -DOSICBC_BUILD -I$(srcdir)/.. $(CBCLIB_CFLAGS)

# Plus the appropriate include directories and libraries for the default solver.
if OSICBC_DFLT_SOLVER_CLP
  AM_CPPFLAGS += $(OSICLP_CFLAGS)
  libOsiCbc_la_LIBADD += $(OSICLP_LFLAGS)
endif

if OSICBC_DFLT_SOLVER_CPX
  AM_CPPFLAGS += $(OSICPX_CFLAGS)
  libOsiCbc_la_LIBADD += $(OSICPX_LFLAGS)
endif

if OSICBC_DFLT_SOLVER_DYLP
  AM_CPPFLAGS += $(OSIDYLP_CFLAGS)
  libOsiCbc_la_LIBADD += $(OSIDYLP_LFLAGS)
endif

if OSICBC_DFLT_SOLVER_GLPK
  AM_CPPFLAGS += $(OSIGLPK_CFLAGS)
  libOsiCbc_la_LIBADD += $(OSIGLPK_LFLAGS)
endif

#if OSICBC_DFLT_SOLVER_GRB
#  AM_CPPFLAGS += ${GRB_CFLAGS}
#endif

if OSICBC_DFLT_SOLVER_MSK
  AM_CPPFLAGS += $(OSIMSK_CFLAGS)
  libOsiCbc_la_LIBADD += $(OSIMSK_LFLAGS)
endif

if OSICBC_DFLT_SOLVER_HIGHS
  AM_CPPFLAGS += $(HIGHS_CFLAGS)
  libOsiCbc_la_LIBADD += $(HIGHS_LFLAGS)
endif

if OSICBC_DFLT_SOLVER_SYM
  AM_CPPFLAGS += $(OSISYM_CFLAGS)
  libOsiCbc_la_LIBADD += $(OSISYM_LFLAGS)
endif

#if OSICBC_DFLT_SOLVER_XPR
#  AM_CPPFLAGS += ${XPR_CFLAGS}
#endif

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin-or'

includecoindir = $(pkgincludedir)
includecoin_HEADERS = OsiCbcSolverInterface.hpp

install-exec-local:
	$(install_sh_DATA) config_osicbc.h $(DESTDIR)$(includecoindir)/OsiCbcConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/OsiCbcConfig.h
