# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                          Subdirectories                              #
########################################################################

SUBDIRS = src

# We don't want to compile the test subdirectory, unless the test target is
# specified.  But we need to list it as subdirectory to make sure that it is
# included in the tarball

if ALWAYS_FALSE
  SUBDIRS += test
endif

########################################################################
#             Additional files to be included in tarball               #
########################################################################

# Here we need include all files that are not mentioned in other Makefiles

EXTRA_DIST = \
	doc/authors.xml \
	doc/cbcexe.xml \
	doc/cbcmodelclass.xml \
	doc/cbcuserguide.xml \
	doc/doxygen.xml \
	doc/faqcontent.xml \
	doc/faq.xml \
	doc/intro.xml \
	doc/legal.xml \
	doc/messages.xml \
	doc/moresamples.xml \
	doc/osibuild.xml \
	doc/otherclasses.xml \
	doc/revhist.xml \
	examples/CbcBranchFollow2.cpp \
	examples/CbcBranchFollow2.hpp \
	examples/CbcBranchLink.cpp \
	examples/CbcBranchLink.hpp \
	examples/CbcBranchUser.cpp \
	examples/CbcBranchUser.hpp \
	examples/CbcCompareUser.cpp \
	examples/CbcCompareUser.hpp \
	examples/CbcSolver2.cpp \
	examples/CbcSolver2.hpp \
	examples/CbcSolver3.cpp \
	examples/CbcSolver3.hpp \
	examples/CbcSolverLongThin.cpp \
	examples/CbcSolverLongThin.hpp \
	examples/CbcSolverLink.hpp \
	examples/CbcSolverLink.cpp \
	examples/ClpQuadInterface.cpp \
	examples/ClpQuadInterface.hpp \
	examples/OsiBranchLink.cpp \
	examples/OsiBranchLink.hpp \
	examples/OsiSolverLink.hpp \
	examples/OsiSolverLink.cpp \
	examples/crew.cpp \
	examples/driver2.cpp \
	examples/driver.cpp \
	examples/fast0507b.cpp \
	examples/fast0507.cpp \
	examples/gear.cpp \
	examples/hotstart.cpp \
	examples/link.cpp \
	examples/longthin.cpp \
	examples/lotsize.cpp \
	examples/Makefile.in \
	examples/minimum.cpp \
	examples/nway.cpp \
	examples/qmip.cpp \
	examples/quad2.mps \
	examples/quad.mps \
	examples/repeat.cpp \
	examples/sample1.cpp \
	examples/sample2.cpp \
	examples/sample3.cpp \
	examples/sample4.cpp \
	examples/sample5.cpp \
	examples/sos.cpp \
	examples/sudoku.cpp \
	examples/sudoku_sample.csv

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd test; $(MAKE) test

unitTest: test

clean-local:
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

distclean-local:
	if test -r test/Makefile; then cd test; $(MAKE) distclean; fi

install-exec-local: install-doc

uninstall-local: uninstall-doc

.PHONY: test unitTest

########################################################################
#                  Installation of the addlibs file                    #
########################################################################

if COIN_HAS_PKGCONFIG
addlibsdir = $(datadir)/coin/doc/Cbc
addlibs_DATA = cbc_addlibs.txt

cbc_addlibs.txt : cbc.pc
	PKG_CONFIG_PATH=.:@COIN_PKG_CONFIG_PATH@ $(PKG_CONFIG) --libs cbc > cbc_addlibs.txt
endif

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = cbc.pc osi-cbc.pc

########################################################################
#                         Maintainer Stuff                             #
########################################################################

CLEANFILES = cbc_addlibs.txt

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =

include ../BuildTools/Makemain.inc
