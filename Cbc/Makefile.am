# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                          Subdirectories                              #
########################################################################

SUBDIRS = src

# Don't do src/OsiCbc recursively, since src/OsiCbc/libOsiCbc depends on src/libCbc
SUBDIRS += src/OsiCbc

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

# Doxygen documentation

doxydoc:
	doxygen doxydoc/doxygen.conf

clean-doxydoc:
	( cd doxydoc ; rm -rf html *.log *.tag )


clean-local: clean-doxydoc
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

distclean-local:
	if test -r test/Makefile; then cd test; $(MAKE) distclean; fi

install-exec-local: install-doc

uninstall-local: uninstall-doc

.PHONY: test unitTest doxydoc

########################################################################
#                  Installation of the addlibs file                    #
########################################################################

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = cbc.pc osi-cbc.pc

addlibsdir = $(DESTDIR)$(datadir)/coin/doc/Cbc

install-data-hook:
	@$(mkdir_p) "$(addlibsdir)"
if COIN_HAS_PKGCONFIG
	PKG_CONFIG_PATH=@COIN_PKG_CONFIG_PATH@:$(DESTDIR)$(pkgconfiglibdir) \
	$(PKG_CONFIG) --libs cbc > $(addlibsdir)/cbc_addlibs.txt
else
if COIN_CXX_IS_CL
	echo "-libpath:`$(CYGPATH_W) @abs_lib_dir@` libCbc.lib @CBCLIB_LIBS_INSTALLED@" > $(addlibsdir)/cbc_addlibs.txt
else
	echo -L@abs_lib_dir@ -lCbcSolver -lCbc @CBCLIB_LIBS_INSTALLED@ > $(addlibsdir)/cbc_addlibs.txt
endif
endif

uninstall-hook:
	rm -f $(addlibsdir)/cbc_addlibs.txt

########################################################################
#                         Maintainer Stuff                             #
########################################################################

CLEANFILES =

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =

include BuildTools/Makemain.inc
