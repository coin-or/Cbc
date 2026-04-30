# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter           IBM    2006-04-13

include BuildTools/Makemain.inc

########################################################################
#                          Subdirectories                              #
########################################################################

SUBDIRS = src src/OsiCbc

pkgconfiglib_DATA = cbc.pc osi-cbc.pc

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd test; $(MAKE) test

unitTest: test

clean-local: clean-doxygen-docs
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

install-exec-local: install-doc install-man1

install-man1:
	@$(MKDIR_P) "$(DESTDIR)$(mandir)/man1"
	$(install_sh_DATA) $(srcdir)/doc/cbc.1 "$(DESTDIR)$(mandir)/man1/cbc.1"

uninstall-local: uninstall-doc uninstall-doxygen-docs uninstall-man1

uninstall-man1:
	rm -f "$(DESTDIR)$(mandir)/man1/cbc.1"

.PHONY: test unitTest doxydoc
