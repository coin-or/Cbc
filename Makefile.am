# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter           IBM    2006-04-13

include BuildTools/Makemain.inc

########################################################################
#                          Subdirectories                              #
########################################################################

SUBDIRS = src/CoinUtils src/Clp src/Cgl src

pkgconfiglib_DATA = mipster.pc

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd test; $(MAKE) test

unitTest: test

clean-local: clean-doxygen-docs
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

install-exec-local: install-doc install-man1 install-param-ref

install-man1:
	@$(MKDIR_P) "$(DESTDIR)$(mandir)/man1"
	$(install_sh_DATA) $(srcdir)/doc/mipster.1 "$(DESTDIR)$(mandir)/man1/mipster.1"

install-param-ref:
	@$(MKDIR_P) "$(DESTDIR)$(docdir)"
	$(install_sh_DATA) $(srcdir)/doc/mipster-parameters.tex "$(DESTDIR)$(docdir)/mipster-parameters.tex"

uninstall-local: uninstall-doc uninstall-doxygen-docs uninstall-man1 uninstall-param-ref

uninstall-man1:
	rm -f "$(DESTDIR)$(mandir)/man1/mipster.1"

uninstall-param-ref:
	rm -f "$(DESTDIR)$(docdir)/mipster-parameters.tex"

.PHONY: test unitTest doxydoc
