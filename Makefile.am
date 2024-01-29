ACLOCAL_AMFLAGS= -I m4
libvna_ladir= $(libdir)

SUBDIRS= src debian package

examplesdir = $(docdir)/examples

dist_doc_DATA = COPYING
dist_examples_DATA = README-EXAMPLES
dist_noinst_SCRIPTS = bootstrap

dist_noinst_DATA = \
	README-BUILD.md \
	README-DEVELOPMENT.md

deb: dist
	(cd package && $(MAKE) $(AM_MAKEFLAGS) $@)

srcrpm rpm rpms: dist
	(cd package && $(MAKE) $(AM_MAKEFLAGS) $@)

pkg: dist
	(cd package && $(MAKE) $(AM_MAKEFLAGS) $@)

plot: all
	(cd src && $(MAKE) $(AM_MAKEFLAGS) $@)

pdfman:
	(cd src && $(MAKE) $(AM_MAKEFLAGS) $@)

#clean-local:
#	rm -f

.PHONY: deb pdfman pkg plot rpm rpms srcrpm
