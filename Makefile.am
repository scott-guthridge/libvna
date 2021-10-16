ACLOCAL_AMFLAGS= -I m4
libvna_ladir= $(libdir)

SUBDIRS= rpm src

dist_doc_DATA= COPYING README-EXAMPLES
dist_noinst_SCRIPTS= bootstrap

srcrpm rpm rpms: dist
	(cd rpm && $(MAKE) $(AM_MAKEFLAGS) $@)

plot: all
	(cd src && $(MAKE) $(AM_MAKEFLAGS) $@)

pdfman:
	(cd src && $(MAKE) $(AM_MAKEFLAGS) $@)

#clean-local:
#	rm -f

.PHONY: srcrpm rpm plot pdfman
