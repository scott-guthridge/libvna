ACLOCAL_AMFLAGS= -I m4
libvna_ladir= $(libdir)

SUBDIRS= rpm src

dist_doc_DATA= COPYING README-EXAMPLES
dist_noinst_SCRIPTS= bootstrap

srcrpm rpms: dist
	(cd rpm && $(MAKE) $(AM_MAKEFLAGS) $@)

plot: all
	(cd src && $(MAKE) $(AM_MAKEFLAGS) $@)

clean-local:
	rm -f archdep.o

.PHONY: srcrpm rpms plot
