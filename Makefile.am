ACLOCAL_AMFLAGS= -I m4
libvna_ladir= $(libdir)

SUBDIRS= src deb rpm

examplesdir = $(docdir)/examples

dist_doc_DATA = COPYING
dist_examples_DATA = README-EXAMPLES
dist_noinst_SCRIPTS = bootstrap

dist_noinst_DATA = \
	debian/changelog \
	debian/control \
	debian/copyright \
	debian/libvna-dev.dirs \
	debian/libvna-dev.docs \
	debian/libvna-dev.examples \
	debian/libvna-dev.install \
	debian/libvna-dev.manpages \
	debian/libvna.dirs \
	debian/libvna.install \
	debian/not-installed \
	debian/postinst_run \
	debian/postrm_run \
	debian/rules \
	debian/source/format \
	debian/source/options

deb: dist
	(cd deb && $(MAKE) $(AM_MAKEFLAGS) $@)

srcrpm rpm rpms: dist
	(cd rpm && $(MAKE) $(AM_MAKEFLAGS) $@)

plot: all
	(cd src && $(MAKE) $(AM_MAKEFLAGS) $@)

pdfman:
	(cd src && $(MAKE) $(AM_MAKEFLAGS) $@)

#clean-local:
#	rm -f

.PHONY: srcrpm rpm plot pdfman
