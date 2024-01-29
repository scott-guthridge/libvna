ACLOCAL_AMFLAGS= -I m4

SUBDIRS= deb pkg rpm

deb:
	(cd deb && $(MAKE) $(AM_MAKEFLAGS) $@)

pkg:
	(cd pkg && $(MAKE) $(AM_MAKEFLAGS) $@)

srcrpm rpm rpms:
	(cd rpm && $(MAKE) $(AM_MAKEFLAGS) $@)


.PHONY: deb pkg rpm rpms srcrpm
