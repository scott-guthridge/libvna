ACLOCAL_AMFLAGS= -I m4

SUBDIRS= deb rpm

deb:
	(cd deb && $(MAKE) $(AM_MAKEFLAGS) $@)

srcrpm rpm rpms:
	(cd rpm && $(MAKE) $(AM_MAKEFLAGS) $@)


.PHONY: deb rpm rpms srcrpm
