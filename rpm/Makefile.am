rpmspec= $(PACKAGE_TARNAME).spec

rpmmacros=\
 --define="_topdir $${PWD}"\
 --define="_sourcedir $${PWD}/.."\
 --define="_specdir $${PWD}"

RPMBUILD= rpmbuild

srcrpm: $(rpmspec)
	$(RPMBUILD) -bs $(rpmmacros) $(rpmspec)

rpm: $(rpmspec)
	$(RPMBUILD) -ba $(rpmmacros) $(rpmspec)

rpms: rpm	# alias

clean-local:
	rm -rf BUILD BUILDROOT RPMS SRPMS

.PHONY: srcrpm rpm rpms
