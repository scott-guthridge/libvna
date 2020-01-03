rpmspec= $(PACKAGE_TARNAME).spec

rpmmacros=\
 --define="_rpmdir $${PWD}"\
 --define="_srcrpmdir $${PWD}"\
 --define="_sourcedir $${PWD}/.."\
 --define="_specdir $${PWD}"\
 --define="_builddir $${PWD}"

RPMBUILD= rpmbuild
RPMFLAGS= --buildroot="$${PWD}/root"

srcrpm: $(rpmspec)
	$(RPMBUILD) $(RPMFLAGS) -bs $(rpmmacros) $(rpmspec)

rpms: $(rpmspec)
	$(RPMBUILD) $(RPMFLAGS) -ba $(rpmmacros) $(rpmspec)

clean-local:
	rm -rf $(PACKAGE_TARNAME)-$(PACKAGE_VERSION) root

.PHONY: srcrpm rpms
