package = @PACKAGE_NAME@
version = @PACKAGE_VERSION@
tarname = @PACKAGE_TARNAME@
tarfile = $(tarname)-$(version).tar.gz

MAKEPKG= makepkg
UPDPKGSUMS= updpkgsums

pkg: PKGBUILD $(tarfile)
	$(UPDPKGSUMS)
	$(MAKEPKG)

$(tarfile): $(top_builddir)/$(tarfile)
	ln -s $(top_builddir)/$(tarfile) $(tarfile)

clean-local:
	rm -rf pkg src libvna-*

.PHONY: pkg
