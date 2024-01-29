ACLOCAL_AMFLAGS= -I m4

package = @PACKAGE_NAME@
version = @PACKAGE_VERSION@
tarname = @PACKAGE_TARNAME@
distdir = $(tarname)-$(version)

ORIGTAR=$(tarname)_$(version).orig.tar.gz

deb: $(distdir)
	( cd $(distdir) && debuild -uc -us --lintian-opts \
	    --suppress-tags file-without-copyright-information )

$(distdir): $(ORIGTAR)
	tar xzf $(ORIGTAR)
	( cd "$(distdir)" && \
	  sed -i -e 's/%%VERSION%%/$(version)/g' debian/changelog )
	dpkg-source --commit $(distdir)

$(ORIGTAR): $(top_builddir)/$(distdir).tar.gz
	ln -s $(top_builddir)/$(distdir).tar.gz "$@"

clean-local:
	rm -rf $(distdir)
	rm -f  $(ORIGTAR) $(tarname)_$(version)-1.debian.tar.xz \
		*.build *.buildinfo *.changes *.deb *.dsc

.PHONY: deb
