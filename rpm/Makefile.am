# To get a spec file w/ the version in it, we copy 'specfile' to the final name
# This is also our chance to substitute in an RPM release number and a date.
RPM_SPECFILE = $(PACKAGE)-$(VERSION).spec
RPM_RELEASE = 1
RPM_DATE = $(shell date +"%a %b %d %Y")

$(RPM_SPECFILE): specfile
	@SED@ -e 's/\@rpm_release\@/$(RPM_RELEASE)/' \
	           -e 's/\@date\@/$(RPM_DATE)/' $< >$@

# Unless told otherwise, try to build for current CPU
RPMBUILD_FLAGS = --target @CR_CPU@
EXTRA_DIST = $(RPM_SPECFILE) license.txt
CLEANFILES = $(RPM_SPECFILE)

if HAVE_RPMBUILD
RPM_TARFILE = $(top_builddir)/$(PACKAGE)-$(VERSION).tar.gz
RPM_CMD = @RPMBUILD@ --define '_topdir '`pwd` $(RPMBUILD_FLAGS)
.PHONY: _rpm _all rpm rpm_prep
_srpm:
	@mkdir -p SRPMS SPECS 
	$(RPM_CMD) -bs $(RPM_SPECFILE)
	@rm -rf SPECS
_rpm:
	@mkdir -p RPMS SPECS BUILD
	$(RPM_CMD) -bb $(RPM_SPECFILE)
	@rm -rf SPECS BUILD
_all:
	@mkdir -p SRPMS RPMS SPECS BUILD
	$(RPM_CMD) -ba $(RPM_SPECFILE)
	@rm -rf SPECS BUILD
rpm_prep:
	@rm -f $(RPM_TARFILE) $(RPM_SPECFILE) # ensure we start fresh, specfile will be rebuilt by dist-gzip
	$(RPM_SUDO) $(MAKE) $(AM_MAKEFLAGS) -C $(top_builddir) dist-gzip TAR_OPTIONS="--owner=0 --group=0"
	$(MAKE) $(AM_MAKEFLAGS) $(RPM_SPECFILE)
	@mkdir -p SOURCES
	cp $(RPM_TARFILE) SOURCES/$(PACKAGE)-$(VERSION).tar.gz

srpm: rpm_prep
	$(MAKE) $(AM_MAKEFLAGS) _srpm
rpm: rpm_prep
	$(MAKE) $(AM_MAKEFLAGS) _rpm
all: rpm_prep
	@true

clean-local:
	@rm -rf RPMS SRPMS SPECS BUILD SOURCES $(RPM_SPECFILE)

distcheck-hook: _all
	$(MAKE) $(AM_MAKEFLAGS) clean
else
all:
	@echo "==========================================="
	@echo "Don't know how to build RPMs on your system"
	@echo "==========================================="
distcheck-hook:
	@true
endif
