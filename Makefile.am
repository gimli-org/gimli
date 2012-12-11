ACLOCAL_AMFLAGS=-I m4

EXTRA_DIST = configure
SUBDIRS = src tests apps rpm

# Never build as part of 'all', 'check' or 'install', but may be built manually.
NOAUTO_DIRS = examples doc rpm

dist-hook:
	echo "dist-hook ###############################################################################"

distcheck-hook:
	echo "distcheck-hook ###############################################################################"

# Target to build source and binary rpms
rpms:
	$(MAKE) $(AM_MAKEFLAGS) -C rpm all

mrproper: maintainer-clean
	rm -f *~
	rm -rf autom4te.cache

DISTCLEANFILES = 
