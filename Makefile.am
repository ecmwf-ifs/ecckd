dist_pkgdata_DATA = README ChangeLog NOTICE
SUBDIRS = src/tools src/lbfgs src/ecckd
EXTRA_DIST = configure_ecmwf.sh src/include test/*.sh test/config.h test/check_configuration.h plot/*.m
ACLOCAL_AMFLAGS = -I m4

clean-local:
	rm -f test/config_*.cfg test/config_*.nam
