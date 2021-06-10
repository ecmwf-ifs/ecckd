dist_pkgdata_DATA = README ChangeLog NOTICE LICENSE
SUBDIRS = src/tools src/lbfgs src/ecckd
EXTRA_DIST = configure_ecmwf.sh src/include/*.h test/*.sh test/config.h test/check_configuration.h plot/*.m doc/Makefile doc/ecckd_documentation.tex doc/flowchart.eps doc/ecckd_documentation.pdf
ACLOCAL_AMFLAGS = -I m4

clean-local:
	rm -f test/config_*.cfg test/config_*.nam
