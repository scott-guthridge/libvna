#
# Vector Network Analyzer Library
# Copyright © 2020, 2021 D Scott Guthridge <scott_guthridge@rompromity.net>
#
ACLOCAL_AMFLAGS= -I m4
GNUPLOT= gnuplot
PDFROFF= pdfroff
SUFFIXES= .3 .pdf

SUBDIRS= . tests

#
# Library
#
include_HEADERS = vnacal.h vnaconv.h vnadata.h vnaerr.h vnafile.h
lib_LTLIBRARIES = libvna.la
libvna_la_SOURCES = archdep.h archdep.c vnacal_internal.h \
	vnacal_delete_parameter.c vnacal_apply.c vnacal_new_add_common.c \
	vnacal_new_parameter.c vnacal_new_set_measurement_error.c \
	vnacal_new_solve.c vnacal_new.c vnacal_create.c vnacal_calibration.c \
	vnacal_make_correlated_parameter.c vnacal_make_scalar_parameter.c \
	vnacal_make_unknown_parameter.c vnacal_make_vector_parameter.c \
	vnacal_add_calibration.c vnacal_delete_calibration.c \
	vnacal_layout.h vnacal_layout.c vnacal_get_parameter_value.c \
	vnacal_error.c vnacal_find_calibration.c vnacal_free.c vnacal_get.c \
	vnacal_load.c vnacal_property.c vnacal_rfi.c vnacal_save.c \
	vnacal_set_dprecision.c vnacal_set_fprecision.c vnacal_parameter.c \
	vnacommon_internal.h vnacommon_lu.c vnacommon_mmultiply.c \
	vnacommon_minverse.c vnacommon_mldivide.c vnacommon_mrdivide.c \
	vnacommon_qrd.c vnacommon_qr.c vnacommon_qrsolve.c \
	vnacommon_qrsolve2.c vnacommon_spline.c \
	vnaerr_internal.h vnaerr_verror.c \
	vnaconv_internal.h vnaconv_atob.c vnaconv_atog.c vnaconv_atoh.c \
	vnaconv_atos.c vnaconv_atot.c vnaconv_atoy.c vnaconv_atoz.c \
	vnaconv_atozi.c vnaconv_btoa.c vnaconv_btog.c vnaconv_btoh.c \
	vnaconv_btos.c vnaconv_btot.c vnaconv_btoy.c vnaconv_btoz.c \
	vnaconv_btozi.c vnaconv_gtoa.c vnaconv_gtob.c vnaconv_gtoh.c \
	vnaconv_gtos.c vnaconv_gtot.c vnaconv_gtoy.c vnaconv_gtoz.c \
	vnaconv_gtozi.c vnaconv_htoa.c vnaconv_htob.c vnaconv_htog.c \
	vnaconv_htos.c vnaconv_htot.c vnaconv_htoy.c vnaconv_htoz.c \
	vnaconv_htozi.c vnaconv_stoa.c vnaconv_stob.c vnaconv_stog.c \
	vnaconv_stoh.c vnaconv_stot.c vnaconv_stoy.c vnaconv_stoyn.c \
	vnaconv_stoz.c vnaconv_stozi.c vnaconv_stozimn.c vnaconv_stozin.c \
	vnaconv_stozn.c vnaconv_ttoa.c vnaconv_ttob.c vnaconv_ttog.c \
	vnaconv_ttoh.c vnaconv_ttos.c vnaconv_ttoy.c vnaconv_ttoz.c \
	vnaconv_ttozi.c vnaconv_ytoa.c vnaconv_ytob.c vnaconv_ytog.c \
	vnaconv_ytoh.c vnaconv_ytos.c vnaconv_ytosn.c vnaconv_ytot.c \
	vnaconv_ytoz.c vnaconv_ytozi.c vnaconv_ytozin.c vnaconv_ytozn.c \
	vnaconv_ztoa.c vnaconv_ztob.c vnaconv_ztog.c vnaconv_ztoh.c \
	vnaconv_ztos.c vnaconv_ztosn.c vnaconv_ztot.c vnaconv_ztoy.c \
	vnaconv_ztoyn.c vnaconv_ztozi.c vnaconv_ztozin.c \
	vnadata_internal.h vnadata_add_frequency.c vnadata_alloc.c \
	vnadata_convert.c vnadata_convert_to_fz0.c vnadata_convert_to_z0.c \
	vnadata_get_fz0.c vnadata_get_fz0_vector.c vnadata_get_typename.c \
	vnadata_get_z0.c vnadata_get_z0_vector.c vnadata_set_all_z0.c \
	vnadata_set_fz0.c vnadata_set_fz0_vector.c vnadata_set_z0.c \
	vnadata_set_z0_vector.c \
	vnafile_internal.h vnafile.c vnafile_load.c vnafile_load_native.c \
	vnafile_load_touchstone.c vnafile_save.c \
	vnaproperty.h vnaproperty_internal.h vnaproperty.c


#
# Man pages
#
dist_man_MANS = vnacal.3 vnacal_new.3 vnacal_parameter.3 vnaconv.3 \
	vnadata.3 vnaerr.3 vnafile.3


#
# Examples
#
noinst_PROGRAMS = vnacal-SOLT-example vnacal-TSD-example \
		  vnacal-TRL-example \
		  vnaconv-example vnadata-example vnafile-example

#vnaproperty_example_SOURCES = vnaproperty-example.c
#vnaproperty_example_LDADD = libvna.la
#vnaproperty_example_LDFLAGS = -static

vnacal_SOLT_example_SOURCES = vnacal-SOLT-example.c
vnacal_SOLT_example_LDADD = libvna.la -lyaml -lm
vnacal_SOLT_example_LDFLAGS = -static

vnacal_TSD_example_SOURCES = vnacal-TSD-example.c
vnacal_TSD_example_LDADD = libvna.la -lyaml -lm
vnacal_TSD_example_LDFLAGS = -static

vnacal_TRL_example_SOURCES = vnacal-TRL-example.c
vnacal_TRL_example_LDADD = libvna.la -lyaml -lm
vnacal_TRL_example_LDFLAGS = -static

vnaconv_example_SOURCES = vnaconv-example.c
vnaconv_example_LDADD = libvna.la -lm
vnaconv_example_LDFLAGS = -static

vnadata_example_SOURCES = vnadata-example.c
vnadata_example_LDADD = libvna.la -lm
vnadata_example_LDFLAGS = -static

vnafile_example_SOURCES = vnafile-example.c
vnafile_example_LDADD = libvna.la -lm
vnafile_example_LDFLAGS = -static

dist_doc_DATA = BFCV-4085+_Plus25DegC.s2p MwT-1F.s2p \
	vnacal-SOLT-example.c vnacal-SOLT-example.gnuplot \
	vnacal-TSD-example.c vnacal-TSD-example.gnuplot \
	vnacal-TRL-example.c vnacal-TRL-example.gnuplot \
	vnaconv-example.c vnadata-example.c vnafile-example.c \
	vnafile-example.ts

vnacal-SOLT-example.out: vnacal-SOLT-example
	./vnacal-SOLT-example > "$@"

plot-SOLT: vnacal-SOLT-example.out vnacal-SOLT-example.gnuplot
	$(GNUPLOT) vnacal-SOLT-example.gnuplot

vnacal-TSD-example.out: vnacal-TSD-example
	./vnacal-TSD-example > "$@"

vnacal-TRL-example.out: vnacal-TRL-example
	./vnacal-TRL-example > "$@"

plot-TSD: vnacal-TSD-example.out vnacal-TSD-example.gnuplot
	$(GNUPLOT) vnacal-TSD-example.gnuplot

plot-TRL: vnacal-TRL-example.out vnacal-TRL-example.gnuplot
	$(GNUPLOT) vnacal-TRL-example.gnuplot

plot: plot-SOLT

.3.pdf:
	$(PDFROFF) -t -e -man --no-toc-relocation $< > "$@"

pdfman: vnacal.pdf vnacal_new.pdf vnacal_parameter.pdf vnaconv.pdf \
	vnadata.pdf vnaerr.pdf vnafile.pdf

clean-local:
	rm -f example.vnacal vnacal-SOLT-example.out vnacal-TSD-example.out \
		vnacal-TRL-example.out vnacal.pdf vnaconv.pdf vnadata.pdf \
		vnafile.pdf
