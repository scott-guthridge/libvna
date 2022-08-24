#
# Vector Network Analyzer Library
# Copyright © 2020-2022 D Scott Guthridge <scott_guthridge@rompromity.net>
#
ACLOCAL_AMFLAGS= -I m4
GNUPLOT= gnuplot
PDFROFF= pdfroff
SUFFIXES= .3 .pdf

SUBDIRS= . tests

#
# Library
#
include_HEADERS = vnacal.h vnaconv.h vnadata.h vnaerr.h vnaproperty.h
lib_LTLIBRARIES = libvna.la
libvna_la_SOURCES = archdep.h archdep.c vnacal_internal.h \
	vnacal_new_internal.h \
	vnacal_add_calibration.c vnacal_apply.c vnacal_calibration.c \
	vnacal_create.c vnacal_delete_calibration.c \
	vnacal_delete_parameter.c vnacal_error.c vnacal_find_calibration.c \
	vnacal_free.c vnacal_get.c vnacal_get_parameter_value.c \
	vnacal_layout.c vnacal_layout.h vnacal_load.c \
	vnacal_make_correlated_parameter.c vnacal_make_scalar_parameter.c \
	vnacal_make_unknown_parameter.c vnacal_make_vector_parameter.c \
	vnacal_name_to_type.c vnacal_new_add_common.c vnacal_new.c \
	vnacal_new_parameter.c vnacal_new_set_m_error.c \
	vnacal_new_set_p_tolerance.c vnacal_new_solve.c \
	vnacal_parameter.c vnacal_property.c vnacal_rfi.c \
	vnacal_save.c vnacal_set_dprecision.c vnacal_set_fprecision.c \
	vnacal_type_to_name.c \
	vnacommon_internal.h vnacommon_lu.c vnacommon_mmultiply.c \
	vnacommon_minverse.c vnacommon_mldivide.c vnacommon_mrdivide.c \
	vnacommon_qrd.c vnacommon_qr.c vnacommon_qrsolve.c \
	vnacommon_qrsolve2.c vnacommon_spline.c \
	vnaerr_internal.h vnaerr_verror.c \
	vnaconv_atob.c vnaconv_atog.c vnaconv_atoh.c vnaconv_atos.c \
	vnaconv_atot.c vnaconv_atou.c vnaconv_atoy.c vnaconv_atoz.c \
	vnaconv_atozi.c vnaconv_btoa.c vnaconv_btog.c vnaconv_btoh.c \
	vnaconv_btos.c vnaconv_btot.c vnaconv_btou.c vnaconv_btoy.c \
	vnaconv_btoz.c vnaconv_btozi.c vnaconv_gtoa.c vnaconv_gtob.c \
	vnaconv_gtoh.c vnaconv_gtos.c vnaconv_gtot.c vnaconv_gtou.c \
	vnaconv_gtoy.c vnaconv_gtoz.c vnaconv_gtozi.c vnaconv_htoa.c \
	vnaconv_htob.c vnaconv_htog.c vnaconv_htos.c vnaconv_htot.c \
	vnaconv_htou.c vnaconv_htoy.c vnaconv_htoz.c vnaconv_htozi.c \
	vnaconv_internal.h vnaconv_stoa.c vnaconv_stob.c vnaconv_stog.c \
	vnaconv_stoh.c vnaconv_stot.c vnaconv_stou.c vnaconv_stoy.c \
	vnaconv_stoyn.c vnaconv_stoz.c vnaconv_stozi.c \
	vnaconv_stozin.c vnaconv_stozn.c vnaconv_ttoa.c vnaconv_ttob.c \
	vnaconv_ttog.c vnaconv_ttoh.c vnaconv_ttos.c vnaconv_ttou.c \
	vnaconv_ttoy.c vnaconv_ttoz.c vnaconv_ttozi.c vnaconv_utoa.c \
	vnaconv_utob.c vnaconv_utog.c vnaconv_utoh.c vnaconv_utos.c \
	vnaconv_utot.c vnaconv_utoy.c vnaconv_utoz.c vnaconv_utozi.c \
	vnaconv_ytoa.c vnaconv_ytob.c vnaconv_ytog.c vnaconv_ytoh.c \
	vnaconv_ytos.c vnaconv_ytosn.c vnaconv_ytot.c vnaconv_ytou.c \
	vnaconv_ytoz.c vnaconv_ytozi.c vnaconv_ytozin.c vnaconv_ytozn.c \
	vnaconv_ztoa.c vnaconv_ztob.c vnaconv_ztog.c vnaconv_ztoh.c \
	vnaconv_ztos.c vnaconv_ztosn.c vnaconv_ztot.c vnaconv_ztou.c \
	vnaconv_ztoy.c vnaconv_ztoyn.c vnaconv_ztozi.c vnaconv_ztozin.c \
	vnadata_add_frequency.c vnadata_alloc.c vnadata_convert.c \
	vnadata_convert_to_fz0.c vnadata_convert_to_z0.c \
	vnadata_format_to_name.c vnadata_get_dprecision.c \
	vnadata_get_filetype.c vnadata_get_format.c vnadata_get_fprecision.c \
	vnadata_get_fz0.c vnadata_get_fz0_vector.c vnadata_get_type_name.c \
	vnadata_get_z0.c vnadata_get_z0_vector.c vnadata_has_fz0.c \
	vnadata_internal.h \
	vnadata_load.c vnadata_load_npd.c vnadata_load_touchstone.c \
	vnadata_parse_filename.c vnadata_save.c vnadata_set_all_z0.c \
	vnadata_set_dprecision.c vnadata_set_filetype.c vnadata_set_format.c \
	vnadata_set_fprecision.c vnadata_set_fz0.c vnadata_set_fz0_vector.c \
	vnadata_set_simple_format.c vnadata_set_z0.c vnadata_set_z0_vector.c \
	vnadata_update_format_string.c vnaproperty_internal.h vnaproperty.c \
	vnaproperty_import_yaml_from_string.c \
	vnaproperty_import_yaml_from_file.c \
	vnaproperty_export_yaml_to_file.c


#
# Man pages
#
dist_man_MANS = vnacal.3 vnacal_new.3 vnacal_parameter.3 vnaconv.3 \
	vnadata.3 vnaerr.3 vnaproperty.3

#
# Examples
#
noinst_PROGRAMS = vnaproperty-example vnacal-SOLT-example vnacal-TSD-example \
		  vnacal-TRL-example \
		  vnaconv-example vnadata-example convert-parameters

vnaproperty_example_SOURCES = vnaproperty-example.c
vnaproperty_example_LDADD = libvna.la -lyaml
vnaproperty_example_LDFLAGS = -static

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

convert_parameters_SOURCES = convert-parameters.c
convert_parameters_LDADD = libvna.la -lm
convert_parameters_LDFLAGS = -static

examplesdir = $(docdir)/examples

dist_examples_DATA = BFCV-4085+_Plus25DegC.s2p MwT-1F.s2p \
	vnacal-SOLT-example.c vnacal-SOLT-example.gnuplot \
	vnacal-TSD-example.c vnacal-TSD-example.gnuplot \
	vnacal-TRL-example.c vnacal-TRL-example.gnuplot \
	vnaconv-example.c vnadata-example.c vnadata-example.npd \
	convert-parameters.c convert-parameters.ts

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
	vnadata.pdf vnaerr.pdf vnaproperty.pdf

clean-local:
	rm -f example.vnacal *.vnacal vnacal-SOLT-example.out \
		vnacal-TSD-example.out vnacal-TRL-example.out \
		vnacal.pdf vnacal_new.pdf vnacal_parameter.pdf \
		vnaconv.pdf vnadata.pdf vnaerr.pdf vnaproperty.pdf
