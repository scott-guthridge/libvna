#
# Vector Network Analyzer Library
# Copyright © 2020, 2021 D Scott Guthridge <scott_guthridge@rompromity.net>
#
ACLOCAL_AMFLAGS= -I m4
GNUPLOT= gnuplot
PDFROFF= pdfroff
SUFFIXES= .3 .pdf

#
# Library
#
include_HEADERS = vnacal.h vnaconv.h vnadata.h vnaerr.h vnafile.h
lib_LTLIBRARIES = libvna.la
libvna_la_SOURCES = archdep.h archdep.c vnacal_internal.h \
	vnacal_delete_parameter.c vnacal_apply.c vnacal_new_add_common.c \
	vnacal_new_parameter.c vnacal_new_solve.c \
	vnacal_new.c vnacal_create.c vnacal_calibration.c \
	vnacal_make_correlated_parameter.c vnacal_make_scalar_parameter.c \
	vnacal_make_unknown_parameter.c vnacal_make_vector_parameter.c \
	vnacal_add_calibration.c vnacal_delete_calibration.c \
	vnacal_layout.h vnacal_layout.c vnacal_get_parameter_value.c \
	vnacal_error.c vnacal_find_calibration.c vnacal_free.c vnacal_get.c \
	vnacal_load.c vnacal_property.c vnacal_rfi.c vnacal_save.c \
	vnacal_set_dprecision.c vnacal_set_fprecision.c vnacal_parameter.c \
	vnacommon_internal.h vnacommon_lu.c vnacommon_minverse.c \
	vnacommon_mldivide.c vnacommon_mrdivide.c \
	vnaerr_internal.h vnaerr_verror.c \
	vnacommon_qrd.c vnacommon_qrsolve.c \
	vnaconv_internal.h vnaconv_a2b.c vnaconv_a2g.c vnaconv_a2h.c \
	vnaconv_a2s.c vnaconv_a2t.c vnaconv_a2y.c vnaconv_a2z.c \
	vnaconv_a2zi.c vnaconv_b2a.c vnaconv_b2g.c vnaconv_b2h.c \
	vnaconv_b2s.c vnaconv_b2t.c vnaconv_b2y.c vnaconv_b2z.c \
	vnaconv_b2zi.c vnaconv_g2a.c vnaconv_g2b.c vnaconv_g2h.c \
	vnaconv_g2s.c vnaconv_g2t.c vnaconv_g2y.c vnaconv_g2z.c \
	vnaconv_g2zi.c vnaconv_h2a.c vnaconv_h2b.c vnaconv_h2g.c \
	vnaconv_h2s.c vnaconv_h2t.c vnaconv_h2y.c vnaconv_h2z.c \
	vnaconv_h2zi.c vnaconv_s2a.c vnaconv_s2b.c vnaconv_s2g.c \
	vnaconv_s2h.c vnaconv_s2t.c vnaconv_s2y.c vnaconv_s2yn.c \
	vnaconv_s2z.c vnaconv_s2zi.c vnaconv_s2zimn.c vnaconv_s2zin.c \
	vnaconv_s2zn.c vnaconv_t2a.c vnaconv_t2b.c vnaconv_t2g.c \
	vnaconv_t2h.c vnaconv_t2s.c vnaconv_t2y.c vnaconv_t2z.c \
	vnaconv_t2zi.c vnaconv_y2a.c vnaconv_y2b.c vnaconv_y2g.c \
	vnaconv_y2h.c vnaconv_y2s.c vnaconv_y2sn.c vnaconv_y2t.c \
	vnaconv_y2z.c vnaconv_y2zi.c vnaconv_y2zin.c vnaconv_y2zn.c \
	vnaconv_z2a.c vnaconv_z2b.c vnaconv_z2g.c vnaconv_z2h.c \
	vnaconv_z2s.c vnaconv_z2sn.c vnaconv_z2t.c vnaconv_z2y.c \
	vnaconv_z2yn.c vnaconv_z2zi.c vnaconv_z2zin.c \
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
# Tests
#
TESTS = vnacommon-test vnaproperty-test vnaconv-test vnacal-test
check_PROGRAMS = vnacommon-test vnaproperty-test vnaconv-test vnacal-test

vnacommon_test_SOURCES = vnacommon-test.c
vnacommon_test_LDADD = libvna.la -lm
vnacommon_test_LDFLAGS = -static

vnaproperty_test_SOURCES = vnaproperty-test.c
vnaproperty_test_LDADD = libvna.la
vnaproperty_test_LDFLAGS = -static

vnaconv_test_SOURCES = vnaconv-test.c
vnaconv_test_LDADD = libvna.la -lm
vnaconv_test_LDFLAGS = -static

#vnadata_test_SOURCES = vnadata-test.c
#vnadata_test_LDADD = libvna.la -lm
#vnadata_test_LDFLAGS = -static

vnacal_test_SOURCES = vnacal-test.c
vnacal_test_LDADD = libvna.la -lyaml -lm
vnacal_test_LDFLAGS = -static

#vnafile_test_SOURCES = vnafile-test.c
#vnafile_test_LDADD = libvna.la -lm
#vnafile_test_LDFLAGS = -static


#
# Examples
#
noinst_PROGRAMS = vnacal-example vnaconv-example vnadata-example \
	vnafile-example

#vnaproperty_example_SOURCES = vnaproperty-example.c
#vnaproperty_example_LDADD = libvna.la
#vnaproperty_example_LDFLAGS = -static

vnacal_example_SOURCES = vnacal-example.c
vnacal_example_LDADD = libvna.la -lyaml -lm
vnacal_example_LDFLAGS = -static

vnaconv_example_SOURCES = vnaconv-example.c
vnaconv_example_LDADD = libvna.la -lm
vnaconv_example_LDFLAGS = -static

vnadata_example_SOURCES = vnadata-example.c
vnadata_example_LDADD = libvna.la -lm
vnadata_example_LDFLAGS = -static

vnafile_example_SOURCES = vnafile-example.c
vnafile_example_LDADD = libvna.la -lm
vnafile_example_LDFLAGS = -static

dist_doc_DATA = vnacal-example.c vnacal-example.gnuplot vnaconv-example.c \
	vnadata-example.c vnafile-example.c vnafile-example.ts

vnacal-example.out: vnacal-example
	./vnacal-example > "$@"

plot: vnacal-example.out vnacal-example.gnuplot
	$(GNUPLOT) vnacal-example.gnuplot

.3.pdf:
	$(PDFROFF) -t -e -man --no-toc-relocation $< > "$@"

pdfman: vnacal.pdf vnacal_new.pdf vnacal_parameter.pdf vnaconv.pdf \
	vnadata.pdf vnaerr.pdf vnafile.pdf

clean-local:
	rm -f example.vnacal vnacal-example.out vnacal-test.vnacal \
		vnacal.pdf vnaconv.pdf vnadata.pdf vnafile.pdf
