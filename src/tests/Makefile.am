#
# Vector Network Analyzer Library
# Copyright © 2020, 2021 D Scott Guthridge <scott_guthridge@rompromity.net>
#
ACLOCAL_AMFLAGS= -I m4

#
# Add -I for src/
#
AM_CFLAGS= -I$(top_srcdir)/src

#
# vnacaltest library
#
noinst_LIBRARIES = libvnacaltest.a
libvnacaltest_a_SOURCES = vnacaltest.h vnacaltest_calibration.c \
	vnacaltest_measurements.c vnacaltest_properties.c \
	vnacaltest_standards.c vnacaltest_terms.c

#
# Test Data Files
#
dist_noinst_DATA = compat-V2.vnacal

#
# Tests
#
TESTS = \
	test-vnacommon-lu test-vnacommon-mldivide test-vnacommon-mrdivide \
	test-vnacommon-minverse test-vnacommon-qr test-vnacommon-qrsolve \
	test-vnacommon-qrsolve2 \
	test-vnaproperty-scalar test-vnaproperty-list test-vnaproperty-map \
	test-vnaproperty-expr \
	test-vnacal-SOLT test-vnacal-Silvonen16 test-vnacal-random \
	test-vnacal-apply test-vnacal-save-load test-vnacal-TRL \
	test-vnacal-compat-V2 \
	test-vnaconv-2x2 test-vnaconv-3x3
check_PROGRAMS = \
	test-vnacommon-lu test-vnacommon-mldivide test-vnacommon-mrdivide \
	test-vnacommon-minverse test-vnacommon-qr test-vnacommon-qrsolve \
	test-vnacommon-qrsolve2 \
	test-vnaproperty-scalar test-vnaproperty-list test-vnaproperty-map \
	test-vnaproperty-expr \
	test-vnacal-SOLT test-vnacal-Silvonen16 test-vnacal-random \
	test-vnacal-apply test-vnacal-save-load test-vnacal-TRL \
	test-vnacal-compat-V2 \
	test-vnaconv-2x2 test-vnaconv-3x3

test_vnacommon_lu_SOURCES = test.c test-vnacommon-lu.c
test_vnacommon_lu_LDADD = $(top_builddir)/src/libvna.la -lm
test_vnacommon_lu_LDFLAGS = -static

test_vnacommon_mldivide_SOURCES = test.h test.c test-vnacommon-mldivide.c
test_vnacommon_mldivide_LDADD = $(top_builddir)/src/libvna.la -lm
test_vnacommon_mldivide_LDFLAGS = -static

test_vnacommon_mrdivide_SOURCES = test.h test.c test-vnacommon-mrdivide.c
test_vnacommon_mrdivide_LDADD = $(top_builddir)/src/libvna.la -lm
test_vnacommon_mrdivide_LDFLAGS = -static

test_vnacommon_minverse_SOURCES = test.h test.c test-vnacommon-minverse.c
test_vnacommon_minverse_LDADD = $(top_builddir)/src/libvna.la -lm
test_vnacommon_minverse_LDFLAGS = -static

test_vnacommon_qr_SOURCES = test.h test.c test-vnacommon-qr.c
test_vnacommon_qr_LDADD = $(top_builddir)/src/libvna.la -lm
test_vnacommon_qr_LDFLAGS = -static

test_vnacommon_qrsolve_SOURCES = test.h test.c test-vnacommon-qrsolve.c
test_vnacommon_qrsolve_LDADD = $(top_builddir)/src/libvna.la -lm
test_vnacommon_qrsolve_LDFLAGS = -static

test_vnacommon_qrsolve2_SOURCES = test.h test.c test-vnacommon-qrsolve2.c
test_vnacommon_qrsolve2_LDADD = $(top_builddir)/src/libvna.la -lm
test_vnacommon_qrsolve2_LDFLAGS = -static

test_vnaproperty_scalar_SOURCES = test.h test.c test-vnaproperty-scalar.c
test_vnaproperty_scalar_LDADD = $(top_builddir)/src/libvna.la
test_vnaproperty_scalar_LDFLAGS = -static

test_vnaproperty_list_SOURCES = test.h test.c test-vnaproperty-list.c
test_vnaproperty_list_LDADD = $(top_builddir)/src/libvna.la
test_vnaproperty_list_LDFLAGS = -static

test_vnaproperty_map_SOURCES = test.h test.c test-vnaproperty-map.c
test_vnaproperty_map_LDADD = $(top_builddir)/src/libvna.la
test_vnaproperty_map_LDFLAGS = -static

test_vnaproperty_expr_SOURCES = test.h test.c test-vnaproperty-expr.c
test_vnaproperty_expr_LDADD = $(top_builddir)/src/libvna.la
test_vnaproperty_expr_LDFLAGS = -static

test_vnaconv_2x2_SOURCES = test.h test.c test-vnaconv-2x2.c
test_vnaconv_2x2_LDADD = $(top_builddir)/src/libvna.la -lm
test_vnaconv_2x2_LDFLAGS = -static

test_vnaconv_3x3_SOURCES = test.h test.c test-vnaconv-3x3.c
test_vnaconv_3x3_LDADD = $(top_builddir)/src/libvna.la -lm
test_vnaconv_3x3_LDFLAGS = -static

#vnadata_test_SOURCES = test.h test.c vnadata-test.c
#vnadata_test_LDADD = $(top_builddir)/src/libvna.la -lm
#vnadata_test_LDFLAGS = -static

test_vnacal_SOLT_SOURCES = test.h test.c test-vnacal-SOLT.c
test_vnacal_SOLT_LDADD = libvnacaltest.a $(top_builddir)/src/libvna.la -lyaml -lm
test_vnacal_SOLT_LDFLAGS = -static

test_vnacal_Silvonen16_SOURCES = test.h test.c test-vnacal-Silvonen16.c
test_vnacal_Silvonen16_LDADD = libvnacaltest.a $(top_builddir)/src/libvna.la -lyaml -lm
test_vnacal_Silvonen16_LDFLAGS = -static

test_vnacal_random_SOURCES = test.h test.c test-vnacal-random.c
test_vnacal_random_LDADD = libvnacaltest.a $(top_builddir)/src/libvna.la -lyaml -lm
test_vnacal_random_LDFLAGS = -static

test_vnacal_apply_SOURCES = test.h test.c test-vnacal-apply.c
test_vnacal_apply_LDADD = libvnacaltest.a $(top_builddir)/src/libvna.la -lyaml -lm
test_vnacal_apply_LDFLAGS = -static

test_vnacal_save_load_SOURCES = test.h test.c test-vnacal-save-load.c
test_vnacal_save_load_LDADD = libvnacaltest.a $(top_builddir)/src/libvna.la \
	-lyaml -lm
test_vnacal_save_load_LDFLAGS = -static

test_vnacal_TRL_SOURCES = test.c test-vnacal-TRL.c
test_vnacal_TRL_LDADD = libvnacaltest.a $(top_builddir)/src/libvna.la \
	-lyaml -lm
test_vnacal_TRL_LDFLAGS = -static

test_vnacal_compat_V2_SOURCES = test.h test.c test-vnacal-compat-V2.c
test_vnacal_compat_V2_LDADD = libvnacaltest.a $(top_builddir)/src/libvna.la \
	-lyaml -lm
test_vnacal_compat_V2_LDFLAGS = -static

#vnafile_test_SOURCES = test.h test.c vnafile-test.c
#vnafile_test_LDADD = $(top_builddir)/src/libvna.la -lm
#vnafile_test_LDFLAGS = -static

clean-local:
	rm -f test-vnacal.vnacal
