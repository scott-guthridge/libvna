#
# Vector Network Analyzer Library
# Copyright © 2020-2024 D Scott Guthridge <scott_guthridge@rompromity.net>
#
ACLOCAL_AMFLAGS= -I m4

#
# Add -I for src/
#
AM_CFLAGS= -I$(top_srcdir)/src

#
# vnacaltest library
#
noinst_LIBRARIES = libt.a
libt_a_SOURCES = libt.h libt.c libt_crand.h libt_crand.c \
	libt_vnacal.h libt_vnacal_calibration.c \
	libt_vnacal_measurements.c libt_vnacal_properties.c \
	libt_vnacal_standards.c libt_vnacal_terms.c \
	libt_vnadata_convert.c libt_vnadata_create.c libt_vnadata.h

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
	test-vnacal-apply test-vnacal-save-load test-vnacal-v-matrices \
	test-vnacal-pvalue test-vnacal-standards \
	test-vnacal-TRL test-vnacal-Van-hamme test-vnacal-compat-V2 \
	test-vnaconv-2x2 test-vnaconv-3x3 \
	test-vnaconv-renormalize-2x2 test-vnaconv-renormalize-NxN \
	test-vnadata-basic test-vnadata-save-load-convert \
	test-vnadata-rconvert
check_PROGRAMS = \
	test-vnacommon-lu test-vnacommon-mldivide test-vnacommon-mrdivide \
	test-vnacommon-minverse test-vnacommon-qr test-vnacommon-qrsolve \
	test-vnacommon-qrsolve2 \
	test-vnaproperty-scalar test-vnaproperty-list test-vnaproperty-map \
	test-vnaproperty-expr \
	test-vnacal-SOLT test-vnacal-Silvonen16 test-vnacal-random \
	test-vnacal-apply test-vnacal-save-load test-vnacal-v-matrices \
	test-vnacal-pvalue test-vnacal-standards \
	test-vnacal-TRL test-vnacal-Van-hamme test-vnacal-compat-V2 \
	test-vnaconv-2x2 test-vnaconv-3x3 \
	test-vnaconv-renormalize-2x2 test-vnaconv-renormalize-NxN \
	test-vnadata-basic test-vnadata-save-load-convert \
	test-vnadata-rconvert

test_vnacommon_lu_SOURCES = test-vnacommon-lu.c
test_vnacommon_lu_LDADD = libt.a $(top_builddir)/src/libvna.la -lm
test_vnacommon_lu_LDFLAGS = -static

test_vnacommon_mldivide_SOURCES = test-vnacommon-mldivide.c
test_vnacommon_mldivide_LDADD = libt.a $(top_builddir)/src/libvna.la -lm
test_vnacommon_mldivide_LDFLAGS = -static

test_vnacommon_mrdivide_SOURCES = test-vnacommon-mrdivide.c
test_vnacommon_mrdivide_LDADD = libt.a $(top_builddir)/src/libvna.la -lm
test_vnacommon_mrdivide_LDFLAGS = -static

test_vnacommon_minverse_SOURCES = test-vnacommon-minverse.c
test_vnacommon_minverse_LDADD = libt.a $(top_builddir)/src/libvna.la -lm
test_vnacommon_minverse_LDFLAGS = -static

test_vnacommon_qr_SOURCES = test-vnacommon-qr.c
test_vnacommon_qr_LDADD = libt.a $(top_builddir)/src/libvna.la -lm
test_vnacommon_qr_LDFLAGS = -static

test_vnacommon_qrsolve_SOURCES = test-vnacommon-qrsolve.c
test_vnacommon_qrsolve_LDADD = libt.a $(top_builddir)/src/libvna.la -lm
test_vnacommon_qrsolve_LDFLAGS = -static

test_vnacommon_qrsolve2_SOURCES = test-vnacommon-qrsolve2.c
test_vnacommon_qrsolve2_LDADD = libt.a $(top_builddir)/src/libvna.la -lm
test_vnacommon_qrsolve2_LDFLAGS = -static

test_vnaproperty_scalar_SOURCES = test-vnaproperty-scalar.c
test_vnaproperty_scalar_LDADD = libt.a $(top_builddir)/src/libvna.la -lyaml
test_vnaproperty_scalar_LDFLAGS = -static

test_vnaproperty_list_SOURCES = test-vnaproperty-list.c
test_vnaproperty_list_LDADD = libt.a $(top_builddir)/src/libvna.la -lyaml
test_vnaproperty_list_LDFLAGS = -static

test_vnaproperty_map_SOURCES = test-vnaproperty-map.c
test_vnaproperty_map_LDADD = libt.a $(top_builddir)/src/libvna.la -lyaml
test_vnaproperty_map_LDFLAGS = -static

test_vnaproperty_expr_SOURCES = test-vnaproperty-expr.c
test_vnaproperty_expr_LDADD = libt.a $(top_builddir)/src/libvna.la -lyaml
test_vnaproperty_expr_LDFLAGS = -static

test_vnaconv_2x2_SOURCES = test-vnaconv-2x2.c
test_vnaconv_2x2_LDADD = libt.a $(top_builddir)/src/libvna.la -lm
test_vnaconv_2x2_LDFLAGS = -static

test_vnaconv_3x3_SOURCES = test-vnaconv-3x3.c
test_vnaconv_3x3_LDADD = libt.a $(top_builddir)/src/libvna.la -lm
test_vnaconv_3x3_LDFLAGS = -static

test_vnaconv_renormalize_2x2_SOURCES = test-vnaconv-renormalize-2x2.c
test_vnaconv_renormalize_2x2_LDADD = libt.a $(top_builddir)/src/libvna.la -lm
test_vnaconv_renormalize_2x2_LDFLAGS = -static

test_vnaconv_renormalize_NxN_SOURCES = test-vnaconv-renormalize-NxN.c
test_vnaconv_renormalize_NxN_LDADD = libt.a $(top_builddir)/src/libvna.la -lm
test_vnaconv_renormalize_NxN_LDFLAGS = -static

#vnadata_test_SOURCES = libt.h libt.c vnadata-libt.c
#vnadata_test_LDADD = $(top_builddir)/src/libvna.la -lm
#vnadata_test_LDFLAGS = -static

test_vnacal_SOLT_SOURCES = test-vnacal-SOLT.c
test_vnacal_SOLT_LDADD = libt.a $(top_builddir)/src/libvna.la -lyaml -lm
test_vnacal_SOLT_LDFLAGS = -static

test_vnacal_Silvonen16_SOURCES = test-vnacal-Silvonen16.c
test_vnacal_Silvonen16_LDADD = libt.a $(top_builddir)/src/libvna.la -lyaml -lm
test_vnacal_Silvonen16_LDFLAGS = -static

test_vnacal_random_SOURCES = test-vnacal-random.c
test_vnacal_random_LDADD = libt.a $(top_builddir)/src/libvna.la -lyaml -lm
test_vnacal_random_LDFLAGS = -static

test_vnacal_apply_SOURCES = test-vnacal-apply.c
test_vnacal_apply_LDADD = libt.a $(top_builddir)/src/libvna.la -lyaml -lm
test_vnacal_apply_LDFLAGS = -static

test_vnacal_save_load_SOURCES = test-vnacal-save-load.c
test_vnacal_save_load_LDADD = libt.a $(top_builddir)/src/libvna.la \
	-lyaml -lm
test_vnacal_save_load_LDFLAGS = -static

test_vnacal_v_matrices_SOURCES = test-vnacal-v-matrices.c
test_vnacal_v_matrices_LDADD = libt.a $(top_builddir)/src/libvna.la \
	-lyaml -lm
test_vnacal_v_matrices_LDFLAGS = -static

test_vnacal_pvalue_SOURCES = test-vnacal-pvalue.c
test_vnacal_pvalue_LDADD = libt.a $(top_builddir)/src/libvna.la \
	-lyaml -lm
test_vnacal_pvalue_LDFLAGS = -static

test_vnacal_TRL_SOURCES = test-vnacal-TRL.c
test_vnacal_TRL_LDADD = libt.a $(top_builddir)/src/libvna.la \
	-lyaml -lm
test_vnacal_TRL_LDFLAGS = -static

test_vnacal_Van_hamme_SOURCES = test-vnacal-Van-hamme.c
test_vnacal_Van_hamme_LDADD = libt.a $(top_builddir)/src/libvna.la \
	-lyaml -lm
test_vnacal_Van_hamme_LDFLAGS = -static

test_vnacal_compat_V2_SOURCES = test-vnacal-compat-V2.c
test_vnacal_compat_V2_LDADD = libt.a $(top_builddir)/src/libvna.la \
	-lyaml -lm
test_vnacal_compat_V2_LDFLAGS = -static

test_vnacal_standards_SOURCES = libt.h libt.c \
	test-vnacal-standards.c
test_vnacal_standards_LDADD = libt.a $(top_builddir)/src/libvna.la \
	-lyaml -lm
test_vnacal_standards_LDFLAGS = -static

test_vnadata_basic_SOURCES = libt.h libt.c test-vnadata-basic.c
test_vnadata_basic_LDADD = libt.a $(top_builddir)/src/libvna.la -lyaml -lm
test_vnadata_basic_LDFLAGS = -static

test_vnadata_save_load_convert_SOURCES = libt.h libt.c \
	test-vnadata-save-load-convert.c
test_vnadata_save_load_convert_LDADD = libt.a $(top_builddir)/src/libvna.la \
	-lyaml -lm
test_vnadata_save_load_convert_LDFLAGS = -static

test_vnadata_rconvert_SOURCES = libt.h libt.c \
	test-vnadata-rconvert.c
test_vnadata_rconvert_LDADD = libt.a $(top_builddir)/src/libvna.la \
	-lyaml -lm
test_vnadata_rconvert_LDFLAGS = -static

clean-local:
	rm -f test-vnacal.vnacal test-vnadata.s1p test-vnadata.s2p \
		test-vnadata.s3p test-vnadata.s4p test-vnadata.ts \
		test-vnadata.npd

