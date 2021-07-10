#include "archdep.h"

#include <assert.h>
#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnaconv.h"
#include "vnadata.h"
#include "libt.h"
#include "libt_vnadata.h"

void libt_vnadata_convert(const double complex *in, double complex *out,
	const double complex *z0, int rows, int columns,
	vnadata_parameter_type_t old_type, vnadata_parameter_type_t new_type)
{
    typedef double complex row2_t[];

    switch (old_type) {
    case VPT_S:
	switch (new_type) {
	case VPT_S:
	    if (out != in) {
		(void)memcpy((void *)out, (void *)in,
			rows * columns * sizeof(double complex));
	    }
	    return;
	case VPT_Z:
	    assert(rows == columns);
	    vnaconv_stozn(in, out, z0, rows);
	    return;
	case VPT_Y:
	    assert(rows == columns);
	    vnaconv_stoyn(in, out, z0, rows);
	    return;
	case VPT_H:
	    assert(rows == 2 && columns == 2);
	    vnaconv_stoh((const row2_t *)in, (row2_t *)out, z0);
	    return;
	case VPT_G:
	    assert(rows == 2 && columns == 2);
	    vnaconv_stog((const row2_t *)in, (row2_t *)out, z0);
	    return;
	case VPT_A:
	    assert(rows == 2 && columns == 2);
	    vnaconv_stoa((const row2_t *)in, (row2_t *)out, z0);
	    return;
	case VPT_B:
	    assert(rows == 2 && columns == 2);
	    vnaconv_stob((const row2_t *)in, (row2_t *)out, z0);
	    return;
	case VPT_T:
	    assert(rows == 2 && columns == 2);
	    vnaconv_stot((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_ZIN:
	    vnaconv_stozimn(in, out, z0, rows, columns);
	    return;
	default:
	    break;
	}
	break;

    case VPT_Z:
	switch (new_type) {
	case VPT_S:
	    assert(rows == columns);
	    vnaconv_ztosn(in, out, z0, rows);
	    return;
	case VPT_Z:
	    assert(rows == columns);
	    if (out != in) {
		(void)memcpy((void *)out, (void *)in,
			rows * columns * sizeof(double complex));
	    }
	    return;
	case VPT_Y:
	    assert(rows == 2 && columns == 2);
	    vnaconv_ztoyn(in, out, rows);
	    return;
	case VPT_H:
	    assert(rows == 2 && columns == 2);
	    vnaconv_ztoh((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_G:
	    assert(rows == 2 && columns == 2);
	    vnaconv_ztog((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_A:
	    assert(rows == 2 && columns == 2);
	    vnaconv_ztoa((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_B:
	    assert(rows == 2 && columns == 2);
	    vnaconv_ztob((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_T:
	    assert(rows == 2 && columns == 2);
	    vnaconv_ztot((const row2_t *)in, (row2_t *)out, z0);
	    return;
	case VPT_ZIN:
	    assert(rows == columns);
	    vnaconv_ztozin(in, out, z0, rows);
	    return;
	default:
	    break;
	}
	break;

    case VPT_Y:
	switch (new_type) {
	case VPT_S:
	    assert(rows == columns);
	    vnaconv_ytosn(in, out, z0, rows);
	    return;
	case VPT_Z:
	    assert(rows == columns);
	    vnaconv_ytozn(in, out, rows);
	    return;
	case VPT_Y:
	    assert(rows == columns);
	    if (out != in) {
		(void)memcpy((void *)out, (void *)in,
			rows * columns * sizeof(double complex));
	    }
	    return;
	case VPT_H:
	    assert(rows == 2 && columns == 2);
	    vnaconv_ytoh((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_G:
	    assert(rows == 2 && columns == 2);
	    vnaconv_ytog((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_A:
	    assert(rows == 2 && columns == 2);
	    vnaconv_ytoa((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_B:
	    assert(rows == 2 && columns == 2);
	    vnaconv_ytob((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_T:
	    assert(rows == 2 && columns == 2);
	    vnaconv_ytot((const row2_t *)in, (row2_t *)out, z0);
	    return;
	case VPT_ZIN:
	    assert(rows == columns);
	    vnaconv_ytozin(in, out, z0, MIN(rows, columns));
	    return;
	default:
	    break;
	}
	break;

    case VPT_H:
	assert(rows == 2 && columns == 2);
	switch (new_type) {
	case VPT_S:
	    vnaconv_htos((const row2_t *)in, (row2_t *)out, z0);
	    return;
	case VPT_Z:
	    vnaconv_htoz((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_Y:
	    vnaconv_htoy((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_H:
	    if (out != in) {
		(void)memcpy((void *)out, (void *)in,
			rows * columns * sizeof(double complex));
	    }
	    return;
	case VPT_G:
	    vnaconv_htog((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_A:
	    vnaconv_htoa((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_B:
	    vnaconv_htob((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_T:
	    vnaconv_htot((const row2_t *)in, (row2_t *)out, z0);
	    return;
	case VPT_ZIN:
	    vnaconv_htozi((const row2_t *)in, out, z0);
	    return;
	default:
	    break;
	}
	break;

    case VPT_G:
	assert(rows == 2 && columns == 2);
	switch (new_type) {
	case VPT_S:
	    vnaconv_gtos((const row2_t *)in, (row2_t *)out, z0);
	    return;
	case VPT_Z:
	    vnaconv_gtoz((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_Y:
	    vnaconv_gtoy((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_H:
	    vnaconv_gtoh((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_G:
	    if (out != in) {
		(void)memcpy((void *)out, (void *)in,
			rows * columns * sizeof(double complex));
	    }
	    return;
	case VPT_A:
	    vnaconv_gtoa((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_B:
	    vnaconv_gtob((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_T:
	    vnaconv_gtot((const row2_t *)in, (row2_t *)out, z0);
	    return;
	case VPT_ZIN:
	    vnaconv_gtozi((const row2_t *)in, out, z0);
	    return;
	default:
	    break;
	}
	break;

    case VPT_A:
	assert(rows == 2 && columns == 2);
	switch (new_type) {
	case VPT_S:
	    vnaconv_atos((const row2_t *)in, (row2_t *)out, z0);
	    return;
	case VPT_Z:
	    vnaconv_atoz((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_Y:
	    vnaconv_atoy((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_H:
	    vnaconv_atoh((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_G:
	    vnaconv_atog((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_A:
	    if (out != in) {
		(void)memcpy((void *)out, (void *)in,
			rows * columns * sizeof(double complex));
	    }
	    return;
	case VPT_B:
	    vnaconv_atob((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_T:
	    vnaconv_atot((const row2_t *)in, (row2_t *)out, z0);
	    return;
	case VPT_ZIN:
	    vnaconv_atozi((const row2_t *)in, out, z0);
	    return;
	default:
	    break;
	}
	break;

    case VPT_B:
	assert(rows == 2 && columns == 2);
	switch (new_type) {
	case VPT_S:
	    vnaconv_btos((const row2_t *)in, (row2_t *)out, z0);
	    return;
	case VPT_Z:
	    vnaconv_btoz((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_Y:
	    vnaconv_btoy((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_H:
	    vnaconv_btoh((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_G:
	    vnaconv_btog((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_A:
	    vnaconv_btoa((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_B:
	    if (out != in) {
		(void)memcpy((void *)out, (void *)in,
			rows * columns * sizeof(double complex));
	    }
	    return;
	case VPT_T:
	    vnaconv_btot((const row2_t *)in, (row2_t *)out, z0);
	    return;
	case VPT_ZIN:
	    vnaconv_btozi((const row2_t *)in, out, z0);
	    return;
	default:
	    break;
	}
	break;

    case VPT_T:
	assert(rows == 2 && columns == 2);
	switch (new_type) {
	case VPT_S:
	    vnaconv_ttos((const row2_t *)in, (row2_t *)out);
	    return;
	case VPT_Z:
	    vnaconv_ttoz((const row2_t *)in, (row2_t *)out, z0);
	    return;
	case VPT_Y:
	    vnaconv_ttoy((const row2_t *)in, (row2_t *)out, z0);
	    return;
	case VPT_H:
	    vnaconv_ttoh((const row2_t *)in, (row2_t *)out, z0);
	    return;
	case VPT_G:
	    vnaconv_ttog((const row2_t *)in, (row2_t *)out, z0);
	    return;
	case VPT_A:
	    vnaconv_ttoa((const row2_t *)in, (row2_t *)out, z0);
	    return;
	case VPT_B:
	    vnaconv_ttob((const row2_t *)in, (row2_t *)out, z0);
	    return;
	case VPT_T:
	    if (out != in) {
		(void)memcpy((void *)out, (void *)in,
			rows * columns * sizeof(double complex));
	    }
	    return;
	case VPT_ZIN:
	    vnaconv_ttozi((const row2_t *)in, out, z0);
	    return;
	default:
	    break;
	}
	break;

    case VPT_ZIN:
	assert(rows == 1);
	switch (new_type) {
	case VPT_ZIN:
	    if (out != in) {
		(void)memcpy((void *)out, (void *)in,
			columns * sizeof(double complex));
	    }
	    return;
	default:
	    break;
	}

    default:
	break;
    }
    libt_error("unexpected conversion: %s -> %s\n",
	    vnadata_get_type_name(old_type),
	    vnadata_get_type_name(new_type));
    return;
}
