/*
 * Vector Network Analyzer Library
 * Copyright © 2020, 2021 D Scott Guthridge <scott_guthridge@rompromity.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "archdep.h"

#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "vnaconv.h"
#include "vnadata_internal.h"


/*
 * conversion_group: a bitwise OR of these values describes a group of
 *		     related conversions
 */
enum conversion_group {
    /*
     * Type of conversion.
     */
    CONV_xtoy	= (1 << 3),
    CONV_xtoI	= (2 << 3),
    CONV_NONE	= (3 << 3),
    CONV_MASK	= (3 << 3),

    /*
     * Dimension of the input matrix.
     */
    DIM_ANY	= (0 << 1),
    DIM_VEC	= (1 << 1),
    DIM_2x2	= (2 << 1),
    DIM_NxN	= (3 << 1),
    DIM_MASK	= (3 << 1),

    /*
     * Does the conversion require Z0?
     */
    Z0_NO	= 0,
    Z0_YES	= 1,
    Z0_MASK	= 1
};

/*
 * Shift and mask values for packing and unpacking the conversion_code.
 */
#define GROUP_SHIFT	8
#define GROUP_MASK	(0xFF00)
#define INDEX_SHIFT	0
#define INDEX_MASK	(0x00FF)

/*
 * MAKE_CODE: make a conversion code from group and index
 *   @group: bitwise or of the values in conversion_group
 *   @index: index into group_* function pointer tables
 */
#define MAKE_CODE(group, index) \
	((group) << GROUP_SHIFT | (index) << INDEX_SHIFT)

/*
 * GET_GROUP: extract the conversion group information the conversion code
 *   @code: a member of the conversion_code enum
 */
#define GET_GROUP(code) \
	(((code) & GROUP_MASK) >> GROUP_SHIFT)

/*
 * GET_INDEX: extract the function pointer index from the conversion code
 *   @code: a member of the conversion_code enum
 */
#define GET_INDEX(code) \
	(((code) & INDEX_MASK) >> INDEX_SHIFT)

/*
 * conversion_code_t: codes describing each possible conversion
 */
typedef enum conversion_code {
    INVAL	= 0x0000,		/* invalid conversion */

    /*
     * Group 2x2_no_xtoy: 2-port to 2-port without Z0
     */
    T0StoT = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy,  0),
    T0TtoS = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy,  1),

    T0ZtoH = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy,  2),
    T0ZtoG = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy,  3),
    T0ZtoA = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy,  4),
    T0ZtoB = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy,  5),

    T0YtoH = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy,  6),
    T0YtoG = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy,  7),
    T0YtoA = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy,  8),
    T0YtoB = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy,  9),

    T0HtoZ = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy, 10),
    T0HtoY = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy, 11),
    T0HtoG = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy, 12),
    T0HtoA = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy, 13),
    T0HtoB = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy, 14),

    T0GtoZ = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy, 15),
    T0GtoY = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy, 16),
    T0GtoH = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy, 17),
    T0GtoA = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy, 18),
    T0GtoB = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy, 19),

    T0AtoZ = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy, 20),
    T0AtoY = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy, 21),
    T0AtoH = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy, 22),
    T0AtoG = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy, 23),
    T0AtoB = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy, 24),

    T0BtoZ = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy, 25),
    T0BtoY = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy, 26),
    T0BtoH = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy, 27),
    T0BtoG = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy, 28),
    T0BtoA = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_xtoy, 29),

    /*
     * Group 2x2_yes_xtoy: 2-port to 2-port with Z0
     */
    T1StoH = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoy,  0),
    T1StoG = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoy,  1),
    T1StoA = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoy,  2),
    T1StoB = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoy,  3),

    T1TtoZ = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoy,  4),
    T1TtoY = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoy,  5),
    T1TtoH = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoy,  6),
    T1TtoG = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoy,  7),
    T1TtoA = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoy,  8),
    T1TtoB = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoy,  9),

    T1ZtoT = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoy, 10),

    T1YtoT = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoy, 11),

    T1HtoS = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoy, 12),
    T1HtoT = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoy, 13),

    T1GtoS = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoy, 14),
    T1GtoT = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoy, 15),

    T1AtoS = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoy, 16),
    T1AtoT = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoy, 17),

    T1BtoS = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoy, 18),
    T1BtoT = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoy, 19),

    /*
     * Group 2x2_yes_x2I: 2-port to Zin vector
     */
    T1TtoI = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoI,  0),
    T1HtoI = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoI,  1),
    T1GtoI = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoI,  2),
    T1AtoI = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoI,  3),
    T1BtoI = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_xtoI,  4),

    /*
     * Group NxN_no_xtoy: N-port to N-port, no Z0
     */
    N0ZtoY = MAKE_CODE(DIM_NxN | Z0_NO  | CONV_xtoy,  0),
    N0YtoZ = MAKE_CODE(DIM_NxN | Z0_NO  | CONV_xtoy,  1),

    /*
     * Group NxN_yes_xtoy: N-port to N-port, with Z0
     */
    N1StoZ = MAKE_CODE(DIM_NxN | Z0_YES | CONV_xtoy,  0),
    N1StoY = MAKE_CODE(DIM_NxN | Z0_YES | CONV_xtoy,  1),
    N1ZtoS = MAKE_CODE(DIM_NxN | Z0_YES | CONV_xtoy,  2),
    N1YtoS = MAKE_CODE(DIM_NxN | Z0_YES | CONV_xtoy,  3),

    /*
     * Group NxN_yes_x2I: N-port to Zin vector
     */
    N1StoI = MAKE_CODE(DIM_NxN | Z0_YES | CONV_xtoI,  0),
    N1ZtoI = MAKE_CODE(DIM_NxN | Z0_YES | CONV_xtoI,  1),
    N1YtoI = MAKE_CODE(DIM_NxN | Z0_YES | CONV_xtoI,  2),

    /*
     * Input and output types the same -- no conversion.
     */
    ASAME = MAKE_CODE(DIM_ANY | Z0_NO  | CONV_NONE,  0),
    TSAME = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_NONE,  0),
    NSAME = MAKE_CODE(DIM_NxN | Z0_NO  | CONV_NONE,  0),
    VSAME = MAKE_CODE(DIM_VEC | Z0_NO  | CONV_NONE,  0),

} conversion_code_t;

/*
 * conversion_table: map a vnadata_parameter_type_t pair to conversion_code_t
 *   Row index is a member of vnadata_parameter_type_t describing the type
 *   of the input matrix.  Column index is the new type.
 *
 *   Name format: [ANTV][01]conversion
 *	A: any dimensions
 *	N: NxN
 *	T: 2x2
 *	V: row vector
 *
 *	0: no    z0 argument
 *	1: needs z0 argument
 *
 *	xtoy	convert x to y, with I == Zin
 *	same	no conversion; just copy
 *
 */
static const conversion_code_t conversion_table[VPT_NTYPES][VPT_NTYPES] = {
     /*  -       S       T       U       Z       Y       H       G       A       B       I  */
/*-*/{ ASAME,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL },
/*S*/{ INVAL,  ASAME, T0StoT,  INVAL, N1StoZ, N1StoY, T1StoH, T1StoG, T1StoA, T1StoB, N1StoI },
/*T*/{ INVAL, T0TtoS,  TSAME,  INVAL, T1TtoZ, T1TtoY, T1TtoH, T1TtoG, T1TtoA, T1TtoB, T1TtoI },
/*U*/{ INVAL,  INVAL,  INVAL,  TSAME,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL },
/*Z*/{ INVAL, N1ZtoS, T1ZtoT,  INVAL,  NSAME, N0ZtoY, T0ZtoH, T0ZtoG, T0ZtoA, T0ZtoB, N1ZtoI },
/*Y*/{ INVAL, N1YtoS, T1YtoT,  INVAL, N0YtoZ,  NSAME, T0YtoH, T0YtoG, T0YtoA, T0YtoB, N1YtoI },
/*H*/{ INVAL, T1HtoS, T1HtoT,  INVAL, T0HtoZ, T0HtoY,  TSAME, T0HtoG, T0HtoA, T0HtoB, T1HtoI },
/*G*/{ INVAL, T1GtoS, T1GtoT,  INVAL, T0GtoZ, T0GtoY, T0GtoH,  TSAME, T0GtoA, T0GtoB, T1GtoI },
/*A*/{ INVAL, T1AtoS, T1AtoT,  INVAL, T0AtoZ, T0AtoY, T0AtoH, T0AtoG,  TSAME, T0AtoB, T1AtoI },
/*B*/{ INVAL, T1BtoS, T1BtoT,  INVAL, T0BtoZ, T0BtoY, T0BtoH, T0BtoG, T0BtoA,  TSAME, T1BtoI },
/*I*/{ INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  VSAME }
};

/*
 * group_2x2_no_xtoy: 2-port to 2-port conversion functions without z0
 */
static void (*group_2x2_no_xtoy[])(const double complex (*in)[2],
	double complex (*out)[2]) = {
    [GET_INDEX(T0StoT)] = vnaconv_stot,

    [GET_INDEX(T0TtoS)] = vnaconv_ttos,

    [GET_INDEX(T0ZtoH)] = vnaconv_ztoh,
    [GET_INDEX(T0ZtoG)] = vnaconv_ztog,
    [GET_INDEX(T0ZtoA)] = vnaconv_ztoa,
    [GET_INDEX(T0ZtoB)] = vnaconv_ztob,

    [GET_INDEX(T0YtoH)] = vnaconv_ytoh,
    [GET_INDEX(T0YtoG)] = vnaconv_ytog,
    [GET_INDEX(T0YtoA)] = vnaconv_ytoa,
    [GET_INDEX(T0YtoB)] = vnaconv_ytob,

    [GET_INDEX(T0HtoZ)] = vnaconv_htoz,
    [GET_INDEX(T0HtoY)] = vnaconv_htoy,
    [GET_INDEX(T0HtoG)] = vnaconv_htog,
    [GET_INDEX(T0HtoA)] = vnaconv_htoa,
    [GET_INDEX(T0HtoB)] = vnaconv_htob,

    [GET_INDEX(T0GtoZ)] = vnaconv_gtoz,
    [GET_INDEX(T0GtoY)] = vnaconv_gtoy,
    [GET_INDEX(T0GtoH)] = vnaconv_gtoh,
    [GET_INDEX(T0GtoA)] = vnaconv_gtoa,
    [GET_INDEX(T0GtoB)] = vnaconv_gtob,

    [GET_INDEX(T0AtoZ)] = vnaconv_atoz,
    [GET_INDEX(T0AtoY)] = vnaconv_atoy,
    [GET_INDEX(T0AtoH)] = vnaconv_atoh,
    [GET_INDEX(T0AtoG)] = vnaconv_atog,
    [GET_INDEX(T0AtoB)] = vnaconv_atob,

    [GET_INDEX(T0BtoZ)] = vnaconv_btoz,
    [GET_INDEX(T0BtoY)] = vnaconv_btoy,
    [GET_INDEX(T0BtoH)] = vnaconv_btoh,
    [GET_INDEX(T0BtoG)] = vnaconv_btog,
    [GET_INDEX(T0BtoA)] = vnaconv_btoa
};

/*
 * group_2x2_yes_xtoy: 2-port to 2-port conversion functions with z0
 */
static void (*group_2x2_yes_xtoy[])(const double complex (*in)[2],
	double complex (*out)[2], const double complex *z0) = {
    [GET_INDEX(T1StoH)] = vnaconv_stoh,
    [GET_INDEX(T1StoG)] = vnaconv_stog,
    [GET_INDEX(T1StoA)] = vnaconv_stoa,
    [GET_INDEX(T1StoB)] = vnaconv_stob,

    [GET_INDEX(T1TtoZ)] = vnaconv_ttoz,
    [GET_INDEX(T1TtoY)] = vnaconv_ttoy,
    [GET_INDEX(T1TtoH)] = vnaconv_ttoh,
    [GET_INDEX(T1TtoG)] = vnaconv_ttog,
    [GET_INDEX(T1TtoA)] = vnaconv_ttoa,
    [GET_INDEX(T1TtoB)] = vnaconv_ttob,

    [GET_INDEX(T1ZtoT)] = vnaconv_ztot,

    [GET_INDEX(T1YtoT)] = vnaconv_ytot,

    [GET_INDEX(T1HtoS)] = vnaconv_htos,
    [GET_INDEX(T1HtoT)] = vnaconv_htot,

    [GET_INDEX(T1GtoS)] = vnaconv_gtos,
    [GET_INDEX(T1GtoT)] = vnaconv_gtot,

    [GET_INDEX(T1AtoT)] = vnaconv_atot,
    [GET_INDEX(T1AtoS)] = vnaconv_atos,

    [GET_INDEX(T1BtoS)] = vnaconv_btos,
    [GET_INDEX(T1BtoT)] = vnaconv_btot
};

/*
 * group_2x2_yes_xtoI: 2-port to Zin vector conversion functions with z0
 */
static void (*group_2x2_yes_xtoI[])(const double complex (*in)[2],
	double complex *out, const double complex *z0) = {
    [GET_INDEX(T1TtoI)] = vnaconv_ttozi,
    [GET_INDEX(T1HtoI)] = vnaconv_htozi,
    [GET_INDEX(T1GtoI)] = vnaconv_gtozi,
    [GET_INDEX(T1AtoI)] = vnaconv_atozi,
    [GET_INDEX(T1BtoI)] = vnaconv_btozi
};

/*
 * group_NxN_no_xtoy: N-port to N-port conversion functions without z0
 */
static void (*group_NxN_no_xtoy[])(const double complex *in,
	double complex *out, int n) = {
    [GET_INDEX(N0ZtoY)] = vnaconv_ztoyn,
    [GET_INDEX(N0YtoZ)] = vnaconv_ytozn
};

/*
 * group_NxN_yes_xtoy: N-port to N-port conversion functions with z0
 */
static void (*group_NxN_yes_xtoy[])(const double complex *in,
	double complex *out, const double complex *z0, int n) = {
    [GET_INDEX(N1StoZ)] = vnaconv_stozn,
    [GET_INDEX(N1StoY)] = vnaconv_stoyn,
    [GET_INDEX(N1ZtoS)] = vnaconv_ztosn,
    [GET_INDEX(N1YtoS)] = vnaconv_ytosn
};

/*
 * group_NxN_yes_xtoI: N-port to Zin vector conversion functions with z0
 */
static void (*group_NxN_yes_xtoI[])(const double complex *in,
	double complex *out, const double complex *z0, int n) = {
    [GET_INDEX(N1StoI)] = vnaconv_stozin,
    [GET_INDEX(N1ZtoI)] = vnaconv_ztozin,
    [GET_INDEX(N1YtoI)] = vnaconv_ytozin
};

/*
 * get_fz0_vector: get the z0 vector (fast)
 *   @vdip:   internal parameter matrix
 *   @findex: frequency index
 */
static double complex *get_fz0_vector(vnadata_internal_t *vdip, int findex)
{
    if (vdip->vdi_flags & VF_PER_F_Z0) {
	return vdip->vdi_z0_vector_vector[findex];
    }
    return vdip->vdi_z0_vector;
}

/*
 * vnadata_convert: convert between parameter types
 *   @vdp_in:  input parameter matrix
 *   @vdp_out: output parameter matrix
 *   @newtype: new type (can be the same as old)
 *
 * Note: vdp_out and vdp_in may be the same.
 */
int vnadata_convert(const vnadata_t *vdp_in, vnadata_t *vdp_out,
	vnadata_parameter_type_t newtype)
{
    vnadata_internal_t *vdip_in;
    conversion_code_t conversion;
    int group, index;

    /*
     * Sanity check the arguments.
     */
    if (vdp_in == NULL) {
	errno = EINVAL;
	return -1;
    }
    vdip_in = VDP_TO_VDIP(vdp_in);
    if (vdip_in->vdi_magic != VDI_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    if (vdp_out == NULL) {
	_vnadata_error(vdip_in, VNAERR_USAGE,
		"vnadata_convert: vdp_out cannot be NULL");
	return -1;
    }
    if (newtype < 0 || newtype >= VPT_NTYPES) {
	_vnadata_error(vdip_in, VNAERR_USAGE,
		"vnadata_convert: invalid new type: %d", (int)newtype);
	return -1;
    }

    /*
     * Look-up the conversion.  Fail if it's invalid.
     */
    conversion = conversion_table[vdp_in->vd_type][newtype];
    if (conversion == INVAL) {
	_vnadata_error(vdip_in, VNAERR_USAGE,
		"vnadata_convert: cannot convert from %s to %s",
		vnadata_get_type_name(vdp_in->vd_type),
		vnadata_get_type_name(newtype));
	return -1;
    }
    group = GET_GROUP(conversion);
    index = GET_INDEX(conversion);

    /*
     * Check the input dimensions.
     */
    switch (group & DIM_MASK) {
    case DIM_ANY:
	break;

    case DIM_VEC:
	if (vdp_in->vd_rows != 1 && vdp_in->vd_columns != 1) {
	    _vnadata_error(vdip_in, VNAERR_USAGE, "vnadata_convert: "
		    "invalid input dimensions: %d x %d: must be vector",
		    vdp_in->vd_rows, vdp_in->vd_columns);
	    return -1;
	}
	break;

    case DIM_2x2:
	if (vdp_in->vd_rows != 2 || vdp_in->vd_columns != 2) {
	    _vnadata_error(vdip_in, VNAERR_USAGE, "vnadata_convert: "
		    "invalid input dimensions: %d x %d: must be 2x2",
		    vdp_in->vd_rows, vdp_in->vd_columns);
	    return -1;
	}
	break;

    case DIM_NxN:
	if (vdp_in->vd_rows != vdp_in->vd_columns) {
	    _vnadata_error(vdip_in, VNAERR_USAGE, "vnadata_convert: "
		    "invalid input dimensions: %d x %d: must be square",
		    vdp_in->vd_rows, vdp_in->vd_columns);
	    return -1;
	}
	break;

    default:
	abort();
	/*NOTREACHED*/
    }

    /*
     * Set up the destination structure if it's not the same as the source.
     * Initially set the type to VPT_UNDEF.
     */
    if (vdp_out != vdp_in) {
	int new_rows    = vdp_in->vd_rows;
	int new_columns = vdp_in->vd_columns;
	int rc;

	/*
	 * If converting from matrix to vector, make it a row vector
	 * of length number of ports.
	 */
	if ((group & CONV_MASK) == CONV_xtoI) {
	    if (new_rows < new_columns) {
		new_columns = new_rows;
	    }
	    new_rows = 1;
	}

	/*
	 * Set up the output matrix.  Transfer everything over except
	 * for error_fn and error_arg.
	 */
	rc = vnadata_init(vdp_out, VPT_UNDEF, new_rows, new_columns,
		vdp_in->vd_frequencies);
	if (rc == -1) {
	    return -1;
	}
	vnadata_set_frequency_vector(vdp_out, vdp_in->vd_frequency_vector);
	if (!(vdip_in->vdi_flags & VF_PER_F_Z0)) {
	    if (vnadata_set_z0_vector(vdp_out, vdip_in->vdi_z0_vector) == -1) {
		return -1;
	    }
	} else {
	    int frequencies = vdp_in->vd_frequencies;

	    for (int findex = 0; findex < frequencies; ++findex) {
		if (vnadata_set_fz0_vector(vdp_out, findex,
			    vdip_in->vdi_z0_vector_vector[findex]) == -1) {
		    return -1;
		}
	    }
	}
	if (vnadata_set_filetype(vdp_out, vdip_in->vdi_filetype) == -1) {
	    return -1;
	}
	if (vnadata_set_format(vdp_out, vdip_in->vdi_format_string) == -1) {
	    return -1;
	}
	if (vnadata_set_fprecision(vdp_out, vdip_in->vdi_fprecision) == -1) {
	    return -1;
	}
	if (vnadata_set_dprecision(vdp_out, vdip_in->vdi_dprecision) == -1) {
	    return -1;
	}
    }

    /*
     * Handle the case of old and new types already the same.  If the
     * output matrix is a different structure, copy the values and type over.
     */
    if (newtype == vdp_in->vd_type) {
	if (vdp_out != vdp_in) {
	    int cells = vdp_in->vd_rows * vdp_in->vd_columns;

	    for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
		(void)memcpy(vdp_out->vd_data[findex], vdp_in->vd_data[findex],
		    cells * sizeof(double complex));
	    }
	    vdp_out->vd_type = newtype;
	}
	return 0;
    }

    /*
     * Do the conversion.
     */
    switch (group) {
    case DIM_2x2 | Z0_NO  | CONV_xtoy:
	{
	    void (*fn)(const double complex (*in)[2], double complex (*out)[2]);

	    fn = group_2x2_no_xtoy[index];
	    for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
		(*fn)((const double complex (*)[2])vdp_in->vd_data[findex],
			(double complex (*)[2])vdp_out->vd_data[findex]);
	    }
	}
	break;

    case DIM_2x2 | Z0_YES | CONV_xtoy:
	{
	    void (*fn)(const double complex (*in)[2], double complex (*out)[2],
		    const double complex *z0);

	    fn = group_2x2_yes_xtoy[index];
	    for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
		(*fn)((const double complex (*)[2])vdp_in->vd_data[findex],
			(double complex (*)[2])vdp_out->vd_data[findex],
			get_fz0_vector(vdip_in, findex));
	    }
	}
	break;

    case DIM_2x2 | Z0_YES | CONV_xtoI:
	{
	    void (*fn)(const double complex (*in)[2], double complex *out,
		    const double complex *z0);

	    fn = group_2x2_yes_xtoI[index];
	    for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
		(*fn)((const double complex (*)[2])vdp_in->vd_data[findex],
			     vdp_out->vd_data[findex],
			     get_fz0_vector(vdip_in, findex));
	    }
	}
	break;

    case DIM_NxN | Z0_NO  | CONV_xtoy:
	{
	    void (*fn)(const double complex *in, double complex *out, int n);

	    fn = group_NxN_no_xtoy[index];
	    for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
		(*fn)(vdp_in->vd_data[findex], vdp_out->vd_data[findex],
		    vdp_in->vd_rows);
	    }
	}
	break;

    case DIM_NxN | Z0_YES | CONV_xtoy:
	{
	    void (*fn)(const double complex *in, double complex *out,
		    const double complex *z0, int n);

	    fn = group_NxN_yes_xtoy[index];
	    for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
		(*fn)(vdp_in->vd_data[findex], vdp_out->vd_data[findex],
		    get_fz0_vector(vdip_in, findex), vdp_in->vd_rows);
	    }
	}
	break;

    case DIM_NxN | Z0_YES | CONV_xtoI:
	{
	    void (*fn)(const double complex *in, double complex *out,
		    const double complex *z0, int n);

	    fn = group_NxN_yes_xtoI[index];
	    for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
		(*fn)(vdp_in->vd_data[findex], vdp_out->vd_data[findex],
		    get_fz0_vector(vdip_in, findex), vdp_in->vd_rows);
	    }
	}
	break;

    default:
	abort();
    }
    vdp_out->vd_type = newtype;

    /*
     * For the case of an in-place conversion of a square matrix to Zin,
     * change the dimensions to a row vector.
     */
    if (vdp_in == vdp_out && (group & CONV_MASK) == CONV_xtoI) {
	if (vdp_out->vd_rows < vdp_out->vd_columns) {
	    vdp_out->vd_columns = vdp_out->vd_rows;
	}
	vdp_out->vd_rows = 1;
    }
    return 0;
}
