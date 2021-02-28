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
    CONV_x2y	= (1 << 3),
    CONV_x2I	= (2 << 3),
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
 * GET_INDEX: extra the function pointer index from the conversion code
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
     * Group 2x2_no_x2y: 2-port to 2-port without Z0
     */
    T0S2T = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y,  0),
    T0Z2H = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y,  1),
    T0Z2G = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y,  2),
    T0Z2A = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y,  3),
    T0Z2B = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y,  4),
    T0Y2H = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y,  5),
    T0Y2G = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y,  6),
    T0Y2A = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y,  7),
    T0Y2B = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y,  8),
    T0T2S = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y,  9),
    T0H2Z = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y, 10),
    T0H2Y = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y, 11),
    T0H2G = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y, 12),
    T0H2A = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y, 13),
    T0H2B = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y, 14),
    T0G2Z = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y, 15),
    T0G2Y = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y, 16),
    T0G2H = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y, 17),
    T0G2A = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y, 18),
    T0G2B = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y, 19),
    T0A2Z = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y, 20),
    T0A2Y = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y, 21),
    T0A2H = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y, 22),
    T0A2G = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y, 23),
    T0A2B = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y, 24),
    T0B2Z = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y, 25),
    T0B2Y = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y, 26),
    T0B2H = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y, 27),
    T0B2G = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y, 28),
    T0B2A = MAKE_CODE(DIM_2x2 | Z0_NO  | CONV_x2y, 29),

    /*
     * Group 2x2_yes_x2y: 2-port to 2-port with Z0
     */
    T1S2H = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2y,  0),
    T1S2G = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2y,  1),
    T1S2A = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2y,  2),
    T1S2B = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2y,  3),
    T1Z2T = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2y,  4),
    T1Y2T = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2y,  5),
    T1T2Z = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2y,  6),
    T1T2Y = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2y,  7),
    T1T2H = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2y,  8),
    T1T2G = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2y,  9),
    T1T2A = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2y, 10),
    T1T2B = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2y, 11),
    T1H2S = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2y, 12),
    T1H2T = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2y, 13),
    T1G2S = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2y, 14),
    T1G2T = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2y, 15),
    T1A2T = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2y, 16),
    T1A2S = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2y, 17),
    T1B2S = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2y, 18),
    T1B2T = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2y, 19),

    /*
     * Group 2x2_yes_x2I: 2-port to Zin vector
     */
    T1T2I = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2I,  0),
    T1H2I = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2I,  1),
    T1G2I = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2I,  2),
    T1A2I = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2I,  3),
    T1B2I = MAKE_CODE(DIM_2x2 | Z0_YES | CONV_x2I,  4),

    /*
     * Group NxN_no_x2y: N-port to N-port, no Z0
     */
    N0Z2Y = MAKE_CODE(DIM_NxN | Z0_NO  | CONV_x2y,  0),
    N0Y2Z = MAKE_CODE(DIM_NxN | Z0_NO  | CONV_x2y,  1),

    /*
     * Group NxN_yes_x2y: N-port to N-port, with Z0
     */
    N1S2Z = MAKE_CODE(DIM_NxN | Z0_YES | CONV_x2y,  0),
    N1S2Y = MAKE_CODE(DIM_NxN | Z0_YES | CONV_x2y,  1),
    N1Z2S = MAKE_CODE(DIM_NxN | Z0_YES | CONV_x2y,  2),
    N1Y2S = MAKE_CODE(DIM_NxN | Z0_YES | CONV_x2y,  3),

    /*
     * Group NxN_yes_x2I: N-port to Zin vector
     */
    N1Z2I = MAKE_CODE(DIM_NxN | Z0_YES | CONV_x2I,  0),
    N1Y2I = MAKE_CODE(DIM_NxN | Z0_YES | CONV_x2I,  1),

    /*
     * Group NxN_yes_x2I: MxN to Zin vector
     */
    M1S2I = MAKE_CODE(DIM_ANY | Z0_YES | CONV_x2I,  0),

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
 */
static const conversion_code_t conversion_table[VPT_NTYPES][VPT_NTYPES] = {
     /*  -      S      Z      Y      T      H      G      A      B      I  */
/*-*/{ ASAME, INVAL, INVAL, INVAL, INVAL, INVAL, INVAL, INVAL, INVAL, INVAL },
/*S*/{ INVAL, ASAME, N1S2Z, N1S2Y, T0S2T, T1S2H, T1S2G, T1S2A, T1S2B, M1S2I },
/*Z*/{ INVAL, N1Z2S, NSAME, N0Z2Y, T1Z2T, T0Z2H, T0Z2G, T0Z2A, T0Z2B, N1Z2I },
/*Y*/{ INVAL, N1Y2S, N1Y2I, NSAME, T1Y2T, T0Y2H, T0Y2G, T0Y2A, T0Y2B, N1Y2I },
/*T*/{ INVAL, T0T2S, T1T2Z, T1T2Y, TSAME, T1T2H, T1T2G, T1T2A, T1T2B, T1T2I },
/*H*/{ INVAL, T1H2S, T0H2Z, T0H2Y, T1H2T, TSAME, T0H2G, T0H2A, T0H2B, T1H2I },
/*G*/{ INVAL, T1G2S, T0G2Z, T0G2Y, T1G2T, T0G2H, TSAME, T0G2A, T0G2B, T1G2I },
/*A*/{ INVAL, T1A2S, T0A2Z, T0A2Y, T1A2T, T0A2H, T0A2G, TSAME, T0A2B, T1A2I },
/*B*/{ INVAL, T1B2S, T0B2Z, T0B2Y, T1B2T, T0B2H, T0B2G, T0B2A, TSAME, T1B2I },
/*I*/{ INVAL, INVAL, INVAL, INVAL, INVAL, INVAL, INVAL, INVAL, INVAL, VSAME }
};

/*
 * group_2x2_no_x2y: 2-port to 2-port conversion functions without z0
 */
static void (*group_2x2_no_x2y[])(const vnaconv_array2_t *in,
	double complex (*out)[2]) = {
    [GET_INDEX(T0S2T)] = vnaconv_s2t,
    [GET_INDEX(T0Z2H)] = vnaconv_z2h,
    [GET_INDEX(T0Z2G)] = vnaconv_z2g,
    [GET_INDEX(T0Z2A)] = vnaconv_z2a,
    [GET_INDEX(T0Z2B)] = vnaconv_z2b,
    [GET_INDEX(T0Y2H)] = vnaconv_y2h,
    [GET_INDEX(T0Y2G)] = vnaconv_y2g,
    [GET_INDEX(T0Y2A)] = vnaconv_y2a,
    [GET_INDEX(T0Y2B)] = vnaconv_y2b,
    [GET_INDEX(T0T2S)] = vnaconv_t2s,
    [GET_INDEX(T0H2Z)] = vnaconv_h2z,
    [GET_INDEX(T0H2Y)] = vnaconv_h2y,
    [GET_INDEX(T0H2G)] = vnaconv_h2g,
    [GET_INDEX(T0H2A)] = vnaconv_h2a,
    [GET_INDEX(T0H2B)] = vnaconv_h2b,
    [GET_INDEX(T0G2Z)] = vnaconv_g2z,
    [GET_INDEX(T0G2Y)] = vnaconv_g2y,
    [GET_INDEX(T0G2H)] = vnaconv_g2h,
    [GET_INDEX(T0G2A)] = vnaconv_g2a,
    [GET_INDEX(T0G2B)] = vnaconv_g2b,
    [GET_INDEX(T0A2Z)] = vnaconv_a2z,
    [GET_INDEX(T0A2Y)] = vnaconv_a2y,
    [GET_INDEX(T0A2H)] = vnaconv_a2h,
    [GET_INDEX(T0A2G)] = vnaconv_a2g,
    [GET_INDEX(T0A2B)] = vnaconv_a2b,
    [GET_INDEX(T0B2Z)] = vnaconv_b2z,
    [GET_INDEX(T0B2Y)] = vnaconv_b2y,
    [GET_INDEX(T0B2H)] = vnaconv_b2h,
    [GET_INDEX(T0B2G)] = vnaconv_b2g,
    [GET_INDEX(T0B2A)] = vnaconv_b2a
};

/*
 * group_2x2_yes_x2y: 2-port to 2-port conversion functions with z0
 */
static void (*group_2x2_yes_x2y[])(const vnaconv_array2_t *in,
	double complex (*out)[2], const double complex *z0) = {
    [GET_INDEX(T1S2H)] = vnaconv_s2h,
    [GET_INDEX(T1S2G)] = vnaconv_s2g,
    [GET_INDEX(T1S2A)] = vnaconv_s2a,
    [GET_INDEX(T1S2B)] = vnaconv_s2b,
    [GET_INDEX(T1Z2T)] = vnaconv_z2t,
    [GET_INDEX(T1Y2T)] = vnaconv_y2t,
    [GET_INDEX(T1T2Z)] = vnaconv_t2z,
    [GET_INDEX(T1T2Y)] = vnaconv_t2y,
    [GET_INDEX(T1T2H)] = vnaconv_t2h,
    [GET_INDEX(T1T2G)] = vnaconv_t2g,
    [GET_INDEX(T1T2A)] = vnaconv_t2a,
    [GET_INDEX(T1T2B)] = vnaconv_t2b,
    [GET_INDEX(T1H2S)] = vnaconv_h2s,
    [GET_INDEX(T1H2T)] = vnaconv_h2t,
    [GET_INDEX(T1G2S)] = vnaconv_g2s,
    [GET_INDEX(T1G2T)] = vnaconv_g2t,
    [GET_INDEX(T1A2T)] = vnaconv_a2t,
    [GET_INDEX(T1A2S)] = vnaconv_a2s,
    [GET_INDEX(T1B2S)] = vnaconv_b2s,
    [GET_INDEX(T1B2T)] = vnaconv_b2t
};

/*
 * group_2x2_yes_x2I: 2-port to Zin vector conversion functions with z0
 */
static void (*group_2x2_yes_x2I[])(const vnaconv_array2_t *in,
	double complex *out, const double complex *z0) = {
    [GET_INDEX(T1T2I)] = vnaconv_t2zi,
    [GET_INDEX(T1H2I)] = vnaconv_h2zi,
    [GET_INDEX(T1G2I)] = vnaconv_g2zi,
    [GET_INDEX(T1A2I)] = vnaconv_a2zi,
    [GET_INDEX(T1B2I)] = vnaconv_b2zi
};

/*
 * group_NxN_no_x2y: N-port to N-port conversion functions without z0
 */
static void (*group_NxN_no_x2y[])(const double complex *in,
	double complex *out, int n) = {
    [GET_INDEX(N0Z2Y)] = vnaconv_z2yn,
    [GET_INDEX(N0Y2Z)] = vnaconv_y2zn
};

/*
 * group_NxN_yes_x2y: N-port to N-port conversion functions with z0
 */
static void (*group_NxN_yes_x2y[])(const double complex *in,
	double complex *out, const double complex *z0, int n) = {
    [GET_INDEX(N1S2Z)] = vnaconv_s2zn,
    [GET_INDEX(N1S2Y)] = vnaconv_s2yn,
    [GET_INDEX(N1Z2S)] = vnaconv_z2sn,
    [GET_INDEX(N1Y2S)] = vnaconv_y2sn
};

/*
 * group_NxN_yes_x2I: N-port to Zin vector conversion functions with z0
 */
static void (*group_NxN_yes_x2I[])(const double complex *in,
	double complex *out, const double complex *z0, int n) = {
    [GET_INDEX(N1Z2I)] = vnaconv_z2zin,
    [GET_INDEX(N1Y2I)] = vnaconv_y2zin
};

/*
 * group_MxN_yes_x2I: MxN to Zin vector conversion functions with z0
 */
static void (*group_MxN_yes_x2I[])(const double complex *in,
	double complex *out, const double complex *z0, int m, int n) = {
    [GET_INDEX(M1S2I)] = vnaconv_s2zimn,
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
    if (vdp_out == NULL || vdp_in == NULL ||
	    vdp_in->vd_type < 0 || vdp_in->vd_type >= VPT_NTYPES ||
	    newtype < 0 || newtype >= VPT_NTYPES) {
	errno = EINVAL;
	return -1;
    }
    vdip_in = VDP_TO_VDIP(vdp_in);
    if (vdip_in->vdi_magic != VDI_MAGIC) {
	errno = EINVAL;
	return -1;
    }

    /*
     * Look-up the conversion.  Fail if it's invalid.
     */
    conversion = conversion_table[vdp_in->vd_type][newtype];
    if (conversion == INVAL) {
	errno = EINVAL;
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
	    errno = EINVAL;
	    return -1;
	}
	break;

    case DIM_2x2:
	if (vdp_in->vd_rows != 2 || vdp_in->vd_columns != 2) {
	    errno = EINVAL;
	    return -1;
	}
	break;

    case DIM_NxN:
	if (vdp_in->vd_rows != vdp_in->vd_columns) {
	    errno = EINVAL;
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
	const double complex *z0_vector;

	/*
	 * If converting from matrix to vector, make it a row vector
	 * of length number of ports.
	 */
	if ((group & CONV_MASK) == CONV_x2I) {
	    if (new_rows < new_columns) {
		new_columns = new_rows;
	    }
	    new_rows = 1;
	}

	/*
	 * Set up the output matrix.
	 */
	rc = vnadata_init(vdp_out, vdp_in->vd_frequencies,
		new_rows, new_columns, VPT_UNDEF);
	if (rc == -1) {
	    return -1;
	}
	vnadata_set_frequency_vector(vdp_out,
		vnadata_get_frequency_vector(vdp_in));
	if ((z0_vector = vnadata_get_z0_vector(vdp_in)) != NULL) {
	    if (vnadata_set_z0_vector(vdp_out, z0_vector) == -1) {
		return -1;
	    }
	} else {
	    int frequencies = vnadata_get_frequencies(vdp_in);

	    for (int findex = 0; findex < frequencies; ++findex) {
		if (vnadata_set_fz0_vector(vdp_out, findex,
			    vnadata_get_fz0_vector(vdp_in, findex)) == -1) {
		    return -1;
		}
	    }
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
    case DIM_2x2 | Z0_NO  | CONV_x2y:
	{
	    void (*fn)(const vnaconv_array2_t *in, double complex (*out)[2]);

	    fn = group_2x2_no_x2y[index];
	    for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
		(*fn)((const vnaconv_array2_t *)vdp_in->vd_data[findex],
			(double complex (*)[2])vdp_out->vd_data[findex]);
	    }
	}
	break;

    case DIM_2x2 | Z0_YES | CONV_x2y:
	{
	    void (*fn)(const vnaconv_array2_t *in, double complex (*out)[2],
		    const double complex *z0);

	    fn = group_2x2_yes_x2y[index];
	    for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
		(*fn)((const vnaconv_array2_t *)vdp_in->vd_data[findex],
			(double complex (*)[2])vdp_out->vd_data[findex],
			get_fz0_vector(vdip_in, findex));
	    }
	}
	break;

    case DIM_2x2 | Z0_YES | CONV_x2I:
	{
	    void (*fn)(const vnaconv_array2_t *in, double complex *out,
		    const double complex *z0);

	    fn = group_2x2_yes_x2I[index];
	    for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
		(*fn)((const vnaconv_array2_t *)vdp_in->vd_data[findex],
			     vdp_out->vd_data[findex],
			     get_fz0_vector(vdip_in, findex));
	    }
	}
	break;

    case DIM_NxN | Z0_NO  | CONV_x2y:
	{
	    void (*fn)(const double complex *in, double complex *out, int n);

	    fn = group_NxN_no_x2y[index];
	    for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
		(*fn)(vdp_in->vd_data[findex], vdp_out->vd_data[findex],
		    vdp_in->vd_rows);
	    }
	}
	break;

    case DIM_NxN | Z0_YES | CONV_x2y:
	{
	    void (*fn)(const double complex *in, double complex *out,
		    const double complex *z0, int n);

	    fn = group_NxN_yes_x2y[index];
	    for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
		(*fn)(vdp_in->vd_data[findex], vdp_out->vd_data[findex],
		    get_fz0_vector(vdip_in, findex), vdp_in->vd_rows);
	    }
	}
	break;

    case DIM_NxN | Z0_YES | CONV_x2I:
	{
	    void (*fn)(const double complex *in, double complex *out,
		    const double complex *z0, int n);

	    fn = group_NxN_yes_x2I[index];
	    for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
		(*fn)(vdp_in->vd_data[findex], vdp_out->vd_data[findex],
		    get_fz0_vector(vdip_in, findex), vdp_in->vd_rows);
	    }
	}
	break;

    case DIM_ANY | Z0_YES | CONV_x2I:
	{
	    void (*fn)(const double complex *in, double complex *out,
		    const double complex *z0, int rows, int columns);

	    fn = group_MxN_yes_x2I[index];
	    for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
		(*fn)(vdp_in->vd_data[findex], vdp_out->vd_data[findex],
		    get_fz0_vector(vdip_in, findex), vdp_in->vd_rows,
		    vdp_in->vd_columns);
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
    if (vdp_in == vdp_out && (group & CONV_MASK) == CONV_x2I) {
	if (vdp_out->vd_rows < vdp_out->vd_columns) {
	    vdp_out->vd_columns = vdp_out->vd_rows;
	}
	vdp_out->vd_rows = 1;
    }
    return 0;
}
