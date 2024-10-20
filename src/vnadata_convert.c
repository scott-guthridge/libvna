/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2023 D Scott Guthridge <scott_guthridge@rompromity.net>
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

#define VNADATA_NO_BOUNDS_CHECK

#include "archdep.h"

#include <assert.h>
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
    CONV_xtoy	= (1 << 4),
    CONV_xtoI	= (2 << 4),	/* to input impedance at each port */
    CONV_NONE	= (3 << 4),
    CONV_MASK	= (3 << 4),

    /*
     * Dimension of the input matrix.
     */
    DIM_ANY	= (0 << 2),
    DIM_VEC	= (1 << 2),
    DIM_2x2	= (2 << 2),
    DIM_NxN	= (3 << 2),
    DIM_MASK	= (3 << 2),

    /*
     * How many z0's does the conversion need?  0, 1 or 2
     */
    Z0_NONE	= 0,
    Z0_ONE	= 1,
    Z0_TWO	= 2,
    Z0_MASK	= 3
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
     * Group 2x2_z0_xtoy: 2-port to 2-port without Z0
     */
    T0StoT = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy,  0),
    T0StoU = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy,  1),

    T0TtoS = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy,  2),
    T0TtoU = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy,  3),

    T0UtoS = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy,  4),
    T0UtoT = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy,  5),

    T0ZtoH = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy,  6),
    T0ZtoG = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy,  7),
    T0ZtoA = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy,  8),
    T0ZtoB = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy,  9),

    T0YtoH = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 10),
    T0YtoG = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 11),
    T0YtoA = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 12),
    T0YtoB = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 13),

    T0HtoZ = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 14),
    T0HtoY = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 15),
    T0HtoG = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 16),
    T0HtoA = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 17),
    T0HtoB = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 18),

    T0GtoZ = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 19),
    T0GtoY = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 20),
    T0GtoH = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 21),
    T0GtoA = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 22),
    T0GtoB = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 23),

    T0AtoZ = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 24),
    T0AtoY = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 25),
    T0AtoH = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 26),
    T0AtoG = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 27),
    T0AtoB = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 28),

    T0BtoZ = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 29),
    T0BtoY = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 30),
    T0BtoH = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 31),
    T0BtoG = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 32),
    T0BtoA = MAKE_CODE(DIM_2x2 | Z0_NONE  | CONV_xtoy, 33),

    /*
     * Group 2x2_z1_xtoy: 2-port to 2-port with Z0
     */
    T1StoH = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy,  0),
    T1StoG = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy,  1),
    T1StoA = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy,  2),
    T1StoB = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy,  3),

    T1TtoZ = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy,  4),
    T1TtoY = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy,  5),
    T1TtoH = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy,  6),
    T1TtoG = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy,  7),
    T1TtoA = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy,  8),
    T1TtoB = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy,  9),

    T1UtoZ = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 10),
    T1UtoY = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 11),
    T1UtoH = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 12),
    T1UtoG = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 13),
    T1UtoA = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 14),
    T1UtoB = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 15),

    T1ZtoT = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 16),
    T1ZtoU = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 17),

    T1YtoT = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 18),
    T1YtoU = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 19),

    T1HtoS = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 20),
    T1HtoT = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 21),
    T1HtoU = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 22),

    T1GtoS = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 23),
    T1GtoT = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 24),
    T1GtoU = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 25),

    T1AtoS = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 26),
    T1AtoT = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 27),
    T1AtoU = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 28),

    T1BtoS = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 29),
    T1BtoT = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 30),
    T1BtoU = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoy, 31),

    /*
     * Group 2x2_z1_x2I: 2-port to Zin vector
     */
    T1TtoI = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoI,  0),
    T1UtoI = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoI,  1),
    T1HtoI = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoI,  2),
    T1GtoI = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoI,  3),
    T1AtoI = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoI,  4),
    T1BtoI = MAKE_CODE(DIM_2x2 | Z0_ONE | CONV_xtoI,  5),

    /*
     * Group NxN_z0_xtoy: N-port to N-port, no Z0
     */
    N0ZtoY = MAKE_CODE(DIM_NxN | Z0_NONE  | CONV_xtoy,  0),
    N0YtoZ = MAKE_CODE(DIM_NxN | Z0_NONE  | CONV_xtoy,  1),

    /*
     * Group NxN_z1_xtoy: N-port to N-port, with Z0
     */
    N1StoZ = MAKE_CODE(DIM_NxN | Z0_ONE | CONV_xtoy,  0),
    N1StoY = MAKE_CODE(DIM_NxN | Z0_ONE | CONV_xtoy,  1),
    N1ZtoS = MAKE_CODE(DIM_NxN | Z0_ONE | CONV_xtoy,  2),
    N1YtoS = MAKE_CODE(DIM_NxN | Z0_ONE | CONV_xtoy,  3),

    /*
     * Group NxN_z1_x2I: N-port to Zin vector
     */
    N1StoI = MAKE_CODE(DIM_NxN | Z0_ONE | CONV_xtoI,  0),
    N1ZtoI = MAKE_CODE(DIM_NxN | Z0_ONE | CONV_xtoI,  1),
    N1YtoI = MAKE_CODE(DIM_NxN | Z0_ONE | CONV_xtoI,  2),

    /*
     * Input and output types the same -- no conversion.
     */
    ASAME = MAKE_CODE(DIM_ANY | Z0_NONE | CONV_NONE,  0),
    TSAME = MAKE_CODE(DIM_2x2 | Z0_NONE | CONV_NONE,  0),
    NSAME = MAKE_CODE(DIM_NxN | Z0_NONE | CONV_NONE,  0),
    VSAME = MAKE_CODE(DIM_VEC | Z0_NONE | CONV_NONE,  0),

    /*
     * Group 2x2_z2_xtoy: 2-port to 2-port, re-normalizing
     */
    T2StoT = MAKE_CODE(DIM_2x2 | Z0_TWO | CONV_xtoy, 0),
    T2StoU = MAKE_CODE(DIM_2x2 | Z0_TWO | CONV_xtoy, 1),
    T2TtoS = MAKE_CODE(DIM_2x2 | Z0_TWO | CONV_xtoy, 2),
    T2TtoT = MAKE_CODE(DIM_2x2 | Z0_TWO | CONV_xtoy, 3),
    T2TtoU = MAKE_CODE(DIM_2x2 | Z0_TWO | CONV_xtoy, 4),
    T2UtoS = MAKE_CODE(DIM_2x2 | Z0_TWO | CONV_xtoy, 5),
    T2UtoT = MAKE_CODE(DIM_2x2 | Z0_TWO | CONV_xtoy, 6),
    T2UtoU = MAKE_CODE(DIM_2x2 | Z0_TWO | CONV_xtoy, 7),

    /*
     * group NxN_z2_xtoy: N-port to N-port, re-normalizing.
     */
    N2StoS = MAKE_CODE(DIM_NxN | Z0_TWO | CONV_xtoy, 0),

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
 *	0: no  z0 arguments
 *	1: one z0 argument
 *
 *	xtoy	convert x to y, with I == Zin
 *	same	no conversion; just copy
 */
static const conversion_code_t conversion_table[VPT_NTYPES][VPT_NTYPES] = {
     /*  -       S       T       U       Z       Y       H       G       A       B       I  */
/*-*/{ ASAME,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL },
/*S*/{ INVAL,  NSAME, T0StoT, T0StoU, N1StoZ, N1StoY, T1StoH, T1StoG, T1StoA, T1StoB, N1StoI },
/*T*/{ INVAL, T0TtoS,  TSAME, T0TtoU, T1TtoZ, T1TtoY, T1TtoH, T1TtoG, T1TtoA, T1TtoB, T1TtoI },
/*U*/{ INVAL, T0UtoS, T0UtoT,  TSAME, T1UtoZ, T1UtoY, T1UtoH, T1UtoG, T1UtoA, T1UtoB, T1UtoI },
/*Z*/{ INVAL, N1ZtoS, T1ZtoT, T1ZtoU,  NSAME, N0ZtoY, T0ZtoH, T0ZtoG, T0ZtoA, T0ZtoB, N1ZtoI },
/*Y*/{ INVAL, N1YtoS, T1YtoT, T1YtoU, N0YtoZ,  NSAME, T0YtoH, T0YtoG, T0YtoA, T0YtoB, N1YtoI },
/*H*/{ INVAL, T1HtoS, T1HtoT, T1HtoU, T0HtoZ, T0HtoY,  TSAME, T0HtoG, T0HtoA, T0HtoB, T1HtoI },
/*G*/{ INVAL, T1GtoS, T1GtoT, T1GtoU, T0GtoZ, T0GtoY, T0GtoH,  TSAME, T0GtoA, T0GtoB, T1GtoI },
/*A*/{ INVAL, T1AtoS, T1AtoT, T1AtoU, T0AtoZ, T0AtoY, T0AtoH, T0AtoG,  TSAME, T0AtoB, T1AtoI },
/*B*/{ INVAL, T1BtoS, T1BtoT, T1BtoU, T0BtoZ, T0BtoY, T0BtoH, T0BtoG, T0BtoA,  TSAME, T1BtoI },
/*I*/{ INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  INVAL,  VSAME }
};

/*
 * renormalizing_table: like above, but for re-normalizing conversions
 *
 *   Name format: [ANTV][01]conversion
 *	A: any dimensions
 *	N: NxN
 *	T: 2x2
 *
 *	2: two z0 arguments
 *
 *	xtoy	convert x to y, with I == Zin
 *	same	no conversion; just copy
 */
static const conversion_code_t renormalizing_table[4][4] = {
	   /* -      S       T       U */
    /*-*/{  ASAME, INVAL,  INVAL,  INVAL  },
    /*S*/{  INVAL, N2StoS, T2StoT, T2StoU },
    /*T*/{  INVAL, T2TtoS, T2TtoT, T2TtoU },
    /*U*/{  INVAL, T2UtoS, T2UtoT, T2UtoU }
};

/*
 * group_2x2_z0_xtoy: 2-port to 2-port conversion functions without z0
 */
static void (*group_2x2_z0_xtoy[])(const double complex (*in)[2],
	double complex (*out)[2]) = {
    [GET_INDEX(T0StoT)] = vnaconv_stot,
    [GET_INDEX(T0StoU)] = vnaconv_stou,

    [GET_INDEX(T0TtoS)] = vnaconv_ttos,
    [GET_INDEX(T0TtoU)] = vnaconv_ttou,

    [GET_INDEX(T0UtoS)] = vnaconv_utos,
    [GET_INDEX(T0UtoT)] = vnaconv_utot,

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
 * group_2x2_z1_xtoy: 2-port to 2-port conversion functions with z0
 */
static void (*group_2x2_z1_xtoy[])(const double complex (*in)[2],
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

    [GET_INDEX(T1UtoZ)] = vnaconv_utoz,
    [GET_INDEX(T1UtoY)] = vnaconv_utoy,
    [GET_INDEX(T1UtoH)] = vnaconv_utoh,
    [GET_INDEX(T1UtoG)] = vnaconv_utog,
    [GET_INDEX(T1UtoA)] = vnaconv_utoa,
    [GET_INDEX(T1UtoB)] = vnaconv_utob,

    [GET_INDEX(T1ZtoT)] = vnaconv_ztot,
    [GET_INDEX(T1ZtoU)] = vnaconv_ztou,

    [GET_INDEX(T1YtoT)] = vnaconv_ytot,
    [GET_INDEX(T1YtoU)] = vnaconv_ytou,

    [GET_INDEX(T1HtoS)] = vnaconv_htos,
    [GET_INDEX(T1HtoT)] = vnaconv_htot,
    [GET_INDEX(T1HtoU)] = vnaconv_htou,

    [GET_INDEX(T1GtoS)] = vnaconv_gtos,
    [GET_INDEX(T1GtoT)] = vnaconv_gtot,
    [GET_INDEX(T1GtoU)] = vnaconv_gtou,

    [GET_INDEX(T1AtoT)] = vnaconv_atot,
    [GET_INDEX(T1AtoS)] = vnaconv_atos,
    [GET_INDEX(T1AtoU)] = vnaconv_atou,

    [GET_INDEX(T1BtoS)] = vnaconv_btos,
    [GET_INDEX(T1BtoT)] = vnaconv_btot,
    [GET_INDEX(T1BtoU)] = vnaconv_btou
};

/*
 * group_2x2_z1_xtoI: 2-port to Zin vector conversion functions with z0
 */
static void (*group_2x2_z1_xtoI[])(const double complex (*in)[2],
	double complex *out, const double complex *z0) = {
    [GET_INDEX(T1TtoI)] = vnaconv_ttozi,
    [GET_INDEX(T1UtoI)] = vnaconv_utozi,
    [GET_INDEX(T1HtoI)] = vnaconv_htozi,
    [GET_INDEX(T1GtoI)] = vnaconv_gtozi,
    [GET_INDEX(T1AtoI)] = vnaconv_atozi,
    [GET_INDEX(T1BtoI)] = vnaconv_btozi
};

/*
 * group_2x2_z2_xtoy: 2-port to 2-port re-normalizing conversions
 */
static void (*group_2x2_z2_xtoy[])(const double complex (*in)[2],
	double complex (*out)[2], const double complex *z1,
	const double complex *z2) = {
    [GET_INDEX(T2StoT)] = vnaconv_stotr,
    [GET_INDEX(T2StoU)] = vnaconv_stour,
    [GET_INDEX(T2TtoS)] = vnaconv_ttosr,
    [GET_INDEX(T2TtoT)] = vnaconv_ttotr,
    [GET_INDEX(T2TtoU)] = vnaconv_ttour,
    [GET_INDEX(T2UtoS)] = vnaconv_utosr,
    [GET_INDEX(T2UtoT)] = vnaconv_utotr,
    [GET_INDEX(T2UtoU)] = vnaconv_utour,
};

/*
 * group_NxN_z0_xtoy: N-port to N-port conversion functions without z0
 */
static void (*group_NxN_z0_xtoy[])(const double complex *in,
	double complex *out, int n) = {
    [GET_INDEX(N0ZtoY)] = vnaconv_ztoyn,
    [GET_INDEX(N0YtoZ)] = vnaconv_ytozn
};

/*
 * group_NxN_z1_xtoy: N-port to N-port conversion functions with z0
 */
static void (*group_NxN_z1_xtoy[])(const double complex *in,
	double complex *out, const double complex *z0, int n) = {
    [GET_INDEX(N1StoZ)] = vnaconv_stozn,
    [GET_INDEX(N1StoY)] = vnaconv_stoyn,
    [GET_INDEX(N1ZtoS)] = vnaconv_ztosn,
    [GET_INDEX(N1YtoS)] = vnaconv_ytosn
};

/*
 * group_NxN_z1_xtoI: N-port to Zin vector conversion functions with z0
 */
static void (*group_NxN_z1_xtoI[])(const double complex *in,
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
 * vnadata_copy_what_t: what to copy
 */
typedef enum {
    COPY_TYPE		= 0x0,	/* copy type and dimensions only */
    COPY_FREQ		= 0x1,	/* copy frequency vector */
    COPY_DATA		= 0x2,	/* copy data array */
    COPY_Z0		= 0x4,	/* copy z0/fz0 */
    COPY_FILETYPE	= 0x8,	/* copy filetype, format, precision */
    COPY_ALL		= 0xF   /* copy everything */
} vnadata_copy_what_t;

/*
 * copy: copy one vnadata structure to another
 *   @vdp_in: source
 *   @vdp_out: destination
 *   @what: bitwise OR of vnadata_copy_what enums
 */
static int vnadata_copy(const vnadata_t *vdp_in, vnadata_t *vdp_out,
			vnadata_copy_what_t what)
{
    vnadata_internal_t *vdip_in = VDP_TO_VDIP(vdp_in);
    const int frequencies = vdp_in->vd_frequencies;
    const int rows        = vdp_in->vd_rows;
    const int columns     = vdp_in->vd_columns;
    int rv;

    if ((rv = vnadata_init(vdp_out, vdp_in->vd_type,
			   rows, columns, frequencies)) == -1) {
	return rv;
    }
    if (what & COPY_FREQ) {
	rv = vnadata_set_frequency_vector(vdp_out,
		vdp_in->vd_frequency_vector);
	assert(rv == 0);
    }
    if (what & COPY_DATA) {
	for (int findex = 0; findex < frequencies; ++findex) {
	    double complex *lfcp;

	    lfcp = vnadata_get_matrix(vdp_in, findex);
	    assert(lfcp != NULL);
	    rv = vnadata_set_matrix(vdp_out, findex, lfcp);
	    assert(rv == 0);
	}
    }
    if (what & COPY_Z0) {
	if (!(vdip_in->vdi_flags & VF_PER_F_Z0)) {
	    rv = vnadata_set_z0_vector(vdp_out, vdip_in->vdi_z0_vector);
	    assert(rv == 0);
	} else {
	    for (int findex = 0; findex < frequencies; ++findex) {
		rv = vnadata_set_fz0_vector(vdp_out, findex,
			vdip_in->vdi_z0_vector_vector[findex]);
		assert(rv == 0);
	    }
	}
    }
    if (what & COPY_FILETYPE) {
	rv = vnadata_set_filetype(vdp_out, vdip_in->vdi_filetype);
	assert(rv == 0);
	rv = vnadata_set_format(vdp_out, vdip_in->vdi_format_string);
	assert(rv == 0);
	rv = vnadata_set_fprecision(vdp_out, vdip_in->vdi_fprecision);
	assert(rv == 0);
	rv = vnadata_set_dprecision(vdp_out, vdip_in->vdi_dprecision);
	assert(rv == 0);
    }
    return 0;
}

/*
 * update_z0: update reference impedances in vdp
 *   @vdp: network parameter data to update
 *   @new_z0: vector/matrix of new z0 values
 *   @stride: number of values per frequency
 */
static void update_z0(vnadata_t *vdp, const double complex *new_z0, int stride)
{
    int rv;

    if (stride == 0) {
	rv = vnadata_set_z0_vector(vdp, new_z0);
	assert(rv == 0);

    } else {
	const int frequencies = vdp->vd_frequencies;
	int rv;

	for (int findex = 0; findex < frequencies; ++findex) {
	    rv = vnadata_set_fz0_vector(vdp, findex, new_z0);
	    assert(rv == 0);
	    new_z0 += stride;
	}
    }
}

/*
 * IS_SCATTERING: test if the given VNA parameter type is S, T or U.
 * 	@type: parameter type to test
 */
#define IS_SCATTERING_TYPE(type) \
	((type) == VPT_S || \
	 (type) == VPT_T || \
	 (type) == VPT_U)

/*
 * z0_update_strategy_t: describes type of z0 conversion
 */
typedef enum {
    Z0U_NONE,	/* no z0 conversion */
    Z0U_CtoC,	/* z, y, h, g, a, b -> z, y, y, g, a, b*/
    Z0U_CtoS,	/* z, y, h, g, a, b -> s, t, u */
    Z0U_StoC,	/* s, t, u -> z, y, h, g, a, b */
    Z0U_StoS	/* s, t, u -> s, t, u */
} z0_update_strategy_t;

/*
 * convert_common: common conversion function
 *   @vdp_in:  input network parameter data
 *   @vdp_out: output network parameter data
 *   @newtype: new type (can be the same as old)
 *   @new_z0:  new reference impedances (can be NULL)
 *   @length:  length of new_z0
 *   @function: name of external function
 */
static int convert_common(const vnadata_t *vdp_in, vnadata_t *vdp_out,
	vnadata_parameter_type_t newtype, const double complex *new_z0,
	int new_z0_length, const char *function)
{
    vnadata_internal_t *vdip_in;
    z0_update_strategy_t z0_update_strategy = Z0U_NONE;
    int z0_stride = 0;
    double complex *temp_z0_vector = NULL;
    conversion_code_t conversion = INVAL;
    int group, index;
    int rv = -1;

    /*
     * Sanity check the arguments.
     */
    if (vdp_in == NULL) {
	errno = EINVAL;
	goto out;
    }
    vdip_in = VDP_TO_VDIP(vdp_in);
    if (vdip_in->vdi_magic != VDI_MAGIC) {
	errno = EINVAL;
	goto out;
    }
    if (vdp_out == NULL) {
	_vnadata_error(vdip_in, VNAERR_USAGE,
		"%s: vdp_out cannot be NULL", function);
	goto out;
    }
    if (newtype < 0 || newtype >= VPT_NTYPES) {
	_vnadata_error(vdip_in, VNAERR_USAGE,
		"%s: invalid new type: %d", function, (int)newtype);
	goto out;
    }

    /*
     * If we're re-normalizing...
     */
    if (new_z0 != NULL) {
	const int frequencies = vdp_in->vd_frequencies;
	const int ports = vdp_in->vd_columns;

	/*
	 * Validate the length argument.  Special-case one z0 value given,
	 * and determine the stride through the z0 matrix by frequency.
	 */
	if (ports > 1 && new_z0_length == 1) {
	    if ((temp_z0_vector = calloc(ports,
			    sizeof(double complex))) == NULL) {
		_vnadata_error(vdip_in, VNAERR_SYSTEM, "calloc: %s",
			strerror(errno));
		goto out;
	    }
	    for (int i = 0; i < ports; ++i) {
		temp_z0_vector[i] = *new_z0;
	    }
	    new_z0 = temp_z0_vector;

	} else if (new_z0_length == ports) {
	    /*NULL*/;

	} else if (new_z0_length == frequencies * ports) {
	    z0_stride = ports;

	} else {
	    _vnadata_error(vdip_in, VNAERR_USAGE,
		    "%s: length must be 1, %d or %d",
		    function, ports, frequencies * ports);
	    goto out;
	}

	/*
	 * Determine how we're converting between types of scattering
	 * parameters (s, t, u) and circuit parameters (z, y, h, g, a, b).
	 */
	if (IS_SCATTERING_TYPE(vdp_in->vd_type)) {
	    if (IS_SCATTERING_TYPE(newtype)) {
		z0_update_strategy = Z0U_StoS;
	    } else {
		z0_update_strategy = Z0U_StoC;
	    }
	} else {
	    if (IS_SCATTERING_TYPE(newtype)) {
		z0_update_strategy = Z0U_CtoS;
	    } else {
		z0_update_strategy = Z0U_CtoC;
	    }
	}
    }

    /*
     * Look-up the conversion.  Fail if it's invalid.
     */
    if (z0_update_strategy != Z0U_StoS) {
	conversion = conversion_table[vdp_in->vd_type][newtype];
    } else {
	if (vdp_in->vd_type <= VPT_U && newtype <= VPT_U) {
	    conversion = renormalizing_table[vdp_in->vd_type][newtype];
	} else {
	    abort();
	}
    }
    if (conversion == INVAL) {
	_vnadata_error(vdip_in, VNAERR_USAGE,
		"%s: cannot convert from %s to %s", function,
		vnadata_get_type_name(vdp_in->vd_type),
		vnadata_get_type_name(newtype));
	goto out;
    }
    group = GET_GROUP(conversion);
    index = GET_INDEX(conversion);

    /*
     * Check the input dimensions.  These should be correct already,
     * but the vnadata_t header is visible in the .h file, so the caller
     * could have broken type and dimensions invariants.
     */
    switch (group & DIM_MASK) {
    case DIM_ANY:
	break;

    case DIM_VEC:
	if (vdp_in->vd_rows != 1) {
	    _vnadata_error(vdip_in, VNAERR_USAGE,
		    "%s: invalid input dimensions: %d x %d: must be row vector",
		    function, vdp_in->vd_rows, vdp_in->vd_columns);
	    goto out;
	}
	break;

    case DIM_2x2:
	if (vdp_in->vd_rows != 2 || vdp_in->vd_columns != 2) {
	    _vnadata_error(vdip_in, VNAERR_USAGE,
		    "%s: invalid input dimensions: %d x %d: must be 2x2",
		    function, vdp_in->vd_rows, vdp_in->vd_columns);
	    goto out;
	}
	break;

    case DIM_NxN:
	if (vdp_in->vd_rows != vdp_in->vd_columns) {
	    _vnadata_error(vdip_in, VNAERR_USAGE,
		    "%s: invalid input dimensions: %d x %d: must be square",
		    function, vdp_in->vd_rows, vdp_in->vd_columns);
	    goto out;
	}
	break;

    default:
	abort();
	/*NOTREACHED*/
    }

    /*
     * If we're re-normalizing and we converting from circuit parameters
     * (z, y, h, g, a, b) to any kind of scattering parameters (s, t, u),
     * update z0 before the conversion.  If we're not converting in-place,
     * then copy the input data structure to the output structure and do
     * an in-place conversion on the later.
     */
    if (z0_update_strategy == Z0U_CtoS) {
	if (vdp_out != vdp_in) {
	    if ((rv = vnadata_copy(vdp_in, vdp_out, COPY_ALL)) == -1) {
		goto out;
	    }
	    vdp_in = vdp_out;
	    vdip_in = VDP_TO_VDIP(vdp_in);
	}
	update_z0(vdp_out, new_z0, z0_stride);
    }

    /*
     * If the output structure isn't the same as the input, copy type
     * and shape, frequency vector and file format.  If no conversion
     * is needed, copy the data.  If not re-normalizing, copy z0.
     */
    if (vdp_out != vdp_in) {
	int what = COPY_FREQ | COPY_FILETYPE;

	if (z0_update_strategy == Z0U_NONE)
	    what |= COPY_Z0;
	if ((group & CONV_MASK) == CONV_NONE)
	    what |= COPY_DATA;
	if ((rv = vnadata_copy(vdp_in, vdp_out, what)) == -1) {
	    goto out;
	}
    }

    /*
     * Do the conversion.
     */
    switch (group) {
    case DIM_2x2 | Z0_NONE | CONV_xtoy:
	{
	    void (*fn)(const double complex (*in)[2], double complex (*out)[2]);

	    fn = group_2x2_z0_xtoy[index];
	    for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
		(*fn)((const double complex (*)[2])vdp_in->vd_data[findex],
			(double complex (*)[2])vdp_out->vd_data[findex]);
	    }
	}
	break;

    case DIM_2x2 | Z0_ONE | CONV_xtoy:
	{
	    void (*fn)(const double complex (*in)[2], double complex (*out)[2],
		    const double complex *z0);

	    fn = group_2x2_z1_xtoy[index];
	    for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
		(*fn)((const double complex (*)[2])vdp_in->vd_data[findex],
			(double complex (*)[2])vdp_out->vd_data[findex],
			get_fz0_vector(vdip_in, findex));
	    }
	}
	break;

    case DIM_2x2 | Z0_ONE | CONV_xtoI:
	{
	    void (*fn)(const double complex (*in)[2], double complex *out,
		    const double complex *z0);

	    fn = group_2x2_z1_xtoI[index];
	    for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
		(*fn)((const double complex (*)[2])vdp_in->vd_data[findex],
			     vdp_out->vd_data[findex],
			     get_fz0_vector(vdip_in, findex));
	    }
	}
	break;

    case DIM_2x2 | Z0_TWO | CONV_xtoy:
	{
	    void (*fn)(const double complex (*in)[2], double complex (*out)[2],
		    const double complex *z1, const double complex *z2);
	    fn = group_2x2_z2_xtoy[index];
	    for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
		const int offset = findex * z0_stride;

		(*fn)((const double complex (*)[2])vdp_in->vd_data[findex],
			(double complex (*)[2])vdp_out->vd_data[findex],
			get_fz0_vector(vdip_in, findex),
			&new_z0[offset]);
	    }
	}
	break;

    case DIM_NxN | Z0_NONE  | CONV_xtoy:
	{
	    void (*fn)(const double complex *in, double complex *out, int n);

	    fn = group_NxN_z0_xtoy[index];
	    for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
		(*fn)(vdp_in->vd_data[findex], vdp_out->vd_data[findex],
		    vdp_in->vd_rows);
	    }
	}
	break;

    case DIM_NxN | Z0_ONE | CONV_xtoy:
	{
	    void (*fn)(const double complex *in, double complex *out,
		    const double complex *z0, int n);

	    fn = group_NxN_z1_xtoy[index];
	    for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
		(*fn)(vdp_in->vd_data[findex], vdp_out->vd_data[findex],
		    get_fz0_vector(vdip_in, findex), vdp_in->vd_rows);
	    }
	}
	break;

    case DIM_NxN | Z0_ONE | CONV_xtoI:
	{
	    void (*fn)(const double complex *in, double complex *out,
		    const double complex *z0, int n);

	    fn = group_NxN_z1_xtoI[index];
	    for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
		(*fn)(vdp_in->vd_data[findex], vdp_out->vd_data[findex],
		    get_fz0_vector(vdip_in, findex), vdp_in->vd_rows);
	    }
	}
	break;

    case DIM_NxN | Z0_TWO | CONV_xtoy:
	assert(vdp_in->vd_type == VPT_S);
	assert(newtype == VPT_S);
	for (int findex = 0; findex < vdp_in->vd_frequencies; ++findex) {
	    const int offset = findex * z0_stride;

	    vnaconv_stosrn(vdp_in->vd_data[findex], vdp_out->vd_data[findex],
			   get_fz0_vector(vdip_in, findex),
			   &new_z0[offset], vdp_in->vd_rows);
	}
	break;

    default:
	/*
	 * Handle no conversion.
	 */
	if ((group & CONV_MASK) == CONV_NONE)
	    break;

	abort();
    }

    /*
     * If we're re-normalizing and we haven't already updated z0, do it now
     */
    if (z0_update_strategy != Z0U_NONE && z0_update_strategy != Z0U_CtoS) {
	update_z0(vdp_out, new_z0, z0_stride);
    }

    /*
     * If we're converting to input impedances, change the dimensions
     * to row vector.
     */
    if ((group & CONV_MASK) == CONV_xtoI) {
	vdp_out->vd_rows = 1;
    }

    /*
     * Set the new type.
     */
    vdp_out->vd_type = newtype;
    rv = 0;

out:
    free((void *)temp_z0_vector);
    return rv;
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
    return convert_common(vdp_in, vdp_out, newtype,
	    NULL, 0, "vnadata_convert");
}

/*
 * vnadata_rconvert: re-normalizing conversion
 *   @vdp_in:  input network parameter data
 *   @vdp_out: output network parameter data
 *   @newtype: new type (can be the same as old)
 *   @new_z0:  new reference impedances (can be NULL)
 *   @length:  length of new_z0
 *
 * Note: vdp_out and vdp_in may be the same.  The length field
 * gives the number of elements in new_z0.  It can be 1, #ports,
 * or #frequencies x #ports.
 */
int vnadata_rconvert(const vnadata_t *vdp_in, vnadata_t *vdp_out,
	vnadata_parameter_type_t newtype, const double complex *new_z0,
	int new_z0_length)
{
    return convert_common(vdp_in, vdp_out, newtype,
	    new_z0, new_z0_length, "vnadata_rconvert");
}
