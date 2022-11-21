/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2022 D Scott Guthridge <scott_guthridge@rompromity.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "archdep.h"

#include <assert.h>
#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_new_internal.h"

typedef enum {
    TRL_T,		/* standard is through */
    TRL_R,		/* standard is reflect */
    TRL_L,		/* standard is line */
    TRL_NONE		/* standard is none of the above */
} trl_which_t;

/*
 * classify_standard: given a standard, determine if it's T, R, L or none
 *   @vnmp: measured standard
 *   @unknown_index: returned unparameter index (or -1 if T)
 */
static trl_which_t classify_standard(const vnacal_new_measurement_t *vnmp,
	int *unknown_index)
{
    const vnacal_new_t *vnp = vnmp->vnm_vnp;
    const vnacal_t *vcp = vnp->vn_vcp;
    const vnacal_parameter_t *vnprp_one;
    vnacal_new_parameter_t **s = vnmp->vnm_s_matrix;

    *unknown_index = -1;
    vnprp_one = _vnacal_get_parameter(vcp, VNACAL_ONE);
    assert(vnprp_one != NULL);
    if (s[1]->vnpr_parameter == vnprp_one) {
	if (s[2]->vnpr_parameter == vnprp_one &&
		s[0] == vnp->vn_zero && s[3] == vnp->vn_zero) {
	    return TRL_T;
	}
	return TRL_NONE;
    }
    if (s[1] == vnp->vn_zero) {
	if (s[0]->vnpr_parameter->vpmr_type == VNACAL_UNKNOWN && s[3] == s[0] &&
		s[2] == vnp->vn_zero) {
	    *unknown_index = s[0]->vnpr_unknown_index;
	    return TRL_R;
	}
	return TRL_NONE;
    }
    if (s[0] == vnp->vn_zero && s[3] == vnp->vn_zero &&
	    s[1]->vnpr_parameter->vpmr_type == VNACAL_UNKNOWN && s[2] == s[1]) {
	*unknown_index = s[1]->vnpr_unknown_index;
	return TRL_L;
    }
    return TRL_NONE;
}

/*
 * _vnacal_new_solve_is_trl: test if the system can be solved using simple TRL
 *   @vnp: pointer to vnacal_new_t structure
 */
bool _vnacal_new_solve_is_trl(const vnacal_new_t *vnp,
	vnacal_new_trl_indices_t *vntip)
{
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    int unknown_index;
    int standard_index;

    /*
     * Init the return structure.
     */
    vntip->vnti_t_standard = -1;
    vntip->vnti_r_standard = -1;
    vntip->vnti_l_standard = -1;
    vntip->vnti_r_unknown  = -1;
    vntip->vnti_l_unknown  = -1;

    /*
     * If the calibration dimensions are not 2x2 or the type is not
     * T8, U8, TE10 or UE10, then it's not simple TRL.
     */
    if (VL_M_ROWS(vlp) != 2 || VL_M_COLUMNS(vlp) != 2 ||
	    (VL_TYPE(vlp) != VNACAL_T8 && VL_TYPE(vlp) != VNACAL_TE10 &&
	     VL_TYPE(vlp) != VNACAL_U8 && VL_TYPE(vlp) != VNACAL_UE10)) {
	return false;
    }

    /*
     * If there are other than 3 standards, other than 2 unknown
     * parameters, or any correlated parameters given, then it's not
     * a match.
     */
    if (vnp->vn_measurement_count != 3) {
	return false;
    }
    if (vnp->vn_unknown_parameters != 2) {
	return false;
    }
    if (vnp->vn_correlated_parameters != 0) {
	return false;
    }

    /*
     * If a measurement error model was given, then use the general
     * autocalibration method instead.
     */
    if (vnp->vn_m_error_vector != NULL) {
	return false;
    }

    /*
     * Classify each standard and collect the unknown indices for R and L.
     */
    standard_index = 0;
    for (vnacal_new_measurement_t *vnmp = vnp->vn_measurement_list;
	    vnmp != NULL; vnmp = vnmp->vnm_next, ++standard_index) {
	switch (classify_standard(vnmp, &unknown_index)) {
	case TRL_T:
	    if (vntip->vnti_t_standard != -1) {	/* check for duplicate T */
		return false;
	    }
	    vntip->vnti_t_standard = standard_index;
	    continue;

	case TRL_R:
	    if (vntip->vnti_r_standard != -1) {	/* check for duplicate R */
		return false;
	    }
	    vntip->vnti_r_unknown = unknown_index;
	    vntip->vnti_r_standard = standard_index;
	    continue;

	case TRL_L:
	    if (vntip->vnti_l_standard != -1) {	/* check for duplicate L */
		return false;
	    }
	    vntip->vnti_l_unknown = unknown_index;
	    vntip->vnti_l_standard = standard_index;
	    continue;

	default:
	    return false;
	}
    }
    return true;
}

#define TRL_EQUATIONS	10
#define TRL_UNKNOWNS	 7

/*
 * _vnacal_new_solve_trl: solve when all s-parameters are known
 *   @vnssp: solve state structure
 *   @x_vector: vector of unknowns filled in by this function
 *   @x_length: length of x_vector
 */
int _vnacal_new_solve_trl(vnacal_new_solve_state_t *vnssp,
	const vnacal_new_trl_indices_t *vntip,
	double complex *x_vector, int x_length)
{
    vnacal_new_t *vnp = vnssp->vnss_vnp;
    vnacal_t *vcp = vnp->vn_vcp;
    const int findex = vnssp->vnss_findex;
    const double complex *mt =
	vnssp->vnss_msv_matrices[vntip->vnti_t_standard].vnmm_m_matrix;
    const double complex *mr =
	vnssp->vnss_msv_matrices[vntip->vnti_r_standard].vnmm_m_matrix;
    const double complex *ml =
	vnssp->vnss_msv_matrices[vntip->vnti_l_standard].vnmm_m_matrix;
    double complex mt11 = mt[0];
    double complex mt12 = mt[1];
    double complex mt21 = mt[2];
    double complex mt22 = mt[3];
    double complex mr11 = mr[0];
    double complex mr22 = mr[3];
    double complex ml11 = ml[0];
    double complex ml12 = ml[1];
    double complex ml21 = ml[2];
    double complex ml22 = ml[3];
    double complex r, l;
    double complex a_matrix[TRL_EQUATIONS][TRL_UNKNOWNS];
    double complex b_vector[TRL_EQUATIONS];
    int eq_count = 0;

    assert(vnp->vn_systems == 1);
    assert(x_length == TRL_UNKNOWNS);

    /*
     * Solve for the unknown line parameter, l.  There are two
     * solutions: choose the one closest to the initial guess.
     */
    {
	double complex a, b, c, u, v;
	double complex guess;
	double d1, d2;

	a = ml12 * mt21;
	b = (ml11 - mt11) * (ml22 - mt22) - ml12 * ml21 - mt12 * mt21;
	c = mt12 * ml21;
	u = -b / (2.0 * a);
	v = csqrt(b*b - 4.0 * a * c) / (2.0 * a);
	guess = vnssp->vnss_p_vector[vntip->vnti_l_unknown][findex];
	d1 = cabs(u + v - guess);
	d2 = cabs(u - v - guess);
	if (d1 <= d2) {
	    l = u + v;
	} else {
	    l = u - v;
	}
    }

    /*
     * Calculate the unknown reflect parameter, r.  There are two
     * solutions: choose the one closest to the intial guess.
     */
    {
	double complex n, d;
	double complex guess;
	double d1, d2;

	n = (-ml12 * mt21 + (mt12 * mt21 - (ml11 - mt11) * (mr22 - mt22)) * l) *
		((mr11 - mt11) * (ml22 - mt22) * l + mt12 * (ml21 - mt21 * l));
	d = (ml21 * (mt22 - mr22) + mt21 * (mr22 - ml22) * l) *
		(-mt11 * ml12 + ml11 * mt12 * l + mr11 * (ml12 - mt12 * l));
	if (d == 0.0) {
	    _vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
		    "solution of unknown reflect parameter is singular");
	    return -1;
	}
	r = csqrt(n / d);
	guess = vnssp->vnss_p_vector[vntip->vnti_r_unknown][findex];
	d1 = cabs(+r - guess);
	d2 = cabs(-r - guess);
	if (d1 > d2) {
	    r = -r;
	}
    }
    vnssp->vnss_p_vector[vntip->vnti_r_unknown][findex] = r;
    vnssp->vnss_p_vector[vntip->vnti_l_unknown][findex] = l;
    vs_update_s_matrices(vnssp);

    /*
     * Build the coefficient matrix (a) and right-hand side vector (b).
     */
    for (int i = 0; i < TRL_EQUATIONS; ++i) {
	for (int j = 0; j < TRL_UNKNOWNS; ++j) {
	    a_matrix[i][j] = 0.0;
	}
	b_vector[i] = 0.0;
    }
    vs_start_system(vnssp, 0);
    while (vs_next_equation(vnssp)) {
	while (vs_next_term(vnssp)) {
	    double complex value = vs_get_negative(vnssp) ? -1.0 : 1.0;
	    const int xindex = vs_get_xindex(vnssp);

	    if (vs_have_m(vnssp)) {
		value *= vs_get_m(vnssp);
	    }
	    if (vs_have_s(vnssp)) {
		value *= vs_get_s(vnssp);
	    }
	    if (vs_have_v(vnssp)) {
		value *= vs_get_v(vnssp);
	    }
	    assert(eq_count < TRL_EQUATIONS);
	    if (xindex == -1) {
		b_vector[eq_count] += value;
	    } else {
		a_matrix[eq_count][xindex] += value;
	    }
	}
	++eq_count;
    }
    assert(eq_count == vnp->vn_equations);
    {
	int rank;

	rank = _vnacommon_qrsolve(x_vector, &a_matrix[0][0],
		b_vector, eq_count, TRL_UNKNOWNS, 1);
	if (rank < TRL_UNKNOWNS) {
	    _vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
		    "singular linear system");
	    return -1;
	}
    }
    return 1;
}
