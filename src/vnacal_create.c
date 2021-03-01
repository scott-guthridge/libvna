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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"


/*
 * _vnacal_cal_diagonal_error_terms: calculate diagonal error terms
 *   @etsp: internal calibration data set
 *   @vcsp: measured calibration data
 *   @column: matrix column index
 */
static void _vnacal_calc_diagonal_error_terms(vnacal_etermset_t *etsp,
	vnacal_calset_t *vcsp, int column)
{
    vnacal_cdata_t *vcdp;
    const vnacal_error_terms_t *etp;

    vcdp = VNACAL_CALIBRATION_DATA(vcsp, column, column);
    etp = VNACAL_ERROR_TERMS(etsp, column, column);
    for (int findex = 0; findex < vcsp->vcs_frequencies; ++findex) {
	double complex a0 = _vnacal_calset_get_reference(vcsp, 0, findex);
	double complex a1 = _vnacal_calset_get_reference(vcsp, 1, findex);
	double complex a2 = _vnacal_calset_get_reference(vcsp, 2, findex);
	double complex m0, m1, m2;
	double complex d;

	/* Compute the diagonal error terms. */
	m0 = _vnacal_calset_get_value(vcdp, VNACAL_Sii_REF0, findex);
	m1 = _vnacal_calset_get_value(vcdp, VNACAL_Sii_REF1, findex);
	m2 = _vnacal_calset_get_value(vcdp, VNACAL_Sii_REF2, findex);
	d = a0 * a1 * (m0 - m1) + a1 * a2 * (m1 - m2) + a2 * a0 * (m2 - m0);
	etp->et_e00[findex] = (a0 * a1 * (m0 - m1) * m2 +
				  a1 * a2 * (m1 - m2) * m0 +
				  a2 * a0 * (m2 - m0) * m1) / d;
	etp->et_e10e01[findex] = (a0 - a1) * (a1 - a2) * (a2 - a0) *
				    (m0 - m1) * (m1 - m2) * (m2 - m0) / (d * d);
	etp->et_e11[findex] = (a0 * (m2 - m1) +
				  a1 * (m0 - m2) +
				  a2 * (m1 - m0)) / d;
    }
}

/*
 * _vnacal_calc_off_diagonal_error_terms: calculate remaining error terms
 *   @etsp: internal calibration data set
 *   @vcsp: measured calibration data
 *   @row: matrix row index
 *   @column: matrix column index
 */
static void _vnacal_calc_off_diagonal_error_terms(vnacal_etermset_t *etsp,
	const vnacal_calset_t *vcsp, int row, int column)
{
    vnacal_cdata_t *vcdp;
    const vnacal_error_terms_t *etp;
    const vnacal_error_terms_t *vcep_diagonal;
    int findex;

    /*
     * Get the needed matrix cells.  The diagonal element, which exists
     * only if column < rows, is used to compute e10e32 and e22 below.
     */
    vcdp = VNACAL_CALIBRATION_DATA(vcsp, row, column);
    etp = VNACAL_ERROR_TERMS(etsp, row, column);
    vcep_diagonal = (column < vcsp->vcs_rows) ?
	VNACAL_ERROR_TERMS(etsp, column, column) : NULL;

    /*
     * For each frequency, calculate the error terms.
     */
    for (findex = 0; findex < vcsp->vcs_frequencies; ++findex) {
	double complex e00    = 0.0;
	double complex e10e01 = 1.0;
	double complex e11    = 0.0;
	double complex sii_through, sij_through, sij_leakage;

	if (vcep_diagonal != NULL) {
	    e00    = vcep_diagonal->et_e00[findex];
	    e10e01 = vcep_diagonal->et_e10e01[findex];
	    e11    = vcep_diagonal->et_e11[findex];
	}
	sii_through = _vnacal_calset_get_value(vcdp, VNACAL_Sjj_THROUGH,
		findex);
	sij_through = _vnacal_calset_get_value(vcdp, VNACAL_Sij_THROUGH,
		findex);
	sij_leakage = _vnacal_calset_get_value(vcdp, VNACAL_Sij_LEAKAGE,
		findex);
	etp->et_e30[findex] = sij_leakage;
	etp->et_e10e32[findex] = e10e01 * (sij_through - sij_leakage) /
	                        (e10e01 + (sii_through - e00) * e11);
	etp->et_e22[findex] = (sii_through - e00) /
	                      (e10e01 + (sii_through - e00) * e11);
    }
}

/*
 * vnacal_create: construct a calibration structure from measured data
 *   @sets: number of calibration sets
 *   @vcspp: vector of pointers to vnacal_calset structures (sets long)
 *   @error_fn: optional error reporting function (NULL if not used)
 *   @error_arg: user data passed through to the error function (or NULL)
 */
vnacal_t *vnacal_create(int sets, vnacal_calset_t **vcspp,
	vnaerr_error_fn_t *error_fn, void *error_arg)
{
    vnacal_t *vcp;

    /*
     * Allocate the vnacal_t.
     */
    if ((vcp = (vnacal_t *)malloc(sizeof(vnacal_t))) == NULL) {
	if (error_fn != NULL) {
	    int saved_errno = errno;
	    char message[80];

	    (void)snprintf(message, sizeof(message), "vnacal_create: %s",
		    strerror(errno));
	    message[sizeof(message)-1] = '\000';
	    errno = saved_errno;
	    (*error_fn)(VNAERR_SYSTEM, message, error_arg);
	    errno = saved_errno;
	}
	return NULL;
    }
    (void)memset((void *)vcp, 0, sizeof(vnacal_t));
    vcp->vc_sets = sets;
    vcp->vc_fprecision = VNACAL_DEFAULT_DATA_PRECISION;
    vcp->vc_dprecision = VNACAL_DEFAULT_FREQUENCY_PRECISION;
    vcp->vc_error_fn = error_fn;
    vcp->vc_error_arg = error_arg;

    /*
     * Allocate vc_set_vector.
     */
    if ((vcp->vc_set_vector = (vnacal_etermset_t **)calloc(sets,
		    sizeof(vnacal_etermset_t *))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"calloc %s", strerror(errno));
	vnacal_free(vcp);
	return NULL;
    }
    for (int set = 0; set < sets; ++set) {
	vnacal_calset_t *vcsp = vcspp[set];
	int min_dimension = MIN(vcsp->vcs_rows, vcsp->vcs_columns);
	double cfmin, cfmax;
	vnacal_etermset_t *etsp;

	/*
	 * Validate dimensions
	 */
	if (vcsp->vcs_rows < 1) {
	    _vnacal_error(vcp, VNAERR_USAGE,
		    "vnacal_create: set %d: rows must be >= 1", set);
	    vnacal_free(vcp);
	    return NULL;
	}
	if (vcsp->vcs_columns < 1) {
	    _vnacal_error(vcp, VNAERR_USAGE,
		    "vnacal_create: set %d: columns must be >= 1", set);
	    vnacal_free(vcp);
	    return NULL;
	}
	if (vcsp->vcs_frequencies < 1) {
	    _vnacal_error(vcp, VNAERR_USAGE,
		    "vnacal_create: set %d: frequencies must be >= 1", set);
	    vnacal_free(vcp);
	    return NULL;
	}

	/*
	 * Make sure the frequency vector was given.
	 */
	if (!vcsp->vcs_frequencies_valid) {
	    _vnacal_error(vcp, VNAERR_USAGE, "vnacal_create: set %d: "
		    "no calibration frequency vector given", set);
	    return NULL;
	}

	/*
	 * Validate the reference frequency ranges.
	 */
	cfmin = vcsp->vcs_frequency_vector[0];
	cfmax = vcsp->vcs_frequency_vector[vcsp->vcs_frequencies - 1];
	for (int i = 0; i < 3; ++i) {
	    vnacal_calset_reference_t *vcdsrp = &vcsp->vcs_references[i];
	    int n;
	    double lower, upper;

	    if (!vcdsrp->vcdsr_is_vector) {
		continue;
	    }
	    n = vcdsrp->u.v.vcdsr_frequencies;
	    lower = (1.0 - VNACAL_F_EXTRAPOLATION) *
		vcdsrp->u.v.vcdsr_frequency_vector[0];
	    upper = (1.0 + VNACAL_F_EXTRAPOLATION) *
		vcdsrp->u.v.vcdsr_frequency_vector[n - 1];
	    if (cfmin < lower || cfmax > upper) {
		_vnacal_error(vcp, VNAERR_USAGE, "vnacal_create: set %d: "
		   "error: frequency range %.3e..%.3e is outside of "
		   "reference %d range %.3e..%.3e",
		   set, cfmin, cfmax, i,
		   vcdsrp->u.v.vcdsr_frequency_vector[0],
		   vcdsrp->u.v.vcdsr_frequency_vector[n - 1]);
		return NULL;
	    }
	}

	/*
	 * Allocate the vnacal_etermset_t and copy in the frequency
	 * vector.
	 */
	if ((vcp->vc_set_vector[set] = _vnacal_etermset_alloc(vcp,
			vcsp->vcs_setname,
			vcsp->vcs_rows, vcsp->vcs_columns,
			vcsp->vcs_frequencies)) == NULL) {
	    vnacal_free(vcp);
	    return NULL;
	}
	etsp = vcp->vc_set_vector[set];
	(void)memcpy((void *)etsp->ets_frequency_vector,
		(void *)vcsp->vcs_frequency_vector,
		vcsp->vcs_frequencies * sizeof(double));

	/*
	 * Copy the system impedance.
	 */
	etsp->ets_z0 = vcsp->vcs_z0;

	/*
	 * Calculate the diagonal error terms.
	 */
	for (int j = 0; j < min_dimension; ++j) {
	    _vnacal_calc_diagonal_error_terms(etsp, vcsp, j);
	}

	/*
	 * Calculate the off-diagonal error terms.
	 */
	for (int i = 0; i < vcsp->vcs_rows; ++i) {
	    for (int j = 0; j < vcsp->vcs_columns; ++j) {
		if (j != i) {
		    _vnacal_calc_off_diagonal_error_terms(etsp, vcsp, i, j);
		}
	    }
	}
    }
    return vcp;
}
