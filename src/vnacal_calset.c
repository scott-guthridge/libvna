/*
 * Vector Network Analyzer Calibration Library
 * Copyright © 2020 D Scott Guthridge <scott_guthridge@rompromity.net>
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

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"


/*
 * vnacal_calset_error: report an error
 *   @vcsp: pointer to vnacal_calset_t
 *   @format: printf format string
 */
static void _vnacal_calset_error(const vnacal_calset_t *vcsp,
	const char *format, ...)
{
    va_list ap;
    char *msg = NULL;
    int rv;

    if (vcsp->vcs_error_fn != NULL) {
	va_start(ap, format);
	rv = vasprintf(&msg, format, ap);
	va_end(ap);
	if (rv == -1) {
	    (*vcsp->vcs_error_fn)(strerror(errno), vcsp->vcs_error_arg);
	    goto out;
	}
	(*vcsp->vcs_error_fn)(msg, vcsp->vcs_error_arg);
    }

out:
    free((void *)msg);
}

/*
 * vnacal_calset_alloc: alloc vnacal_calset
 *   @rows: number of VNA ports where signal is detected
 *   @columns: number of VNA ports where signal is generated
 *   @frequencies: number of frequency points
 *   @error_fn: optional error reporting function (NULL if not used)
 *   @error_arg: user data passed through to the error function (or NULL)
 *
 *   Initialize the caller-supplied vnacal_calset
 *   object and allocate and init the contained vcs_frequency_vector,
 *   vcs_matrix, vdc_sii_reference, vdc_sii_through, vdc_sji_through
 *   and vdc_sdj_leakage vectors.
 *
 *   Set errno and return NULL on error.
 */
vnacal_calset_t *vnacal_calset_alloc(const char *setname,
	int rows, int columns, int frequencies,
	vnacal_error_fn_t *error_fn, void *error_arg)
{
    vnacal_calset_t *vcsp;
    int ncells = rows *columns;

    /*
     * Allocate the vnacal_calset_t object.
     */
    vcsp = (vnacal_calset_t *)malloc(sizeof(vnacal_calset_t));
    if (vcsp == NULL) {
	if (error_fn != NULL) {
	    char buf[80];

	    (void)snprintf(buf, sizeof(buf), "vnacal_create: %s",
                    strerror(errno));
            buf[sizeof(buf)-1] = '\000';
            (*error_fn)(error_arg, buf);
	}
	return NULL;
    }
    (void)memset((void *)vcsp, 0, sizeof(vnacal_calset_t));
    vcsp->vcs_error_fn  = error_fn;
    vcsp->vcs_error_arg = error_arg;

    if (setname == NULL) {
	_vnacal_calset_error(vcsp, "invalid NULL setname");
	errno = EINVAL;
	vnacal_calset_free(vcsp);
	return NULL;
    }
    if (rows < 1 || columns < 1) {
	_vnacal_calset_error(vcsp, "vnacal_calset_alloc: "
		"invalid dimension (%d x %d)", rows, columns);
	errno = EINVAL;
	return NULL;
    }
    if (frequencies < 0) {
	_vnacal_calset_error(vcsp, "vnacal_calset_alloc: "
		"invalid frequency count (%d)", frequencies);
	errno = EINVAL;
	return NULL;
    }
    if ((vcsp->vcs_setname = strdup(setname)) == NULL) {
	vnacal_calset_free(vcsp);
	return NULL;
    }
    vcsp->vcs_rows = rows;
    vcsp->vcs_columns = columns;
    vcsp->vcs_frequencies = frequencies;
    vcsp->vcs_references[0].u.vcdsr_gamma = -1.0;
    vcsp->vcs_references[1].u.vcdsr_gamma =  1.0;
    vcsp->vcs_references[2].u.vcdsr_gamma =  0.0;
    if ((vcsp->vcs_frequency_vector = calloc(frequencies,
		    sizeof(double))) == NULL) {
	_vnacal_calset_error(vcsp, "calloc: %s", strerror(errno));
	vnacal_calset_free(vcsp);
	return NULL;
    }
    vcsp->vcs_z0 = VNADATA_DEFAULT_Z0;
    if ((vcsp->vcs_matrix = calloc(ncells,
		    sizeof(vnacal_cdata_t))) == NULL) {
	_vnacal_calset_error(vcsp, "calloc: %s", strerror(errno));
	vnacal_calset_free(vcsp);
	return NULL;
    }
    for (int i = 0; i < ncells; ++i) {
	vnacal_cdata_t *vcdp = &vcsp->vcs_matrix[i];

	for (int j = 0; j < 3; ++j) {
	    if ((vcdp->vcd_data_vectors[j] = calloc(frequencies,
			    sizeof(double complex))) == NULL) {
		_vnacal_calset_error(vcsp, "calloc: %s", strerror(errno));
		vnacal_calset_free(vcsp);
		return NULL;
	    }
	}
    }
    return vcsp;
}

/*
 * vnacal_calset_set_frequency_vector: set the frequency vector
 *   @vcsp: pointer to vnacal_calset_t
 *   @frequency_vector: vector of increasing frequencies
 */
int vnacal_calset_set_frequency_vector(vnacal_calset_t *vcsp,
	const double *frequency_vector)
{
    if (vcsp == NULL) {
	errno = EINVAL;
	return -1;
    }
    if (frequency_vector == NULL) {
	_vnacal_calset_error(vcsp, "vnacal_calset_set_frequency_vector: "
		"invalid NULL frequency_vector");
	errno = EINVAL;
	return -1;
    }
    for (int i = 0; i < vcsp->vcs_frequencies; ++i) {
	if (isnan(frequency_vector[i]) || frequency_vector[i] < 0.0) {
	    _vnacal_calset_error(vcsp, "vnacal_calset_set_frequency_vector: "
		    "invalid frequency, %f", frequency_vector[i]);
	    errno = EINVAL;
	    return -1;
	}
    }
    for (int i = 0; i < vcsp->vcs_frequencies - 1; ++i) {
	if (frequency_vector[i] >= frequency_vector[i + 1]) {
	    _vnacal_calset_error(vcsp, "vnacal_calset_set_frequency_vector: "
		    "error: frequencies not ascending");
	    errno = EINVAL;
	    return -1;
	}
    }
    (void)memcpy((void *)vcsp->vcs_frequency_vector, (void *)frequency_vector,
	    vcsp->vcs_frequencies * sizeof(double));
    vcsp->vcs_frequencies_valid = true;
    return 0;
}

/*
 * vnacal_calset_set_z0: set the system impedance for all VNA ports
 *   @vcsp: pointer to vnacal_calset_t
 *   @z0:   nominal impedance looking into a VNA port
 *
 * Note:
 *   We currently assume all VNA ports have the same system impedance.
 *   To change this, we'd probably first want to be able to set the
 *   reference gamma values on a per port basis.  Also, some changes
 *   would be needed in the calibration calculations for the "through"
 *   impedance tests which would see an impedance mismatch.
 *
 *   If not set, the default is 50 ohms.
 */
int vnacal_calset_set_z0(vnacal_calset_t *vcsp, double complex z0)
{
    if (vcsp == NULL) {
	errno = EINVAL;
	return -1;
    }
    vcsp->vcs_z0 = z0;
    return 0;
}

/*
 * vnacal_calset_add_vector: add a data vector to the calibration set
 *   @vcsp: pointer to vnacal_calset_t
 *   @row: VNA detector port (zero-based)
 *   @column: VNA driving port (zero-based)
 *   @term: calibration term
 *   @vector: vector of measured complex voltages
 */
int vnacal_calset_add_vector(vnacal_calset_t *vcsp, int row, int column,
	int term, const double complex *vector)
{
    vnacal_cdata_t *vcdp;
    double complex *destination;

    if (vcsp == NULL) {
	errno = EINVAL;
	return -1;
    }
    if (row < 0 || row >= vcsp->vcs_rows) {
	_vnacal_calset_error(vcsp, "vnacal_calset_add_vector: "
		"invalid row: %d", row);
	errno = EINVAL;
	return -1;
    }
    if (column < 0 || column >= vcsp->vcs_columns) {
	_vnacal_calset_error(vcsp, "vnacal_calset_add_vector: "
		"invalid column: %d", column);
	errno = EINVAL;
	return -1;
    }
    if (vector == NULL) {
	_vnacal_calset_error(vcsp, "vnacal_calset_add_vector: "
		"invalid NULL vector");
	errno = EINVAL;
	return -1;
    }
    if ((term & _VNACAL_DIAGONAL) && row != column) {
	_vnacal_calset_error(vcsp, "vnacal_calset_add_vector: "
		"error: diagonal term given on off-diagonal cell");
	errno = EINVAL;
	return -1;
    }
    if ((term & _VNACAL_OFF_DIAGONAL) && row == column) {
	_vnacal_calset_error(vcsp, "vnacal_calset_add_vector: "
		"error: off-diagonal term given on diagonal cell");
	errno = EINVAL;
	return -1;
    }
    term &= ~(_VNACAL_DIAGONAL | _VNACAL_OFF_DIAGONAL);
    if (term < 0 || term > 2) {
	_vnacal_calset_error(vcsp, "vnacal_calset_add_vector: "
		"invalid term");
	errno = EINVAL;
	return -1;
    }
    vcdp = VNACAL_CALIBRATION_DATA(vcsp, row, column);
    destination = vcdp->vcd_data_vectors[term];
    for (int findex = 0; findex < vcsp->vcs_frequencies; ++findex) {
	destination[findex] += vector[findex];
    }
    ++vcdp->vcd_counts[term];
    return 0;
}

/*
 * vnacal_calset_set_reference: set a scalar reference gamma value
 *   @vcsp: pointer to vnacal_calset_t
 *   @reference: reference index 0-3
 *   @gamma: gamma value, e.g. -1 short, 1 open, 0 load
 */
int vnacal_calset_set_reference(vnacal_calset_t *vcsp, int reference,
	double complex gamma)
{
    vnacal_calset_reference_t *vcdsrp;

    if (vcsp == NULL) {
	errno = EINVAL;
	return -1;
    }
    if (reference < 0 || reference > 2) {
	_vnacal_calset_error(vcsp, "vnacal_calset_set_reference: "
		"error: reference index %d not in 0..2", reference);
	errno = EINVAL;
	return -1;
    }
    vcdsrp = &vcsp->vcs_references[reference];
    if (vcdsrp->vcdsr_is_vector) {
	free((void *)vcdsrp->u.v.vcdsr_frequency_vector);
	free((void *)vcdsrp->u.v.vcdsr_gamma_vector);
	vcdsrp->u.v.vcdsr_frequency_vector = NULL;
	vcdsrp->u.v.vcdsr_gamma_vector = NULL;
    }
    vcdsrp->vcdsr_is_vector = false;
    vcdsrp->u.vcdsr_gamma = gamma;
    return 0;
}

/*
 * vnacal_calset_set_reference_vector: set a vector of reference gamma values
 *   @vcsp: pointer to vnacal_calset_t
 *   @reference: reference index 0-3
 *   @frequencies: length of frequency_vector and gamma_vector
 *   @frequency_vector: vector of increasing frequency values
 *   @gamma_vector: vector of gamma values
 *
 *   The frequency vector given to this function must span the full
 *   range of that given to vnacal_calset_set_frequency_vector, but doesn't
 *   have to be identical.
 */
int vnacal_calset_set_reference_vector(vnacal_calset_t *vcsp,
	int reference, int frequencies, const double *frequency_vector,
	const double complex *gamma_vector)
{
    double *new_frequency_vector = NULL;
    double complex *new_gamma_vector = NULL;
    vnacal_calset_reference_t *vcdsrp;

    /*
     * Validate parameters.
     */
    if (vcsp == NULL) {
	errno = EINVAL;
	return -1;
    }
    if (reference < 0 || reference > 2) {
	_vnacal_calset_error(vcsp, "vnacal_calset_set_reference_vector: "
		"error: reference index %d not in 0..2", reference);
	errno = EINVAL;
	return -1;
    }
    if (frequencies < 1) {
	_vnacal_calset_error(vcsp, "vnacal_calset_set_reference_vector: "
		"invalid number of frequencies: %d", frequencies);
	errno = EINVAL;
	return -1;
    }

    /*
     * Make sure the frequencies are ascending.
     */
    for (int i = 1; i < frequencies; ++i) {
	if (frequency_vector[i - 1] >= frequency_vector[i]) {
	    _vnacal_calset_error(vcsp,
		    "vnacal_calset_set_reference_vector: error: "
		    "frequencies not ascending");
	    errno = EINVAL;
	    return -1;
	}
    }

    /*
     * Allocate new vectors.
     */
    if ((new_frequency_vector = calloc(frequencies, sizeof(double))) == NULL) {
	_vnacal_calset_error(vcsp, "calloc: %s", strerror(errno));
	return -1;
    }
    if ((new_gamma_vector = calloc(frequencies,
		    sizeof(double complex))) == NULL) {
	_vnacal_calset_error(vcsp, "calloc: %s", strerror(errno));
	free((void *)new_frequency_vector);
	return -1;
    }
    (void)memcpy((void *)new_frequency_vector, (void *)frequency_vector,
	    frequencies * sizeof(double));
    (void)memcpy((void *)new_gamma_vector, (void *)gamma_vector,
	    frequencies * sizeof(double complex));

    /*
     * If we there are existing vectors, free them.
     */
    vcdsrp = &vcsp->vcs_references[reference];
    if (vcdsrp->vcdsr_is_vector) {
	free((void *)vcdsrp->u.v.vcdsr_frequency_vector);
	free((void *)vcdsrp->u.v.vcdsr_gamma_vector);
    }

    /*
     * Install the new vectors.
     */
    vcdsrp->u.v.vcdsr_frequencies = frequencies;
    vcdsrp->u.v.vcdsr_frequency_vector = new_frequency_vector;
    vcdsrp->u.v.vcdsr_gamma_vector = new_gamma_vector;
    vcdsrp->vcdsr_is_vector = true;
    return 0;
}

/*
 * _vnacal_calset_get_value: return the requested calibration value
 *   @vcdp: pointer to vnacal_cdata_t
 *   @term: calibration term
 *   @findex: frequency index
 */
double complex _vnacal_calset_get_value(const vnacal_cdata_t *vcdp,
	int term, int findex)
{
    term &= ~(_VNACAL_DIAGONAL | _VNACAL_OFF_DIAGONAL);
    assert(term >= 0 && term <= 2);

    if (vcdp->vcd_counts[term] == 0) {
        return 0.0;
    }
    return vcdp->vcd_data_vectors[term][findex] /
           (double)vcdp->vcd_counts[term];
}

/*
 * _vnacal_calset_get_reference: get the given reference value
 *   @vcsp: pointer to vnacal_calset_t
 *   @reference: reference index 0-3
 *   @findex: index into frequency vector
 *
 *   If the frequencies given to vnaset_calset_set_reference_vector and
 *   vnacal_calset_set_findex don't align, this function uses rational
 *   function interpolation to find the desired value.
 */
double complex _vnacal_calset_get_reference(const vnacal_calset_t *vcsp,
	int reference, int findex)
{
    const vnacal_calset_reference_t *vcdsrp;
    int segment = 0;

    assert(reference >= 0 && reference <= 2);
    assert(findex >= 0 && findex < vcsp->vcs_frequencies);
    vcdsrp = &vcsp->vcs_references[reference];

    if (!vcdsrp->vcdsr_is_vector) {
	return vcdsrp->u.vcdsr_gamma;
    }
    return _vnacal_rfi(vcdsrp->u.v.vcdsr_frequency_vector,
	    vcdsrp->u.v.vcdsr_gamma_vector, vcdsrp->u.v.vcdsr_frequencies,
	    MIN(vcdsrp->u.v.vcdsr_frequencies, VNACAL_MAX_M), &segment,
	    vcsp->vcs_frequency_vector[findex]);
}

/*
 * vnacal_calset_free: free buffers allocated by
 *      vnacal_calset_alloc
 *   @vcsp: pointer to vnacal_calset_t
 *
 *   Free vcs_frequency_vector, vcs_matrix, vdc_sii_reference, vdc_sii_through,
 *   vdc_sji_through and vdc_sdj_leakage.
 */
void vnacal_calset_free(vnacal_calset_t *vcsp)
{
    int ncells = vcsp->vcs_rows * vcsp->vcs_columns;

    if (vcsp != NULL) {
	for (int i = 0; i < 3; ++i) {
	    vnacal_calset_reference_t *vcdsrp = &vcsp->vcs_references[i];

	    if (vcdsrp->vcdsr_is_vector) {
		free((void *)vcdsrp->u.v.vcdsr_gamma_vector);
		free((void *)vcdsrp->u.v.vcdsr_frequency_vector);
	    }
	}
	if (vcsp->vcs_matrix != NULL) {
	    for (int i = 0; i < ncells; ++i) {
		vnacal_cdata_t *vcdp = &vcsp->vcs_matrix[i];

		free((void *)vcdp->vcd_data_vectors[2]);
		free((void *)vcdp->vcd_data_vectors[1]);
		free((void *)vcdp->vcd_data_vectors[0]);
	    }
	    free((void *)vcsp->vcs_matrix);
	}
	free((void *)vcsp->vcs_frequency_vector);
	free((void *)vcsp->vcs_setname);
	(void)memset((void *)vcsp, 0, sizeof(vnacal_calset_t));
    }
}
