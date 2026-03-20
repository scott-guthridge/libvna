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

#include "archdep.h"

#include <assert.h>
#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"

/*
 * vnacal_eval_parameter_matrix: evaluate parameter matrix at given frequency
 *   @function: name of function user called
 *   @vpmmp: vnacal_parameter_matrix_map structure
 *   @frequency: frequency at which to evaluate
 *   @z0_vector: reference impedances results should be returned in
 *   @result_matrix: caller-allocated matrix to hold result
 */
int _vnacal_eval_parameter_matrix_i(const char *function,
        const vnacal_parameter_matrix_map_t *vpmmp, double frequency,
	const double complex *z0_vector, double complex *result_matrix)
{
    int rows = vpmmp->vpmm_rows;
    int columns = vpmmp->vpmm_columns;

    /*
     * Init result matrix to all zeros.
     */
#if BINARY_ZERO_IS_DOUBLE_ZERO
    (void)memset((void *)result_matrix, 0,
	    rows * columns * sizeof(double complex));
#else
    for (int cell = 0; cell < rows * columns; ++cell) {
	result_matrix[cell] = 0.0;
    }
#endif

    /*
     * Evaluate standards.
     */
    assert(vpmmp->vpmm_standard_rmap == NULL || z0_vector != NULL);
    for (const vnacal_standard_rmap_t *vsrmp = vpmmp->vpmm_standard_rmap;
	    vsrmp != NULL; vsrmp = vsrmp->vsrm_next) {
	vnacal_standard_t *stdp = vsrmp->vsrm_stdp;
	const int *port_map = vsrmp->vsrm_rmap_vector;
	const int std_ports = stdp->std_ports;
	double complex std_z0_vector[std_ports];
	double complex std_result_matrix[std_ports * std_ports];

	/*
	 * Fill std_z0_vector.
	 */
	for (int port = 0; port < std_ports; ++port) {
	    std_z0_vector[port] = z0_vector[port_map[port]];
	}

	/*
	 * Evaluate the standard into std_result_matrix.
	 */
	if ((*stdp->std_ops->stdo_eval)(stdp, function, std_z0_vector,
		    frequency, std_result_matrix) == -1) {
	    return -1;
	}

	/*
	 * Copy the result back to the full matrix being careful to
	 * skip rows or columns of the standard that are missing if
	 * the result matrix is rectangular.
	 */
	for (int std_row = 0; std_row < std_ports; ++std_row) {
	    int row = port_map[std_row];

	    assert(row >= 0);
	    if (row >= rows)
		continue;

	    for (int std_column = 0; std_column < std_ports; ++std_column) {
		int column = port_map[std_column];
		int std_cell;
		int cell;

		assert(column >= 0);
		if (column >= columns)
		    continue;

		std_cell = std_row * std_ports + std_column;
		cell = row * columns + column;
		result_matrix[cell] = std_result_matrix[std_cell];
	    }
	}
    }

    /*
     * Evaluate regular parameters.
     */
    for (const vnacal_parameter_rmap_t *vprmp = vpmmp->vpmm_parameter_rmap;
	    vprmp != NULL; vprmp = vprmp->vprm_next) {
	vnacal_parameter_t *vpmrp = vprmp->vprm_parameter;;
	double complex value;

	switch (vpmrp->vpmr_type) {
	case VNACAL_SCALAR:
	    value = vpmrp->vpmr_coefficient;
	    break;

	case VNACAL_VECTOR:
	case VNACAL_UNKNOWN:	/* always solved values here */
	case VNACAL_CORRELATED:
	    {
		double fmin, fmax;
		double lower, upper;

		assert(vpmrp->vpmr_frequency_vector != NULL);
		fmin = vpmrp->vpmr_frequency_vector[0];
		fmax = vpmrp->vpmr_frequency_vector[vpmrp->vpmr_frequencies-1];
		lower = (1.0 - VNACAL_F_EXTRAPOLATION) * fmin;
		upper = (1.0 + VNACAL_F_EXTRAPOLATION) * fmax;
		if (frequency < lower || frequency > upper) {
		    _vnacal_error(vpmmp->vpmm_vcp, VNAERR_USAGE,
			    "%s: frequency %e must be between %e and %e\n",
			    function, frequency, fmin, fmax);
		    return -1;
		}
		value = _vnacal_rfi(vpmrp->vpmr_frequency_vector,
			vpmrp->vpmr_coefficient_vector,
			vpmrp->vpmr_frequencies,
			MIN(vpmrp->vpmr_frequencies, VNACAL_MAX_M),
			&vpmrp->vpmr_segment,
			frequency);
	    }
	    break;

	default:
	    abort();
	}
	result_matrix[vprmp->vprm_cell] = value;
    }
    return 0;
}
