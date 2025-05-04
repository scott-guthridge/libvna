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
 * vnacal_eval_parameter_matrix: evaluate parameter matrix at a given frequency
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @parameter_matrix: index of parameter
 *   @ports: number of rows and nuumber of columns in paramter_matrix
 *   @frequency: frequency at which to evaluate the parameter
 *   @z0_vector: reference impedance for each port of the result
 *   @result_matrix: caller-supplied matrix to hold the result
 */
int vnacal_eval_parameter_matrix(vnacal_t *vcp,
	const int *parameter_matrix, int rows, int columns, double frequency,
	const double complex *z0_vector, double complex *result_matrix)
{
    vnacal_parameter_t **matrix = NULL;
    vnacal_parameter_matrix_map_t *vpmmp = NULL;
    int rc = -1;

    matrix = calloc(rows * columns, sizeof(vnacal_parameter_t *));
    if (matrix == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s", strerror(errno));
	goto out;
    }
    for (int cell = 0; cell < rows * columns; ++cell) {
	vnacal_parameter_t *vpmrp;

	vpmrp = _vnacal_get_parameter(vcp, parameter_matrix[cell]);
	if (vpmrp == NULL) {
	    goto out;
	}
	matrix[cell] = vpmrp;
    }
    if ((vpmmp = _vnacal_analyze_parameter_matrix(__func__, vcp,
		    matrix, rows, columns, /*initial=*/false)) == NULL) {
	goto out;
    }
    rc = _vnacal_eval_parameter_matrix_i(__func__, vpmmp, frequency,
	    z0_vector, result_matrix);

out:
    _vnacal_free_parameter_matrix_map(vpmmp);
    free((void *)matrix);
    return rc;
}
