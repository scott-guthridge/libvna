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
 * vnacal_delete_parameter_matrix: delete the parameters in the given matrix
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @parameter_matrix: matrix of parameter indices
 *   @rows: rows in parameter_matrix
 *   @columns: columns in parameter_matrix
 *
 * Note: this function does not free the parameter matrix itself.  Use
 * free to return the memory.
 */
void vnacal_delete_parameter_matrix(vnacal_t *vcp,
	const int *parameter_matrix, int rows, int columns)
{
    for (int row = 0; row < rows; ++row) {
	for (int column = 0; column < columns; ++column) {
	    const int cell = row * columns + column;
	    int parameter = parameter_matrix[cell];

	    if (parameter >= VNACAL_PREDEFINED_PARAMETERS) {
		vnacal_delete_parameter(vcp, parameter);
	    }
	}
    }
}
