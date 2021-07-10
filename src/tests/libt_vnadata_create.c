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

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include "libt.h"
#include "libt_vnadata.h"


/*
 * libt_vnadata_z0_names
 */
const char *libt_vnadata_z0_names[] = {
    "single",
    "real-vector",
    "complex-vector",
    "per-f"
};

/*
 * libt_vnadata_fill_names
 */
const char *libt_vnadata_fill_names[] = {
    "cell",
    "matrix",
    "vector"
};

/*
 * libt_vnadata_create: create test network parameter data
 *   @type: type of parameters to create
 *   @rows: number of rows in the matrix
 *   @columns: number of columns in the matrix
 *   @frequencies: number of frequencies
 *   @z0_type: type of reference impedances
 */
libt_vnadata_t *libt_vnadata_create(vnadata_parameter_type_t type, int rows,
	int columns, int frequencies, libt_vnadata_z0_type_t z0_type)
{
    libt_vnadata_t *tdp;
    int ports = MAX(rows, columns);

    if ((tdp = malloc(sizeof(libt_vnadata_t))) == NULL) {
	libt_error("malloc: %s", strerror(errno));
    }
    (void)memset((void *)tdp, 0, sizeof(*tdp));
    tdp->td_type = type;
    tdp->td_rows = rows;
    tdp->td_columns = columns;
    tdp->td_frequencies = frequencies;
    if ((tdp->td_vector = calloc(frequencies,
		    sizeof(double complex *))) == NULL) {
	libt_error("calloc: %s", strerror(errno));
    }
    for (int findex = 0; findex < frequencies; ++findex) {
	if ((tdp->td_vector[findex] = calloc(rows * columns,
			sizeof(double complex))) == NULL) {
	    libt_error("calloc: %s", strerror(errno));
	}
    }
    tdp->td_z0_type = z0_type;
    for (int findex = 0; findex < frequencies; ++findex) {
	for (int cell = 0; cell < rows * columns; ++cell) {
	    tdp->td_vector[findex][cell] = libt_crandn();
	}
    }
    if (frequencies > 0) {
	if ((tdp->td_frequency_vector = calloc(frequencies,
			sizeof(double))) == NULL) {
	    libt_error("calloc: %s", strerror(errno));
	}
	if (frequencies == 1) {
	    tdp->td_frequency_vector[0] = 1.0e+9;
	} else {
	    for (int i = 0; i < frequencies; ++i) {
		tdp->td_frequency_vector[i] = 1.0e+6 * pow(1.0e+3,
			(double)i / (double)(frequencies - 1));
	    }
	}
    }
    if (ports > 0) {
	switch (z0_type) {
	case Z0_SINGLE:
	case Z0_REAL_VECTOR:
	case Z0_COMPLEX_VECTOR:
	    if ((tdp->td_z0_vector = calloc(ports,
			    sizeof(double complex))) == NULL) {
		libt_error("calloc: %s", strerror(errno));
	    }
	    switch (z0_type) {
	    case Z0_SINGLE:
		for (int i = 0; i < ports; ++i) {
		    tdp->td_z0_vector[i] = 75.0;
		}
		break;

	    case Z0_REAL_VECTOR:
		for (int i = 0; i < ports; ++i) {
		    tdp->td_z0_vector[i] = (i + 1) * 10.0;
		}
		break;

	    case Z0_COMPLEX_VECTOR:
		for (int i = 0; i < ports; ++i) {
		    tdp->td_z0_vector[i] = libt_crandn();
		}
		break;

	    default:
		assert(!"unhandled case in switch");
	    }
	    break;

	case Z0_PER_F:
	    if (frequencies > 0) {
		if ((tdp->td_fz0_vector = calloc(frequencies,
				sizeof(double complex))) == NULL) {
		    libt_error("calloc: %s", strerror(errno));
		}
		for (int findex = 0; findex < frequencies; ++findex) {
		    if ((tdp->td_fz0_vector[findex] = calloc(ports,
				    sizeof(double complex))) == NULL) {
			libt_error("calloc: %s", strerror(errno));
		    }
		    for (int port = 0; port < ports; ++port) {
			tdp->td_fz0_vector[findex][port] = libt_crandn_nz();
		    }
		}
	    }
	    break;

	default:
	    assert(!"unhandled case in switch");
	}
    }
    if (opt_v >= 2) {
	(void)printf("Test data: %s %d %d %d %s\n",
		vnadata_get_type_name(type), rows, columns, frequencies,
		libt_vnadata_z0_names[z0_type]);
	for (int findex = 0; findex < frequencies; ++findex) {
	    (void)printf("f %d: %f Hz\n",
		    findex, tdp->td_frequency_vector[findex]);
	    if (tdp->td_z0_type == Z0_PER_F) {
		(void)printf("  z0:");
		for (int port = 0; port < ports; ++port) {
		    (void)printf(" %9.6f%+9.6fj",
			    creal(tdp->td_fz0_vector[findex][port]),
			    cimag(tdp->td_fz0_vector[findex][port]));
		}
		(void)printf("\n");
	    }
	    for (int row = 0; row < rows; ++row) {
		for (int column = 0; column < columns; ++column) {
		    int cell = row * columns + column;
		    double complex value = tdp->td_vector[findex][cell];

		    (void)printf("  %9.6f%+9.6fj", creal(value), cimag(value));
		}
		(void)printf("\n");
	    }
	    (void)printf("\n");
	}
	if (tdp->td_z0_type != Z0_PER_F) {
	    (void)printf("z0:");
	    for (int port = 0; port < ports; ++port) {
		(void)printf(" %9.6f%+9.6fj",
			creal(tdp->td_z0_vector[port]),
			cimag(tdp->td_z0_vector[port]));
	    }
	}
	(void)printf("\n\n");
	(void)fflush(stdout);
    }
    return tdp;
}

/*
 * libt_vnadata_free: free test data allocated in libt_vnadata_create
 *   @tdp: pointer to libt_vnadata_t structure
 */
void libt_vnadata_free(libt_vnadata_t *tdp)
{
    if (tdp != NULL) {
	if (tdp->td_z0_type == Z0_PER_F && tdp->td_fz0_vector != NULL) {
	    for (int findex = 0; findex < tdp->td_frequencies; ++findex) {
		free((void *)tdp->td_fz0_vector[findex]);
	    }
	}
	if (tdp->td_z0_type == Z0_PER_F) {
	    free((void *)tdp->td_fz0_vector);
	} else {
	    free((void *)tdp->td_z0_vector);
	}
	free((void *)tdp->td_frequency_vector);
	if (tdp->td_vector != NULL) {
	    for (int findex = 0; findex < tdp->td_frequencies; ++findex) {
		free((void *)tdp->td_vector[findex]);
	    }
	    free((void *)tdp->td_vector);
	}
	free((void *)tdp);
    }
}

/*
 * libt_vnadata_validate: check a vnadata_t structure against the test data
 *   @tdp: test parameters
 *   @vdp: vnadata structure to validate
 */
libt_result_t libt_vnadata_validate(const libt_vnadata_t *tdp,
	const vnadata_t *vdp)
{
    const vnadata_parameter_type_t type = vnadata_get_type(vdp);
    const int rows = tdp->td_rows;
    const int columns = tdp->td_columns;
    const int frequencies = tdp->td_frequencies;
    const int ports = MAX(rows, columns);
    vnadata_parameter_type_t temp_type;
    int i_temp;

    /*
     * Check type and dimensions.
     */
    temp_type = vnadata_get_type(vdp);

    if (temp_type != type) {
	libt_fail("vnadata_get_type: returned %d; expected %d\n",
		temp_type, type);
	return T_FAIL;
    }
    i_temp = vnadata_get_rows(vdp);
    if (i_temp != rows) {
	libt_fail("vnadata_get_rows: returned %d; expected %d\n",
		i_temp, rows);
	return T_FAIL;
    }
    i_temp = vnadata_get_columns(vdp);
    if (i_temp != columns) {
	libt_fail("vnadata_get_columns: returned %d; expected %d\n",
		i_temp, columns);
	return T_FAIL;
    }
    i_temp = vnadata_get_frequencies(vdp);
    if (i_temp != frequencies) {
	libt_fail("vnadata_get_frequencies: returned %d; expected %d\n",
		i_temp, frequencies);
	return T_FAIL;
    }

    /*
     * Check frequencies and data using the single cell accessors.
     */
    for (int findex = 0; findex < tdp->td_frequencies; ++findex) {
	double f_actual = vnadata_get_frequency(vdp, findex);
	double f_expected = vdp->vd_frequency_vector[findex];

	if (!libt_isequal_d_rpt("vnadata_get_frequency",
		    f_actual, f_expected)) {
	    libt_fail(": findex %d\n", findex);
	    return T_FAIL;
	}
	for (int row = 0; row < rows; ++row) {
	    for (int column = 0; column < columns; ++column) {
		int cell = row * columns + column;
		double complex value;

		value = vnadata_get_cell(vdp, findex, row, column);
		if (!libt_isequal_c_rpt("vnadata_get_cell", value,
			    tdp->td_vector[findex][cell])) {
		    libt_fail(": findex %d row %d column %d\n",
			    findex, row, column);
		    return T_FAIL;
		}
	    }
	}
    }

    /*
     * Check frequencies using the vector interface.
     */
    {
	const double *vector = vnadata_get_frequency_vector(vdp);

	for (int findex = 0; findex < tdp->td_frequencies; ++findex) {
	    if (!libt_isequal_d_rpt("vnadata_get_frequency_vector",
			vector[findex], vdp->vd_frequency_vector[findex])) {
		libt_fail(": port %d\n", findex);
		return T_FAIL;
	    }
	}
    }

    /*
     * Check data using vnadata_get_matrix.
     */
    for (int findex = 0; findex < frequencies; ++findex) {
	const double complex *matrix;

	if ((matrix = vnadata_get_matrix(vdp, findex)) == NULL &&
		rows != 0 && columns != 0) {
	    libt_fail("vnadata_get_matrix: findex %d\n", findex);
	    return T_FAIL;
	}
	for (int row = 0; row < rows; ++row) {
	    for (int column = 0; column < columns; ++column) {
		const int cell = row * columns + column;

		if (!libt_isequal_c_rpt("vnadata_get_matrix",
			    matrix[cell], tdp->td_vector[findex][cell])) {
		    libt_fail(": findex %d row %d column %d\n",
			    findex, row, column);
		    return T_FAIL;
		}
	    }
	}
    }

    /*
     * Check data using vnadata_get_to_vector.
     */
    {
	for (int row = 0; row < rows; ++row) {
	    for (int column = 0; column < columns; ++column) {
		const int cell = row * columns + column;
		double complex vector[frequencies];

		if (vnadata_get_to_vector(vdp, row, column, vector) == -1) {
		    libt_fail("vnadata_get_to_vector: row %d columns %d\n",
			    row, column);
		    return T_FAIL;
		}
		for (int findex = 0; findex < frequencies; ++findex) {
		    if (!libt_isequal_c_rpt("vnadata_get_to_vector",
				vector[findex], tdp->td_vector[findex][cell])) {
			libt_fail(": row %d column %d findex %d\n",
				row, column, findex);
			return T_FAIL;
		    }
		}
	    }
	}
    }

    /*
     * Check z0, no z0 per F.
     */
    if (tdp->td_z0_type != Z0_PER_F) {
	const double complex *vector;

	/*
	 * Check single vnadata_get_z0 and vnadata_fz0.
	 */
	for (int port = 0; port < ports; ++port) {
	    double complex value;

	    value = vnadata_get_z0(vdp, port);
	    if (!libt_isequal_c_rpt("vnadata_get_z0", value,
			tdp->td_z0_vector[port])) {
		libt_fail(": port %d\n", port);
		return T_FAIL;
	    }
	    value = vnadata_get_fz0(vdp, 0, port);
	    if (!libt_isequal_c_rpt("vnadata_get_fz0", value,
			tdp->td_z0_vector[port])) {
		libt_fail(": port %d no per-f-z0\n", port);
		return T_FAIL;
	    }
	}

	/*
	 * Check vnadata_get_z0_vector.
	 */
	if ((vector = vnadata_get_z0_vector(vdp)) == NULL && ports != 0) {
	    libt_fail("vnadata_get_z0_vector: returned NULL\n");
	    return T_FAIL;
	}
	for (int port = 0; port < ports; ++port) {
	    if (!libt_isequal_c_rpt("vnadata_get_z0_vector", vector[port],
			tdp->td_z0_vector[port])) {
		libt_fail(": port %d\n", port);
		return T_FAIL;
	    }
	}

    /*
     * Check z0 per F.
     */
    } else {
	for (int findex = 0; findex < frequencies; ++findex) {
	    const double complex *vector;

	    /*
	     * Check vnadata_get_fz0.
	     */
	    for (int port = 0; port < ports; ++port) {
		double complex value;

		value = vnadata_get_fz0(vdp, findex, port);
		if (!libt_isequal_c_rpt("vnadata_get_fz0", value,
			    tdp->td_fz0_vector[findex][port])) {
		    libt_fail(": findex %d port %d\n", findex, port);
		    return T_FAIL;
		}
	    }

	    /*
	     * Check vnadata_get_fz0_vector.
	     */
	    if ((vector = vnadata_get_fz0_vector(vdp, findex)) == NULL &&
		    ports != 0) {
		libt_fail("vnadata_get_fz0_vector: findex %d\n", findex);
		return T_FAIL;
	    }
	    for (int port = 0; port < ports; ++port) {
		if (!libt_isequal_c_rpt("vnadata_get_fz0_vector", vector[port],
			    tdp->td_fz0_vector[findex][port])) {
		    libt_fail(": findex %d port %d\n", findex, port);
		    return T_FAIL;
		}
	    }
	}
    }
    return T_PASS;
}

/*
 * libt_vnadata_fill: fill the vnadata_t structure
 *   @tdp: test parameters
 *   @vdp: vnadata structure to validate
 *   @fill_method: which method to use
 */
libt_result_t libt_vnadata_fill(const libt_vnadata_t *tdp,
	vnadata_t *vdp, libt_vnadata_fill_method_t fill_method)
{
    const vnadata_parameter_type_t type = tdp->td_type;
    const int rows = tdp->td_rows;
    const int columns = tdp->td_columns;
    const int frequencies = tdp->td_frequencies;
    const int ports = MAX(rows, columns);

    /*
     * Init the vnadata_t structure.
     */
    if (vnadata_init(vdp, type, rows, columns, frequencies) == -1) {
	libt_fail("vnadata_init: type %s rows %d columns %d frequencies %d\n",
		vnadata_get_type_name(type), rows, columns, frequencies);
	return T_FAIL;
    }

    /*
     * Load frequencies.
     */
    if (fill_method == FM_CELL) {
	for (int findex = 0; findex < frequencies; ++findex) {
	    if (vnadata_set_frequency(vdp, findex,
			tdp->td_frequency_vector[findex]) == -1) {
		libt_fail("vnadata_set_frequency: findex %d value %e\n",
			findex, tdp->td_frequency_vector[findex]);
		return T_FAIL;
	    }
	}
    } else {
	if (vnadata_set_frequency_vector(vdp, tdp->td_frequency_vector) == -1 &&
		frequencies != 0) {
	    libt_fail("vnadata_set_frequency_vector: returned -1\n");
	    return T_FAIL;
	}
    }

    /*
     * Load data
     */
    switch (fill_method) {
    case FM_CELL:
	for (int findex = 0; findex < frequencies; ++findex) {
	    for (int row = 0; row < rows; ++row) {
		for (int column = 0; column < columns; ++column) {
		    const int cell = row * columns + column;

		    if (vnadata_set_cell(vdp, findex, row, column,
				tdp->td_vector[findex][cell]) == -1) {
			libt_fail("vnadata_set_cell: "
				"findex %d row %d column %d\n",
				findex, row, column);
			return T_FAIL;
		    }
		}
	    }
	}
	break;

    case FM_MATRIX:
	for (int findex = 0; findex < frequencies; ++findex) {
	    if (vnadata_set_matrix(vdp, findex, tdp->td_vector[findex]) == -1) {
		libt_fail("vnadata_set_matrix: findex %d\n", findex);
		return T_FAIL;
	    }
	}
	break;

    case FM_VECTOR:
	{
	    double complex vector[frequencies];

	    for (int row = 0; row < rows; ++row) {
		for (int column = 0; column < columns; ++column) {
		    for (int findex = 0; findex < frequencies; ++findex) {
			const int cell = row * columns + column;

			vector[findex] = tdp->td_vector[findex][cell];
		    }
		    if (vnadata_set_from_vector(vdp, row, column,
				vector) == -1) {
			libt_fail("vnadata_set_from_vector: "
				"row %d column %d\n", row, column);
			return T_FAIL;
		    }
		}
	    }
	}
	break;

    default:
	assert(!"unhandled case in switch");
    }

    /*
     * Load z0
     */
    switch (tdp->td_z0_type) {
    case Z0_SINGLE:
	if (ports != 0) {
	    if (vnadata_set_all_z0(vdp, tdp->td_z0_vector[0]) == -1) {
		libt_fail("vnadata_set_all_z0: returned -1\n");
		return T_FAIL;
	    }
	}
	break;

    case Z0_REAL_VECTOR:
    case Z0_COMPLEX_VECTOR:
	if (fill_method == FM_CELL) {
	    for (int port = 0; port < ports; ++port) {
		if (vnadata_set_z0(vdp, port, tdp->td_z0_vector[port]) == -1) {
		    libt_fail("vnadata_set_z0: returned -1\n");
		    return T_FAIL;
		}
	    }
	} else {
	    if (vnadata_set_z0_vector(vdp, tdp->td_z0_vector) == -1) {
		libt_fail("vnadata_set_z0_vector: returned -1\n");
		return T_FAIL;
	    }
	}
	break;

    case Z0_PER_F:
	for (int findex = 0; findex < frequencies; ++findex) {
	    if (fill_method == FM_CELL) {
		for (int port = 0; port < ports; ++port) {
		    if (vnadata_set_fz0(vdp, findex, port,
				tdp->td_fz0_vector[findex][port]) == -1) {
			libt_fail("vnadata_set_fz0: returned -1\n");
			return T_FAIL;
		    }
		}
	    } else if (ports != 0) {
		if (vnadata_set_fz0_vector(vdp, findex,
			    tdp->td_fz0_vector[findex]) == -1) {
		    libt_fail("vnadata_set_fz0_vector: returned -1\n");
		    return T_FAIL;
		}
	    }
	}
	break;

    default:
	assert(!"unhandled case in switch");
    }
    return T_PASS;
}
