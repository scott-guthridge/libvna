/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2026 D Scott Guthridge <scott_guthridge@rompromity.net>
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
 * vnacal_parameter_forw_map_t: forward map for standard parameters
 */
typedef struct vnacal_parameter_forw_map {
    /* standard reverse map entry for this port */
    vnacal_standard_rmap_t *vpfm_vsrmp;

    /* port of standard corresponding to the parameter matrix port, -1 for
       not set, or -2 for regular parameter found in this row or column */
    int vpfm_standard_port;

    /* cell of parameter matrix that created this mapping */
    int vpfm_origin_cell;

} vnacal_parameter_forw_map_t;

/*
 * report_port_conflict: report inconsistent mapping
 *   @function: name of user-called function
 *   @vcp: vnacal_t structure
 *   @matrix: parameter matrix
 *   @columns: number of columns in matrix
 *   @row: row of entry found in conflict
 *   @column: column of entry found in conflict
 *   @origin_cell: cell of parameter matrix with which this entry conficts
 */
static void report_port_conflict(const char *function, vnacal_t *vcp,
    vnacal_parameter_t **matrix, int columns, int row, int column,
    int origin_cell)
{
    char parameter_name1[PARAMETER_BUFFER_ALLOC];
    char parameter_name2[PARAMETER_BUFFER_ALLOC];
    char buf1[SXX_BUFFER_ALLOC];
    char buf2[SXX_BUFFER_ALLOC];

    _vnacal_get_parameter_name(matrix[row * columns + column],
	    /*with_sxx=*/true, parameter_name1);
    _vnacal_get_parameter_name(matrix[origin_cell],
	    /*with_sxx=*/true, parameter_name2);
    _vnacal_format_sxx(buf1, &buf1[SXX_BUFFER_ALLOC], row, column);
    _vnacal_format_sxx(buf2, &buf2[SXX_BUFFER_ALLOC],
	    origin_cell / columns, origin_cell % columns);
    _vnacal_error(vcp, VNAERR_USAGE,
	    "%s: %s at %s in parameter matrix conflicts with %s at %s",
	    function, parameter_name1, buf1, parameter_name2, buf2);
}

/*
 * _vnacal_get_calkit_name: return calkit name and number of ports
 *   @vcdp: vnacal_calkit_data_t structure
 *   @ip_ports: address of int to receive number of ports
 */
const char *_vnacal_get_calkit_name(const vnacal_calkit_data_t *vcdp,
	int *ip_ports)
{
    switch (vcdp->vcd_type) {
    case VNACAL_CALKIT_SHORT:
	*ip_ports = 1;
	return "calkit short";

    case VNACAL_CALKIT_OPEN:
	*ip_ports = 1;
	return "calkit open";

    case VNACAL_CALKIT_LOAD:
	*ip_ports = 1;
	return "calkit load";

    case VNACAL_CALKIT_THROUGH:
	*ip_ports = 2;
	return "calkit through";

    default:
	break;
    }
    return NULL;
}

/*
 * _vnacal_fill_standard_parameter_matrix: fill parameter matrix
 *   @function: user-called function
 *   @stdp: pointer to standard
 *   @parameter_matrix: caller-supplied parameter matrix
 *
 *   Note: parameter_matrix must be initialized to -1's.
 */
int _vnacal_fill_standard_parameter_matrix(const char *function,
	vnacal_standard_t *stdp, int *parameter_matrix)
{
    const vnacal_parameter_type_t type = stdp->std_ops->stdo_type;
    vnacal_t *vcp = stdp->std_vcp;
    const int ports = stdp->std_ports;

    for (int row = 0; row < ports; ++row) {
	for (int column = 0; column < ports; ++column) {
	    const int cell = ports * row + column;
	    vnacal_parameter_t *vpmrp;

	    vpmrp = _vnacal_alloc_parameter(function, vcp);
	    if (vpmrp == NULL) {
		goto error;
	    }
	    ++stdp->std_refcount;
	    vpmrp->vpmr_type = type;
	    vpmrp->vpmr_stdp = stdp;
	    vpmrp->vpmr_row = row;
	    vpmrp->vpmr_column = column;
	    parameter_matrix[cell] = vpmrp->vpmr_index;
	}
    }
    return 0;

error:
    vnacal_delete_parameter_matrix(stdp->std_vcp, parameter_matrix,
	    ports, ports);
    return -1;
}

/*
 * _vnacal_analyze_parameter_matrix: check matrix and build standard maps
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @matrix: parameter matrix (caller must maintain holds)
 *   @rows: number of rows in parameter matrix
 *   @columns: number of columns in parameter matrix
 *   @initial: for unknown parameters, use initial rather than solved value
 */
vnacal_parameter_matrix_map_t *_vnacal_analyze_parameter_matrix(
	const char *function, vnacal_t *vcp, vnacal_parameter_t **matrix,
	int rows, int columns, bool initial)
{
    const int ports = MAX(rows, columns);
    vnacal_parameter_forw_map_t *vpfmp_vector = NULL;
    vnacal_parameter_matrix_map_t *vpmmp_result = NULL;
    vnacal_standard_rmap_t **standard_rmap_anchor;
    vnacal_parameter_rmap_t **parameter_rmap_anchor;


    /*
     * If the initial flag is set, resolve all unknown parameters to
     * their initial values.  If not set, test that the parameter has
     * a solved value.
     */
    if (initial) {
	for (int cell = 0; cell < rows * columns; ++cell) {
	    vnacal_parameter_t **vpmrpp = &matrix[cell];

	    if (*vpmrpp == NULL)
		continue;

	    while ((*vpmrpp)->vpmr_type == VNACAL_UNKNOWN ||
	           (*vpmrpp)->vpmr_type == VNACAL_CORRELATED) {
		*vpmrpp = (*vpmrpp)->vpmr_other;
	    }
	}
    } else {
	for (int cell = 0; cell < rows * columns; ++cell) {
	    vnacal_parameter_t **vpmrpp = &matrix[cell];

	    if (*vpmrpp == NULL)
		continue;

	    if ((*vpmrpp)->vpmr_type != VNACAL_UNKNOWN &&
		(*vpmrpp)->vpmr_type != VNACAL_CORRELATED)
		continue;

	    if (matrix[cell]->vpmr_frequency_vector == NULL) {
		int row = cell / columns;
		int column = cell % columns;

		_vnacal_error(vcp, VNAERR_USAGE,
		    "%s: unknown parameter at s%d%s%d has no solved value",
		    function, row + 1, (row >= 9 || column >= 9) ? "_" : "",
		    column + 1);
		return NULL;
	    }
	}
    }

    /*
     * Allocate and init the forward port map.
     */
    if ((vpfmp_vector = calloc(ports,
		    sizeof(vnacal_parameter_forw_map_t))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"calloc: %s", strerror(errno));
	return NULL;
    }
    for (int i = 0; i < ports; ++i) {
	vpfmp_vector[i].vpfm_standard_port = -1;
	vpfmp_vector[i].vpfm_origin_cell = -1;
    }

    /*
     * Allocate and init the parameter matrix port map structure.
     */
    vpmmp_result = malloc(sizeof(vnacal_parameter_matrix_map_t));
    if (vpmmp_result == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "malloc: %s", strerror(errno));
	goto error;
    }
    (void)memset((void *)vpmmp_result, 0, sizeof(*vpmmp_result));
    vpmmp_result->vpmm_vcp = vcp;
    vpmmp_result->vpmm_rows = rows;
    vpmmp_result->vpmm_columns = columns;
    standard_rmap_anchor = &vpmmp_result->vpmm_standard_rmap;
    parameter_rmap_anchor = &vpmmp_result->vpmm_parameter_rmap;

    /*
     * Analyze each cell of the parameter matrix.
     */
    for (int row = 0; row < rows; ++row) {
	vnacal_parameter_forw_map_t *row_fmap = &vpfmp_vector[row];

	for (int column = 0; column < columns; ++column) {
	    int cell = columns * row + column;
	    vnacal_parameter_forw_map_t *column_fmap = &vpfmp_vector[column];
	    vnacal_parameter_t *vpmrp = matrix[cell];

	    /*
	     * Ignore NULL and zero entries.
	     */
	    if (vpmrp == NULL) {
		continue;
	    }
	    if (vpmrp->vpmr_index == VNACAL_ZERO) {
		continue;
	    }
	    if (VNACAL_IS_STANDARD_PARAMETER(vpmrp)) {
		vnacal_standard_t *stdp = vpmrp->vpmr_stdp;
		int r = vpmrp->vpmr_row;
		int c = vpmrp->vpmr_column;
		int origin_cell = -1;
		vnacal_standard_rmap_t *vsrmp = NULL;

                /*
		 * Check for off-diagonal standard elements in the major
		 * diagonal of the parameter matrix:
		 *   s12 ...
		 *   ... ...
		 */
		if (r != c && row == column) {
		    char buf1[SXX_BUFFER_ALLOC];
		    char buf2[SXX_BUFFER_ALLOC];

		    _vnacal_format_sxx(buf1, &buf1[SXX_BUFFER_ALLOC], r, c);
		    _vnacal_format_sxx(buf2, &buf2[SXX_BUFFER_ALLOC],
			    row, column);
		    _vnacal_error(vcp, VNAERR_USAGE,
			     "%s: off-diagonal element %s of %s standard "
			     "cannot appear in diagonal element %s of "
			     "parameter matrix",
			     function, buf1, stdp->std_name,
			     buf2);
		    goto error;
		}
		/*
		 * Check for diagonal standard elements outside of the
		 * major diagonal of the parameter matrix:
		 *   ... s11
		 *   ... ...
		 */
		if (r == c && row != column) {
		    char buf1[SXX_BUFFER_ALLOC];
		    char buf2[SXX_BUFFER_ALLOC];

		    _vnacal_format_sxx(buf1, &buf1[SXX_BUFFER_ALLOC], r, c);
		    _vnacal_format_sxx(buf2, &buf2[SXX_BUFFER_ALLOC],
			    row, column);
		    _vnacal_error(vcp, VNAERR_USAGE,
			    "%s: diagonal element %s of %s standard "
			    "cannot appear in off-diagonal element %s of "
			    "parameter matrix",
			    function, buf1, stdp->std_name, buf2);
		    goto error;
		}
		/*
		 * Test if we already have a vnacal_standard_rmap_t
		 * for this row or column.  If not, create one.
		 */
		if (row_fmap->vpfm_vsrmp != NULL) {
		    vsrmp = row_fmap->vpfm_vsrmp;
		    origin_cell = row_fmap->vpfm_origin_cell;
		} else if (column_fmap->vpfm_vsrmp != NULL) {
		    vsrmp = column_fmap->vpfm_vsrmp;
		    origin_cell = column_fmap->vpfm_origin_cell;
		} else {
		    int std_ports = stdp->std_ports;

		    vsrmp = malloc(sizeof(vnacal_standard_rmap_t));
		    if (vsrmp == NULL) {
			_vnacal_error(vcp, VNAERR_SYSTEM,
				"malloc: %s", strerror(errno));
			goto error;
		    }
		    (void)memset((void *)vsrmp, 0, sizeof(*vsrmp));
		    vsrmp->vsrm_stdp = stdp;
		    vsrmp->vsrm_rmap_vector = malloc(std_ports * sizeof(int));
		    if (vsrmp->vsrm_rmap_vector == NULL) {
			free((void *)vsrmp);
			_vnacal_error(vcp, VNAERR_SYSTEM,
				"malloc: %s", strerror(errno));
			goto error;
		    }
		    vsrmp->vsrm_cell_vector = malloc(std_ports * sizeof(int));
		    if (vsrmp->vsrm_cell_vector == NULL) {
			free((void *)vsrmp);
			_vnacal_error(vcp, VNAERR_SYSTEM,
				"malloc: %s", strerror(errno));
			goto error;
		    }
		    for (int i = 0; i < std_ports; ++i) {
			vsrmp->vsrm_rmap_vector[i] = -1;
			vsrmp->vsrm_cell_vector[i] = -1;
		    }
		    *standard_rmap_anchor = vsrmp;
		    standard_rmap_anchor = &vsrmp->vsrm_next;
		}
		/*
		 * Check for elements of more than one standard in the
		 * same row or column of the parameter matrix.
		 */
		if (vsrmp->vsrm_stdp != stdp) {
		    report_port_conflict(function, vcp, matrix, columns,
			    row, column, origin_cell);
		    goto error;
		}
		/*
		 * Fill in map from parameter port to standard port for
		 * this row.  If aleady set, check for inconsistency:
		 *   s11 s21
		 *        ^
		 */
		if (row_fmap->vpfm_standard_port == -1) {
		    row_fmap->vpfm_vsrmp = vsrmp;
		    row_fmap->vpfm_standard_port = r;
		    row_fmap->vpfm_origin_cell = cell;
		} else if (row_fmap->vpfm_standard_port != r) {
		    report_port_conflict(function, vcp, matrix, columns,
			    row, column, row_fmap->vpfm_origin_cell);
		    goto error;
		}
		/*
		 * Fill in map from parameter port to standard port for
		 * this column.  If aleady set, check for inconsistency:
		 *   s11
		 *   s12
		 *     ^
		 */
		if (column_fmap->vpfm_standard_port == -1) {
		    column_fmap->vpfm_vsrmp = vsrmp;
		    column_fmap->vpfm_standard_port = c;
		    column_fmap->vpfm_origin_cell = cell;
		} else if (column_fmap->vpfm_standard_port != c) {
		    report_port_conflict(function, vcp, matrix, columns,
			    row, column, column_fmap->vpfm_origin_cell);
		    goto error;
		}
		/*
		 * Fill in reverse map from standard to parameter matrix
		 * for this row.  If already set, check for inconsistency:
		 *   s11	or s11 s12 s12
		 *       s11                 ^
		 */
		if (vsrmp->vsrm_rmap_vector[r] == -1) {
		    vsrmp->vsrm_rmap_vector[r] = row;
		    vsrmp->vsrm_cell_vector[r] = cell;
		} else if (vsrmp->vsrm_rmap_vector[r] != row) {
		    report_port_conflict(function, vcp, matrix, columns,
			    row, column, vsrmp->vsrm_cell_vector[r]);
		    goto error;
		}
		/*
		 * Fill in reverse map from standard to parameter
		 * matrix for this column.  If already set, check for
		 * inconsistency:
		 *   s11
		 *   s21
		 *   s21
		 *    ^
		 */
		if (vsrmp->vsrm_rmap_vector[c] == -1) {
		    vsrmp->vsrm_rmap_vector[c] = column;
		    vsrmp->vsrm_cell_vector[c] = cell;
		} else if (vsrmp->vsrm_rmap_vector[c] != column) {
		    report_port_conflict(function, vcp, matrix, columns,
			    row, column, vsrmp->vsrm_cell_vector[c]);
		    goto error;
		}
	    } else {
		vnacal_parameter_rmap_t *vprmp = NULL;

 		/*
		 * Check for a row of the parameter matrix containing
		 * both standard parameters and regular parameters:
		 *   s11 s12 r11<---
		 *   s21 s22
		 */
		if (row_fmap->vpfm_vsrmp != NULL) {
		    report_port_conflict(function, vcp, matrix, columns,
			    row, column, row_fmap->vpfm_origin_cell);
		    goto error;
		}
		/*
		 * Check for a column of the parameter matrix containing
		 * both standard parameters and regular parameters:
		 *   s11 s12
		 *   s21 s22
		 *   r11
		 *   ^^^
		 */
		if (column_fmap->vpfm_vsrmp != NULL) {
		    report_port_conflict(function, vcp, matrix, columns,
			    row, column, column_fmap->vpfm_origin_cell);
		    goto error;
		}
		/*
		 * Set vpfm_standard_port for this row and column to -2 to
		 * indicate that the port contains regular parameters.
		 */
		if (row_fmap->vpfm_origin_cell == -1) {
		    row_fmap->vpfm_standard_port = -2;
		    row_fmap->vpfm_origin_cell = cell;
		}
		if (column_fmap->vpfm_origin_cell == -1) {
		    column_fmap->vpfm_standard_port = -2;
		    column_fmap->vpfm_origin_cell = cell;
		}
		/*
		 * Create a vnacal_parameter_rmap_t entry for this
		 * parameter.
		 */
		vprmp = malloc(sizeof(vnacal_parameter_rmap_t));
		if (vprmp == NULL) {
		    _vnacal_error(vcp, VNAERR_SYSTEM,
			    "malloc: %s", strerror(errno));
		    goto error;
		}
		(void)memset((void *)vprmp, 0, sizeof(*vprmp));
		vprmp->vprm_parameter = vpmrp;
		vprmp->vprm_cell = cell;
		*parameter_rmap_anchor = vprmp;
		parameter_rmap_anchor = &vprmp->vprm_next;
	    }
	}
    }

    /*
     * Check for ports of the standard that don't appear in the
     * parameter matrix:
     *   s11 0
     *   0   0
     *
     * This check is necessary because we need to know the reference
     * impedances of all ports.  We do allow missing rows, missing
     * columns and missing individual cells, though, as long as all
     * ports are covered:
     *   s11 s12    s11 0    s11 s12    s11 s12
     *   0   0      s21 0    0   s22    s21 0
     */
    for (vnacal_standard_rmap_t *vsrmp = vpmmp_result->vpmm_standard_rmap;
	    vsrmp != NULL; vsrmp = vsrmp->vsrm_next) {
	int std_ports = vsrmp->vsrm_stdp->std_ports;
	for (int i = 0; i < std_ports; ++i) {
	    if (vsrmp->vsrm_rmap_vector[i] < 0) {
		int origin = -1;
		char buf[SXX_BUFFER_ALLOC];

		for (int j = 0; j < std_ports; ++j) {
		    int temp;

		    if ((temp = vsrmp->vsrm_rmap_vector[j]) >= 0 &&
			    (origin == -1 || temp < origin)) {
			origin = temp;
		    }
		}
		assert(origin >= 0);
		_vnacal_format_sxx(buf, &buf[SXX_BUFFER_ALLOC],
			origin / std_ports, origin % std_ports);
		_vnacal_error(vcp, VNAERR_USAGE,
			"%s: no elements of port %d of the %s standard "
			"at %s appear in the parameter matrix",
			function, i + 1, vsrmp->vsrm_stdp->std_name,
			buf);
		goto error;
	    }
	}
    }

out:
    free((void *)vpfmp_vector);
    return vpmmp_result;

error:
    _vnacal_free_parameter_matrix_map(vpmmp_result);
    vpmmp_result = NULL;
    goto out;
}

/*
 * _vnacal_free_parameter_matrix_map: free a vnacal_parameter_matrix_map_t
 *   @vpmmp: struct returned by _vnacal_analyze_parameter_matrix
 */
void _vnacal_free_parameter_matrix_map(vnacal_parameter_matrix_map_t *vpmmp)
{
    vnacal_standard_rmap_t *vsrmp;
    vnacal_parameter_rmap_t *vprmp;

    if (vpmmp == NULL)
	return;

    while ((vsrmp = vpmmp->vpmm_standard_rmap) != NULL) {
	vpmmp->vpmm_standard_rmap = vsrmp->vsrm_next;

	free((void *)vsrmp->vsrm_cell_vector);
	free((void *)vsrmp->vsrm_rmap_vector);
	free((void *)vsrmp);
    }
    while ((vprmp = vpmmp->vpmm_parameter_rmap) != NULL) {
	vpmmp->vpmm_parameter_rmap = vprmp->vprm_next;

	free((void *)vprmp);
    }
    free((void *)vpmmp);
}

/*
 * _vnacal_alloc_standard: allocate a vnacal_standard_t structure
 *   @function: name of user-called function
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ops: type and class-specific operations
 *   @ports: number of ports, assumed square
 *   @size: size of derived subtype in bytes
 *
 * Return:
 *   Derived standard structure with reference count of 1.  When
 *   the reference count becomes zero, stdo_free is called to free
 *   subclass-specific data, then the structure is freed.
 */
void *_vnacal_alloc_standard(const char *function, vnacal_t *vcp,
	const vnacal_standard_ops_t *ops, int ports, size_t size)
{
    vnacal_standard_t *stdp;

    assert(size >= sizeof(vnacal_standard_t));
    if ((stdp = malloc(size)) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "malloc: %s", strerror(errno));
	return NULL;
    }
    (void)memset((void *)stdp, 0, size);
    stdp->std_ops = ops;
    stdp->std_name = NULL;
    stdp->std_ports = ports;
    stdp->std_refcount = 1;
    stdp->std_vcp = vcp;

    return stdp;
}

/*
 * _vnacal_free_standard: free a vnacal_standard_t structure
 *   @stdp: standard
 */
static void _vnacal_free_standard(vnacal_standard_t *stdp)
{
    if (stdp != NULL) {
	const vnacal_standard_ops_t *stdop = stdp->std_ops;

	assert(stdp->std_refcount == 0);
	if (stdop->stdo_free != NULL) {
	    stdop->stdo_free(stdp);
	}
	(void)free((void *)stdp->std_name);
	free((void *)stdp);
    }
}

/*
 * _vnacal_release_standard: decrement the reference count on a standard
 *   @stdpp: address of pointer to standard
 *
 *  Frees standard and sets *stdpp to NULL if the reference count goes
 *  to zero.
 */
void _vnacal_release_standard(vnacal_standard_t **stdpp)
{
    vnacal_standard_t *stdp = *stdpp;

    assert(stdp->std_refcount > 0);
    if (--stdp->std_refcount == 0) {
	_vnacal_free_standard(stdp);
	*stdpp = stdp = NULL;
    }
}
