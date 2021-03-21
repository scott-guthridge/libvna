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

#include "src/archdep.h"

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include "test.h"
#include "vnacaltest.h"

/*
 * DIVROUND: divide and round up
 */
#ifndef DIVROUND
#define DIVROUND(k, n)	(((k) + (n) - 1) / (n))
#endif /* DIVROUND */

/*
 * test_vnacal_generate_random_parameters: generate n random scalar parameters
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @vector: caller-supplied vector to fill with parameter indices
 *   @n: number of parameters to generate
 */
int test_vnacal_generate_random_parameters(vnacal_t *vcp, int *vector, int n)
{
    for (int i = 0; i < n; ++i) {
	if ((vector[i] = vnacal_make_scalar_parameter(vcp,
			test_crandn())) == -1) {
	    while (--i >= 0) {
		(void)vnacal_delete_parameter(vcp, vector[i]);
	    }
	    return -1;
	}
    }
    return 0;
}

/*
 * test_vnacal_calc_needed_standards: calculate the number standards needed
 *   @type: error term type
 *   @m_rows: rows in the measurement matrix
 *   @m_columns: columns in the measurement matrix
 *   @add_all_match: set to true if all-match standard needed
 *
 * This function may sometimes overestimate T8, U8, T16 and U16 where
 * we add an extra standard.
 */
int test_vnacal_calc_needed_standards(vnacal_type_t type,
	int m_rows, int m_columns, bool *add_all_match)
{
    const int ports = MAX(m_rows, m_columns);
    int standards = 0;

    *add_all_match = false;
    if (ports == 1) {			/* special-case */
	return 3;
    }
    switch (type) {
    case VNACAL_T8:
    case VNACAL_U8:
	{
	    int terms = 2 * (m_rows + m_columns) - 1;

	    standards = DIVROUND(terms, m_rows * m_columns) + 1;
	}
	break;

    case VNACAL_TE10:
    case VNACAL_UE10:
	{
	    int terms = 2 * (m_rows + m_columns) - 1;

	    *add_all_match = true;
	    standards = DIVROUND(terms, m_rows * m_columns);
	}
	break;

    case VNACAL_T16:
    case VNACAL_U16:
	{
	    int terms = (m_rows + m_columns) * 2 * ports - 1;

	    standards = DIVROUND(terms, m_rows * m_columns) + 1;
	}
	break;

    case VNACAL_UE14:
    case _VNACAL_E12_UE14:
    case VNACAL_E12:
	{
	    int terms = m_columns * (2 * m_rows + 1);

	    *add_all_match = true;
	    standards = DIVROUND(terms, m_rows * m_columns);
	}
	break;

    default:
	abort();
    }
    return standards;
}

/*
 * make_random_calibration: make a random calibration
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @type: error term type
 *   @m_rows: number of VNA ports that detect signal
 *   @m_columns: number of VNA ports that generate signal
 *   @frequencies: number of test frequenciens
 *   @ab: true: use a, b matrices; false: use m matrix
 */
test_vnacal_terms_t *make_random_calibration(vnacal_t *vcp, vnacal_type_t type,
	int m_rows, int m_columns, int frequencies, bool ab)
{
    test_vnacal_terms_t *ttp = NULL;
    const int ports = MAX(m_rows, m_columns);
    bool add_all_match;
    int standards;
    test_vnacal_measurements_t *tmp = NULL;

    /*
     * Generate random error parameters.
     */
    if ((ttp = test_vnacal_generate_error_terms(vcp, type, m_rows, m_columns,
		    frequencies, NULL, 1.0, false)) == NULL) {
	(void)fprintf(stderr, "%s: test_vnacal_generate_error_terms: %s\n",
		progname, strerror(errno));
	goto error;
    }

    /*
     * Calculate the number of standards needed.
     */
    standards = test_vnacal_calc_needed_standards(type, m_rows, m_columns,
	    &add_all_match);

    /*
     * Allocate the measurements matrices.
     */
    if ((tmp = test_vnacal_alloc_measurements(type, m_rows, m_columns,
		    frequencies, ab)) == NULL) {
	goto error;
    }

    /*
     * If needed, add an all match matrix.
     */
    if (add_all_match) {
	int s[ports * ports];

	for (int i = 0; i < ports * ports; ++i) {
	    s[i] = VNACAL_MATCH;
	}
	if (test_vnacal_calculate_measurements(ttp, tmp, s, ports, ports,
		    NULL) == -1) {
	    goto error;
	}
	if (ab) {
	    if (vnacal_new_add_mapped_matrix(ttp->tt_vnp,
			tmp->tm_a_matrix, tmp->tm_a_rows, tmp->tm_a_columns,
			tmp->tm_b_matrix, m_rows, m_columns,
			s, ports, ports, NULL) == -1) {

		goto error;
	    }
	} else {
	    if (vnacal_new_add_mapped_matrix_m(ttp->tt_vnp,
			tmp->tm_b_matrix, m_rows, m_columns,
			s, ports, ports, NULL) == -1) {

		goto error;
	    }
	}
    }

    /*
     * Add random standards.
     */
    for (int standard = 0; standard < standards; ++standard) {
	int s[ports * ports];

	if (test_vnacal_generate_random_parameters(vcp, s, ports * ports) == -1) {
	    goto error;
	}
	if (test_vnacal_calculate_measurements(ttp, tmp, s, ports, ports,
		    NULL) == -1) {
	    goto error;
	}
	if (ab) {
	    if (vnacal_new_add_mapped_matrix(ttp->tt_vnp,
			tmp->tm_a_matrix, tmp->tm_a_rows, tmp->tm_a_columns,
			tmp->tm_b_matrix, m_rows, m_columns,
			s, ports, ports, NULL) == -1) {

		goto error;
	    }
	} else {
	    if (vnacal_new_add_mapped_matrix_m(ttp->tt_vnp,
			tmp->tm_b_matrix, m_rows, m_columns,
			s, ports, ports, NULL) == -1) {

		goto error;
	    }
	}
	for (int i = 0; i < ports * ports; ++i) {
	    if (vnacal_delete_parameter(vcp, s[i]) == -1) {
		goto error;
	    }
	}
    }
    test_vnacal_free_measurements(tmp);
    tmp = NULL;

    /*
     * Solve for the error parameters and check.
     */
    if (vnacal_new_solve(ttp->tt_vnp) == -1) {
	(void)fprintf(stderr, "%s: vnacal_solve: %s\n",
		progname, strerror(errno));
	goto error;
    }
    if (test_vnacal_validate_calibration(ttp, NULL) == -1) {
	goto error;
    }
    return ttp;

error:
    test_vnacal_free_measurements(tmp);
    test_vnacal_free_error_terms(ttp);
    return NULL;
}

/*
 * test_vnacal_print_standard: show a calibration standard
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @s: s-parameter indices matrix describing the standard
 *   @s_rows: rows in s_matrix
 *   @s_columns: columns in s_matrix
 *   @frequencies: number of calibration frequencies
 *   @frequency_vector: vector of frequencies
 *   @port_map: map from standard port to VNA port
 */
void test_vnacal_print_standard(vnacal_t *vcp, const int *s,
	int s_rows, int s_columns,
	int frequencies, const double *frequency_vector,
	const int *port_map)
{
    bool has_vector = false;

    /*
     * First scan to determine if any standards are of vector type.
     */
    for (int row = 0; row < s_rows; ++row) {
	for (int column = 0; column < s_columns; ++column) {
	    int cell = row * s_columns + column;
	    vnacal_parameter_t *vpmrp;

	    vpmrp = _vnacal_get_parameter(vcp, s[cell]);
	    for (;;) {
		switch (vpmrp->vpmr_type) {
		case VNACAL_NEW:
		    abort();
		case VNACAL_SCALAR:
		    break;
		case VNACAL_VECTOR:
		    has_vector = true;
		    break;
		case VNACAL_UNKNOWN:
		case VNACAL_CORRELATED:
		    vpmrp = vpmrp->vpmr_other;
		    continue;
		}
		break;
	    }
	}
    }

    /*
     * Print
     */
    (void)printf("standard %d x %d:\n", s_rows, s_columns);
    if (has_vector) {
	for (int findex = 0; findex < frequencies; ++findex) {
	    double f = frequency_vector[findex];

	    (void)printf("f %e\n", frequency_vector[findex]);
	    for (int row = 0; row < s_rows; ++row) {
		for (int column = 0; column < s_columns; ++column) {
		    int cell = row * s_columns + column;
		    vnacal_parameter_t *vpmrp;
		    double complex value;

		    vpmrp = _vnacal_get_parameter(vcp, s[cell]);
		    value = _vnacal_get_parameter_value_i(vpmrp, f);
		    (void)printf("  s%d%d: %8.5f%+8.5fj\n",
			row + 1, column + 1, creal(value), cimag(value));
		}
	    }
	}
    } else {
	for (int row = 0; row < s_rows; ++row) {
	    for (int column = 0; column < s_columns; ++column) {
		int cell = row * s_columns + column;
		vnacal_parameter_t *vpmrp;
		double complex value;

		vpmrp = _vnacal_get_parameter(vcp, s[cell]);
		value = _vnacal_get_parameter_value_i(vpmrp, 0.0);
		(void)printf("  s%d%d: %8.5f%+8.5fj\n",
			row + 1, column + 1, creal(value), cimag(value));
	    }
	}
    }
    if (port_map != NULL) {
	int ports = MAX(s_rows, s_columns);

	(void)printf("map:");
	for (int port = 0; port < ports; ++port) {
	    (void)printf(" %d", port_map[port]);
	}
	(void)printf("\n");
    }
    (void)printf("\n");
}

/*
 * test_vnacal_add_single_reflect: measure a single reflect standard
 *   @ttp: pointer to test error terms structure
 *   @tmp: test measurements structure
 *   @s11: reflection parameter
 *   @port: port to measure
 */
int test_vnacal_add_single_reflect(const test_vnacal_terms_t *ttp,
	test_vnacal_measurements_t *tmp, int s11, int port)
{
    const vnacal_layout_t *vlp = &ttp->tt_layout;
    const int m_rows = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    vnacal_new_t *vnp = ttp->tt_vnp;
    int a_rows = VL_HAS_COLUMN_SYSTEMS(vlp) ? 1 : m_columns;

    if (test_vnacal_calculate_measurements(ttp, tmp, &s11, 1, 1, &port) == -1) {
	return -1;
    }
    if (tmp->tm_a_matrix != NULL) {
	if (vnacal_new_add_single_reflect(vnp,
		    tmp->tm_a_matrix, a_rows, m_columns,
		    tmp->tm_b_matrix, m_rows, m_columns,
		    s11, port) == -1) {
	    return -1;
	}
    } else {
	if (vnacal_new_add_single_reflect_m(vnp,
		    tmp->tm_b_matrix, m_rows, m_columns,
		    s11, port) == -1) {
	    return -1;
	}
    }
    return 0;
}

/*
 * test_vnacal_add_double_reflect: measure a double reflect standard
 *   @ttp: pointer to test error terms structure
 *   @tmp: test measurements structure
 *   @s11: first reflection parameter
 *   @s22: second reflection parameter
 *   @port1: first port
 *   @port2: second port
 */
int test_vnacal_add_double_reflect(const test_vnacal_terms_t *ttp,
	test_vnacal_measurements_t *tmp, int s11, int s22, int port1, int port2)
{
    const vnacal_layout_t *vlp = &ttp->tt_layout;
    const int m_rows = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    vnacal_new_t *vnp = ttp->tt_vnp;
    int a_rows = VL_HAS_COLUMN_SYSTEMS(vlp) ? 1 : m_columns;
    int s_matrix[2][2] = {
	{ s11, VNACAL_ZERO },
	{ VNACAL_ZERO, s22 }
    };
    int port_map[2] = { port1, port2 };

    if (test_vnacal_calculate_measurements(ttp, tmp, &s_matrix[0][0], 2, 2,
		port_map) == -1) {
	return -1;
    }
    if (tmp->tm_a_matrix != NULL) {
	if (vnacal_new_add_double_reflect(vnp,
		    tmp->tm_a_matrix, a_rows, m_columns,
		    tmp->tm_b_matrix, m_rows, m_columns,
		    s11, s22, port1, port2) == -1) {
	    return -1;
	}
    } else {
	if (vnacal_new_add_double_reflect_m(vnp,
		    tmp->tm_b_matrix, m_rows, m_columns,
		    s11, s22, port1, port2) == -1) {
	    return -1;
	}
    }
    return 0;
}

/*
 * test_vnacal_add_through: measure a through standard between the given ports
 *   @ttp: pointer to test error terms structure
 *   @tmp: test measurements structure
 *   @port1: first port
 *   @port2: second port
 */
int test_vnacal_add_through(const test_vnacal_terms_t *ttp,
	test_vnacal_measurements_t *tmp, int port1, int port2)
{
    const vnacal_layout_t *vlp = &ttp->tt_layout;
    const int m_rows = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    vnacal_new_t *vnp = ttp->tt_vnp;
    int a_rows = VL_HAS_COLUMN_SYSTEMS(vlp) ? 1 : m_columns;
    static int s_matrix[2][2] = { { VNACAL_MATCH, VNACAL_ONE },
				  { VNACAL_ONE, VNACAL_MATCH } };
    int port_map[2] = { port1, port2 };

    if (test_vnacal_calculate_measurements(ttp, tmp,
		&s_matrix[0][0], 2, 2, port_map) == -1) {
	return -1;
    }
    if (tmp->tm_a_matrix != NULL) {
	if (vnacal_new_add_through(vnp,
		    tmp->tm_a_matrix, a_rows, m_columns,
		    tmp->tm_b_matrix, m_rows, m_columns,
		    port1, port2) == -1) {
	    return -1;
	}
    } else {
	if (vnacal_new_add_through_m(vnp,
		    tmp->tm_b_matrix, m_rows, m_columns,
		    port1, port2) == -1) {
	    return -1;
	}
    }
    return 0;
}

/*
 * test_vnacal_add_line: measure a line standard between the given ports
 *   @ttp: pointer to test error terms structure
 *   @s_2x2: 2x2 parameter matrix
 *   @port1: first port
 *   @port2: second port
 */
int test_vnacal_add_line(const test_vnacal_terms_t *ttp,
	test_vnacal_measurements_t *tmp, const int *s_2x2,
	int port1, int port2)
{
    const vnacal_layout_t *vlp = &ttp->tt_layout;
    const int m_rows = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    vnacal_new_t *vnp = ttp->tt_vnp;
    int a_rows = VL_HAS_COLUMN_SYSTEMS(vlp) ? 1 : m_columns;
    int port_map[2] = { port1, port2 };

    if (test_vnacal_calculate_measurements(ttp, tmp, s_2x2, 2, 2,
		port_map) == -1) {
	return -1;
    }
    if (tmp->tm_a_matrix != NULL) {
	if (vnacal_new_add_line(vnp,
		    tmp->tm_a_matrix, a_rows, m_columns,
		    tmp->tm_b_matrix, m_rows, m_columns,
		    s_2x2, port1, port2) == -1) {
	    return -1;
	}
    } else {
	if (vnacal_new_add_line_m(vnp,
		    tmp->tm_b_matrix, m_rows, m_columns,
		    s_2x2, port1, port2) == -1) {
	    return -1;
	}
    }
    return 0;
}
