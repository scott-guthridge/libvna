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

#ifndef _TEST_VNACAL_H
#define _TEST_VNACAL_H

#include <stdbool.h>
#include "vnacal_internal.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * test_vnacal_terms_t: generated error terms
 */
typedef struct test_vnacal_terms {
    /* error term type and layout */
    vnacal_layout_t tt_layout;

    /* vector of test frequencies */
    double *tt_frequency_vector;

    /* number test frequencies */
    int tt_frequencies;

    /* vector (one per frequency) of vectors of error terms */
    double complex **tt_error_term_vector;

    /* associated vnacal_new_t structure, if not NULL */
    vnacal_new_t *tt_vnp;

} test_vnacal_terms_t;

/*
 * test_vnacal_measurements_t: simulated measurements of a standard
 */
typedef struct test_vnacal_measurements {
    double complex **tm_a_matrix;
    double complex **tm_b_matrix;
    int tm_a_rows;
    int tm_a_columns;
    int tm_b_rows;
    int tm_b_columns;
} test_vnacal_measurements_t;

/* opt_a: assert on test failure flag */
extern bool opt_a;

/* opt_v: test output verbosity */
extern int opt_v;

/* test_vnacal_sigma_n: standard deviation of noise to add to measurements */
extern double test_vnacal_sigma_n;

/* test_vnacal_sigma_t: standard deviation of tracking error to add to meas. */
extern double test_vnacal_sigma_t;

/* test_vnacal_generate_error_terms: generate random error terms */
extern test_vnacal_terms_t *test_vnacal_generate_error_terms(
	vnacal_t *vcp, vnacal_type_t type, int m_rows, int m_columns,
	int frequencies, const double *frequency_vector, double sigma,
	bool ab);

/* test_vnacal_print_error_terms: show the generated error terms */
extern void test_vnacal_print_error_terms(const test_vnacal_terms_t *ttp);

/* test_vnacal_free_error_terms: free test error terms */
extern void test_vnacal_free_error_terms(test_vnacal_terms_t *ttp);

/* test_vnacal_generate_random_parameters: generate random scalar parameters */
extern int test_vnacal_generate_random_parameters(vnacal_t *vcp, int *vector,
	int n);

/* test_vnacal_calc_needed_standards: calculate the number standards needed */
extern int test_vnacal_calc_needed_standards(vnacal_type_t type,
	int m_rows, int m_columns, bool *add_all_match);

/* make_random_calibration: make a random calibration */
extern test_vnacal_terms_t *make_random_calibration(vnacal_t *vcp,
	vnacal_type_t type, int m_rows, int m_columns, int frequencies,
	bool ab);

/* test_vnacal_print_standard: show a calibration standard */
extern void test_vnacal_print_standard(vnacal_t *vcp, const int *s,
        int s_rows, int s_columns, int frequencies,
	const double *frequency_vector, const int *port_map);

/* test_vnacal_add_single_reflect: measure a single reflect standard */
extern int test_vnacal_add_single_reflect(const test_vnacal_terms_t *ttp,
        test_vnacal_measurements_t *tmp, int s11, int port);

/* test_vnacal_add_double_reflect: measure a double reflect standard */
extern int test_vnacal_add_double_reflect(const test_vnacal_terms_t *ttp,
        test_vnacal_measurements_t *tmp, int s11, int s22,
	int port1, int port2);

/* test_vnacal_add_through: measure a through standard */
extern int test_vnacal_add_through(const test_vnacal_terms_t *ttp,
        test_vnacal_measurements_t *tmp, int port1, int port2);

/* test_vnacal_add_line: measure a line standard between the given ports */
extern int test_vnacal_add_line(const test_vnacal_terms_t *ttp,
        test_vnacal_measurements_t *tmp, const int *s_2x2,
        int port1, int port2);

/* test_vnacal_alloc_measurements: allocate test measurements */
extern test_vnacal_measurements_t *test_vnacal_alloc_measurements(
	vnacal_type_t type, int m_rows, int m_columns, int frequencies,
	bool ab);

/* test_vnacal_calculate_measurements: calculate measurements of standard */
extern int test_vnacal_calculate_measurements(
	const test_vnacal_terms_t *ttp,
        test_vnacal_measurements_t *tmp,
        const int *s_matrix, int s_matrix_rows, int s_matrix_columns,
        const int *port_map);

/* test_vnacal_print_measurements: print the measured values */
extern void test_vnacal_print_measurements(test_vnacal_measurements_t *tmp,
        int frequencies);

/* test_vnacal_free_measurements: free a test_vnacal_measurements_t structure */
extern void test_vnacal_free_measurements(test_vnacal_measurements_t *tmp);

/* test_vnacal_print_properties: print a property list */
extern void test_vnacal_print_properties(const vnaproperty_t *vprp, int indent);

/* test_vnacal_print_calibration: print solved calibration error terms */
extern void test_vnacal_print_calibration(vnacal_calibration_t *calp);

/* test_vnacal_validate_calibration: compare calculated error terms to actual */
extern int test_vnacal_validate_calibration(
	const test_vnacal_terms_t *ttp, vnacal_calibration_t *calp);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _TEST_VNACAL_H */
