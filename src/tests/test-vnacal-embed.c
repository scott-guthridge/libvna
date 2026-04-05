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
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include "vnacal_internal.h"
#include "vnaconv.h"
#include "libt_crand.h"
#include "libt.h"
#include "libt_crand.h"

#define NTRIALS		90
#define FREQUENCIES	3

/*
 * Command Line Options
 */
char *progname;
static const char options[] = "av";
static const char *const usage[] = {
    "[-av]",
    NULL
};
static const char *const help[] = {
    "-a	 abort on data miscompare",
    "-v	 show verbose output",
    NULL
};
bool opt_a = false;
int  opt_v = 0;

/*
 * TEST_EQUAL: fail the test if x and y are not equal
 *   Assumes local variable "result" and label "out" are defined.
 */
#define TEST_EQUAL(x, y, label) \
    if (opt_a) { \
        assert(libt_isequal_label((x), (y), (label))); \
    } else { \
        if (!libt_isequal_label((x), (y), (label))) { \
            result = T_FAIL; \
            goto out; \
        } \
    }

/*
 * error_fn: error reporting function
 *   @message: error message
 *   @arg: (unused)
 *   @category: error category (unused)
 */
static void error_fn(const char *message, void *arg, vnaerr_category_t category)
{
    (void)printf("%s: %s\n", progname, message);
}

/* Random number generator parameters. */
#define R_NU		0.857148
#define R_SIGMA		0.5
#define R_MIN		0.1
#define R_MAX		4.0

/*
 * random_z0: return a random reference impedance
 *     Ensure that the real part isn't too close to zero.
 */
static double complex random_z0()
{
    return libt_crand_nsmmra(R_NU, R_SIGMA, R_MIN, R_MAX, 0.0, 140.0);
}

/*
 * random_fixture_cell: generate a random fixture cell.
 *	Ensure that there are no values close to zero in the through
 *	component.
 */
static double complex random_fixture_cell(int fixture_ports,
	int row, int column)
{
    if (row + fixture_ports / 2 == column ||
	    row == column + fixture_ports / 2) {
	return libt_crand_nsmm(R_NU, R_SIGMA, R_MIN, R_MAX);
    }
    return libt_crandn();
}

/*
 * test_scalar: test the scalar embed/de-embed functions
 *   @trial: trial number
 *   @z0_type: use scalar, vector or matrix reference impedances
 */
static libt_result_t test_scalar(int trial, vnacal_z0_type_t z0_type)
{
    vnacal_t *vcp = NULL;
    double frequency_vector[FREQUENCIES];
    double complex target_z0[FREQUENCIES];
    double complex fixture_z0[FREQUENCIES][2];
    double complex result_z0[FREQUENCIES];
    double complex target_matrix[FREQUENCIES];
    double complex target_r_matrix[FREQUENCIES];
    double complex fixture_matrix[FREQUENCIES][2][2];
    double complex fixture_r_matrix[FREQUENCIES][2][2];
    double complex result_matrix[FREQUENCIES];
    double complex new_fixture_z0[2];
    int target_parameter = -1;
    int fixture_parameter_matrix[2][2];
    int embed_parameter = -1;
    int deembed_parameter = -1;
    libt_result_t result = T_FAIL;

    /*
     * Print the test header.
     */
    if (opt_v != 0) {
	(void)printf("Test vnacal embed: trial %3d scalar: z0_type %d\n",
		trial, (int)z0_type);
    }

    /*
     * Fill in random data.
     */
    for (int findex = 0; findex < FREQUENCIES; ++findex) {
	double f = 10e+9 * findex / (double)FREQUENCIES;

	frequency_vector[findex] = f;
	if (z0_type == VNACAL_Z0_MATRIX) {
	    target_z0[findex] = random_z0();
	    for (int port = 0; port < 2; ++port) {
		fixture_z0[findex][port] = random_z0();
	    }
	    result_z0[findex] = libt_crand_nsmmra(R_NU, R_SIGMA,
		    R_MIN, R_MAX, 0.0, 140.0);
	}
	target_matrix[findex] = libt_crandn();
	for (int row = 0; row < 2; ++row) {
	    for (int column = 0; column < 2; ++column) {
		fixture_matrix[findex][row][column] =
		    random_fixture_cell(2, row, column);
	    }
	}
    }

    /*
     * Renormalize as needed.
     */
    switch (z0_type) {
    case VNACAL_Z0_SCALAR:
	for (int findex = 0; findex < FREQUENCIES; ++findex) {
	    target_z0[findex] = 1.0;
	    fixture_z0[findex][0] = 1.0;
	    fixture_z0[findex][1] = 1.0;
	    result_z0[findex] = 1.0;
	}
	(void)memcpy((void *)&target_r_matrix[0],
		(void *)&target_matrix[0], sizeof(target_matrix));
	(void)memcpy((void *)&fixture_r_matrix[0][0][0],
		(void *)&fixture_matrix[0][0][0], sizeof(fixture_matrix));
	break;

    case VNACAL_Z0_VECTOR:
	target_z0[0] = random_z0();
	fixture_z0[0][0] = random_z0();
	fixture_z0[0][1] = random_z0();
	result_z0[0] = random_z0();
	for (int findex = 1; findex < FREQUENCIES; ++findex) {
	    target_z0[findex] = target_z0[0];
	    fixture_z0[findex][0] = fixture_z0[0][0];
	    fixture_z0[findex][1] = fixture_z0[0][1];
	    result_z0[findex] = result_z0[0];
	}
	new_fixture_z0[0] = result_z0[0];
	new_fixture_z0[1] = conj(result_z0[0]);
	for (int findex = 0; findex < FREQUENCIES; ++findex) {
	    vnaconv_stosrn(&target_matrix[findex], &target_r_matrix[findex],
		    &target_z0[0], &result_z0[0], 1);
	    vnaconv_stosrn(&fixture_matrix[findex][0][0],
		    &fixture_r_matrix[findex][0][0],
		    fixture_z0[0], new_fixture_z0, 2);
	}
	break;

    case VNACAL_Z0_MATRIX:
	for (int findex = 0; findex < FREQUENCIES; ++findex) {
	    new_fixture_z0[0] = result_z0[findex];
	    new_fixture_z0[1] = conj(result_z0[findex]);
	    vnaconv_stosrn(&target_matrix[findex], &target_r_matrix[findex],
		    &target_z0[findex], &result_z0[findex], 1);
	    vnaconv_stosrn(&fixture_matrix[findex][0][0],
		    &fixture_r_matrix[findex][0][0],
		    &fixture_z0[findex][0], new_fixture_z0, 2);
	}
	break;

    default:
	abort();
    }

    /*
     * Make parameter matrices, embed and evaluate.
     */
    if ((vcp = vnacal_create(error_fn, NULL)) == NULL) {
	goto out;
    }
    if ((target_parameter = vnacal_make_vector_parameter(vcp,
		frequency_vector, FREQUENCIES, target_r_matrix)) == -1) {
	goto out;
    }
    for (int row = 0; row < 2; ++row) {
	for (int column = 0; column < 2; ++column) {
	    double complex values[FREQUENCIES];

	    for (int findex = 0; findex < FREQUENCIES; ++findex) {
		values[findex] = fixture_r_matrix[findex][row][column];
	    }
	    if ((fixture_parameter_matrix[row][column] =
			vnacal_make_vector_parameter(vcp,
			    frequency_vector, FREQUENCIES, values)) == -1) {
		goto out;
	    }
	}
    }
    if ((embed_parameter = vnacal_embed_parameter(vcp, target_parameter,
		    fixture_parameter_matrix)) == -1) {
	goto out;
    }
    for (int findex = 0; findex < FREQUENCIES; ++findex) {
	if ((result_matrix[findex] = vnacal_eval_parameter(vcp,
			embed_parameter,
			frequency_vector[findex],
			target_z0[findex])) == HUGE_VAL) {
	    goto out;
	}
    }

    /*
     * Compute and check.
     */
    for (int findex = 0; findex < FREQUENCIES; ++findex) {
	double complex el;
	double complex er;
	double complex et;
	double complex em;
	double complex ts;
	double complex ti;
	double complex tx;
	double complex tm;
	double complex u;
	double complex v;
	double complex m;

	if (opt_v > 1) {
	    (void)printf("findex %d\n", findex);
	    libt_print_cmatrix("target_z0", &target_z0[findex], 1, 1);
	    libt_print_cmatrix("target", &target_matrix[findex], 1, 1);
	    libt_print_cmatrix("fixture_z0", &fixture_z0[findex][0], 2, 1);
	    libt_print_cmatrix("fixture", &fixture_matrix[findex][0][0], 2, 2);
	    if (z0_type != VNACAL_Z0_SCALAR) {
		libt_print_cmatrix("target_r", &target_r_matrix[findex], 1, 1);
		libt_print_cmatrix("fixture_r",
			&fixture_r_matrix[findex][0][0], 2, 2);
	    }
	    libt_print_cmatrix("result_z0", &result_z0[findex], 1, 1);
	    libt_print_cmatrix("result", &result_matrix[findex], 1, 1);
	}

	/*
	 * Split fixture_matrix into el, er, et and em.
	 */
	el = fixture_r_matrix[findex][0][0];
	er = fixture_r_matrix[findex][0][1];
	et = fixture_r_matrix[findex][1][0];
	em = fixture_r_matrix[findex][1][1];
	if (opt_v > 1) {
	    libt_print_cmatrix("el", &el, 1, 1);
	    libt_print_cmatrix("er", &er, 1, 1);
	    libt_print_cmatrix("et", &et, 1, 1);
	    libt_print_cmatrix("em", &em, 1, 1);
	}

	/*
	 * Convert to T parameters.
	 *   ts = er - el * et^-1 em
	 *   ti = el * et^-1
	 *   tx = -et^-1 * em
	 *   tm = et^-1
	 */
	tm = 1.0 / et;
	ti = el * tm;
	ts = er - ti * em;
	tx = -tm * em;
	if (opt_v > 1) {
	    libt_print_cmatrix("ts", &ts, 1, 1);
	    libt_print_cmatrix("ti", &ti, 1, 1);
	    libt_print_cmatrix("tx", &tx, 1, 1);
	    libt_print_cmatrix("tm", &tm, 1, 1);
	}

	/*
	 * Embed target_matrix in standard:
	 *   m = (ts * s + ti) * (tx * s + tm)^-1
	 *
	 *   u = ts * target_matrix + ti
	 *   v = tx * target_matrix + tm
	 *   m = u / v
	 */
	u = ts * target_r_matrix[findex] + ti;
	v = tx * target_r_matrix[findex] + tm;
	m = u / v;
	if (opt_v > 1) {
	    libt_print_cmatrix("m", &m, 1, 1);
	}

	/*
	 * Check
	 */
	{
	    char label[3 * sizeof(int) + 3];

	    (void)snprintf(label, sizeof(label), "E %d", findex);
	    TEST_EQUAL(result_matrix[findex], m, label);
	}
    }

    /*
     * De-embed the fixture_matrix.
     */
    if ((deembed_parameter = vnacal_deembed_parameter(vcp,
		    embed_parameter, fixture_parameter_matrix)) == -1) {
	goto out;
    }
    for (int findex = 0; findex < FREQUENCIES; ++findex) {
	{
	    char label[3 * sizeof(int) + 3];
	    double complex check;

	    if ((check = vnacal_eval_parameter(vcp, deembed_parameter,
			    frequency_vector[findex],
			    target_z0[findex])) == HUGE_VAL) {
		goto out;
	    }
	    vnaconv_stosrn(&check, &check,
		    &result_z0[findex], &target_z0[findex], 1);
	    (void)snprintf(label, sizeof(label), "D %d", findex);
	    TEST_EQUAL(check, target_matrix[findex], label);
	}
    }
    result = T_PASS;

out:
    vnacal_free(vcp);
    return result;
}

/*
 * test_matrix: test the matrix versions
 *   @trial: trial number
 *   @target_ports: number of standard/DUT ports
 *   @fixture_ports: number of ports in the test fixture
 *   @frequencies: number of frequency points
 *   @z0_type: use scalar, vector or matrix reference impedances
 */
static libt_result_t test_matrix(int trial, int target_ports, int fixture_ports,
	int frequencies, vnacal_z0_type_t z0_type)
{
    vnacal_t *vcp = NULL;
    vnadata_t *vdp_target = NULL;
    vnadata_t *vdp_target_r = NULL;
    vnadata_t *vdp_fixture = NULL;
    vnadata_t *vdp_fixture_r = NULL;
    vnadata_t *vdp_result = NULL;
    vnadata_t *vdp_check = NULL;
    int target_parameter_matrix[target_ports][target_ports];
    int fixture_parameter_matrix[fixture_ports][fixture_ports];
    int full_fixture_parameter_matrix[2 * target_ports][2 * target_ports];
    int port_map[target_ports]; /* zero-based */
    libt_result_t result = T_FAIL;

    assert(target_ports > 0);
    assert(fixture_ports > 0);
    assert(fixture_ports % 2 == 0);
    assert(fixture_ports <= target_ports * 2);

    /*
     * Print the test header.
     */
    if (opt_v != 0) {
	(void)printf("Test vnacal embed: trial %3d target_ports %d "
		"fixture_ports %d z0_type %d\n",
		trial, target_ports, fixture_ports, (int)z0_type);
    }

    /*
     * Allocate network parameter data containers.
     */
    if ((vdp_target = vnadata_alloc_and_init(error_fn, NULL, VPT_S,
		    target_ports, target_ports, frequencies)) == NULL) {
	goto out;
    }
    if ((vdp_target_r = vnadata_alloc_and_init(error_fn, NULL, VPT_S,
		    target_ports, target_ports, frequencies)) == NULL) {
	goto out;
    }
    if ((vdp_fixture = vnadata_alloc_and_init(error_fn, NULL, VPT_S,
		    fixture_ports, fixture_ports, frequencies)) == NULL) {
	goto out;
    }
    if ((vdp_fixture_r = vnadata_alloc_and_init(error_fn, NULL, VPT_S,
		    fixture_ports, fixture_ports, frequencies)) == NULL) {
	goto out;
    }
    if ((vdp_result = vnadata_alloc_and_init(error_fn, NULL, VPT_S,
		    target_ports, target_ports, frequencies)) == NULL) {
	goto out;
    }
    if ((vdp_check = vnadata_alloc_and_init(error_fn, NULL, VPT_S,
		    target_ports, target_ports, frequencies)) == NULL) {
	goto out;
    }
    for (int findex = 0; findex < frequencies; ++findex) {
	double f = 10e+9 * findex / (double)frequencies;

	if (vnadata_set_frequency(vdp_target, findex, f) == -1) {
	    goto out;
	}
	if (vnadata_set_frequency(vdp_fixture, findex, f) == -1) {
	    goto out;
	}
	if (vnadata_set_frequency(vdp_result, findex, f) == -1) {
	    goto out;
	}
	if (vnadata_set_frequency(vdp_check, findex, f) == -1) {
	    goto out;
	}
	if (z0_type == VNACAL_Z0_MATRIX) {
	    for (int port = 0; port < target_ports; ++port) {
		vnadata_set_fz0(vdp_target, findex, port, random_z0());
	    }
	    for (int port = 0; port < fixture_ports; ++port) {
		vnadata_set_fz0(vdp_fixture, findex, port, random_z0());
	    }
	    for (int port = 0; port < target_ports; ++port) {
		vnadata_set_fz0(vdp_result, findex, port, random_z0());
	    }
	}
	for (int row = 0; row < target_ports; ++row) {
	    for (int column = 0; column < target_ports; ++column) {
		vnadata_set_cell(vdp_target, findex, row, column,
			libt_crandn());
	    }
	}
	for (int row = 0; row < fixture_ports; ++row) {
	    for (int column = 0; column < fixture_ports; ++column) {
		vnadata_set_cell(vdp_fixture, findex, row, column,
			random_fixture_cell(fixture_ports, row, column));
	    }
	}
    }
    if (vnadata_set_frequency_vector(vdp_check,
		vnadata_get_frequency_vector(vdp_target)) == -1) {
	goto out;
    }
    switch (z0_type) {
    case VNACAL_Z0_SCALAR:
	vnadata_set_all_z0(vdp_target, 1.0);
	vnadata_set_all_z0(vdp_result, 1.0);
	vnadata_set_all_z0(vdp_fixture, 1.0);
	vnadata_set_all_z0(vdp_check, 1.0);
	break;

    case VNACAL_Z0_VECTOR:
	for (int port = 0; port < target_ports; ++port) {
	    vnadata_set_z0(vdp_target, port, random_z0());
	    vnadata_set_z0(vdp_result, port, random_z0());
	}
	for (int port = 0; port < fixture_ports; ++port) {
	    vnadata_set_z0(vdp_fixture, port, random_z0());
	}
	if (vnadata_set_z0_vector(vdp_check,
		vnadata_get_z0_vector(vdp_target)) == -1) {
	    goto out;
	}
	break;

    case VNACAL_Z0_MATRIX:
	for (int findex = 0; findex < frequencies; ++findex) {
	    if (vnadata_set_fz0_vector(vdp_check, findex,
		    vnadata_get_fz0_vector(vdp_target, findex)) == -1) {
		goto out;
	    }
	}
	break;

    default:
	abort();

    }

    /*
     * Half the time, create random permutation of a subset of the DUT
     * ports in port_map.  This is to demonstrate using permutations of
     * the fixture matrix and subsets of the fixture matrix, i.e. where
     * some DUT ports bypass the fixture.
     */
    for (int i = 0; i < target_ports; ++i) {
	port_map[i] = i;
    }
    if (fixture_ports < 2 * target_ports || (random() & 1)) {
	int i;

	for (i = 0; i < fixture_ports / 2; ++i) {
	    int n = target_ports - i;
	    int r;

	    if (n < 2) {
		continue;
	    }
	    r = i + random() % (target_ports - i);
	    if (r != i) {
		int temp = port_map[i];
		port_map[i] = port_map[r];
		port_map[r] = temp;
	    }
	}
	for (; i < target_ports; ++i) {	/* valid only to fixture_ports / 2 */
	    port_map[i] = -1;
	}
	if (opt_v > 1) {
	    (void)printf("port map:");
	    for (int i = 0; i < fixture_ports / 2; ++i) {
		(void)printf(" %d", port_map[i]);
	    }
	    (void)printf("\n");
	}
    } else if (opt_v > 1) {
	(void)printf("no port map.\n");
    }

    /*
     * Make the target and fixture parameter matrices.
     */
    if ((vcp = vnacal_create(error_fn, NULL)) == NULL) {
	goto out;
    }
    if (vnacal_make_data_parameter_matrix(vcp, vdp_target,
		&target_parameter_matrix[0][0],
		sizeof(target_parameter_matrix)) == -1) {
	goto out;
    }
    if (vnacal_make_data_parameter_matrix(vcp, vdp_fixture,
		&fixture_parameter_matrix[0][0],
		sizeof(fixture_parameter_matrix)) == -1) {
	goto out;
    }

    /*
     * Initialize the full fixture matrix to a no-op fixture.  Note
     * that scalar VNACAL_ONE and VNACAL_ZERO don't go through port
     * re-normalization in vnacal_eval_parameter_matrix.  We're relying
     * both on that and on the fact that the embed functions evaluate the
     * DUT in the same impedance reference context as the combined fixture
     * and DUT, so bypassed ports are already correctly normalized.
     */
    for (int row = 0; row < 2 * target_ports; ++row) {
	for (int column = 0; column < 2 * target_ports; ++column) {
	    if (row == column + target_ports ||
		    row + target_ports == column) {
		full_fixture_parameter_matrix[row][column] = VNACAL_ONE;
	    } else {
		full_fixture_parameter_matrix[row][column] = VNACAL_ZERO;
	    }
	}
    }

    /*
     * Layer in the (possibly smaller) fixture matrix.
     */
    for (int r = 0; r < fixture_ports / 2; ++r) {
	int row = port_map[r];

	for (int c = 0; c < fixture_ports / 2; ++c) {
	    int column = port_map[c];

	    full_fixture_parameter_matrix[row][column] =
		fixture_parameter_matrix[r][c];
	    full_fixture_parameter_matrix[row][column + target_ports] =
		fixture_parameter_matrix[r][c + fixture_ports / 2];
	    full_fixture_parameter_matrix[row + target_ports][column] =
		fixture_parameter_matrix[r + fixture_ports / 2][c];
	    full_fixture_parameter_matrix[row + target_ports][column +
		target_ports] = fixture_parameter_matrix[r +
		fixture_ports / 2][c + fixture_ports / 2];
	}
    }

    /*
     * Embed and export to data.
     */
    if (random() & 1) {
	int embed_parameter_matrix[target_ports][target_ports];

	if (vnacal_embed_parameter_matrix(vcp,
		    &target_parameter_matrix[0][0], target_ports,
		    &full_fixture_parameter_matrix[0][0],
		    &embed_parameter_matrix[0][0],
		    sizeof(embed_parameter_matrix)) == -1) {
	    goto out;
	}
	if (vnacal_parameter_matrix_to_data(vcp, &embed_parameter_matrix[0][0],
		    target_ports, target_ports, vdp_result) == -1) {
	    goto out;
	}
    } else {
	if (vnacal_embed(vcp, vdp_target, vdp_result,
		    &full_fixture_parameter_matrix[0][0],
		    2 * target_ports) == -1) {
	    goto out;
	}
    }

    /*
     * Compute and check.
     */
    switch (z0_type) {
    case VNACAL_Z0_SCALAR:
	vdp_target_r = vdp_target;
	vdp_fixture_r = vdp_fixture;
	break;

    case VNACAL_Z0_VECTOR:
	{
	    double complex fixture_z0_vector[fixture_ports];
	    int i;

	    if (vnadata_rconvert(vdp_target, vdp_target_r, VPT_S,
			vnadata_get_z0_vector(vdp_result),
			target_ports) == -1) {
		goto out;
	    }
	    for (i = 0; i < fixture_ports / 2; ++i) {
		const int j = port_map[i];
		double complex z0 = vnadata_get_z0(vdp_result, j);

		fixture_z0_vector[i] = z0;
		fixture_z0_vector[i + fixture_ports / 2] = conj(z0);
	    }
	    if (vnadata_rconvert(vdp_fixture, vdp_fixture_r, VPT_S,
			fixture_z0_vector, fixture_ports) == -1) {
		goto out;
	    }
	}
	break;

    case VNACAL_Z0_MATRIX:
	{
	    double complex target_fz0_matrix[frequencies][target_ports];
	    double complex fixture_fz0_matrix[frequencies][fixture_ports];

	    for (int findex = 0; findex < frequencies; ++findex) {
		int i;

		(void *)memcpy((void *)target_fz0_matrix[findex],
			(void *)vnadata_get_fz0_vector(vdp_result, findex),
			target_ports * sizeof(double complex));
		for (i = 0; i < fixture_ports / 2; ++i) {
		    const int j = port_map[i];
		    double complex z0 = vnadata_get_fz0(vdp_result, findex, j);

		    fixture_fz0_matrix[findex][i] = z0;
		    fixture_fz0_matrix[findex][i + fixture_ports / 2] =
			conj(z0);
		}
	    }
	    if (vnadata_rconvert(vdp_target, vdp_target_r, VPT_S,
			&target_fz0_matrix[0][0],
			frequencies * target_ports) == -1) {
		goto out;
	    }
	    if (vnadata_rconvert(vdp_fixture, vdp_fixture_r, VPT_S,
			&fixture_fz0_matrix[0][0],
			frequencies * fixture_ports) == -1) {
		goto out;
	    }
	}
	break;

    default:
	abort();
    }
    if (opt_v > 1) {
	libt_print_vnadata("target", vdp_target);
	libt_print_vnadata("fixture", vdp_fixture);
	if (vdp_target_r != vdp_target) {
	    libt_print_vnadata("target_r", vdp_target_r);
	}
	if (vdp_fixture_r != vdp_fixture) {
	    libt_print_vnadata("fixture_r", vdp_fixture_r);
	}
	libt_print_vnadata("result", vdp_result);
    }

    for (int findex = 0; findex < frequencies; ++findex) {
	double complex el[target_ports][target_ports];
	double complex er[target_ports][target_ports];
	double complex et[target_ports][target_ports];
	double complex em[target_ports][target_ports];
	double complex ts[target_ports][target_ports];
	double complex ti[target_ports][target_ports];
	double complex tx[target_ports][target_ports];
	double complex tm[target_ports][target_ports];
	double complex u[target_ports][target_ports];
	double complex v[target_ports][target_ports];
	double complex m[target_ports][target_ports];

	if (opt_v > 1) {
	    (void)printf("findex %d\n", findex);
	}

	/*
	 * Init the fixture S parameters to a null fixture, then layer
	 * in the actual fixture data.
	 */
	for (int row = 0; row < target_ports; ++row) {
	    for (int column = 0; column < target_ports; ++column) {
		double complex trc = row == column ? 1.0 : 0.0;

		el[row][column] = 0.0;
		er[row][column] = trc;
		et[row][column] = trc;
		em[row][column] = 0.0;
	    }
	}
	for (int row = 0; row < fixture_ports / 2; ++row) {
	    int r = port_map[row];

	    for (int column = 0; column < fixture_ports / 2; ++column) {
		int c = port_map[column];

		el[r][c] = vnadata_get_cell(vdp_fixture_r, findex,
			row, column);
		er[r][c] = vnadata_get_cell(vdp_fixture_r, findex,
			row, column + fixture_ports / 2);
		et[r][c] = vnadata_get_cell(vdp_fixture_r, findex,
			row + fixture_ports / 2, column);
		em[r][c] = vnadata_get_cell(vdp_fixture_r, findex,
			row + fixture_ports / 2, column + fixture_ports / 2);
	    }
	}
	if (opt_v > 1) {
	    libt_print_cmatrix("el", *el, target_ports, target_ports);
	    libt_print_cmatrix("er", *er, target_ports, target_ports);
	    libt_print_cmatrix("et", *et, target_ports, target_ports);
	    libt_print_cmatrix("em", *em, target_ports, target_ports);
	}

	/*
	 * Convert to T parameters.
	 *   ts = er - el * et^-1 em
	 *   ti = el * et^-1
	 *   tx = -et^-1 * em
	 *   tm = et^-1
	 */
	_vnacommon_minverse(&tm[0][0], &et[0][0], target_ports);
	_vnacommon_mmultiply(&ti[0][0], &el[0][0], &tm[0][0],
		target_ports, target_ports, target_ports);
	_vnacommon_mmultiply(&ts[0][0], &ti[0][0], &em[0][0],
		target_ports, target_ports, target_ports);
	for (int row = 0; row < target_ports; ++row) {
	    for (int column = 0; column < target_ports; ++column) {
		ts[row][column] = er[row][column] - ts[row][column];
	    }
	}
	_vnacommon_mmultiply(&tx[0][0], &tm[0][0], &em[0][0],
		target_ports, target_ports, target_ports);
	for (int row = 0; row < target_ports; ++row) {
	    for (int column = 0; column < target_ports; ++column) {
		tx[row][column] *= -1;
	    }
	}
	if (opt_v > 1) {
	    libt_print_cmatrix("ts", *ts, target_ports, target_ports);
	    libt_print_cmatrix("ti", *ti, target_ports, target_ports);
	    libt_print_cmatrix("tx", *tx, target_ports, target_ports);
	    libt_print_cmatrix("tm", *tm, target_ports, target_ports);
	}

	/*
	 * Embed target in standard:
	 *   m = (ts * s + ti) * (tx * s + tm)^-1
	 *
	 *   u = ts * target + ti;
	 *   v = tx * target + tm;
	 *   m = u / v
	 */
	_vnacommon_mmultiply(&u[0][0], &ts[0][0],
		vnadata_get_matrix(vdp_target_r, findex),
		target_ports, target_ports, target_ports);
	for (int row = 0; row < target_ports; ++row) {
	    for (int column = 0; column < target_ports; ++column) {
		u[row][column] += ti[row][column];
	    }
	}
	_vnacommon_mmultiply(&v[0][0], &tx[0][0],
		vnadata_get_matrix(vdp_target_r, findex),
		target_ports, target_ports, target_ports);
	for (int row = 0; row < target_ports; ++row) {
	    for (int column = 0; column < target_ports; ++column) {
		v[row][column] += tm[row][column];
	    }
	}
	_vnacommon_mrdivide(&m[0][0], &u[0][0], &v[0][0],
		target_ports, target_ports);
	if (opt_v > 1) {
	    libt_print_cmatrix("m", *m, target_ports, target_ports);
	}

	/*
	 * Check
	 */
	for (int row = 0; row < target_ports; ++row) {
	    for (int column = 0; column < target_ports; ++column) {
		char label[9 * sizeof(int) + 10];

		(void)snprintf(label, sizeof(label), "E %d %d %d",
			findex, row, column);
		TEST_EQUAL(m[row][column],
			vnadata_get_cell(vdp_result, findex, row, column),
			label);
	    }
	}
    }

    /*
     * De-embed the fixture.
     */
    if (random() & 1) {
	int embed_parameter_matrix[target_ports][target_ports];
	int deembed_parameter_matrix[target_ports][target_ports];

	if (vnacal_make_data_parameter_matrix(vcp, vdp_result,
		    &embed_parameter_matrix[0][0],
		    sizeof(embed_parameter_matrix)) == -1) {
	    goto out;
	}
	if (vnacal_deembed_parameter_matrix(vcp,
		    &embed_parameter_matrix[0][0], target_ports,
		    &full_fixture_parameter_matrix[0][0],
		    &deembed_parameter_matrix[0][0],
		    sizeof(deembed_parameter_matrix)) == -1) {
	    goto out;
	}
	if (vnacal_parameter_matrix_to_data(vcp,
		    &deembed_parameter_matrix[0][0],
		    target_ports, target_ports, vdp_check) == -1) {
	    goto out;
	}
    } else {
	if (vnacal_deembed(vcp, vdp_result, vdp_check,
		    &full_fixture_parameter_matrix[0][0],
		    2 * target_ports) == -1) {
	    goto out;
	}
    }
    if (opt_v > 1) {
	libt_print_vnadata("check", vdp_check);
    }
    for (int findex = 0; findex < frequencies; ++findex) {
	for (int row = 0; row < target_ports; ++row) {
	    for (int column = 0; column < target_ports; ++column) {
		char label[9 * sizeof(int) + 10];

		(void)snprintf(label, sizeof(label), "D %d %d %d",
			findex, row, column);
		TEST_EQUAL(vnadata_get_cell(vdp_target, findex, row, column),
			vnadata_get_cell(vdp_check, findex, row, column),
			label);
	    }
	}
    }
    result = T_PASS;

out:
    vnacal_free(vcp);
    vnadata_free(vdp_check);
    vnadata_free(vdp_result);
    if (vdp_fixture_r != vdp_fixture) {
	vnadata_free(vdp_fixture_r);
    }
    vnadata_free(vdp_fixture);
    if (vdp_target_r != vdp_target) {
	vnadata_free(vdp_target_r);
    }
    vnadata_free(vdp_target);
    return result;
}

/*
 * run_trials: run all test trials
 */
static libt_result_t run_trials()
{
    libt_result_t result = T_FAIL;

    for (int trial = 1; trial <= NTRIALS; ++trial) {
	for (vnacal_z0_type_t z0_type = VNACAL_Z0_SCALAR;
		z0_type <= VNACAL_Z0_MATRIX; ++z0_type) {
	    result = test_scalar(trial, z0_type);
	    if (result != T_PASS)
		goto out;
	    for (int target_ports = 1; target_ports <= 5; ++target_ports) {
		for (int fixture_ports = 2; fixture_ports <= 2 * target_ports;
			fixture_ports += 2) {
		    result = test_matrix(trial, target_ports, fixture_ports,
			    FREQUENCIES, z0_type);
		    if (result != T_PASS)
			goto out;
		}
	    }
	}
    }
    result = T_PASS;

out:
    libt_report(result);
    return result;
}

/*
 * print_usage: print a usage message and exit
 */
static void print_usage()
{
    const char *const *cpp;

    for (cpp = usage; *cpp != NULL; ++cpp) {
	(void)fprintf(stderr, "%s: usage %s\n", progname, *cpp);
    }
    for (cpp = help; *cpp != NULL; ++cpp) {
	(void)fprintf(stderr, "%s\n", *cpp);
    }
    exit(99);
}

/*
 * main
 */
int
main(int argc, char **argv)
{
    /*
     * Parse Options
     */
    if ((char *)NULL == (progname = strrchr(argv[0], '/'))) {
	progname = argv[0];
    } else {
	++progname;
    }
    for (;;) {
	switch (getopt(argc, argv, options)) {
	case -1:
	    break;

	case 'a':
	    opt_a = true;
	    continue;

	case 'v':
	    ++opt_v;
	    continue;

	default:
	    print_usage();
	}
	break;
    }
    libt_isequal_init();
    exit(run_trials());
}
