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
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include "../vnacal_internal.h"
#include "libt.h"
#include "libt_vnacal.h"


#define NTRIALS		600

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
 * error_fn: error reporting function
 *   @category: error category
 *   @message: error message
 *   @arg: (unused)
 */
static void error_fn(vnaerr_category_t category, const char *message, void *arg)
{
    (void)fprintf(stderr, "%s: %s\n", progname, message);
}

/*
 * FREQUENCIES: number of frequency points to test
 */
#define FREQUENCIES	1/*ZZ:5*/

/*
 * standards_t
 */
typedef enum {
    GAMMA	= 0,
    P1		= 1,
    P2		= 2,
    P3		= 3,
    P4		= 4,
    P5		= 5,
    P6		= 6,
    P7		= 7,
    P8		= 8,
    P9		= 9,
    P10		= 10,
    N		= 11
} standards_t;

/*
 * Quadratic approximations of the sigma curves given in the paper.
 * With t in GHz, the sigma value is:
 *
 *    10^(coef[0] + coef[1] * t + coef[2] * t^2 - 6)
 */
static const double sigma_r_coef[] = {
    3.0143e+00,  0.0,         0.0
};
static const double sigma_st_coef[] = {
    3.2714286,   0.0778061,  -0.0029337
};
static const double sigma_l_coef[] = {
    2.59071429,  0.02817602,  0.00041454
};
static const double sigma_o_coef[] = {
    3.5842857,   0.0967602,  -0.0032526
};
static const double sigma_a_coef[] = {
    3.5442857,  -0.0113010,   0.0015944
};
static const double sigma_fl_coef[] = {
    0.67928571,  0.01437500,  0.00086097
};
static const double sigma_tr_coef[] = {
    2.9671429,   0.0486990,  -0.0017219
};

/*
 * SIGMA: calculation sigma from coefficient array and frequency, f, in GHz
 */
#define SIGMA(coef, g) \
	pow(10.0, (coef)[0] + (coef)[1] * (g) + (coef)[2] * (g)*(g) - 6.0)

/*
 * run_vnacal_van_hamme_trial: run a through-reflect-line calibration trial
 */
static libt_result_t run_vnacal_van_hamme_trial(int trial, vnacal_type_t type)
{
    libt_result_t result = T_SKIPPED;
    vnacal_t *vcp = NULL;
    libt_vnacal_terms_t *ttp = NULL;
    vnacal_new_t *vnp = NULL;
    libt_vnacal_measurements_t *tmp = NULL;
    double sigma_r[FREQUENCIES];
    double sigma_st[FREQUENCIES];
    double sigma_l[FREQUENCIES];
    double sigma_o[FREQUENCIES];
    double sigma_a[FREQUENCIES];
    double sigma_fl[FREQUENCIES];
    double sigma_tr[FREQUENCIES];
    double complex actual_values[N][FREQUENCIES];
    int actual[N];
    int unknown[N];
    int s_matrix[2][2];

    /*
     * If -v, print the test header.
     */
    if (opt_v != 0) {
	(void)printf("Test vnacal Van-hamme calibration trial %d "
		"type %-4s Van hamme\n", trial, _vnacal_type_to_name(type));
    }

    /*
     * Init actual and unknown.
     */
    for (int i = 0; i < N; ++i) {
	unknown[i] = -1;
	actual[i] = -1;
    }

    /*
     * Create the calibration structure.
     */
    if ((vcp = vnacal_create(error_fn, NULL)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_create: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }

    /*
     * Generate random error parameters.
     */
    if ((ttp = libt_vnacal_generate_error_terms(vcp, type, 2, 2,
		    FREQUENCIES, /*frequency_vector*/NULL,
		    0.1, LIBT_GET_2_10_GHZ)) == NULL) {
	(void)fprintf(stderr, "%s: libt_vnacal_generate_error_terms: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    vnp = ttp->tt_vnp;

    /*
     * Generate the sigma values and actual standard parameters.
     */
    for (int findex = 0; findex < FREQUENCIES; ++findex) {
	double g = ttp->tt_frequency_vector[findex] / 1.0e+9;

	sigma_r[findex]  = SIGMA(sigma_r_coef, g);
	sigma_st[findex] = SIGMA(sigma_st_coef, g);
	sigma_l[findex]  = SIGMA(sigma_l_coef, g);
	sigma_o[findex]  = SIGMA(sigma_o_coef, g);
	sigma_a[findex]  = SIGMA(sigma_a_coef, g);
	sigma_fl[findex] = SIGMA(sigma_fl_coef, g);
	sigma_tr[findex] = SIGMA(sigma_tr_coef, g);

	actual_values[GAMMA][findex] =  0.0 + sigma_a[findex] * libt_crandn();
	actual_values[P1][findex]    =  0.0 + sigma_r[findex] * libt_crandn();
	actual_values[P2][findex]    =  0.0 + sigma_r[findex] * libt_crandn();
	actual_values[P3][findex]    =  1.0 + sigma_st[findex] * libt_crandn();
	actual_values[P4][findex]    =  actual_values[GAMMA][findex] +
					      sigma_l[findex] * libt_crandn();
	actual_values[P5][findex]    = -1.0 + sigma_a[findex] * libt_crandn();
	actual_values[P6][findex]    =  1.0 + sigma_a[findex] * libt_crandn();
	actual_values[P7][findex]    =  actual_values[GAMMA][findex] +
					      sigma_l[findex] * libt_crandn();
	actual_values[P8][findex]    =  actual_values[P5][findex] +
					      sigma_st[findex] * libt_crandn();
	actual_values[P9][findex]    =  actual_values[P6][findex] +
					      sigma_o[findex] * libt_crandn();
	actual_values[P10][findex]   =  actual_values[GAMMA][findex] +
					      sigma_l[findex] * libt_crandn();
    }
    if (opt_v > 1) {
	(void)printf("actual:\n");
	for (int findex = 0; findex < FREQUENCIES; ++findex) {
	    (void)printf("%7.1e Hz\n", ttp->tt_frequency_vector[findex]);
	    (void)printf("  sigma_r  %e\n", sigma_r[findex]);
	    (void)printf("  sigma_st %e\n", sigma_st[findex]);
	    (void)printf("  sigma_l  %e\n", sigma_l[findex]);
	    (void)printf("  sigma_o  %e\n", sigma_o[findex]);
	    (void)printf("  sigma_a  %e\n", sigma_a[findex]);
	    (void)printf("  sigma_fl %e\n", sigma_fl[findex]);
	    (void)printf("  sigma_tr %e\n", sigma_tr[findex]);
	    for (int p = 0; p < N; ++p) {
		if (p == 0) {
		    (void)printf("  G  ");
		} else {
		    (void)printf("  P%d%s", p, p < 10 ? " " : "");
		}
		(void)printf(" %12.5e %+12.5ej    %3.5f <%8.3f\n",
			creal(actual_values[p][findex]),
			cimag(actual_values[p][findex]),
			cabs(actual_values[p][findex]),
			carg(actual_values[p][findex]) * 180.0 / M_PI);
	    }
	    (void)printf("\n");
	}
    }

    /*
     * Create the actual parameters.  These are never shown to
     * the vnacal_new_t structure; they're used only internally
     * in libt_vnacal_calculate_measurements.
     */
    for (int i = 0; i < N; ++i) {
	if ((actual[i] = vnacal_make_vector_parameter(vcp, FREQUENCIES,
		ttp->tt_frequency_vector, actual_values[i])) == -1) {
	    result = T_FAIL;
	    goto out;
	}
    }

    /*
     * Create the unknown parameters.
     */
    if ((unknown[GAMMA] = vnacal_make_correlated_parameter(vcp, VNACAL_MATCH,
		    FREQUENCIES, ttp->tt_frequency_vector,
		    sigma_a)) == -1) {
	result = T_FAIL;
	goto out;
    }
    if ((unknown[P1] = vnacal_make_correlated_parameter(vcp, VNACAL_MATCH,
		    FREQUENCIES, ttp->tt_frequency_vector,
		    sigma_r)) == -1) {
	result = T_FAIL;
	goto out;
    }
    if ((unknown[P2] = vnacal_make_correlated_parameter(vcp, VNACAL_MATCH,
		    FREQUENCIES, ttp->tt_frequency_vector,
		    sigma_r)) == -1) {
	result = T_FAIL;
	goto out;
    }
    if ((unknown[P3] = vnacal_make_correlated_parameter(vcp, VNACAL_ONE,
		    FREQUENCIES, ttp->tt_frequency_vector,
		    sigma_st)) == -1) {
	result = T_FAIL;
	goto out;
    }
    if ((unknown[P4] = vnacal_make_correlated_parameter(vcp, unknown[GAMMA],
		    FREQUENCIES, ttp->tt_frequency_vector,
		    sigma_l)) == -1) {
	result = T_FAIL;
	goto out;
    }
    if ((unknown[P5] = vnacal_make_unknown_parameter(vcp,
		    VNACAL_SHORT)) == -1) {
	result = T_FAIL;
	goto out;
    }
    if ((unknown[P6] = vnacal_make_unknown_parameter(vcp, VNACAL_OPEN)) == -1) {
	result = T_FAIL;
	goto out;
    }
    if ((unknown[P7] = vnacal_make_correlated_parameter(vcp, unknown[GAMMA],
		    FREQUENCIES, ttp->tt_frequency_vector,
		    sigma_l)) == -1) {
	result = T_FAIL;
	goto out;
    }
    if ((unknown[P8] = vnacal_make_correlated_parameter(vcp, unknown[P5],
		    FREQUENCIES, ttp->tt_frequency_vector,
		    sigma_st)) == -1) {
	result = T_FAIL;
	goto out;
    }
    if ((unknown[P9] = vnacal_make_correlated_parameter(vcp, unknown[P6],
		    FREQUENCIES, ttp->tt_frequency_vector,
		    sigma_o)) == -1) {
	result = T_FAIL;
	goto out;
    }
    if ((unknown[P10] = vnacal_make_correlated_parameter(vcp, unknown[GAMMA],
		    FREQUENCIES, ttp->tt_frequency_vector,
		    sigma_l)) == -1) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Set the measurement error.
     */
    if (vnacal_new_set_m_error(vnp, ttp->tt_frequency_vector,
		FREQUENCIES, sigma_fl, sigma_tr) == -1) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Allocate the measurements matrices.
     */
    if ((tmp = libt_vnacal_alloc_measurements(type, /*m_rows*/2, /*m_columns*/2,
		    FREQUENCIES, /*ab*/false)) == NULL) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Add the unknown through standard.
     */
    s_matrix[0][0] = actual[P1];
    s_matrix[0][1] = actual[P3];
    s_matrix[1][0] = actual[P3];
    s_matrix[1][1] = actual[P2];
    if (libt_vnacal_calculate_measurements(ttp, tmp, &s_matrix[0][0],
		2, 2, NULL) == -1) {
	result = T_FAIL;
	goto out;
    }
    s_matrix[0][0] = unknown[P1];
    s_matrix[0][1] = unknown[P3];
    s_matrix[1][0] = unknown[P3];
    s_matrix[1][1] = unknown[P2];
    if (vnacal_new_add_line_m(vnp,
		tmp->tm_b_matrix, /*m_rows*/2, /*m_columns*/2,
		&s_matrix[0][0], 1, 2) == -1) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Add load-short standard.
     */
    s_matrix[0][0] = actual[P4];
    s_matrix[0][1] = VNACAL_ZERO;
    s_matrix[1][0] = VNACAL_ZERO;
    s_matrix[1][1] = actual[P5];
    if (libt_vnacal_calculate_measurements(ttp, tmp, &s_matrix[0][0],
		2, 2, NULL) == -1) {
	result = T_FAIL;
	goto out;
    }
    if (vnacal_new_add_double_reflect_m(vnp,
		tmp->tm_b_matrix, /*m_rows*/2, /*m_columns*/2,
		unknown[P4], unknown[P5], 1, 2) == -1) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Add open-short standard.
     */
    s_matrix[0][0] = actual[P6];
    s_matrix[0][1] = VNACAL_ZERO;
    s_matrix[1][0] = VNACAL_ZERO;
    s_matrix[1][1] = actual[P5];
    if (libt_vnacal_calculate_measurements(ttp, tmp, &s_matrix[0][0],
		2, 2, NULL) == -1) {
	result = T_FAIL;
	goto out;
    }
    if (vnacal_new_add_double_reflect_m(vnp,
		tmp->tm_b_matrix, /*m_rows*/2, /*m_columns*/2,
		unknown[P6], unknown[P5], 1, 2) == -1) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Add open-load standard.
     */
    s_matrix[0][0] = actual[P6];
    s_matrix[0][1] = VNACAL_ZERO;
    s_matrix[1][0] = VNACAL_ZERO;
    s_matrix[1][1] = actual[P7];
    if (libt_vnacal_calculate_measurements(ttp, tmp, &s_matrix[0][0],
		2, 2, NULL) == -1) {
	result = T_FAIL;
	goto out;
    }
    if (vnacal_new_add_double_reflect_m(vnp,
		tmp->tm_b_matrix, /*m_rows*/2, /*m_columns*/2,
		unknown[P6], unknown[P7], 1, 2) == -1) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Add short-load standard.
     */
    s_matrix[0][0] = actual[P8];
    s_matrix[0][1] = VNACAL_ZERO;
    s_matrix[1][0] = VNACAL_ZERO;
    s_matrix[1][1] = actual[P7];
    if (libt_vnacal_calculate_measurements(ttp, tmp, &s_matrix[0][0],
		2, 2, NULL) == -1) {
	result = T_FAIL;
	goto out;
    }
    if (vnacal_new_add_double_reflect_m(vnp,
		tmp->tm_b_matrix, /*m_rows*/2, /*m_columns*/2,
		unknown[P8], unknown[P7], 1, 2) == -1) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Add short-open standard.
     */
    s_matrix[0][0] = actual[P8];
    s_matrix[0][1] = VNACAL_ZERO;
    s_matrix[1][0] = VNACAL_ZERO;
    s_matrix[1][1] = actual[P9];
    if (libt_vnacal_calculate_measurements(ttp, tmp, &s_matrix[0][0],
		2, 2, NULL) == -1) {
	result = T_FAIL;
	goto out;
    }
    if (vnacal_new_add_double_reflect_m(vnp,
		tmp->tm_b_matrix, /*m_rows*/2, /*m_columns*/2,
		unknown[P8], unknown[P9], 1, 2) == -1) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Add load-open standard.
     */
    s_matrix[0][0] = actual[P10];
    s_matrix[0][1] = VNACAL_ZERO;
    s_matrix[1][0] = VNACAL_ZERO;
    s_matrix[1][1] = actual[P9];
    if (libt_vnacal_calculate_measurements(ttp, tmp, &s_matrix[0][0],
		2, 2, NULL) == -1) {
	result = T_FAIL;
	goto out;
    }
    if (vnacal_new_add_double_reflect_m(vnp,
		tmp->tm_b_matrix, /*m_rows*/2, /*m_columns*/2,
		unknown[P10], unknown[P9], 1, 2) == -1) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Solve for the error parameters and check.
     */
    if (vnacal_new_solve(ttp->tt_vnp) == -1) {
	(void)fprintf(stderr, "%s: vnacal_solve: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (libt_vnacal_validate_calibration(ttp, NULL) == -1) {
	result = T_FAIL;
	goto out;
    }
#if 1
    for (int findex = 0; findex < FREQUENCIES; ++findex) {
	double frequency = ttp->tt_frequency_vector[findex];

	if (opt_v > 1) {
	    (void)printf("findex %d frequency %e:\n", findex, frequency);
	}
	for (int i = 0; i < N; ++i) {
	    double complex actual_value;
	    double complex solved_value;

	    actual_value = actual_values[i][findex];
	    if ((solved_value = vnacal_get_parameter_value(vcp, unknown[i],
			    frequency)) == HUGE_VAL) {
		result = T_FAIL;
		goto out;
	    }
	    if (opt_v > 1) {
		if (i == 0) {
		    (void)printf("  G  ");
		} else {
		    (void)printf("  P%d%s", i, i < 10 ? " " : "");
		}
		(void)printf(" actual % e%+ej\n",
			creal(actual_value), cimag(actual_value));
		(void)printf("      solved % e%+ej\n",
			creal(solved_value), cimag(solved_value));
		(void)printf("      delta  % e\n",
			cabs(solved_value - actual_value));
	    }
	    if (!libt_isequal(solved_value, actual_value)) {
		if (opt_a) {
		    assert(!"data miscompare");
		}
		result = T_FAIL;
		goto out;
	    }
	}
    }
#endif
    result = T_PASS;

out:
    for (int i = 0; i < N; ++i) {
	if (actual[i] != -1) {
	    vnacal_delete_parameter(vcp, actual[i]);
	}
	if (unknown[i] != -1) {
	    vnacal_delete_parameter(vcp, unknown[i]);
	}
    }
    libt_vnacal_free_measurements(tmp);
    libt_vnacal_free_error_terms(ttp);
    vnacal_free(vcp);
    return result;
}

/*
 * test_vnacal_van_hamme: through-reflect-line calibration
 */
static libt_result_t test_vnacal_van_hamme()
{
    libt_result_t result = T_SKIPPED;
    static vnacal_type_t type_array[] = {
	VNACAL_T16, VNACAL_U16
    };

    for (int trial = 0; trial < NTRIALS; ++trial) {
	for (int t_index = 0; t_index < 2; ++t_index) {
	    result = run_vnacal_van_hamme_trial(trial, type_array[t_index]);
	    if (result != T_PASS)
		goto out;
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
    libt_isequal_eps = 0.1;
    exit(test_vnacal_van_hamme());
}
