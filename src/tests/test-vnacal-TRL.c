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


#define NTRIALS		1000

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
 * TRL_FREQUENCIES: number of frequency points to test
 */
#define TRL_FREQUENCIES	2

/*
 * make_random_parameters: make random actual and guess
 *   @r_actual: address to receive actual reflect
 *   @l_actual: address to receive actual line
 *   @r_guess:  address to receive guess reflect
 *   @l_guess:  address to receive guess line
 */
static void make_random_parameters(double complex *r_actual,
	double complex *l_actual, double complex *r_guess,
	double complex *l_guess)
{
    double actual_magnitude;
    double actual_angle;
    double distance2;

    /*
     * Find actual reflect.  Magnitude is constrained to 0.25 .. 1.0;
     * angle is not constrained.
     */
    actual_magnitude = 0.25 + 0.75 * random() / RAND_MAX;
    actual_angle = 360.0 * random() / RAND_MAX - 180.0;
    *r_actual = actual_magnitude * cexp(I * M_PI / 180.0 * actual_angle);

    /*
     * Find the actual line.  Magnitude is constrained from 0.25 .. 1.0;
     * angle is constrained to 20..160 or 200..350 degrees to prevent
     * it from being * too close to through.
     */
    actual_magnitude = 0.25 + 0.75 * random() / RAND_MAX;
    actual_angle = 140.0 * (2.0 * random() / RAND_MAX - 1.0);
    if (actual_angle >= 0.0) {
	actual_angle += 20.0;
    } else {
	actual_angle -= 20.0;
    }
    *l_actual = actual_magnitude * cexp(I * M_PI / 180.0 * actual_angle);

    /*
     * There are four solutions to TRL:
     *     R,   L
     *    -R,   L
     *    1/R, 1/L
     *   -1/R, 1/L
     *
     * We need an initial guess that's always closer to the first solution
     * than to the others.  Stategy: find the minimum distance from
     * the first solution to each of the others in 4-space.  Choose a
     * gaussian random vector to serve as the uniform direction in
     * 4-space between the actual and guess.  Scale the vector to be
     * less than half of the distance calculated in the first step.
     * Add the actual answer to produce the guess values.
     */
    {
	double complex ri = 1.0 / *r_actual, li = 1.0 / *l_actual;
	double complex ctemp1, ctemp2;
	double dtemp;

	/*
	 * Find squared distance to -R, L
	 */
	ctemp1 = -*r_actual - *r_actual;
	distance2 = creal(ctemp1 * conj(ctemp1));

	/*
	 * Find squared distance to 1/R, 1/L
	 */
	ctemp1 = ri - *r_actual;
	ctemp2 = li - *l_actual;
	dtemp = creal(ctemp1 * conj(ctemp1)) + creal(ctemp2 * conj(ctemp2));
	if (dtemp < distance2) {
	    distance2 = dtemp;
	}

	/*
	 * Find squared distance to -1/R, 1/L
	 */
	ctemp1 = -ri - *r_actual;
	dtemp = creal(ctemp1 * conj(ctemp1)) + creal(ctemp2 * conj(ctemp2));
	if (dtemp < distance2) {
	    distance2 = dtemp;
	}
    }
    assert(distance2 > 0.0);

    /*
     * Choose gaussian random guess values and normalize to be less
     * than half of sqrt(distance2).
     */
    {
	double temp;

	*r_guess = libt_crandn();
	*l_guess = libt_crandn();
	temp = distance2 / (creal(*r_guess * conj(*r_guess)) +
	                    creal(*l_guess * conj(*l_guess)));
	temp = 0.45 * sqrt(temp);	/* maximum magnitude */
	temp *= random() / RAND_MAX;
	*r_guess = temp * *r_guess + *r_actual;
	*l_guess = temp * *l_guess + *l_actual;
    }
}

/*
 * run_vnacal_trl_trial: run a through-reflect-line calibration trial
 */
static libt_result_t run_vnacal_trl_trial(int trial, vnacal_type_t type)
{
    libt_result_t result = T_SKIPPED;
    vnacal_t *vcp = NULL;
    libt_vnacal_terms_t *ttp = NULL;
    vnacal_new_t *vnp = NULL;
    libt_vnacal_measurements_t *tmp = NULL;
    double complex r_actual[TRL_FREQUENCIES];	/* actual reflection */
    double complex r_guess[TRL_FREQUENCIES];	/* guess for reflection */
    double complex l_actual[TRL_FREQUENCIES];	/* actual line */
    double complex l_guess[TRL_FREQUENCIES];	/* guess for line */
    int r_unknown = -1, l_unknown = -1;

    /*
     * If -v, print the test header.
     */
    if (opt_v != 0) {
	(void)printf("Test vnacal TRL calibration trial %d type %-4s TRL\n",
		trial, _vnacal_type_to_name(type));
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
		    TRL_FREQUENCIES, /*frequency_vector*/NULL,
		    0.10, 0)) == NULL) {
	(void)fprintf(stderr, "%s: libt_vnacal_generate_error_terms: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    vnp = ttp->tt_vnp;

    /*
     * Generate random reflect and line parameters.
     */
    for (int findex = 0; findex < TRL_FREQUENCIES; ++findex) {
	make_random_parameters(&r_actual[findex], &l_actual[findex],
		               &r_guess[findex],  &l_guess[findex]);
    }
    if (opt_v > 1) {
	(void)printf("actual:\n");
	for (int findex = 0; findex < TRL_FREQUENCIES; ++findex) {
	    (void)printf("R %8.5f%+8.5f  L %8.5f%+8.5fj\n",
		    creal(r_actual[findex]), cimag(r_actual[findex]),
		    creal(l_actual[findex]), cimag(l_actual[findex]));
	    (void)printf("    %8.5f <%8.3f  %8.5f <%8.3f\n",
		    cabs(r_actual[findex]),
		    carg(r_actual[findex]) * 180.0 / M_PI,
		    cabs(l_actual[findex]),
		    carg(l_actual[findex]) * 180.0 / M_PI);
	}
	(void)printf("\n");
	(void)printf("guess:\n");
	for (int findex = 0; findex < TRL_FREQUENCIES; ++findex) {
	    (void)printf("R %8.5f%+8.5f  L %8.5f%+8.5fj\n",
		    creal(r_guess[findex]), cimag(r_guess[findex]),
		    creal(l_guess[findex]), cimag(l_guess[findex]));
	    (void)printf("    %8.5f <%8.3f %8.5f <%8.3f\n",
		    cabs(r_guess[findex]),
		    carg(r_guess[findex]) * 180.0 / M_PI,
		    cabs(l_guess[findex]),
		    carg(l_guess[findex]) * 180.0 / M_PI);
	}
	(void)printf("\n");
    }

    /*
     * Allocate the measurements matrices.
     */
    if ((tmp = libt_vnacal_alloc_measurements(type, /*m_rows*/2, /*m_columns*/2,
		    TRL_FREQUENCIES, /*ab*/false)) == NULL) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Add the through standard.
     */
    if (libt_vnacal_add_through(ttp, tmp, 1, 2) == -1) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Add the reflect standard.
     */
    {
	int p_actual = -1, p_guess = -1;
	int s_matrix[2][2];

	if ((p_actual = vnacal_make_vector_parameter(vcp, TRL_FREQUENCIES,
			ttp->tt_frequency_vector, r_actual)) == -1) {
	    result = T_FAIL;
	    goto out;
	}
	if ((p_guess = vnacal_make_vector_parameter(vcp, TRL_FREQUENCIES,
			ttp->tt_frequency_vector, r_guess)) == -1) {
	    result = T_FAIL;
	    goto out;
	}
	if ((r_unknown = vnacal_make_unknown_parameter(vcp, p_guess)) == -1) {
	    result = T_FAIL;
	    goto out;
	}
	s_matrix[0][0] = p_actual;
	s_matrix[0][1] = VNACAL_ZERO;
	s_matrix[1][0] = VNACAL_ZERO;
	s_matrix[1][1] = p_actual;
	if (libt_vnacal_calculate_measurements(ttp, tmp, &s_matrix[0][0],
		    2, 2, NULL) == -1) {
	    result = T_FAIL;
	    goto out;
	}
	if (vnacal_new_add_double_reflect_m(vnp,
		    tmp->tm_b_matrix, /*m_rows*/2, /*m_columns*/2,
		    r_unknown, r_unknown, 1, 2) == -1) {
	    result = T_FAIL;
	    goto out;
	}
	vnacal_delete_parameter(vcp, p_guess);
	vnacal_delete_parameter(vcp, p_actual);
    }

    /*
     * Add the line standard.
     */
    {
	int p_actual = -1, p_guess = -1;
	int s_matrix[2][2];

	if ((p_actual = vnacal_make_vector_parameter(vcp, TRL_FREQUENCIES,
			ttp->tt_frequency_vector, l_actual)) == -1) {
	    result = T_FAIL;
	    goto out;
	}
	if ((p_guess = vnacal_make_vector_parameter(vcp, TRL_FREQUENCIES,
			ttp->tt_frequency_vector, l_guess)) == -1) {
	    result = T_FAIL;
	    goto out;
	}
	if ((l_unknown = vnacal_make_unknown_parameter(vcp, p_guess)) == -1) {
	    result = T_FAIL;
	    goto out;
	}
	s_matrix[0][0] = VNACAL_MATCH;
	s_matrix[0][1] = p_actual;
	s_matrix[1][0] = p_actual;
	s_matrix[1][1] = VNACAL_MATCH;
	if (libt_vnacal_calculate_measurements(ttp, tmp, &s_matrix[0][0],
		    2, 2, NULL) == -1) {
	    result = T_FAIL;
	    goto out;
	}
	s_matrix[0][1] = l_unknown;
	s_matrix[1][0] = l_unknown;
	if (vnacal_new_add_line_m(vnp, tmp->tm_b_matrix, 2, 2,
		    &s_matrix[0][0], 1, 2) == -1) {
	    result = T_FAIL;
	    goto out;
	}
	vnacal_delete_parameter(vcp, p_guess);
	vnacal_delete_parameter(vcp, p_actual);
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
    for (int findex = 0; findex < TRL_FREQUENCIES; ++findex) {
	double frequency = ttp->tt_frequency_vector[findex];
	double complex r_solved, l_solved;

	if ((r_solved = vnacal_get_parameter_value(vcp, r_unknown,
			frequency)) == HUGE_VAL) {
	    result = T_FAIL;
	    goto out;
	}
	if ((l_solved = vnacal_get_parameter_value(vcp, l_unknown,
			frequency)) == HUGE_VAL) {
	    result = T_FAIL;
	    goto out;
	}
	if (opt_v > 1) {
	    (void)printf("findex %d:\n", findex);
	    (void)printf("  r_actual %e%+ej\n",
		    creal(r_actual[findex]), cimag(r_actual[findex]));
	    (void)printf("  r_solved %e%+ej\n",
		    creal(r_solved), cimag(r_solved));
	    (void)printf("  delta %e\n",
		    cabs(r_solved - r_actual[findex]));
	    (void)printf("  l_actual %e%+ej\n",
		    creal(l_actual[findex]), cimag(l_actual[findex]));
	    (void)printf("  l_solved %e%+ej\n",
		    creal(l_solved), cimag(l_solved));
	    (void)printf("  delta %e\n",
		    cabs(l_solved - l_actual[findex]));
	    (void)printf("\n");
	}
	if (!libt_isequal(r_solved, r_actual[findex])) {
	    if (opt_a) {
		assert(!"data miscompare");
	    }
	    result = T_FAIL;
	    goto out;
	}
	if (!libt_isequal(l_solved, l_actual[findex])) {
	    if (opt_a) {
		assert(!"data miscompare");
	    }
	    result = T_FAIL;
	    goto out;
	}
    }
    result = T_PASS;

out:
    if (r_unknown >= 0) {
	vnacal_delete_parameter(vcp, r_unknown);
    }
    if (l_unknown >= 0) {
	vnacal_delete_parameter(vcp, l_unknown);
    }
    libt_vnacal_free_measurements(tmp);
    libt_vnacal_free_error_terms(ttp);
    vnacal_free(vcp);
    return result;
}

/*
 * test_vnacal_trl: through-reflect-line calibration
 */
static libt_result_t test_vnacal_trl()
{
    libt_result_t result = T_SKIPPED;
    static vnacal_type_t type_array[] = {
	VNACAL_T8, VNACAL_U8, VNACAL_TE10, VNACAL_UE10
    };

    for (int trial = 0; trial < NTRIALS; ++trial) {
	for (int t_index = 0; t_index < 4; ++t_index) {
	    result = run_vnacal_trl_trial(trial, type_array[t_index]);
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
    libt_isequal_eps = 0.01;
    exit(test_vnacal_trl());
}
