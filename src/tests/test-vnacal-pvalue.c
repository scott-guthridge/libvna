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
#include "libt.h"
#include "libt_vnacal.h"


/*
 * Number of test trials to run
 */
#define NTRIALS		5

/*
 * N is the number of calibrations we solve on each trial and the
 * number of points in the empiricle CDF.  KS_THRESHOLD is the test
 * threshold for N=1000, p=0.001.  These constants must go together.
 *
 * To find the constants, first numerically invert kolmogorov_smirnov_cdf
 * in octave to find P(x < 1.9495) = 0.999.  The threshold is then:
 *
 *     1.9495 / sqrt(N)
 */
#define N		1000
#define KS_THRESHOLD	0.0616486

/*
 * Vector of pvalues, one per experiment.  Number of frequencies
 * is fixed at 1.
 */
static double pvalues[N];

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
 * Allow a small number of vnacal_new_solve calls to fail in each trial
 * due to random error.
 */
#define ALLOWED_SOLVE_FAILURES	3

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

/*
 * f_cmp: compare doubles for qsort
 */
static int f_cmp(const void *pv1, const void *pv2)
{
    double v1 = *(const double *)pv1;
    double v2 = *(const double *)pv2;

    if (v1 < v2) {
	return -1;
    }
    if (v1 > v2) {
	return +1;
    }
    return 0;
}


/*
 * run_solt_trial_helper: add short, open and load calibrations on port
 *   @ttp: pointer to test error terms structure
 *   @tmp: test measurements structure
 *   @port: port to calibrate
 */
static libt_result_t run_solt_trial_helper(libt_vnacal_terms_t *ttp,
	libt_vnacal_measurements_t *tmp, int port)
{
    if (libt_vnacal_add_single_reflect(ttp, tmp, VNACAL_SHORT, port) == -1) {
	return T_FAIL;
    }
    if (libt_vnacal_add_single_reflect(ttp, tmp, VNACAL_OPEN, port) == -1) {
	return T_FAIL;
    }
    if (libt_vnacal_add_single_reflect(ttp, tmp, VNACAL_MATCH, port) == -1) {
	return T_FAIL;
    }

    /*
     * If UE14 or E12, one row or one column, add another standard to keep
     * the system overdetermined.
     */
    if (ttp->tt_layout.vl_type == VNACAL_UE14 ||
	    ttp->tt_layout.vl_type == VNACAL_E12 ||
	    VL_M_ROWS(&ttp->tt_layout) == 1 ||
	    VL_M_COLUMNS(&ttp->tt_layout) == 1) {
	vnacal_t *vcp = ttp->tt_vnp->vn_vcp;
	int p;

	if ((p = vnacal_make_scalar_parameter(vcp, I)) == -1) {
	    return T_FAIL;
	}
	if (libt_vnacal_add_single_reflect(ttp, tmp, p, port) == -1) {
	    return T_FAIL;
	}
	if (vnacal_delete_parameter(vcp, p) == -1) {
	    return T_FAIL;
	}
    }
    return T_PASS;
}

/*
 * run_one_experiment: test if pvalue distribution is linear
 *   @experiment: experiment count (starting on zero)
 *   @type: error term type
 *   @m_rows: number of VNA ports that detect signal
 *   @m_columns: number of VNA ports that generate signal
 */
static libt_result_t run_one_experiment(int experiment,
	vnacal_type_t type, int m_rows, int m_columns)
{
    const int ports     = MAX(m_rows, m_columns);
    const int diagonals = MIN(m_rows, m_columns);
    vnacal_t *vcp = NULL;
    libt_vnacal_terms_t *ttp = NULL;
    vnacal_new_t *vnp = NULL;
    libt_vnacal_measurements_t *tmp = NULL;
    const double sigma_fl = 1.0e-3;
    libt_result_t result = T_FAIL;

    /*
     * If -vv, print a header.
     */
    if (opt_v >= 2) {
	(void)printf("experiment %3d size %d x %d "
		"type %-4s\n",
		experiment, m_rows, m_columns, vnacal_type_to_name(type));
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
    if ((ttp = libt_vnacal_generate_error_terms(vcp, type, m_rows, m_columns,
		    /*frequencies*/1, NULL, 0)) == NULL) {
	(void)fprintf(stderr, "%s: libt_vnacal_generate_error_terms: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    vnp = ttp->tt_vnp;

    /*
     * Allocate the measurements matrices.
     */
    if ((tmp = libt_vnacal_alloc_measurements(type, m_rows, m_columns,
		    /*frequencies*/1, /*ab*/false)) == NULL) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Because we're generating many experiments and we expect the pvalue
     * to be uniformly distributed between 0 and 1, it's not unlikely
     * that we'll hit a low value along the way and fail the pvalue test.
     * Set the threshold very low to avoid false positives: we're testing
     * the distribution; not than that all trials succeed.
     */
    if (vnacal_new_set_pvalue_limit(vnp, 1.0e-8) == -1) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Set the measurement error.
     */
    libt_vnacal_sigma_n = &sigma_fl;
    if (vnacal_new_set_m_error(vnp, NULL, /*frequencies*/1,
		&sigma_fl, /*sigma_tr*/NULL) == -1) {
	result = T_FAIL;
	goto out;
    }

    /*
     * For two-port T16 and U16, use T-MM-SS-SM-MS.
     */
    if ((type == VNACAL_T16 || type == VNACAL_U16) && ports > 1) {
	assert(ports == 2);	/* for > 2, need to use mapped matrix */

	if (libt_vnacal_add_through(ttp, tmp, 1, 2) == -1) {
	    result = T_FAIL;
	    goto out;
	}
	if (libt_vnacal_add_double_reflect(ttp, tmp,
		    VNACAL_MATCH, VNACAL_MATCH, 1, 2) == -1) {
	    result = T_FAIL;
	    goto out;
	}
	if (libt_vnacal_add_double_reflect(ttp, tmp,
		    VNACAL_SHORT, VNACAL_SHORT, 1, 2) == -1) {
	    result = T_FAIL;
	    goto out;
	}
	if (libt_vnacal_add_double_reflect(ttp, tmp,
		    VNACAL_SHORT, VNACAL_MATCH, 1, 2) == -1) {
	    result = T_FAIL;
	    goto out;
	}
	if (libt_vnacal_add_double_reflect(ttp, tmp,
		    VNACAL_MATCH, VNACAL_SHORT, 1, 2) == -1) {
	    result = T_FAIL;
	    goto out;
	}

    /*
     * Otherwise, use SOLT.
     */
    } else {
	/*
	 * Run short, open, and load test on every diagonal port.
	 */
	for (int port = 1; port <= diagonals; ++port) {
	    result = run_solt_trial_helper(ttp, tmp, port);
	    if (result != T_PASS)
		goto out;
	}

	/*
	 * Run through tests between every diagonal port and every other port.
	 */
	for (int port1 = 1; port1 <= diagonals; ++port1) {
	    for (int port2 = port1 + 1; port2 <= ports; ++port2) {
		if (libt_vnacal_add_through(ttp, tmp, port1, port2) == -1) {
		    result = T_FAIL;
		    goto out;
		}
	    }
	}
    }
    libt_vnacal_free_measurements(tmp);
    tmp = NULL;

    /*
     * Use hidden API to receive pvalue back from vnacal_new_solve.
     */
    assert(experiment < N);
    vnp->vn_pvalue_vector = &pvalues[experiment];

    /*
     * Solve for the error parameters and check.
     */
    if (vnacal_new_solve(ttp->tt_vnp) == -1) {
	(void)printf("%s: vnacal_solve: %s\n", progname, strerror(errno));
	result = T_SKIPPED;
	goto out;
    }
    /*
     * We skip this because there's a non-neglible chance that out our
     * measurement errors cause it to occasionally fail.  This isn't
     * what we're testing.
     */
#if 0
    if (libt_vnacal_validate_calibration(ttp, NULL) == -1) {
	goto out;
    }
#endif
    result = T_PASS;

out:
    libt_vnacal_free_error_terms(ttp);
    vnacal_free(vcp);
    return result;
}

/*
 * run_trial: run a test trial
 *   @trial: trial number (starting on 1)
 *   @type: error term type
 *   @m_rows: number of VNA ports that detect signal
 *   @m_columns: number of VNA ports that generate signal
 */
static libt_result_t run_trial(int trial,
	vnacal_type_t type, int m_rows, int m_columns)
{
    double max_deviation = 0.0;
    int n_solve_failures = 0;
    libt_result_t result = T_FAIL;

    /*
     * If -v, print the test header.
     */
    if (opt_v != 0) {
	(void)printf("Test vnacal-new-pvalue: trial %3d size %d x %d "
		"type %-4s\n",
		trial, m_rows, m_columns, vnacal_type_to_name(type));
    }

    /*
     * Run N experiments.  Sort the resulting pvalues to create the
     * empiricle CDF.
     */
    for (int experiment = 0; experiment < N; ++experiment) {
	result = run_one_experiment(experiment,
		type, m_rows, m_columns);
	if (result == T_SKIPPED) {
	    if (++n_solve_failures > ALLOWED_SOLVE_FAILURES) {
		result = T_FAIL;
		goto out;
	    }
	    --experiment;
	    continue;
	}
	if (result != T_PASS)
	    goto out;
    }
    qsort((void *)pvalues, N, sizeof(double), f_cmp);

    /*
     * Find maximum deviation from a uniform distribution and
     * apply KS test for N points at the chosen confidence.
     */
    if (opt_v >= 2) {
	(void)printf("pvalues = [\n");
    }
    for (int i = 0; i < N; ++i) {
	double dy = fabs(pvalues[i] - (double)i / (N - 1));

	if (opt_v >= 2) {
	    (void)printf("    %f\n", pvalues[i]);
	}
	if (dy > max_deviation) {
	    max_deviation = dy;
	}
    }
    if (opt_v >= 2) {
	(void)printf("];\n");
    }
    if (opt_v > 0) {
	(void)printf("%% max_deviation: %f\n", max_deviation);
    }
    if (max_deviation > KS_THRESHOLD) {
	(void)printf("max_deviation %f failed KS test\n", max_deviation);
	result = T_FAIL;
	goto out;
    }
    result = T_PASS;

out:
    return result;
}

/*
 * test_vnacal_new_pvalue: test vnacal_new_* with random multi-port standards
 */
static libt_result_t test_vnacal_new_pvalue()
{
    /* Note: we're testing only 2x2 mainly because the test runs too
       long when we include the other sizes.  Also, for non 2x2 T16/U16,
       we need appropriate standards above. */
    //static const int sizes[] = { 1, 2, 3, 4 };
    static const int sizes[] = { 2 };
    static const vnacal_type_t types[] = {
	VNACAL_T8, VNACAL_U8, VNACAL_TE10, VNACAL_UE10,
	VNACAL_T16, VNACAL_U16, VNACAL_UE14, VNACAL_E12
    };
    libt_result_t result = T_FAIL;

    /*
     * For each trial...
     */
    for (int trial = 1; trial <= NTRIALS; ++trial) {
	for (int si = 0; si < sizeof(sizes) / sizeof(int); ++si) {
	    for (int sj = 0; sj < sizeof(sizes) / sizeof(int); ++sj) {
		int rows = sizes[si];
		int columns = sizes[sj];

		for (int ti = 0; ti < sizeof(types) / sizeof(types[0]); ++ti) {
		    vnacal_type_t type = types[ti];

		    if (type == VNACAL_T8 || type == VNACAL_TE10 ||
			    type == VNACAL_T16) {
			if (rows > columns) {
			    continue;
			}
		    } else {
			if (rows < columns) {
			    continue;
			}
		    }
		    result = run_trial(trial, type, rows, columns);
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
    libt_isequal_eps = 0.1;	/* we're not testing this: be tolerant */
    exit(test_vnacal_new_pvalue());
}
