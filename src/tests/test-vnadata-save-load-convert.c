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
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include "vnadata.h"
#include "libt.h"
#include "libt_vnadata.h"


#define N_TRIALS	4

/*
 * Options
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
int opt_v = 0;

/*
 * error_fn: error reporting function
 *   @message: error message
 *   @arg: (unused)
 *   @category: error category (unused)
 */
static void error_fn(const char *message, void *arg, vnaerr_category_t category)
{
    (void)printf("error: %s: %s\n", progname, message);
}

#define FREQUENCIES	5


/*
 * run_trial
 */
static int run_trial(int trial, vnadata_filetype_t filetype,
	vnadata_parameter_type_t type, int rows, int columns,
	const char *format, vnadata_parameter_type_t load_type,
	libt_vnadata_z0_type_t z0_type)
{
    libt_vnadata_t *tdp = NULL;
    vnadata_t *vdp = NULL;
    int load_rows, load_columns;
    const char *filename;
    FILE *fp = NULL;
    libt_result_t result = T_FAIL;

    /*
     * Get the test filename.
     */
    switch (filetype) {
    case VNADATA_FILETYPE_TOUCHSTONE1:
	switch (rows) {
	case 1:
	    filename = "test-vnadata.s1p";
	    break;
	case 2:
	    filename = "test-vnadata.s2p";
	    break;
	case 3:
	    filename = "test-vnadata.s3p";
	    break;
	case 4:
	    filename = "test-vnadata.s4p";
	    break;
	default:
	    assert(!"unhandled case in switch");
	}
	break;
    case VNADATA_FILETYPE_TOUCHSTONE2:
	filename = "test-vnadata.ts";
	break;
    case VNADATA_FILETYPE_NPD:
	filename = "test-vnadata.npd";
	break;
    default:
	assert(!"unhandled case in switch");
    }

    /*
     * If verbose, report the test case.
     */
    if (opt_v >= 1) {
	const char *type_name = vnadata_get_type_name(type);

	(void)printf("Test SLC: trial %2d type %-3s size %d x %d "
		"%s %s %s\n",
		trial, (type == VPT_UNDEF) ? "-" : type_name,
		rows, columns, filename, format,
		libt_vnadata_z0_names[z0_type]);
	(void)fflush(stdout);
    }

    /*
     * Create test values.  Note that libt_vnadata_create simply fills
     * the matrix with small random numbers centered around zero.  For S
     * parameters, this is OK, but for parameter types with impedances
     * or addmittances, using small numbers far from the Z0 leads to
     * poorly conditioned matrices.  Fix this by pretending the parameters
     * are in S and converting them to the desired type.
     */
    if ((tdp = libt_vnadata_create(type, rows, columns, FREQUENCIES,
		    z0_type)) == NULL) {
	libt_error("libt_vnadata_create: %s\n", strerror(errno));
	/*NOTREACHED*/
    }
    if (type != VPT_S && type != VPT_ZIN) {
	const double complex *z0_vector = NULL;

	if (z0_type != Z0_PER_F) {
	    z0_vector = tdp->td_z0_vector;
	}
	for (int findex = 0; findex < FREQUENCIES; ++findex) {
	    if (z0_type == Z0_PER_F) {
		z0_vector = tdp->td_fz0_vector[findex];
	    }
	    libt_vnadata_convert(tdp->td_vector[findex],
		    tdp->td_vector[findex], z0_vector, tdp->td_rows,
		    tdp->td_columns, VPT_S, type);
	}
	if (opt_v >= 2) {
	    (void)printf("After conversion to %s:\n",
		    vnadata_get_type_name(type));
	    for (int findex = 0; findex < FREQUENCIES; ++findex) {
		for (int row = 0; row < rows; ++row) {
		    for (int column = 0; column < columns; ++column) {
			int cell = row * columns + column;
			double complex value = tdp->td_vector[findex][cell];

			(void)printf("  %9.6f%+9.6fj",
				creal(value), cimag(value));
		    }
		    (void)printf("\n");
		}
		(void)printf("\n");
	    }
	}
    }

    /*
     * Allocate the vnadata_t structure, fill it from the test values
     * and verify that they match.
     */
    if ((vdp = vnadata_alloc(error_fn, NULL)) == NULL) {
	libt_fail("vnadata_alloc: returned NULL\n");
	result = T_FAIL;
	goto out;
    }
    if ((result = libt_vnadata_fill(tdp, vdp, FM_MATRIX)) != T_PASS) {
	goto out;
    }
    if ((result = libt_vnadata_validate(tdp, vdp)) != T_PASS) {
	goto out;
    }

#ifdef TEST_FULL_PRECISION
    /*
     * Set both frequency and data precision to maximum, which
     * uses hexadecimal floating point notation to avoid losing
     * any precision in save and load.
     *
     * This is ifdefed out by default to avoid test failures on
     * architectures that don't support hexadecimal floating point.
     */
    if (vnadata_set_fprecision(vdp, VNADATA_MAX_PRECISION) == -1) {
	libt_fail("vnadata_set_fprecision: returned -1\n");
	result = T_FAIL;
	goto out;
    }
    if (vnadata_set_dprecision(vdp, VNADATA_MAX_PRECISION) == -1) {
	libt_fail("vnadata_set_dprecision: returned -1\n");
	result = T_FAIL;
	goto out;
    }
#endif

    /*
     * Set type file format.
     */
    if (vnadata_set_format(vdp, format) == -1) {
	libt_fail("vnadata_set_format: returned -1\n");
	result = T_FAIL;
	goto out;
    }

    /*
     * If trial is one less than doulby even, run vnadata_cksave().
     */
    if (((trial + 1) & 2) == 0) {
	if (vnadata_cksave(vdp, filename) == -1) {
	    libt_fail("vnadata_cksave: returned -1\n");
	    result = T_FAIL;
	    goto out;
	}
    }

    /*
     * Save the parameters to a file.  If trial is even, use save;
     * otherwise use fsave.
     */
    if ((trial & 1) == 0) {
	if (vnadata_save(vdp, filename) == -1) {
	    libt_fail("vnadata_save: returned -1\n");
	    result = T_FAIL;
	    goto out;
	}
    } else {
	if (vnadata_set_filetype(vdp, filetype) == -1) {
	    libt_fail("vnadata_set_filetype: returned -1\n");
	    result = T_FAIL;
	    goto out;
	}
	if ((fp = fopen(filename, "w")) == NULL) {
	    libt_error("fopen: %s: %s\n", filename, strerror(errno));
	    /*NOTREACHED*/
	}
	if (vnadata_fsave(vdp, fp, "-") == -1) {
	    libt_fail("vnadata_fsave: returned -1\n");
	    result = T_FAIL;
	    goto out;
	}
	if (fclose(fp) == -1) {
	    libt_error("fclose: %s: %s\n", filename, strerror(errno));
	}
	fp = NULL;
    }
    vnadata_free(vdp);
    vdp = NULL;

    /*
     * A few formats, IL, RL and VSWR don't give enough information to
     * reconstruct any parameter type, so are not loadable.  If we saved
     * one of those, this is as far as we can go.  That we can't load
     * those formats back leaves a small gap in the testing in that we
     * can't verify that we saved them correctly.  Another test could
     * be written to parse the save file and check.
     */
    if (load_type == VPT_UNDEF) {
	result = T_PASS;
	goto out;
    }

    /*
     * Create a new vnadata_t structure and load from the file.  If trial
     * is doubly even, use load; otherwise use use fload.
     */
    if ((vdp = vnadata_alloc(error_fn, NULL)) == NULL) {
	libt_fail("vnadata_alloc: returned NULL\n");
	result = T_FAIL;
	goto out;
    }
    if ((trial & 2) == 0) {
	if (vnadata_load(vdp, filename) == -1) {
	    libt_fail("vnadata_load: returned -1\n");
	    result = T_FAIL;
	    goto out;
	}
    } else {
	if (vnadata_set_filetype(vdp, filetype) == -1) {
	    libt_fail("vnadata_set_filetype: returned -1\n");
	    result = T_FAIL;
	    goto out;
	}
	if ((fp = fopen(filename, "r")) == NULL) {
	    libt_error("fopen: %s: %s\n", filename, strerror(errno));
	    /*NOTREACHED*/
	}
	if (vnadata_fload(vdp, fp, "-") == -1) {
	    libt_fail("vnadata_fload: returned -1\n");
	    result = T_FAIL;
	    goto out;
	}
	if (fclose(fp) == -1) {
	    libt_error("fclose: %s: %s\n", filename, strerror(errno));
	}
	fp = NULL;
    }

    /*
     * Test that the loaded parameter type is what we expected.
     */
    if (vnadata_get_type(vdp) != load_type) {
	libt_fail("expected load to return type %s but found type %s\n",
		vnadata_get_type_name(load_type),
		vnadata_get_type_name(vnadata_get_type(vdp)));
	result = T_FAIL;
	goto out;
    }

    /*
     * Check the loaded dimensions.
     */
    if (type != VPT_ZIN && load_type == VPT_ZIN) {
	load_rows    = 1;
	load_columns = MIN(tdp->td_rows, tdp->td_columns);
    } else {
	load_rows    = tdp->td_rows;
	load_columns = tdp->td_columns;
    }
    if (vnadata_get_rows(vdp) != load_rows) {
	libt_fail("expected %d columns from load; found %d\n",
		rows, vnadata_get_rows(vdp));
	result = T_FAIL;
	goto out;
    }
    if (vnadata_get_columns(vdp) != load_columns) {
	libt_fail("expected %d columns from load; found %d\n",
		columns, vnadata_get_columns(vdp));
	result = T_FAIL;
	goto out;
    }
    if (vnadata_get_frequencies(vdp) != FREQUENCIES) {
	libt_fail("expected %d frequencies from load; found %d\n",
		FREQUENCIES, vnadata_get_frequencies(vdp));
	result = T_FAIL;
	goto out;
    }

    /*
     * Check the loaded frequencies.
     */
    {
	const double *frequency_vector = vnadata_get_frequency_vector(vdp);

	for (int findex = 0; findex < FREQUENCIES; ++findex) {
	    if (!libt_isequal_d_rpt("frequency",
			frequency_vector[findex],
			tdp->td_frequency_vector[findex])) {
		libt_fail(": findex %d\n", findex);
		result = T_FAIL;
		goto out;
	    }
	}
    }

    /*
     * Check the loaded z0's.
     */
    {
	int ports = MAX(tdp->td_rows, tdp->td_columns);

	/*
	 * There's a small oddity here where loading Zin can discard
	 * some z0 entries.  Example: we start with a 2x1 S parameter
	 * matrix and save in Zinri format.  In the save file, we have
	 * z0 entries for both ports.  But when we load the file back,
	 * we have only the z0 for the port that had an input impedance.
	 * Fix the port count for this case.
	 */
	if (load_type == VPT_ZIN) {
	    ports = load_columns;
	}
	if (tdp->td_z0_type != Z0_PER_F) {
	    const double complex *z0_vector = vnadata_get_z0_vector(vdp);

	    for (int port = 0; port < ports; ++port) {
		if (!libt_isequal_c_rpt("z0_vector",
			    z0_vector[port],
			    tdp->td_z0_vector[port])) {
		    libt_fail(": port %d\n", port);
		    result = T_FAIL;
		    goto out;
		}
	    }
	} else {
	    for (int findex = 0; findex < FREQUENCIES; ++findex) {
		const double complex *z0_vector;

		z0_vector = vnadata_get_fz0_vector(vdp, findex);
		for (int port = 0; port < ports; ++port) {
		    if (!libt_isequal_c_rpt("fz0_vector",
				z0_vector[port],
				tdp->td_fz0_vector[findex][port])) {
			libt_fail(": port %d\n", port);
			result = T_FAIL;
			goto out;
		    }
		}
	    }
	}
    }

    /*
     * Check the loaded data.
     */
    for (int findex = 0; findex < FREQUENCIES; ++findex) {
	const double complex *actual = vnadata_get_matrix(vdp, findex);
	const double complex *z0_vector = vnadata_get_fz0_vector(vdp, findex);
	double complex expected[load_rows][load_columns];

	(void)memset((void *)expected, 0, sizeof(expected));
	libt_vnadata_convert(tdp->td_vector[findex], &expected[0][0],
	    z0_vector, tdp->td_rows, tdp->td_columns, tdp->td_type, load_type);
	if (opt_v >= 2) {
	    (void)printf("findex %d\n", findex);
	    (void)printf("expected load values:\n");
	    for (int row = 0; row < load_rows; ++row) {
		for (int column = 0; column < load_columns; ++column) {
		    (void)printf(" %+f%+fj",
			    creal(expected[row][column]),
			    cimag(expected[row][column]));
		}
		(void)printf("\n");
	    }
	    (void)printf("\n");
	    (void)printf("actual load values:\n");
	    for (int row = 0; row < load_rows; ++row) {
		for (int column = 0; column < load_columns; ++column) {
		    int cell = load_columns * row + column;

		    (void)printf(" %+f%+fj",
			    creal(actual[cell]), cimag(actual[cell]));
		}
		(void)printf("\n");
	    }
	    (void)printf("\n");
	}
	for (int row = 0; row < load_rows; ++row) {
	    for (int column = 0; column < load_columns; ++column) {
		int cell = load_columns * row + column;

		if (!libt_isequal_c_rpt("data", actual[cell],
			    expected[row][column])) {
		    libt_fail(": findex %d row %d columns %d\n",
			    findex, row, column);
		    result = T_FAIL;
		    goto out;
		}
	    }
	}
    }

    /*
     * If possible, convert back to the original type and validate.
     */
    if (vnadata_get_type(vdp) != VPT_ZIN || tdp->td_type == VPT_ZIN) {
	if (vnadata_convert(vdp, vdp, tdp->td_type) == -1) {
	    libt_fail("vnadata_convert: returned -1\n");
	    result = T_FAIL;
	    goto out;
	}
	if ((result = libt_vnadata_validate(tdp, vdp)) != T_PASS) {
	    goto out;
	}
    }
    result = T_PASS;

out:
    if (fp != NULL) {
	(void)fclose(fp);
    }
    vnadata_free(vdp);
    libt_vnadata_free(tdp);
    return result;
}

/*
 * file_format_t: describes a file format
 */
typedef struct file_format {
    const char		       *ff_name;	/* format */
    vnadata_parameter_type_t	ff_type;	/* loaded type */
    bool			ff_touchstone;	/* touchstone compatible */
    bool			ff_2x2_only;	/* 2x2 matrix only */
} file_format_t;

/*
 * Table of Save Formats
 */
static const file_format_t format_table[] = {
    /*
     * Touchstone-compatible formats
     */
    { "SRI",		VPT_S,		true,	false },
    { "SMA",		VPT_S,		true,	false },
    { "SDB",		VPT_S,		true,	false },
    { "ZRI",		VPT_Z,		true,	true  },
    { "ZMA",		VPT_Z,		true,	true  },
    { "YRI",		VPT_Y,		true,	true  },
    { "YMA",		VPT_Y,		true,	true  },
    { "HRI",		VPT_H,		true,	true  },
    { "HMA",		VPT_H,		true,	true  },
    { "GRI",		VPT_G,		true,	true  },
    { "GMA",		VPT_G,		true,	true  },

    /*
     * Other matrix types (NPD only)
     */
    { "TRI",		VPT_T,		false,	true  },
    { "TMA",		VPT_T,		false,	true  },
    { "TDB",		VPT_T,		false,	true  },
    { "URI",		VPT_U,		false,	true  },
    { "UMA",		VPT_U,		false,	true  },
    { "UDB",		VPT_U,		false,	true  },
    { "ARI",		VPT_A,		false,  true  },
    { "AMA",		VPT_A,		false,  true  },
    { "BRI",		VPT_B,		false,  true  },
    { "BMA",		VPT_B,		false,  true  },

    /*
     * Multiple types (NPD only)
     */
    { "SDB,ARI,ZMA",	VPT_A,		false,	true  },
    { "ZRI,YMA",	VPT_Z,		false,	true  },
    { "IL,VSWR,SRC",	VPT_ZIN,	false,	true  },

    /*
     * RL, IL & VSWR: saveable but not loadable (NPD only)
     */
    { "IL",		VPT_UNDEF,	false,  false },
    { "RL",		VPT_UNDEF,	false,  false },
    { "VSWR",		VPT_UNDEF,	false,  false },

    /*
     * Input impedance types (NPD only)
     */
    { "ZINRI",		VPT_ZIN,	false,	false },
    { "ZINMA",		VPT_ZIN,	false,	false },
    { "PRC",		VPT_ZIN,	false,	false },
    { "PRL",		VPT_ZIN,	false,	false },
    { "SRC",		VPT_ZIN,	false,	false },
    { "SRL",		VPT_ZIN,	false,	false },

    { NULL,		VPT_UNDEF,	false,	false },
};


/*
 * test_vnadata_slc_helper: choose format
 */
static libt_result_t test_vnadata_slc_helper(int trial,
	vnadata_filetype_t filetype, vnadata_parameter_type_t type,
	int rows, int columns)
{
    libt_result_t result = T_SKIPPED;
    const file_format_t *ffp;

    for (ffp = format_table; ffp->ff_name != NULL; ++ffp) {
	/*
	 * When we reach the end of the Touchstone formats with a
	 * Touchstone filetype, skip the rest.
	 */
	if (filetype != VNADATA_FILETYPE_NPD && !ffp->ff_touchstone) {
	    break;
	}
	/*
	 * If the input matrix type is Zin, then we can only save
	 * in formats that load as Zin.
	 */
	if (type == VPT_ZIN && ffp->ff_type != VPT_ZIN) {
	    continue;
	}
	/*
	 * Return loss requires at least one off-diagonal element.
	 */
	if (strcasecmp(ffp->ff_name, "RL") == 0 && columns < 2) {
	    continue;
	}
	/*
	 * If we're given a rectangular S matrix, we can save only
	 * in S formats.
	 */
	if (type == VPT_S && rows != columns && ffp->ff_name[0] != 'S') {
	    continue;
	}
	/*
	 * If we're given a matrix of dimension other than 2x2, make
	 * sure we don't try to convert to a 2x2 format.
	 */
	if ((rows != 2 || columns != 2) && ffp->ff_2x2_only) {
	    continue;
	}
	/*
	 * Run with a single Z0 for all ports.
	 */
	result = run_trial(trial, filetype, type, rows, columns,
		ffp->ff_name, ffp->ff_type, Z0_SINGLE);
	if (result != T_PASS) {
	    return result;
	}
	if (filetype == VNADATA_FILETYPE_TOUCHSTONE1) {
	    continue;
	}
	/*
	 * If not Touchstone 1, run with a vector of real Z0's,
	 * different per port.
	 */
	result = run_trial(trial, filetype, type, rows, columns,
		ffp->ff_name, ffp->ff_type, Z0_REAL_VECTOR);
	if (result != T_PASS) {
	    return result;
	}
	if (filetype == VNADATA_FILETYPE_TOUCHSTONE2) {
	    continue;
	}
	/*
	 * If not Touchstone, run with complex system impedances
	 * and per-frequency complex system impedances.
	 */
	result = run_trial(trial, filetype, type, rows, columns,
		ffp->ff_name, ffp->ff_type, Z0_COMPLEX_VECTOR);
	if (result != T_PASS) {
	    return result;
	}
	result = run_trial(trial, filetype, type, rows, columns,
		ffp->ff_name, ffp->ff_type, Z0_PER_F);
	if (result != T_PASS) {
	    return result;
	}
    }
    return T_PASS;
}

/*
 * test_vnadata_slc: run save/load/convert tests on vnadata
 */
static libt_result_t test_vnadata_slc()
{
    libt_result_t result = T_SKIPPED;
    vnadata_parameter_type_t type;
    vnadata_filetype_t filetype;

    for (int trial = 0; trial < N_TRIALS; ++trial) {
	for (filetype = VNADATA_FILETYPE_TOUCHSTONE1;
		filetype <= VNADATA_FILETYPE_NPD; ++filetype) {
	    for (type = 0; type < VPT_NTYPES; ++type) {
		switch (type) {
		case VPT_UNDEF:
		    continue;

		case VPT_S:
		case VPT_Z:
		case VPT_Y:
		    {
			int max_ports;

			if (filetype == VNADATA_FILETYPE_TOUCHSTONE1) {
			    max_ports = 4;
			} else {
			    max_ports = 7;
			}
			for (int ports = 1; ports <= max_ports; ++ports) {
			    result = test_vnadata_slc_helper(trial,
				    filetype, type, ports, ports);
			    if (result != T_PASS) {
				goto out;
			    }
			}
		    }
		    break;

		case VPT_T:
		case VPT_U:
		case VPT_H:
		case VPT_G:
		case VPT_A:
		case VPT_B:
		    result = test_vnadata_slc_helper(trial, filetype,
			    type, 2, 2);
		    if (result != T_PASS) {
			goto out;
		    }
		    break;

		case VPT_ZIN:
		    if (filetype != VNADATA_FILETYPE_NPD) {
			continue;
		    }
		    for (int ports = 1; ports <= 7; ++ports) {
			result = test_vnadata_slc_helper(trial, filetype,
				type, 1, ports);
			if (result != T_PASS) {
			    goto out;
			}
		    }
		    break;

		default:
		    assert(!"unhandled case in switch");
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
    exit(2);
}

/*
 * main: test program
 */
int main(int argc, char **argv)
{
    if ((char *)NULL == (progname = strrchr(argv[0], '/'))) {
	progname = argv[0];
    } else {
	++progname;
    }
    for (;;) {
	switch (getopt(argc, argv, options)) {
	case 'a':
	    opt_a = true;
	    continue;

	case 'v':
	    ++opt_v;
	    continue;

	case -1:
	    break;

	default:
	    print_usage();
	}
	break;
    }
    argc -= optind;
    argv += optind;
    if (argc != 0) {
	print_usage();
    }

    /*
     * Set the error limit for numeric comparisons.  When not using
     * hexadecimal floating point, the numeric error due to using only
     * 6 digits of precision in the save files accumulates pretty high
     * in the 5x5 case, so we have to be a little lenient here.
     */
#ifdef TEST_FULL_PRECISION
    libt_isequal_init();
#else
    libt_isequal_eps = 0.1;
#endif

    exit(test_vnadata_slc());
}
