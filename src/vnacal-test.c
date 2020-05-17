/*
 * Vector Network Analyzer Calibration Library
 * Copyright © 2020 D Scott Guthridge <scott_guthridge@rompromity.net>
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
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include "vnaconv_internal.h"
#include "vnacal_internal.h"


char *progname;

#define PI	3.1415926535897932384626433832795
#define EPS	1.0e-4

#define NTRIALS		100


/*
 * Command Line Options
 */
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
static bool opt_a = false;
static bool opt_v = false;

/*
 * crandn: generate a random complex number where real and imaginary parts
 *	are normally distributed with zero mean and unit standard deviation
 */
static double complex crandn()
{
    double u1 = (random() + 1.0) / RAND_MAX;
    double u2 = (double)random() / RAND_MAX;
    double r = sqrt(-2.0 * log(u1));
    double a = 2 * PI * u2;

    return r * (cos(a) + I * sin(a));
}

/*
 * isequal: test if x and y are approximately equal
 */
static int isequal(double complex x, double complex y)
{
    int rv;
    double d = cabs(csqrt(x * y));

    if (d < 1.0) {
	d = 1.0;
    }
    rv = cabs(x - y) / d < EPS;
    if (!rv) {
	printf("|x-y| = %f\n", cabs(x - y));
	printf("%f%+fi != %f%+fi\n", creal(x), cimag(x), creal(y), cimag(y));
    }
    return rv;
}

/*
 * cmatrix_print: print an m by n serialized complex matrix
 */
static void cmatrix_print(double complex *a, int m, int n)
{
    for (int i = 0; i < m; ++i) {
	for (int j = 0; j < n; ++j) {
	    double complex v = a[i * n + j];

	    (void)printf(" %8.5f%+8.5fj", creal(v), cimag(v));
	}
	(void)printf("\n");
    }
    (void)printf("\n");
}

/*
 * test_result_type
 */
typedef enum test_result {
    T_PASS,
    T_FAIL,
    T_SKIPPED
} test_result_type;

/*
 * test counters
 */
static int test_count = 0;
static int fail_count = 0;

/*
 * report_test_result: report a test result
 */
static void report_test_result(const char *test_name, test_result_type result)
{
    const char *result_name;

    switch (result) {
    case T_PASS:
	result_name = "PASS";
	break;
    case T_FAIL:
	result_name = "FAIL";
	break;
    case T_SKIPPED:
	result_name = "SKIPPED";
	break;
    default:
	abort();
	/*NOTREACHED*/
    }
    (void)printf("Test %2d: %-58s %s\n", ++test_count, test_name, result_name);
    (void)fflush(stdout);
    if (result == T_FAIL) {
	++fail_count;
    }
}

static const char *error_term_names[] = {
    "e00", "e10e01", "e11"
};

/*
 * gen_error_terms: fill in the vnacal_calset_t with calibration
 *	values and return the error term matrix
 *   @vcsp: pointer to newly allocated vnacal_calset_t
 */
static double complex *(*gen_error_terms(vnacal_calset_t *vcsp))[3]
{
    int rows = vcsp->vcs_rows;
    int columns = vcsp->vcs_columns;
    int frequencies = vcsp->vcs_frequencies;
    double frequency_vector[frequencies];
    double complex references[3][frequencies];
    double complex cdata[rows][columns][3][frequencies];
    double complex *(*error_terms)[3];
    int ndiagonal = (rows <= columns) ? rows : columns;

    /*
     * Allocate the error terms matrix and contained frequency vectors.
     */
    if ((error_terms = (double complex *(*)[3])calloc(rows * columns,
		    sizeof(double complex *[3]))) == NULL) {
	(void)fprintf(stderr, "%s: calloc: %s\n", progname, strerror(errno));
	return NULL;
    }
    for (int i = 0; i < rows * columns; ++i) {
	for (int k = 0; k < 3; ++k) {
	    if ((error_terms[i][k] = calloc(frequencies,
			    sizeof(double complex))) == NULL) {
		(void)fprintf(stderr, "%s: calloc: %s\n",
			progname, strerror(errno));
		return NULL;
	    }
	}
    }

    /*
     * Generate the frequency vector.
     */
    if (frequencies == 1) {
	frequency_vector[0] = 1.0e6;
    } else if (frequencies == 2) {
	frequency_vector[0] = 0.0;
	frequency_vector[1] = 1.0e6;
    } else {
	frequency_vector[0] = 0.0;
	for (int i = 1; i < frequencies; ++i) {
	    frequency_vector[i] = pow(1.0e6,
		    (double)(i - 1) / (double)(frequencies - 2));
	}
    }
    if (vnacal_calset_set_frequency_vector(vcsp, frequency_vector) == -1) {
	(void)fprintf(stderr, "%s: vnacal_calset_set_frequency_vector: %s\n",
		progname, strerror(errno));
	return NULL;
    }

    /*
     * Generate the reference gamma values.
     */
    for (int findex = 0; findex < frequencies; ++findex) {
	for (int reference = 0; reference < 3; ++reference) {
	    int singular;

	    do {
		references[reference][findex] = crandn();
		singular = 0;
		for (int i = 0; i < reference; ++i) {
		    if (cabs(references[reference][findex] -
		             references[i][findex]) < EPS) {
			singular = 1;
			break;
		    }
		}
	    } while (singular);
	}
    }
    if (vnacal_calset_set_reference_vector(vcsp, 0, frequencies,
		frequency_vector, references[0]) == -1) {
	(void)fprintf(stderr, "%s: vnacal_calset_set_reference_vector: %s\n",
		progname, strerror(errno));
	return NULL;
    }
    if (vnacal_calset_set_reference_vector(vcsp, 1, frequencies,
		frequency_vector, references[1]) == -1) {
	(void)fprintf(stderr, "%s: vnacal_calset_set_reference_vector: %s\n",
		progname, strerror(errno));
	return NULL;
    }
    if (vnacal_calset_set_reference_vector(vcsp, 2, frequencies,
		frequency_vector, references[2]) == -1) {
	(void)fprintf(stderr, "%s: vnacal_calset_set_reference_vector: %s\n",
		progname, strerror(errno));
	return NULL;
    }

    /*
     * For each frequency...
     */
    for (int findex = 0; findex < frequencies; ++findex) {

	/*
	 * Generate the diagonal terms.
	 */
	for (int column = 0; column < ndiagonal; ++column) {
	    double complex e00, e10e01, e11;
	    double complex **epp = error_terms[column * columns + column];

	    /*
	     * Generate e00, e10e01 and e11.
	     */
	    e00 = crandn();
	    do {
		e10e01 = crandn();
	    } while (cabs(e10e01) <= EPS);
	    e11 = crandn();

	    /*
	     * Compute vcd_data_vectors for each reference gamma.
	     */
	    for (int reference = 0; reference < 3; ++reference) {
		double complex gamma;

		gamma = _vnacal_calset_get_reference(vcsp, reference, findex);
		cdata[column][column][reference][findex] = e00 +
		    e10e01 * gamma / (1.0 - e11 * gamma);
	    }
	    epp[0][findex] = e00;
	    epp[1][findex] = e10e01;
	    epp[2][findex] = e11;
	}

	/*
	 * Generate the off-diagonal terms.
	 */
	for (int row = 0; row < rows; ++row) {
	    for (int column = 0; column < columns; ++column) {
		double complex **epp;
		double complex e30, e10e32, e22 = 0.0;

		/*
		 * Skip diagonal entries.
		 */
		if (column == row) {
		    continue;
		}
		epp = error_terms[row * columns + column];

		/*
		 * Generate e30 and e10e32.
		 */
		e30 = crandn();
		do {
		    e10e32 = crandn();
		} while (cabs(e10e32) <= EPS);
		e22 = 0.0;

		/*
		 * If this column has a diagonal entry, generate e22 and use
		 * the diagonal terms to calculate vcd_sii_through_vector and
		 * vcd_sji_through_vector for full 6 error terms.  Otherwise,
		 * the VNA cannot calculate e22 and we can only supply 5 terms.
		 */
		if (column < rows) {
		    double complex **temp =
			error_terms[column * columns + column];
		    double complex e00	  = temp[0][findex];
		    double complex e10e01 = temp[1][findex];
		    double complex e11	  = temp[2][findex];

		    e22 = crandn();
		    cdata[row][column][0][findex] = e00 + e10e01 * e22 /
			(1.0 - e11 * e22);
		    cdata[row][column][1][findex] = e30 + e10e32 /
			(1.0 - e11 * e22);
		    cdata[row][column][2][findex] = e30;

		} else {
		    cdata[row][column][0][findex] = 0.0;
		    cdata[row][column][1][findex] = e30 + e10e32;
		    cdata[row][column][2][findex] = e30;
		}
		epp[0][findex] = e30;
		epp[1][findex] = e10e32;
		epp[2][findex] = e22;
	    }
	}
    }
    for (int row = 0; row < rows; ++row) {
	for (int column = 0; column < columns; ++column) {
	    for (int term = 0; term < 3; ++term) {
		if (vnacal_calset_add_vector(vcsp, row, column,
			    term, cdata[row][column][term]) == -1) {
		    (void)fprintf(stderr, "%s: vnacal_calset_add_vector: "
			    "%s\n", progname, strerror(errno));
		    return NULL;
		}
	    }
	}
    }
    return error_terms;
}

/*
 * free_error_terms: free the error terms matrix
 *   @error_terms: serialized error term matrix
 *   @vcsp: vnacal_calset_t for matrix dimensions
 */
static void free_error_terms(double complex *(*error_terms)[3],
			     vnacal_calset_t *vcsp)
{
    int ncells = vcsp->vcs_rows * vcsp->vcs_columns;

    if (error_terms != NULL) {
	for (int i = 0; i < ncells; ++i) {
	    for (int k = 0; k < 3; ++k) {
		free((void *)error_terms[i][k]);
	    }
	}
	free((void *)error_terms);
    }
}

/*
 * alloc_matrix_of_vectors: allocate a matrix of per-frequency vectors
 *   @ncells: rows x columns
 *   @frequencies: length of frequency vectors
 */
static double complex **alloc_matrix_of_vectors(int ncells, int frequencies)
{
    double complex **result;

    if ((result = (double complex **)calloc(ncells,
		    sizeof(double complex *))) == NULL) {
	(void)fprintf(stderr, "%s: calloc: %s\n",
		progname, strerror(errno));
	return NULL;
    }
    for (int cell = 0; cell < ncells; ++cell) {
	if ((result[cell] = (double complex *)calloc(frequencies,
			sizeof(double complex))) == NULL) {
	    (void)fprintf(stderr, "%s: calloc: %s\n",
		    progname, strerror(errno));
	    return NULL;
	}
    }
    return result;
}

/*
 * free_matrix_of_vectors: allocate a matrix of per-frequency vectors
 *   @matrix: serialized matrix of pointers to complex
 *   @ncells: rows x columns
 */
static void free_matrix_of_vectors(double complex **matrix, int ncells)
{
    for (int cell = 0; cell < ncells; ++cell) {
	free((void *)matrix[cell]);
    }
    free((void *)matrix);
}

/*
 * error_fn: error reporting function
 *   @msg: error message
 *   @arg: (unused)
 */
static void error_fn(const char *msg, void *arg)
{
    (void)fprintf(stderr, "%s: %s\n", progname, msg);
}

/*
 * test_vnacal_new_helper
 *   @trial: trial number
 *   @rows: rows in calibration matrix
 *   @columns: columns in calibration matrix
 *   @frequencies: number of frequency points
 */
static test_result_type test_vnacal_new_helper(int trial, int rows,
	int columns, int frequencies)
{
    vnacal_calset_t *vcsp = NULL;
    double complex *(*error_terms)[3] = NULL;
    vnacal_t *vcp = NULL;
    vnacal_etermset_t *etsp;
    test_result_type result = T_SKIPPED;

    /*
     * If -v, print the test header.
     */
    if (opt_v) {
	(void)printf("Test vnacal_create: trial %3d size %d x %d\n",
		trial, rows, columns);
    }

    /*
     * Generate the error terms and calibration measurements.
     */
    if ((vcsp = vnacal_calset_alloc("test", rows, columns,
		    frequencies, error_fn, NULL)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_calset_alloc: "
		"%s\n", progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if ((error_terms = gen_error_terms(vcsp)) == NULL) {
	result = T_FAIL;
	goto out;
    }


    /*
     * Create a new vnacal_t based on the calibration measurements.
     */
    if ((vcp = vnacal_create(1, &vcsp, error_fn, NULL)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_create: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }

    /*
     * Verify the error terms.
     */
    if (opt_v) {
	(void)printf("R C F ET\n");
    }
    etsp = vcp->vc_set_vector[0];
    for (int findex = 0; findex < vcsp->vcs_frequencies;
	    ++findex) {
	for (int row = 0; row < vcsp->vcs_rows; ++row) {
	    for (int column = 0; column < vcsp->vcs_columns; ++column) {
		int cell = row * vcsp->vcs_columns + column;
		double complex **epp = error_terms[cell];
		vnacal_error_terms_t *etp = &etsp->ets_error_term_matrix[cell];

		for (int k = 0; k < 3; ++k) {
		    if (opt_v) {
			(void)printf("%d %d %d %-6s %+e%+ei %+e%+ei\n",
				row, column, findex, error_term_names[k],
				creal(etp->et_data_vectors[k][findex]),
				cimag(etp->et_data_vectors[k][findex]),
				creal(epp[k][findex]),
				cimag(epp[k][findex]));
		    }
		    if (!isequal(etp->et_data_vectors[k][findex],
				epp[k][findex])) {
			if (opt_a) {
			    assert(!"data miscompare");
			}
			result = T_FAIL;
			goto out;
		    }
		}
	    }
	}
    }
    if (opt_v) {
	(void)printf("\n");
    }
    result = T_PASS;

out:
    if (vcp != NULL) {
	if (error_terms != NULL) {
	    free_error_terms(error_terms, vcsp);
	    error_terms = NULL;
	}
	vnacal_free(vcp);
	vcp = NULL;
    }
    if (vcsp != NULL) {
	vnacal_calset_free(vcsp);
	vcsp = NULL;
    }
    return result;
}

/*
 * test_vnacal_new
 */
static void test_vnacal_new()
{
    static const int sizes[] = { 1, 2, 3, 4 };
    test_result_type result = T_SKIPPED;

    for (int trial = 1; trial <= NTRIALS; ++trial) {
	for (int si = 0; si < sizeof(sizes) / sizeof(int); ++si) {
	    for (int sj = 0; sj < sizeof(sizes) / sizeof(int); ++sj) {
		int m = sizes[si];
		int n = sizes[sj];

		result = test_vnacal_new_helper(trial, m, n, 2);
		if (result != T_PASS)
		    goto out;
	    }
	}
    }
    result = T_PASS;

out:
    report_test_result("vnacal_create", result);
}

/*
 * test_vnacal_apply_helper
 *   @trial: trial number
 *   @vrows: rows in calibration matrix
 *   @vcolumns: columns in calibration matrix
 *   @drows: rows in DUT matrix
 *   @dcolumns: columns in DUT matrix
 *   @frequencies: number of frequency points
 *   @map_flag: map DUT ports to VNA ports (required if DUT larger than cal)
 */
static test_result_type test_vnacal_apply_helper(int trial,
	int vrows, int vcolumns, int drows, int dcolumns,
	int frequencies, bool map_flag)
{
    vnacal_calset_t *vcsp = NULL;
    vnacal_input_t *vip = NULL;
    double complex *(*error_terms)[3] = NULL;
    vnacal_t *vcp = NULL;
    int map[drows * dcolumns];
    double complex **actual_matrix   = NULL;
    double complex **measured_matrix = NULL;
    vnadata_t *output_matrix = NULL;
    test_result_type result = T_SKIPPED;

    /*
     * If -v, print the test header.
     */
    if (opt_v) {
	(void)printf("Test vnacal_input_apply: trial %3d cal size (%d x %d) "
		"S size (%d x %d) map %d\n", trial, vrows, vcolumns,
		drows, dcolumns, (int)map_flag);
    }

    /*
     * Generate the error terms and calibration measurements.
     */
    if ((vcsp = vnacal_calset_alloc("test", vrows, vcolumns,
		    frequencies, error_fn, NULL)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_calset_alloc: "
		"%s\n", progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if ((error_terms = gen_error_terms(vcsp)) == NULL) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Create a new vnacal_t based on the calibration measurements.
     */
    if ((vcp = vnacal_create(1, &vcsp, error_fn, NULL)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_create: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (opt_v) {
	vnacal_etermset_t *etsp = vcp->vc_set_vector[0];

	(void)printf("error terms:\n");
	(void)printf("R C F ET\n");
	for (int findex = 0; findex < frequencies; ++findex) {
	    for (int ci = 0; ci < vrows; ++ci) {
		for (int cj = 0; cj < vcolumns; ++cj) {
		    int c_cell = ci * vcolumns + cj;
		    double complex **epp = error_terms[c_cell];
		    vnacal_error_terms_t *etp;

		    etp = &etsp->ets_error_term_matrix[c_cell];
		    for (int k = 0; k < 3; ++k) {
			(void)printf("%d %d %d %-6s %+e%+ei %+e%+ei\n",
				ci, cj, findex, error_term_names[k],
				creal(etp->et_data_vectors[k][findex]),
				cimag(etp->et_data_vectors[k][findex]),
				creal(epp[k][findex]),
				cimag(epp[k][findex]));
		    }
		}
	    }
	}
	(void)printf("\n");
    }

    /*
     * If map_flag, generate a random map between S parameter
     * ports and VNA ports.
     */
    if (map_flag) {
	int c_ndiagonal = MIN(vrows, vcolumns);

	for (int row = 0; row < drows; ++row) {
	    for (int column = 0; column < dcolumns; ++column) {
		int cell = row * dcolumns + column;

		if (row == column) {
		    int c_diagonal = random() % c_ndiagonal;

		    map[cell] = c_diagonal * vcolumns + c_diagonal;

		} else if (vcolumns > 1) {
		    int vrow = random() % vrows;
		    int vcolumn = random() % (vcolumns - 1);

		    if (vcolumn >= vrow) {
			++vcolumn;
		    }
		    assert(vrow != vcolumn);
		    map[cell] = vrow * vcolumns + vcolumn;

		} else {
		    int vrow;

		    assert(vrows > 1);
		    vrow = random() % (vrows - 1) + 1;

		    map[cell] = vrow * vcolumns;
		}
	    }
	}
	if (opt_v) {
	    (void)printf("map:\n");
	    for (int i = 0; i < drows; ++i) {
		for (int j = 0; j < dcolumns; ++j) {
		    int cell = map[i * dcolumns + j];
		    int mrow = cell / dcolumns;
		    int mcol = cell % dcolumns;

		    (void)printf("   %2d %2d", mrow, mcol);
		}
		(void)printf("\n");
	    }
	    (void)printf("\n");
	}
    }

    /*
     * Allocate S parameter matrices.
     */
    actual_matrix   = alloc_matrix_of_vectors(drows * dcolumns, frequencies);
    measured_matrix = alloc_matrix_of_vectors(drows * dcolumns, frequencies);
    if (actual_matrix == NULL || measured_matrix == NULL) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Generate the "actual" S-parameters.
     */
    for (int row = 0; row < drows; ++row) {
	for (int column = 0; column < dcolumns; ++column) {
	    for (int findex = 0; findex < vcsp->vcs_frequencies; ++findex) {
		actual_matrix[row * dcolumns + column][findex] = crandn();
	    }
	}
    }
    if (opt_v) {
	(void)printf("actual_matrix:\n");
	(void)printf("R C F\n");
	for (int findex = 0; findex < vcsp->vcs_frequencies; ++findex) {
	    for (int row = 0; row < drows; ++row) {
		for (int column = 0; column < dcolumns; ++column) {
		    double complex v;

		    v = actual_matrix[row * dcolumns + column][findex];
		    (void)printf("%d %d %d %+e%+ei\n",
			row, column, findex, creal(v), cimag(v));
		}
	    }
	}
	(void)printf("\n");
    }

    /*
     * Generate the "measured" S-parameters given actual and error terms.
     */
    for (int findex = 0; findex < vcsp->vcs_frequencies; ++findex) {

	/*
	 * Init measured to zero.
	 */
	for (int i = 0; i < drows * dcolumns; ++i) {
	    measured_matrix[i][findex] = 0.0;
	}

	/*
	 * For each column, (each driven port) find the corresponding
	 * column in the measured_matrix.
	 */
	for (int dcolumn = 0; dcolumn < dcolumns; ++dcolumn) {
	    double complex a[drows * drows];
	    double complex x[drows];
	    double complex b[drows];
	    double complex d;
#define A(i, j) (a[(i) * drows + (j)])

	    /*
	     * We start by forming a (drows x drows) matrix A and column
	     * vector b which we'll use to solve for column vector x.
	     * The A matrix corresponds to the square portion of the
	     * S parameter matrix.  The columns of A correspond to the
	     * elements of the b and x vectors and to the elements of the
	     * current column of the "measured" matrix.  Column vector
	     * x that we're solving for represents the voltage out of
	     * each DUT port, taking into account the reflections due to
	     * errors at each of the other ports.  It's easy to calculate
	     * the column of the "measured" matrix once we know x.
	     *
	     * Initialize A to the identity matrix and b to the current
	     * column in the actual_matrix.
	     */
	    for (int i = 0; i < drows; ++i) {
		for (int j = 0; j < drows; ++j) {
		    A(i, j) = (i == j) ? 1.0 : 0.0;
		}
		b[i] = actual_matrix[i * dcolumns + dcolumn][findex];
	    }
	    /*
	     * Make A = I - S E, where E is a diagonal matrix made of
	     * the e11/e22 error terms for this column.
	     */
	    for (int j = 0; j < drows; ++j) {	/* column in A */
		double complex e11;

		if (j < dcolumns) {
		    int cell;

		    /*
		     * Find the error parameters for row j in this
		     * DUT column.  Note that j is the index of a row in
		     * the S matrix and in the X column vector, but it's
		     * the index of a column in the A matrix.  We start
		     * with the assumption that cells of the DUT matrix
		     * are 1:1 with cells in the calibration matrix,
		     * then if mapping is enabled, we apply the map.
		     */
		    cell = j * dcolumns + dcolumn;
		    if (map_flag) {
			cell = map[cell];
		    }
		    assert(cell < vrows * vcolumns);
		    e11 = error_terms[cell][2][findex];
		    for (int i = 0; i < drows; ++i) {
			A(i, j) -= actual_matrix[i*dcolumns + j][findex] * e11;
		    }
		}
	    }
	    if (opt_v) {
		(void)printf("findex %d column %d:\n", findex, dcolumn);
		(void)printf("a:\n");
		cmatrix_print(a, drows, drows);
		(void)printf("b:\n");
		cmatrix_print(b, drows, 1);
	    }
	    /*
	     * Find x = A^-1 b.
	     */
	    d = _vnacommon_mldivide(x, a, b, drows, 1);
	    if (cabs(d) <= EPS)	 {
		(void)fprintf(stderr, "%s: test_vnacal_mrdivide: warning: "
			"skipping nearly singular test matrix\n",
			progname);
		result = T_SKIPPED;
		goto out;
	    }
	    if (opt_v) {
		(void)printf("x:\n");
		cmatrix_print(x, drows, 1);
	    }
	    /*
	     * From x, calculate the "measured" S-parameters for this
	     * column
	     */
	    for (int drow = 0; drow < drows; ++drow) {
		int cell = drow * dcolumns + dcolumn;
		double complex e00;
		double complex e10e01;

		if (map_flag) {
		    cell = map[cell];
		}
		assert(cell < vrows * vcolumns);
		e00    = error_terms[cell][0][findex];
		e10e01 = error_terms[cell][1][findex];
		measured_matrix[drow * dcolumns + dcolumn][findex] =
		    e00 + e10e01 * x[drow];
	    }
	}
    }
    if (opt_v) {
	(void)printf("measured_matrix:\n");
	(void)printf("R C F\n");
	for (int findex = 0; findex < vcsp->vcs_frequencies; ++findex) {
	    for (int row = 0; row < drows; ++row) {
		for (int column = 0; column < dcolumns; ++column) {
		    double complex v = measured_matrix[row * dcolumns +
			column][findex];

		    (void)printf("%d %d %d %+e%+ei\n",
			row, column, findex, creal(v), cimag(v));
		}
	    }
	}
	(void)printf("\n");
    }

    /*
     * Create the vnacal_input_t.
     */
    vip = vnacal_input_alloc(vcp, /*set=*/0, drows, dcolumns, frequencies);
    if (vip == NULL) {
	(void)fprintf(stderr, "%s: vnacal_input_alloc: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_input_set_frequency_vector(vip,
		vcsp->vcs_frequency_vector) == -1) {
	(void)fprintf(stderr, "%s: vnacal_input_set_frequency_vector: "
		"%s\n", progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    for (int drow = 0; drow < drows; ++drow) {
	for (int dcolumn = 0; dcolumn < dcolumns; ++dcolumn) {
	    int cell = dcolumns * drow + dcolumn;

	    if (!map_flag) {
		if (vnacal_input_add_vector(vip, drow, dcolumn,
			    measured_matrix[cell]) == -1) {
		    (void)fprintf(stderr, "%s: vnacal_input_add_vector: "
			    "drow %d dcolumn %d: %s\n",
			    progname, drow, dcolumn, strerror(errno));
		    result = T_FAIL;
		    goto out;
		}
	    } else {
		int vrow    = map[cell] / vcsp->vcs_columns;
		int vcolumn = map[cell] % vcsp->vcs_columns;

		if (vnacal_input_add_mapped_vector(vip, vrow, vcolumn,
			    drow, dcolumn, measured_matrix[cell]) == -1) {
		    (void)fprintf(stderr, "%s: vnacal_input_add_vector: "
			    "drow %d dcolumn %d: %s\n",
			    progname, drow, dcolumn, strerror(errno));
		    result = T_FAIL;
		    goto out;
		}
	    }
	}
    }

    /*
     * Get the computed S-parameters.
     */
    if ((output_matrix = vnadata_alloc()) == NULL)  {
	(void)fprintf(stderr, "%s: vnadata_alloc: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_input_apply(vip, output_matrix) == -1) {
	(void)fprintf(stderr, "%s: vnacal_input_apply: "
		"%s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (opt_v) {
	(void)printf("computed_vector:\n");
	(void)printf("R C F\n");
	for (int findex = 0; findex < vcsp->vcs_frequencies; ++findex) {
	    for (int row = 0; row < drows; ++row) {
		for (int column = 0; column < dcolumns; ++column) {
		    double complex v;

		    v = vnadata_get_cell(output_matrix, findex, row, column);
		    (void)printf("%d %d %d %+e%+ei\n",
			row, column, findex, creal(v), cimag(v));
		}
	    }
	}
	(void)printf("\n");
    }

    /*
     * Check the result.
     */
    for (int i = 0; i < drows; ++i) {
	for (int j = 0; j < dcolumns; ++j) {
	    for (int findex = 0; findex < vcsp->vcs_frequencies; ++findex) {
		double complex v;
		double dy;

		v = vnadata_get_cell(output_matrix, findex, i, j);
		dy = cabs(v - actual_matrix[i * dcolumns + j][findex]);
		if (dy >= EPS) {
		    if (opt_a) {
			assert(!"data miscompare");
		    }
		    result = T_FAIL;
		    goto out;
		}
	    }
	}
    }
    result = T_PASS;

out:
    if (vip != NULL) {
	vnacal_input_free(vip);
	vip = NULL;
    }
    free_matrix_of_vectors(measured_matrix, drows * dcolumns);
    free_matrix_of_vectors(actual_matrix,   drows * dcolumns);
    if (vcp != NULL) {
	if (error_terms != NULL) {
	    free_error_terms(error_terms, vcsp);
	    error_terms = NULL;
	}
	vnacal_free(vcp);
	vcp = NULL;
    }
    if (vcsp != NULL) {
	vnacal_calset_free(vcsp);
	vcsp = NULL;
    }
    return result;
}
#undef B
#undef S
#undef U

/*
 * test_vnacal_input_apply: test vnacal_input_apply and vnacal_apply_with_map
 */
static void test_vnacal_input_apply()
{
    static const int sizes[] = { 1, 2, 3, 4 };
    test_result_type result = T_SKIPPED;
    bool pass = false;

    for (int trial = 1; trial <= NTRIALS; ++trial) {
	for (int si = 0; si < sizeof(sizes) / sizeof(int); ++si) {
	    for (int sj = 0; sj < sizeof(sizes) / sizeof(int); ++sj) {
		int drows = sizes[si];
		int dcolumns = sizes[sj];

		result = test_vnacal_apply_helper(trial,
			drows, dcolumns, drows, dcolumns, 2, false);
		switch (result) {
		case T_PASS:
		    pass = true;
		    break;
		case T_SKIPPED:
		    break;
		case T_FAIL:
		default:
		    goto out;
		}
		result = test_vnacal_apply_helper(trial,
			2, 1, drows, dcolumns, 2, true);
		switch (result) {
		case T_PASS:
		    pass = true;
		    break;
		case T_SKIPPED:
		    break;
		case T_FAIL:
		default:
		    goto out;
		}
		result = test_vnacal_apply_helper(trial,
			1, 2, drows, dcolumns, 2, true);
		switch (result) {
		case T_PASS:
		    pass = true;
		    break;
		case T_SKIPPED:
		    break;
		case T_FAIL:
		default:
		    goto out;
		}
		result = test_vnacal_apply_helper(trial,
			2, 2, drows, dcolumns, 2, true);
		switch (result) {
		case T_PASS:
		    pass = true;
		    break;
		case T_SKIPPED:
		    break;
		case T_FAIL:
		default:
		    goto out;
		}
	    }
	}
    }
    result = pass ? T_PASS : T_SKIPPED;

out:
    report_test_result("vnacal_input_apply", result);
}

/*
 * Test Strings for vnacal_property_set
 */
static const char property_foo_value[] = "1234567890";
static const char property_bar_value[] = "abcdefghijkl\nmnopqrstuvwxyz";
static const char property3_value[] = "αβγδεζηθικλμνξοπρστυφχψω";

/*
 * test_vnacal_save
 */
static void test_vnacal_save()
{
    vnacal_calset_t *cal_sets[2] = { NULL, NULL };
    double complex *(*error_terms0)[3] = NULL;
    double complex *(*error_terms1)[3] = NULL;
    vnacal_t *vcp = NULL;
    const char *cp_temp;
    test_result_type result = T_SKIPPED;

    /*
     * If -v, print the test header.
     */
    if (opt_v) {
	(void)printf("Test vnacal_save, vnacal_load\n");
    }

    /*
     * Generate the first calibration set.
     */
    if ((cal_sets[0] = vnacal_calset_alloc("first-set", 2, 1, 20,
		    error_fn, NULL)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_calset_alloc: "
		"%s\n", progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if ((error_terms0 = gen_error_terms(cal_sets[0])) == NULL) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Generate the second calibration set.
     */
    if ((cal_sets[1] = vnacal_calset_alloc("second-set", 3, 5, 10,
		    error_fn, NULL)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_calset_alloc: "
		"%s\n", progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if ((error_terms1 = gen_error_terms(cal_sets[1])) == NULL) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Create a new vnacal_t based on the calibration measurements.
     */
    if ((vcp = vnacal_create(2, cal_sets, error_fn, NULL)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_create: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }

    /*
     * Set test properties.
     */
    if (vnacal_property_set(vcp, 0, "foo=999999999999") == -1) {
	(void)fprintf(stderr, "%s: vnacal_property_set: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_property_set(vcp, 0, "bar=%s", property_bar_value) == -1) {
	(void)fprintf(stderr, "%s: vnacal_property_set: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_property_set(vcp, 0, "foo=%s", property_foo_value) == -1) {
	(void)fprintf(stderr, "%s: vnacal_property_set: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_property_set(vcp, 1, "baz=!!!") == -1) {
	(void)fprintf(stderr, "%s: vnacal_property_set: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_property_set(vcp, 1, "property3=%s", property3_value) == -1) {
	(void)fprintf(stderr, "%s: vnacal_property_set: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_property_delete(vcp, 1, "baz") == -1) {
	(void)fprintf(stderr, "%s: vnacal_property_delete: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    for (int row = 0; row < cal_sets[0]->vcs_rows; ++row) {
	for (int column = 0; column < cal_sets[0]->vcs_columns; ++column) {
	    int cell = row * cal_sets[0]->vcs_columns + column;
	    int value = (cell + 1) %
		(cal_sets[0]->vcs_rows * cal_sets[0]->vcs_columns);

	    if (vnacal_property_set(vcp, 0, "switches[%d][%d]=%d",
			row, column, value) == -1) {
		(void)fprintf(stderr, "%s: vnacal_property_set: %s\n",
			progname, strerror(errno));
		result = T_FAIL;
		goto out;
	    }
	}
    }
    for (int row = 0; row < cal_sets[1]->vcs_rows; ++row) {
	for (int column = 0; column < cal_sets[1]->vcs_columns; ++column) {
	    int cell = row * cal_sets[1]->vcs_columns + column;
	    int value = (cell + 3) %
		(cal_sets[1]->vcs_rows * cal_sets[1]->vcs_columns);

	    if (vnacal_property_set(vcp, 1, "switches[%d][%d]=%d",
			row, column, value) == -1) {
		(void)fprintf(stderr, "%s: vnacal_property_set: %s\n",
			progname, strerror(errno));
		result = T_FAIL;
		goto out;
	    }
	}
    }

    /*
     * Save and free.
     */
    if (vnacal_set_dprecision(vcp, 7) == -1) {
	(void)fprintf(stderr, "%s: vnacal_set_precision: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_save(vcp, "vnacal-test.vnacal", ".testcal") == -1) {
	(void)fprintf(stderr, "%s: vnacal_save: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    vnacal_free(vcp);
    vcp = NULL;

    /*
     * Load
     */
    if ((vcp = vnacal_load("vnacal-test.vnacal", ".testcal",
		    error_fn, NULL)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_load: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_get_sets(vcp) != 2) {
	(void)printf("expected 2 sets; found %d\n", vnacal_get_sets(vcp));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_get_rows(vcp, 0) != 2) {
	(void)printf("expected 2 rows in set 0; found %d\n",
	    vnacal_get_rows(vcp, 0));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_get_columns(vcp, 0) != 1) {
	(void)printf("expected 1 column in set 0; found %d\n",
	    vnacal_get_rows(vcp, 0));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_get_frequencies(vcp, 0) != 20) {
	(void)printf("expected 20 frequencies in set 0; found %d\n",
	    vnacal_get_frequencies(vcp, 0));
	result = T_FAIL;
	goto out;
    }
    for (int i = 0; i < 2; ++i) {
	vnacal_etermset_t *etsp = vcp->vc_set_vector[0];
	vnacal_error_terms_t *etp = &etsp->ets_error_term_matrix[i];

	for (int j = 0; j < 3; ++j) {
	    for (int k = 0; k < 20; ++k) {
		if (!isequal(etp->et_data_vectors[j][k],
			    error_terms0[i][j][k])) {
		    result = T_FAIL;
		    goto out;
		}
	    }
	}
    }
    if (vnacal_get_rows(vcp, 1) != 3) {
	(void)printf("expected 3 rows in set 1; found %d\n",
	    vnacal_get_rows(vcp, 1));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_get_columns(vcp, 1) != 5) {
	(void)printf("expected 5 columns in set 1; found %d\n",
	    vnacal_get_rows(vcp, 1));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_get_frequencies(vcp, 1) != 10) {
	(void)printf("expected 10 frequencies in set 1; found %d\n",
	    vnacal_get_frequencies(vcp, 1));
	result = T_FAIL;
	goto out;
    }
    for (int i = 0; i < 3*5; ++i) {
	vnacal_etermset_t *etsp = vcp->vc_set_vector[1];
	vnacal_error_terms_t *etp = &etsp->ets_error_term_matrix[i];

	for (int j = 0; j < 3; ++j) {
	    for (int k = 0; k < 10; ++k) {
		if (!isequal(etp->et_data_vectors[j][k],
			    error_terms1[i][j][k])) {
		    result = T_FAIL;
		    goto out;
		}
	    }
	}
    }
    if ((cp_temp = vnacal_property_get(vcp, 0, "foo")) == NULL) {
	(void)printf("property \"foo\" in set 0 not found\n");
	result = T_FAIL;
	goto out;
    }
    if (strcmp(cp_temp, property_foo_value) != 0) {
	(void)printf("expected \"%s\" for property \"foo\"; found \"%s\"\n",
		property_foo_value, cp_temp);
	result = T_FAIL;
	goto out;
    }
    if ((cp_temp = vnacal_property_get(vcp, 0, "bar")) == NULL) {
	(void)printf("property \"bar\" in set 0 not found\n");
	result = T_FAIL;
	goto out;
    }
    if (strcmp(cp_temp, property_bar_value) != 0) {
	(void)printf("expected \"%s\" for property \"bar\"; found \"%s\"\n",
		property_bar_value, cp_temp);
	result = T_FAIL;
	goto out;
    }
    if ((cp_temp = vnacal_property_get(vcp, 0, "baz")) != NULL) {
	(void)printf("property \"baz\" not expected in set 0; "
		"found it with value \"%s\"\n", cp_temp);
	result = T_FAIL;
	goto out;
    }
    if ((cp_temp = vnacal_property_get(vcp, 1, "property3")) == NULL) {
	(void)printf("property \"property3\" in set 1 not found\n");
	result = T_FAIL;
	goto out;
    }
    if (strcmp(cp_temp, property3_value) != 0) {
	(void)printf("expected \"%s\" for property \"property3\"; "
		"found \"%s\"\n", property3_value, cp_temp);
	result = T_FAIL;
	goto out;
    }
    vnacal_free(vcp);
    vcp = NULL;
    result = T_PASS;

out:
    report_test_result("vnacal_save/vnacal_load", result);
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
    /*NOTREACHED*/
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
	    opt_v = true;
	    continue;
	default:
	    print_usage();
	    /*NOTREACHED*/
	}
	break;
    }
    test_vnacal_new();
    test_vnacal_input_apply();
    test_vnacal_save();

    exit(fail_count != 0);
    /*NOTREACHED*/
}
