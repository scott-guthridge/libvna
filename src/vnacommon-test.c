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
#include <unistd.h>
#include "vnacommon_internal.h"


#define PI	3.1415926535897932384626433832795
#define EPS	1.0e-4

#define N_MATRIX_TRIALS	100

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
static bool opt_a = false;
static bool opt_v = false;

/*
 * crandn: generate normally distributed complex numbers
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
 * cmatrix_multiply: find C = A x B
 *   @c: serialized result matrix, m x o
 *   @a: serialized A matrix, m x n
 *   @b: serialized B matrix, n x o
 *   @m: first dimension of C and A
 *   @n: second dimension of A, first dimension of B
 *   @o: second dimension of C and B
 */
static void cmatrix_multiply(double complex *c, const double complex *a,
	const double complex *b, int m, int n, int o)
{
#define A(i, j) (a[(i) * n + (j)])
#define B(i, j) (b[(i) * o + (j)])
#define C(i, j) (c[(i) * o + (j)])

    for (int i = 0; i < m; ++i) {
	for (int k = 0; k < o; ++k) {
	    double complex s = 0.0;

	    for (int j = 0; j < n; ++j) {
		s += A(i, j) * B(j, k);
	    }
	    C(i, k) = s;
	}
    }
}
#undef A
#undef B
#undef C

/*
 * cmatrix_print: print an m by n serialized complex matrix
 */
static void cmatrix_print(const char *tag, double complex *a, int m, int n)
{
    (void)printf("%s:\n", tag);
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
    (void)(void)fflush(stdout);
    if (result == T_FAIL) {
	++fail_count;
    }
}

/*
 * test_vnacommon_lu: test LU factorization
 */
static void test_vnacommon_lu()
{
    static const int sizes[] = { 1, 2, 3, 10 };
    double complex a[10 * 10];
    double complex t[10 * 10];
    int row_index[10];
    test_result_type result = T_SKIPPED;

#define A(i, j) (a[(i) * n + (j)])
#define T(i, j) (t[(i) * n + (j)])
    for (int trial = 1; trial <= N_MATRIX_TRIALS; ++trial) {
	for (int size_index = 0; size_index < sizeof(sizes) / sizeof(int);
		++size_index) {
	    int n = sizes[size_index];
	    double complex d;

	    /*
	     * If -v, print the test header.
	     */
	    if (opt_v) {
		(void)printf("Test vnacommon_lu: trial %3d size %d x %d\n",
			trial, n, n);
		(void)fflush(stdout);
	    }

	    /*
	     * Generate a random complex matrix in T and copy to A.
	     */
	    for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
		    A(i, j) = T(i, j) = crandn();
		}
	    }
	    if (opt_v) {
		cmatrix_print("a", a, n, n);
		(void)fflush(stdout);
	    }

	    /*
	     * Compute the in-place L U factorization.
	     */
	    errno = 0;
	    d = _vnacommon_lu(a, row_index, n);
	    if (opt_v) {
		cmatrix_print("LU factorization", a, n, n);
		(void)printf("determinant %8.5f%+8.5fj\n", creal(d), cimag(d));
		(void)printf("\n");
		(void)fflush(stdout);
	    }
	    if (cabs(d) < EPS) {
		(void)fprintf(stderr, "%s: test_vnacommon_lu: warning: "
			"skipping nearly singular test matrix\n", progname);
		continue;
	    }

	    /*
	     * Check the result.
	     */
	    for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
		    double complex s = (i <= j) ? A(i, j) : 0.0;
		    double dy;

		    for (int k = 0; k < i; ++k) {
			if (k <= j) {
			    s += A(i, k) * A(k, j);
			}
		    }
		    dy = cabs(s - T(row_index[i], j));
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
    }
    result = T_PASS;

out:
    report_test_result("LU Factorization", result);
}
#undef A
#undef T

/*
 * test_vnacommon_mldivide: test matrix left division
 */
static void test_vnacommon_mldivide()
{
    static const int sizes[] = { 1, 2, 3, 5 };
    double complex x[5 * 5];
    double complex a[5 * 5];
    double complex b[5 * 5];
    double complex t[5 * 5];
    test_result_type result = T_SKIPPED;

#define A(i, j) (a[(i) * m + (j)])
#define B(i, j) (b[(i) * n + (j)])
#define X(i, j) (x[(i) * n + (j)])
#define T(i, j) (t[(i) * n + (j)])
    for (int trial = 1; trial <= N_MATRIX_TRIALS; ++trial) {
	for (int si = 0; si < sizeof(sizes) / sizeof(int); ++si) {
	    for (int sj = 0; sj < sizeof(sizes) / sizeof(int); ++sj) {
		int m = sizes[si];
		int n = sizes[sj];
		double complex d;

		/*
		 * If -v, print the test header.
		 */
		if (opt_v) {
		    (void)printf("Test vnacommon_mldivide: trial %3d size "
			    "%d x %d\n", trial, m, n);
		    (void)fflush(stdout);
		}

		/*
		 * Generate A and T.  Multiply to find B.
		 */
		for (int i = 0; i < m; ++i) {
		    for (int j = 0; j < m; ++j) {
			A(i, j) = crandn();	/* m x m */
		    }
		    for (int j = 0; j < n; ++j) {
			T(i, j) = crandn();	/* m x n */
		    }
		}
		cmatrix_multiply(b, a, t, m, m, n);
		if (opt_v) {
		    cmatrix_print("a", a, m, m);
		    cmatrix_print("b", b, m, n);
		    cmatrix_print("t", t, m, n);
		    (void)fflush(stdout);
		}

		/*
		 * Solve for X.
		 */
		d = _vnacommon_mldivide(x, a, b, m, n);
		if (opt_v) {
		    cmatrix_print("x", x, m, n);
		    (void)printf("determinant %8.5f%+8.5fj\n",
			    creal(d), cimag(d));
		    (void)printf("\n");
		    (void)fflush(stdout);
		}
		if (cabs(d) < EPS) {
		    (void)fprintf(stderr, "%s: test_vnacommon_mldivide: "
			    "warning: skipping nearly singular test matrix\n",
			    progname);
		    continue;
		}

		/*
		 * Check the result.
		 */
		for (int i = 0; i < m; ++i) {
		    for (int j = 0; j < n; ++j) {
			double dy = cabs(X(i, j) - T(i, j));

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
	}
    }
    result = T_PASS;

out:
    report_test_result("Matrix Left Division", result);
}
#undef A
#undef B
#undef X
#undef T

/*
 * test_vnacommon_mrdivide: test matrix right division
 */
static void test_vnacommon_mrdivide()
{
    static const int sizes[] = { 1, 2, 3, 5 };
    double complex x[5 * 5];
    double complex a[5 * 5];
    double complex b[5 * 5];
    double complex t[5 * 5];
    test_result_type result = T_SKIPPED;

#define A(i, j) (a[(i) * n + (j)])
#define B(i, j) (b[(i) * n + (j)])
#define X(i, j) (x[(i) * n + (j)])
#define T(i, j) (t[(i) * n + (j)])
    for (int trial = 1; trial <= N_MATRIX_TRIALS; ++trial) {
	for (int si = 0; si < sizeof(sizes) / sizeof(int); ++si) {
	    for (int sj = 0; sj < sizeof(sizes) / sizeof(int); ++sj) {
		int m = sizes[si];
		int n = sizes[sj];
		double complex d;

		/*
		 * If -v, print the test header.
		 */
		if (opt_v) {
		    (void)printf("Test vnacommon_mrdivide: trial %3d size "
			    "%d x %d\n", trial, m, n);
		    (void)fflush(stdout);
		}

		/*
		 * Generate A and T.  Multiply to find B.
		 */
		for (int i = 0; i < n; ++i) {
		    for (int j = 0; j < n; ++j) {
			A(i, j) = crandn();	/* n x n */
		    }
		}
		for (int i = 0; i < m; ++i) {
		    for (int j = 0; j < n; ++j) {
			T(i, j) = crandn();	/* m x n */
		    }
		}
		cmatrix_multiply(b, t, a, m, n, n);
		if (opt_v) {
		    cmatrix_print("a", a, n, n);
		    cmatrix_print("b", b, m, n);
		    cmatrix_print("t", t, m, n);
		    (void)fflush(stdout);
		}

		/*
		 * Solve for X.
		 */
		d = _vnacommon_mrdivide(x, b, a, m, n);
		if (opt_v) {
		    cmatrix_print("x", x, m, n);
		    (void)printf("determinant %8.5f%+8.5fj\n",
			    creal(d), cimag(d));
		    (void)printf("\n");
		    (void)fflush(stdout);
		}
		if (cabs(d) < EPS) {
		    (void)fprintf(stderr, "%s: test_vnacommon_mrdivide: "
			    "warning: skipping nearly singular test matrix\n",
			    progname);
		    continue;
		}

		/*
		 * Check the result.
		 */
		for (int i = 0; i < m; ++i) {
		    for (int j = 0; j < n; ++j) {
			double dy = cabs(X(i, j) - T(i, j));

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
	}
    }
    result = T_PASS;

out:
    report_test_result("Matrix Right Division", result);
}
#undef A
#undef B
#undef X
#undef T

/*
 * test_vnacommon_minverse: test matrix inverse
 */
static void test_vnacommon_minverse()
{
    static const int sizes[] = { 1, 2, 3, 5 };
    double complex a[5 * 5];
    double complex t[5 * 5];
    double complex x[5 * 5];
    test_result_type result = T_SKIPPED;

#define A(i, j) (a[(i) * n + (j)])
#define T(i, j) (t[(i) * n + (j)])
#define X(i, j) (x[(i) * n + (j)])
    for (int trial = 1; trial <= N_MATRIX_TRIALS; ++trial) {
	for (int si = 0; si < sizeof(sizes) / sizeof(int); ++si) {
	    int n = sizes[si];
	    double complex d;

	    /*
	     * If -v, print the test header.
	     */
	    if (opt_v) {
		(void)printf("Test vnacommon_minverse: trial %3d size "
			"%d x %d\n", trial, n, n);
		(void)fflush(stdout);
	    }

	    /*
	     * Generate A and copy to T.
	     */
	    for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
		    A(i, j) = crandn();	/* n x n */
		    T(i, j) = A(i, j);
		}
	    }
	    if (opt_v) {
		cmatrix_print("a", a, n, n);
		(void)fflush(stdout);
	    }

	    /*
	     * Find X = T^-1
	     */
	    d = _vnacommon_minverse(x, t, n);
	    if (opt_v) {
		cmatrix_print("x", x, n, n);
		(void)printf("determinant %8.5f%+8.5fj\n",
			creal(d), cimag(d));
		(void)printf("\n");
		(void)fflush(stdout);
	    }
	    if (cabs(d) < EPS) {
		(void)fprintf(stderr, "%s: test_vnacommon_mldivide: warning: "
			"skipping nearly singular test matrix\n",
			progname);
		continue;
	    }

	    /*
	     * Find T = A * X and check that the result is the
	     * identity matrix.
	     */
	    cmatrix_multiply(t, a, x, n, n, n);
	    for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
		    double dy = cabs((i == j ? 1.0 : 0.0) - T(i, j));

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
    }
    result = T_PASS;

out:
    report_test_result("Matrix Inverse", result);
}
#undef X
#undef T
#undef A

/*
 * test_qrd: test QR decomposition
 */
static void test_vnacommon_qrd()
{
    test_result_type result = T_SKIPPED;

    for (int trial = 1; trial <= N_MATRIX_TRIALS; ++trial) {
	for (int m = 1; m <= 5; ++m) {
	    for (int n = 1; n <= 5; ++n) {
		int diagonals = MIN(m, n);
		double complex a[m][n];
		double complex t[m][n];
		double complex d[diagonals];
		double complex q[m][m];
		double complex r[m][n];

		/*
		 * If -v, print the test header.
		 */
		if (opt_v) {
		    (void)printf("Test vnacommon_qrd: trial %3d "
			    "size %d x %d\n", trial, m, n);
		    (void)fflush(stdout);
		}

		/*
		 * Fill A with random numbers and copy to R.
		 * Decompose into factored Q_n and R components.
		 */
		for (int i = 0; i < m; ++i) {
		    for (int j = 0; j < n; ++j) {
			r[i][j] = a[i][j] = crandn();
		    }
		}
		_vnacommon_qrd(*r, d, m, n);
		if (opt_v) {
		    cmatrix_print("a",  *a, m, n);
		    cmatrix_print("qr", *r, m, n);
		    cmatrix_print("d",  d, 1, diagonals);
		    (void)fflush(stdout);
		}

		/*
		 * Initialize q to the identity matrix.
		 */
		for (int i = 0; i < m; ++i) {
		    for (int j = 0; j < m; ++j) {
			q[i][j] = (i == j) ? 1.0 : 0.0;
		    }
		}

		/*
		 * Form Q.
		 */
		for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
		    for (int i = 0; i < m; ++i) {
			double complex s = 0.0;

			for (int j = diagonal; j < m; ++j) {
			    s += q[i][j] * r[j][diagonal];
			}
			for (int j = diagonal; j < m; ++j) {
			    q[i][j] -= 2.0 * s * conj(r[j][diagonal]);
			}
		    }
		}
		if (opt_v) {
		    cmatrix_print("q", *q, m, m);
		    (void)fflush(stdout);
		}

		/*
		 * Form R.
		 */
		for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
		    r[diagonal][diagonal] = d[diagonal];
		}
		for (int i = 0; i < diagonals; ++i) {
		    for (int j = 0; j < i; ++j) {
			r[i][j] = 0.0;
		    }
		}
		for (int i = n; i < m; ++i) {
		    for (int j = 0; j < n; ++j) {
			r[i][j] = 0.0;
		    }
		}
		if (opt_v) {
		    cmatrix_print("r", *r, m, n);
		    (void)fflush(stdout);
		}

		/*
		 * Test that Q Q' is the identity matrix.
		 */
		for (int i = 0; i < m; ++i) {
		    for (int j = 0; j < m; ++j) {
			double complex s = 0.0;

			for (int k = 0; k < m; ++k) {
			    s += q[i][k] * conj(q[j][k]);
			}
			s -= (i == j) ? 1.0 : 0.0;
			if (cabs(s) > EPS) {
			    if (opt_a) {
				assert(!"data miscompare");
			    }
			    result = T_FAIL;
			    goto out;
			}
		    }
		}

		/*
		 * Test that Q R == A
		 */
		cmatrix_multiply(*t, *q, *r, m, m, n);
		for (int i = 0; i < m; ++i) {
		    for (int j = 0; j < n; ++j) {
			double dy;

			dy = cabs(t[i][j] - a[i][j]);
			if (dy > EPS) {
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
    }
    result = T_PASS;

out:
    report_test_result("QR Decomposition", result);
}

/*
 * qrsolve_helper: generate a random system of equations and solve
 */
static int qrsolve_helper(double complex *x, double complex *a,
	double complex *b, int m, int n, int o)
{
    double complex u[m][n];
    double complex v[m][o];

    /*
     * Generate random matrices A and B, and make copies
     * in U and V, respectively.
     */
    for (int i = 0; i < m; ++i) {
	for (int j = 0; j < n; ++j) {
	    a[i * n + j] = u[i][j] = crandn();
	}
	for (int k = 0; k < o; ++k) {
	    b[i * o + k] = v[i][k] = crandn();
	}
    }

    /*
     * Solve the system.  This call destroys both u and v.
     */
    return _vnacommon_qrsolve(x, &u[0][0], &v[0][0], m, n, o);
}

/*
 * find_axb_error: find the squared error in A * X = B
 */
static double find_axb_error(const double complex *a, const double complex *x,
	const double complex *b, int m, int n, int o)
{
    double squared_error = 0.0;

    for (int k = 0; k < o; ++k) {
	for (int i = 0; i < m; ++i) {
	    double complex s = 0.0;
	    double e;

	    for (int j = 0; j < n; ++j) {
		s += a[i * n + j] * x[j * o + k];
	    }
	    e = cabs(s - b[i * o + k]);
	    squared_error += e * e;
	}
    }
    return squared_error;
}

/*
 * Test vnacommon_qrsolve.
 */
static void test_vnacommon_qrsolve()
{
    test_result_type result = T_SKIPPED;

    for (int trial = 1; trial <= N_MATRIX_TRIALS; ++trial) {
	/*
	 * Test square coefficient matrices.
	 */
	for (int n = 1; n <= 10; ++n) {
	    const int o = 3;
	    double complex a[n][n];
	    double complex b[n][o];
	    double complex x[n][o];
	    int rank;

	    /*
	     * If -v, print the test header.
	     */
	    if (opt_v) {
		(void)printf("Test vnacommon_qrsolve: trial %3d "
			"size %d x %d\n", trial, n, n);
		(void)fflush(stdout);
	    }

	    /*
	     * Generate random matrices A and B, and solve for X.
	     */
	    rank = qrsolve_helper(&x[0][0], &a[0][0], &b[0][0], n, n, o);
	    if (opt_v) {
		cmatrix_print("a", *a, n, n);
		cmatrix_print("b", *b, n, o);
		cmatrix_print("x", *x, n, o);
		(void)printf("rank %d\n", rank);
		(void)fflush(stdout);
	    }

	    /*
	     * Verify A X == B.
	     */
	    for (int k = 0; k < o; ++k) {
		for (int i = 0; i < n; ++i) {
		    double complex s = 0.0;

		    for (int j = 0; j < n; ++j) {
			s += a[i][j] * x[j][k];
		    }
		    if (cabs(s - b[i][k]) > EPS) {
			if (opt_a) {
			    assert(!"data miscompare");
			}
			result = T_FAIL;
			goto out;
		    }
		}
	    }

	    /*
	     * Check rank.
	     */
	    if (rank != n) {
		if (opt_a) {
		    assert(!"incorrect rank");
		}
		result = T_FAIL;
		goto out;
	    }
	}

	/*
	 * Test more columns than rows (underdetermined case).
	 */
	for (int m = 1; m <= 4; ++m) {
	    for (int n = m + 1; n <= 5; ++n) {
		for (int o = 1; o <= 2; ++o) {
		    double complex a[m][n];
		    double complex b[m][o];
		    double complex x[n][o];
		    int rank;

		    /*
		     * If -v, print the test header.
		     */
		    if (opt_v) {
			(void)printf("Test vnacommon_qrsolve: trial %3d "
				"A size %d x %d, B size %d x %d\n",
				trial, m, n, n, o);
			(void)fflush(stdout);
		    }

		    /*
		     * Generate random matrices A and B, and solve for X.
		     */
		    rank = qrsolve_helper(&x[0][0], &a[0][0], &b[0][0], m, n, o);
		    if (opt_v) {
			cmatrix_print("a", *a, m, n);
			cmatrix_print("b", *b, n, o);
			cmatrix_print("x", *x, n, o);
			(void)printf("rank %d\n", rank);
			(void)fflush(stdout);
		    }

		    /*
		     * Verify A X == B.
		     */
		    for (int k = 0; k < o; ++k) {
			for (int i = 0; i < m; ++i) {
			    double complex s = 0.0;

			    for (int j = 0; j < n; ++j) {
				s += a[i][j] * x[j][k];
			    }
			    if (cabs(s - b[i][k]) > EPS) {
				if (opt_a) {
				    assert(!"data miscompare");
				}
				result = T_FAIL;
				goto out;
			    }
			}
		    }

		    /*
		     * Check rank.
		     */
		    if (rank != m) {
			if (opt_a) {
			    assert(!"incorrect rank");
			}
			result = T_FAIL;
			goto out;
		    }
		}
	    }
	}

	/*
	 * Test more rows than columns (overdetermined case).
	 */
	for (int n = 1; n <= 4; ++n) {
	    for (int m = n + 1; m <= 5; ++m) {
		for (int o = 1; o <= 2; ++o) {
		    double complex a[m][n];
		    double complex b[m][o];
		    double complex x[n][o];
		    int rank;
		    double error0;

		    /*
		     * If -v, print the test header.
		     */
		    if (opt_v) {
			(void)printf("Test vnacommon_qrsolve: trial %3d "
				"A size %d x %d, B size %d x %d\n",
				trial, m, n, n, o);
			(void)fflush(stdout);
		    }

		    /*
		     * Generate random matrices A and B, and solve for X.
		     */
		    rank = qrsolve_helper(&x[0][0], &a[0][0], &b[0][0], m, n, o);
		    if (opt_v) {
			cmatrix_print("a", *a, m, n);
			cmatrix_print("b", *b, m, o);
			cmatrix_print("x", *x, n, o);
			(void)printf("rank %d\n", rank);
			(void)fflush(stdout);
		    }

		    /*
		     * Get the squared error of the result, then perturb
		     * each X_{j,k} value and verify that the error
		     * doesn't decrease when moving away from the result.
		     */
		    error0 = find_axb_error(&a[0][0], &x[0][0], &b[0][0],
			    m, n, o);
		    for (int k = 0; k < o; ++k) {
			for (int j = 0; j < n; ++j) {
			    static const double complex deltas[] =
				{ 0.001, 0.001 * I, -0.001, -0.001 * I };
			    double complex x0 = x[j][k];

			    for (int i = 0; i < 4; ++i) {
				double e;

				x[j][k] = x0 + deltas[i];
				e = find_axb_error(&a[0][0], &x[0][0], &b[0][0],
					m, n, o);
				if (e < error0) {
				    if (opt_a) {
					assert(!"bad result");
				    }
				    result = T_FAIL;
				    goto out;
				}
			    }
			    x[j][k] = x0;	/* restore x[j][k] */
			}
		    }

		    /*
		     * Check rank.
		     */
		    if (rank != n) {
			if (opt_a) {
			    assert(!"incorrect rank");
			}
			result = T_FAIL;
			goto out;
		    }
		}
	    }
	}
    }
    result = T_PASS;

out:
    report_test_result("QR Solve", result);
}

/*
 * qrsolve_q_helper: generate a random system of equations and solve
 */
static int qrsolve_q_helper(double complex *x, double complex *a,
	double complex *b, double complex *q, int m, int n, int o)
{
    double complex u[m][n];
    double complex v[m][o];

    /*
     * Generate random matrices A and B, and make copies
     * in U and V, respectively.
     */
    for (int i = 0; i < m; ++i) {
	for (int j = 0; j < n; ++j) {
	    a[i * n + j] = u[i][j] = crandn();
	}
	for (int k = 0; k < o; ++k) {
	    b[i * o + k] = v[i][k] = crandn();
	}
    }

    /*
     * Solve the system.  This call destroys both u and v.
     */
    return _vnacommon_qrsolve_q(x, &u[0][0], &v[0][0], q, m, n, o);
}

/*
 * test_vnacommon_qrsolve_q: test QR decomposition with Q
 */
static void test_vnacommon_qrsolve_q()
{
    test_result_type result = T_SKIPPED;

    for (int trial = 1; trial <= N_MATRIX_TRIALS; ++trial) {
	for (int m = 1; m <= 5; ++m) {
	    for (int n = 1; n <= 5; ++n) {
		for (int o = 1; o <= 2; ++o) {
		    double complex a[m][n];
		    double complex b[m][o];
		    double complex x[n][o];
		    double complex q[m][m];
		    int diagonals = MIN(m, n);
		    int rank;

		    /*
		     * If -v, print the test header.
		     */
		    if (opt_v) {
			(void)printf("Test vnacommon_qrsolve_q: trial %3d "
				"size %d x %d\n", trial, m, n);
			(void)fflush(stdout);
		    }

		    /*
		     * Generate random matrices A and B, and solve for X.
		     */
		    rank = qrsolve_q_helper(*x, *a, *b, *q, m, n, o);
		    if (opt_v) {
			cmatrix_print("a", *a, m, n);
			cmatrix_print("b", *b, m, o);
			cmatrix_print("x", *x, n, o);
			cmatrix_print("q", *q, m, m);
			(void)printf("rank %d\n", rank);
			(void)fflush(stdout);
		    }

		    /*
		     * If m <= n, verify A X == B.  Otherwise, the system
		     * is overdetermined and the equality won't hold.  We
		     * checked that case already anyway.
		     */
		    if (m <= n) {
			for (int k = 0; k < o; ++k) {
			    for (int i = 0; i < m; ++i) {
				double complex s = 0.0;

				for (int j = 0; j < n; ++j) {
				    s += a[i][j] * x[j][k];
				}
				if (cabs(s - b[i][k]) > EPS) {
				    if (opt_a) {
					assert(!"data miscompare");
				    }
				    result = T_FAIL;
				    goto out;
				}
			    }
			}
		    }

		    /*
		     * Check rank.
		     */
		    if (rank != diagonals) {
			if (opt_a) {
			    assert(!"incorrect rank");
			}
			result = T_FAIL;
			goto out;
		    }

#if 0
		    /*
		     * Test that Q Q' is the identity matrix.
		     */
		    for (int i = 0; i < m; ++i) {
			for (int j = 0; j < m; ++j) {
			    double complex s = 0.0;

			    for (int k = 0; k < m; ++k) {
				s += q[i][k] * conj(q[j][k]);
			    }
			    s -= (i == j) ? 1.0 : 0.0;
			    if (cabs(s) > EPS) {
				if (opt_a) {
				    assert(!"data miscompare");
				}
				result = T_FAIL;
				goto out;
			    }
			}
		    }
#endif
		}
	    }
	}
    }
    result = T_PASS;

out:
    report_test_result("QR Solve Q", result);
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
 * main: test program
 */
int
main(int argc, char **argv)
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
	    opt_v = true;
	    continue;

	case -1:
	    break;

	default:
	    print_usage();
	    /*NOTREACHED*/
	}
	break;
    }
    argc -= optind;
    argv += optind;
    if (argc != 0) {
	print_usage();
	/*NOTREACHED*/
    }

    test_vnacommon_lu();
    test_vnacommon_mldivide();
    test_vnacommon_mrdivide();
    test_vnacommon_minverse();
    test_vnacommon_qrd();
    test_vnacommon_qrsolve();
    test_vnacommon_qrsolve_q();

    exit(fail_count != 0);
    /*NOTREACHED*/
}
