/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2022 D Scott Guthridge <scott_guthridge@rompromity.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_new_internal.h"

/* #define DEBUG 1 */

/*
 * PHI_INV: inverse of the golden ratio
 * PHI_INV2: inverse of the golden ratio squared
 */
#define PHI_INV		0.61803398874989484820
#define PHI_INV2	0.38196601125010515180

#ifdef DEBUG
/*
 * print_rmatrix: print an m by n serialized real matrix in octave form
 *   @name: name of matrix
 *   @a: pointer to first element of matrix
 *   @m: number of rows
 *   @n: number of columns
 */
static void print_rmatrix(const char *name, double *a, int m, int n)
{
    (void)printf("%s = [\n", name);
    for (int i = 0; i < m; ++i) {
	for (int j = 0; j < n; ++j) {
	    (void)printf(" %9.5f", a[i * n + j]);
	}
	(void)printf("\n");
    }
    (void)printf("]\n");
}

/*
 * print_cmatrix: print an m by n serialized complex matrix in octave form
 *   @name: name of matrix
 *   @a: pointer to first element of matrix
 *   @m: number of rows
 *   @n: number of columns
 */
static void print_cmatrix(const char *name, double complex *a, int m, int n)
{
    (void)printf("%s = [\n", name);
    for (int i = 0; i < m; ++i) {
	for (int j = 0; j < n; ++j) {
	    double complex v = a[i * n + j];

	    (void)printf(" %9.5f%+9.5fj", creal(v), cimag(v));
	}
	(void)printf("\n");
    }
    (void)printf("]\n");
}
#endif

/*
 * calc_weights: compute w_vector from x_vector, p_vector
 *   @vnssp: solve state structure
 *   @x_vector: current error parameter solution
 *   @w_vector: new weight vector
 *
 * TODO: we're not calculating weights according to the more elegant
 * solution given in the Van hamme paper.  We need to do some rework of
 * the way we represent equations in order to do it right. This is kind
 * of close, though.
 */
static void calc_weights(vnacal_new_solve_state_t *vnssp,
	const double complex *x_vector, double *w_vector)
{
    vnacal_new_t *vnp = vnssp->vnss_vnp;
    const int findex = vnssp->vnss_findex;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int m_cells = VL_M_COLUMNS(vlp) * VL_M_ROWS(vlp);
    vnacal_new_m_error_t *m_error_vector = vnp->vn_m_error_vector;
    const double noise = m_error_vector[findex].vnme_noise;
    const double tracking = m_error_vector[findex].vnme_tracking;
    double complex m_weight_vector[m_cells];
    int equation = 0;

    for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
	int offset = sindex * (vlp->vl_t_terms - 1);

	vs_start_system(vnssp, sindex);
	while (vs_next_equation(vnssp)) {
	    vnacal_new_measurement_t *vnmp;
	    vnacal_new_msv_matrices_t *vnmmp;
	    double u = 0.0;

	    vnmp = vnssp->vnss_vnep->vne_vnmp;
	    vnmmp = &vnssp->vnss_msv_matrices[vnmp->vnm_index];
	    for (int m_cell = 0; m_cell < m_cells; ++m_cell) {
		m_weight_vector[m_cell] = 0.0;
	    }
	    while (vs_next_term(vnssp)) {
		int xindex = vs_get_xindex(vnssp);
		int m_cell = vs_get_m_cell(vnssp);

		if (m_cell >= 0) {
		    double complex v = 1.0;

		    if (vs_get_negative(vnssp)) {
			v *= -1.0;
		    }
		    if (vs_have_s(vnssp)) {
			v *= vs_get_s(vnssp);
		    }
		    if (xindex >= 0) {
			v *= x_vector[offset + xindex];
		    } else {
			v *= -1.0;
		    }
		    assert(m_cell < m_cells);
		    m_weight_vector[m_cell] += v;
		}
	    }

	    /*
	     * Calculate the new weight.
	     */
	    for (int m_cell = 0; m_cell < m_cells; ++m_cell) {
		double complex v = m_weight_vector[m_cell];

		if (v != 0.0) {
		    double complex m;
		    double temp;

		    /*
		     * Add the squared error contributed by m_cell.
		     */
		    assert(vnmp->vnm_m_matrix[m_cell] != NULL);
		    m = vnmmp->vnmm_m_matrix[m_cell];
		    temp = noise * noise +
			tracking * tracking * creal(m * conj(m));
		    u += creal(v * conj(v)) * temp;
		}
	    }
	    u = sqrt(u);
	    if (u < noise) {	/* avoid divide by zero */
		u = noise;
	    }
	    w_vector[equation] = 1.0 / u;
	    ++equation;
	}
    }
}

/*
 * _vnacal_new_solve_auto: solve for both error terms and unknown s-parameters
 *   @vnssp: pointer to state structure
 *   @x_vector: vector of unknowns filled in by this function
 *   @x_length: length of x_vector
 *
 * This implementation is based on the algorithm described in H. Van Hamme
 * and M. Vanden Bossche, "Flexible vector network analyzer calibration
 * with accuracy bounds using an 8-term or a 16-term error correction
 * model," in IEEE Transactions on Microwave Theory and Techniques,
 * vol. 42, no. 6, pp. 976-987, June 1994, doi: 10.1109/22.293566.  There
 * are a few differences, however.  For example, instead of calculating
 * the error bounds on the error parameters, we simply test if the data
 * are consistent with the given linear model and error model.  Also,
 * we do the equation weighting wrong instead of using the "V" matrices.
 */
int _vnacal_new_solve_auto(vnacal_new_solve_state_t *vnssp,
	double complex *x_vector, int x_length)
{
    /* pointer to vnacal_new_t structore */
    vnacal_new_t *vnp = vnssp->vnss_vnp;

    /* current frequency index */
    const int findex = vnssp->vnss_findex;

    /* current frequency */
    const double frequency = vnp->vn_frequency_vector[findex];

    /* pointer to vnacal_t structure */
    vnacal_t *vcp = vnp->vn_vcp;

    /* number of unknown (including correlated) parameters */
    const int p_length = vnp->vn_unknown_parameters;

    /* number of correlated parameters */
    const int correlated = vnp->vn_correlated_parameters;

    /* pointer to vnacal_layout_t structure */
    const vnacal_layout_t *vlp = &vnp->vn_layout;

    /* map indicating which of the unknown parameters are correlated type */
    bool is_correlated_vector[p_length];

    /* vector of weights for each measurement */
    double *w_vector = NULL;

    /* weight vector that was used to create best solution */
    double *best_w_vector = NULL;

    /* best error parameters */
    double complex best_x_vector[x_length];

    /* best unknown parameters */
    double complex best_p_vector[p_length];

    /* Gauss-Newton correction vector generated from the best solution */
    double complex best_d_vector[p_length];

    /* sum of squares of best_d_vector, initially infinite */
    double best_sum_d_squared = INFINITY;

    /* count of equations in the linear error term system */
    int equations = 0;

    /* count of "excess" equations used to solve for the unknown standards */
    int p_equations;

    /* count of rows in the Jacobian matrix */
    int j_rows;

    /* current number of iterations in the backtracking line search */
    int backtrack_count = 0;

    /* return status of this function */
    int rv = -1;

    /*
     * Test that we have at least as many equations as unknowns.
     */
    equations = vnp->vn_equations;
    assert(x_length == vnp->vn_systems * (vlp->vl_t_terms - 1));
    if (equations + correlated < x_length + p_length) {
	_vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: not enough "
		"standards given to solve the system");
	return -1;
    }
    p_equations = equations - x_length;
    j_rows = p_equations + correlated;

    /*
     * Determine which unknown parameters are of correlated type.
     */
    (void)memset((void *)is_correlated_vector, 0, sizeof(is_correlated_vector));
    for (vnacal_new_parameter_t *vnprp = vnp->vn_unknown_parameter_list;
	    vnprp != NULL; vnprp = vnprp->vnpr_next_unknown) {
	assert(vnprp->vnpr_unknown_index >= 0);
	assert(vnprp->vnpr_unknown_index < p_length);
	if (vnprp->vnpr_parameter->vpmr_type == VNACAL_CORRELATED) {
	    is_correlated_vector[vnprp->vnpr_unknown_index] = true;
	}
    }

    /*
     * Iterate using Gauss-Newton to find the unknown parameters, vnss_p_vector.
     */
    for (int iteration = 0; /*EMPTY*/; ++iteration) {
	/* coefficient matrix of the linear error term system */
	double complex a_matrix[equations][x_length];

	/* right-hand side of the linear error term system */
	double complex b_vector[equations];

	/* orthogonal matrix from QR decomposition of a_matrix */
	double complex q_matrix[equations][equations];

	/* upper-triangular matrix from QR decompositin of a_matrix */
	double complex r_matrix[equations][x_length];

	/* Jacobian matrix for Gauss-Newton */
	double complex j_matrix[j_rows][p_length];

	/* right-hand side residual vector for Gauss-Newton */
	double complex k_vector[j_rows];

	/* difference vector from Gauss-Newton */
	double complex d_vector[p_length];

	/* current equation index */
	int equation;

	/* rank of a_matrix */
	int rank;

	/* sum of squared magnitudes of the elements of d_vector */
	double sum_d_squared;

#if DEBUG >= 2
	/* Jacobian of a_matrix with respect to vnss_p_vector */
	double complex aprimex_matrix[equations][p_length];
#endif

#if DEBUG >= 2
	(void)printf("# iteration %d\n", iteration);
#endif
	/*
	 * Build a_matrix and right-hand-side b_vector.  This linear
	 * system is used to solve for the error parameters (x_vector).
	 * It's built from the measurements of the calibration standards
	 * added to the vnacal_new_t structure via the vnacal_new_add_*
	 * functions.
	 *
	 * Note that in calibration types other than T16 and U16, the
	 * leakage equations are handled outside of the system and will
	 * have already been subtracted out.  For example, a double
	 * reflect standard in 2x2 T8 contributes only two equations
	 * intead of four.  In TE10 and UE10, the other two are used to
	 * compute leakage terms -- that's done outside of this function.
	 */
	for (int i = 0; i < equations; ++i) {
	    for (int j = 0; j < x_length; ++j) {
		a_matrix[i][j] = 0.0;
	    }
	    b_vector[i] = 0.0;
	}
	equation = 0;
	for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
	    int offset = sindex * (vlp->vl_t_terms - 1);

	    /*
	     * The vs_start_system, vs_next_equation and vs_next_term
	     * functions form an abstract iterator that systematically
	     * walks through the equations added via vnacal_new_add_*.
	     *
	     * In the case of UE14 (used to solve classic E12 SOLT), each
	     * column of the measurement matrix forms an independent
	     * linear system with its own separate error terms.
	     * These independent systems, however, share the same unknown
	     * calibration parameters (vnss_p_vector), and for simplicity
	     * of solving them, we create one big block-diagonal matrix
	     * equation to solve all systems at once.
	     */
	    vs_start_system(vnssp, sindex);
	    while (vs_next_equation(vnssp)) {
		while (vs_next_term(vnssp)) {
		    int xindex = vs_get_xindex(vnssp);
		    double complex v = vs_get_negative(vnssp) ? -1.0 : 1.0;

		    if (vs_have_m(vnssp)) {
			v *= vs_get_m(vnssp);
		    }
		    if (vs_have_s(vnssp)) {
			v *= vs_get_s(vnssp);
		    }
		    if (w_vector != NULL) {
			v *= w_vector[equation];
		    }
		    if (xindex == -1) {
			b_vector[equation] = v;
		    } else {
			a_matrix[equation][offset + xindex] = v;
		    }
		}
		++equation;
	    }
	}
	assert(equation == equations);
#ifdef DEBUG
	print_cmatrix("a", &a_matrix[0][0], equations, x_length);
	print_cmatrix("b", b_vector, equations, 1);
#endif /* DEBUG */

	/*
	 * Find the QR decomposition of a_matrix, creating q_matrix and
	 * r_matrix, destroying a_matrix.
	 *
	 * Conceptually, Q and R are partitioned as follows:
	 *
	 *   [ Q1 Q2 ] [ R1
	 *               0 ]
	 *
	 * with dimensions:
	 *   Q1: equations x x_length
	 *   Q2: equations x (equations - x_length)
	 *   R1: x_length  x x_length
	 */
	rank = _vnacommon_qr(*a_matrix, *q_matrix, *r_matrix,
		equations, x_length);
	if (rank < x_length) {
	    _vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
		    "singular linear system");
	    goto out;
	}
#if DEBUG >= 2
	print_cmatrix("q", &q_matrix[0][0], equations, equations);
	print_cmatrix("r", &r_matrix[0][0], equations, x_length);
#endif /* DEBUG >= 2 */

	/*
	 * Solve for x_vector.
	 *   R x = Q^H b, where Q^H is the conjugate transpose of Q
	 */
	_vnacommon_qrsolve2(x_vector, *q_matrix, *r_matrix, b_vector,
		equations, x_length, 1);
#ifdef DEBUG
	print_cmatrix("x", x_vector, x_length, 1);
#endif /* DEBUG */

	/*
	 * If measurement error was given (via vnacal_new_set_m_error),
	 * then we weight the equations in the system to compensate for
	 * uncertainty in the measurements.
	 *
	 * If we haven't already done so, allocate w_vector and
	 * best_w_vector, compute weights based on the initial guesses
	 * of vnss_p_vector, and restart the loop from the top.
	 */
	if (vnp->vn_m_error_vector != NULL && w_vector == NULL) {
	    if ((w_vector = malloc(equations * sizeof(double))) == NULL) {
		_vnacal_error(vcp, VNAERR_SYSTEM,
			"malloc: %s", strerror(errno));
		goto out;
	    }
	    if (p_length != 0) {
		if ((best_w_vector = calloc(equations,
				sizeof(double))) == NULL) {
		    _vnacal_error(vcp, VNAERR_SYSTEM,
			    "calloc: %s", strerror(errno));
		    goto out;
		}
	    }
	    for (int i = 0; i < equations; ++i) {
		w_vector[i] = 1.0;
	    }
	    calc_weights(vnssp, x_vector, w_vector);
#ifdef DEBUG
	    print_rmatrix("w", w_vector, equations, 1);
#endif /* DEBUG */
	    --iteration;
	    continue;
	}

#ifdef DEBUG
    (void)printf("p = [\n");
    for (int i = 0; i < p_length; ++i) {
	(void)printf("  %f%+fj\n",
		creal(vnssp->vnss_p_vector[i][findex]),
		cimag(vnssp->vnss_p_vector[i][findex]));
    }
    (void)printf("]\n\n");
#endif /* DEBUG */

	/*
	 * If there are no unknown parameters, we're done.
	 */
	if (p_length == 0)
	    goto success;

	/*
	 * At this point, we know that the system is nonlinear.
	 * We have two sets of variables to solve: the error terms,
	 * x_vector, and the unknown calibration parameters, p_vector.
	 * The a_matrix depends on p_vector; consequently, the
	 * system A x = b contains products of p and x variables,
	 * thus is quadratic.  It is, however, a separable nonlinear
	 * least squares problem that can be solved using the variable
	 * projection method as described by Golub and LeVeque, 1979
	 * http://faculty.washington.edu/rjl/pubs/GolubLeVeque1979/
	 * GolubLeVeque1979.pdf
	 *
	 * Using this method, we make an initial guess for p, solve x as
	 * a linear system, project the remaining equations into a new
	 * space that lets us construct the Jacobian in terms of p only,
	 * use Gauss-Newton to improve our estimate of p and repeat from
	 * the solve for x step until we have suitable convergence.
	 *
	 * Following is a brief derivation of the variable projection method.
	 *
	 * Our goal is to minimize the system A(p) x = b in a least-squares
	 * sense, where A(p) is matrix valued function of vector p, b is a
	 * known vector, and x and p are the unknown vectors we need to find
	 * in order to to minimize:
	 *
	 *     || b - A(p) x ||^2
	 *
	 * There must be an orthogonal matrix Q that diagonalizes A to R.
	 * Both of the new resulting matrices still depend on p.
	 *
	 *   A(p) = Q(p) R(p)
	 *
	 * Partition Q(p) and R(p) as follows:
	 *
	 *   A(p) = [ Q1(p) Q2(p) ] [ R1(p) ]
	 *                          [   0   ]
	 *        = Q1(p) R1(p)
	 *
	 * It follows also that:
	 *
	 *   Q2(p)^H A(p) = 0
	 *
	 * where ^H is the conjugate transpose.
	 *
	 * Solve A(p) x = b for x:
	 *
	 *   A(p)        x = b
	 *   Q1(p) R1(p) x = b
	 *   R1(p)       x = Q1(p)^H b
	 *               x = R1(p)^-1 Q1(p)^H b
	 *
	 * which minimizes:
	 *
	 *     || b - A(p) x ||^2
	 *
	 * with our current guess for p.
	 *
	 * From the invariance of the 2-norm under orthogonal
	 * transformations, we can multiply the inside of the
	 * above by Q^H * without changing the norm:
	 *
	 *   = || Q(p)^H (b - A(p) x) ||^2
	 *
	 *   = || Q1(p)^H b - Q1(p)^H A(p) x ||^2
	 *     || Q2(p)^H b - Q2(p)^H A(p) x ||^2
	 *
	 * But Q1(p)^H A(p) = R1(p), and Q2(p)^H A(p) = 0, so
	 *
	 *   = || Q1(p)^H b - R1(p) x ||^2
	 *     || Q2(p)^H b - 0       ||^2
	 *
	 * and because R1(p) x = Q1(p)^H b from above, Q1(p)^H b - R1(p) x = 0
	 *
	 *   = || 0         ||^2
	 *     || Q2(p)^H b ||^2
	 *
	 * so we simply need to minimize:
	 *
	 *     || Q2(p)^H b ||^2
	 *
	 * We will improve p using Gauss-Newton.  We need the Jacobian
	 * matrix of the residuals above with respect to each p_k.
	 *
	 * Recall from above that Q2(p)^H A(p) = 0.  If we take the
	 * partial derivative of each side with respect respect to each
	 * p_k, then from the product rule, we get:
	 *
	 *   Q2'(p)^H A(p) +  Q2(p)^H A'(p) = 0
	 *
	 * Re-arranging:
	 *
	 *   Q2'(p)^H A(p) = -Q2(p)^H A'(p)
	 *
	 * Using A(p) = Q1(p) R1(p):
	 *
	 *   Q2'(p)^H Q1(p) R1(p) = -Q2(p)^H A'(p)
	 *
	 * Multiply on the right by R1(p)^-1 Q1(p)^H b:
	 *
	 *   Q2'(p)^H Q1(p) Q1(p)^H b = -Q2(p)^H A'(p) R1(p)^-1 Q1(p)^H b
	 *
	 * We'd really like Q1(p) Q1(p)^H to cancel, but they don't in
	 * this direction.  But Kaufman 1975 suggests the approximation:
	 *
	 *   Q2'(p)^H ≈ -Q2(p)^H A'(p) A(p)^+
	 *   where A(p)^+ is the pseudoinverse of A(p), or R1(p)^-1 Q(p)^H
	 *
	 * Using the approximation, we can treat Q1(p) Q1(p)^H as if they
	 * do cancel:
	 *
	 *   Q2'(p)^H b ≈ -Q2(p)^H A'(p) R1(p)^-1 Q1(p)^H b
	 *
	 * From above, R1(p)^-1 Q1(p)^H b = x:
	 *
	 *   Q2'(p)^H b ≈ -Q2(p)^H A'(p) x
	 *
	 * We can easily find A'(p) since it's just the coefficients
	 * of A that contain the given p.  Thus our Jacobian matrix
	 * (j_matrix) is:
	 *
	 *   J(p) ≈ -Q2(p)^H A'(p) x
	 *
	 * And the right hand side residual for Gauss-Newton is:
	 *
	 *   k(p) =  Q2(p)^H b
	 *
	 * To find the correction in p, we use QR decomposition to solve:
	 *
	 *   J(p) d = k(p)
	 *
	 * and apply the correction:
	 *
	 *   p -= d
	 */
	for (int i = 0; i < j_rows; ++i) {
	    for (int j = 0; j < p_length; ++j) {
		j_matrix[i][j] = 0.0;
	    }
	    k_vector[i] = 0.0;
	}
#if DEBUG >= 2
	for (int i = 0; i < equations; ++i) {
	    for (int j = 0; j < p_length; ++j) {
		aprimex_matrix[i][j] = 0.0;
	    }
	}
#endif /* DEBUG >= 2 */
	equation = 0;
	for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
	    int offset = sindex * (vlp->vl_t_terms - 1);

	    vs_start_system(vnssp, sindex);
	    while (vs_next_equation(vnssp)) {
		vnacal_new_equation_t *vnep = vnssp->vnss_vnep;
		vnacal_new_measurement_t *vnmp = vnep->vne_vnmp;

		while (vs_next_term(vnssp)) {
		    int xindex = vs_get_xindex(vnssp);
		    int s_cell = vs_get_s_cell(vnssp);
		    vnacal_new_parameter_t *vnprp = NULL;

		    /*
		     * Apply this term's contribution to the current
		     * row of the Jacobian matrix.  We're computing
		     * -Q2(p)^H A'(p) x, but doing the the first matrix
		     * multiplication with loop nesting inverted from the
		     * usual order so that we can go row by row through A.
		     */
		    if (s_cell >= 0 && (vnprp =
				vnmp->vnm_s_matrix[s_cell])->vnpr_unknown) {
			double complex v = vs_get_negative(vnssp) ? -1.0 : 1.0;
			int unknown = vnprp->vnpr_unknown_index;

			if (vs_have_m(vnssp)) {
			    v *= vs_get_m(vnssp);
			}
			if (w_vector != NULL) {
			    v *= w_vector[equation];
			}
			assert(xindex >= 0);
			v *= x_vector[offset + xindex];
#if DEBUG >= 2
			aprimex_matrix[equation][unknown] += v;
#endif /* DEBUG >= 2 */
			for (int k = 0; k < p_equations; ++k) {
			    j_matrix[k][unknown] -=
				conj(q_matrix[equation][x_length + k]) * v;
			}
		    }
		}

		/*
		 * Build the right-hand-side vector of residuals, k_vector:
		 *     k(p) = Q2(p)^H b
		 */
		for (int k = 0; k < p_equations; ++k) {
		    k_vector[k] += conj(q_matrix[equation][x_length + k]) *
			b_vector[equation];
		}
		++equation;
	    }
	}
	assert(equation == equations);
#if DEBUG >= 2
	print_cmatrix("aprimex", &aprimex_matrix[0][0], equations, p_length);
#endif /* DEBUG >= 2 */

	/*
	 * Add an additional row to j_matrix and k_vector for each
	 * correlated parameter.
	 */
	if (correlated != 0) {
	    int j_row = p_equations;
	    vnacal_new_parameter_t *vnprp1;

	    for (vnprp1 = vnp->vn_unknown_parameter_list; vnprp1 != NULL;
		    vnprp1 = vnprp1->vnpr_next_unknown) {
		vnacal_parameter_t *vpmrp1 = vnprp1->vnpr_parameter;
		vnacal_new_parameter_t *vnprp2;
		int pindex1;
		double weight;

		/*
		 * Skip if not a correlated parameter.
		 */
		if (vpmrp1->vpmr_type != VNACAL_CORRELATED)
		    continue;

		/*
		 * When the correlated parameter is correlated with
		 * another unknown parameter, we can describe it with
		 * an equation of the form:
		 *
		 *   weight p[i] - weight p[j] = 0
		 *
		 * where weight is one over the sigmal value (standard
		 * deviation) associated with the correlated parameter.
		 * When the correlated parameter is correlated with
		 * a constant parameter, we can describe it with an
		 * equation of the form:
		 *
		 *   weight p[i] = weight K
		 *
		 * We represent these equations using a matrix, E, and
		 * column vector, f, such that:
		 *
		 *   E p = f
		 *
		 * In the first case, we store weight and -weight into the
		 * columns of E corresponding to p[i] and p[j], with zero
		 * in the corresponding row of f, thus setting the two
		 * parameters equal under the weight.  In the second case,
		 * we store weight into the column of E corresponding
		 * to p[i], and the constant parameter info f.
		 *
		 * In the J k system, however, we're not computing p,
		 * but rather the error in p0 that leads us to a better
		 * prediction, p1:
		 *
		 *   E d = E p0 - f
		 *   p1 = p0 - d
		 *
		 * Thus we store E into the lower rows of j_matrix and
		 * (E p0 - f) into the lower rows of k_vector.	We do
		 * the mulplication E*p0 by row.
		 */
		weight = 1.0 / _vnacal_get_correlated_sigma(vpmrp1,
			frequency);
		vnprp2 = vnprp1->vnpr_correlate;
		pindex1 = vnprp1->vnpr_unknown_index;
		j_matrix[j_row][pindex1] = weight;	/* partial derivative */
		k_vector[j_row] += weight *
		    vnssp->vnss_p_vector[pindex1][findex];/*contr. to residual*/
		if (vnprp2->vnpr_unknown) {
		    int pindex2 = vnprp2->vnpr_unknown_index;
		    j_matrix[j_row][pindex2] = -weight; /* partial drvtv */
		    k_vector[j_row] -= weight *
			vnssp->vnss_p_vector[pindex2][findex];   /* resid */
		} else { /* known parameter value */
		    k_vector[j_row] -= weight *	/* contr. to resid */
			_vnacal_get_parameter_value_i(vnprp2->vnpr_parameter,
				frequency);
		}
		++j_row;
	    }
	    assert(j_row == j_rows);
	}
#ifdef DEBUG
	print_cmatrix("j", &j_matrix[0][0], j_rows, p_length);
	print_cmatrix("k", k_vector, j_rows, 1);
#endif /* DEBUG */

	/*
	 * Solve the j_matrix, k_vector system to create d_vector, the
	 * Gauss-Newton correction to vnss_p_vector.
	 */
	if (j_rows == p_length) {
	    double complex determinant;

	    determinant = _vnacommon_mldivide(d_vector,
		    &j_matrix[0][0], k_vector, p_length, 1);
	    if (determinant == 0.0 || !isnormal(cabs(determinant))) {
		_vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			"singular linear system");
		return -1;
	    }
	} else {
	    int rank;

	    rank = _vnacommon_qrsolve(d_vector, &j_matrix[0][0],
		    k_vector, j_rows, p_length, 1);
	    if (rank < p_length) {
		_vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			"singular linear system");
		return -1;
	    }
	}
#ifdef DEBUG
	(void)printf("# findex %d iteration %d\n", findex, iteration);
	print_cmatrix("d", d_vector, p_length, 1);
#endif /* DEBUG */

	/*
	 * Calculate the squared magnitude of d_vector.
	 */
	sum_d_squared = 0.0;
	for (int i = 0; i < p_length; ++i) {
	    sum_d_squared += creal(d_vector[i] * conj(d_vector[i]));
	}
#ifdef DEBUG
	(void)printf("# sum_d_squared      %13.6e\n", sum_d_squared);
	(void)printf("# best_sum_d_squared %13.6e\n", best_sum_d_squared);
#endif /* DEBUG */

	/*
	 * If the error is within the target tolerance, stop.
	 */
	if (sum_d_squared / (double)p_length <= vnp->vn_p_tolerance *
						vnp->vn_p_tolerance) {
#ifdef DEBUG
	    (void)printf("# stop: converged\n");
#endif /* DEBUG */
	    break;
	}

	/*
	 * If we have the best solution so far (or the first), remember
	 * this solution and apply the new correction to vnss_p_vector.
	 */
	if (sum_d_squared < best_sum_d_squared) {
	    double sum_p_squared = 0.0;

#ifdef DEBUG
	    (void)printf("# best\n");
#endif /* DEBUG */
	    /*
	     * Limit the magnitude of d_vector to keep it smaller than
	     * the magnitude of vnss_p_vector (or smaller than one if
	     * vnss_p_vector is less than one).  This improves stability
	     * at the cost of slowing convergence.  It makes it less
	     * likely that we jump into an adjacent basin of attraction.
	     * We use one over the golden ratio as the maximum norm of
	     * d_vector relative to vnss_p_vector (or one).
	     */
	    for (int i = 0; i < p_length; ++i) {
		double complex p = vnssp->vnss_p_vector[i][findex];

		sum_p_squared += creal(p * conj(p));
	    }
	    if (sum_p_squared < 1.0) {	/* if less than 1, make it 1 */
		sum_p_squared = 1.0;
	    }
	    if (sum_d_squared > sum_p_squared * PHI_INV2) {
		double scale = sqrt(sum_p_squared / sum_d_squared) * PHI_INV;

#ifdef DEBUG
		(void)printf("# scaling d_vector by: %f\n", scale);
#endif /* DEBUG */
		for (int i = 0; i < p_length; ++i) {
		    d_vector[i] *= scale;
		}
	    }

	    /*
	     * Remember this solution.
	     */
	    (void)memcpy((void *)best_x_vector, (void *)x_vector,
		    x_length * sizeof(double complex));
	    for (int i = 0; i < p_length; ++i) {
		best_p_vector[i] = vnssp->vnss_p_vector[i][findex];
		best_d_vector[i] = d_vector[i];
	    }
	    if (w_vector != NULL) {
		(void)memcpy((void *)best_w_vector, (void *)w_vector,
			equations * sizeof(double));
	    }
	    best_sum_d_squared = sum_d_squared;

	    /*
	     * Update the weight vector.
	     */
	    if (w_vector != NULL) {
		calc_weights(vnssp, x_vector, w_vector);
#ifdef DEBUG
		for (int i = 0; i < equations; ++i) {
		    (void)printf("# w[%2d] = %f\n", i, w_vector[i]);
		}
#endif /* DEBUG */
	    }

	    /*
	     * Apply d_vector to vnss_p_vector.
	     */
	    for (int i = 0; i < p_length; ++i) {
		vnssp->vnss_p_vector[i][findex] -= d_vector[i];
	    }
#ifdef DEBUG
	    (void)printf("p = [\n");
	    for (int i = 0; i < p_length; ++i) {
		(void)printf("  %f%+fj\n",
			creal(vnssp->vnss_p_vector[i][findex]),
			cimag(vnssp->vnss_p_vector[i][findex]));
	    }
	    (void)printf("]\n\n");
#endif /* DEBUG */
	    backtrack_count = 0;

	/*
	 * If the new solution is worse: we must have over-corrected.
	 * Use a backtracking line search that keeps dividing d_vector
	 * in half and retrying from the best solution.
	 */
	} else {
	    if (++backtrack_count > 6) {
#ifdef DEBUG
	    (void)printf("# stop: backtrack count\n");
#endif /* DEBUG */
		break;
	    }
#ifdef DEBUG
	    (void)printf("# retry with half d_vector\n");
#endif /* DEBUG */
	    for (int i = 0; i < p_length; ++i) {
		best_d_vector[i] *= 0.5;
#ifdef DEBUG
	    (void)printf("# half-d[%2d] = %13.6e %+13.6ej\n", i,
		    creal(best_d_vector[i]),
		    cimag(best_d_vector[i]));
#endif /* DEBUG */
		vnssp->vnss_p_vector[i][findex] =
		    best_p_vector[i] + best_d_vector[i];
	    }
#ifdef DEBUG
	    for (int i = 0; i < p_length; ++i) {
		(void)printf("# p[%2d] = %13.6e %+13.6ej\n", i,
			creal(vnssp->vnss_p_vector[i][findex]),
			cimag(vnssp->vnss_p_vector[i][findex]));
	    }
#endif /* DEBUG */
	    if (w_vector != NULL) {
		(void)memcpy((void *)w_vector, (void *)best_w_vector,
			equations * sizeof(double));
		calc_weights(vnssp, x_vector, w_vector);
#ifdef DEBUG
		print_rmatrix("w", w_vector, equations, 1);
#endif /* DEBUG */
	    }
	}
	vs_update_s_matrices(vnssp);

	/*
	 * Limit the number of iterations.
	 *
	 * TODO: instead of failing here, just return what we have so
	 * far and let calc_rms_error check if it's close enough.
	 */
	if (iteration >= 50) {
	    _vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
		    "system failed to converge at %e Hz", frequency);
	    goto out;
	}
    }

    /*
     * Load the best solution.
     */
    (void)memcpy((void *)x_vector, (void *)best_x_vector,
	    x_length * sizeof(double complex));
    for (int i = 0; i < p_length; ++i) {
	vnssp->vnss_p_vector[i][findex] = best_p_vector[i];
    }
    vs_update_s_matrices(vnssp);

success:
    rv = 0;
    /*FALLTHROUGH*/

out:
    free((void *)w_vector);
    free((void *)best_w_vector);
    return rv;
}
