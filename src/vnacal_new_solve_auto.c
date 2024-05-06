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

#ifdef DEBUG
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

	    (void)printf(" %+.6f%+.6fj", creal(v), cimag(v));
	}
	(void)printf("\n");
    }
    (void)printf("]\n");
}
#endif

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
 * are consistent with the given linear model and error model.
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

    /* best error parameters */
    double complex best_x_vector[x_length];

    /* best unknown parameters */
    double complex best_p_vector[p_length];

    /* lowest seen sum of squares in k_vector */
    double best_sum_k_squared = INFINITY;

    /* scales the Marquardt parameter */
    double marquardt_multiplier = 1.0;

    /* count of equations in the linear error term system */
    int equations = 0;

    /* count of "excess" equations used to solve for the unknown standards */
    int p_equations;

    /* count of rows in the Jacobian matrix */
    int j_rows;

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
     * If a measurement error vector was given, calculate weights
     * for each measurement.
     */
    if (vnp->vn_m_error_vector != NULL) {
	if ((w_vector = vs_calc_weights(vnssp)) == NULL) {
	    goto out;
	}
    }

    /*
     * Init best_x_vector.
     */
    _vnacal_new_solve_init_x_vector(vnssp, best_x_vector, x_length);

#ifdef DEBUG
    (void)printf("p = [\n");
    for (int i = 0; i < p_length; ++i) {
	(void)printf("  %9.6f%+9.6fj\n",
		creal(vnssp->vnss_p_vector[i][findex]),
		cimag(vnssp->vnss_p_vector[i][findex]));
    }
    (void)printf("]\n\n");
#endif /* DEBUG */

    /*
     * Iterate using Levenberg-Marquardt to find the unknown parameters,
     * vnss_p_vector.
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

	/* Jacobian matrix */
	double complex j_matrix[j_rows][p_length];

	/* right-hand side residual vector */
	double complex k_vector[j_rows];

	/* difference vector */
	double complex d_vector[p_length];

	/* current equation index */
	int equation;

	/* rank of a_matrix */
	int rank;

	/* true if current solution is best so far */
	bool best;

	/* Marquardt parameter */
	double lambda = 0.0;

	/* sum of squared mangnitudes of the elements of k_vector */
	double sum_k_squared = 0.0;

	/* root mean squared difference between x_vector and best_x_vector */
	double rms_delta_x = 0.0;

	/* root mean squared difference to be applied to p_vector */
	double rms_delta_p = 0.0;

#if DEBUG >= 3
	/* Jacobian matrix with respect to vnss_p_vector */
	double complex aprimex_matrix[equations][p_length];
#endif

#if DEBUG
	(void)printf("# iteration %d\n", iteration);
#endif /* DEBUG */

	/*
	 * If this is not the first iteration, update the V matrices
	 * from best_x_vector.
	 */
	if (iteration > 0) {
	    if (vs_update_all_v_matrices("vnacal_new_solve",
			vnssp, best_x_vector, x_length) == -1) {
		goto out;
	    }
	}
#ifdef DEBUG
	if (vs_have_v(vnssp)) {
	    int standard = 0;
	    int v_rows, v_columns;
	    vnacal_new_measurement_t *vnmp;

	    if (VL_IS_T(vlp)) {
		v_rows    = VL_S_COLUMNS(vlp);
		v_columns = VL_M_COLUMNS(vlp);
	    } else {
		v_rows    = VL_M_ROWS(vlp);
		v_columns = VL_S_ROWS(vlp);
	    }
	    for (vnmp = vnp->vn_measurement_list; vnmp != NULL;
		    vnmp = vnmp->vnm_next) {
		vnacal_new_msv_matrices_t *vnmmp;

		vnmmp = &vnssp->vnss_msv_matrices[vnmp->vnm_index];
		for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
		    char name[24];

		    (void)sprintf(name, "v%d_%d", standard + 1, sindex + 1);
		    print_cmatrix(name, vnmmp->vnsm_v_matrices[sindex],
			    v_rows, v_columns);
		}
		++standard;
	    }
	}
#endif /* DEBUG */
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
		    if (vs_have_v(vnssp)) {
			v *= vs_get_v(vnssp);
		    }
		    if (w_vector != NULL) {
			v *= w_vector[equation];
		    }
		    if (xindex == -1) {
			b_vector[equation] += v;
		    } else {
			a_matrix[equation][offset + xindex] += v;
		    }
		}
		++equation;
	    }
	}
	assert(equation == equations);
#if DEBUG >= 2
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
#if DEBUG >= 3
	print_cmatrix("q", &q_matrix[0][0], equations, equations);
	print_cmatrix("r", &r_matrix[0][0], equations, x_length);
#endif /* DEBUG >= 3 */

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
	 * If there are no unknown parameters, we're done.
	 */
	if (p_length == 0)
	    goto success;

	/*
	 * Calculate the root mean squared difference between x_vector
	 * and best_x_vector.
	 */
	{
	    double sum = 0.0;

	    for (int i = 0; i < x_length; ++i) {
		sum += _vnacommon_cabs2(x_vector[i] - best_x_vector[i]);
	    }
	    rms_delta_x = sqrt(sum / (double)x_length);
	}

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
	 * Using this method, we make an initial guess for p, solve x as a
	 * linear system, project the remaining equations into a new space
	 * that lets us construct the Jacobian matrix in terms of p only,
	 * use Levenberg-Marquardt to improve our estimate of p and repeat
	 * from the solve for x step until we have suitable convergence.
	 *
	 * The following comments describe the variable projection method.
	 *
	 * Our goal is to minimize the system A(p) x = b in a least-squares
	 * sense, where A(p) is matrix valued function of vector p, b is a
	 * known vector, and x and p are the unknown vectors we need to find
	 * in order to to minimize:
	 *
	 *     || A(p) x - b ||^2
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
	 *          A(p) x = b
	 *   Q1(p) R1(p) x = b
	 *         R1(p) x = Q1(p)^H b
	 *               x = R1(p)^-1 Q1(p)^H b
	 *
	 * which minimizes:
	 *
	 *     || A(p) x - b ||^2
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
	 *     || Q2(p)^H b - Q2(p)^H A(p) x ||
	 *
	 * But Q1(p)^H A(p) = R1(p), and Q2(p)^H A(p) = 0, so
	 *
	 *   = || Q1(p)^H b - R1(p) x ||^2
	 *     || Q2(p)^H b - 0       ||
	 *
	 * and because R1(p) x = Q1(p)^H b from above, Q1(p)^H b - R1(p) x = 0
	 *
	 *   = || 0         ||^2
	 *     || Q2(p)^H b ||
	 *
	 * so we simply need to minimize:
	 *
	 *     || Q2(p)^H b ||^2
	 *
	 * We will improve p using Levenberg-Marquardt	We need the
	 * Jacobian matrix for the residuals in the new system with
	 * respect to each p_k, which we'll now work toward.
	 *
	 * In the equations below, a prime (') symbol on a matrix
	 * represents the element by element partial derivative with
	 * respect to p[k].  We'll use the notication, A'(p)_k to
	 * represent the parital derivative of A with respect to p[k].
	 * We'll consider each k separately, one at a time.
	 *
	 * Recall from above that Q2(p)^H A(p) = 0.  If we take the
	 * partial derivative of each side with respect respect to each
	 * p_k, then from the product rule, we get:
	 *
	 *   Q2'(p)^H_k A(p) +  Q2(p)^H A'(p)_k = 0
	 *
	 * Re-arranging:
	 *
	 *   Q2'(p)^H_k A(p) = -Q2(p)^H A'(p)_k
	 *
	 * Using A(p) = Q1(p) R1(p):
	 *
	 *   Q2'(p)^H_k Q1(p) R1(p) = -Q2(p)^H A'(p)_k
	 *
	 * Multiply on the right by R1(p)^-1 Q1(p)^H b:
	 *
	 *   Q2'(p)^H_k Q1(p) Q1(p)^H b = -Q2(p)^H A'(p)_k R1(p)^-1 Q1(p)^H b
	 *
	 * From above, R1(p)^-1 Q1(p)^H b = x:
	 *
	 *   Q2'(p)^H_k Q1(p) Q1(p)^H b = -Q2(p)^H A'(p)_k x
	 *
	 * Note that Q1(p) Q1(p)^H don't cancel in this direction.
	 *
	 * We can easily find A'(p) because it's simply the coefficients
	 * of the elements of A that contain the given p, but we have no
	 * obvious way of finding Q'(p).  However, Kaufman "A variable
	 * projection method for solving separable nonlinear least squares
	 * problems", BIT 15(1975), pp 49-57, suggests the approximation:
	 *
	 *   Q2'(p)^H_k ≈ -Q2(p)^H A'(p)_k A(p)^+
	 *   where A(p)^+ is the pseudoinverse of A(p), or R1(p)^-1 Q(p)^H
	 *
	 * so:
	 *
	 *   Q2'(p)^H_k b ≈ -Q2(p)^H A'(p)_k R1(p)^-1 Q1(p)^H b
	 *
	 * Again, substituting: R1(p)^-1 Q1(p)^H b = x:
	 *
	 *   Q2'(p)^H_k b ≈ -Q2(p)^H A'(p)_k x
	 *
	 * Thus, we form each column, k (and dummy i), of our Jacobian
	 * matrix (j_matrix) from:
	 *
	 *   J(p)_ik ≈ -Q2(p)^H A'(p)_k x
	 *
	 * And the right hand side residual is:
	 *
	 *   k(p) =  Q2(p)^H b
	 *
	 * To find the correction in p, we can solve:
	 *
	 *   J(p) d = k(p)
	 *
	 * Which would be the Gauss-Newton solution.  But Gauss-Newton
	 * may not converge if the initial guesses aren't very close.
	 * Instead, we create a modified system, J1 d = k1, that
	 * introduces the Marquardt parameter.  From here on, we'll
	 * drop the (p) argument from the equations.
	 *
	 *   J1 = J^H J + lambda I
	 *   k1 = J^H k
	 *
	 * There are many suggestions in the literature for how to choose
	 * lambda, some more practical than others.   N. Yamashita
	 * and M. Fukushima, “On the rate of convergence of the
	 * levenberg-marquardt method,” in Topics in Numerical Analysis,
	 * pp. 239–249, Springer, Vienna, AS, USA, 2001, shows that the
	 * choice lamba = ||j||^2 provides quadratic convergence.  We use
	 * a variation on this: lambda = marquardt_multiplier * ||j||^2.
	 * Where marquardt_multiplier is initially 1.  If the system
	 * diverges, then we double marquart_multiplier and try again
	 * until we get a better solution.  When we get a better solution,
	 * we shrink marquardt_multiplier such that it's the greater of
	 * 1 and the previous value scaled by the improvement in ||j||^2.
	 *
	 * Finally, we use LU decomposition to solve:
	 *
	 *   J1 d = k1
	 *
	 * and apply the correction:
	 *
	 *   p -= d
	 *
	 * until the changes in p_vector and x_vector are sufficiently
	 * small.
	 */
	for (int i = 0; i < j_rows; ++i) {
	    for (int j = 0; j < p_length; ++j) {
		j_matrix[i][j] = 0.0;
	    }
	    k_vector[i] = 0.0;
	}
#if DEBUG >= 3
	for (int i = 0; i < equations; ++i) {
	    for (int j = 0; j < p_length; ++j) {
		aprimex_matrix[i][j] = 0.0;
	    }
	}
#endif /* DEBUG >= 3 */
	equation = 0;
	for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
	    int offset = sindex * (vlp->vl_t_terms - 1);

	    vs_start_system(vnssp, sindex);
	    while (vs_next_equation(vnssp)) {
		vnacal_new_equation_t *vnep = vnssp->vnss_vnep;
		vnacal_new_measurement_t *vnmp = vnep->vne_vnmp;

		while (vs_next_term(vnssp)) {
		    const int s_cell = vs_get_s_cell(vnssp);
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
			int unknown = vnprp->vnpr_unknown_index;
			const int xindex = vs_get_xindex(vnssp);
			double complex v = vs_get_negative(vnssp) ? -1.0 : 1.0;

			if (vs_have_m(vnssp)) {
			    v *= vs_get_m(vnssp);
			}
			if (vs_have_v(vnssp)) {
			    v *= vs_get_v(vnssp);
			}
			assert(xindex >= 0);
			if (w_vector != NULL) {
			    v *= w_vector[equation];
			}
			v *= x_vector[offset + xindex];
#if DEBUG >= 3
			aprimex_matrix[equation][unknown] += v;
#endif /* DEBUG >= 3 */
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
#if DEBUG >= 3
	print_cmatrix("aprimex", &aprimex_matrix[0][0], equations, p_length);
#endif /* DEBUG >= 3 */

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
		weight = 1.0 / _vnacal_get_correlated_sigma(vpmrp1, frequency);
		pindex1 = vnprp1->vnpr_unknown_index;
		vnprp2 = vnprp1->vnpr_correlate;
		j_matrix[j_row][pindex1] = weight;
		k_vector[j_row] = weight *
		    vnssp->vnss_p_vector[pindex1][findex];
		if (vnprp2->vnpr_unknown) {
		    int pindex2 = vnprp2->vnpr_unknown_index;
		    j_matrix[j_row][pindex2] = -weight;
		    k_vector[j_row] -=
			weight * vnssp->vnss_p_vector[pindex2][findex];
		} else { /* known parameter value */
		    k_vector[j_row] -= weight *
			_vnacal_get_parameter_value_i(vnprp2->vnpr_parameter,
				frequency);
		}
		++j_row;
	    }
	    assert(j_row == j_rows);
	}
#if DEBUG >= 2
	print_cmatrix("j", &j_matrix[0][0], j_rows, p_length);
	print_cmatrix("k", k_vector, j_rows, 1);
#endif /* DEBUG */

	/*
	 * Calculate the squared magnitude of k_vector.  From it,
	 * calculate the marquardt parameter as marquardt_multiplier
	 * times the square of the norm.
	 */
	sum_k_squared = 0.0;
	for (int i = 0; i < j_rows; ++i) {
	    sum_k_squared += _vnacommon_cabs2(k_vector[i]);
	}
	lambda = marquardt_multiplier * sum_k_squared;
#if DEBUG >= 2
	(void)printf("# marquardt_multiplier %13.6e\n", marquardt_multiplier);
	(void)printf("# lambda               %13.6e\n", lambda);
#endif /* DEBUG */

	/*
	 * Solve the j_matrix, k_vector system with Marquardt parameter
	 * to create d_vector, the Levenburg-Marquardt correction to
	 * vnss_p_vector.
	 */
	{
	    double complex j1_matrix[p_length][p_length];
	    double complex k1_vector[p_length];
	    double complex determinant;

	    /*
	     * Find J1 = (J^H J + lambda I)
	     */
	    for (int i = 0; i < p_length; ++i) {
		for (int j = 0; j < p_length; ++j) {
		    double complex s = 0.0;

		    for (int k = 0; k < j_rows; ++k) {
			s += conj(j_matrix[k][i]) * j_matrix[k][j];
		    }
		    if (i == j) {
			s += lambda;
		    }
		    j1_matrix[i][j] = s;
		}
	    }

	    /*
	     * Find k1 = J^H k
	     */
	    for (int i = 0; i < p_length; ++i) {
		double complex s = 0.0;

		for (int k = 0; k < j_rows; ++k) {
		    s += conj(j_matrix[k][i]) * k_vector[k];
		}
		k1_vector[i] = s;
	    }
#if DEBUG >= 2
	    print_cmatrix("j1", *j1_matrix, p_length, p_length);
	    print_cmatrix("k1", k1_vector, p_length, 1);
#endif /* DEBUG */

	    /*
	     * Find d = J1 \ k1
	     */
	    determinant = _vnacommon_mldivide(d_vector,
		    &j1_matrix[0][0], k1_vector, p_length, 1);
	    if (determinant == 0.0 || !isnormal(cabs(determinant))) {
		_vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			"singular linear system");
		return -1;
	    }
	}
#ifdef DEBUG
	print_cmatrix("d", d_vector, p_length, 1);
#endif /* DEBUG */

	/*
	 * Calculate the root mean squared of d_vector.
	 */
	{
	    double sum = 0.0;

	    for (int i = 0; i < p_length; ++i) {
		sum += _vnacommon_cabs2(d_vector[i]);
	    }
	    rms_delta_p = sqrt(sum / (double)p_length);
	}

	/*
	 * If we have the best solution so far (or the first), save
	 * the solution.
	 */
#ifdef DEBUG
	(void)printf("# sum_k_squared      %13.6e\n", sum_k_squared);
	(void)printf("# best_sum_k_squared %13.6e\n", best_sum_k_squared);
#endif /* DEBUG */
	best = sum_k_squared < best_sum_k_squared;
	if (best) {
#ifdef DEBUG
	    (void)printf("# best\n");
#endif /* DEBUG */

	    /*
	     * Save the best p_vector and x_vector.
	     */
	    for (int i = 0; i < p_length; ++i) {
		best_p_vector[i] = vnssp->vnss_p_vector[i][findex];
	    }
	    (void)memcpy((void *)best_x_vector, (void *)x_vector,
		    x_length * sizeof(double complex));
#if DEBUG >= 2
	    print_cmatrix("best_p_vector", best_p_vector, p_length, 1);
	    print_cmatrix("best_x_vector", best_x_vector, x_length, 1);
#endif /* DEBUG 2 */

	    /*
	     * Reduce the marquardt multiplier based on the improvement
	     * in sum_k_squared.
	     */
	    marquardt_multiplier *= sum_k_squared / best_sum_k_squared;
	    if (marquardt_multiplier < 1.0) {
		marquardt_multiplier = 1.0;
	    }
	    best_sum_k_squared = sum_k_squared;

	/*
	 * If the new solution is worse: we must have over-corrected.
	 * Restore state to the best solution and increase Marquardt
	 * multiplier.
	 */
	} else {
#if DEBUG >= 2
	    (void)printf("# increasing marquardt parameter\n");
#endif

	    /*
	     * Restore best solution.
	     */
	    for (int i = 0; i < p_length; ++i) {
		vnssp->vnss_p_vector[i][findex] = best_p_vector[i];
	    }
	    vs_update_s_matrices(vnssp);
	    (void)memcpy((void *)x_vector, (void *)best_x_vector,
		    x_length * sizeof(double complex));
	    marquardt_multiplier *= 2.0;
	}

	/*
	 * Stop if the changes in x_vector and p_vector are within
	 * tolerance.
	 */
#ifdef DEBUG
	(void)printf("# rms_delta_x        %13.6e\n", rms_delta_x);
	(void)printf("# rms_delta_p        %13.6e\n", rms_delta_p);
	(void)printf("# vn_et_tolerance    %13.6e\n", vnp->vn_et_tolerance);
	(void)printf("# vn_p_tolerance     %13.6e\n", vnp->vn_p_tolerance);
#endif /* DEBUG */
	if (rms_delta_p <= vnp->vn_p_tolerance &&
	    rms_delta_x <= vnp->vn_et_tolerance) {
#ifdef DEBUG
	    (void)printf("# stop: converged (iteration %d)\n", iteration);
#endif /* DEBUG */
	    break;
	}

	/*
	 * If we have the best solution, apply d_vector to vnss_p_vector.
	 * Otherwise, repeat the previous best solution with the increased
	 * Marquardt parameter.
	 */
	if (best) {
	    for (int i = 0; i < p_length; ++i) {
		vnssp->vnss_p_vector[i][findex] -= d_vector[i];
	    }
	    vs_update_s_matrices(vnssp);
#if DEBUG
	    (void)printf("p = [\n");
	    for (int i = 0; i < p_length; ++i) {
		(void)printf("  %9.6f%+9.6fj\n",
			creal(vnssp->vnss_p_vector[i][findex]),
			cimag(vnssp->vnss_p_vector[i][findex]));
	    }
	    (void)printf("]\n\n");
#endif /* DEBUG */
	}

	/*
	 * Limit the number of iterations.
	 */
	if (iteration >= vnp->vn_iteration_limit) {
	    _vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
		    "system failed to converge at %e Hz", frequency);
	    goto out;
	}
    }
    /*
     * The code above breaks the loop only in the converged case.
     * Both x_vector and vnss_p_vector are guaranteed to be up-to-date
     * here because either the last solution found was the best, or we
     * restored the best solution.
     */

success:
    rv = 0;
    /*FALLTHROUGH*/

out:
    free((void *)w_vector);
    return rv;
}
