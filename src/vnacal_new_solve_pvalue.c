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
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_new_internal.h"


/*
 * chisq_pvalue: return 1 - CDF of the chi-squared distribution with k d.f.
 *   @n:  degrees of freedom
 *   @x2: chi-squared statistic
 *
 *   Returns the probability of finding a chi-squared statistic equal
 *   to or greater than x2 given that x2 is distributed according to
 *   a chi-squared distribution with k degrees of freedom.
 *
 *   The implementation follows a much simplified version of the method
 *   described in "Computation of the Incomplete Gamma Function Ratios
 *   and their Inverse", Didonato, Morris, 1986. Instead of selecting
 *   one of ten different methods for computing the regularized upper
 *   incomplete gamma function depending on the values of the parameters,
 *   we use use only erfc for the special case of n == 1, and the two
 *   finite sums (equations 14) for everything else.  While the paper
 *   suggests different methods for some parameter combinations such
 *   as very large n, this simplified version achieves a maximum error
 *   of about 10^-16 over the range of values we expect to be used in
 *   this application.  Further, it doesn't require computation of the
 *   complete gamma function.
 */
static double chisq_pvalue(int n, double x2)
{
    double x = x2 / 2.0;
    double result;

    assert(n >= 1);

    /*
     * For zero x, the result is 1.
     */
    if (x <= 0.0) {
	result = 1.0;

    /*
     * For the special case of one degree of freedom, use erfc.
     */
    } else if (n == 1) {
	result = erfc(sqrt(x));

    /*
     * For n even,
     *     Q(1, x) = exp(-x)
     *     Q(a + 1, x) = Q(a, x) + R(a, x) / a
     *
     *     where:
     *         Q is one minus the regularized upper incomplete gamma function,
     *         a = n/2, and
     *         R(a, x) = e^(-x) x^a / Gamma(a), which results from an
     *         integration by parts on the defintion of the regularized
     *         upper incomplete gamma function.
     */
    } else if ((n & 1) == 0) {
	double c = exp(-x);
	double f = 1.0;
	double s = 0.0;

	n >>= 1;
	for (int i = 0; i < n; ++i) {
	    if (i != 0) {
		f *= x / (double)i;
	    }
	    s += f;
	}
	result = c * s;

    /*
     * For n odd,
     *     Q(1/2, x) = erfc(sqrt(x))
     *     Q(a + 1, x) = Q(a, x) + R(a, x) / a
     *
     *     with the same conditions as above
     */
    } else {
	double c1 = erfc(sqrt(x));
	double c2 = exp(-x) / sqrt(M_PI * x);
	double f = 1.0;
	double s = 0.0;

	n >>= 1;
	for (int i = 1; i <= n; ++i) {
	    f *= x / (i - 0.5);
	    s += f;
	}
	result = c1 + c2 * s;
    }

    return result;
}

/*
 * _vnacal_new_solve_calc_pvalue: calculate probability system is consistent
 *   @vnssp: solve state structure
 *   @x_vector: vector of solved error terms
 *   @x_length: length of x_vector
 */
double _vnacal_new_solve_calc_pvalue(vnacal_new_solve_state_t *vnssp,
	const double complex *x_vector, int x_length)
{
    vnacal_new_t *vnp = vnssp->vnss_vnp;
    const int findex = vnssp->vnss_findex;
    vnacal_new_leakage_term_t **leakage_matrix = vnssp->vnss_leakage_matrix;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int m_rows    = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    const vnacal_new_m_error_t *m_error_vector = vnp->vn_m_error_vector;
    double noise, tracking;
    double chisq = 0.0;
    int df = 0;

    /*
     * Get the expected measurement error.
     */
    assert(m_error_vector != NULL);
    noise = m_error_vector[findex].vnme_noise;
    tracking = m_error_vector[findex].vnme_tracking;

    /*
     * Accumulate the squared magnitudes of the residuals of the
     * linear system, all normalized to one standard deviation.
     */
    for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
	int offset = sindex * (vlp->vl_t_terms - 1);

	vs_start_system(vnssp, sindex);
	while (vs_next_equation(vnssp)) {
	    vnacal_new_equation_t *vnep = vnssp->vnss_vnep;
	    vnacal_new_measurement_t *vnmp = vnep->vne_vnmp;
	    vnacal_new_msv_matrices_t *vnsmp =
		&vnssp->vnss_msv_matrices[vnmp->vnm_index];
	    const int eq_row     = vnep->vne_row;
	    const int eq_column  = vnep->vne_column;
	    const int eq_cell    = eq_row * m_columns + eq_column;
	    double complex residual = 0.0;
	    double complex m_value;
	    double squared_residual;
	    double divisor;

	    while (vs_next_term(vnssp)) {
		double complex value = vs_get_negative(vnssp) ? -1.0 : 1.0;
		int xindex = vs_get_xindex(vnssp);

		if (vs_have_m(vnssp)) {
		    value *= vs_get_m(vnssp);
		}
		if (vs_have_s(vnssp)) {
		    value *= vs_get_s(vnssp);
		}
		if (vs_have_v(vnssp)) {
		    value *= vs_get_v(vnssp);
		}
		if (xindex >= 0) {
		    assert(offset + xindex < x_length);
		    value *= x_vector[offset + xindex];
		} else {
		    value = -value;
		}
		residual += value;
	    }
	    squared_residual = creal(residual * conj(residual));

	    /*
	     * Normlize the residual to 1 standard deviation.
	     */
	    m_value = vnsmp->vnsm_m_matrix[eq_cell];
	    divisor = creal(m_value * conj(m_value));
	    divisor *= tracking * tracking;
	    divisor += noise * noise;
	    squared_residual /= divisor;

	    /*
	     * Because the residuals are complex, each contributes
	     * two degress of freedom.  It's also necessary to multiply
	     * the squared residual by 2 because the complex residual
	     * normalized to 1 standard deviation is really a real and
	     * an imaginary part, each with only 1 / sqrt(2) standard
	     * deviations.  Normalize the components to 1.
	     */
	    chisq += 2.0 * squared_residual;
	    df += 2;
	}
	/*
	 * Subtract two degrees of freedom for each dependent complex
	 * variable.
	 */
	df -= 2 * (vlp->vl_t_terms - 1);
    }

    /*
     * Accumulate variance from leakage parameters outside of
     * the linear system.
     */
    if (leakage_matrix != NULL) {
	for (int row = 0; row < m_rows; ++row) {
	    for (int column = 0; column < m_columns; ++column) {
		if (row != column) {
		    const int m_cell = row * m_columns + column;
		    const vnacal_new_leakage_term_t *ltp =
			leakage_matrix[m_cell];
		    double value;

		    if (ltp->vnlt_count > 1) {
			double complex sum_x = ltp->vnlt_sum;
			const int n = ltp->vnlt_count;
			double n_mean_squared;
			double weight;

			n_mean_squared = creal(sum_x * conj(sum_x)) / n;
			value = ltp->vnlt_sumsq - n_mean_squared;
			weight = 1.0 / (noise * noise +
				n_mean_squared / n * tracking * tracking);
			value *= weight;
			chisq += 2.0 * value;
			df += 2 * (n - 1);
		    }
		}
	    }
	}
    }

    /*
     * Note that we don't collect residuals from correlated parameters.
     * The reason is that these are already accounted for in the linear
     * system.  If we know the p's, we can find the x's by solution of
     * a linear system, and while it's less obvious, the converse is
     * also true: given the x's, we can find the p's as a linear system.
     * Thus the x's and p's are dependent.  Residuals from the correlated
     * parameters apply pressure on the p values which are then reflected
     * in the x values.
     */

    /*
     * If there are no degrees of freedom, then the p-value is zero.
     */
    if (df < 1) {
	return 0.0;
    }

    /*
     * Calculate the probabilty that the chi square statistic in our
     * assumed statistical model is greater than or equal to chisq.
     * If the result is small, we can reject the null hypothesis that
     * the data are consistent with the model.
     */
    assert(!isnan(chisq));
    assert(chisq >= 0.0);
    return chisq_pvalue(df, chisq);
}
