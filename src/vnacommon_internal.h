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

#ifndef VNACOMMON_INTERNAL_H
#define VNACOMMON_INTERNAL_H

#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* _vnacommon_lu: find replace A11 with its LU decomposition */
extern double complex _vnacommon_lu(complex double *a, int *row_index, int n);

/* _vnacommon_mmultiply: find C = A x B */
extern void _vnacommon_mmultiply(double complex *c, const double complex *a,
        const double complex *b, int m, int n, int o);

/* _vnacommon_mldivide: find X = A^-1 * B, X mxn, A mxm, B mxn (A destroyed) */
double complex _vnacommon_mldivide(complex double *x, complex double *a,
	const double complex *b, int m, int n);

/* _vnacommon_mrdivide: find X = B * A^-1, X mxn, B mxn, A nxn (A destroyed) */
extern double complex _vnacommon_mrdivide(double complex *x,
	const double complex *b, double complex *a, int m, int n);

/* _vnacommon_minverse: find X = A^-1, X nxn, A nxn (A destroyed) */
extern double complex _vnacommon_minverse(complex double *x, complex double *a,
	int n);

/* _vnacommon_qrd: find the QR decomposition of A, destroying A */
extern void _vnacommon_qrd(complex double *a, complex double *d,
	int rows, int columns);

/* _vnacommon_qr: find the QR decomposition of A */
extern int _vnacommon_qr(complex double *a, complex double *q,
	complex double *r, int m, int n);

/* _vnacommon_qrsolve: solve the system A X = B, destroying A and B */
extern int _vnacommon_qrsolve(complex double *x, complex double *a,
	complex double *b, int m, int n, int o);

/* _vnacommon_qrsolve2: solve the system Q R X = B */
extern void _vnacommon_qrsolve2(double complex *x, const double complex *q,
	const double complex *r, const double complex *b,
	const int m, const int n, const int o);

/* _vnacommon_spline_calc: find natural cubic spline coefficients */
extern int _vnacommon_spline_calc(int n, const double *x_vector,
	const double *y_vector, double (*c_vector)[3]);

/* _vnacommon_spline_eval: evaluate the spline at x */
extern double _vnacommon_spline_eval(int n, const double *x_vector,
	const double *y_vector, const double (*c_vector)[3], double x);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* VNACOMMON_INTERNAL_H */
