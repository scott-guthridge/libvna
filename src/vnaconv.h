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

#ifndef VNACONV_H
#define VNACONV_H

#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * vnaconv_array2_t: array of two double complex
 *   This is a work around a syntactic limitation of const
 *   in declarators.
 */
typedef double complex vnaconv_array2_t[2];

/*
 * 2x2 Conversions
 */
extern void vnaconv_atob(const vnaconv_array2_t *a, double complex (*b)[2]);
extern void vnaconv_atog(const vnaconv_array2_t *a, double complex (*g)[2]);
extern void vnaconv_atoh(const vnaconv_array2_t *a, double complex (*h)[2]);
extern void vnaconv_atos(const vnaconv_array2_t *a, double complex (*s)[2],
	const double complex *z0);
extern void vnaconv_atot(const vnaconv_array2_t *a, double complex (*t)[2],
	const double complex *z0);
extern void vnaconv_atoy(const vnaconv_array2_t *a, double complex (*y)[2]);
extern void vnaconv_atoz(const vnaconv_array2_t *a, double complex (*z)[2]);
extern void vnaconv_atozi(const vnaconv_array2_t *a, double complex *zi,
	const double complex *z0);
extern void vnaconv_btoa(const vnaconv_array2_t *b, double complex (*a)[2]);
extern void vnaconv_btog(const vnaconv_array2_t *b, double complex (*g)[2]);
extern void vnaconv_btoh(const vnaconv_array2_t *b, double complex (*h)[2]);
extern void vnaconv_btos(const vnaconv_array2_t *b, double complex (*s)[2],
	const double complex *z0);
extern void vnaconv_btot(const vnaconv_array2_t *b, double complex (*t)[2],
	const double complex *z0);
extern void vnaconv_btoy(const vnaconv_array2_t *b, double complex (*y)[2]);
extern void vnaconv_btoz(const vnaconv_array2_t *b, double complex (*z)[2]);
extern void vnaconv_btozi(const vnaconv_array2_t *b, double complex *zi,
	const double complex *z0);
extern void vnaconv_gtoa(const vnaconv_array2_t *g, double complex (*a)[2]);
extern void vnaconv_gtob(const vnaconv_array2_t *g, double complex (*b)[2]);
extern void vnaconv_gtoh(const vnaconv_array2_t *g, double complex (*h)[2]);
extern void vnaconv_gtos(const vnaconv_array2_t *g, double complex (*s)[2],
	const double complex *z0);
extern void vnaconv_gtot(const vnaconv_array2_t *g, double complex (*t)[2],
	const double complex *z0);
extern void vnaconv_gtoy(const vnaconv_array2_t *g, double complex (*y)[2]);
extern void vnaconv_gtoz(const vnaconv_array2_t *g, double complex (*z)[2]);
extern void vnaconv_gtozi(const vnaconv_array2_t *g, double complex *zi,
	const double complex *z0);
extern void vnaconv_htoa(const vnaconv_array2_t *h, double complex (*a)[2]);
extern void vnaconv_htob(const vnaconv_array2_t *h, double complex (*b)[2]);
extern void vnaconv_htog(const vnaconv_array2_t *h, double complex (*g)[2]);
extern void vnaconv_htos(const vnaconv_array2_t *h, double complex (*s)[2],
	const double complex *z0);
extern void vnaconv_htot(const vnaconv_array2_t *h, double complex (*t)[2],
	const double complex *z0);
extern void vnaconv_htoy(const vnaconv_array2_t *h, double complex (*y)[2]);
extern void vnaconv_htoz(const vnaconv_array2_t *h, double complex (*z)[2]);
extern void vnaconv_htozi(const vnaconv_array2_t *h, double complex *zi,
	const double complex *z0);
extern void vnaconv_stoa(const vnaconv_array2_t *s, double complex (*a)[2],
	const double complex *z0);
extern void vnaconv_stob(const vnaconv_array2_t *s, double complex (*b)[2],
	const double complex *z0);
extern void vnaconv_stog(const vnaconv_array2_t *s, double complex (*g)[2],
	const double complex *z0);
extern void vnaconv_stoh(const vnaconv_array2_t *s, double complex (*h)[2],
	const double complex *z0);
extern void vnaconv_stot(const vnaconv_array2_t *s, double complex (*t)[2]);
extern void vnaconv_stoy(const vnaconv_array2_t *s, double complex (*y)[2],
	const double complex *z0);
extern void vnaconv_stoz(const vnaconv_array2_t *s, double complex (*z)[2],
	const double complex *z0);
extern void vnaconv_stozi(const vnaconv_array2_t *s, double complex *zi,
	const double complex *z0);
extern void vnaconv_ttoa(const vnaconv_array2_t *t, double complex (*a)[2],
	const double complex *z0);
extern void vnaconv_ttob(const vnaconv_array2_t *t, double complex (*b)[2],
	const double complex *z0);
extern void vnaconv_ttog(const vnaconv_array2_t *t, double complex (*g)[2],
	const double complex *z0);
extern void vnaconv_ttoh(const vnaconv_array2_t *t, double complex (*h)[2],
	const double complex *z0);
extern void vnaconv_ttos(const vnaconv_array2_t *t, double complex (*s)[2]);
extern void vnaconv_ttoy(const vnaconv_array2_t *t, double complex (*y)[2],
	const double complex *z0);
extern void vnaconv_ttoz(const vnaconv_array2_t *t, double complex (*z)[2],
	const double complex *z0);
extern void vnaconv_ttozi(const vnaconv_array2_t *t, double complex *zi,
	const double complex *z0);
extern void vnaconv_ytoa(const vnaconv_array2_t *y, double complex (*a)[2]);
extern void vnaconv_ytob(const vnaconv_array2_t *y, double complex (*b)[2]);
extern void vnaconv_ytog(const vnaconv_array2_t *y, double complex (*g)[2]);
extern void vnaconv_ytoh(const vnaconv_array2_t *y, double complex (*h)[2]);
extern void vnaconv_ytos(const vnaconv_array2_t *y, double complex (*s)[2],
	const double complex *z0);
extern void vnaconv_ytot(const vnaconv_array2_t *y, double complex (*t)[2],
	const double complex *z0);
extern void vnaconv_ytoz(const vnaconv_array2_t *y, double complex (*z)[2]);
extern void vnaconv_ytozi(const vnaconv_array2_t *y, double complex *zi,
	const double complex *z0);
extern void vnaconv_ztoa(const vnaconv_array2_t *z, double complex (*a)[2]);
extern void vnaconv_ztob(const vnaconv_array2_t *z, double complex (*b)[2]);
extern void vnaconv_ztog(const vnaconv_array2_t *z, double complex (*g)[2]);
extern void vnaconv_ztoh(const vnaconv_array2_t *z, double complex (*h)[2]);
extern void vnaconv_ztos(const vnaconv_array2_t *z, double complex (*s)[2],
	const double complex *z0);
extern void vnaconv_ztot(const vnaconv_array2_t *z, double complex (*t)[2],
	const double complex *z0);
extern void vnaconv_ztoy(const vnaconv_array2_t *z, double complex (*y)[2]);
extern void vnaconv_ztozi(const vnaconv_array2_t *z, double complex *zi,
	const double complex *z0);

/*
 * NxN conversions
 */
extern void vnaconv_stoyn(const double complex *s, double complex *y,
	const double complex *z0, int n);
extern void vnaconv_stozin(const double complex *s, double complex *zi,
	const double complex *z0, int n);
extern void vnaconv_stozn(const double complex *s, double complex *z,
	const double complex *z0, int n);
extern void vnaconv_ytosn(const double complex *y, double complex *s,
	const double complex *z0, int n);
extern void vnaconv_ytozin(const double complex *y, double complex *zi,
	const double complex *z0, int n);
extern void vnaconv_ytozn(const double complex *y, double complex *s, int n);
extern void vnaconv_ztosn(const double complex *z, double complex *s,
	const double complex *z0, int n);
extern void vnaconv_ztoyn(const double complex *z, double complex *y, int n);
extern void vnaconv_ztozin(const double complex *z, double complex *zi,
	const double complex *z0, int n);

/*
 * MxN conversions
 */
extern void vnaconv_stozimn(const double complex *s, double complex *zi,
	const double complex *z0, int rows, int columns);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* VNACONV_H */
