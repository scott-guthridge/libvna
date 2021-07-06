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
 * 2x2 Conversions
 */
extern void vnaconv_atob(const double complex (*a)[2], double complex (*b)[2]);
extern void vnaconv_atog(const double complex (*a)[2], double complex (*g)[2]);
extern void vnaconv_atoh(const double complex (*a)[2], double complex (*h)[2]);
extern void vnaconv_atos(const double complex (*a)[2], double complex (*s)[2],
	const double complex *z0);
extern void vnaconv_atot(const double complex (*a)[2], double complex (*t)[2],
	const double complex *z0);
extern void vnaconv_atoy(const double complex (*a)[2], double complex (*y)[2]);
extern void vnaconv_atoz(const double complex (*a)[2], double complex (*z)[2]);
extern void vnaconv_atozi(const double complex (*a)[2], double complex *zi,
	const double complex *z0);
extern void vnaconv_btoa(const double complex (*b)[2], double complex (*a)[2]);
extern void vnaconv_btog(const double complex (*b)[2], double complex (*g)[2]);
extern void vnaconv_btoh(const double complex (*b)[2], double complex (*h)[2]);
extern void vnaconv_btos(const double complex (*b)[2], double complex (*s)[2],
	const double complex *z0);
extern void vnaconv_btot(const double complex (*b)[2], double complex (*t)[2],
	const double complex *z0);
extern void vnaconv_btoy(const double complex (*b)[2], double complex (*y)[2]);
extern void vnaconv_btoz(const double complex (*b)[2], double complex (*z)[2]);
extern void vnaconv_btozi(const double complex (*b)[2], double complex *zi,
	const double complex *z0);
extern void vnaconv_gtoa(const double complex (*g)[2], double complex (*a)[2]);
extern void vnaconv_gtob(const double complex (*g)[2], double complex (*b)[2]);
extern void vnaconv_gtoh(const double complex (*g)[2], double complex (*h)[2]);
extern void vnaconv_gtos(const double complex (*g)[2], double complex (*s)[2],
	const double complex *z0);
extern void vnaconv_gtot(const double complex (*g)[2], double complex (*t)[2],
	const double complex *z0);
extern void vnaconv_gtoy(const double complex (*g)[2], double complex (*y)[2]);
extern void vnaconv_gtoz(const double complex (*g)[2], double complex (*z)[2]);
extern void vnaconv_gtozi(const double complex (*g)[2], double complex *zi,
	const double complex *z0);
extern void vnaconv_htoa(const double complex (*h)[2], double complex (*a)[2]);
extern void vnaconv_htob(const double complex (*h)[2], double complex (*b)[2]);
extern void vnaconv_htog(const double complex (*h)[2], double complex (*g)[2]);
extern void vnaconv_htos(const double complex (*h)[2], double complex (*s)[2],
	const double complex *z0);
extern void vnaconv_htot(const double complex (*h)[2], double complex (*t)[2],
	const double complex *z0);
extern void vnaconv_htoy(const double complex (*h)[2], double complex (*y)[2]);
extern void vnaconv_htoz(const double complex (*h)[2], double complex (*z)[2]);
extern void vnaconv_htozi(const double complex (*h)[2], double complex *zi,
	const double complex *z0);
extern void vnaconv_stoa(const double complex (*s)[2], double complex (*a)[2],
	const double complex *z0);
extern void vnaconv_stob(const double complex (*s)[2], double complex (*b)[2],
	const double complex *z0);
extern void vnaconv_stog(const double complex (*s)[2], double complex (*g)[2],
	const double complex *z0);
extern void vnaconv_stoh(const double complex (*s)[2], double complex (*h)[2],
	const double complex *z0);
extern void vnaconv_stot(const double complex (*s)[2], double complex (*t)[2]);
extern void vnaconv_stoy(const double complex (*s)[2], double complex (*y)[2],
	const double complex *z0);
extern void vnaconv_stoz(const double complex (*s)[2], double complex (*z)[2],
	const double complex *z0);
extern void vnaconv_stozi(const double complex (*s)[2], double complex *zi,
	const double complex *z0);
extern void vnaconv_ttoa(const double complex (*t)[2], double complex (*a)[2],
	const double complex *z0);
extern void vnaconv_ttob(const double complex (*t)[2], double complex (*b)[2],
	const double complex *z0);
extern void vnaconv_ttog(const double complex (*t)[2], double complex (*g)[2],
	const double complex *z0);
extern void vnaconv_ttoh(const double complex (*t)[2], double complex (*h)[2],
	const double complex *z0);
extern void vnaconv_ttos(const double complex (*t)[2], double complex (*s)[2]);
extern void vnaconv_ttoy(const double complex (*t)[2], double complex (*y)[2],
	const double complex *z0);
extern void vnaconv_ttoz(const double complex (*t)[2], double complex (*z)[2],
	const double complex *z0);
extern void vnaconv_ttozi(const double complex (*t)[2], double complex *zi,
	const double complex *z0);
extern void vnaconv_ytoa(const double complex (*y)[2], double complex (*a)[2]);
extern void vnaconv_ytob(const double complex (*y)[2], double complex (*b)[2]);
extern void vnaconv_ytog(const double complex (*y)[2], double complex (*g)[2]);
extern void vnaconv_ytoh(const double complex (*y)[2], double complex (*h)[2]);
extern void vnaconv_ytos(const double complex (*y)[2], double complex (*s)[2],
	const double complex *z0);
extern void vnaconv_ytot(const double complex (*y)[2], double complex (*t)[2],
	const double complex *z0);
extern void vnaconv_ytoz(const double complex (*y)[2], double complex (*z)[2]);
extern void vnaconv_ytozi(const double complex (*y)[2], double complex *zi,
	const double complex *z0);
extern void vnaconv_ztoa(const double complex (*z)[2], double complex (*a)[2]);
extern void vnaconv_ztob(const double complex (*z)[2], double complex (*b)[2]);
extern void vnaconv_ztog(const double complex (*z)[2], double complex (*g)[2]);
extern void vnaconv_ztoh(const double complex (*z)[2], double complex (*h)[2]);
extern void vnaconv_ztos(const double complex (*z)[2], double complex (*s)[2],
	const double complex *z0);
extern void vnaconv_ztot(const double complex (*z)[2], double complex (*t)[2],
	const double complex *z0);
extern void vnaconv_ztoy(const double complex (*z)[2], double complex (*y)[2]);
extern void vnaconv_ztozi(const double complex (*z)[2], double complex *zi,
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
