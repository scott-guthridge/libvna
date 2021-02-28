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
extern void vnaconv_a2b(const vnaconv_array2_t *a, double complex (*b)[2]);
extern void vnaconv_a2g(const vnaconv_array2_t *a, double complex (*g)[2]);
extern void vnaconv_a2h(const vnaconv_array2_t *a, double complex (*h)[2]);
extern void vnaconv_a2s(const vnaconv_array2_t *a, double complex (*s)[2],
	const double complex *z0);
extern void vnaconv_a2t(const vnaconv_array2_t *a, double complex (*t)[2],
	const double complex *z0);
extern void vnaconv_a2y(const vnaconv_array2_t *a, double complex (*y)[2]);
extern void vnaconv_a2z(const vnaconv_array2_t *a, double complex (*z)[2]);
extern void vnaconv_a2zi(const vnaconv_array2_t *a, double complex *zi,
	const double complex *z0);
extern void vnaconv_b2a(const vnaconv_array2_t *b, double complex (*a)[2]);
extern void vnaconv_b2g(const vnaconv_array2_t *b, double complex (*g)[2]);
extern void vnaconv_b2h(const vnaconv_array2_t *b, double complex (*h)[2]);
extern void vnaconv_b2s(const vnaconv_array2_t *b, double complex (*s)[2],
	const double complex *z0);
extern void vnaconv_b2t(const vnaconv_array2_t *b, double complex (*t)[2], const double complex *z0);
extern void vnaconv_b2y(const vnaconv_array2_t *b, double complex (*y)[2]);
extern void vnaconv_b2z(const vnaconv_array2_t *b, double complex (*z)[2]);
extern void vnaconv_b2zi(const vnaconv_array2_t *b, double complex *zi,
	const double complex *z0);
extern void vnaconv_g2a(const vnaconv_array2_t *g, double complex (*a)[2]);
extern void vnaconv_g2b(const vnaconv_array2_t *g, double complex (*b)[2]);
extern void vnaconv_g2h(const vnaconv_array2_t *g, double complex (*h)[2]);
extern void vnaconv_g2s(const vnaconv_array2_t *g, double complex (*s)[2],
	const double complex *z0);
extern void vnaconv_g2t(const vnaconv_array2_t *g, double complex (*t)[2],
	const double complex *z0);
extern void vnaconv_g2y(const vnaconv_array2_t *g, double complex (*y)[2]);
extern void vnaconv_g2z(const vnaconv_array2_t *g, double complex (*z)[2]);
extern void vnaconv_g2zi(const vnaconv_array2_t *g, double complex *zi,
	const double complex *z0);
extern void vnaconv_h2a(const vnaconv_array2_t *h, double complex (*a)[2]);
extern void vnaconv_h2b(const vnaconv_array2_t *h, double complex (*b)[2]);
extern void vnaconv_h2g(const vnaconv_array2_t *h, double complex (*g)[2]);
extern void vnaconv_h2s(const vnaconv_array2_t *h, double complex (*s)[2],
	const double complex *z0);
extern void vnaconv_h2t(const vnaconv_array2_t *h, double complex (*t)[2],
	const double complex *z0);
extern void vnaconv_h2y(const vnaconv_array2_t *h, double complex (*y)[2]);
extern void vnaconv_h2z(const vnaconv_array2_t *h, double complex (*z)[2]);
extern void vnaconv_h2zi(const vnaconv_array2_t *h, double complex *zi,
	const double complex *z0);
extern void vnaconv_s2a(const vnaconv_array2_t *s, double complex (*a)[2],
	const double complex *z0);
extern void vnaconv_s2b(const vnaconv_array2_t *s, double complex (*b)[2],
	const double complex *z0);
extern void vnaconv_s2g(const vnaconv_array2_t *s, double complex (*g)[2],
	const double complex *z0);
extern void vnaconv_s2h(const vnaconv_array2_t *s, double complex (*h)[2],
	const double complex *z0);
extern void vnaconv_s2t(const vnaconv_array2_t *s, double complex (*t)[2]);
extern void vnaconv_s2y(const vnaconv_array2_t *s, double complex (*y)[2],
	const double complex *z0);
extern void vnaconv_s2z(const vnaconv_array2_t *s, double complex (*z)[2],
	const double complex *z0);
extern void vnaconv_s2zi(const vnaconv_array2_t *s, double complex *zi,
	const double complex *z0);
extern void vnaconv_t2a(const vnaconv_array2_t *t, double complex (*a)[2],
	const double complex *z0);
extern void vnaconv_t2b(const vnaconv_array2_t *t, double complex (*b)[2],
	const double complex *z0);
extern void vnaconv_t2g(const vnaconv_array2_t *t, double complex (*g)[2],
	const double complex *z0);
extern void vnaconv_t2h(const vnaconv_array2_t *t, double complex (*h)[2],
	const double complex *z0);
extern void vnaconv_t2s(const vnaconv_array2_t *t, double complex (*s)[2]);
extern void vnaconv_t2y(const vnaconv_array2_t *t, double complex (*y)[2],
	const double complex *z0);
extern void vnaconv_t2z(const vnaconv_array2_t *t, double complex (*z)[2],
	const double complex *z0);
extern void vnaconv_t2zi(const vnaconv_array2_t *t, double complex *zi,
	const double complex *z0);
extern void vnaconv_y2a(const vnaconv_array2_t *y, double complex (*a)[2]);
extern void vnaconv_y2b(const vnaconv_array2_t *y, double complex (*b)[2]);
extern void vnaconv_y2g(const vnaconv_array2_t *y, double complex (*g)[2]);
extern void vnaconv_y2h(const vnaconv_array2_t *y, double complex (*h)[2]);
extern void vnaconv_y2s(const vnaconv_array2_t *y, double complex (*s)[2],
	const double complex *z0);
extern void vnaconv_y2t(const vnaconv_array2_t *y, double complex (*t)[2],
	const double complex *z0);
extern void vnaconv_y2z(const vnaconv_array2_t *y, double complex (*z)[2]);
extern void vnaconv_y2zi(const vnaconv_array2_t *y, double complex *zi,
	const double complex *z0);
extern void vnaconv_z2a(const vnaconv_array2_t *z, double complex (*a)[2]);
extern void vnaconv_z2b(const vnaconv_array2_t *z, double complex (*b)[2]);
extern void vnaconv_z2g(const vnaconv_array2_t *z, double complex (*g)[2]);
extern void vnaconv_z2h(const vnaconv_array2_t *z, double complex (*h)[2]);
extern void vnaconv_z2s(const vnaconv_array2_t *z, double complex (*s)[2],
	const double complex *z0);
extern void vnaconv_z2t(const vnaconv_array2_t *z, double complex (*t)[2],
	const double complex *z0);
extern void vnaconv_z2y(const vnaconv_array2_t *z, double complex (*y)[2]);
extern void vnaconv_z2zi(const vnaconv_array2_t *z, double complex *zi,
	const double complex *z0);

/*
 * NxN conversions
 */
extern void vnaconv_s2yn(const double complex *s, double complex *y,
	const double complex *z0, int n);
extern void vnaconv_s2zin(const double complex *s, double complex *zi,
	const double complex *z0, int n);
extern void vnaconv_s2zn(const double complex *s, double complex *z,
	const double complex *z0, int n);
extern void vnaconv_y2sn(const double complex *y, double complex *s,
	const double complex *z0, int n);
extern void vnaconv_y2zin(const double complex *y, double complex *zi,
	const double complex *z0, int n);
extern void vnaconv_y2zn(const double complex *y, double complex *s, int n);
extern void vnaconv_z2sn(const double complex *z, double complex *s,
	const double complex *z0, int n);
extern void vnaconv_z2yn(const double complex *z, double complex *y, int n);
extern void vnaconv_z2zin(const double complex *z, double complex *zi,
	const double complex *z0, int n);

/*
 * MxN conversions
 */
extern void vnaconv_s2zimn(const double complex *s, double complex *zi,
	const double complex *z0, int rows, int columns);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* VNACONV_H */
