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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _LIBT_VNACAL_H
#define _LIBT_VNACAL_H

#include <stdbool.h>
#include "vnacal_internal.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum libt_vnadata_z0_type {
    Z0_SINGLE,
    Z0_REAL_VECTOR,
    Z0_COMPLEX_VECTOR,
    Z0_PER_F,
    Z0_NTYPES
} libt_vnadata_z0_type_t;

typedef enum libt_vnadata_fill_method {
    FM_CELL,
    FM_MATRIX,
    FM_VECTOR,
    FM_NMETHODS
} libt_vnadata_fill_method_t;

/*
 * test_data
 */
typedef struct libt_vnadata {
    vnadata_parameter_type_t td_type;
    int td_rows;
    int td_columns;
    int td_frequencies;
    double complex **td_vector;
    double *td_frequency_vector;
    libt_vnadata_z0_type_t td_z0_type;
    union {
	double complex *td_z0_vector;
	double complex **td_fz0_vector;
    } u;
} libt_vnadata_t;
#define td_z0_vector	u.td_z0_vector
#define td_fz0_vector	u.td_fz0_vector

extern const char *libt_vnadata_z0_names[];
extern const char *libt_vnadata_fill_names[];
extern libt_vnadata_t *libt_vnadata_create(vnadata_parameter_type_t type,
	int rows, int columns, int frequencies,
	libt_vnadata_z0_type_t z0_type);
extern void libt_vnadata_free(libt_vnadata_t *tdp);
extern libt_result_t libt_vnadata_validate(const libt_vnadata_t *tdp,
	const vnadata_t *vdp);
extern libt_result_t libt_vnadata_fill(const libt_vnadata_t *tdp,
	vnadata_t *vdp, libt_vnadata_fill_method_t fill_method);
extern void libt_vnadata_convert(const double complex *in, double complex *out,
	const double complex *z0, int rows, int columns,
	vnadata_parameter_type_t old_type, vnadata_parameter_type_t new_type);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _LIBT_VNACAL_H */
