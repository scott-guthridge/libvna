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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef VNAERR_INTERNAL_H
#define VNAERR_INTERNAL_H

#include <stdarg.h>
#include "vnaerr.h"

#ifdef __cplusplus
extern "C" {
#endif

/* _vnaerr_error: report an error */
extern void _vnaerr_verror(vnaerr_error_fn_t *error_fn, void *error_arg,
	vnaerr_category_t category, const char *format, va_list ap);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* VNAERR_INTERNAL_H */
