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

#ifndef _VNAERR_H
#define _VNAERR_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * vnaerr_category_t: category of error
 */
typedef enum vnaerr_category {
    VNAERR_SYSTEM,	/* system error, e.g. ENOMEM, EPERM */
    VNAERR_USAGE,	/* invalid argument to library function */
    VNAERR_VERSION,	/* unsupported file version */
    VNAERR_SYNTAX,	/* syntax error in file */
    VNAERR_WARNING,	/* non-fatal warning */
    VNAERR_MATH,	/* singular matrix or convergence failure */
    VNAERR_INTERNAL	/* unexpected internal inconsistency */
} vnaerr_category_t;

/*
 * vnaerr_error_fn_t: error reporting function type
 *   @message: error message without newline
 *   @error_arg: user-supplied argument passed through to error function
 *   @category: category of error
 */
typedef void vnaerr_error_fn_t(const char *message, void *error_arg,
	vnaerr_category_t category);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _VNAERR_H */
