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

#include "archdep.h"

#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnaerr_internal.h"


/*
 * vnaerr_verror: report an error
 *   @error_fn: user-supplied callback function (can be NULL)
 *   @error_arg: user-supplied data to callback function (can be NULL)
 *   @category: category of error
 *   @format: printf format string
 *   @va_list: argument list
 */
void _vnaerr_verror(vnaerr_error_fn_t *error_fn, void *error_arg,
	vnaerr_category_t category, const char *format, va_list ap)
{
    int new_errno;
    char *message = NULL;

    /*
     * Select the new errno value based on the category.  For system
     * errors, errno should be already be set on entry.
     */
    switch (category) {
    case VNAERR_SYSTEM:
	new_errno = errno;
	break;

    case VNAERR_USAGE:
	new_errno = EINVAL;
	break;

    case VNAERR_VERSION:
	new_errno = ENOPROTOOPT;
	break;

    case VNAERR_SYNTAX:
	new_errno = EBADMSG;
	break;

    case VNAERR_WARNING:
	new_errno = 0;
	break;

    case VNAERR_MATH:
	new_errno = EDOM;
	break;

    case VNAERR_INTERNAL:
    default:
	new_errno = ENOSYS;
	break;
    }

    /*
     * If the user supplied a callback function, call it.  If the
     * vasprintf fails, still try to report the original error.
     */
    if (error_fn != NULL) {
	if (vasprintf(&message, format, ap) == -1) {
	    errno = new_errno;
	    (*error_fn)(strerror(new_errno), error_arg, category);
	    goto out;
	}
	errno = new_errno;
	(*error_fn)(message, error_arg, category);
    }

out:
    free((void *)message);
    errno = new_errno;
}
