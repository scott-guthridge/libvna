/*
 * Vector Network Analyzer Calibration Library
 * Copyright © 2020 D Scott Guthridge <scott_guthridge@rompromity.net>
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
#include "vnacal_internal.h"


/*
 * vnacal_error: report an error
 *   @vccp: configuration
 *   @format: printf format string
 */
void _vnacal_error(const vnacal_t *vcp, const char *format, ...)
{
    va_list ap;
    char *msg = NULL;
    int rv;

    if (vcp->vc_error_fn != NULL) {
	va_start(ap, format);
	rv = vasprintf(&msg, format, ap);
	va_end(ap);
	if (rv == -1) {
	    (*vcp->vc_error_fn)(strerror(errno), vcp->vc_error_arg);
	    goto out;
	}
	(*vcp->vc_error_fn)(msg, vcp->vc_error_arg);
    }

out:
    free((void *)msg);
}
