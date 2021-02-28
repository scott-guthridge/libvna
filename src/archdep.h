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
#include "config.h"

#define _GNU_SOURCE	/* for declaration of vasprintf on Linux */

#ifndef MAX
#define MAX(a, b)       ((a) >= (b) ? (a) : (b))
#endif /* MAX */

#ifndef MIN
#define MIN(a, b)       ((a) <= (b) ? (a) : (b))
#endif /* MIN */

#ifndef HAVE_ISASCII	/* C89, but not C99 */
#define isascii(c)	((c) >= 0 && (c) <= 0x7f)
#endif /* HAVE_ISASCII */

#ifndef HAVE_RANDOM	/* POSIX.1-2001, POSIX.1-2008, 4.3BSD */
#define random rand
#define srandom srand
#endif

#ifndef HAVE_STRCASECMP /* 4.4BSD, POSIX.1-2001, POSIX.1-2008 */
extern int _vna_strcasecmp(const char *s1, const char *s2);
#define strcasecmp _vna_strcasecmp
#endif /* HAVE_STRCASECMP */

#ifndef HAVE_STRDUP	/* SVr4, 4.3BSD, POSIX.1-2001 */
extern char *_vna_strdup(const char *s);
#define strdup _vna_strdup
#endif /* HAVE_STRDUP */

#ifndef HAVE_VASPRINTF	/* GNU */
#include <stdarg.h>
extern int _vna_vasprintf(char **strp, const char *fmt, va_list ap);
#define vasprintf _vna_vasprintf
#endif /* HAVE_VASPRINTF */
