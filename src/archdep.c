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

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef HAVE_INSQUE

/*
 * _vna_insque: insert an element into a circular list
 *   @elem: element to insert
 *   @prev: element to insert after
 */
void _vna_insque(void *elem, void *prev)
{
    list_t *lp_elem = elem;
    list_t *lp_prev = prev;
    list_t *lp_next = lp_prev->l_forw;

    lp_elem->l_forw = lp_next;
    lp_elem->l_back = lp_prev;
    lp_prev->l_forw = lp_elem;
    lp_next->l_back = lp_elem;
}

#endif /* HAVE_INSQUE */
#ifndef HAVE_REMQUE

/*
 * _vna_remque: remove an element from a circular list
 *   @elem: element to remove
 */
void _vna_remque(void *elem)
{
    list_t *lp_elem = elem;
    list_t *lp_prev = lp_elem->l_back;
    list_t *lp_next = lp_elem->l_forw;

    lp_prev->l_forw = lp_next;
    lp_next->l_back = lp_prev;
    lp_elem->l_forw = NULL;
    lp_elem->l_back = NULL;
}

#endif /* HAVE_REMQUE */
#ifndef HAVE_STRCASECMP

/*
 * _vna_strcasecmp: case-insensitive string compare (4.4BSD, POSIX.1-2001)
 */
int _vna_strcasecmp(const char *s1, const char *s2)
{
    unsigned char c1, c2;

    do {
	c1 = *s1++;
	c2 = *s2++;

	if (isupper(c1)) {
	    c1 = tolower(c1);
	}
	if (isupper(c2)) {
	    c2 = tolower(c2);
	}
	if (c1 < c2) {
	    return -1;
	}
	if (c1 > c2) {
	    return  1;
	}
    } while (c1 != '\000' && c2 != '\000');
    return 0;
}
#endif /* HAVE_STRCASECMP */
#ifndef HAVE_STRDUP

/*
 * _vna_strdup: duplicate a string (XOPEN >= 500, POSIX >= 200809, BSD, SVID)
 */
char *_vna_strdup(const char *s)
{
    size_t length = strlen(s);
    char *buf;

    if ((buf = malloc(length + 1)) == NULL) {
	return NULL;
    }
    (void)strcpy(buf, s);
    return buf;
}
#endif /* HAVE_STRDUP */
#ifndef HAVE_VASPRINTF

/*
 * _vna_vasprintf: print to allocated string (GNU)
 */
int _vna_vasprintf(char **strp, const char *fmt, va_list ap)
{
    char *buf;
    int rv;

    if (strp == NULL) {
	errno = EINVAL;
	return -1;
    }
#define MAX_BUF	4095
    if ((buf = malloc(MAX_BUF + 1)) == NULL) {
	return -1;
    }
    (void)memset((void *)buf, '\000', MAX_BUF + 1);
    if ((rv = vsnprintf(buf, MAX_BUF + 1, fmt, ap)) == -1) {
	free((void *)buf);
	return -1;
    }
    if (rv < (MAX_BUF + 1) / 2) {
	char *cp;

	if ((cp = realloc(buf, rv + 1)) != NULL) {
	    buf = cp;
	}
    }
    *strp = buf;
    return rv;
}
#endif /* HAVE_VASPRINTF */
