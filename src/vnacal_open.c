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

#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"


/*
 * _make_dirs: make directories as needed after &filename[offset]
 *   @filename: calibration file name
 *   @offset: offset of first component after $HOME
 */
static int _make_dirs(vnacal_t *vcp, char *filename, int offset)
{
    char *start = &filename[offset];
    char *cp;
    struct stat ss;

    while ((cp = strchr(start, '/')) != NULL) {
	*cp = '\000';
	if (stat(filename, &ss) == -1) {
	    if (errno != ENOENT) {
		_vnacal_error(vcp, "stat: %s: %s", filename, strerror(errno));
		*cp = '/';
		return -1;
	    }
	    if (mkdir(filename, 0777) == -1) {
		_vnacal_error(vcp, "mkdir: %s: %s", filename, strerror(errno));
		*cp = '/';
		return -1;
	    }
	}
	*cp = '/';
	start = cp + 1;
    }
    return 0;
}

/*
 * _vnacal_open: open a calibration file with given fopen mode
 *   @vcp: a pointer to the structure returned by vnacal_create or vnacal_load
 *   @pathname: calibration file name
 *   @dotdir: directory under $HOME or NULL
 *   @mode: mode (see fopen(3))
 */
FILE *_vnacal_open(vnacal_t *vcp, const char *pathname, const char *dotdir,
	const char *mode)
{
    size_t pathname_length = strlen(pathname);
    bool has_extension;
    char *filename;
    const char *home = NULL;
    FILE *fp;
    static const char extension[] = ".vnacal";

    /*
     * Test if the given filename contains a .vnacal extension.
     */
    has_extension = (pathname_length >= sizeof(extension)-1) &&
	    strcmp(&pathname[pathname_length - (sizeof(extension)-1)],
		    extension) == 0;

    /*
     * If opening for read, or if given an absolute path, or if an extension
     * was given, try pathname unmodified.
     */
    if (mode[0] == 'r' || pathname[0] == '/' || has_extension) {
	if ((filename = malloc(pathname_length + 1)) == NULL) {
	    _vnacal_error(vcp, "malloc: %s", strerror(errno));
	    return NULL;
	}
	(void)memcpy((void *)filename, (void *)pathname, pathname_length + 1);
	if ((fp = fopen(filename, mode)) != NULL) {
	    free((void *)vcp->vc_filename);
	    vcp->vc_filename = filename;
	    return fp;
	}
	if (mode[0] != 'r' || errno != ENOENT || has_extension) {
	    _vnacal_error(vcp, "%s: %s", filename, strerror(errno));
	    goto error;
	}
	free((void *)filename);
	filename = NULL;
    }
    assert(!has_extension);

    /*
     * If read, try adding the extension.
     */
    if (mode[0] == 'r') {
	char *cp;

	if ((filename = malloc(pathname_length + sizeof(extension))) == NULL) {
	    _vnacal_error(vcp, "malloc: %s", strerror(errno));
	    return NULL;
	}
	cp = filename;
	(void)memcpy((void *)cp, (void *)pathname, pathname_length);
	cp += pathname_length;
	(void)memcpy((void *)cp, (void *)extension, sizeof(extension));
	if ((fp = fopen(filename, mode)) != NULL) {
	    free((void *)vcp->vc_filename);
	    vcp->vc_filename = filename;
	    return fp;
	}
	if (errno != ENOENT) {
	    _vnacal_error(vcp, "%s: %s", filename, strerror(errno));
	    goto error;
	}
	free((void *)filename);
	filename = NULL;
    }

    /*
     * If pathname is relative, dotdir is given and $HOME is set,
     * try $HOME/{dotdir}/{pathname}.vnacal
     */
    if (pathname[0] != '/' && dotdir != NULL &&
	    (home = getenv("HOME")) != NULL) {
	size_t home_length = strlen(home);
	size_t dotdir_length = strlen(dotdir);
	char *cp;

	if ((filename = malloc(home_length + 1 + dotdir_length + 1 +
			pathname_length + sizeof(extension))) == NULL) {
	    _vnacal_error(vcp, "malloc: %s", strerror(errno));
	    return NULL;
	}
	cp = filename;
	(void)memcpy((void *)cp, (void *)home, home_length);
	cp += home_length;
	*cp++ = '/';
	(void)memcpy((void *)cp, (void *)dotdir, dotdir_length);
	cp += dotdir_length;
	*cp++ = '/';
	(void)memcpy((void *)cp, (void *)pathname, pathname_length + 1);
	cp += pathname_length;
	(void)memcpy((void *)cp, (void *)extension, sizeof(extension));
	if (mode[0] == 'w' || mode[0] == 'a') {
	    if (_make_dirs(vcp, filename, home_length + 1) == -1) {
		free((void *)filename);
		filename = NULL;
		return NULL;
	    }
	}
	if ((fp = fopen(filename, mode)) == NULL) {
	    _vnacal_error(vcp, "%s: %s", filename, strerror(errno));
	    goto error;
	}
	free((void *)vcp->vc_filename);
	vcp->vc_filename = filename;
	return fp;
    }
    _vnacal_error(vcp, "%s: %s", pathname, strerror(errno));
    /*FALLTHROUGH*/

error:
    free((void *)filename);
    filename = NULL;
    return NULL;
}
