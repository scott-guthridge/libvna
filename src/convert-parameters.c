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

/*
 * Network parameter converter: converts between network parameter
 * types and between Touchstone 1, Touchstone 2 and NPD file format.
 * The file type is based on filename extension using ".s1p", ".s2p",
 * ".s3p", etc.  for Touchstone 1, ".ts" for Touchstone 2, and ".npd"
 * or other for NPD format.
 *
 * Example:
 *     Convert 4x4 network data from a Touchstone 1 file to Z parameters
 *     in magnitude/angle format, saving as Touchstone 2.
 *
 *     npd-convert -f zma data.s4p data.ts
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <vnadata.h>


static char *progname;

/*
 * usage: usage text
 */
static const char usage[] =
    "%s [-f format] input-file output-file\n"
    "where format is a comma-separated list of:\n"
    "  s[ri|ma|dB]  scattering parameters\n"
    "  t[ri|ma|dB]  scattering-transfer parameters\n"
    "  z[ri|ma]     impedance parameters\n"
    "  y[ri|ma]     admittance parameters\n"
    "  h[ri|ma]     hybrid parameters\n"
    "  g[ri|ma]     inverse-hybrid parameters\n"
    "  a[ri|ma]     ABCD parameters\n"
    "  b[ri|ma]     inverse ABCD parameters\n"
    "  Zin[ri|ma]   input impedances\n"
    "  PRC          Zin as parallel RC\n"
    "  PRL          Zin as parallel RL\n"
    "  SRC          Zin as series RC\n"
    "  SRL          Zin as series RL\n"
    "  IL           insertion loss\n"
    "  RL           return loss\n"
    "  VSWR         voltage standing wave ratio\n"
    "\n"
    "Coordinates\n"
    "  ri  real, imaginary\n"
    "  ma  magnitude, angle\n"
    "  dB  decibels, angle\n"
    "\n"
    "Specifiers are case-insensitive.\n";

/*
 * error_fn: error printing function for the library
 *   @message: single line error message without a newline
 *   @error_arg: passed through to the error function (unused here)
 *   @category: category of error (ignored here)
 */
static void error_fn(const char *message, void *error_arg,
	vnaerr_category_t category)
{
    (void)fprintf(stderr, "%s: %s\n", progname, message);
}

/*
 * main
 */
int main(int argc, char **argv)
{
    vnadata_t *vdp;
    const char *f_opt = NULL;

    if ((char *)NULL == (progname = strrchr(argv[0], '/'))) {
	progname = argv[0];
    } else {
	++progname;
    }
    for (;;) {
	switch (getopt(argc, argv, "f:")) {
	case -1:
	    break;

	case 'f':
	    f_opt = optarg;
	    continue;

	default:
	    (void)fprintf(stderr, usage, progname);
	    exit(2);
	}
	break;
    }
    argc -= optind;
    argv += optind;
    if (argc != 2) {
	(void)fprintf(stderr, usage, progname);
	exit(2);
    }
    if ((vdp = vnadata_alloc(error_fn, NULL)) == NULL) {
	exit(3);
    }
    if (vnadata_load(vdp, argv[0]) == -1) {
	exit(4);
    }
    /*
     * Set the filetype back to auto so that saving to a .ts file forces
     * Touchstone 2 format.
     */
    if (vnadata_set_filetype(vdp, VNADATA_FILETYPE_AUTO) == -1) {
	exit(5);
    }
    if (f_opt != NULL) {
	if (vnadata_set_format(vdp, f_opt) == -1) {
	    exit(6);
	}
    }
    if (vnadata_save(vdp, argv[1]) == -1) {
	exit(7);
    }
    vnadata_free(vdp);
    exit(0);
}
