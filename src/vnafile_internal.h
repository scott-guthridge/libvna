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

#ifndef _VNAFILE_INTERNAL_H
#define _VNAFILE_INTERNAL_H

#include <stdint.h>
#include <stdio.h>
#include "vnafile.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * PI: constant
 */
#ifndef PI
#define PI	3.14159265358979323846264338327950288419716939937508
#endif

/*
 * vnafile_format_t: describes whether result is real or complex, and
 *      whether it should be printed as scalar, decibels, rectangular,
 *      polar, RC, RL, or VSWR
 */
typedef enum vnafile_format_type {
    VNAFILE_FORMAT_DB_ANGLE,	/* dB and angle */
    VNAFILE_FORMAT_MAG_ANGLE,	/* magnitude, angle */
    VNAFILE_FORMAT_REAL_IMAG,	/* real, imaginary */
    VNAFILE_FORMAT_PRC,		/* parallel R-C (VPT_ZIN only) */
    VNAFILE_FORMAT_PRL,		/* parallel R-L (VPT_ZIN only) */
    VNAFILE_FORMAT_SRC,		/* series   R-C (VPT_ZIN only) */
    VNAFILE_FORMAT_SRL,		/* series   R-L (VPT_ZIN only) */
    VNAFILE_FORMAT_IL,		/* insertion loss (VPT_S only) */
    VNAFILE_FORMAT_RL,		/* return loss    (VPT_S only) */
    VNAFILE_FORMAT_VSWR		/* voltage standing wave ratio (VPT_S only) */
} vnafile_format_type_t;

/*
 * vnafile_format_t: parsed format descriptor
 */
typedef struct vnafile_format {
    vnadata_parameter_type_t vff_parameter;
    vnafile_format_type_t    vff_format;
} vnafile_format_t;

/*
 * vnafile_t: format information for loading/saving network parameters
 */
struct vnafile {
    vnafile_error_fn_t *vf_error_fn;
    void *vf_error_arg;
    vnafile_type_t vf_type;
    vnafile_format_t *vf_format_vector;
    int vf_format_count;
    char *vf_format_string;
    int vf_fprecision;
    int vf_dprecision;
};

/* vnafile_error: report an error */
extern void _vnafile_error(const vnafile_t *vfp, const char *format, ...);

/* _vnafile_format_to_name: return the given format descriptor as a string */
extern const char *_vnafile_format_to_name(const vnafile_format_t *vffp);

/* _vnafile_update_format_string: recompute vf_format_string */
extern int _vnafile_update_format_string(vnafile_t *vfp);

/* _vnafile_find_type: try to determine the file type from the filename */
extern vnafile_type_t _vnafile_find_type(const char *filename, int *ports);

/* _vnafile_set_simple_format: set a single parameter */
extern int _vnafile_set_simple_format(vnafile_t *vfp,
	vnadata_parameter_type_t parameter, vnafile_format_type_t format);

/* _vnafile_load_native: load a native format file */
extern int _vnafile_load_native(vnafile_t *vfp, FILE *fp,
	const char *filename, vnadata_t *vdp);

/* _vnafile_load_touchstone: load a touchstone file */
extern int _vnafile_load_touchstone(vnafile_t *vfp, FILE *fp,
	const char *filename, vnadata_t *vdp);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _VNAFILE_INTERNAL_H */
