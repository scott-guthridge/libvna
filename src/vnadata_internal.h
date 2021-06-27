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

#ifndef VNADATA_INTERNAL_H
#define VNADATA_INTERNAL_H

#include <stdbool.h>
#include <stddef.h>
#include "vnaerr_internal.h"
#include "vnadata.h"

/*
 * VDI_MAGIC: magic number
 */
#define VDI_MAGIC	0x56444930	/* VDI0 */

/*
 * vdi_flags values: optional feature flags (see vnadata_init)
 */
#define VF_PER_F_Z0	0x0001

/*
 * _VNADATA_IS_POWER: return true if parameter represents power or log-power
 */
#define _VNADATA_IS_POWER(parameter) \
	((parameter) == VPT_S || (parameter) == VPT_T)

/*
 * _VNADATA_IS_MATRIX: return true if parameter represents a convertible matrix
 */
#define _VNADATA_IS_MATRIX(parameter) \
	((parameter) != VPT_UNDEF && (parameter) != VPT_ZIN)

/*
 * vnadata_format_descriptor_t: describes whether a load/save field is real or complex,
 *	and whether it should be printed as scalar, decibels, rectangular,
 *      polar, RC, RL, or VSWR
 */
typedef enum vnadata_format {
    VNADATA_FORMAT_DB_ANGLE,	/* dB and angle */
    VNADATA_FORMAT_MAG_ANGLE,	/* magnitude, angle */
    VNADATA_FORMAT_REAL_IMAG,	/* real, imaginary */
    VNADATA_FORMAT_PRC,		/* parallel R-C (VPT_ZIN only) */
    VNADATA_FORMAT_PRL,		/* parallel R-L (VPT_ZIN only) */
    VNADATA_FORMAT_SRC,		/* series   R-C (VPT_ZIN only) */
    VNADATA_FORMAT_SRL,		/* series   R-L (VPT_ZIN only) */
    VNADATA_FORMAT_IL,		/* insertion loss (VPT_S only) */
    VNADATA_FORMAT_RL,		/* return loss    (VPT_S only) */
    VNADATA_FORMAT_VSWR		/* voltage standing wave ratio (VPT_S only) */
} vnadata_format_t;

/*
 * vnadata_format_descriptor_t: parsed load/save format descriptor
 */
typedef struct vnadata_format_descriptor {
    vnadata_parameter_type_t	vfd_parameter;
    vnadata_format_t		vfd_format;
} vnadata_format_descriptor_t;

/*
 * vnadata_internal_t: internal version of vnadata_t
 */
typedef struct vnadata_internal {
    /* magic number used to validate struct */
    uint32_t vdi_magic;

    /* flags: bitwise-OR of vnadata_flags */
    uint32_t vdi_flags;

    /* user-visible portion of this structure */
    vnadata_t vdi_vd;

    /* user-supplied error callback or NULL */
    vnaerr_error_fn_t *vdi_error_fn;

    /* user-supplied error callback argument or NULL */
    void *vdi_error_arg;

    /* allocation of vdi_z0_vector_vector */
    int vdi_p_allocation;

    /* allocation of vd_frequency_vector */
    int vdi_f_allocation;

    /* allocation of each vd_data[findex] matrix */
    int vdi_m_allocation;

    /* per-port and optionally also per-frequency system impedances */
    union {
	double complex *vdi_z0_vector;		/* frequency-independent */
	double complex **vdi_z0_vector_vector;	/* frequency-dependent */
    } vdi_z0;

    /* file format for vnadata_load */
    vnadata_filetype_t vdi_filetype;

    /* vector of field formats for load/save */
    vnadata_format_descriptor_t *vdi_format_vector;

    /* length of vdi_format_vector */
    int vdi_format_count;

    /* string version of vdi_format_vector */
    char *vdi_format_string;

    /* numeric precision for frequency values */
    int vdi_fprecision;

    /* numeric precision for data values */
    int vdi_dprecision;

} vnadata_internal_t;

/*
 * Aliases to hide union
 */
#define vdi_z0_vector		vdi_z0.vdi_z0_vector
#define vdi_z0_vector_vector	vdi_z0.vdi_z0_vector_vector

/*
 * VDP_TO_VDIP: convert from pointer to vnadata_t to vnadata_internal_t
 */
#define VDP_TO_VDIP(vdp) \
    (vnadata_internal_t *)((char *)(vdp) - offsetof(vnadata_internal_t, vdi_vd))

/* vnadata_error: report an error */
extern void _vnadata_error(const vnadata_internal_t *vdip,
	vnaerr_category_t category, const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 3, 4)));
#else
    ;
#endif

/* _vnadata_extend_p: extend the port allocation for Z0 */
extern int _vnadata_extend_p(vnadata_internal_t *vdip, int new_p_allocation);

/* _vnadata_extend_m: extend the matrix allocation */
extern int _vnadata_extend_m(vnadata_internal_t *vdip, int new_m_allocation);

/* _vnadata_extend_f: extend the frequency allocation */
extern int _vnadata_extend_f(vnadata_internal_t *vdip, int new_f_allocation);

/* _vnadata_convert_to_z0: convert from simple z0 to frequency-dependent z0 */
extern int _vnadata_convert_to_fz0(vnadata_internal_t *vdip);

/* _vnadata_convert_to_z0: convert from frequency-dependent z0 to simple z0 */
extern int _vnadata_convert_to_z0(vnadata_internal_t *vdip);

/* _vnadata_format_to_name: return the given format descriptor as a string */
extern const char *_vnadata_format_to_name(const vnadata_format_descriptor_t
	*format);

/* _vnadata_update_format_string: recompute vdi_format_string */
extern int _vnadata_update_format_string(vnadata_internal_t *vdip);

/* _vnadata_parse_filename: try to determine the filetype from the filename */
extern vnadata_filetype_t _vnadata_parse_filename(const char *filename,
	int *ports);

/* _vnadata_set_simple_format: set a single parameter */
extern int _vnadata_set_simple_format(vnadata_internal_t *vdip,
	vnadata_parameter_type_t type, vnadata_format_t format);

/* _vnadata_load_npd: load a NPD format file */
extern int _vnadata_load_npd(vnadata_internal_t *vdip, FILE *fp,
	const char *filename);

/* _vnadata_load_touchstone: load a touchstone file */
extern int _vnadata_load_touchstone(vnadata_internal_t *vdip, FILE *fp,
	const char *filename);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* VNADATA_INTERNAL_H */
