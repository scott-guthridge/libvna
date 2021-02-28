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

#include <stddef.h>
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
 * vnadata_internal_t: internal version of vnadata_t
 */
typedef struct vnadata_internal {
    uint32_t vdi_magic;
    uint32_t vdi_flags;			/* bitwize-OR of vnadata_flags */
    vnadata_t vdi_vd;
    int vdi_p_allocation;
    int vdi_f_allocation;
    int vdi_m_allocation;
    union {
	double complex *vdi_z0_vector;		/* frequency-independent */
	double complex **vdi_z0_vector_vector;	/* frequency-dependent */
    } vdi_z0;
} vnadata_internal_t;

#define vdi_z0_vector		vdi_z0.vdi_z0_vector
#define vdi_z0_vector_vector	vdi_z0.vdi_z0_vector_vector

/*
 * VDP_TO_VDIP: convert from pointer to vnadata_t to vnadata_internal_t
 */
#define VDP_TO_VDIP(vdp) \
    (vnadata_internal_t *)((char *)(vdp) - offsetof(vnadata_internal_t, vdi_vd))

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

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* VNADATA_INTERNAL_H */
