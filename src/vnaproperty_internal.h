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

#ifndef _VNAPROPERTY_INTERNAL_H
#define _VNAPROPERTY_INTERNAL_H

#include <stdint.h>
#include "vnaproperty.h"

#ifdef __cplusplus
extern "C" {
#endif

#define VNAPROPERTY_MAP_PAIR_ELEMENT_MAGIC	0x564D5045	/* "VMPE" */

/* forward */
struct vnaproperty_map;

/*
 * vnaproperty_t: base class for VNA properties
 */
struct vnaproperty {
    uint32_t vpr_type;
    uint32_t vpr_refcount;
};

/*
 * vnaproperty_scalar_t: scalar (string) object
 */
typedef struct vnaproperty_scalar {
    vnaproperty_t vps_base;
    char *vps_value;
} vnaproperty_scalar_t;

/*
 * vnaproperty_list_t: list object
 */
typedef struct vnaproperty_list {
    vnaproperty_t vpl_base;
    size_t vpl_length;
    size_t vpl_allocation;
    vnaproperty_t **vpl_vector;
} vnaproperty_list_t;

/*
 * vnaproperty_map_element_t: internal element of a map object
 */
typedef struct vnaproperty_map_element {
    vnaproperty_map_pair_t vme_pair;
    uint32_t vme_magic;
    uint32_t vme_hashval;
    struct vnaproperty_map_element *vme_hash_next;
    struct vnaproperty_map_element *vme_order_next;
    struct vnaproperty_map_element *vme_order_prev;
} vnaproperty_map_element_t;

/*
 * vnaproperty_map_t: map object
 *
 * Note:
 *   For cosmetic reasons, we remember the insertion order, and the
 *   iterator returns elements in this order so that related fields
 *   may remain grouped together.  However, there is no intended
 *   semantic significance to the map order.
 */
typedef struct vnaproperty_map {
    vnaproperty_t vpm_base;
    size_t vpm_count;
    size_t vpm_hash_size;
    vnaproperty_map_element_t **vpm_hash_table;
    vnaproperty_map_element_t *vpm_order_head;
    vnaproperty_map_element_t *vpm_order_tail;
} vnaproperty_map_t;


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _VNAPROPERTY_INTERNAL_H */
