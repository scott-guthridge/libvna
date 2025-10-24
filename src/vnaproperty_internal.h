/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2023 D Scott Guthridge <scott_guthridge@rompromity.net>
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

typedef enum {
    VNAPROPERTY_ERROR	= -1,
    VNAPROPERTY_SCALAR	= 0x56505253,	/* "VPRS" */
    VNAPROPERTY_MAP	= 0x5650524D,	/* "VPRM" */
    VNAPROPERTY_LIST	= 0x5650524C	/* "VPRL" */
} vnaproperty_type_t;

#define VNAPROPERTY_MAP_PAIR_ELEMENT_MAGIC	0x564D5045	/* "VMPE" */

/* forward */
struct vnaproperty_map;

/*
 * vnaproperty_t: base class for VNA properties
 */
struct vnaproperty {
    uint32_t vpr_type;
    int vpr_line;	/* line number if imported from file */
};

/*
 * vnaproperty_scalar_t: scalar (string) structure
 */
typedef struct vnaproperty_scalar {
    vnaproperty_t vps_base;
    char *vps_value;
} vnaproperty_scalar_t;

/*
 * vnaproperty_list_t: property list structure
 */
typedef struct vnaproperty_list {
    vnaproperty_t vpl_base;
    size_t vpl_length;
    size_t vpl_allocation;
    vnaproperty_t **vpl_vector;
} vnaproperty_list_t;

/*
 * vnaproperty_map_pair: key-value pair for _vnaproperty_map_get_pairs
 */
typedef struct vnaproperty_map_pair {
    const char *vmpr_key;
    vnaproperty_t *vmpr_value;
} vnaproperty_map_pair_t;

/*
 * vnaproperty_map_element_t: internal element of a map structure
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
 * vnaproperty_map_t: property map structure
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

/*
 * vnaproperty_yaml_t: common argument structure for yaml import/export
 */
typedef struct vnaproperty_yaml {
    void	       *vyml_document;	/* yaml_document_t */
    const char         *vyml_filename;	/* filename for error messages */
    vnaerr_error_fn_t  *vyml_error_fn;	/* error reporting function */
    void	       *vyml_error_arg;	/* argument to error function */
} vnaproperty_yaml_t;

/* _vnaproperty_yaml_import: import properties from a YAML document */
extern int _vnaproperty_yaml_import(vnaproperty_yaml_t *vymlp,
	vnaproperty_t **rootptr, void *yaml_node);

/* _vnaproperty_yaml_export: export properties to a YAML document */
extern int _vnaproperty_yaml_export(vnaproperty_yaml_t *vymlp,
	const vnaproperty_t *root);

/* _vnaproperty_yaml_error: report an error */
extern void _vnaproperty_yaml_error(const vnaproperty_yaml_t *vymlp,
	vnaerr_category_t category, const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 3, 4)))
#endif
;

/* _vnaproperty_get_line: return the line number where a node was parsed */
extern int _vnaproperty_get_line(const vnaproperty_t *vprp);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _VNAPROPERTY_INTERNAL_H */
