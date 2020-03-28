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

#ifndef _VNAPROPERTY_H
#define _VNAPROPERTY_H

#include <stdarg.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum vnaproperty_type {
    VNAPROPERTY_ERROR	= -1,
    VNAPROPERTY_SCALAR	= 0x56505253,	/* "VPRS" */
    VNAPROPERTY_MAP	= 0x5650524D,	/* "VPRM" */
    VNAPROPERTY_LIST	= 0x5650524C	/* "VPRL" */
} vnaproperty_type_t;

/*
 * vnaproperty_t: opaque scalar, map or list object
 */
typedef struct vnaproperty vnaproperty_t;

/*
 * vnaproperty_map_item_t: opaque
 */
typedef struct vnaproperty_map_item vnaproperty_map_item_t;

/*
 * vnaproperty_map_pair: key-value pair for vnaproperty_map_get_pairs
 */
typedef struct vnaproperty_map_pair {
    const char *vmpr_key;
    vnaproperty_t *vmpr_value;
} vnaproperty_map_pair_t;

/*
 * vnaproperty_type: return the type of the given element
 *   @element: object to test
 */
extern vnaproperty_type_t vnaproperty_type(const vnaproperty_t *element);

/*
 * vnaproperty_hold: increment the reference count on element
 *   @element: object to reference
 */
void vnaproperty_hold(vnaproperty_t *element);

/*
 * vnaproperty_free: decrement the reference count on element and free if zero
 *   @element: object to release
 */
void vnaproperty_free(vnaproperty_t *element);

/*
 * vnaproperty_alloc_scalar: allocate a scalar object
 *   @value: value of the scalar (string)
 */
extern vnaproperty_t *vnaproperty_scalar_alloc(const char *value);

/*
 * vnaproperty_scalar_get: return the value of a scalar
 *   @scalar: scalar object
 */
extern const char *vnaproperty_scalar_get(const vnaproperty_t *scalar);

/*
 * vmaproperty_scalar_set: change the value of a scalar object
 *   @scalar: scalar object
 *   @value: new value
 */
extern int vmaproperty_scalar_set(vnaproperty_t *scalar,
	const char *value);

/*
 * vnaproperty_alloc_list: allocate a new list object
 */
extern vnaproperty_t *vnaproperty_list_alloc();

/*
 * vnaproperty_list_count: return the number of items in the list
 *   @list: list object
 */
extern int vnaproperty_list_count(const vnaproperty_t *list);

/*
 * vnaproperty_list_get: get the element at the given index
 *   @list: list object
 *   @index: index of element to get
 *
 * Note: this function doesn't increment the reference count on the
 * retrieved element.
 */
extern vnaproperty_t *vnaproperty_list_get(const vnaproperty_t *list,
	int index);

/*
 * vnaproperty_list_set: replace the element at the given index
 *   @list: list object
 *   @index: index where element should be placed
 *   @element: element to insert (reference is transferred to list)
 *
 * Index must be in 0..N where N is the number of elements in the list.
 */
extern int vnaproperty_list_set(vnaproperty_t *list, int index,
	vnaproperty_t *element);

/*
 * vnaproperty_list_append: append element to the list
 *   @list: list object
 *   @element: object to append (reference is transferred to list)
 */
extern int vnaproperty_list_append(vnaproperty_t *list, vnaproperty_t *element);

/*
 * vnaproperty_list_insert: insert element at the given index
 *   @list: list object
 *   @index: index where element should be inserted
 *   @element: element to insert (reference is transferred to list)
 *
 * Index must be in 0..N where N is the number of elements in the list.
 */
extern int vnaproperty_list_insert(vnaproperty_t *list, int index,
	vnaproperty_t *element);

/*
 * vnaproperty_list_delete: delete the element at the given index
 *   @list: list object
 *   @index: index of element to delete
 */
extern int vnaproperty_list_delete(vnaproperty_t *list, int index);

/*
 * vnaproperty_alloc_map: allocate a new map object
 */
extern vnaproperty_t *vnaproperty_map_alloc();

/*
 * vnaproperty_map_count: return the number of items in the map
 *   @map: map object
 */
extern int vnaproperty_map_count(const vnaproperty_t *map);

/*
 * vnaproperty_map_get: get the element with given key
 */
extern vnaproperty_t *vnaproperty_map_get(const vnaproperty_t *map,
	const char *key);

/*
 * vnaproperty_map_set: add an element to the map (replacing if key exists)
 *   @map: map object
 *   @key: search key
 *   @element: element to add (reference is transferred to list)
 */
extern int vnaproperty_map_set(vnaproperty_t *map, const char *key,
	vnaproperty_t *element);

/*
 * vnaproperty_map_delete: delete the element with given key
 *   @map: map object
 *   @key: search key
 */
extern int vnaproperty_map_delete(vnaproperty_t *map, const char *key);

/*
 * vnaproperty_map_begin: begin iteration
 *   @map: map object
 */
extern const vnaproperty_map_pair_t *
vnaproperty_map_begin(const vnaproperty_t *map);

/*
 * vnaproperty_map_next: return the next key-value pair
 *   @cur: current pair
 */
extern const vnaproperty_map_pair_t *vnaproperty_map_next(
    const vnaproperty_map_pair_t *cur);

/*
 * vnaproperty_expr
 */

extern vnaproperty_type_t vnaproperty_expr_vtype(const vnaproperty_t *root,
	const char *format, va_list ap);
extern int vnaproperty_expr_vcount(const vnaproperty_t *root,
	const char *format, va_list ap);
extern const char **vnaproperty_expr_vkeys(const vnaproperty_t *root,
	const char *format, va_list ap);
extern const char *vnaproperty_expr_vget(const vnaproperty_t *root,
	const char *format, va_list ap);
extern int vnaproperty_expr_vset(vnaproperty_t **anchor,
	const char *format, va_list ap);
extern int vnaproperty_expr_vdelete(vnaproperty_t **anchor,
	const char *format, va_list ap);

/*
 * vnaproperty_expr_type: get the type of the given property expression
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 */
extern vnaproperty_type_t vnaproperty_expr_type(const vnaproperty_t *root,
	const char *format, ...);

/*
 * vnaproperty_expr_count: return count of objects in given collection
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 */
extern int vnaproperty_expr_count(const vnaproperty_t *root,
	const char *format, ...);

/*
 * vnaproperty_expr_keys: return a vector of keys for the given map expr
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 *
 * Caller can free the vector by a call to free.
 */
extern const char **vnaproperty_expr_keys(const vnaproperty_t *root,
	const char *format, ...);

/*
 * vnaproperty_expr_get: get a property value from a property expression
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 */
extern const char *vnaproperty_expr_get(const vnaproperty_t *root,
	const char *format, ...);

/*
 * vnaproperty_expr_set: set a property value from a property expression
 *   @anchor: address of root property pointer
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 */
extern int vnaproperty_expr_set(vnaproperty_t **anchor,
	const char *format, ...);

/*
 * vnaproperty_expr_delete: delete the value described by format
 *   @anchor: address of root property pointer
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 */
extern int vnaproperty_expr_delete(vnaproperty_t **anchor,
	const char *format, ...);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _VNAPROPERTY_H */
