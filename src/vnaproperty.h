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

#ifndef _VNAPROPERTY_H
#define _VNAPROPERTY_H

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <vnaerr.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * vnaproperty_t: opaque scalar, map or list structure
 */
typedef struct vnaproperty vnaproperty_t;

/*
 * vnaproperty_vtype: get the type of the given property expression
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @ap:     variable arguments
 *
 * Return:
 *   'm' map
 *   'l' list
 *   's' scalar
 *   -1  error
 */
extern int vnaproperty_vtype(const vnaproperty_t *root,
	const char *format, va_list ap)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 2, 0)))
#endif
;

/*
 * vnaproperty_vcount: return count of elements in given collection
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @ap:     variable arguments
 */
extern int vnaproperty_vcount(const vnaproperty_t *root,
	const char *format, va_list ap)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 2, 0)))
#endif
;

/*
 * vnaproperty_vkeys: return a vector of keys for the given map expr
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @ap:     variable arguments
 *
 * Caller can free the vector by a call to free.
 */
extern const char **vnaproperty_vkeys(const vnaproperty_t *root,
	const char *format, va_list ap)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 2, 0)))
#endif
;

/*
 * vnaproperty_vget: get a property value from a property expression
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @ap:     variable arguments
 */
extern const char *vnaproperty_vget(const vnaproperty_t *root,
	const char *format, va_list ap)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 2, 0)))
#endif
;

/*
 * vnaproperty_vset: set a property value from a property expression
 *   @rootptr: address of root property pointer
 *   @format:  printf-like format string forming the property expression
 *   @ap:      variable arguments
 */
extern int vnaproperty_vset(vnaproperty_t **rootptr,
	const char *format, va_list ap)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 2, 0)))
#endif
;

/*
 * vnaproperty_delete: vdelete the value described by format
 *   @rootptr: address of root property pointer
 *   @format:  printf-like format string forming the property expression
 *   @ap:      variable arguments
 */
extern int vnaproperty_vdelete(vnaproperty_t **rootptr,
	const char *format, va_list ap)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 2, 0)))
#endif
;

/*
 * vnaproperty_vget_subtree: get the subtree described by format
 *   @root:    property data root (can be NULL)
 *   @format:  printf-like format string forming the property expression
 *   @ap:      variable arguments
 */
extern vnaproperty_t *vnaproperty_vget_subtree(const vnaproperty_t *root,
	const char *format, va_list ap)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 2, 0)))
#endif
;

/*
 * vnaproperty_vset_subtree: make the tree conform and return address of subtree
 *   @rootptr: address of root property pointer
 *   @format:  printf-like format string forming the property expression
 *   @ap:      variable arguments
 */
extern vnaproperty_t **vnaproperty_vset_subtree(vnaproperty_t **rootptr,
	const char *format, va_list ap)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 2, 0)))
#endif
;

/*
 * vnaproperty_type: get the type of the given property expression
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 *
 * Return:
 *   'm' map
 *   'l' list
 *   's' scalar
 *   -1  error
 */
extern int vnaproperty_type(const vnaproperty_t *root,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 2, 3)))
#endif
;

/*
 * vnaproperty_count: return count of elements in given collection
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 */
extern int vnaproperty_count(const vnaproperty_t *root,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 2, 3)))
#endif
;

/*
 * vnaproperty_keys: return a vector of keys for the given map expr
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 *
 * Caller can free the vector by a call to free.
 */
extern const char **vnaproperty_keys(const vnaproperty_t *root,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 2, 3)))
#endif
;

/*
 * vnaproperty_get: get a property value from a property expression
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 */
extern const char *vnaproperty_get(const vnaproperty_t *root,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 2, 3)))
#endif
;

/*
 * vnaproperty_set: set a property value from a property expression
 *   @rootptr: address of root property pointer
 *   @format:  printf-like format string forming the property expression
 *   @...:     optional variable arguments
 */
extern int vnaproperty_set(vnaproperty_t **rootptr,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 2, 3)))
#endif
;

/*
 * vnaproperty_delete: delete the value described by format
 *   @rootptr: address of root property pointer
 *   @format:  printf-like format string forming the property expression
 *   @...:     optional variable arguments
 */
extern int vnaproperty_delete(vnaproperty_t **rootptr,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 2, 3)))
#endif
;

/*
 * vnaproperty_get_subtree: get the subtree described by format
 *   @root:    address of root property pointer
 *   @format:  printf-like format string forming the property expression
 *   @...:     optional variable arguments
 */
extern vnaproperty_t *vnaproperty_get_subtree(const vnaproperty_t *root,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 2, 3)))
#endif
;

/*
 * vnaproperty_set_subtree: make the tree conform and return address of subtree
 *   @rootptr: address of root property pointer
 *   @format:  printf-like format string forming the property expression
 *   @...:     optional variable arguments
 */
extern vnaproperty_t **vnaproperty_set_subtree(vnaproperty_t **rootptr,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 2, 3)))
#endif
;

/*
 * vnaproperty_copy: copy a subtree
 *   @destination: subtree to be replaced by copy
 *   @source:      subtree to copy
 */
extern int vnaproperty_copy(vnaproperty_t **destination,
	const vnaproperty_t *source);

/*
 * vnaproperty_quote_key: quote a map ID that contains spaces or reserved chars
 *   @key: map key to quote
 * 
 * Caller must free the returned memory by a call to free().
 */
extern char *vnaproperty_quote_key(const char *key);

/*
 * vnaproperty_import_yaml_from_string: import YAML from a string
 *   @rootptr:   address of vnaproperty root
 *   @input:     string to parse
 *   @error_fn:  optional error reporting function
 *   @error_arg: optional argument to error reporting function
 */
extern int vnaproperty_import_yaml_from_string(vnaproperty_t **rootptr,
	const char *input, vnaerr_error_fn_t *error_fn, void *error_arg);

/*
 * vnaproperty_import_yaml_from_file: import YAML from a file pointer
 *   @rootptr:   address of vnaproperty root
 *   @fp:        open file pointer
 *   @filename:  file name for error messages
 *   @error_fn:  optional error reporting function
 *   @error_arg: optional argument to error reporting function
 */
extern int vnaproperty_import_yaml_from_file(vnaproperty_t **rootptr, FILE *fp,
	const char *filename, vnaerr_error_fn_t *error_fn, void *error_arg);

/*
 * vnaproperty_export_yaml_to_file: import YAML from a file pointer
 *   @rootptr:   address of vnaproperty root
 *   @fp:        open file pointer
 *   @filename:  file name for error messages
 *   @error_fn:  optional error reporting function
 *   @error_arg: optional argument to error reporting function
 */
extern int vnaproperty_export_yaml_to_file(const vnaproperty_t *root, FILE *fp,
	const char *filename, vnaerr_error_fn_t *error_fn, void *error_arg);

/***********************************************************************
 * ZZ: delete me
 * Undocumented YAML Import / Export
 *
 * The _vnaproperty_yaml_import and _vnaproperty_yaml_export functions
 * build property trees from YAML subtrees, and build YAML subtrees from
 * a property trees, respectively.  These are public so that they can
 * be used in libraries that are friendly with this one.  They remain
 * undocumented, however, to avoid creating dependencies on libyaml in
 * our external API, making it possible to switch to newer versions of
 * libyaml or to other YAML libraries without breaking compatibility.
 * Note that we pass pointers to the libyaml types as pointer to void
 * so that this file doesn't have to include libyaml.h.
 *
 **********************************************************************/

/*
 * vnaproperty_yaml_t: common argument structure for yaml import/export
 */
typedef struct vnaproperty_yaml {
    void	       *vyml_document;	/* yaml_document_t */
    const char         *vyml_filename;	/* filename for error messages */
    vnaerr_error_fn_t  *vyml_error_fn;	/* error reporting function */
    void	       *vyml_error_arg;	/* argument to error function */
} vnaproperty_yaml_t;

/*
 * _vnaproperty_yaml_import: import properties from a YAML document
 *   @vymlp:     common argument structure
 *   @rootptr:   address of vnaproperty root
 *   @yaml_node: YAML node type cast to void pointer
 */
extern int _vnaproperty_yaml_import(vnaproperty_yaml_t *vymlp,
	vnaproperty_t **rootptr, void *yaml_node);

/*
 * _vnaproperty_yaml_export: export properties to a YAML document
 *   @vymlp:     common argument structure
 *   @root:      address of vnaproperty root
 *
 * Returns the index of a YAML node, or -1 on error.
 */
extern int _vnaproperty_yaml_export(vnaproperty_yaml_t *vymlp,
	const vnaproperty_t *root);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _VNAPROPERTY_H */
