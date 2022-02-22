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

#include <assert.h>
#include <errno.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"


/*
 * _get_property_root: return the address of the appropriate property root
 */
static vnaproperty_t **_get_property_root(vnacal_t *vcp, int ci)
{
    vnacal_calibration_t *calp;

    if (vcp == NULL || vcp->vc_magic != VC_MAGIC) {
	errno = EINVAL;
	return NULL;
    }
    if (ci == -1) {
	return &vcp->vc_properties;
    }
    if ((calp = _vnacal_get_calibration(vcp, ci)) != NULL) {
	return &calp->cal_properties;
    }
    errno = EINVAL;
    return NULL;
}

/*
 * vnacal_property_type: get the type of the given property expression
 *   @root:   property data root (can be NULL)
 *   @format: printf format string forming the property expression
 *   @...:    optional variable arguments
 */
int vnacal_property_type(vnacal_t *vcp, int ci,
	const char *format, ...)
{
    va_list ap;
    vnaproperty_t **anchor;
    int rv;

    if ((anchor = _get_property_root(vcp, ci)) == NULL) {
	return -1;
    }
    va_start(ap, format);
    rv = vnaproperty_vtype(*anchor, format, ap);
    va_end(ap);

    return rv;
}

/*
 * vnacal_property_count: return count of elements in given collection
 *   @root:   property data root (can be NULL)
 *   @format: printf format string forming the property expression
 *   @...:    optional variable arguments
 */
int vnacal_property_count(vnacal_t *vcp, int ci,
	const char *format, ...)
{
    va_list ap;
    vnaproperty_t **anchor;
    int count;

    if ((anchor = _get_property_root(vcp, ci)) == NULL) {
	return -1;
    }
    va_start(ap, format);
    count = vnaproperty_vcount(*anchor, format, ap);
    va_end(ap);

    return count;
}

/*
 * vnacal_property_keys: return a vector of keys for the given map expr
 *   @root:   property data root (can be NULL)
 *   @format: printf format string forming the property expression
 *   @...:    optional variable arguments
 *
 * Caller can free the vector by a call to free.
 */
const char **vnacal_property_keys(vnacal_t *vcp, int ci,
	const char *format, ...)
{
    va_list ap;
    vnaproperty_t **anchor;
    const char **keys;

    if ((anchor = _get_property_root(vcp, ci)) == NULL) {
	return NULL;
    }
    va_start(ap, format);
    keys = vnaproperty_vkeys(*anchor, format, ap);
    va_end(ap);

    return keys;
}

/*
 * vnacal_property_get: get a property value from a property expression
 *   @root:   property data root (can be NULL)
 *   @format: printf format string forming the property expression
 *   @...:    optional variable arguments
 */
const char *vnacal_property_get(vnacal_t *vcp, int ci,
	const char *format, ...)
{
    va_list ap;
    vnaproperty_t **anchor;
    const char *value;

    if ((anchor = _get_property_root(vcp, ci)) == NULL) {
	return NULL;
    }
    va_start(ap, format);
    value = vnaproperty_vget(*anchor, format, ap);
    va_end(ap);

    return value;
}

/*
 * vnacal_property_set: set a property value from a property expression
 *   @anchor: address of root property pointer
 *   @format: printf format string forming the property expression
 *   @...:    optional variable arguments
 */
int vnacal_property_set(vnacal_t *vcp, int ci,
	const char *format, ...)
{
    va_list ap;
    vnaproperty_t **anchor;
    int rv;

    if ((anchor = _get_property_root(vcp, ci)) == NULL) {
	return -1;
    }
    va_start(ap, format);
    rv = vnaproperty_vset(anchor, format, ap);
    va_end(ap);

    return rv;
}

/*
 * vnacal_property_delete: delete the value described by format
 *   @anchor: address of root property pointer
 *   @format: printf format string forming the property expression
 *   @...:    optional variable arguments
 */
int vnacal_property_delete(vnacal_t *vcp, int ci,
	const char *format, ...)
{
    va_list ap;
    vnaproperty_t **anchor;
    int rv;

    if ((anchor = _get_property_root(vcp, ci)) == NULL) {
	return -1;
    }
    va_start(ap, format);
    rv = vnaproperty_vdelete(anchor, format, ap);
    va_end(ap);

    return rv;
}

/*
 * vnacal_property_get_subtree: get the subtree described by format
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 *   @format:  printf-like format string forming the property expression
 *   @...:     optional variable arguments
 */
vnaproperty_t *vnacal_property_get_subtree(vnacal_t *vcp, int ci,
	const char *format, ...)
{
    va_list ap;
    vnaproperty_t **anchor;
    vnaproperty_t *subtree;

    if ((anchor = _get_property_root(vcp, ci)) == NULL) {
	return NULL;
    }
    va_start(ap, format);
    subtree = vnaproperty_vget_subtree(*anchor, format, ap);
    va_end(ap);

    return subtree;
}

/*
 * vnacal_property_set_subtree: for subtree and return address
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 *   @format:  printf-like format string forming the property expression
 *   @...:     optional variable arguments
 */
vnaproperty_t **vnacal_property_set_subtree(vnacal_t *vcp, int ci,
	const char *format, ...)
{
    va_list ap;
    vnaproperty_t **anchor;
    vnaproperty_t **subtree;

    if ((anchor = _get_property_root(vcp, ci)) == NULL) {
	return NULL;
    }
    va_start(ap, format);
    subtree = vnaproperty_vset_subtree(anchor, format, ap);
    va_end(ap);

    return subtree;
}
