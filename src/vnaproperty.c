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
#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnaproperty_internal.h"


/*
 * Castagnoli CRC-32 Table
 */
static uint32_t crc32c_table[] = {
	0x00000000, 0x1edc6f41, 0x3db8de82, 0x2364b1c3,
	0x7b71bd04, 0x65add245, 0x46c96386, 0x58150cc7,
	0xf6e37a08, 0xe83f1549, 0xcb5ba48a, 0xd587cbcb,
	0x8d92c70c, 0x934ea84d, 0xb02a198e, 0xaef676cf,
	0xf31a9b51, 0xedc6f410, 0xcea245d3, 0xd07e2a92,
	0x886b2655, 0x96b74914, 0xb5d3f8d7, 0xab0f9796,
	0x05f9e159, 0x1b258e18, 0x38413fdb, 0x269d509a,
	0x7e885c5d, 0x6054331c, 0x433082df, 0x5deced9e,
	0xf8e959e3, 0xe63536a2, 0xc5518761, 0xdb8de820,
	0x8398e4e7, 0x9d448ba6, 0xbe203a65, 0xa0fc5524,
	0x0e0a23eb, 0x10d64caa, 0x33b2fd69, 0x2d6e9228,
	0x757b9eef, 0x6ba7f1ae, 0x48c3406d, 0x561f2f2c,
	0x0bf3c2b2, 0x152fadf3, 0x364b1c30, 0x28977371,
	0x70827fb6, 0x6e5e10f7, 0x4d3aa134, 0x53e6ce75,
	0xfd10b8ba, 0xe3ccd7fb, 0xc0a86638, 0xde740979,
	0x866105be, 0x98bd6aff, 0xbbd9db3c, 0xa505b47d,
	0xef0edc87, 0xf1d2b3c6, 0xd2b60205, 0xcc6a6d44,
	0x947f6183, 0x8aa30ec2, 0xa9c7bf01, 0xb71bd040,
	0x19eda68f, 0x0731c9ce, 0x2455780d, 0x3a89174c,
	0x629c1b8b, 0x7c4074ca, 0x5f24c509, 0x41f8aa48,
	0x1c1447d6, 0x02c82897, 0x21ac9954, 0x3f70f615,
	0x6765fad2, 0x79b99593, 0x5add2450, 0x44014b11,
	0xeaf73dde, 0xf42b529f, 0xd74fe35c, 0xc9938c1d,
	0x918680da, 0x8f5aef9b, 0xac3e5e58, 0xb2e23119,
	0x17e78564, 0x093bea25, 0x2a5f5be6, 0x348334a7,
	0x6c963860, 0x724a5721, 0x512ee6e2, 0x4ff289a3,
	0xe104ff6c, 0xffd8902d, 0xdcbc21ee, 0xc2604eaf,
	0x9a754268, 0x84a92d29, 0xa7cd9cea, 0xb911f3ab,
	0xe4fd1e35, 0xfa217174, 0xd945c0b7, 0xc799aff6,
	0x9f8ca331, 0x8150cc70, 0xa2347db3, 0xbce812f2,
	0x121e643d, 0x0cc20b7c, 0x2fa6babf, 0x317ad5fe,
	0x696fd939, 0x77b3b678, 0x54d707bb, 0x4a0b68fa,
	0xc0c1d64f, 0xde1db90e, 0xfd7908cd, 0xe3a5678c,
	0xbbb06b4b, 0xa56c040a, 0x8608b5c9, 0x98d4da88,
	0x3622ac47, 0x28fec306, 0x0b9a72c5, 0x15461d84,
	0x4d531143, 0x538f7e02, 0x70ebcfc1, 0x6e37a080,
	0x33db4d1e, 0x2d07225f, 0x0e63939c, 0x10bffcdd,
	0x48aaf01a, 0x56769f5b, 0x75122e98, 0x6bce41d9,
	0xc5383716, 0xdbe45857, 0xf880e994, 0xe65c86d5,
	0xbe498a12, 0xa095e553, 0x83f15490, 0x9d2d3bd1,
	0x38288fac, 0x26f4e0ed, 0x0590512e, 0x1b4c3e6f,
	0x435932a8, 0x5d855de9, 0x7ee1ec2a, 0x603d836b,
	0xcecbf5a4, 0xd0179ae5, 0xf3732b26, 0xedaf4467,
	0xb5ba48a0, 0xab6627e1, 0x88029622, 0x96def963,
	0xcb3214fd, 0xd5ee7bbc, 0xf68aca7f, 0xe856a53e,
	0xb043a9f9, 0xae9fc6b8, 0x8dfb777b, 0x9327183a,
	0x3dd16ef5, 0x230d01b4, 0x0069b077, 0x1eb5df36,
	0x46a0d3f1, 0x587cbcb0, 0x7b180d73, 0x65c46232,
	0x2fcf0ac8, 0x31136589, 0x1277d44a, 0x0cabbb0b,
	0x54beb7cc, 0x4a62d88d, 0x6906694e, 0x77da060f,
	0xd92c70c0, 0xc7f01f81, 0xe494ae42, 0xfa48c103,
	0xa25dcdc4, 0xbc81a285, 0x9fe51346, 0x81397c07,
	0xdcd59199, 0xc209fed8, 0xe16d4f1b, 0xffb1205a,
	0xa7a42c9d, 0xb97843dc, 0x9a1cf21f, 0x84c09d5e,
	0x2a36eb91, 0x34ea84d0, 0x178e3513, 0x09525a52,
	0x51475695, 0x4f9b39d4, 0x6cff8817, 0x7223e756,
	0xd726532b, 0xc9fa3c6a, 0xea9e8da9, 0xf442e2e8,
	0xac57ee2f, 0xb28b816e, 0x91ef30ad, 0x8f335fec,
	0x21c52923, 0x3f194662, 0x1c7df7a1, 0x02a198e0,
	0x5ab49427, 0x4468fb66, 0x670c4aa5, 0x79d025e4,
	0x243cc87a, 0x3ae0a73b, 0x198416f8, 0x075879b9,
	0x5f4d757e, 0x41911a3f, 0x62f5abfc, 0x7c29c4bd,
	0xd2dfb272, 0xcc03dd33, 0xef676cf0, 0xf1bb03b1,
	0xa9ae0f76, 0xb7726037, 0x9416d1f4, 0x8acabeb5
};

/*
 * crc32c: compute Castagnoli CRC-32 on data
 *   @value: initial value (typically -1)
 *   @data: data on which to compute the CRC
 *   @length: length of data in bytes
 */
static uint32_t crc32c(uint32_t value, const void *data, size_t length)
{
	const uint8_t *bp;

	for (bp = data; bp < &((uint8_t *)data)[length]; ++bp) {
		value = (value << 8) ^ crc32c_table[(value >> 24) ^ *bp];
	}
	return value;
}


/*
 * vnaproperty_type: return the type of the given element
 *   @element: object to test
 */
vnaproperty_type_t vnaproperty_type(const vnaproperty_t *element)
{
    return element->vpr_type;
}

/*
 * vnaproperty_hold: increment the reference count on element
 *   @element: object to reference
 */
void vnaproperty_hold(vnaproperty_t *element)
{
    assert(element->vpr_type == VNAPROPERTY_SCALAR ||
           element->vpr_type == VNAPROPERTY_LIST   ||
	   element->vpr_type == VNAPROPERTY_MAP);
    ++element->vpr_refcount;
}

/*
 * vnaproperty_free: decrement the reference count on element and free if zero
 *   @element: object to release
 */
void vnaproperty_free(vnaproperty_t *element)
{
    if (element == NULL)
	return;

    assert(element->vpr_type == VNAPROPERTY_SCALAR ||
           element->vpr_type == VNAPROPERTY_LIST   ||
	   element->vpr_type == VNAPROPERTY_MAP);
    assert(element->vpr_refcount > 0);
    if (--element->vpr_refcount == 0) {
	switch (element->vpr_type) {
	case VNAPROPERTY_SCALAR:
	    {
		vnaproperty_scalar_t *vpsp = (vnaproperty_scalar_t *)element;
		free((void *)vpsp->vps_value);
		(void)memset((void *)vpsp, 'X', sizeof(*vpsp));
	    }
	    break;

	case VNAPROPERTY_LIST:
	    {
		vnaproperty_list_t *vplp  = (vnaproperty_list_t *)element;
		for (size_t s = 0; s < vplp->vpl_length; ++s) {
		    vnaproperty_free((vnaproperty_t *)vplp->vpl_vector[s]);
		}
		free((void *)vplp->vpl_vector);
		(void)memset((void *)vplp, 'X', sizeof(*vplp));
	    }
	    break;

	case VNAPROPERTY_MAP:
	    {
		vnaproperty_map_t *vpmp = (vnaproperty_map_t *)element;
		vnaproperty_map_element_t *vmep, *next;

		for (vmep = vpmp->vpm_order_head; vmep != NULL; vmep = next) {
		    next = vmep->vme_order_next;
		    free((void *)vmep->vme_pair.vmpr_key);
		    vnaproperty_free((vnaproperty_t *)vmep->vme_pair.
			    vmpr_value);
		    memset((void *)vmep, 'X', sizeof(*vmep));
		    free((void *)vmep);
		    vmep = NULL;
		}
		free((void *)vpmp->vpm_hash_table);
		(void)memset((void *)vpmp, 'X', sizeof(*vpmp));
	    }
	    break;
	}
	free((void *)element);
    }
}

/*
 * vnaproperty_alloc_scalar: allocate a scalar object
 *   @value: value of the scalar (string)
 */
vnaproperty_t *vnaproperty_scalar_alloc(const char *value)
{
    char *copy;
    vnaproperty_scalar_t *vpsp;

    if (strcasecmp(value, "null") == 0) {
	value = "~";
    }
    if ((copy = strdup(value)) == NULL) {
	return NULL;
    }
    if ((vpsp = malloc(sizeof(vnaproperty_scalar_t))) == NULL) {
	free((void *)copy);
	return NULL;
    }
    (void)memset((void *)vpsp, 0, sizeof(*vpsp));
    vpsp->vps_base.vpr_type = VNAPROPERTY_SCALAR;
    vpsp->vps_base.vpr_refcount = 1;
    vpsp->vps_value = copy;

    return ((vnaproperty_t *)vpsp);
}

/*
 * vnaproperty_scalar_get: return the value of a scalar
 *   @scalar: scalar object
 */
const char *vnaproperty_scalar_get(const vnaproperty_t *scalar)
{
    vnaproperty_scalar_t *vpsp;

    if (scalar->vpr_type != VNAPROPERTY_SCALAR) {
	errno = EINVAL;
	return NULL;
    }
    vpsp = (vnaproperty_scalar_t *)scalar;
    return vpsp->vps_value;
}

/*
 * vmaproperty_scalar_set: change the value of a scalar object
 *   @scalar: scalar object
 *   @value: new value
 */
int vmaproperty_scalar_set(vnaproperty_t *scalar, const char *value)
{
    char *copy;
    vnaproperty_scalar_t *vpsp;

    if (scalar->vpr_type != VNAPROPERTY_SCALAR) {
	errno = EINVAL;
	return -1;
    }
    vpsp = (vnaproperty_scalar_t *)scalar;

    if (strcasecmp(value, "null") == 0) {
	value = "~";
    }
    if ((copy = strdup(value)) == NULL) {
	return -1;
    }
    free((void *)vpsp->vps_value);
    vpsp->vps_value = copy;
    return 0;
}

/*
 * _vnaproperty_list_check_allocation: extend allocation to >= size
 *   @vplp: list object
 *   @size: needed size
 */
int _vnaproperty_list_check_allocation(vnaproperty_list_t *vplp, size_t size)
{
    size_t new_allocation;
    vnaproperty_t **new_vector = NULL;

    /* If there's already enough, return. */
    if (size <= vplp->vpl_allocation)
	return 0;

    /*
     * Find the next power of 2 allocation higher than size.
     */
    new_allocation = MAX(8, vplp->vpl_allocation);
    while (new_allocation <= size) {
	new_allocation <<= 1;
	if (new_allocation == 0) {
	    errno = EINVAL;
	    return -1;
	}
    }

    /* Realloc */
    if ((new_vector = realloc(vplp->vpl_vector, new_allocation *
		    sizeof(vnaproperty_t *))) == NULL) {
	return -1;
    }
    (void)memset((void *)&new_vector[vplp->vpl_allocation], 0,
                 (new_allocation - vplp->vpl_allocation) *
		 sizeof(vnaproperty_t *));
    vplp->vpl_vector = new_vector;
    vplp->vpl_allocation = new_allocation;
    return 0;
}

/*
 * vnaproperty_alloc_list: allocate a new list object
 */
vnaproperty_t *vnaproperty_list_alloc()
{
    vnaproperty_list_t *vplp;

    if ((vplp = malloc(sizeof(vnaproperty_list_t))) == NULL) {
	return NULL;
    }
    (void)memset((void *)vplp, 0, sizeof(*vplp));
    vplp->vpl_base.vpr_type = VNAPROPERTY_LIST;
    vplp->vpl_base.vpr_refcount = 1;

    return ((vnaproperty_t *)vplp);
}

/*
 * vnaproperty_list_count: return the number of items in the list
 *   @list: list object
 */
int vnaproperty_list_count(const vnaproperty_t *list)
{
    vnaproperty_list_t *vplp;

    if (list->vpr_type != VNAPROPERTY_LIST) {
	errno = EINVAL;
	return -1;
    }
    vplp = (vnaproperty_list_t *)list;
    return vplp->vpl_length;
}

/*
 * vnaproperty_list_get: get the element at the given index
 *   @list: list object
 *   @index: index of element to get
 *
 * Note: this function doesn't increment the reference count on the
 * retrieved element.
 */
vnaproperty_t *vnaproperty_list_get(const vnaproperty_t *list, int index)
{
    vnaproperty_list_t *vplp;

    if (list->vpr_type != VNAPROPERTY_LIST) {
	errno = EINVAL;
	return NULL;
    }
    vplp = (vnaproperty_list_t *)list;
    if (index < 0 || index >= vplp->vpl_length) {
	errno = EDOM;
	return NULL;
    }
    return vplp->vpl_vector[index];
}

/*
 * vnaproperty_list_set: replace the element at the given index
 *   @list: list object
 *   @index: index where element should be placed
 *   @element: element to insert (reference is transferred to list)
 */
int vnaproperty_list_set(vnaproperty_t *list, int index,
	vnaproperty_t *element)
{
    vnaproperty_list_t *vplp;

    if (list->vpr_type != VNAPROPERTY_LIST) {
	errno = EINVAL;
	return -1;
    }
    vplp = (vnaproperty_list_t *)list;
    if (index < 0) {
	errno = EDOM;
	return -1;
    }
    if (index >= vplp->vpl_length) {	/* extend case */
	if (_vnaproperty_list_check_allocation(vplp, index + 1) == -1) {
	    return -1;
	}
	for (int i = vplp->vpl_length; i < index; ++i) { /* pad with nulls */
	    if ((vplp->vpl_vector[i] = vnaproperty_scalar_alloc("~")) == NULL) {
		return -1;
	    }
	    ++vplp->vpl_length;
	}
	assert(vplp->vpl_length == index);
	vplp->vpl_length = index + 1;
    } else {				/* replace case */
	vnaproperty_free(vplp->vpl_vector[index]);
	vplp->vpl_vector[index] = NULL;
    }
    vplp->vpl_vector[index] = element;
    return 0;
}

/*
 * vnaproperty_list_append: append element to the list
 *   @list: list object
 *   @element: object to append (reference is transferred to list)
 */
int vnaproperty_list_append(vnaproperty_t *list, vnaproperty_t *element)
{
    vnaproperty_list_t *vplp;

    if (list->vpr_type != VNAPROPERTY_LIST) {
	errno = EINVAL;
	return -1;
    }
    vplp = (vnaproperty_list_t *)list;
    if (_vnaproperty_list_check_allocation(vplp, vplp->vpl_length + 1) == -1) {
	return -1;
    }
    vplp->vpl_vector[vplp->vpl_length++] = element;
    return 0;
}

/*
 * vnaproperty_list_insert: insert element at the given index
 *   @list: list object
 *   @index: index where element should be inserted
 *   @element: element to insert (reference is transferred to list)
 *
 * Index must be in 0..N where N is the number of elements in the list.
 */
int vnaproperty_list_insert(vnaproperty_t *list, int index,
	vnaproperty_t *element)
{
    vnaproperty_list_t *vplp;

    if (list->vpr_type != VNAPROPERTY_LIST) {
	errno = EINVAL;
	return -1;
    }
    vplp = (vnaproperty_list_t *)list;
    if (index < 0) {
	errno = EDOM;
	return -1;
    }
    if (index >= vplp->vpl_length) {
	return vnaproperty_list_set(list, index, element);
    }
    if (_vnaproperty_list_check_allocation(vplp, vplp->vpl_length + 1) == -1) {
	return -1;
    }
    (void)memmove((void *)&vplp->vpl_vector[index + 1],
                  (void *)&vplp->vpl_vector[index],
		  (vplp->vpl_length - index) * sizeof(vnaproperty_t *));
    vplp->vpl_vector[index] = element;
    ++vplp->vpl_length;
    return 0;
}

/*
 * vnaproperty_list_delete: delete the element at the given index
 *   @list: list object
 *   @index: index of element to delete
 */
int vnaproperty_list_delete(vnaproperty_t *list, int index)
{
    vnaproperty_list_t *vplp;

    if (list->vpr_type != VNAPROPERTY_LIST) {
	errno = EINVAL;
	return -1;
    }
    vplp = (vnaproperty_list_t *)list;
    if (index < 0 || index >= vplp->vpl_length) {
	errno = EDOM;
	return -1;
    }
    (void)memmove((void *)&vplp->vpl_vector[index],
                  (void *)&vplp->vpl_vector[index + 1],
		  (vplp->vpl_length - index) * sizeof(vnaproperty_t *));
    vplp->vpl_vector[--vplp->vpl_length] = NULL;
    return 0;
}

/*
 * _vnaproperty_map_compare_keys: compare key and hashval to the given element
 *   @key: key string
 *   @hashval: full 32-bit hashval
 *   @vmep: element with which to compare
 *
 * compares by hashval first, then by string compare
 */
static int _vnaproperty_map_compare_keys(const char *key, uint32_t hashval,
	vnaproperty_map_element_t *vmep)
{
    if (hashval < vmep->vme_hashval)
	return -1;

    if (hashval > vmep->vme_hashval)
	return +1;

    return strcmp(key, vmep->vme_pair.vmpr_key);
}

/*
 * _vnaproperty_map_find_anchor: find insertion point
 *   @vpmp: map object
 *   @anchor:  address of address of the insertion point
 *   @key: search key
 *   @hashval: full 32 bit hash value
 *
 * If the key is found return true, else return false.  In either case,
 * *anchor is set to the insertion point, i.e. if the element is found,
 * **anchor is the found element; if not found, *anchor is the address
 * of the pointer where the element should be inserted.
 */
static bool _vnaproperty_map_find_anchor(
	vnaproperty_map_t *vpmp, vnaproperty_map_element_t ***anchor,
	const char *key, uint32_t hashval)
{
    vnaproperty_map_element_t *vmep, **vmepp;
    int bucket;
    int cmp = 1;

    assert(vpmp->vpm_hash_size != 0);
    bucket = hashval % vpmp->vpm_hash_size;
    vmepp = &vpmp->vpm_hash_table[bucket];
    for (; (vmep = *vmepp) != NULL; vmepp = &vmep->vme_hash_next) {
	if ((cmp = _vnaproperty_map_compare_keys(key, hashval, vmep)) <= 0)
	    break;
    }
    *anchor = vmepp;
    return cmp == 0;
}

static int _vnaproperty_map_expand(vnaproperty_map_t *vpmp)
{
    size_t old_allocation = vpmp->vpm_hash_size;
    size_t new_allocation = vpmp->vpm_count + 1;
    vnaproperty_map_element_t **new_table;

    /*
     * Calculate the new allocation and alloc/resize the table.
     */
    new_allocation += (new_allocation + 1) / 2;
    if (new_allocation < 11) {
	new_allocation = 11;
    }
    new_table = (vnaproperty_map_element_t **)realloc(
	    (void *)vpmp->vpm_hash_table,
	    new_allocation * sizeof(vnaproperty_map_element_t *));
    if (new_table == NULL) {
	return -1;
    }
    (void)memset((void *)&new_table[old_allocation], 0,
	    (new_allocation - old_allocation) *
	    sizeof(vnaproperty_map_element_t *));
    vpmp->vpm_hash_table = new_table;

    /*
     * Walk the old chains and rehash all elements.
     */
    for (size_t s = 0; s < old_allocation; ++s) {
	vnaproperty_map_element_t *head, *cur;

	/*
	 * Move the elements from the current hash chain to a local list
	 * at head.
	 */
	head = new_table[s];
	new_table[s] = NULL;

	/*
	 * For each element on the local list...
	 */
	while ((cur = head) != NULL) {
	    vnaproperty_map_element_t *next, **anchor;

	    /*
	     * Remove the first element as "cur".
	     */
	    head = cur->vme_hash_next;
	    cur->vme_hash_next = NULL;

	    /*
	     * Re-insert "cur" into its new hash chain.  As we go through
	     * the buckets, we will re-process elements we've already moved.
	     * But since a reprocessed element is always already in the
	     * correct location, it can be reprocessed only once.
	     */
	    anchor = &new_table[cur->vme_hashval % new_allocation];
	    for (; (next = *anchor) != NULL; anchor = &next->vme_hash_next) {
		if (_vnaproperty_map_compare_keys(cur->vme_pair.vmpr_key,
			    cur->vme_hashval, next) <= 0)
		    break;
	    }
	    cur->vme_hash_next = next;
	    *anchor = cur;
	}
    }
    vpmp->vpm_hash_size = new_allocation;
    return 0;
}

/*
 * _vnaproperty_map_append_order_element: append to order list
 *   @vpmp: map object
 *   @vmep: element to append
 */
void _vnaproperty_map_append_order_element(vnaproperty_map_t *vpmp,
					   vnaproperty_map_element_t *vmep)
{
    assert(vmep->vme_order_next == NULL);
    assert(vmep->vme_order_prev == NULL);
    if (vpmp->vpm_order_head == NULL) {
	assert(vpmp->vpm_order_tail == NULL);
	assert(vpmp->vpm_count == 0);
	vpmp->vpm_order_head = vmep;
    } else {
	assert(vpmp->vpm_order_tail != NULL);
	assert(vpmp->vpm_count > 0);
	vmep->vme_order_prev = vpmp->vpm_order_tail;
	vpmp->vpm_order_tail->vme_order_next = vmep;
    }
    vpmp->vpm_order_tail = vmep;
}

/*
 * _vnaproperty_map_delete_order_element: delete from order list
 *   @vpmp: map object
 *   @vmep: element to delete
 */
void _vnaproperty_map_delete_order_element(vnaproperty_map_t *vpmp,
					   vnaproperty_map_element_t *vmep)
{
    vnaproperty_map_element_t *prev = vmep->vme_order_prev;
    vnaproperty_map_element_t *next = vmep->vme_order_next;

    if (prev != NULL) {
	prev->vme_order_next = next;
    } else {
	vpmp->vpm_order_head = next;
    }
    if (next != NULL) {
	next->vme_order_prev = prev;
    } else {
	vpmp->vpm_order_tail = prev;
    }
}

/*
 * vnaproperty_alloc_map: allocate a new map object
 */
vnaproperty_t *vnaproperty_map_alloc()
{
    vnaproperty_map_t *vpmp;

    if ((vpmp = malloc(sizeof(vnaproperty_map_t))) == NULL) {
	return NULL;
    }
    (void)memset((void *)vpmp, 0, sizeof(*vpmp));
    vpmp->vpm_base.vpr_type = VNAPROPERTY_MAP;
    vpmp->vpm_base.vpr_refcount = 1;
    return ((vnaproperty_t *)vpmp);
}

/*
 * vnaproperty_map_count: return the number of items in the map
 *   @map: map object
 */
int vnaproperty_map_count(const vnaproperty_t *map)
{
    vnaproperty_map_t *vpmp;

    if (map->vpr_type != VNAPROPERTY_MAP) {
	errno = EINVAL;
	return -1;
    }
    vpmp = (vnaproperty_map_t *)map;
    return vpmp->vpm_count;
}

/*
 * vnaproperty_map_get: get the element with given key
 */
vnaproperty_t *vnaproperty_map_get(const vnaproperty_t *map, const char *key)
{
    vnaproperty_map_t *vpmp;
    vnaproperty_map_element_t **anchor;
    uint32_t hashval;

    if (map->vpr_type != VNAPROPERTY_MAP) {
	errno = EINVAL;
	return NULL;
    }
    vpmp = (vnaproperty_map_t *)map;
    if (vpmp->vpm_count == 0) {
	errno = ENOENT;
	return NULL;
    }
    hashval = crc32c(-1, (void *)key, strlen(key));
    if (!_vnaproperty_map_find_anchor(vpmp, &anchor, key, hashval)) {
	errno = ENOENT;
	return NULL;
    }
    return (vnaproperty_t *)(*anchor)->vme_pair.vmpr_value;
}

/*
 * vnaproperty_map_set: add an element to the map (replacing if key exists)
 *   @map: map object
 *   @key: search key
 *   @element: element to add (reference is transferred to list)
 */
int vnaproperty_map_set(vnaproperty_t *map, const char *key,
	vnaproperty_t *element)
{
    vnaproperty_map_t *vpmp;
    vnaproperty_map_element_t **anchor, *vmep;
    uint32_t hashval;

    if (map->vpr_type != VNAPROPERTY_MAP) {
	errno = EINVAL;
	return -1;
    }
    vpmp = (vnaproperty_map_t *)map;
    if (vpmp->vpm_count + 1 >= 2 * vpmp->vpm_hash_size) {
	if (_vnaproperty_map_expand(vpmp) == -1) {
	    return -1;
	}
    }
    hashval = crc32c(-1, (void *)key, strlen(key));
    if (_vnaproperty_map_find_anchor(vpmp, &anchor, key, hashval)) {
	vmep = *anchor;
	vnaproperty_free(vmep->vme_pair.vmpr_value);
	vmep->vme_pair.vmpr_value = element;
	return 0;
    }
    if ((vmep = malloc(sizeof(vnaproperty_map_element_t))) == NULL) {
	return -1;
    }
    (void)memset((void *)vmep, 0, sizeof(*vmep));
    if ((vmep->vme_pair.vmpr_key = strdup(key)) == NULL) {
	free((void *)vmep);
	return -1;
    }
    vmep->vme_pair.vmpr_value = element;
    vmep->vme_magic = VNAPROPERTY_MAP_PAIR_ELEMENT_MAGIC;
    vmep->vme_hashval = hashval;
    vmep->vme_hash_next = *anchor;
    *anchor = vmep;
    _vnaproperty_map_append_order_element(vpmp, vmep);
    ++vpmp->vpm_count;
    return 0;
}

/*
 * vnaproperty_map_delete: delete the element with given key
 *   @map: map object
 *   @key: search key
 */
int vnaproperty_map_delete(vnaproperty_t *map, const char *key)
{
    vnaproperty_map_t *vpmp;
    vnaproperty_map_element_t **anchor, *vmep;
    uint32_t hashval;

    if (map->vpr_type != VNAPROPERTY_MAP) {
	errno = EINVAL;
	return -1;
    }
    vpmp = (vnaproperty_map_t *)map;
    if (vpmp->vpm_count == 0) {
	errno = ENOENT;
	return -1;
    }
    hashval = crc32c(-1, (void *)key, strlen(key));
    if (!_vnaproperty_map_find_anchor(vpmp, &anchor, key, hashval)) {
	errno = ENOENT;
	return -1;
    }
    vmep = *anchor;
    assert(vpmp->vpm_count > 0);
    _vnaproperty_map_delete_order_element(vpmp, vmep);
    *anchor = vmep->vme_hash_next;
    free((void *)vmep->vme_pair.vmpr_key);
    vnaproperty_free((vnaproperty_t *)vmep->vme_pair.vmpr_value);
    memset((void *)vmep, 'X', sizeof(*vmep));
    free((void *)vmep);
    --vpmp->vpm_count;
    if (vpmp->vpm_count == 0) {
	assert(vpmp->vpm_order_head == NULL);
	assert(vpmp->vpm_order_tail == NULL);
    } else {
	assert(vpmp->vpm_order_head != NULL);
	assert(vpmp->vpm_order_tail != NULL);
    }
    return 0;
}

/*
 * vnaproperty_map_begin: return the first key-value pair
 *   @map: map object
 */
const vnaproperty_map_pair_t *vnaproperty_map_begin(const vnaproperty_t *map)
{
    vnaproperty_map_t *vpmp;

    if (map->vpr_type != VNAPROPERTY_MAP) {
	errno = EINVAL;
	return NULL;
    }
    vpmp = (vnaproperty_map_t *)map;
    return (vnaproperty_map_pair_t *)vpmp->vpm_order_head;
}

/*
 * vnaproperty_map_next: return the next key-value pair
 *   @cur: current position
 */
const vnaproperty_map_pair_t *vnaproperty_map_next(
	const vnaproperty_map_pair_t *cur)
{
    vnaproperty_map_element_t *vmep = (vnaproperty_map_element_t *)cur;

    if (vmep->vme_magic != VNAPROPERTY_MAP_PAIR_ELEMENT_MAGIC) {
	errno = EINVAL;
	return NULL;
    }
    return (const vnaproperty_map_pair_t *)vmep->vme_order_next;
}

typedef enum vnaproperty_token {
    T_ERROR	= -1,
    T_EOF,
    T_DOT,
    T_ASSIGN,
    T_LBRACKET,
    T_RBRACKET,
    T_ID,
    T_INT
} vnaproperty_token_t;

/*
 * vnaproperty_scanner: lexical analyzer state
 */
typedef struct vnaproperty_scanner {
    char               *vs_input;
    char               *vs_position;
    char                vs_cur;
    char               *vs_text;
    vnaproperty_token_t vs_token;
} vnaproperty_scanner_t;

/*
 * VNAPROPERTY_GETCHAR: advance to the next input character
 */
#define VNAPROPERTY_GETCHAR(vsp) \
	((vsp)->vs_cur = *++(vsp)->vs_position)

/*
 * VNAPROPERTY_ISIDCHAR: return true if c is a valid identifier character
 */
#define VNAPROPERTY_ISIDCHAR(c) \
	(isalpha(c) || isdigit(c) || !isascii(c) || (c) == '_' || \
	 (c) == '-' || (c) == '+')

/*
 * _vnaproperty_scan: return the next input token
 *   @vsp: scanner state
 */
static void _vnaproperty_scan(vnaproperty_scanner_t *vsp)
{
    for (;;) {
	switch (vsp->vs_cur) {
	case '\000':
	    vsp->vs_token = T_EOF;
	    return;

	/*
	 * Ignore whitespace.
	 */
	case '\f':
	case '\n':
	case '\r':
	case '\t':
	case '\v':
	case ' ':
	    VNAPROPERTY_GETCHAR(vsp);
	    continue;

	/*
	 * Handle simple tokens.
	 */
	case '.':
	    VNAPROPERTY_GETCHAR(vsp);
	    vsp->vs_token = T_DOT;
	    return;

	case '=':
	    VNAPROPERTY_GETCHAR(vsp);
	    vsp->vs_token = T_ASSIGN;
	    return;

	case '[':
	    VNAPROPERTY_GETCHAR(vsp);
	    vsp->vs_token = T_LBRACKET;
	    return;

	case ']':
	    VNAPROPERTY_GETCHAR(vsp);
	    vsp->vs_token = T_RBRACKET;
	    return;

	default:
	    break;
	}

	/*
	 * Integers and Identifiers
	 */
	if (VNAPROPERTY_ISIDCHAR(vsp->vs_cur)) {
	    vsp->vs_text = vsp->vs_position;
	    do {
		VNAPROPERTY_GETCHAR(vsp);
	    } while (VNAPROPERTY_ISIDCHAR(vsp->vs_cur));
	    *vsp->vs_position = '\000';
	    vsp->vs_token = T_ID;
	    return;
	}

	/*
	 * All other characters are reserved.
	 */
	vsp->vs_token = T_ERROR;
	return;
    }
}

/*
 * vnaproperty_expr_t: expression node
 */
typedef struct vnaproperty_expr {
    vnaproperty_token_t vex_token;	/* T_DOT=map, T_LBRACKET=list */
    struct vnaproperty_expr *vex_next;	/* next expression node */
    vnaproperty_t *vex_collection;	/* existing collection if exists */
    union {
	char   *vex_name;		/* key for map */
	int     vex_index;		/* index for list */
    } u;
} vnaproperty_expr_t;

/*
 * _vnaproperty_expr_free: free an expression list
 *   @head: expression list returned from _vnaproperty_expr_parse
 */
static void _vnaproperty_expr_free(vnaproperty_expr_t *head)
{

    while (head != NULL) {
	vnaproperty_expr_t *vexp = head;

	head = vexp->vex_next;
	free((void *)vexp);
    }
}

/*
 * _vnaproperty_parse: parse a property expression and return expression list
 *   @vsp:    caller allocated scanner state initialized here
 *   @anchor: location to return expression list
 *   @root:   property root (can be NULL)
 *   @format: sprintf format
 *   @ap:     variable argument pointer
 *
 *   On success, *anchor is a linked list of expression nodes
 *   representing the given expression in right-to-left order.
 *   For example, given "foo.bar[5]", the returned list would
 *   be T_LIST/index=5, T_MAP/name="bar", T_MAP/name="foo".
 *   If the indicated collections exist (starting from root), the
 *   vex_collection member points to the existing object; otherwise
 *   if the object doesn't exist or is not the requested type,
 *   vex_collection is NULL.  For get operations, the first
 *   element in the returned list is the bottom most collection
 *   containing the requested element, if one exists.  For set
 *   operations, one builds new property nodes from the bottom up
 *   until finding an expression node with a non-NULL vex_collection
 *   pointer.  A set operation on that collection completes the
 *   update.  Thus, in both cases, the reversed list is the most
 *   useful order.
 */
static int _vnaproperty_parse(vnaproperty_scanner_t *vsp,
	vnaproperty_expr_t **anchor, vnaproperty_t *root,
	const char *format, va_list ap)
{
    vnaproperty_expr_t *head = NULL, *vexp;

    /*
     * Format the VNA property expression and init the scanner.
     */
    (void)memset((void *)vsp, 0, sizeof(vnaproperty_scanner_t));
    if (vasprintf(&vsp->vs_input, format, ap) == -1) {
	return -1;
    }
    vsp->vs_position = vsp->vs_input;
    vsp->vs_cur = vsp->vs_input[0];
    _vnaproperty_scan(vsp);

    /*
     * expr    : T_DOT					     // root
     *         | T_DOT T_ID tail_opt			     // absolute path
     *         | T_DOT T_LBRACKET T_ID T_BRACKET tail_opt    // absolute path
     *         | T_ID tail_opt				     // relative path
     *         | T_LBRACKET T_ID T_BRACKET tail_opt	     // relative path
     *         ;
     */
    switch (vsp->vs_token) {
    case T_DOT:
	    _vnaproperty_scan(vsp);
	    switch (vsp->vs_token) {
	    case T_ID:
		goto map;

	    case T_LBRACKET:
		break;

	    default:
		goto accept;
	    }
	    break;

    case T_ID:
	goto map;

    case T_LBRACKET:
	break;

    default:
	goto error;
    }

    /*
     * tail_opt: %empty
     *         | tail_opt T_DOT T_ID
     *         | tail_opt T_LBRACKET T_ID T_BRACKET
     *         ;
     */
    for (;;) {
	switch (vsp->vs_token) {
	case T_DOT:
	    /*
             * tail_opt : tail_opt <here> T_DOT T_ID
	     */
	    _vnaproperty_scan(vsp);
	    if (vsp->vs_token != T_ID) {
		goto error;
	    }
	    /*FALLTHROUGH*/

	map:
	    if ((vexp = malloc(sizeof(vnaproperty_expr_t))) == NULL) {
		goto error;
	    }
	    (void)memset((void *)vexp, 0, sizeof(vnaproperty_expr_t));
	    vexp->vex_token = T_DOT;
	    vexp->u.vex_name = vsp->vs_text;
	    if (root != NULL) {
		if (root->vpr_type == VNAPROPERTY_MAP) {
		    vexp->vex_collection = root;
		    root = vnaproperty_map_get(root, vexp->u.vex_name);
		} else {
		    root = NULL;
		}
	    }
	    vexp->vex_next = head;
	    head = vexp;
	    _vnaproperty_scan(vsp);
	    continue;

	case T_LBRACKET:
	    /*
	     * tail_opt : tail_opt <here> T_LBRACKET T_ID T_RBRACKET
	     */
	    _vnaproperty_scan(vsp);
	    if (vsp->vs_token != T_ID) {
		errno = EINVAL;
		goto error;
	    }
	    {
		char *end;
		int index;

		index = strtol(vsp->vs_text, &end, 0);
		if (end == vsp->vs_text) {
		    errno = EINVAL;
		    goto error;
		}
		while (isascii(*end) && isspace(*end)) {
		    ++end;
		}
		if (*end != '\000') {
		    errno = EINVAL;
		    goto error;
		}
		if ((vexp = malloc(sizeof(vnaproperty_expr_t))) == NULL) {
		    goto error;
		}
		(void)memset((void *)vexp, 0, sizeof(vnaproperty_expr_t));
		vexp->vex_token = T_LBRACKET;
		vexp->u.vex_index = index;
		if (root != NULL) {
		    if (root->vpr_type == VNAPROPERTY_LIST) {
			vexp->vex_collection = root;
			root = vnaproperty_list_get(root, vexp->u.vex_index);
		    } else {
			root = NULL;
		    }
		}
		vexp->vex_next = head;
		head = vexp;
	    }
	    _vnaproperty_scan(vsp);
	    if (vsp->vs_token != T_RBRACKET) {
		errno = EINVAL;
		goto error;
	    }
	    _vnaproperty_scan(vsp);
	    continue;

	default:
	    break;
	}
	break;
    }
accept:
    *anchor = head;
    return 0;

error:
    _vnaproperty_expr_free(head);
    return -1;
}

/*
 * _vnaproperty_expr_getobj: return the node referended by the given expression
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @ap      variable argument pointer
 */
const vnaproperty_t *_vnaproperty_expr_getobj(const vnaproperty_t *root,
	const char *format, va_list ap)
{
    vnaproperty_expr_t *head = NULL;
    vnaproperty_scanner_t vs;
    const vnaproperty_t *property = NULL;

    /*
     * Parse the expression.
     */
    if (_vnaproperty_parse(&vs, &head, (vnaproperty_t *)root,
		format, ap) == -1) {
	goto out;
    }

    /*
     * Fail if there are extra tokens at the end.
     */
    if (vs.vs_token != T_EOF) {
	errno = EINVAL;
	goto out;
    }

    /*
     * If the expression list is '.' or empty, return the root element.
     */
    if (head == NULL) {
	if (root == NULL) {
	    errno = ENOENT;
	    goto out;
	}
	property = root;
	goto out;
    }

    /*
     * If vex_collection is NULL, it means the path leading down to
     * the property either doesn't exist or the object types don't
     * match; errno will already be set.
     */
    if (head->vex_collection == NULL) {
	goto out;
    }

    /*
     * Get the value from the bottom-most collection.
     */
    switch (head->vex_token) {
    case T_DOT:
	property = vnaproperty_map_get(head->vex_collection,
		head->u.vex_name);
	break;

    case T_LBRACKET:
	property = vnaproperty_list_get(head->vex_collection,
		head->u.vex_index);
	break;

    default:
	abort();
    }

out:
    _vnaproperty_expr_free(head);
    free((void *)vs.vs_input);
    return property;
}

/*
 * vnaproperty_expr_vtype: get the type of the given property expression
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @ap:     variable argument pointer
 */
vnaproperty_type_t vnaproperty_expr_vtype(const vnaproperty_t *root,
	const char *format, va_list ap)
{
    const vnaproperty_t *property;

    if ((property = _vnaproperty_expr_getobj(root, format, ap)) == NULL) {
	return VNAPROPERTY_ERROR;
    }
    return property->vpr_type;
}

/*
 * vnaproperty_expr_vcount: return count of objects in given collection
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @ap:     variable argument pointer
 */
int vnaproperty_expr_vcount(const vnaproperty_t *root,
	const char *format, va_list ap)
{
    const vnaproperty_t *property;

    if ((property = _vnaproperty_expr_getobj(root, format, ap)) == NULL) {
	return -1;
    }
    switch (property->vpr_type) {
    case VNAPROPERTY_MAP:
	return vnaproperty_map_count(property);

    case VNAPROPERTY_LIST:
	return vnaproperty_list_count(property);

    default:
	break;
    }
    errno = EINVAL;
    return -1;
}

/*
 * vnaproperty_expr_vkeys: return a vector of keys for the given map expr
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @ap:     variable argument pointer
 */
const char **vnaproperty_expr_vkeys(const vnaproperty_t *root,
	const char *format, va_list ap)
{
    const vnaproperty_t *property;
    const vnaproperty_map_t *map;
    const vnaproperty_map_element_t *vmep;
    const char **vector, **cpp;

    if ((property = _vnaproperty_expr_getobj(root, format, ap)) == NULL) {
	return NULL;
    }
    if (property->vpr_type != VNAPROPERTY_MAP) {
	errno = EINVAL;
	return NULL;
    }
    map = (vnaproperty_map_t *)property;
    if ((vector = calloc(map->vpm_count + 1, sizeof(char *))) == NULL) {
	return NULL;
    }
    cpp = vector;
    for (vmep = map->vpm_order_head; vmep != NULL;
	    vmep = vmep->vme_order_next) {
	*cpp++ = vmep->vme_pair.vmpr_key;
    }
    *cpp = NULL;
    return vector;
}

/*
 * vnaproperty_expr_vget: get a property value from a property expression
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @ap:     variable argument pointer
 */
const char *vnaproperty_expr_vget(const vnaproperty_t *root,
	const char *format, va_list ap)
{
    const vnaproperty_t *scalar;

    if ((scalar = _vnaproperty_expr_getobj(root, format, ap)) == NULL) {
	return NULL;
    }
    if (scalar->vpr_type != VNAPROPERTY_SCALAR) {
	errno = EINVAL;
	return NULL;
    }
    return vnaproperty_scalar_get(scalar);
}

/*
 * vnaproperty_expr_vset: set a property value from a property expression
 *   @anchor: address of root property pointer
 *   @format: printf-like format string forming the property expression
 *   @ap:     variable argument pointer
 */
int vnaproperty_expr_vset(vnaproperty_t **anchor, const char *format,
	va_list ap)
{
    vnaproperty_expr_t *head = NULL;
    vnaproperty_scanner_t vs;
    vnaproperty_t *value = NULL;
    int rv, result = -1;

    /*
     * Parse the expression.
     */
    if ((rv = _vnaproperty_parse(&vs, &head, *anchor, format, ap)) == -1)
	goto out;

    /*
     * Make sure the current token is T_ASSIGN.  Allocate a new scalar
     * property from the remaining input text.
     */
    if (vs.vs_token != T_ASSIGN) {
	errno = EINVAL;
	goto out;
    }
    if ((value = vnaproperty_scalar_alloc(vs.vs_position)) == NULL) {
	goto out;
    }

    /*
     * Go through the expression list (working from leaf to root).
     */
    for (vnaproperty_expr_t *vexp = head; vexp != NULL; vexp = vexp->vex_next) {
	vnaproperty_t *collection;

	/*
	 * If the current collection exists, use it; otherwise,
	 * allocate a new one.
	 */
	if ((collection = vexp->vex_collection) == NULL) {
	    switch (vexp->vex_token) {
	    case T_DOT:
		if ((collection = vnaproperty_map_alloc()) == NULL) {
		    goto out;
		}
		break;

	    case T_LBRACKET:
		if ((collection = vnaproperty_list_alloc()) == NULL) {
		    goto out;
		}
		break;

	    default:
		abort();
	    }
	}

	/*
	 * Add the value to the collection.
	 */
	switch (vexp->vex_token) {
	case T_DOT:
	    if (vnaproperty_map_set(collection, vexp->u.vex_name,
			value) == -1) {
		goto out;
	    }
	    break;

	case T_LBRACKET:
	    if (vnaproperty_list_set(collection, vexp->u.vex_index,
			value) == -1) {
		goto out;
	    }
	    break;

	default:
	    abort();
	}

	/*
	 * If the current collection already existed, then we're done.
	 * Otherwise, set value to our newly allocated collection and
	 * keep going.
	 */
	if (vexp->vex_collection != NULL) {
	    assert(collection == vexp->vex_collection);
	    value = NULL;
	    result = 0;
	    goto out;
	}
	assert(collection != NULL);
	value = collection;
    }

    /*
     * Replace *anchor with value.
     */
    assert(value != NULL);
    vnaproperty_free(*anchor);
    *anchor = value;
    value = NULL;
    result = 0;

out:
    vnaproperty_free(value);
    _vnaproperty_expr_free(head);
    free((void *)vs.vs_input);
    return result;
}

/*
 * vnaproperty_expr_vdelete: delete the value described by format
 *   @anchor: address of root property pointer
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 */
int vnaproperty_expr_vdelete(vnaproperty_t **anchor, const char *format,
	va_list ap)
{
    vnaproperty_expr_t *head = NULL;
    vnaproperty_scanner_t vs;
    vnaproperty_t *property;
    int rv, result = -1;

    /*
     * Parse the expression.
     */
    if ((rv = _vnaproperty_parse(&vs, &head, *anchor, format, ap)) == -1)
	goto out;

    /*
     * Make sure the current token is T_EOF.
     */
    if (vs.vs_token != T_EOF) {
	errno = EINVAL;
	goto out;
    }

    /*
     * If the expression list is '.' or empty, then delete the root.
     */
    if (head == NULL) {
	vnaproperty_free(*anchor);
	*anchor = NULL;
	result = 0;
	goto out;
    }

    /*
     * If vex_collection is NULL, then the path leading down either
     * doesn't all exist, or something along the way is the wrong
     * type.  Errno is already set.
     */
    if ((property = head->vex_collection) == NULL) {
	goto out;
    }

    /*
     * Delete the indicated key.
     */
    switch (head->vex_token) {
    case T_DOT:
	result = vnaproperty_map_delete(property, head->u.vex_name);
	goto out;

    case T_LBRACKET:
	result = vnaproperty_list_delete(property, head->u.vex_index);
	goto out;

    default:
	abort();
    }

out:
    _vnaproperty_expr_free(head);
    free((void *)vs.vs_input);
    return result;
}

/*
 * vnaproperty_expr_type: get the type of the given property expression
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the expression
 *   @...:    optional variable arguments
 */
vnaproperty_type_t vnaproperty_expr_type(const vnaproperty_t *root,
	const char *format, ...)
{
    va_list ap;
    vnaproperty_type_t type;

    va_start(ap, format);
    type = vnaproperty_expr_vtype(root, format, ap);
    va_end(ap);

    return type;
}

/*
 * vnaproperty_expr_count: return count of objects in given collection
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 */
int vnaproperty_expr_count(const vnaproperty_t *root,
	const char *format, ...)
{
    va_list ap;
    int count;

    va_start(ap, format);
    count = vnaproperty_expr_vcount(root, format, ap);
    va_end(ap);

    return count;
}

/*
 * vnaproperty_expr_keys: return a vector of keys for the given map expr
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 *
 * Caller can free the vector by a call to free.
 */
const char **vnaproperty_expr_keys(const vnaproperty_t *root,
	const char *format, ...)
{
    va_list ap;
    const char **keys;

    va_start(ap, format);
    keys = vnaproperty_expr_vkeys(root, format, ap);
    va_end(ap);

    return keys;
}

/*
 * vnaproperty_expr_get: get a property value from a property expression
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 */
const char *vnaproperty_expr_get(const vnaproperty_t *root,
	const char *format, ...)
{
    va_list ap;
    const char *value;

    va_start(ap, format);
    value = vnaproperty_expr_vget(root, format, ap);
    va_end(ap);

    return value;
}

/*
 * vnaproperty_expr_set: set a property value from a property expression
 *   @anchor: address of root property pointer
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 */
int vnaproperty_expr_set(vnaproperty_t **anchor, const char *format, ...)
{
    va_list ap;
    int rv;

    va_start(ap, format);
    rv = vnaproperty_expr_vset(anchor, format, ap);
    va_end(ap);

    return rv;
}

/*
 * vnaproperty_expr_delete: delete the value described by format
 *   @anchor: address of root property pointer
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 */
int vnaproperty_expr_delete(vnaproperty_t **anchor, const char *format, ...)
{
    va_list ap;
    int rv;

    va_start(ap, format);
    rv = vnaproperty_expr_vdelete(anchor, format, ap);
    va_end(ap);

    return rv;
}
