/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2022 D Scott Guthridge <scott_guthridge@rompromity.net>
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
#include <yaml.h>
#include "vnaerr_internal.h"
#include "vnaproperty_internal.h"

static void vnaproperty_free(vnaproperty_t *root);	/* forward */


/***********************************************************************
 * Scalars
 **********************************************************************/

/*
 * scalar_alloc: allocate a new scalar element
 *   @value: value of the scalar (string)
 */
static vnaproperty_t *scalar_alloc(const char *value)
{
    char *copy;
    vnaproperty_scalar_t *vpsp;

    if ((copy = strdup(value)) == NULL) {
	return NULL;
    }
    if ((vpsp = malloc(sizeof(vnaproperty_scalar_t))) == NULL) {
	free((void *)copy);
	return NULL;
    }
    (void)memset((void *)vpsp, 0, sizeof(*vpsp));
    vpsp->vps_base.vpr_type = VNAPROPERTY_SCALAR;
    vpsp->vps_value = copy;

    return ((vnaproperty_t *)vpsp);
}

/*
 * scalar_get: return the value of a scalar
 *   @scalar: pointer to scalar element
 */
static const char *scalar_get(const vnaproperty_t *scalar)
{
    vnaproperty_scalar_t *vpsp;

    if (scalar->vpr_type != VNAPROPERTY_SCALAR) {
	errno = EINVAL;
	return NULL;
    }
    vpsp = (vnaproperty_scalar_t *)scalar;
    return vpsp->vps_value;
}


/***********************************************************************
 * Maps
 **********************************************************************/

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
 * map_compare_keys: compare key and hashval to the given element
 *   @key: key string
 *   @hashval: full 32-bit hashval
 *   @vmep: element with which to compare
 *
 * compares by hashval first, then by string compare
 */
static int map_compare_keys(const char *key, uint32_t hashval,
	vnaproperty_map_element_t *vmep)
{
    if (hashval < vmep->vme_hashval)
	return -1;

    if (hashval > vmep->vme_hashval)
	return +1;

    return strcmp(key, vmep->vme_pair.vmpr_key);
}

/*
 * map_find_anchor: find insertion point
 *   @vplp:    pointer to property map structure
 *   @anchor:  address of address of the insertion point
 *   @key:     search key
 *   @hashval: full 32 bit hash value
 *
 * If the key is found return true, else return false.  In either case,
 * *anchor is set to the insertion point, i.e. if the element is found,
 * **anchor is the found element; if not found, *anchor is the address
 * of the pointer where the element should be inserted.
 */
static bool map_find_anchor(vnaproperty_map_t *vpmp,
	vnaproperty_map_element_t ***anchor,
	const char *key, uint32_t hashval)
{
    vnaproperty_map_element_t *vmep, **vmepp;
    int bucket;
    int cmp = 1;

    assert(vpmp->vpm_hash_size != 0);
    bucket = hashval % vpmp->vpm_hash_size;
    vmepp = &vpmp->vpm_hash_table[bucket];
    for (; (vmep = *vmepp) != NULL; vmepp = &vmep->vme_hash_next) {
	if ((cmp = map_compare_keys(key, hashval, vmep)) <= 0)
	    break;
    }
    *anchor = vmepp;
    return cmp == 0;
}

/*
 * map_expand: increase the size of the hash table and rehash all elements
 *   @vpmp: map structure
 */
static int map_expand(vnaproperty_map_t *vpmp)
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
	     * correct location, it has to be reprocessed at most once.
	     */
	    anchor = &new_table[cur->vme_hashval % new_allocation];
	    for (; (next = *anchor) != NULL; anchor = &next->vme_hash_next) {
		if (map_compare_keys(cur->vme_pair.vmpr_key,
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
 * map_append_order_element: append to order list
 *   @vplp: pointer to property map structure
 *   @vmep: element to append
 */
static void map_append_order_element(vnaproperty_map_t *vpmp,
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
 * map_delete_order_element: delete from order list
 *   @vplp: pointer to property map structure
 *   @vmep: element to delete
 */
static void map_delete_order_element(vnaproperty_map_t *vpmp,
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
 * map_alloc: allocate a new map element
 */
static vnaproperty_t *map_alloc()
{
    vnaproperty_map_t *vpmp;

    if ((vpmp = malloc(sizeof(vnaproperty_map_t))) == NULL) {
	return NULL;
    }
    (void)memset((void *)vpmp, 0, sizeof(*vpmp));
    vpmp->vpm_base.vpr_type = VNAPROPERTY_MAP;
    return ((vnaproperty_t *)vpmp);
}

/*
 * map_count: return the number of items in the map
 *   @list: pointer to map element
 */
static int map_count(const vnaproperty_t *map)
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
 * map_subtree: return the address of the subtree of key
 *   @list: pointer to map element
 *   @add:  insert key if not found
 *   @key:  search key
 */
static vnaproperty_t **map_subtree(vnaproperty_t *map, bool add,
	const char *key)
{
    vnaproperty_map_t *vpmp;
    vnaproperty_map_element_t **anchor, *vmep;
    uint32_t hashval;

    if (map->vpr_type != VNAPROPERTY_MAP) {
	errno = EINVAL;
	return NULL;
    }
    vpmp = (vnaproperty_map_t *)map;
    if (vpmp->vpm_count + 1 >= 2 * vpmp->vpm_hash_size) {
	if (map_expand(vpmp) == -1) {
	    return NULL;
	}
    }
    hashval = crc32c(-1, (void *)key, strlen(key));
    if (map_find_anchor(vpmp, &anchor, key, hashval)) {
	vmep = *anchor;
	return &vmep->vme_pair.vmpr_value;
    }
    if ((vmep = malloc(sizeof(vnaproperty_map_element_t))) == NULL) {
	return NULL;
    }
    if (!add) {
	errno = ENOENT;
	return NULL;
    }
    (void)memset((void *)vmep, 0, sizeof(*vmep));
    if ((vmep->vme_pair.vmpr_key = strdup(key)) == NULL) {
	free((void *)vmep);
	return NULL;
    }
    vmep->vme_magic = VNAPROPERTY_MAP_PAIR_ELEMENT_MAGIC;
    vmep->vme_hashval = hashval;
    vmep->vme_hash_next = *anchor;
    *anchor = vmep;
    map_append_order_element(vpmp, vmep);
    ++vpmp->vpm_count;
    return &vmep->vme_pair.vmpr_value;
}

/*
 * map_delete: delete the element with given key
 *   @list: pointer to map element
 *   @key: search key
 */
static int map_delete(vnaproperty_t *map, const char *key)
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
    if (!map_find_anchor(vpmp, &anchor, key, hashval)) {
	errno = ENOENT;
	return -1;
    }
    vmep = *anchor;
    assert(vpmp->vpm_count > 0);
    map_delete_order_element(vpmp, vmep);
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


/***********************************************************************
 * Lists
 **********************************************************************/

/*
 * list_check_allocation: extend allocation to >= size
 *   @vplp: pointer to property list structure
 *   @size: needed size
 */
static int list_check_allocation(vnaproperty_list_t *vplp, size_t size)
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
 * list_alloc: allocate a new list element
 */
static vnaproperty_t *list_alloc()
{
    vnaproperty_list_t *vplp;

    if ((vplp = malloc(sizeof(vnaproperty_list_t))) == NULL) {
	return NULL;
    }
    (void)memset((void *)vplp, 0, sizeof(*vplp));
    vplp->vpl_base.vpr_type = VNAPROPERTY_LIST;

    return ((vnaproperty_t *)vplp);
}

/*
 * list_count: return the number of items in the list
 *   @list: pointer to list element
 */
static int list_count(const vnaproperty_t *list)
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
 * list_subtree: return address of subtree at given index
 *   @list:  pointer to list element
 *   @add:   add elements to list if needed
 *   @index: index into list
 */
static vnaproperty_t **list_subtree(vnaproperty_t *list,
	bool add, int index)
{
    vnaproperty_list_t *vplp;

    if (list->vpr_type != VNAPROPERTY_LIST || index < 0) {
	errno = EINVAL;
	return NULL;
    }
    vplp = (vnaproperty_list_t *)list;
    if (index >= vplp->vpl_length) {	/* extend case */
	if (!add) {
	    errno = ENOENT;
	    return NULL;
	}
	if (list_check_allocation(vplp, index + 1) == -1) {
	    return NULL;
	}
	vplp->vpl_length = index + 1;
    }
    return &vplp->vpl_vector[index];
}

/*
 * list_insert: insert null element at index and return address
 *   @list: pointer to list element
 *   @index: index where element should be inserted
 *
 * Index must be in 0..N where N is the number of elements in the list.
 */
static vnaproperty_t **list_insert(vnaproperty_t *list, int index)
{
    vnaproperty_list_t *vplp;

    if (list->vpr_type != VNAPROPERTY_LIST || index < 0) {
	errno = EINVAL;
	return NULL;
    }
    vplp = (vnaproperty_list_t *)list;
    if (index >= vplp->vpl_length) {
	return list_subtree(list, /*add*/true, index);
    }
    if (list_check_allocation(vplp, vplp->vpl_length + 1) == -1) {
	return NULL;
    }
    (void)memmove((void *)&vplp->vpl_vector[index + 1],
		  (void *)&vplp->vpl_vector[index],
		  (vplp->vpl_length - index) * sizeof(vnaproperty_t *));
    vplp->vpl_vector[index] = NULL;
    ++vplp->vpl_length;

    return &vplp->vpl_vector[index];
}

/*
 * list_append: append null element and return address
 *   @list: pointer to list element
 */
static vnaproperty_t **list_append(vnaproperty_t *list)
{
    vnaproperty_list_t *vplp;
    vnaproperty_t **result;

    if (list->vpr_type != VNAPROPERTY_LIST) {
	errno = EINVAL;
	return NULL;
    }
    vplp = (vnaproperty_list_t *)list;
    if (list_check_allocation(vplp, vplp->vpl_length + 1) == -1) {
	return NULL;
    }
    *(result = &vplp->vpl_vector[vplp->vpl_length++]) = NULL;

    return result;
}

/*
 * list_delete: delete the element at the given index
 *   @list: pointer to list element
 *   @index: index of element to delete
 */
static int list_delete(vnaproperty_t *list, int index)
{
    vnaproperty_list_t *vplp;

    if (list->vpr_type != VNAPROPERTY_LIST || index < 0) {
	errno = EINVAL;
	return -1;
    }
    vplp = (vnaproperty_list_t *)list;
    if (index >= vplp->vpl_length) {
	errno = ENOENT;
	return -1;
    }
    (void)memmove((void *)&vplp->vpl_vector[index],
		  (void *)&vplp->vpl_vector[index + 1],
		  (vplp->vpl_length - index) * sizeof(vnaproperty_t *));
    vplp->vpl_vector[--vplp->vpl_length] = NULL;
    return 0;
}


/***********************************************************************
 * Property Scanner & Parser
 **********************************************************************/

/*
 * vnaproperty_token_t: tokens returned by the scanner
 */
typedef enum vnaproperty_token {
    T_ERROR	= -1,
    T_EOF,
    T_HASH,
    T_PLUS,
    T_DOT,
    T_ASSIGN,
    T_LBRACKET,
    T_RBRACKET,
    T_LCURLY,
    T_RCURLY,
    T_ID,
    T_INT
} vnaproperty_token_t;

/*
 * scanner_t: scanner state
 */
typedef struct scanner {
    char               *scn_input;
    char               *scn_position;
    char                scn_cur;
    char               *scn_text;
    union {
	int		scn_int;
    } u;
    vnaproperty_token_t scn_token;
} scanner_t;

/*
 * VNAPROPERTY_GETCHAR: advance to the next input character
 */
#define VNAPROPERTY_GETCHAR(scnp) \
	((scnp)->scn_cur = *++(scnp)->scn_position)

/*
 * ISIDCHAR1: return true if c can start an identifier
 */
#define ISIDCHAR1(c) \
	(isalpha(c) || !isascii(c) || (c) == '_' || (c) == '\\')

/*
 * ISIDCHAR: return true if c can continue an identifier
 */
#define ISIDCHAR(c) \
	(isalpha(c) || isdigit(c) || !isascii(c) || (c) == ' ' || \
	 (c) == '_' || (c) == '-' || (c) == '\\')

/*
 * scan: return the next input token
 *   @scnp: scanner state
 */
static void scan(scanner_t *scnp)
{
    for (;;) {
	switch (scnp->scn_cur) {
	case '\000':
	    scnp->scn_token = T_EOF;
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
	    VNAPROPERTY_GETCHAR(scnp);
	    continue;

	/*
	 * Handle simple tokens.
	 */
	case '#':
	    VNAPROPERTY_GETCHAR(scnp);
	    scnp->scn_token = T_HASH;
	    return;

	case '+':
	    VNAPROPERTY_GETCHAR(scnp);
	    scnp->scn_token = T_PLUS;
	    return;

	case '.':
	    VNAPROPERTY_GETCHAR(scnp);
	    scnp->scn_token = T_DOT;
	    return;

	case '=':
	    VNAPROPERTY_GETCHAR(scnp);
	    scnp->scn_token = T_ASSIGN;
	    return;

	case '[':
	    VNAPROPERTY_GETCHAR(scnp);
	    scnp->scn_token = T_LBRACKET;
	    return;

	case ']':
	    VNAPROPERTY_GETCHAR(scnp);
	    scnp->scn_token = T_RBRACKET;
	    return;

	case '{':
	    VNAPROPERTY_GETCHAR(scnp);
	    scnp->scn_token = T_LCURLY;
	    return;

	case '}':
	    VNAPROPERTY_GETCHAR(scnp);
	    scnp->scn_token = T_RCURLY;
	    return;

	default:
	    break;
	}

	/*
	 * Integers
	 */
	if (isdigit(scnp->scn_cur)) {
	    scnp->scn_text = scnp->scn_position;
	    do {
		VNAPROPERTY_GETCHAR(scnp);
	    } while (isdigit(scnp->scn_cur));
	    *scnp->scn_position = '\000';
	    scnp->u.scn_int = strtol(scnp->scn_text, NULL, 10);
	    scnp->scn_token = T_INT;
	    return;
	}

	/*
	 * Identifiers
	 */
	if (ISIDCHAR1(scnp->scn_cur)) {
	    char *protected = scnp->scn_position;
	    char *destination = scnp->scn_position;

	    scnp->scn_text = scnp->scn_position;
	    do {
		if (scnp->scn_cur == '\\') {
		    VNAPROPERTY_GETCHAR(scnp);
		    if (scnp->scn_cur == '\000') {
			scnp->scn_token = T_ERROR;
			return;
		    }
		    protected = scnp->scn_position;
		}
		*destination++ = scnp->scn_cur;
		VNAPROPERTY_GETCHAR(scnp);
	    } while (ISIDCHAR(scnp->scn_cur));

	    /*
	     * Trim trailing spaces.
	     */
	    while (&destination[-1] > protected &&
		    destination[-1] == ' ') {
		--destination;
	    }
	    *destination = '\000';
	    scnp->scn_token = T_ID;
	    return;
	}

	/*
	 * All other characters are reserved.
	 */
	scnp->scn_token = T_ERROR;
	return;
    }
}

/*
 * expr_type_t: expression node type
 */
typedef enum expr_type {
    E_MAP,				/* {}   */
    E_MAP_ELEMENT,			/* foo  */
    E_LIST,				/* []   */
    E_LIST_ELEMENT,			/* [5]  */
    E_LIST_INSERT,			/* [5+] */
    E_LIST_APPEND,			/* [+] */
    E_DOT				/* .    */
} expr_type_t;

/*
 * expr_t: expression node
 */
typedef struct expr {
    expr_type_t ex_type;		/* expression node type */
    struct expr *ex_next;		/* next expression node */
    union {
	char   *ex_name;		/* key for map element */
	int     ex_index;		/* index for list element */
    } u;
} expr_t;

/*
 * parser_t: parser state
 */
typedef struct parser {
    scanner_t		prs_scn;		/* scanner state */
    expr_t	       *prs_head;	/* expression list head */
    expr_t	       *prs_tail;	/* expression list tail */
    vnaproperty_t      *prs_collection;	/* if last elem is map/list */
} parser_t;

/*
 * expr_free: free scanner and parser resources
 *   @parser: pointer to parser state structure
 */
static void parser_free(parser_t *parser)
{
    while (parser->prs_head != NULL) {
	expr_t *exp = parser->prs_head;

	parser->prs_head = exp->ex_next;
	free((void *)exp);
    }
    free((void *)parser->prs_scn.scn_input);
    parser->prs_scn.scn_input = NULL;
}

/*
 * parse: parse a property expression and return expression list
 *   @parser: address of caller-allocated parser state structure
 *   @result: address of vnaproperty_t pointer to receive result
 *   @format: sprintf format
 *   @ap:     variable argument pointer
 *
 *   On success, *exp_anchor is a linked list of expression nodes
 *   representing the given expression in left-to-right order.
 */
static int parse(parser_t *parser,
	const char *format, va_list ap)
{
    scanner_t *scanner = &parser->prs_scn;
    expr_t *exp, **expp;
    int state;

    /*
     * Init the parser and format the user's arguments.
     */
    (void)memset((void *)parser, 0, sizeof(*parser));
    if (vasprintf(&scanner->scn_input, format, ap) == -1) {
	return -1;
    }
    scanner->scn_position = scanner->scn_input;
    scanner->scn_cur = scanner->scn_input[0];
    expp = &parser->prs_head;
    exp = NULL;
    scan(scanner);

    /*
     * LL(1) Grammar
     *
     * expr is the starting state
     * [token] means the token is left on the input
     *
     * expr		: T_DOT dot			// .
     *			| [T_ID] map_element		// abc
     *			| T_LBRACKET subscript		// [5], [5+], [+], []
     *			| [T_LCURLY] abstract_map	// {}
     *			;
     *
     * dot		: [T_ID] map_element		// .abc
     *			| T_LBRACKET subscript		// .[5], .[5+], etc
     *			| [T_LCURLY] abstract_map	// .{}
     *			| final_dot			// accepting
     *			;
     *
     * chain		: T_DOT dot			// foo.
     *			| T_LBRACKET subscript		// foo.bar[5]
     *			| [T_LCURLY] abstract_map	// foo.bar{}
     *			| λ				// accepting
     *			;
     *
     * subscript	: T_INT list_element		// 5], 5+]
     *			| T_PLUS list_append		// +]
     *			| [T_RBRACKET] abstract_list	// ]
     *			;
     *
     * map_element	: T_ID chain			// foo
     *			;
     *
     * list_element	: T_PLUS list_insert		// 5+]
     *			| T_RBRACKET chain		// 5]
     *			;
     *
     * list_insert	: T_RBRACKET chain		// 5+]
     *			;
     *
     * list_append	: T_RBRACKET chain		// +]
     *			;
     *
     * final_dot	: λ
     *			;
     *
     * abstract_map	: T_LCURLY T_RCURLY
     *			;
     *
     * abstract_list	: T_RBRACKET
     *			;
     *
     * The language is regular, so we parse it simply using Duff's device.
     */
    state = 0;
    for (;;) {
	switch (state) {
	case 0:
	    /*
	     * expr	: T_DOT dot
	     *		| [T_ID] map_element
	     *		| T_LBRACKET subscript
	     *		| [T_LCURLY] abstract_map
	     *		;
	     */
	    switch (scanner->scn_token) {
	    case T_DOT:
		scan(scanner);
		state = 1;
		continue;

	    case T_ID:
		goto map_element;

	    case T_LBRACKET:
		scan(scanner);
		state = 3;
		continue;

	    case T_LCURLY:
		goto abstract_map;

	    default:
		goto error;				/* NOT accepting */
	    }
	    break;

	case 1:
	    /*
	     * dot	: [T_ID] map_element
	     *		| T_LBRACKET subscript
	     *		| [T_LCURLY] abstract_map
	     *		| [default]  final_dot
	     *		;
	     */
	    switch (scanner->scn_token) {
	    case T_ID:
		goto map_element;

	    case T_LBRACKET:
		scan(scanner);
		state = 3;
		continue;

	    case T_LCURLY:
		goto abstract_map;

	    default:
		goto final_dot;				/* accepting state */
	    }
	    break;

	case 2:
	    /*
	     * chain	: T_DOT dot
	     *		| T_LBRACKET subscript
	     *		| [T_LCURLY] abstract_map
	     *		| λ
	     *		;
	     */
	    switch (scanner->scn_token) {
	    case T_DOT:
		scan(scanner);
		state = 1;
		continue;

	    case T_LBRACKET:
		scan(scanner);
		state = 3;
		continue;

	    case T_LCURLY:
		goto abstract_map;

	    default:
		break;					/* accepting state */
	    }
	    break;

	case 3:
	    /*
	     * subscript	: T_INT list_element
	     *			| T_PLUS list_append
	     *			| [T_RBRACKET] abstract_list
	     *			;
	     */
	    switch (scanner->scn_token) {
	    case T_INT:
		goto list_element;

	    case T_PLUS:
		goto list_append;

	    case T_RBRACKET:
		goto abstract_list;

	    default:
		goto error;				/* NOT accepting */
	    }
	    break;

	map_element:
	    /*
	     * map_element	: T_ID chain ;
	     */
	    assert(scanner->scn_token == T_ID);
	    if ((exp = malloc(sizeof(expr_t))) == NULL) {
		goto error;
	    }
	    (void)memset((void *)exp, 0, sizeof(expr_t));
	    exp->ex_type = E_MAP_ELEMENT;
	    exp->u.ex_name = scanner->scn_text;
	    exp->ex_next = NULL;
	    *expp = exp;
	    expp = &exp->ex_next;
	    scan(scanner);
	    state = 2;
	    continue;

	list_element:	/* and list_insert */
	    /*
	     * list_element	: T_PLUS list_insert
	     *			| T_RBRACKET chain
	     *			;
	     *
	     * list_insert	: T_RBRACKET chain
	     *			;
	     */
	    assert(scanner->scn_token == T_INT);
	    if ((exp = malloc(sizeof(expr_t))) == NULL) {
		goto error;
	    }
	    (void)memset((void *)exp, 0, sizeof(expr_t));
	    exp->ex_type = E_LIST_ELEMENT;
	    exp->u.ex_index = scanner->u.scn_int;
	    scan(scanner);
	    if (scanner->scn_token == T_PLUS) {
		exp->ex_type = E_LIST_INSERT;
		scan(scanner);
	    }
	    *expp = exp;
	    expp = &exp->ex_next;
	    if (scanner->scn_token != T_RBRACKET) {
		errno = EINVAL;
		goto error;
	    }
	    scan(scanner);
	    state = 2;
	    continue;

	list_append:
	    /*
	     * list_append	: T_RBRACKET chain
	     */
	    assert(scanner->scn_token == T_PLUS);
	    if ((exp = malloc(sizeof(expr_t))) == NULL) {
		goto error;
	    }
	    (void)memset((void *)exp, 0, sizeof(expr_t));
	    exp->ex_type = E_LIST_APPEND;
	    *expp = exp;
	    expp = &exp->ex_next;
	    scan(scanner);
	    if (scanner->scn_token != T_RBRACKET) {
		errno = EINVAL;
		goto error;
	    }
	    scan(scanner);
	    state = 2;
	    continue;

	final_dot:
	    /*
	     * final_dot	: λ ;
	     */
	    if ((exp = malloc(sizeof(expr_t))) == NULL) {
		goto error;
	    }
	    (void)memset((void *)exp, 0, sizeof(expr_t));
	    exp->ex_type = E_DOT;
	    exp->ex_next = NULL;
	    *expp = exp;
	    expp = &exp->ex_next;
	    break;

	abstract_map:
	    /*
	     * abstract_map	: T_LCURLY T_RCURLY ;
	     */
	    assert(scanner->scn_token == T_LCURLY);
	    scan(scanner);
	    if (scanner->scn_token != T_RCURLY) {
		goto error;
	    }
	    scan(scanner);
	    if ((exp = malloc(sizeof(expr_t))) == NULL) {
		goto error;
	    }
	    (void)memset((void *)exp, 0, sizeof(expr_t));
	    exp->ex_type = E_MAP;
	    exp->ex_next = NULL;
	    *expp = exp;
	    expp = &exp->ex_next;
	    break;

	abstract_list:
	    /*
	     * abstract_list	: T_RBRACKET ;
	     */
	    assert(scanner->scn_token == T_RBRACKET);
	    scan(scanner);
	    if ((exp = malloc(sizeof(expr_t))) == NULL) {
		goto error;
	    }
	    (void)memset((void *)exp, 0, sizeof(expr_t));
	    exp->ex_type = E_LIST;
	    exp->ex_next = NULL;
	    *expp = exp;
	    expp = &exp->ex_next;
	    break;
	}
	break;
    }
    parser->prs_tail = exp;
    return 0;

error:
    parser_free(parser);
    errno = EINVAL;
    return -1;
}


/***********************************************************************
 * Internal Tree Operations
 **********************************************************************/

/*
 * parse_and_descend: parse the expression and descend down the tree
 *   @parser:  address of caller-allocated parser state structure
 *   @rootptr: address of property data root
 *   @set:     force the tree to conform to the indicated expression
 *   @format:  printf-like format string forming the property expression
 *   @ap       variable argument pointer
 */
static vnaproperty_t **parse_and_descend(parser_t *parser,
	vnaproperty_t **rootptr, bool set, const char *format, va_list ap)
{
    vnaproperty_t **anchor = rootptr;
    vnaproperty_t *node = *anchor;
    vnaproperty_t *collection = NULL;

    /*
     * Parse the expression.
     */
    if (parse(parser, format, ap) == -1) {
	return NULL;
    }

    /*
     * Following the expression list, walk down the tree.
     */
    for (expr_t *exp = parser->prs_head; exp != NULL; exp = exp->ex_next) {
	collection = NULL;
	if (exp->ex_type == E_DOT) {
	    assert(exp->ex_next == NULL);
	    break;
	}
	switch (exp->ex_type) {
	case E_MAP:
	case E_MAP_ELEMENT:
	    if (node == NULL) {
		if (!set) {
		    errno = ENOENT;
		    goto error;
		}
		if ((node = map_alloc()) == NULL) {
		    goto error;
		}
		*anchor = node;

	    } else if (node->vpr_type != VNAPROPERTY_MAP) {
		if (!set) {
		    node = NULL;
		    errno = EINVAL;
		    goto error;
		}
		vnaproperty_free(node);
		*anchor = node = NULL;
		if ((node = map_alloc()) == NULL) {
		    goto error;
		}
		*anchor = node;
	    }
	    if (exp->ex_type == E_MAP) {
		assert(exp->ex_next == NULL);
		break;
	    }
	    collection = node;
	    if ((anchor = map_subtree(node, set,
			    exp->u.ex_name)) == NULL) {
		goto error;
	    }
	    node = *anchor;
	    continue;

	case E_LIST:
	case E_LIST_ELEMENT:
	case E_LIST_INSERT:
	case E_LIST_APPEND:
	    if (node == NULL) {
		if (!set) {
		    errno = ENOENT;
		    goto error;
		}
		if ((node = list_alloc()) == NULL) {
		    goto error;
		}
		*anchor = node;

	    } else if (node->vpr_type != VNAPROPERTY_LIST) {
		if (!set) {
		    node = NULL;
		    errno = EINVAL;
		    goto error;
		}
		vnaproperty_free(node);
		*anchor = node = NULL;
		if ((node = list_alloc()) == NULL) {
		    goto error;
		}
		*anchor = node;
	    }
	    switch (exp->ex_type) {
	    case E_LIST:
		assert(exp->ex_next == NULL);
		break;

	    case E_LIST_ELEMENT:
		collection = node;
		if ((anchor = list_subtree(node, set,
				exp->u.ex_index)) == NULL) {
		    goto error;
		}
		node = *anchor;
		continue;

	    case E_LIST_INSERT:
		if (!set) {
		    errno = EINVAL;
		    goto error;
		}
		collection = node;
		if ((anchor = list_insert(node, exp->u.ex_index)) == NULL) {
		    goto error;
		}
		node = *anchor;
		continue;

	    case E_LIST_APPEND:
		if (!set) {
		    errno = EINVAL;
		    goto error;
		}
		collection = node;
		if ((anchor = list_append(node)) == NULL) {
		    goto error;
		}
		node = *anchor;
		continue;

	    default:
		abort();
	    }
	    break;

	    if (exp->ex_type == E_LIST) {
		assert(exp->ex_next == NULL);
		break;
	    }
	    collection = node;
	    if ((anchor = list_subtree(node, set,
			    exp->u.ex_index)) == NULL) {
		goto error;
	    }
	    node = *anchor;
	    continue;

	default:
	    abort();
	}
    }
    parser->prs_collection = collection;
    return anchor;

error:
    parser_free(parser);
    return NULL;
}

/*
 * get_node: parse the expression and return the indicated node
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @ap:     variable argument pointer
 */
static const vnaproperty_t *get_node(const vnaproperty_t *root,
	const char *format, va_list ap)
{
    parser_t parser;
    scanner_t *scanner = &parser.prs_scn;
    vnaproperty_t **anchor;

    /*
     * Parse the expression and descend to the requested node.
     */
    if ((anchor = parse_and_descend(&parser, (vnaproperty_t **)&root,
		    /*set*/false, format, ap)) == NULL) {
	return NULL;
    }

    /*
     * Make sure there are no unexpected tokens at the end of line.
     */
    if (scanner->scn_token != T_EOF) {
	errno = EINVAL;
	parser_free(&parser);
	return NULL;
    }

    parser_free(&parser);
    return *anchor;
}

/*
 * vnaproperty_free: free the given tree
 *   @root: pointer to tree to free
 */
static void vnaproperty_free(vnaproperty_t *root)
{
    if (root == NULL)
	return;

    assert(root->vpr_type == VNAPROPERTY_SCALAR ||
	   root->vpr_type == VNAPROPERTY_LIST   ||
	   root->vpr_type == VNAPROPERTY_MAP);
    switch (root->vpr_type) {
    case VNAPROPERTY_SCALAR:
	{
	    vnaproperty_scalar_t *vpsp = (vnaproperty_scalar_t *)root;
	    free((void *)vpsp->vps_value);
	    (void)memset((void *)vpsp, 'X', sizeof(*vpsp));
	}
	break;

    case VNAPROPERTY_LIST:
	{
	    vnaproperty_list_t *vplp  = (vnaproperty_list_t *)root;
	    for (size_t s = 0; s < vplp->vpl_length; ++s) {
		vnaproperty_free((vnaproperty_t *)vplp->vpl_vector[s]);
	    }
	    free((void *)vplp->vpl_vector);
	    (void)memset((void *)vplp, 'X', sizeof(*vplp));
	}
	break;

    case VNAPROPERTY_MAP:
	{
	    vnaproperty_map_t *vpmp = (vnaproperty_map_t *)root;
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
    free((void *)root);
}


/***********************************************************************
 * External API
 **********************************************************************/

/*
 * vnaproperty_vtype: get the type of the given property expression
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @ap:     variable argument pointer
 *
 * Return:
 *   'm' map
 *   'l' list
 *   's' scalar
 *   -1  error
 */
int vnaproperty_vtype(const vnaproperty_t *root, const char *format, va_list ap)
{
    const vnaproperty_t *node;

    if ((node = get_node(root, format, ap)) == NULL) {
	return VNAPROPERTY_ERROR;
    }

    /*
     * Convert to character.
     */
    switch (node->vpr_type) {
    case VNAPROPERTY_SCALAR:
	return 's';

    case VNAPROPERTY_MAP:
	return 'm';

    case VNAPROPERTY_LIST:
	return 'l';

    default:
	return -1;
    }
}

/*
 * vnaproperty_vcount: return count of elements in given collection
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @ap:     variable argument pointer
 */
int vnaproperty_vcount(const vnaproperty_t *root,
	const char *format, va_list ap)
{
    const vnaproperty_t *node;

    if ((node = get_node(root, format, ap)) == NULL) {
	return -1;
    }
    switch (node->vpr_type) {
    case VNAPROPERTY_MAP:
	return map_count(node);

    case VNAPROPERTY_LIST:
	return list_count(node);

    default:
	break;
    }
    errno = EINVAL;
    return -1;
}

/*
 * vnaproperty_vkeys: return a vector of keys for the given map expr
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @ap:     variable argument pointer
 */
const char **vnaproperty_vkeys(const vnaproperty_t *root,
	const char *format, va_list ap)
{
    const vnaproperty_t *node;
    const vnaproperty_map_t *map;
    const vnaproperty_map_element_t *vmep;
    const char **vector, **cpp;

    if ((node = get_node(root, format, ap)) == NULL) {
	return NULL;
    }
    if (node->vpr_type != VNAPROPERTY_MAP) {
	errno = EINVAL;
	return NULL;
    }
    map = (vnaproperty_map_t *)node;
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
 * vnaproperty_vget: get a property value from a property expression
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @ap:     variable argument pointer
 */
const char *vnaproperty_vget(const vnaproperty_t *root,
	const char *format, va_list ap)
{
    const vnaproperty_t *node;

    if ((node = get_node(root, format, ap)) == NULL) {
	return NULL;
    }
    if (node->vpr_type != VNAPROPERTY_SCALAR) {
	errno = EINVAL;
	return NULL;
    }
    return scalar_get(node);
}

/*
 * vnaproperty_vset: set a property value from a property expression
 *   @rootptr: address of root property pointer
 *   @format:  printf-like format string forming the property expression
 *   @ap:      variable argument pointer
 */
int vnaproperty_vset(vnaproperty_t **rootptr, const char *format, va_list ap)
{
    parser_t parser;
    scanner_t *scanner = &parser.prs_scn;
    vnaproperty_t **anchor;
    vnaproperty_t *value = NULL;
    int rv = -1;

    if ((anchor = parse_and_descend(&parser, rootptr, /*set*/true,
		    format, ap)) == NULL) {
	return -1;
    }

    /*
     * Make sure we're not trying to assign to a map or list.
     */
    switch (parser.prs_tail->ex_type) {
    case E_MAP_ELEMENT:
    case E_LIST_ELEMENT:
    case E_LIST_INSERT:
    case E_LIST_APPEND:
    case E_DOT:
	break;

    case E_MAP:
    case E_LIST:
    default:
	errno = EINVAL;
	goto out;
    }

    /*
     * Get the value to assign.
     *   foo=text	set foo to a scalar with given text
     *   foo#		set foo to null
     */
    switch (scanner->scn_token) {
    case T_ASSIGN:
	value = scalar_alloc(scanner->scn_position);
	if (value == NULL) {
	    goto out;
	}
	break;

    case T_HASH:
	break;

    default:
	errno = EINVAL;
	goto out;
    }

    /*
     * Free any old value and install the new value.
     */
    vnaproperty_free(*anchor);
    *anchor = value;
    rv = 0;

out:
    parser_free(&parser);
    return rv;
}

/*
 * vnaproperty_vdelete: delete the value described by format
 *   @rootptr:   address of root property pointer
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 */
int vnaproperty_vdelete(vnaproperty_t **rootptr, const char *format,
	va_list ap)
{
    parser_t parser;
    scanner_t *scanner = &parser.prs_scn;
    expr_t *tail;
    vnaproperty_t **anchor, *collection;
    int rv = -1;

    /*
     * Parse the expression and descend to the requested node.
     */
    if ((anchor = parse_and_descend(&parser, rootptr, /*set*/false,
		    format, ap)) == NULL) {
	return -1;
    }

    /*
     * Make sure there are no unexpected trailing tokens.
     */
    if (scanner->scn_token != T_EOF) {
	errno = EINVAL;
	goto out;
    }

    /*
     * Delete the indicated key.
     */
    tail = parser.prs_tail;
    collection = parser.prs_collection;
    switch (tail->ex_type) {
    case E_MAP_ELEMENT:
	assert(collection != NULL);
	rv = map_delete(collection, tail->u.ex_name);
	goto out;

    case E_LIST_ELEMENT:
	assert(collection != NULL);
	rv = list_delete(collection, tail->u.ex_index);
	goto out;

    default:
	vnaproperty_free(*anchor);
	*anchor = NULL;
    }
    rv = 0;

out:
    parser_free(&parser);
    return rv;
}

/*
 * vnaproperty_get_subtree: get subtree described by format
 *   @root:   address of root property pointer
 *   @format: printf-like format string forming the property expression
 *   @ap:     variable arguments
 */
vnaproperty_t *vnaproperty_vget_subtree(const vnaproperty_t *root,
	const char *format, va_list ap)
{
    parser_t parser;
    scanner_t *scanner = &parser.prs_scn;
    vnaproperty_t **anchor;

    if ((anchor = parse_and_descend(&parser, (vnaproperty_t **)&root,
		    /*set*/false, format, ap)) == NULL) {
	return NULL;
    }

    /*
     * Make sure there are no unexpected trailing tokens.
     */
    if (scanner->scn_token != T_EOF) {
	errno = EINVAL;
	anchor = NULL;
	goto out;
    }

out:
    parser_free(&parser);
    return *anchor;
}

/*
 * vnaproperty_set_subtree: make the tree conform and return address of subtree
 *   @rootptr: address of root property pointer
 *   @format:  printf-like format string forming the property expression
 *   @ap:      variable arguments
 */
vnaproperty_t **vnaproperty_vset_subtree(vnaproperty_t **rootptr,
	const char *format, va_list ap)
{
    parser_t parser;
    scanner_t *scanner = &parser.prs_scn;
    vnaproperty_t **anchor;

    if ((anchor = parse_and_descend(&parser, rootptr,
		    /*set*/true, format, ap)) == NULL) {
	return NULL;
    }

    /*
     * Make sure there are no unexpected trailing tokens.
     */
    if (scanner->scn_token != T_EOF) {
	errno = EINVAL;
	anchor = NULL;
	goto out;
    }

out:
    parser_free(&parser);
    return anchor;
}

/*
 * vnaproperty_type: get the type of the given property expression
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the expression
 *   @...:    optional variable arguments
 */
int vnaproperty_type(const vnaproperty_t *root, const char *format, ...)
{
    va_list ap;
    int type;

    va_start(ap, format);
    type = vnaproperty_vtype(root, format, ap);
    va_end(ap);

    return type;
}

/*
 * vnaproperty_count: return count of elements in given collection
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 */
int vnaproperty_count(const vnaproperty_t *root,
	const char *format, ...)
{
    va_list ap;
    int count;

    va_start(ap, format);
    count = vnaproperty_vcount(root, format, ap);
    va_end(ap);

    return count;
}

/*
 * vnaproperty_keys: return a vector of keys for the given map expr
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 *
 * Caller can free the vector by a call to free.
 */
const char **vnaproperty_keys(const vnaproperty_t *root,
	const char *format, ...)
{
    va_list ap;
    const char **keys;

    va_start(ap, format);
    keys = vnaproperty_vkeys(root, format, ap);
    va_end(ap);

    return keys;
}

/*
 * vnaproperty_get: get a property value from a property expression
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 */
const char *vnaproperty_get(const vnaproperty_t *root,
	const char *format, ...)
{
    va_list ap;
    const char *value;

    va_start(ap, format);
    value = vnaproperty_vget(root, format, ap);
    va_end(ap);

    return value;
}

/*
 * vnaproperty_set: set a property value from a property expression
 *   @rootptr: address of root property pointer
 *   @format:  printf-like format string forming the property expression
 *   @...:     optional variable arguments
 */
int vnaproperty_set(vnaproperty_t **rootptr, const char *format, ...)
{
    va_list ap;
    int rv;

    va_start(ap, format);
    rv = vnaproperty_vset(rootptr, format, ap);
    va_end(ap);

    return rv;
}

/*
 * vnaproperty_delete: delete the value described by format
 *   @rootptr: address of root property pointer
 *   @format:  printf-like format string forming the property expression
 *   @...:     optional variable arguments
 */
int vnaproperty_delete(vnaproperty_t **rootptr, const char *format, ...)
{
    va_list ap;
    int rv;

    va_start(ap, format);
    rv = vnaproperty_vdelete(rootptr, format, ap);
    va_end(ap);

    return rv;
}

/*
 * vnaproperty_get_subtree: get subtree described by format
 *   @root:    address of root property pointer
 *   @format:  printf-like format string forming the property expression
 *   @...:     optional variable arguments
 */
vnaproperty_t *vnaproperty_get_subtree(const vnaproperty_t *root,
	const char *format, ...)
{
    va_list ap;
    vnaproperty_t *subtree;

    va_start(ap, format);
    subtree = vnaproperty_vget_subtree(root, format, ap);
    va_end(ap);

    return subtree;
}

/*
 * vnaproperty_set_subtree: make the tree conform and return address of subtree
 *   @rootptr: address of root property pointer
 *   @format:  printf-like format string forming the property expression
 *   @...:     optional variable arguments
 */
vnaproperty_t **vnaproperty_set_subtree(vnaproperty_t **rootptr,
	const char *format, ...)
{
    va_list ap;
    vnaproperty_t **subtree;

    va_start(ap, format);
    subtree = vnaproperty_vset_subtree(rootptr, format, ap);
    va_end(ap);

    return subtree;
}

/*
 * vnaproperty_quote_key: quote a map ID that contains reserved chars
 *   @key: map key to quote
 *
 * Caller must free the returned memory by a call to free().
 */
char *vnaproperty_quote_key(const char *key)
{
    int length;
    bool *map = NULL;
    int n_special = 0;
    char *result = NULL;
    char *cur;

    /*
     * Handle NULL.
     */
    if (key == NULL) {
	return NULL;
    }

    /*
     * Handle empty string.
     */
    if (key[0] == '\000') {
	result = strdup("");
	goto out;
    }

    /*
     * Allocate a vector of bool to record positions of special characters.
     */
    length = strlen(key);
    if ((map = calloc(length, sizeof(bool))) == NULL) {
	goto out;
    }

    /*
     * Find the positions that require quotes.
     */
    if (!ISIDCHAR1(key[0]) || key[0] == '\\') {
	map[0] = true;
	++n_special;
    }
    {
	int i;

	for (i = 1; i < length; ++i) {
	    if (!ISIDCHAR(key[i]) || key[i] == '\\') {
		map[i] = true;
		++n_special;
	    }
	}

	/*
	 * Trailing spaces are special.
	 */
	while (i > 1 && key[i - 1] == ' ') {
	    map[--i] = true;
	    ++n_special;
	}
    }

    /*
     * If no quotes required, return a copy of the original string.
     */
    if (n_special == 0) {
	result = strdup(key);
	goto out;
    }

    /*
     * Form the new string.
     */
    if ((result = malloc(length + n_special + 1)) == NULL) {
	goto out;
    }
    cur = result;
    for (int i = 0; i < length; ++i) {
	if (map[i]) {
	    *cur++ = '\\';
	}
	*cur++ = key[i];
    }
    *cur = '\000';

out:
    free((void *)map);
    return result;
}

/*
 * dfs_copy: recursely copy properties
 *   @destination: address of root of the destination tree
 *   @source: root of the source tree
 */
static int dfs_copy(vnaproperty_t **destination, const vnaproperty_t *source)
{
    const char **keys = NULL;
    const char *cp = NULL;
    int count;

    if (source == NULL) {
	return 0;
    }
    switch (source->vpr_type) {
    case VNAPROPERTY_SCALAR:
	if ((cp = scalar_get(source)) != NULL) {
	    if (vnaproperty_set(destination, ".=%s", cp) == -1) {
		return -1;
	    }
	} else {
	    if (vnaproperty_set(destination, ".#") == -1) {
		return -1;
	    }
	}
	break;

    case VNAPROPERTY_MAP:
	if ((keys = vnaproperty_keys(source, ".")) == NULL) {
	    return -1;
	}
	for (const char **cpp = keys; *cpp != NULL; ++cpp) {
	    char *key = NULL;
	    vnaproperty_t **new_destination, *new_source;

	    if ((key = vnaproperty_quote_key(*cpp)) == NULL) {
		free((void *)keys);
		return -1;
	    }
	    new_destination = vnaproperty_set_subtree(destination, ".%s", key);
	    if (new_destination == NULL) {
		free((void *)keys);
		free((void *)key);
		return -1;
	    }
	    new_source = vnaproperty_get_subtree(source, "%s", key);
	    if (dfs_copy(new_destination, new_source) == -1) {
		free((void *)keys);
		free((void *)key);
		return -1;
	    }
	    free((void *)key);
	}
	free((void *)keys);
	break;

    case VNAPROPERTY_LIST:
	count = vnaproperty_count(source, ".");
	for (int i = 0; i < count; ++i) {
	    vnaproperty_t **new_destination, *new_source;

	    new_destination = vnaproperty_set_subtree(destination, "[%d]", i);
	    if (new_destination == NULL) {
		return -1;
	    }
	    new_source = vnaproperty_get_subtree(source, "[%d]", i);
	    if (dfs_copy(new_destination, new_source) == -1) {
		return -1;
	    }
	}
	break;

    default:
	abort();
    }
    return 0;
}

/*
 * vnaproperty_copy: copy a subtree
 *   @destination: address of node where copy is placed
 *   @source: subtree to copy
 */
int vnaproperty_copy(vnaproperty_t **destination, const vnaproperty_t *source)
{
    (void)vnaproperty_delete(destination, ".");
    return dfs_copy(destination, source);
}


/***********************************************************************
 * Undocumented YAML Import / Export
 **********************************************************************/

/*
 * _vnaproperty_yaml_error: report an error
 *   @vymlp:    common argument structure
 *   @category: category of error
 *   @format:   printf format string
 */
void _vnaproperty_yaml_error(const vnaproperty_yaml_t *vymlp,
	vnaerr_category_t category, const char *format, ...)
{
    va_list ap;

    va_start(ap, format);
    _vnaerr_verror(vymlp->vyml_error_fn, vymlp->vyml_error_arg,
	    category, format, ap);
    va_end(ap);
}

/*
 * is_yaml_null_value: test if a string is a yaml null value
 */
static bool is_yaml_null_value(const char *s)
{
    return (s[0] == '~' && s[1] == '\000') ||
	strcmp(s, "null") == 0 ||
	strcmp(s, "Null") == 0 ||
	strcmp(s, "NULL") == 0;
}

/*
 * add_mapping_entry: add a simple scalar tag to value mapping entry
 *   @vymlp:    common argument structure
 *   @t_map:    map into which we're adding
 *   @key:      string valued key
 *   @value:    tag of value to insert
 */
static int add_mapping_entry(vnaproperty_yaml_t *vymlp, int t_map,
	const char *key, int t_value)
{
    yaml_document_t *document = vymlp->vyml_document;
    int t_key;

    errno = 0;
    if ((t_key = yaml_document_add_scalar(document, NULL,
		    (yaml_char_t *)key, strlen(key),
		    YAML_ANY_SCALAR_STYLE)) == 0) {
	if (errno == 0) {
	    errno = EINVAL;
	}
	_vnaproperty_yaml_error(vymlp, VNAERR_SYSTEM,
		"yaml_document_add_scalar: %s: %s",
		vymlp->vyml_filename, strerror(errno));
	return -1;
    }
    if (yaml_document_append_mapping_pair(document, t_map,
		t_key, t_value) == 0) {
	if (errno == 0) {
	    errno = EINVAL;
	}
	_vnaproperty_yaml_error(vymlp, VNAERR_SYSTEM,
		"yaml_document_add_mapping_pair: %s: %s",
		vymlp->vyml_filename, strerror(errno));
	return -1;
    }
    return 0;
}

/*
 * _vnaproperty_yaml_import: import properties from the given YAML document
 *   @vymlp:    common argument structure
 *   @rootptr:  address of property tree root
 *   @vp_node:  yaml node cast to void pointer
 */
int _vnaproperty_yaml_import(vnaproperty_yaml_t *vymlp,
	vnaproperty_t **rootptr, void *vp_node)
{
    yaml_document_t *document = vymlp->vyml_document;
    yaml_node_t *node = vp_node;

    switch (node->type) {
    case YAML_SCALAR_NODE:
	/*
	 * Handle NULL.
	 *
	 * Libyaml 0.2.2 parses untagged null values as !!str.	The scalar
	 * style seems to be the most reliable way to distinguish ~ from
	 * "~", so we use this instead of the tag.
	 */
	if (is_yaml_null_value((const char *)node->data.scalar.value)  &&
		    node->data.scalar.style == YAML_PLAIN_SCALAR_STYLE) {
	    return 0;	/* root is already NULL */
	}

	/*
	 * Handle scalars.
	 */
	if (vnaproperty_set(rootptr, ".=%s",
		    (const char *)node->data.scalar.value) == -1) {
	    _vnaproperty_yaml_error(vymlp, VNAERR_SYSTEM,
		    "_vnaproperty_set: %s: %s",
		    vymlp->vyml_filename, strerror(errno));
	    goto out;
	}
	return 0;

    case YAML_MAPPING_NODE:
	{
	    yaml_node_pair_t *pair;

	    if (vnaproperty_set_subtree(rootptr, "{}") == NULL) {
		_vnaproperty_yaml_error(vymlp, VNAERR_SYSTEM,
			"_vnaproperty_set_subtree: %s: %s",
			vymlp->vyml_filename, strerror(errno));
		goto out;
	    }
	    for (pair = node->data.mapping.pairs.start;
		    pair < node->data.mapping.pairs.top; ++pair) {
		yaml_node_t *key, *value;
		vnaproperty_t **subtree;

		key = yaml_document_get_node(document, pair->key);
		if (key->type != YAML_SCALAR_NODE) {
		    _vnaproperty_yaml_error(vymlp, VNAERR_WARNING,
			    "%s (line %ld) warning: "
			    "non-scalar property key ignored\n",
			    vymlp->vyml_filename, key->start_mark.line + 1);
		    continue;
		}
		value = yaml_document_get_node(document, pair->value);
		if ((subtree = vnaproperty_set_subtree(rootptr, "%s",
			    (const char *)key->data.scalar.value)) == NULL) {
		    _vnaproperty_yaml_error(vymlp, VNAERR_SYSTEM,
			    "_vnaproperty_set_subtree: %s: %s",
			    vymlp->vyml_filename, strerror(errno));
		    goto out;
		}
		if (_vnaproperty_yaml_import(vymlp, subtree, value) == -1) {
		    goto out;
		}
	    }
	    return 0;
	}

    case YAML_SEQUENCE_NODE:
	{
	    yaml_node_item_t *item;

	    if (vnaproperty_set_subtree(rootptr, "[]") == NULL) {
		_vnaproperty_yaml_error(vymlp, VNAERR_SYSTEM,
			"_vnaproperty_set_subtree: %s: %s",
			vymlp->vyml_filename, strerror(errno));
		goto out;
	    }
	    for (item = node->data.sequence.items.start;
		    item < node->data.sequence.items.top; ++item) {
		int index = item - node->data.sequence.items.start;
		yaml_node_t *value;
		vnaproperty_t **subtree;

		value = yaml_document_get_node(document, *item);
		if ((subtree = vnaproperty_set_subtree(rootptr, "[%d]",
				index)) == NULL) {
		    _vnaproperty_yaml_error(vymlp, VNAERR_SYSTEM,
			    "_vnaproperty_set_subtree: %s: %s",
			    vymlp->vyml_filename, strerror(errno));
		    goto out;
		}
		if (_vnaproperty_yaml_import(vymlp, subtree, value) == -1) {
		    goto out;
		}
	    }
	    return 0;
	}

    default:
	abort();
    }

out:
    return -1;
}

/*
 * _vnaproperty_yaml_export: add a property list to the YAML document
 *   @vymlp:    common argument structure
 *   @document: yaml document
 *   @root:     property list root
 */
int _vnaproperty_yaml_export(vnaproperty_yaml_t *vymlp,
	const vnaproperty_t *root)
{
    yaml_document_t *document = vymlp->vyml_document;

    /*
     * Handle NULL
     */
    if (root == NULL) {
	int scalar;
	static const char *value = "~";

	errno = 0;
	if ((scalar = yaml_document_add_scalar(document, NULL,
			(yaml_char_t *)value, strlen(value),
			YAML_PLAIN_SCALAR_STYLE)) == 0) {
	    if (errno == 0) {
		errno = EINVAL;
	    }
	    _vnaproperty_yaml_error(vymlp, VNAERR_SYSTEM,
		    "yaml_document_add_scalar: %s: %s",
		    vymlp->vyml_filename, strerror(errno));
	    return -1;
	}
	return scalar;
    }

    /*
     * Handle scalars, maps and lists.
     */
    switch (root->vpr_type) {
    case VNAPROPERTY_SCALAR:
	{
	    int scalar;
	    const char *value;
	    yaml_scalar_style_t style = YAML_ANY_SCALAR_STYLE;

	    if ((value = vnaproperty_get(root, ".")) == NULL) {
		_vnaproperty_yaml_error(vymlp, VNAERR_INTERNAL,
			"%s: _vnaproperty_get: %s: %s",
			__func__, vymlp->vyml_filename, strerror(errno));
		return -1;
	    }
	    if (strchr(value, '\n') != NULL) {
		style = YAML_LITERAL_SCALAR_STYLE;
	    } else if (is_yaml_null_value(value)) {
		style = YAML_DOUBLE_QUOTED_SCALAR_STYLE;
	    }
	    errno = 0;
	    if ((scalar = yaml_document_add_scalar(document, NULL,
			    (yaml_char_t *)value, strlen(value), style)) == 0) {
		if (errno == 0) {
		    errno = EINVAL;
		}
		_vnaproperty_yaml_error(vymlp, VNAERR_SYSTEM,
			"yaml_document_add_scalar: %s: %s",
			vymlp->vyml_filename, strerror(errno));
		return -1;
	    }
	    return scalar;
	}

    case VNAPROPERTY_MAP:
	{
	    int map;
	    const char **keys = NULL;

	    if ((keys = vnaproperty_keys(root, "{}")) == NULL) {
		_vnaproperty_yaml_error(vymlp, VNAERR_SYSTEM,
			"vnaproperty_keys: %s", strerror(errno));
		return -1;
	    }
	    errno = 0;
	    if ((map = yaml_document_add_mapping(document, NULL,
			    YAML_ANY_MAPPING_STYLE)) == 0) {
		if (errno == 0) {
		    errno = EINVAL;
		}
		_vnaproperty_yaml_error(vymlp, VNAERR_SYSTEM,
			"yaml_document_add_mapping: %s: %s",
			vymlp->vyml_filename, strerror(errno));
		return -1;
	    }
	    for (const char **cpp = keys; *cpp != NULL; ++cpp) {
		char *key = NULL;
		vnaproperty_t *subtree;
		int value;

		if ((key = vnaproperty_quote_key(*cpp)) == NULL) {
		    free((void *)keys);
		    return -1;
		}
		subtree = vnaproperty_get_subtree(root, "%s", key);
		if ((value = _vnaproperty_yaml_export(vymlp, subtree)) == -1) {
		    free((void *)keys);
		    free((void *)key);
		    return -1;
		}
		if (add_mapping_entry(vymlp, map, key, value) == -1) {
		    free((void *)keys);
		    free((void *)key);
		    return -1;
		}
		free((void *)key);
	    }
	    free((void *)keys);
	    return map;
	}

    case VNAPROPERTY_LIST:
	{
	    int sequence;
	    int count = vnaproperty_count(root, "[]");

	    errno = 0;
	    if ((sequence = yaml_document_add_sequence(document, NULL,
			    YAML_BLOCK_SEQUENCE_STYLE)) == 0) {
		if (errno == 0) {
		    errno = EINVAL;
		}
		_vnaproperty_yaml_error(vymlp, VNAERR_SYSTEM,
			"yaml_document_add_sequence: %s: %s",
			vymlp->vyml_filename, strerror(errno));
		return -1;
	    }
	    for (int i = 0; i < count; ++i) {
		vnaproperty_t *subtree;
		int value;

		subtree = vnaproperty_get_subtree(root, "[%d]", i);
		if ((value = _vnaproperty_yaml_export(vymlp, subtree)) == -1) {
		    return -1;
		}
		if (yaml_document_append_sequence_item(document, sequence,
			    value) == 0) {
		    if (errno == 0) {
			errno = EINVAL;
		    }
		    _vnaproperty_yaml_error(vymlp, VNAERR_SYSTEM,
			    "yaml_document_append_sequence_item: %s: %s",
			    vymlp->vyml_filename, strerror(errno));
		    return -1;
		}
	    }
	    return sequence;
	}

    default:
	break;
    }
    abort();
}
