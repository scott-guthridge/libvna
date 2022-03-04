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

#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <vnaproperty.h>

/*
 * errfn: report YAML errors
 */
static void errfn(const char *message, void *error_arg,
	vnaerr_category_t category)
{
    fprintf(stderr, "%s\n", message);
}

/*
 * indent: indent level steps
 *   @level: number of indents
 */
static void indent(int level)
{
    for (int i = 0; i < level; ++i) {
	printf("  ");
    }
}

/*
 * print_subtree: print the subtree at root
 *   @root:  root of subtree to print
 *   @level: current indent level
 */
static void print_subtree(vnaproperty_t *root, int level)
{
    /*
     * Handle NULL subtree.
     */
    if (root == NULL) {
	printf("~");
	return;
    }

    /*
     * Handle each node type...
     */
    switch (vnaproperty_type(root, ".")) {
    case 's': /* scalar */
	printf("\"%s\"", vnaproperty_get(root, "."));
	break;

    case 'm': /* map */
	{
	    const char **keys;

	    /*
	     * Get the list of keys.
	     */
	    if ((keys = vnaproperty_keys(root, "{}")) == NULL) {
		fprintf(stderr, "vnaproperty_keys: %s\n",
			strerror(errno));
		exit(5);
	    }

	    /*
	     * For each key, recurse.
	     */
	    printf("{\n");
	    ++level;
	    for (const char **cpp = keys; *cpp != NULL; ++cpp) {
		char *quoted;
		vnaproperty_t *subtree;

		if ((quoted = vnaproperty_quote_key(*cpp)) == NULL) {
		    fprintf(stderr, "vnaproperty_quote_key: %s\n",
			    strerror(errno));
		    exit(6);
		}
		subtree = vnaproperty_get_subtree(root, "%s", quoted);
		free((void *)quoted);
		indent(level);
		printf("%s: ", *cpp);
		print_subtree(subtree, level);
		if (cpp[1] != NULL) {
		    printf(",");
		}
		printf("\n");
	    }
	    --level;
	    indent(level);
	    printf("}");
	    free((void *)keys);
	}
	break;

    case 'l': /* list */
	{
	    int count;

	    /*
	     * Get the count of elements in the list.
	     */
	    if ((count = vnaproperty_count(root, "[]")) == -1) {
		fprintf(stderr, "vnaproperty_count: %s\n",
			strerror(errno));
		exit(8);
	    }

	    /*
	     * For each element, recurse.
	     */
	    printf("[\n");
	    ++level;
	    for (int i = 0; i < count; ++i) {
		vnaproperty_t *subtree;

		subtree = vnaproperty_get_subtree(root, "[%d]", i);
		indent(level);
		print_subtree(subtree, level);
		if (i + 1 < count) {
		    printf(",");
		}
		printf("\n");
	    }
	    --level;
	    indent(level);
	    printf("]");
	}
	break;
    }
}

/*
 * main
 */
int main(int argc, char **argv)
{
    vnaproperty_t *root = NULL;
    FILE *fp;

    /*
     * Check usage.
     */
    if (argc != 2) {
	fprintf(stderr, "usage: yaml-file\n");
	exit(2);
    }

    /*
     * Build the property tree from the input file.
     */
    if ((fp = fopen(argv[1], "r")) == NULL) {
	fprintf(stderr, "fopen: %s: %s\n", argv[1], strerror(errno));
	exit(3);
    }
    if (vnaproperty_import_yaml_from_file(&root, fp, argv[1],
		errfn, NULL) == -1) {
	fprintf(stderr, "vnaproperty_import_yaml_from_file: %s: %s\n",
		argv[1], strerror(errno));
	exit(4);
    }
    fclose(fp);

    /*
     * Print the tree.
     */
    print_subtree(root, 0);
    printf("\n");

    /*
     * Free all objects.
     */
    (void)vnaproperty_delete(&root, ".");

    exit(0);
}
