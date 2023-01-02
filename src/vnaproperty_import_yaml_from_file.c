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

/*
 * vnaproperty_import_yaml_from_file: import YAML from a file pointer
 *   @rootptr:   address of vnaproperty root
 *   @fp:        open file pointer
 *   @filename:  file name for error messages
 *   @error_fn:  optional error reporting function
 *   @error_arg: optional argument to error reporting function
 */
int vnaproperty_import_yaml_from_file(vnaproperty_t **rootptr, FILE *fp,
	const char *filename, vnaerr_error_fn_t *error_fn, void *error_arg)
{
    vnaproperty_yaml_t vyml;
    yaml_parser_t parser;
    yaml_document_t document;
    bool delete_document = false;
    yaml_node_t *root;

    (void)memset((void *)&vyml, 0, sizeof(vyml));
    vyml.vyml_filename = filename;
    vyml.vyml_error_fn = error_fn;
    vyml.vyml_error_arg = error_arg;

    yaml_parser_initialize(&parser);
    yaml_parser_set_input_file(&parser, fp);
    if (!yaml_parser_load(&parser, &document)) {
	_vnaproperty_yaml_error(&vyml, VNAERR_SYNTAX, "%s (line %ld) error: %s",
		vyml.vyml_filename, (long)parser.problem_mark.line + 1,
		parser.problem);
	goto error;
    }
    delete_document = true;
    vyml.vyml_document = &document;
    if ((root = yaml_document_get_root_node(&document)) == NULL) {
	_vnaproperty_yaml_error(&vyml, VNAERR_SYNTAX,
		"%s error: empty YAML document", vyml.vyml_filename);
	goto error;
    }
    if (_vnaproperty_yaml_import(&vyml, rootptr, (void *)root) == -1) {
	goto error;
    }
    yaml_document_delete(&document);
    return 0;

error:
    if (delete_document) {
	yaml_document_delete(&document);
    }
    return -1;
}
