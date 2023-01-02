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
 * vnaproperty_export_yaml_to_file: import YAML from a file pointer
 *   @rootptr:   address of vnaproperty root
 *   @fp:        open file pointer
 *   @filename:  file name for error messages
 *   @error_fn:  optional error reporting function
 *   @error_arg: optional argument to error reporting function
 */
int vnaproperty_export_yaml_to_file(const vnaproperty_t *root, FILE *fp,
	const char *filename, vnaerr_error_fn_t *error_fn, void *error_arg)
{
    vnaproperty_yaml_t vyml;
    yaml_tag_directive_t tags[1];
    yaml_document_t document;
    bool delete_document = false;
    yaml_emitter_t emitter;

    /*
     * Init the vnaproperty_yaml_t structure.
     */
    (void)memset((void *)&vyml, 0, sizeof(vyml));
    vyml.vyml_filename = filename;
    vyml.vyml_error_fn = error_fn;
    vyml.vyml_error_arg = error_arg;

    /*
     * Init the yaml document.
     */
    if (!yaml_document_initialize(&document, NULL, &tags[0], &tags[0], 0, 0)) {
	if (errno == 0) {
	    errno = EINVAL;
	}
	_vnaproperty_yaml_error(&vyml, VNAERR_SYSTEM,
		"yaml_document_initialize: %s: %s",
		vyml.vyml_filename, strerror(errno));
	goto error;
    }
    delete_document = true;
    vyml.vyml_document = &document;

    /*
     * Export properties to YAML.
     */
    if (_vnaproperty_yaml_export(&vyml, root) == -1) {
	goto error;
    }

    /*
     * Write
     */
    if (!yaml_emitter_initialize(&emitter)) {
	if (errno == 0) {
	    errno = EINVAL;
	}
	_vnaproperty_yaml_error(&vyml, VNAERR_SYSTEM,
		"yaml_emitter_initialize: %s: %s",
		vyml.vyml_filename, strerror(errno));
	goto error;
    }
    yaml_emitter_set_output_file(&emitter, fp);
    yaml_emitter_set_encoding(&emitter, YAML_UTF8_ENCODING);
    yaml_emitter_set_canonical(&emitter, 0);
    //yaml_emitter_set_indent(&emitter, 2);
    yaml_emitter_set_width(&emitter, 80);
    yaml_emitter_set_unicode(&emitter, 1);
    yaml_emitter_set_break(&emitter, YAML_ANY_BREAK);
    if (!yaml_emitter_open(&emitter)) {
	if (errno == 0) {
	    errno = EINVAL;
	}
	_vnaproperty_yaml_error(&vyml, VNAERR_SYSTEM,
		"yaml_emitter_open: %s: %s",
		vyml.vyml_filename, emitter.problem);
	goto error;
    }
    delete_document = false;	/* deleted by yaml_emitter_dump */
    if (!yaml_emitter_dump(&emitter, &document)) {
	if (errno == 0) {
	    errno = EINVAL;
	}
	_vnaproperty_yaml_error(&vyml, VNAERR_SYSTEM,
		"yaml_emitter_dump: %s: %s",
		vyml.vyml_filename, emitter.problem);
	goto error;
    }
    if (!yaml_emitter_close(&emitter)) {
	if (errno == 0) {
	    errno = EINVAL;
	}
	_vnaproperty_yaml_error(&vyml, VNAERR_SYSTEM,
		"yaml_emitter_close: %s: %s",
		vyml.vyml_filename, emitter.problem);
	goto error;
    }
    (void)yaml_emitter_delete(&emitter);
    return 0;

error:
    if (delete_document) {
	yaml_document_delete(&document);
    }
    return -1;
}
