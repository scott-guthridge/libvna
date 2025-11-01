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
#include <complex.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <yaml.h>
#include "vnacal_internal.h"
#include "vnaproperty_internal.h"

/*
 * format_complex: format a complex number into the given buffer
 *   @value: complex value to format
 *   @precision: digits of precision
 *   @buf: buffer
 *   @size: number of bytes in buffer
 */
static void format_complex(double complex value, int precision,
	char *buf, size_t size)
{
    double real = creal(value);
    double imag = cimag(value);
    size_t s;

    assert(precision >= 1);
    if (precision == VNACAL_MAX_PRECISION) {
	s = snprintf(buf, size, "%+a %+aj", real, imag);
    } else {
	s = snprintf(buf, size, "%+.*e %+.*ej",
		precision - 1, real,
		precision - 1, imag) == -1;
    }
    assert(s < size); /* otherwise result was truncated */
}

/*
 * add_properties: add user-properties
 *   @root: where to place the properties entry
 *   @properties: properties to place
 */
static int add_properties(vnaproperty_t **root, const vnaproperty_t *properties)
{
    vnaproperty_t **anchor;

    if ((anchor = vnaproperty_set_subtree(root, "properties")) == NULL) {
	return -1;
    }
    if (vnaproperty_copy(anchor, properties) == -1) {
	return -1;
    }
    return 0;
}

/*
 * vnacal_save: create or overwrite a calibration file with new data
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @pathname: calibration file name
 *
 *   The pathname and basename parameters work as in vnacal_load except
 *   that the $HOME/{pathname} directory is created if necessary.
 */
int vnacal_save(vnacal_t *vcp, const char *pathname)
{
    FILE *fp = NULL;
    vnaproperty_t *vprpp_root = NULL;
    vnaproperty_t **vprpp_calibrations;
    vnacal_error_term_matrix_t *matrix_list = NULL;
    const int dprecision = vcp->vc_dprecision;
    char buf[81];
    int minor_version = 0;
    int rc = -1;

    if ((fp = fopen(pathname, "w")) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "fopen: %s: %s",
		vcp->vc_filename, strerror(errno));
	return -1;
    }
    free((void *)vcp->vc_filename);
    if ((vcp->vc_filename = strdup(pathname)) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "strdup: %s", strerror(errno));
	goto out;
    }
    if (add_properties(&vprpp_root, vcp->vc_properties) == -1) {
	goto vnaproperty_set_error;
    }
    vprpp_calibrations = vnaproperty_set_subtree(&vprpp_root, "calibrations[]");
    if (vprpp_calibrations == NULL) {
	goto vnaproperty_set_error;
    }
    for (int ci = 0; ci < vcp->vc_calibration_allocation; ++ci) {
	vnacal_calibration_t *calp = vcp->vc_calibration_vector[ci];
	vnacal_layout_t vl;
	vnaproperty_t **vprpp_calibration;
	vnaproperty_t **vprpp_data;

	if (calp == NULL) {
	    continue;
	}
	if (calp->cal_z0_type != VNACAL_Z0_SCALAR && minor_version < 1) {
	    minor_version = 1;
	}
	_vnacal_layout(&vl, calp->cal_type, calp->cal_rows, calp->cal_columns);
	if (_vnacal_build_error_term_list(calp, &vl, &matrix_list) == -1) {
	    goto out;
	}
	if ((vprpp_calibration = vnaproperty_set_subtree(
		vprpp_calibrations, "[+]{}")) == NULL) {
	    goto vnaproperty_set_error;
	}
	if (vnaproperty_set(vprpp_calibration, "name=%s",
		    calp->cal_name) == -1) {
	    goto vnaproperty_set_error;
	}
	if (vnaproperty_set(vprpp_calibration, "type=%s",
		    vnacal_type_to_name(calp->cal_type)) == -1) {
	    goto vnaproperty_set_error;
	}
	if (vnaproperty_set(vprpp_calibration, "rows=%d",
		    calp->cal_rows) == -1) {
	    goto vnaproperty_set_error;
	}
	if (vnaproperty_set(vprpp_calibration, "columns=%d",
		    calp->cal_columns) == -1) {
	    goto vnaproperty_set_error;
	}
	if (vnaproperty_set(vprpp_calibration, "frequencies=%d",
		    calp->cal_frequencies) == -1) {
	    goto vnaproperty_set_error;
	}
	switch (calp->cal_z0_type) {
	case VNACAL_Z0_SCALAR:
	    format_complex(calp->cal_z0, dprecision, buf, sizeof(buf));
	    if (vnaproperty_set(vprpp_calibration, "z0=%s", buf) == -1) {
		goto vnaproperty_set_error;
	    }
	    break;

	case VNACAL_Z0_VECTOR:
	    {
		vnaproperty_t **vprpp;
		int ports = MAX(calp->cal_rows, calp->cal_columns);

		vprpp = vnaproperty_set_subtree(vprpp_calibration, "z0[]");
		if (vprpp == NULL) {
		    goto vnaproperty_set_error;
		}
		for (int port = 0; port < ports; ++port) {
		    format_complex(calp->cal_z0_vector[port],
			    dprecision, buf, sizeof(buf));
		    if (vnaproperty_set(vprpp, "z0[%d]=%s", port, buf) == -1) {
			goto vnaproperty_set_error;
		    }
		}
	    }
	    break;

	case VNACAL_Z0_MATRIX:
	    break;

	default:
	    abort();
	}
	if (add_properties(vprpp_calibration, calp->cal_properties) == -1) {
	    goto vnaproperty_set_error;
	}
	vprpp_data = vnaproperty_set_subtree(vprpp_calibration, "data[]");
	if (vprpp_data == NULL) {
	    goto vnaproperty_set_error;
	}
	for (int findex = 0; findex < calp->cal_frequencies; ++findex) {
	    double f = calp->cal_frequency_vector[findex];
	    vnaproperty_t **vprpp_frequency;

	    if ((vprpp_frequency = vnaproperty_set_subtree(vprpp_data, "[%d]{}",
			    findex)) == NULL) {
		goto vnaproperty_set_error;
	    }
	    if (vnaproperty_set(vprpp_frequency, "f=%.*e",
			vcp->vc_fprecision, f) == -1) {
		goto vnaproperty_set_error;
	    }
	    if (calp->cal_z0_type == VNACAL_Z0_MATRIX) {
		vnaproperty_t **vprpp;
		int ports = MAX(calp->cal_rows, calp->cal_columns);

		vprpp = vnaproperty_set_subtree(vprpp_calibration, "z0[]");
		if (vprpp == NULL) {
		    goto vnaproperty_set_error;
		}
		for (int port = 0; port < ports; ++port) {
		    format_complex(calp->cal_z0_matrix[port][findex],
			    dprecision, buf, sizeof(buf));
		    if (vnaproperty_set(vprpp, "z0[%d]=%s", port, buf) == -1) {
			goto vnaproperty_set_error;
		    }
		}
	    }
	    for (vnacal_error_term_matrix_t *vetmp = matrix_list;
		    vetmp != NULL; vetmp = vetmp->vetm_next) {
		vnaproperty_t **vprpp_term;
		vnacal_error_term_matrix_type_t type = vetmp->vetm_type;
		double complex **matrix = vetmp->vetm_matrix;
		const int rows = vetmp->vetm_rows;
		const int columns = vetmp->vetm_columns;

		if ((vprpp_term = vnaproperty_set_subtree(
			vprpp_frequency, "%s[]", vetmp->vetm_name)) == NULL) {
		    goto vnaproperty_set_error;
		}
		switch (type) {
		case VETM_VECTOR:
		    assert(vetmp->vetm_rows == 1);
		    for (int i = 0; i < columns; ++i) {
			format_complex(matrix[i][findex],
				dprecision, buf, sizeof(buf));
			if (vnaproperty_set(vprpp_term, "[%d]=%s",
				    i, buf) == -1) {
			    goto vnaproperty_set_error;
			}
		    }
		    break;

		case VETM_MATRIX_ND:
		case VETM_MATRIX:
		    for (int row = 0; row < rows; ++row) {
			vnaproperty_t **vprpp_row;

			if ((vprpp_row = vnaproperty_set_subtree(
				vprpp_term, "[%d][]", row)) == NULL) {
			    goto vnaproperty_set_error;
			}
			for (int column = 0; column < columns; ++column) {
			    if (row != column || type != VETM_MATRIX_ND) {
				format_complex((*matrix++)[findex],
					dprecision, buf, sizeof(buf));
				if (vnaproperty_set(vprpp_row, "[%d]=%s",
					    column, buf) == -1) {
				    goto vnaproperty_set_error;
				}
			    } else {
				if (vnaproperty_set(vprpp_row, "[%d]#",
					    column) == -1) {
				    goto vnaproperty_set_error;
				}
			    }
			}
		    }
		    break;

		default:
		    abort();
		}
	    }
	}
	_vnacal_free_error_term_matrices(&matrix_list);
    }
    (void)fprintf(fp, "#VNACal 1.%d\n", minor_version);
    if (vnaproperty_export_yaml_to_file(vprpp_root, fp, pathname,
		vcp->vc_error_fn, vcp->vc_error_arg) == -1) {
	goto out;
    }
    if (fclose(fp) == -1) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "fclose: %s: %s",
		vcp->vc_filename, strerror(errno));
	fp = NULL;
	goto out;
    }
    fp = NULL;
    rc = 0;
    goto out;

vnaproperty_set_error:
    switch (errno) {
    case EINVAL:
    default:
	_vnacal_error(vcp, VNAERR_INTERNAL, "vnacal_save: internal error: %s",
		strerror(errno));
	break;
    case ENOMEM:
	_vnacal_error(vcp, VNAERR_SYSTEM, "malloc: %s", strerror(errno));
	break;
    }
out:
    (void)vnaproperty_delete(&vprpp_root, ".");
    if (fp != NULL) {
	(void)fclose(fp);
    }
    if (rc == -1) {
	free((void *)vcp->vc_filename);
	vcp->vc_filename = NULL;
    }
    return rc;
}
