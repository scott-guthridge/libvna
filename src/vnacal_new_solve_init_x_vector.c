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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "archdep.h"

#include <assert.h>
#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_new_internal.h"

/*
 * _vnacal_new_init_x_vector: initialize the error terms to perfect values
 *   @vnssp: solve state structure
 *   @x_vector: vector of unknowns to init
 *   @x_length: length of x_vector (for defensive code)
 */
void _vnacal_new_solve_init_x_vector(vnacal_new_solve_state_t *vnssp,
	double complex *x_vector, int x_length)
{
    vnacal_new_t *vnp = vnssp->vnss_vnp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    int idx = 0;

    switch (VL_TYPE(vlp)) {
    case VNACAL_T8:
    case VNACAL_TE10:
	for (int d = 0; d < VL_TS_TERMS(vlp); ++d) {
	    x_vector[idx++] = 1.0;
	}
	for (int d = 0; d < VL_TI_TERMS(vlp); ++d) {
	    x_vector[idx++] = 0.0;
	}
	for (int d = 0; d < VL_TX_TERMS(vlp); ++d) {
	    x_vector[idx++] = 0.0;
	}
	for (int d = 0; d < VL_TM_TERMS(vlp) - 1; ++d) {
	    x_vector[idx++] = 1.0;
	}
	break;

    case VNACAL_U8:
    case VNACAL_UE10:
	for (int d = 0; d < VL_UM_TERMS(vlp) - 1; ++d) {
	    x_vector[idx++] = 1.0;
	}
	for (int d = 0; d < VL_UI_TERMS(vlp); ++d) {
	    x_vector[idx++] = 0.0;
	}
	for (int d = 0; d < VL_UX_TERMS(vlp); ++d) {
	    x_vector[idx++] = 0.0;
	}
	for (int d = 0; d < VL_US_TERMS(vlp); ++d) {
	    x_vector[idx++] = 1.0;
	}
	break;

    case VNACAL_T16:
	for (int r = 0; r < VL_TS_ROWS(vlp); ++r) {
	    for (int c = 0; c < VL_TS_COLUMNS(vlp); ++c) {
		x_vector[idx++] = (r == c) ? 1.0 : 0.0;
	    }
	}
	for (int r = 0; r < VL_TI_ROWS(vlp); ++r) {
	    for (int c = 0; c < VL_TI_COLUMNS(vlp); ++c) {
		x_vector[idx++] = 0.0;
	    }
	}
	for (int r = 0; r < VL_TX_ROWS(vlp); ++r) {
	    for (int c = 0; c < VL_TX_COLUMNS(vlp); ++c) {
		x_vector[idx++] = 0.0;
	    }
	}
	for (int r = 0; r < VL_TM_ROWS(vlp); ++r) {
	    for (int c = 0; c < VL_TM_COLUMNS(vlp); ++c) {
		if (!(r == 0 && c == 0)) {
		    x_vector[idx++] = (r == c) ? 1.0 : 0.0;
		}
	    }
	}
	break;

    case VNACAL_U16:
	for (int r = 0; r < VL_UM_ROWS(vlp); ++r) {
	    for (int c = 0; c < VL_UM_COLUMNS(vlp); ++c) {
		if (!(r == 0 && c == 0)) {
		    x_vector[idx++] = (r == c) ? 1.0 : 0.0;
		}
	    }
	}
	for (int r = 0; r < VL_UI_ROWS(vlp); ++r) {
	    for (int c = 0; c < VL_UI_COLUMNS(vlp); ++c) {
		x_vector[idx++] = 0.0;
	    }
	}
	for (int r = 0; r < VL_UX_ROWS(vlp); ++r) {
	    for (int c = 0; c < VL_UX_COLUMNS(vlp); ++c) {
		x_vector[idx++] = 0.0;
	    }
	}
	for (int r = 0; r < VL_US_ROWS(vlp); ++r) {
	    for (int c = 0; c < VL_US_COLUMNS(vlp); ++c) {
		x_vector[idx++] = (r == c) ? 1.0 : 0.0;
	    }
	}
	break;

    case VNACAL_UE14:
    case _VNACAL_E12_UE14:
	for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
	    for (int d = 0; d < VL_UM_TERMS(vlp) - 1; ++d) {
		x_vector[idx++] = 1.0;
	    }
	    for (int d = 0; d < VL_UI_TERMS(vlp); ++d) {
		x_vector[idx++] = 0.0;
	    }
	    for (int d = 0; d < VL_UX_TERMS(vlp); ++d) {
		x_vector[idx++] = 0.0;
	    }
	    for (int d = 0; d < VL_US_TERMS(vlp); ++d) {
		x_vector[idx++] = 1.0;
	    }
	}
	break;

    default:
	abort();
    }
    assert(idx == x_length);
}
