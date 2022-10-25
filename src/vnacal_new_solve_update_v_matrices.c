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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_new_internal.h"

/*
 * update_v_t8: update the v matrices for error term type T8/TE10
 *   @vnssp:    solve state structure
 *   @idx:      index of measured standard
 *   @x_vector: current vector of error terms
 *   @x_length: length of x_vector
 */
static int update_v_t8(vnacal_new_solve_state_t *vnssp, int idx,
	const double complex *x_vector, int x_length)
{
    vnacal_new_t *vnp = vnssp->vnss_vnp;
    vnacal_new_msv_matrices_t *vnmmp = &vnssp->vnss_msv_matrices[idx];
    vnacal_new_measurement_t *vnmp = vnmmp->vnmm_vnmp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int m_rows    = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    double complex vi_matrix[m_columns * m_columns];
    int base = 0;
    int tx_base, tm_base;

    /*
     * Dimensions:
     *   ts: m_rows    x m_columns (diagonal)
     *   ti: m_rows    x m_columns (diagonal)
     *   tx: m_columns x m_columns (diagonal)
     *   tm: m_columns x m_columns (diagonal)
     *   s:  m_columns x m_columns
     */
    assert(m_rows <= m_columns);
    base += m_rows;			/* skip ts diagonal matrix */
    base += m_rows;			/* skip ti diagonal matrix */
    tx_base = base;
    base += m_columns;			/* skip tx diagonal matrix */
    tm_base = base;
    base += m_columns - 1;		/* skip tm diagonal matrix */
    assert(base == x_length);

    /*
     * Find v = (tx s + tm)^-1
     */
    for (int i = 0; i < m_columns; ++i) {
	for (int j = 0; j < m_columns; ++j) {
	    const int vi_cell = i * m_columns + j;
	    const int s_cell = vi_cell;

	    /*
	     * First find vi = tm
	     */
	    if (i != j) {
		vi_matrix[vi_cell] = 0.0;
	    } else if (i == 0) {
		vi_matrix[vi_cell] = 1.0;	/* tm11 */
		--tm_base;
	    } else {
		vi_matrix[vi_cell] = x_vector[tm_base + i];
	    }

	    /*
	     * Add tx * s, treating unknown s values as zero.
	     */
	    if (vnmp->vnm_s_matrix[s_cell] != NULL) {
		vi_matrix[vi_cell] += x_vector[tx_base + i] *
		    vnmmp->vnmm_s_matrix[s_cell];
	    }
	}
    }
    if (_vnacommon_minverse(vnmmp->vnsm_v_matrices[0],
	    vi_matrix, m_columns) == 0.0) {
	return -1;
    }
    return 0;
}

/*
 * update_v_u8: update the v matrices for error term type U8/UE10
 *   @vnssp:    solve state structure
 *   @idx:      index of measured standard
 *   @x_vector: current vector of error terms
 *   @x_length: length of x_vector
 */
static int update_v_u8(vnacal_new_solve_state_t *vnssp, int idx,
	const double complex *x_vector, int x_length)
{
    vnacal_new_t *vnp = vnssp->vnss_vnp;
    vnacal_new_msv_matrices_t *vnmmp = &vnssp->vnss_msv_matrices[idx];
    vnacal_new_measurement_t *vnmp = vnmmp->vnmm_vnmp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int m_rows    = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    double complex vi_matrix[m_rows * m_rows];
    int base = 0;
    int um_base, ux_base;

    /*
     * Dimensions:
     *   um: m_rows x m_rows    (diagonal)
     *   ui: m_rows x m_columns (diagonal)
     *   ux: m_rows x m_rows    (diagonal)
     *   us: m_rows x m_columns (diagonal)
     *   s:  m_rows x m_rows
     */
    assert(m_rows >= m_columns);
    um_base = base;
    base += m_rows - 1;			/* skip um diagonal matrix */
    base += m_columns;			/* skip ui diagonal matrix */
    ux_base = base;
    base += m_rows;			/* skip ux diagonal matrix */
    base += m_columns;			/* skip us diagonal matrix */
    assert(base == x_length);

    /*
     * Find v = (um - s ux)^-1
     */
    for (int i = 0; i < m_rows; ++i) {
	for (int j = 0; j < m_rows; ++j) {
	    const int vi_cell = i * m_rows + j;
	    const int s_cell = vi_cell;

	    /*
	     * First find vi = um
	     */
	    if (i != j) {
		vi_matrix[vi_cell] = 0.0;
	    } else if (i == 0) {
		vi_matrix[vi_cell] = 1.0;	/* tm11 */
		--um_base;
	    } else {
		vi_matrix[vi_cell] = x_vector[um_base + i];
	    }

	    /*
	     * Subract s * ux, treating unknown s values as zero.
	     */
	    if (vnmp->vnm_s_matrix[s_cell] != NULL) {
		vi_matrix[vi_cell] -= vnmmp->vnmm_s_matrix[s_cell] *
		    x_vector[ux_base + j];
	    }
	}
    }
    if (_vnacommon_minverse(vnmmp->vnsm_v_matrices[0],
	    vi_matrix, m_rows) == 0.0) {
	return -1;
    }
    return 0;
}

/*
 * update_v_t16: update the v matrices for error term type T16
 *   @vnssp:    solve state structure
 *   @idx:      index of measured standard
 *   @x_vector: current vector of error terms
 *   @x_length: length of x_vector
 */
static int update_v_t16(vnacal_new_solve_state_t *vnssp, int idx,
	const double complex *x_vector, int x_length)
{
    vnacal_new_t *vnp = vnssp->vnss_vnp;
    vnacal_new_msv_matrices_t *vnmmp = &vnssp->vnss_msv_matrices[idx];
    vnacal_new_measurement_t *vnmp = vnmmp->vnmm_vnmp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int m_rows    = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    double complex vi_matrix[m_columns * m_columns];
    int base = 0;
    int tx_base, tm_base;

    /*
     * Dimensions:
     *   ts: m_rows    x m_columns
     *   ti: m_rows    x m_columns
     *   tx: m_columns x m_columns
     *   tm: m_columns x m_columns
     *   s:  m_columns x m_columns
     */
    assert(m_rows <= m_columns);
    base += m_rows * m_columns;			/* skip ts matrix */
    base += m_rows * m_columns;			/* skip ti matrix */
    tx_base = base;
    base += m_columns * m_columns;		/* skip tx matrix */
    tm_base = base;
    base += m_columns * m_columns - 1;		/* skip tm matrix */
    assert(base == x_length);

    /*
     * Find v = (tx s + tm)^-1
     */
    for (int i = 0; i < m_columns; ++i) {
	for (int j = 0; j < m_columns; ++j) {
	    const int vi_cell = i * m_columns + j;
	    const int tm_cell = vi_cell;

	    /*
	     * First add the contribution of tm.
	     */
	    if (tm_cell == 0) {
		vi_matrix[vi_cell] = 1.0;	/* tm11 */
		--tm_base;
	    } else {
		vi_matrix[vi_cell] = x_vector[tm_base + tm_cell];
	    }

	    /*
	     * Now add the contribution of tx * s.
	     */
	    for (int k = 0; k < m_columns; ++k) {
		const int tx_cell = i * m_columns + k;
		const int s_cell  = k * m_columns + j;

		/*
		 * Add tx * s, treating unknown s values as zero.
		 */
		if (vnmp->vnm_s_matrix[s_cell] != NULL) {
		    vi_matrix[vi_cell] += x_vector[tx_base + tx_cell] *
			vnmmp->vnmm_s_matrix[s_cell];
		}
	    }
	}
    }
    if (_vnacommon_minverse(vnmmp->vnsm_v_matrices[0],
	    vi_matrix, m_columns) == 0.0) {
	return -1;
    }
    return 0;
}

/*
 * update_v_u16: update the v matrices for error term type U16
 *   @vnssp:    solve state structure
 *   @idx:      index of measured standard
 *   @x_vector: current vector of error terms
 *   @x_length: length of x_vector
 */
static int update_v_u16(vnacal_new_solve_state_t *vnssp, int idx,
	const double complex *x_vector, int x_length)
{
    vnacal_new_t *vnp = vnssp->vnss_vnp;
    vnacal_new_msv_matrices_t *vnmmp = &vnssp->vnss_msv_matrices[idx];
    vnacal_new_measurement_t *vnmp = vnmmp->vnmm_vnmp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int m_rows    = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    double complex vi_matrix[m_rows * m_rows];
    int base = 0;
    int um_base, ux_base;

    /*
     * Dimensions:
     *   um: m_rows x m_rows
     *   ui: m_rows x m_columns
     *   ux: m_rows x m_rows
     *   us: m_rows x m_columns
     *   s:  m_rows x m_rows
     */
    assert(m_rows >= m_columns);
    um_base = base;
    base += m_rows * m_rows - 1;		/* skip um matrix */
    base += m_rows * m_columns;			/* skip ui matrix */
    ux_base = base;
    base += m_rows * m_rows;			/* skip ux matrix */
    base += m_rows * m_columns;			/* skip us matrix */
    assert(base == x_length);

    /*
     * Find v = (um - s ux)^-1
     */
    for (int i = 0; i < m_rows; ++i) {
	for (int j = 0; j < m_rows; ++j) {
	    const int vi_cell = i * m_rows + j;
	    const int um_cell = vi_cell;

	    /*
	     * First add the contribution of um.
	     */
	    if (um_cell == 0) {
		vi_matrix[vi_cell] = 1.0;	/* tm11 */
		--um_base;
	    } else {
		vi_matrix[vi_cell] = x_vector[um_base + um_cell];
	    }

	    /*
	     * Next, add the contribution of s * ux.
	     */
	    for (int k = 0; k < m_rows; ++k) {
		const int s_cell  = i * m_rows + k;
		const int ux_cell = k * m_rows + j;

		if (vnmp->vnm_s_matrix[s_cell] != NULL) {
		    vi_matrix[vi_cell] -= vnmmp->vnmm_s_matrix[s_cell] *
			x_vector[ux_base + ux_cell];
		}
	    }

	}
    }
    if (_vnacommon_minverse(vnmmp->vnsm_v_matrices[0],
	    vi_matrix, m_rows) == 0.0) {
	return -1;
    }
    return 0;
}

/*
 * update_v_ue14: update the v matrices for error term type UE14/E12
 *   @vnssp:    solve state structure
 *   @idx:      index of measured standard
 *   @sindex:   index of linear system
 *   @x_vector: current vector of error terms
 *   @x_length: length of x_vector
 */
static double complex update_v_ue14(vnacal_new_solve_state_t *vnssp, int idx,
	int sindex, const double complex *x_vector, int x_length)
{
    vnacal_new_t *vnp = vnssp->vnss_vnp;
    vnacal_new_msv_matrices_t *vnmmp = &vnssp->vnss_msv_matrices[idx];
    vnacal_new_measurement_t *vnmp = vnmmp->vnmm_vnmp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int m_rows = VL_M_ROWS(vlp);
    double complex vi_matrix[m_rows * m_rows];
    int base = 0;
    int um_base, ux_base;

    /*
     * Dimensions:
     *   um: m_rows x m_rows	(diagonal)
     *   ui: m_rows x 1		(diagonal)
     *   ux: m_rows x m_rows	(diagonal)
     *   us: m_rows x 1		(diagonal)
     *   s:  m_rows x m_rows
     */
    um_base = base;
    base += m_rows - 1;			/* skip um diagonal matrix */
    base += 1;				/* skip ui diagonal matrix */
    ux_base = base;
    base += m_rows;			/* skip ux diagonal matrix */
    base += 1;				/* skip us diagonal matrix */

    /*
     * Skip if no v matrix here.
     */
    if (vnmmp->vnsm_v_matrices[sindex] == NULL) {
	return 0;
    }

    /*
     * Find v = (um - s ux)^-1
     */
    for (int i = 0; i < m_rows; ++i) {
	for (int j = 0; j < m_rows; ++j) {
	    const int vi_cell = i * m_rows + j;
	    const int s_cell = vi_cell;

	    /*
	     * First find vi = um
	     */
	    if (i != j) {
		vi_matrix[vi_cell] = 0.0;
	    } else if (i == sindex) {
		vi_matrix[vi_cell] = 1.0;	/* tm11 */
		--um_base;
	    } else {
		vi_matrix[vi_cell] = x_vector[um_base + i];
	    }

	    /*
	     * Subract s * ux, treating unknown s values as zero.
	     */
	    if (vnmp->vnm_s_matrix[s_cell] != NULL) {
		vi_matrix[vi_cell] -= vnmmp->vnmm_s_matrix[s_cell] *
		    x_vector[ux_base + j];
	    }
	}
    }
    if (_vnacommon_minverse(vnmmp->vnsm_v_matrices[sindex],
		vi_matrix, m_rows) == 0.0) {
	return -1;
    }
    assert(base == x_length);
    return 0;
}

/*
 * _vnacal_new_solve_update_v_matrices: update vnsm_v_matrices for 1 system
 *   @function: name of user-called function
 *   @vnssp:    solve state structure
 *   @sindex:   index of linear system
 *   @x_vector: current vector of error terms
 *   @x_lenth:  length of x_vector (for sanity check)
 *
 * Use this version when x_vector holds the error terms for a
 * single linear system.
 *
 * The V matrices convert the residuals of the linear error term systems
 * to errors in their associated measurement.
 *	-Ts S V - Ti V + M (Tx S + Tm) V = 0
 *	V Um M + V Ui - V S (Ux M + Us) = 0
 */
int _vnacal_new_solve_update_v_matrices(const char *function,
	vnacal_new_solve_state_t *vnssp, int sindex,
	const double complex *x_vector, int x_length)
{
    vnacal_new_t *vnp = vnssp->vnss_vnp;
    vnacal_t *vcp = vnp->vn_vcp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;

    switch (VL_TYPE(vlp)) {
    case VNACAL_T8:
    case VNACAL_TE10:
	for (int idx = 0; idx < vnp->vn_measurement_count; ++idx) {
	    if (vnssp->vnss_msv_matrices[idx].vnsm_v_matrices == NULL) {
		continue;
	    }
	    if (vnssp->vnss_msv_matrices[idx].vnsm_v_matrices[0] == NULL) {
		continue;
	    }
	    if (update_v_t8(vnssp, idx, x_vector, x_length) == -1) {
		_vnacal_error(vcp, VNAERR_MATH, "%s: "
			"singular matrix", function);
		return -1;
	    }
	}
	break;

    case VNACAL_U8:
    case VNACAL_UE10:
	for (int idx = 0; idx < vnp->vn_measurement_count; ++idx) {
	    if (vnssp->vnss_msv_matrices[idx].vnsm_v_matrices == NULL) {
		continue;
	    }
	    if (vnssp->vnss_msv_matrices[idx].vnsm_v_matrices[0] == NULL) {
		continue;
	    }
	    if (update_v_u8(vnssp, idx, x_vector, x_length) == -1) {
		_vnacal_error(vcp, VNAERR_MATH, "%s: "
			"singular matrix", function);
		return -1;
	    }
	}
	break;

    case VNACAL_T16:
	for (int idx = 0; idx < vnp->vn_measurement_count; ++idx) {
	    if (vnssp->vnss_msv_matrices[idx].vnsm_v_matrices == NULL) {
		continue;
	    }
	    if (vnssp->vnss_msv_matrices[idx].vnsm_v_matrices[0] == NULL) {
		continue;
	    }
	    if (update_v_t16(vnssp, idx, x_vector, x_length) == -1) {
		_vnacal_error(vcp, VNAERR_MATH, "%s: "
			"singular matrix", function);
		return -1;
	    }
	}
	break;

    case VNACAL_U16:
	for (int idx = 0; idx < vnp->vn_measurement_count; ++idx) {
	    if (vnssp->vnss_msv_matrices[idx].vnsm_v_matrices == NULL) {
		continue;
	    }
	    if (vnssp->vnss_msv_matrices[idx].vnsm_v_matrices[0] == NULL) {
		continue;
	    }
	    if (update_v_u16(vnssp, idx, x_vector, x_length) == -1) {
		_vnacal_error(vcp, VNAERR_MATH, "%s: "
			"singular matrix", function);
		return -1;
	    }
	}
	break;

    case VNACAL_UE14:
    case _VNACAL_E12_UE14:
	for (int idx = 0; idx < vnp->vn_measurement_count; ++idx) {
	    if (vnssp->vnss_msv_matrices[idx].vnsm_v_matrices == NULL) {
		continue;
	    }
	    if (update_v_ue14(vnssp, idx, sindex, x_vector, x_length) == -1) {
		_vnacal_error(vcp, VNAERR_MATH, "%s: "
			"singular matrix", function);
		return -1;
	    }
	}
	break;

    case VNACAL_E12:
    default:
	abort();
    }
    return 0;
}

/*
 * _vnacal_new_solve_update_all_v_matrices: update v_matrices for all systems
 *   @function: name of user-called function
 *   @vnssp:    solve state structure
 *   @x_vector: current vector of error terms
 *   @x_lenth:  length of x_vector (for sanity check)
 *
 * Use this version when x_vector is the concatenation of error terms
 * from one or more linear systems.
 */
int _vnacal_new_solve_update_all_v_matrices(const char *function,
	vnacal_new_solve_state_t *vnssp,
	const double complex *x_vector, int x_length)
{
    vnacal_new_t *vnp = vnssp->vnss_vnp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;

    assert(x_length == vnp->vn_systems * (vlp->vl_t_terms - 1));
    for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
	int base = (vlp->vl_t_terms - 1) * sindex;

	if (_vnacal_new_solve_update_v_matrices(function, vnssp, sindex,
		&x_vector[base], vlp->vl_t_terms - 1) == -1) {
	    return -1;
	}
    }
    return 0;
}
