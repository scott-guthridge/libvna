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

#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include "vnadata_internal.h"


/*
 * _vnadata_format_to_name: return the given format descriptor as a string
 *   @vfdp: format descriptor
 */
const char *_vnadata_format_to_name(const vnadata_format_descriptor_t *vfdp)
{
    switch (vfdp->vfd_format) {
    case VNADATA_FORMAT_PRC:
	return "PRC";

    case VNADATA_FORMAT_PRL:
	return "PRL";

    case VNADATA_FORMAT_SRC:
	return "SRC";

    case VNADATA_FORMAT_SRL:
	return "SRL";

    case VNADATA_FORMAT_IL:
	return "IL";

    case VNADATA_FORMAT_RL:
	return "RL";

    case VNADATA_FORMAT_VSWR:
	return "VSWR";

    default:
	switch (vfdp->vfd_parameter) {
	case VPT_UNDEF:
	    switch (vfdp->vfd_format) {
	    case VNADATA_FORMAT_REAL_IMAG:
		return "ri";
	    case VNADATA_FORMAT_MAG_ANGLE:
		return "ma";
	    case VNADATA_FORMAT_DB_ANGLE:
		return "dB";
	    default:
		break;
	    }
	    break;

	case VPT_S:
	    switch (vfdp->vfd_format) {
	    case VNADATA_FORMAT_REAL_IMAG:
		return "Sri";
	    case VNADATA_FORMAT_MAG_ANGLE:
		return "Sma";
	    case VNADATA_FORMAT_DB_ANGLE:
		return "SdB";
	    default:
		break;
	    }
	    break;

	case VPT_T:
	    switch (vfdp->vfd_format) {
	    case VNADATA_FORMAT_REAL_IMAG:
		return "Tri";
	    case VNADATA_FORMAT_MAG_ANGLE:
		return "Tma";
	    case VNADATA_FORMAT_DB_ANGLE:
		return "TdB";
	    default:
		break;
	    }
	    break;

	case VPT_Z:
	    switch (vfdp->vfd_format) {
	    case VNADATA_FORMAT_REAL_IMAG:
		return "Zri";
	    case VNADATA_FORMAT_MAG_ANGLE:
		return "Zma";
	    case VNADATA_FORMAT_DB_ANGLE:
		return "ZdB";
	    default:
		break;
	    }
	    break;

	case VPT_Y:
	    switch (vfdp->vfd_format) {
	    case VNADATA_FORMAT_REAL_IMAG:
		return "Yri";
	    case VNADATA_FORMAT_MAG_ANGLE:
		return "Yma";
	    case VNADATA_FORMAT_DB_ANGLE:
		return "YdB";
	    default:
		break;
	    }
	    break;

	case VPT_H:
	    switch (vfdp->vfd_format) {
	    case VNADATA_FORMAT_REAL_IMAG:
		return "Hri";
	    case VNADATA_FORMAT_MAG_ANGLE:
		return "Hma";
	    case VNADATA_FORMAT_DB_ANGLE:
		return "HdB";
	    default:
		break;
	    }
	    break;

	case VPT_G:
	    switch (vfdp->vfd_format) {
	    case VNADATA_FORMAT_REAL_IMAG:
		return "Gri";
	    case VNADATA_FORMAT_MAG_ANGLE:
		return "Gma";
	    case VNADATA_FORMAT_DB_ANGLE:
		return "GdB";
	    default:
		break;
	    }
	    break;

	case VPT_A:
	    switch (vfdp->vfd_format) {
	    case VNADATA_FORMAT_REAL_IMAG:
		return "Ari";
	    case VNADATA_FORMAT_MAG_ANGLE:
		return "Ama";
	    case VNADATA_FORMAT_DB_ANGLE:
		return "AdB";
	    default:
		break;
	    }
	    break;

	case VPT_B:
	    switch (vfdp->vfd_format) {
	    case VNADATA_FORMAT_REAL_IMAG:
		return "Bri";
	    case VNADATA_FORMAT_MAG_ANGLE:
		return "Bma";
	    case VNADATA_FORMAT_DB_ANGLE:
		return "BdB";
	    default:
		break;
	    }
	    break;

	case VPT_ZIN:
	    switch (vfdp->vfd_format) {
	    case VNADATA_FORMAT_REAL_IMAG:
		return "Zinri";
	    case VNADATA_FORMAT_MAG_ANGLE:
		return "Zinma";
	    default:
		break;
	    }
	    break;

	default:
	    break;
	}
    }
    abort();
    /*NOTREACHED*/
}
