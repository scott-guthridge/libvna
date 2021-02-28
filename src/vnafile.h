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

#ifndef VNAFILE_H
#define VNAFILE_H

#include <errno.h>
#include <vnadata.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * VNAFILE_MAX_PRECISION: argument to vnafile_set_fprecision and
 *      vnafile_set_dprecision for hexadecimal floating point
 *
 * Note: must be the same as VNACAL_MAX_PRECISION
 */
#define VNAFILE_MAX_PRECISION	1000

/* vnafile_type_t: file type */
typedef enum vnafile_type {
	/* automatically determine format from the filename */
	VNAFILE_AUTO		= 0,
	/* simple delimited fields */
	VNAFILE_NATIVE		= 1,
	/* touchstone v1 format */
	VNAFILE_TOUCHSTONE1	= 2,
	/* touchstone v2 format */
	VNAFILE_TOUCHSTONE2	= 3,
} vnafile_type_t;

/*
 * vnafile_t: format information for loading/saving network parameters
 */
typedef struct vnafile vnafile_t;

/*
 * vnafile_error_fn_t: error reporting function type
 *   @message: error message without newline
 *   @error_arg: user-supplied argument passed through to error function
 */
typedef void vnafile_error_fn_t(const char *message, void *error_arg);

/*
 * vnafile_alloc: allocate the vnafile_t parameter structure
 *   @error_fn:  optional error reporting function
 *   @error_arg: optional opaque argument passed through to the error function
 */
extern vnafile_t *vnafile_alloc(vnafile_error_fn_t *error_fn, void *error_arg);

/*
 * vnafile_get_file_type: return the file type
 *   @vfp: pointer to the object returned from vnafile_alloc
 */
extern vnafile_type_t vnafile_get_file_type(const vnafile_t *vfp);

/*
 * vnafile_set_file_type: set the file type
 *   @vfp: pointer to the object returned from vnafile_alloc
 *   @type: file type
 *
 *   The default type is VNAFILE_AUTO where the library tries to intuit
 *   the type from the filename.
 */
extern int vnafile_set_file_type(vnafile_t *vfp, vnafile_type_t type);

/*
 * vnafile_get_format: current the format string
 *   @vfp: pointer to the object returned from vnafile_alloc
 */
extern const char *vnafile_get_format(const vnafile_t *vfp);

/*
 * vnafile_set_format: set the format string
 *   @vfp: pointer to the object returned from vnafile_alloc
 *   @format: a comma-separated case-insensitive list of the following:
 *     {S,Z,Y,T,H,G,A,B}[{ri,ma,dB}]
 *     {il,rl}
 *     zin[{ri,ma}]
 *     {prc,prl,src,srl}
 *     vswr
 *
 *   If not set, the default is "SdB".
 */
extern int vnafile_set_format(vnafile_t *vfp, const char *format);

/*
 * vnafile_get_fprecision: get the frequency value precision
 *   @vfp: pointer to the object returned from vnafile_alloc
 */
extern int vnafile_get_fprecision(const vnafile_t *vcp);

/*
 * vnafile_set_fprecision: set the frequency value precision
 *   @vfp: pointer to the object returned from vnafile_alloc
 *   @precision: precision in decimal places (1..n) or VNAFILE_MAX_PRECISION
 */
extern int vnafile_set_fprecision(vnafile_t *vcp, int precision);

/*
 * vnafile_get_dprecision: set the data value precision
 *   @vfp: pointer to the object returned from vnafile_alloc
 */
extern int vnafile_get_dprecision(const vnafile_t *vcp);

/*
 * vnafile_set_dprecision: set the data value precision for vnafile_save
 *   @vfp: pointer to the object returned from vnafile_alloc
 *   @precision: precision in decimal places (1..n) or VNAFILE_MAX_PRECISION
 */
extern int vnafile_set_dprecision(vnafile_t *vcp, int precision);

/*
 * vnafile_load: load network parameters from filename
 *   @vfp: pointer to the object returned from vnafile_alloc
 *   @filename: file to load
 *   @vdp: output data (reshaped as needed)
 */
extern int vnafile_load(vnafile_t *vfp, const char *filename, vnadata_t *vdp);

/*
 * vnafile_load: load network parameters from a file pointer
 *   @vfp: pointer to the object returned from vnafile_alloc
 *   @fp: file pointer
 *   @filename: filename used in error messages and to intuit the file type
 *   @vdp: output data (reshaped as needed)
 */
extern int vnafile_fload(vnafile_t *vfp, FILE *fp, const char *filename,
	vnadata_t *vdp);

/*
 * vnafile_check: test if parameters are valid for save
 *   @vfp: pointer to the object returned from vnafile_alloc
 *   @filename: file to save
 *   @vdp: input data
 */
extern int vnafile_check(vnafile_t *vfp, const char *filename,
	const vnadata_t *vdp);

/*
 * vnafile_save: save network parameters to filename
 *   @vfp: pointer to the object returned from vnafile_alloc
 *   @filename: file to save
 *   @vdp: input data
 */
extern int vnafile_save(vnafile_t *vfp, const char *filename,
	const vnadata_t *vdp);

/*
 * vnafile_fsave: save network parameters to a file pointer
 *   @vfp: pointer to the object returned from vnafile_alloc
 *   @fp: file pointer
 *   @filename: filename used in error messages and to intuit the file type
 *   @vdp: input data
 */
extern int vnafile_fsave(vnafile_t *vfp, FILE *fp, const char *filename,
	const vnadata_t *vdp);

/*
 * vnafile_free: free the object obtained from vnafile_alloc
 *   @vfp: object to free
 */
extern void vnafile_free(vnafile_t *vfp);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* VNAFILE_H */
