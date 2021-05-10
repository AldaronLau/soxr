/* SoX Resampler Library      Copyright (c) 2007-18 robs@users.sourceforge.net
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */



/* -------------------------------- Gubbins --------------------------------- */

#if !defined soxr_included
#define soxr_included


  #include <stddef.h>

typedef struct soxr_io_spec soxr_io_spec_t;

/* ---------------------------- API conventions --------------------------------

Buffer lengths (and occupancies) are expressed as the number of contained
samples per channel.

Parameter names for buffer lengths have the suffix `len'.

A single-character `i' or 'o' is often used in names to give context as
input or output (e.g. ilen, olen).                                            */



/* --------------------------- Version management --------------------------- */

/* E.g. #if SOXR_THIS_VERSION >= SOXR_VERSION(0,1,1) ...                      */

#define SOXR_VERSION(x,y,z)     (((x)<<16)|((y)<<8)|(z))
#define SOXR_THIS_VERSION       SOXR_VERSION(0,1,3)
#define SOXR_THIS_VERSION_STR               "0.1.3"



/* --------------------------- Type declarations ---------------------------- */

typedef struct soxr * soxr_t;          /* A resampler for 1 or more channels. */
typedef char const * soxr_error_t;                /* 0:no-error; non-0:error. */

typedef void       * soxr_buf_t;  /* 1 buffer of channel-interleaved samples. */
typedef void const * soxr_cbuf_t;                        /* Ditto; read-only. */

typedef soxr_buf_t const  * soxr_bufs_t;/* Or, a separate buffer for each ch. */
typedef soxr_cbuf_t const * soxr_cbufs_t;                /* Ditto; read-only. */

typedef void const * soxr_in_t;      /* Either a soxr_cbuf_t or soxr_cbufs_t,
                                        depending on itype in soxr_io_spec_t. */
typedef void       * soxr_out_t;     /* Either a soxr_buf_t or soxr_bufs_t,
                                        depending on otype in soxr_io_spec_t. */



/* --------------------------- API main functions --------------------------- */

#define soxr_strerror(e)               /* Soxr counterpart to strerror. */     \
    ((e)?(e):"no error")


/* Create a stream resampler: */

soxr_t soxr_create(
    double      input_rate,      /* Input sample-rate. */
    double      output_rate,     /* Output sample-rate. */
        /* All following arguments are optional (may be set to NULL). */
    soxr_error_t *              /* To report any error during creation. */
);

/* If not using an app-supplied input function, after creating a stream
 * resampler, repeatedly call: */

soxr_error_t soxr_process(
    soxr_t      resampler,      /* As returned by soxr_create. */
                            /* Input (to be resampled): */
    soxr_in_t   in,             /* Input buffer(s); may be NULL (see below). */
    size_t      ilen,           /* Input buf. length (samples per channel). */
    size_t      * idone,        /* To return actual # samples used (<= ilen). */
                            /* Output (resampled): */
    soxr_out_t  out,            /* Output buffer(s).*/
    size_t      olen,           /* Output buf. length (samples per channel). */
    size_t      * odone);       /* To return actual # samples out (<= olen).

    Note that no special meaning is associated with ilen or olen equal to
    zero.  End-of-input (i.e. no data is available nor shall be available)
    may be indicated by seting `in' to NULL.                                  */


/* then repeatedly call: */

size_t /*odone*/ soxr_output(/* Resample and output a block of data.*/
    soxr_t resampler,            /* As returned by soxr_create. */
    soxr_out_t data,             /* App-supplied buffer(s) for resampled data.*/
    size_t olen);                /* Amount of data to output; >= odone. */



/* Common stream resampler operations: */

soxr_error_t soxr_error(soxr_t);   /* Query error status. */
size_t   * soxr_num_clips(soxr_t); /* Query int. clip counter (for R/W). */
double     soxr_delay(soxr_t);  /* Query current delay in output samples.*/
char const * soxr_engine(soxr_t);  /* Query resampling engine name. */

soxr_error_t soxr_clear(soxr_t); /* Ready for fresh signal, same config. */
void         soxr_delete(soxr_t);  /* Free resources. */


/* For variable-rate resampling. See example # 5 for how to create a
 * variable-rate resampler and how to use this function. */

soxr_error_t soxr_set_io_ratio(soxr_t, double io_ratio, size_t slew_len);

#undef SOXR

#endif
