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


#if defined __cplusplus
  #include <cstddef>
  extern "C" {
#else
  #include <stddef.h>
#endif

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



/* -------------------------- API type definitions -------------------------- */


// FIXME
#define soxr_datatype_size(x)  /* Returns `sizeof' a soxr_datatype_t sample. */\
  ((unsigned char *)"\4\10\4\2")[(x)&3]




#define SOXR_TPDF              0     /* Applicable only if otype is INT16. */
#define SOXR_NO_DITHER         8u    /* Disable the above. */

#define SOXR_ROLLOFF_SMALL     0u    /* <= 0.01 dB */
#define SOXR_ROLLOFF_MEDIUM    1u    /* <= 0.35 dB */
#define SOXR_ROLLOFF_NONE      2u    /* For Chebyshev bandwidth. */

#define SOXR_HI_PREC_CLOCK     8u  /* Increase `irrational' ratio accuracy. */
#define SOXR_DOUBLE_PRECISION 16u  /* Use D.P. calcs even if precision <= 20. */
#define SOXR_VR               32u  /* Variable-rate resampling. */



                                   /* For `irrational' ratios only: */
#define SOXR_COEF_INTERP_AUTO  0u    /* Auto select coef. interpolation. */
#define SOXR_COEF_INTERP_LOW   2u    /* Man. select: less CPU, more memory. */
#define SOXR_COEF_INTERP_HIGH  3u    /* Man. select: more CPU, less memory. */



/* -------------------------- API type constructors ------------------------- */

/* These functions allow setting of the most commonly-used structure
 * parameters, with other parameters being given default values.  The default
 * values may then be overridden, directly in the structure, if needed.  */

                                  /* The 5 standard qualities found in SoX: */
#define SOXR_QQ                 0   /* 'Quick' cubic interpolation. */
#define SOXR_LQ                 1   /* 'Low' 16-bit with larger rolloff. */
#define SOXR_MQ                 2   /* 'Medium' 16-bit with medium rolloff. */
#define SOXR_HQ                 SOXR_20_BITQ /* 'High quality'. */
#define SOXR_VHQ                SOXR_28_BITQ /* 'Very high quality'. */

#define SOXR_16_BITQ            3
#define SOXR_20_BITQ            4
#define SOXR_24_BITQ            5
#define SOXR_28_BITQ            6
#define SOXR_32_BITQ            7
                                /* Reserved for internal use (to be removed): */
#define SOXR_LSR0Q              8     /* 'Best sinc'. */
#define SOXR_LSR1Q              9     /* 'Medium sinc'. */
#define SOXR_LSR2Q              10    /* 'Fast sinc'. */

#define SOXR_LINEAR_PHASE       0x00
#define SOXR_INTERMEDIATE_PHASE 0x10
#define SOXR_MINIMUM_PHASE      0x30

#define SOXR_STEEP_FILTER       0x40

soxr_io_spec_t soxr_io_spec(void);

/* --------------------------- Advanced use only ---------------------------- */

/* For new designs, the following functions/usage will probably not be needed.
 * They might be useful when adding soxr into an existing design where values
 * for the resampling-rate and/or number-of-channels parameters to soxr_create
 * are not available when that function will be called.  In such cases, the
 * relevant soxr_create parameter(s) can be given as 0, then one or both of the
 * following (as appropriate) later invoked (but prior to calling soxr_process
 * or soxr_output):
 *
 * soxr_set_error(soxr, soxr_set_io_ratio(soxr, io_ratio, 0));
 * soxr_set_error(soxr, soxr_set_num_channels(soxr, num_channels));
 */

soxr_error_t soxr_set_error(soxr_t, soxr_error_t);


#undef SOXR

#if defined __cplusplus
}
#endif

#endif
