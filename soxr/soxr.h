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

/* --------------------------- Type declarations ---------------------------- */

typedef struct soxr * soxr_t;          /* A resampler for 1 or more channels. */
typedef char const * soxr_error_t;                /* 0:no-error; non-0:error. */

typedef void       * soxr_buf_t;  /* 1 buffer of channel-interleaved samples. */


/* --------------------------- API main functions --------------------------- */

/* Create a stream resampler: */

soxr_t soxr_create(
    double      io_rate      /* Input / Output sample-rate. */
);

/* If not using an app-supplied input function, after creating a stream
 * resampler, repeatedly call: */

void soxr_process(
    soxr_t      resampler,      /* As returned by soxr_create. */
                            /* Input (to be resampled): */
    void const* in,             /* Input buffer(s); may be NULL (see below). */
    size_t      ilen,           /* Input buf. length (samples per channel). */
    size_t      * idone,        /* To return actual # samples used (<= ilen). */
                            /* Output (resampled): */
    void*       out,            /* Output buffer(s).*/
    size_t      olen,           /* Output buf. length (samples per channel). */
    size_t      * odone);       /* To return actual # samples out (<= olen).

    Note that no special meaning is associated with ilen or olen equal to
    zero.  End-of-input (i.e. no data is available nor shall be available)
    may be indicated by seting `in' to NULL.                                  */

#endif
