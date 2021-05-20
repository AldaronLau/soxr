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

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#include "cr.h"

/* --------------------------- Type declarations ---------------------------- */

// Resampler for a single channel.
typedef struct {
    rate_t resampler; /* For one channel. */
    rate_shared_t shared; /* Between channels. */
} resampler_t;

typedef char const * soxr_error_t;                /* 0:no-error; non-0:error. */

typedef void       * soxr_buf_t;  /* 1 buffer of channel-interleaved samples. */


/* --------------------------- API main functions --------------------------- */

/* Create a stream resampler: */

resampler_t* soxr_create(
    double      io_rate      /* Input / Output sample-rate. */
);
