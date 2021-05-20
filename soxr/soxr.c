/* SoX Resampler Library      Copyright (c) 2007-18 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details. */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include "soxr.h"
#include "internal.h"
#include "cr.h"

// Prototypes
float* resampler_input(rate_t * p, float const * samples, size_t n);
void resampler_process(struct rate * p, size_t olen);
float const * resampler_output(struct rate * p, float * samples, size_t * n0);
void resampler_flush(struct rate * p);

char const * resampler_create(void * channel, void * shared, double io_ratio);

resampler_t* soxr_create(double io_ratio) {
    resampler_t* p = calloc(sizeof(*p), 1);

    resampler_create(
        &p->resampler,
        &p->shared,
        io_ratio
    );

    return p;
}
