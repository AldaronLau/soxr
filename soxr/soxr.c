/* SoX Resampler Library      Copyright (c) 2007-18 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details. */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include "soxr.h"
#include "internal.h"

// Prototypes
typedef struct rate rate_t;
float* resampler_input(rate_t * p, float const * samples, size_t n);
void resampler_process(struct rate * p, size_t olen);
float const * resampler_output(struct rate * p, float * samples, size_t * n0);
void resampler_flush(struct rate * p);
void resampler_sizes(size_t * shared, size_t * channel);
char const * resampler_create(void * channel, void * shared, double io_ratio);

struct soxr {
    double io_ratio;

    void* resampler; /* For one channel. */
    void* shared; /* Between channels. */

    int flushing;
};

#include "filter.h"

soxr_t soxr_create(double io_ratio) {
    soxr_t p = calloc(sizeof(*p), 1);

    p->io_ratio = io_ratio;

    size_t shared_size, channel_size;

    resampler_sizes(&shared_size, &channel_size);
    p->shared = calloc(shared_size, 1);
    p->resampler = calloc(channel_size, 1);

    resampler_create(
        p->resampler,
        p->shared,
        p->io_ratio
    );

    return p;
}

static size_t soxr_input(soxr_t p, void const * in, size_t len) {
    printf("SOXR_INPUT\n");

    if (!len) {
        p->flushing = true;
        return 0;
    }

    void* chan = resampler_input(p->resampler, NULL, len);
    memcpy(chan, in, len * sizeof(float));

    return len;
}

static size_t soxr_output(soxr_t p, void * out, size_t len) {
    printf("SOXR_OUTPUT\n");

    float const * src;
    if (p->flushing) {
        resampler_flush(p->resampler);
    }
    resampler_process(p->resampler, len);
    src = resampler_output(p->resampler, NULL, &len);

    memcpy(out, (float const *) src, len * sizeof(float));

    return len;
}

void soxr_process(soxr_t p,
    void const * in , size_t ilen0, size_t * idone0,
    void       * out, size_t olen , size_t * odone0)
{
    size_t ilen, idone, odone = 0;
    bool flush_requested = false;

    if (!in) {
        flush_requested = true;
        ilen = 0;
        ilen0 = 0;
    } else {
        ilen = ilen0;
    }
    p->flushing |= ilen == ilen0 && flush_requested;

    idone = ilen? soxr_input (p, in , ilen) : 0;
    odone = soxr_output(p, out, olen);

    if (idone0) *idone0 = idone;
    if (odone0) *odone0 = odone;
}
