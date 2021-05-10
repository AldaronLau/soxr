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

typedef void (* fn_t)(void);

typedef void * resampler_t; /* For one channel. */
typedef void * resampler_shared_t; /* Between channels. */

struct soxr {
  double io_ratio;

  size_t max_ilen;

  resampler_shared_t shared;
  resampler_t resampler;

  void * * channel_ptrs;
  int flushing;
};

#include "filter.h"

soxr_t soxr_create(double io_ratio) {
    soxr_t p = calloc(sizeof(*p), 1);

    p->io_ratio = io_ratio;

    size_t shared_size, channel_size;

    resampler_sizes(&shared_size, &channel_size);
    p->channel_ptrs = calloc(sizeof(*p->channel_ptrs), 1);
    p->shared = calloc(shared_size, 1);
    p->resampler = calloc(channel_size, 1);

    resampler_create(
        p->resampler,
        p->shared,
        p->io_ratio
    );

    return p;
}

static size_t soxr_input(soxr_t p, void const * in, size_t len)
{
    printf("EFI\n");

    if (!len) {
        p->flushing = true;
        return 0;
    }

    p->channel_ptrs[0] = resampler_input(p->resampler, NULL, len);
    memcpy(p->channel_ptrs[0], in, len * sizeof(float));

    return len;
}

static size_t soxr_output_1ch(soxr_t p, unsigned i, size_t len)
{
  printf("YOYO\n");

  float const * src;
  if (p->flushing)
    resampler_flush(p->resampler);
  resampler_process(p->resampler, len);
  src = resampler_output(p->resampler, NULL, &len);
  p->channel_ptrs[i] = (void /* const */ *)src;
  return len;
}

static size_t soxr_output_no_callback(soxr_t p, soxr_buf_t out, size_t len)
{
  printf("CZSF\n");
  
  unsigned u;
  size_t done = 0;

  for (u = 0; u < 1; ++u) {
    done = soxr_output_1ch(p, u, len);
  }

    memcpy(out, (float const *)p->channel_ptrs[0], done * sizeof(float));

  return done;
}

static size_t soxr_output(soxr_t p, void * out, size_t len0) {
    size_t odone, odone0 = 0, olen = len0;

    odone = soxr_output_no_callback(p, out, olen);
    odone0 += odone;

    return odone0;
}

static size_t soxr_i_for_o(soxr_t p, size_t olen, size_t ilen)
{
  size_t result = (size_t)ceil((double)olen * p->io_ratio);
  return min(result, ilen);
}

void soxr_process(soxr_t p,
    void const * in , size_t ilen0, size_t * idone0,
    void       * out, size_t olen , size_t * odone0)
{
  size_t ilen, idone, odone = 0;
  bool flush_requested = false;

  if (!in) {
    flush_requested = true, ilen = ilen0 = 0;
  } else {
    if ((ptrdiff_t)ilen0 < 0)
      flush_requested = true, ilen0 = ~ilen0;
    if (idone0 && (1 || flush_requested))
      ilen = soxr_i_for_o(p, olen, ilen0);
    else
      ilen = ilen0/*, olen = soxr_o_for_i(p, ilen, olen)*/;
  }
  p->flushing |= ilen == ilen0 && flush_requested;

    idone = ilen? soxr_input (p, in , ilen) : 0;
    odone = soxr_output(p, out, olen);

  if (idone0) *idone0 = idone;
  if (odone0) *odone0 = odone;
}
