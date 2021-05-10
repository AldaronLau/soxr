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
void resampler_close(struct rate * p);
double resampler_delay(struct rate * p);
void resampler_sizes(size_t * shared, size_t * channel);

static void _soxr_deinterleave_f(float * * dest, void const * * src0, size_t n) {
    printf("deinterlieva\n");

    float const * src = *src0;

    memcpy(dest[0], src, n * sizeof(float));

    src = &src[n];
    
    *src0 = src;
}

static void _soxr_interleave_f(void * * dest0,
  float const * const * src, size_t n)
{
    float * dest = *dest0;

    memcpy(dest, src[0], n * sizeof(float));

    dest = &dest[n];

    *dest0 = dest;
}

typedef void (* fn_t)(void);
typedef fn_t control_block_t[9];

#define resampler_create       (*(char const * (*)(void * channel, void * shared, double io_ratio, double scale))p->control_block[7])

typedef void * resampler_t; /* For one channel. */
typedef void * resampler_shared_t; /* Between channels. */

struct soxr {
  double io_ratio;
  soxr_error_t error;

  size_t max_ilen;

  resampler_shared_t shared;
  resampler_t resampler;
  control_block_t control_block;

  void * * channel_ptrs;
  int flushing;
};

#include "filter.h"

soxr_error_t soxr_error(soxr_t p)
{
  return p->error;
}

extern control_block_t _soxr_rate32_cb;

static void soxr_delete0(soxr_t p)
{
    if (p->resampler)
        resampler_close(p->resampler);
    free(p->resampler);

    free(p->channel_ptrs);
    free(p->shared);

    memset(p, 0, sizeof(*p));
}

soxr_t soxr_create(double io_ratio) {
    soxr_t p = calloc(sizeof(*p), 1);

    control_block_t * control_block;

    p->io_ratio = io_ratio;

    control_block = &_soxr_rate32_cb;

    memcpy(&p->control_block, control_block, sizeof(p->control_block));

    p->io_ratio = io_ratio;

    size_t shared_size, channel_size;

    resampler_sizes(&shared_size, &channel_size);
    p->channel_ptrs = calloc(sizeof(*p->channel_ptrs), 1);
    p->shared = calloc(shared_size, 1);
    p->resampler = calloc(channel_size, 1);

    resampler_create(
        p->resampler,
        p->shared,
        p->io_ratio,
        1.0
    );

    return p;
}

double soxr_delay(soxr_t p)
{
  return (p && !p->error)? resampler_delay(p->resampler) : 0;
}

void soxr_delete(soxr_t p)
{
  if (p)
    soxr_delete0(p), free(p);
}

static size_t soxr_input(soxr_t p, void const * in, size_t len)
{
    printf("EFI\n");

    unsigned i;
    if (!p || p->error) return 0;
    if (!in && len) {p->error = "null input buffer pointer"; return 0;}
    if (!len) {
        p->flushing = true;
        return 0;
    }

    for (i = 0; i < 1; ++i) {
        p->channel_ptrs[i] = resampler_input(p->resampler, NULL, len);
    }

    _soxr_deinterleave_f((float **)p->channel_ptrs, &in, len);

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

  _soxr_interleave_f(&out, (float const * const *)p->channel_ptrs, done);

  return done;
}



size_t soxr_output(soxr_t p, void * out, size_t len0)
{
  size_t odone, odone0 = 0, olen = len0, idone;
  bool was_flushing;

  if (!p || p->error) return 0;
  if (!out && len0) {p->error = "null output buffer pointer"; return 0;}

  do {
    odone = soxr_output_no_callback(p, out, olen);
    odone0 += odone;
    break;
  } while (odone || idone || (!was_flushing && p->flushing));
  return odone0;
}



static size_t soxr_i_for_o(soxr_t p, size_t olen, size_t ilen)
{
  size_t result = (size_t)ceil((double)olen * p->io_ratio);
  return min(result, ilen);
}

soxr_error_t soxr_process(soxr_t p,
    void const * in , size_t ilen0, size_t * idone0,
    void       * out, size_t olen , size_t * odone0)
{
  size_t ilen, idone, odone = 0;
  bool flush_requested = false;

  if (!p) return "null pointer";

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
  return p->error;
}


