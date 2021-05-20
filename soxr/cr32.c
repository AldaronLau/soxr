/* SoX Resampler Library      Copyright (c) 2007-16 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details. */

#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "filter.h"

#include "internal.h"
#include "cr.h"

#include "half-coefs.h"

static half_fir_info_t const half_firs[] = {
    (half_fir_info_t) {
        .num_coefs = 7,
        .coefs = half_fir_coefs_7,
        .att = 120.65f
    },
    (half_fir_info_t) {
        .num_coefs = 8,
        .coefs = half_fir_coefs_8,
        .att = 136.51f
    },
    (half_fir_info_t) {
        .num_coefs = 9,
        .coefs = half_fir_coefs_9,
        .att = 152.32f
    }
};

static void vpoly0(stage_t * p, fifo_t * output_fifo) {
  int num_in = min(stage_occupancy(p), p->input_size);
  if (num_in) {
    float const * input = stage_read_p(p);
    int at = p->at.integer, step = p->step.integer;
    int i, num_out = (num_in * p->L - at + step - 1) / step;
    float * __restrict output = fifo_reserve(output_fifo, num_out);

    for (i = 0; at < num_in * p->L; ++i, at += step) {
        int const div = at / p->L;
        int rem = at % p->L;
        float const * const __restrict at = input + div;
        int j = 0; float sum = 0;
        float const * const __restrict coefs =
        (float *)(float * __restrict)p->shared->poly_fir_coefs + p->n * rem;
        while (j < (p->n)) {
            sum += coefs[j] * at[j], ++j;
        }
        output[i] = sum;
    }

    assert(i == num_out);
    fifo_read(&p->fifo, at / p->L, NULL);
    p->at.integer = at % p->L;
  }
}

static void U100_0(stage_t * p, fifo_t * output_fifo) {
  int num_in = min(stage_occupancy(p), p->input_size);
  if (num_in) {
    int at = p->at.integer, step = p->step.integer;
    int i, num_out = (num_in * p->L - at + step - 1) / step;
    float * __restrict output = fifo_reserve(output_fifo, num_out);

    for (i = 0; at < num_in * p->L; ++i, at += step) {
        float sum = 0;
        output[i] = sum;
    }

    assert(i == num_out);
    fifo_read(&p->fifo, at / p->L, NULL);
    p->at.integer = at % p->L;
  }
}

static void u100_0(stage_t * p, fifo_t * output_fifo) {
  int num_in = min(stage_occupancy(p), p->input_size);
  if (num_in) {
    int at = p->at.integer, step = p->step.integer;
    int i, num_out = (num_in * p->L - at + step - 1) / step;
    float * __restrict output = fifo_reserve(output_fifo, num_out);

    for (i = 0; at < num_in * p->L; ++i, at += step) {
        float sum = 0;
        output[i] = sum;
    }

    assert(i == num_out);
    fifo_read(&p->fifo, at / p->L, NULL);
    p->at.integer = at % p->L;
  }
}

// Only vpoly0 is used it looks like
static poly_fir_t const poly_firs[] = {
    (poly_fir_t) {-1, {{0, vpoly0}}},
    (poly_fir_t) {-1, {{0, vpoly0}}},
    (poly_fir_t) {-1, {{0, vpoly0}}},
    (poly_fir_t) {-1, {{0, vpoly0}}},
    (poly_fir_t) {-1, {{0, vpoly0}}},
    (poly_fir_t) {-1, {{0, vpoly0}}},

    (poly_fir_t) {-1, {{0, vpoly0}}},
    (poly_fir_t) {-1, {{0, vpoly0}}},
    (poly_fir_t) {-1, {{0, vpoly0}}},
    (poly_fir_t) {-1, {{0, vpoly0}}},
    (poly_fir_t) {-1, {{0, vpoly0}}},
    (poly_fir_t) {-1, {{0, vpoly0}}},

    (poly_fir_t) {10.62f, {{U100_l, U100_0}}},
    (poly_fir_t) {11.28f, {{11, u100_0}}},

    (poly_fir_t) {-1, {{0, vpoly0}}},
    (poly_fir_t) {-1, {{0, vpoly0}}},
    (poly_fir_t) {-1, {{0, vpoly0}}},
    (poly_fir_t) {-1, {{0, vpoly0}}},
    (poly_fir_t) {-1, {{0, vpoly0}}},
};

static cr_core_t const cr_core = (cr_core_t) {
    .half_firs = half_firs,
    .half_firs_len = array_length(half_firs),
    .poly_firs = poly_firs
};

char const* resampler_create(void * channel, void * shared, double io_ratio) {
    return resampler_init(channel, shared, io_ratio, &cr_core);
}
