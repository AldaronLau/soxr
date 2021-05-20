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

static const float half_fir_coefs_7[] = {
     3.1062656496657370e-01, -8.4998810699955796e-02,  3.4007044621123500e-02,
    -1.2839903789829387e-02,  3.9899380181723145e-03, -8.9355202017945374e-04,
     1.0918292424806546e-04,
};

static const float half_fir_coefs_8[] = {
     3.1154652365332069e-01, -8.7344917685739543e-02,  3.6814458353637280e-02,
    -1.5189204581464479e-02,  5.4540855610738801e-03, -1.5643862626630416e-03,
     3.1816575906323303e-04, -3.4799449225005688e-05,
};

static const float half_fir_coefs_9[] = {
     3.1227034755311189e-01, -8.9221517147969526e-02,  3.9139704015071934e-02,
    -1.7250558515852023e-02,  6.8589440230476112e-03, -2.3045049636430419e-03,
     6.0963740543348963e-04, -1.1323803957431231e-04,  1.1197769991000046e-05,
};

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

    (poly_fir_t) {10.62f, {{42, U100_0}}},
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
