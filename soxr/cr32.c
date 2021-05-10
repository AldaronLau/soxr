/* SoX Resampler Library      Copyright (c) 2007-16 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details. */

#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "filter.h"

#include "internal.h"
#include "cr.h"

extern fn_t _soxr_rdft32_cb[];

static void cubic_stage_fn(stage_t * p, fifo_t * output_fifo)
{
  float const * input = stage_read_p(p);
  int num_in = min(stage_occupancy(p), p->input_size);
  int i, max_num_out = 1 + (int)(num_in * p->out_in_ratio);
  float * output = fifo_reserve(output_fifo, max_num_out);

  for (i = 0; p->at.integer < num_in; ++i, p->at.whole += p->step.whole) {
    float const * s = input + p->at.integer;
    double x = p->at.fraction * (1 / MULT32);
    double b = .5*(s[1]+s[-1])-*s, a = (1/6.)*(s[2]-s[1]+s[-1]-*s-4*b);
    double c = s[1]-*s-a-b;
    output[i] = (float)(p->mult * (((a*x + b)*x + c)*x + *s));
  }
  assert(max_num_out - i >= 0);
  fifo_trim_by(output_fifo, max_num_out - i);
  fifo_read(&p->fifo, p->at.integer, NULL);
  p->at.integer = 0;
}



#if defined __AVX__
  #define DEFINED_AVX 1
#else
  #define DEFINED_AVX 0
#endif

#if defined __x86_64__ || defined _M_X64 || defined i386 || defined _M_IX86
  #define DEFINED_X86 1
#else
  #define DEFINED_X86 0
#endif

#if defined __arm__
  #define DEFINED_ARM 1
#else
  #define DEFINED_ARM 0
#endif

#include "half-coefs.h"

#define FUNCTION_H h7
#define CONVOLVE ____ __ _
#include "half-fir.h"

#define FUNCTION_H h8
#define CONVOLVE ____ ____
#include "half-fir.h"

#define FUNCTION_H h9
#define CONVOLVE ____ ____ _
#include "half-fir.h"

static half_fir_info_t const half_firs[] = {
  { 7, half_fir_coefs_7 , h7 , 0  , 120.65f},
  { 8, half_fir_coefs_8 , h8 , 0  , 136.51f},
  { 9, half_fir_coefs_9 , h9 , 0  , 152.32f},
};

#define COEFS (float * __restrict)p->shared->poly_fir_coefs
#define VAR_LENGTH p->n
#define VAR_CONVOLVE(n) while (j < (n)) _
#define VAR_POLY_PHASE_BITS p->phase_bits



#define FUNCTION vpoly0
#define FIR_LENGTH VAR_LENGTH
#define CONVOLVE(n) VAR_CONVOLVE(n)
#include "poly-fir0.h"

#define FUNCTION vpoly1
#define COEF_INTERP 1
#define PHASE_BITS VAR_POLY_PHASE_BITS
#define FIR_LENGTH VAR_LENGTH
#define CONVOLVE(n) VAR_CONVOLVE(n)
#include "poly-fir.h"

#define FUNCTION vpoly2
#define COEF_INTERP 2
#define PHASE_BITS VAR_POLY_PHASE_BITS
#define FIR_LENGTH VAR_LENGTH
#define CONVOLVE(n) VAR_CONVOLVE(n)
#include "poly-fir.h"

#define FUNCTION vpoly3
#define COEF_INTERP 3
#define PHASE_BITS VAR_POLY_PHASE_BITS
#define FIR_LENGTH VAR_LENGTH
#define CONVOLVE(n) VAR_CONVOLVE(n)
#include "poly-fir.h"




#define poly_fir_convolve_U100 _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
#define FUNCTION U100_0
#define FIR_LENGTH U100_l
#define CONVOLVE(n) poly_fir_convolve_U100
#include "poly-fir0.h"

#define u100_l 11
#define poly_fir_convolve_u100 _ _ _ _ _ _ _ _ _ _ _
#define FUNCTION u100_0
#define FIR_LENGTH u100_l
#define CONVOLVE(n) poly_fir_convolve_u100
#include "poly-fir0.h"

#define FUNCTION u100_1
#define COEF_INTERP 1
#define PHASE_BITS 8
#define FIR_LENGTH u100_l
#define CONVOLVE(n) poly_fir_convolve_u100
#include "poly-fir.h"

#define FUNCTION u100_2
#define COEF_INTERP 2
#define PHASE_BITS 6
#define FIR_LENGTH u100_l
#define CONVOLVE(n) poly_fir_convolve_u100
#include "poly-fir.h"

static poly_fir_t const poly_firs[] = {
  {-1, {{0, vpoly0}, { 7.2f, vpoly1}, {5.0f, vpoly2}}},
  {-1, {{0, vpoly0}, { 9.4f, vpoly1}, {6.7f, vpoly2}}},
  {-1, {{0, vpoly0}, {12.4f, vpoly1}, {7.8f, vpoly2}}},
  {-1, {{0, vpoly0}, {13.6f, vpoly1}, {9.3f, vpoly2}}},
  {-1, {{0, vpoly0}, {10.5f, vpoly2}, {8.4f, vpoly3}}},
  {-1, {{0, vpoly0}, {11.85f,vpoly2}, {9.0f, vpoly3}}},

  {-1, {{0, vpoly0}, { 8.0f, vpoly1}, {5.3f, vpoly2}}},
  {-1, {{0, vpoly0}, { 8.6f, vpoly1}, {5.7f, vpoly2}}},
  {-1, {{0, vpoly0}, {10.6f, vpoly1}, {6.75f,vpoly2}}},
  {-1, {{0, vpoly0}, {12.6f, vpoly1}, {8.6f, vpoly2}}},
  {-1, {{0, vpoly0}, { 9.6f, vpoly2}, {7.6f, vpoly3}}},
  {-1, {{0, vpoly0}, {11.4f, vpoly2}, {8.65f,vpoly3}}},

  {10.62f, {{U100_l, U100_0}, {0, 0}, {0, 0}}},
  {11.28f, {{u100_l, u100_0}, {8, u100_1}, {6, u100_2}}},

  {-1, {{0, vpoly0}, {   9, vpoly1}, {  6, vpoly2}}},
  {-1, {{0, vpoly0}, {  11, vpoly1}, {  7, vpoly2}}},
  {-1, {{0, vpoly0}, {  13, vpoly1}, {  8, vpoly2}}},
  {-1, {{0, vpoly0}, {  10, vpoly2}, {  8, vpoly3}}},
  {-1, {{0, vpoly0}, {  12, vpoly2}, {  9, vpoly3}}},
};

static cr_core_t const cr_core = {
  {malloc, calloc, free},
  half_firs, array_length(half_firs),
  0, 0,
  cubic_stage_fn,
  poly_firs, _soxr_rdft32_cb
};

char const * rate_create(void * channel, void * shared, double io_ratio,
    double scale)
{
  return resampler_init(channel, shared, io_ratio, scale, &cr_core);
}

fn_t _soxr_rate32_cb[] = {
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  (fn_t)rate_create,
  (fn_t)0,
};

