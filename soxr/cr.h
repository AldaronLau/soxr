/* SoX Resampler Library      Copyright (c) 2007-16 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details. */

#if !defined soxr_cr_included
#define soxr_cr_included

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

/* --------------------------- Type declarations ---------------------------- */

typedef char const * soxr_error_t;                /* 0:no-error; non-0:error. */

typedef void       * soxr_buf_t;  /* 1 buffer of channel-interleaved samples. */


/* --------------------------- API main functions --------------------------- */

#define  FIFO_SIZE_T int
#include "fifo.h"

struct stage;
typedef void (* stage_fn_t)(struct stage * input, fifo_t * output);
typedef struct half_fir_info {
  int num_coefs;
  float const * coefs;
  float att;
} half_fir_info_t;

typedef struct {
    float scalar;
    stage_fn_t fn;
} poly_fir1_t;

typedef struct {
    float beta;
    poly_fir1_t interp[1];
} poly_fir_t;

/* Conceptually: coef_p is &coefs[num_phases][fir_len][interp_order+1]: */
#define coef(coef_p, interp_order, fir_len, phase_num, coef_interp_num, fir_coef_num) (coef_p)[\
  (fir_len) * ((interp_order) + 1) * (phase_num) + \
  ((interp_order) + 1) * (fir_coef_num) + \
  ((interp_order) - (coef_interp_num))]

typedef union { /* Int64 in parts */
  #if HAVE_BIGENDIAN
  struct {int32_t ms; uint32_t ls;} parts;
  #else
  struct {uint32_t ls; int32_t ms;} parts;
  #endif
  int64_t all;
} int64p_t;

typedef union { /* Uint64 in parts */
  #if HAVE_BIGENDIAN
  struct {uint32_t ms, ls;} parts;
  #else
  struct {uint32_t ls, ms;} parts;
  #endif
  uint64_t all;
} uint64p_t;

typedef struct {
  int dft_length;
  int num_taps;
  int post_peak;
  float* coefs;
} dft_filter_t;

typedef struct { /* So generated filter coefs may be shared between channels */
  float   * poly_fir_coefs;
  dft_filter_t dft_filter[2];
} rate_shared_t;

typedef double float_step_t; /* Or long double or __float128. */

typedef union { /* Fixed point arithmetic */
  struct {uint64p_t ls; int64p_t ms;} fix;  /* Hi-prec has ~96 bits. */
  float_step_t flt;
} step_t;

#define integer  fix.ms.parts.ms
#define fraction fix.ms.parts.ls
#define whole    fix.ms.all

typedef int core_flags_t;

typedef struct stage {
  int        num;

  /* Common to all stage types: */
  core_flags_t   core_flags;
  stage_fn_t fn;
  fifo_t     fifo;
  int        pre;       /* Number of past samples to store */
  int        pre_post;  /* pre + number of future samples to store */
  int        preload;   /* Number of zero samples to pre-load the fifo */
  double     out_in_ratio; /* For buffer management. */
  int        input_size;
  bool       is_input;

  /* For a stage with variable (run-time generated) filter coefs: */
  rate_shared_t * shared;
  unsigned   dft_filter_num; /* Which, if any, of the 2 DFT filters to use */
  float       * dft_scratch;
  float      * dft_out;
  float const * coefs;

  /* For a stage with variable L/M: */
  step_t     at, step;
  int        L, remM;
  int        n, phase_bits, block_len;
  double     mult, phase0;
} stage_t;

#define stage_occupancy(s) max(0, fifo_occupancy(&(s)->fifo) - (s)->pre_post)
#define stage_read_p(s) ((float *)fifo_read_ptr(&(s)->fifo) + (s)->pre)

typedef enum {rolloff_small, rolloff_medium, rolloff_none} rolloff_t;

typedef struct {
  half_fir_info_t  const * half_firs;
  size_t half_firs_len;
  poly_fir_t const * poly_firs;
} cr_core_t;

// Resampler for a single channel.
typedef struct {
    cr_core_t const * core;
    double            io_ratio;
    int64_t           samples_in;
    int64_t           samples_out;
    int               num_stages;
    int               flushing;
    stage_t*          stages;

    rate_shared_t shared; /* Between channels. */
} resampler_t;

/* Create a stream resampler: */

resampler_t* soxr_create(
    double      io_rate      /* Input / Output sample-rate. */
);

char const * resampler_init(
  resampler_t * const p,                /* Per audio channel.                            */
  rate_shared_t * const shared,    /* Between channels (undergoing same rate change)*/
  double const io_ratio,           /* Input rate divided by output rate.            */
  cr_core_t const * const core);

#endif
