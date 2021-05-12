/* SoX Resampler Library      Copyright (c) 2007-16 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details. */

/* Resample using an interpolated poly-phase FIR with length LEN. */
/* Input must be followed by FIR_LENGTH-1 samples. */

#if COEF_INTERP != 1 && COEF_INTERP != 2 && COEF_INTERP != 3
  #error COEF_INTERP
#endif

#define N FIR_LENGTH

#if COEF_INTERP == 1
    #define _ sum += (b*x + a)*in[j], ++j;
#elif COEF_INTERP == 2
    #define _ sum += ((c*x + b)*x + a)*in[j], ++j;
#else
    #define _ sum += (((d*x + c)*x + b)*x + a)*in[j], ++j;
#endif

#define a (coef(COEFS, COEF_INTERP, N, phase, 0,j))
#define b (coef(COEFS, COEF_INTERP, N, phase, 1,j))
#define c (coef(COEFS, COEF_INTERP, N, phase, 2,j))
#define d (coef(COEFS, COEF_INTERP, N, phase, 3,j))

#define BEGINNING float sum = 0
#define END output[i] = sum
#define CORE(n) core(n)



#define floatPrecCore(n) { \
  float_step_t at = p->at.flt; \
  for (i = 0; (int)at < num_in; ++i, at += p->step.flt) { \
    float const * const __restrict in = input + (int)at; \
    float_step_t frac = at - (int)at; \
    int phase = (int)(frac * (1 << PHASE_BITS)); \
    float x = (float)(frac * (1 << PHASE_BITS) - phase); \
    int j = 0; \
    BEGINNING; CONVOLVE(n); END; \
  } \
  fifo_read(&p->fifo, (int)at, NULL); \
  p->at.flt = at - (int)at; } /* Could round to 1 in some cirmcumstances. */

#define stdPrecCore(n) { \
  int64p_t at; at.all = p->at.whole; \
  for (i = 0; at.parts.ms < num_in; ++i, at.all += p->step.whole) { \
    float const * const __restrict in = input + at.parts.ms; \
    uint32_t const frac = at.parts.ls; \
    int phase = (int)(frac >> (32 - PHASE_BITS)); /* high-order bits */ \
    /* Low-order bits, scaled to [0,1): */ \
    float x = (float)((frac << PHASE_BITS) * (1 / MULT32)); \
    int j = 0; \
    BEGINNING; CONVOLVE(n); END; \
  } \
  fifo_read(&p->fifo, at.parts.ms, NULL); \
  p->at.whole = at.parts.ls; }

#define SPCORE stdPrecCore

#define core(n) SPCORE(n)

static void FUNCTION(stage_t * p, fifo_t * output_fifo)
{
  float const * input = stage_read_p(p);
  int num_in = min(stage_occupancy(p), p->input_size);
  int i, max_num_out = 1 + (int)(num_in * p->out_in_ratio);
  float * const __restrict output = fifo_reserve(output_fifo, max_num_out);

  CORE(N);
  assert(max_num_out - i >= 0);
  fifo_trim_by(output_fifo, max_num_out - i);
}



#undef _
#undef a
#undef b
#undef c
#undef d
#undef CORE
#undef cc
#undef core
#undef COEF_INTERP
#undef N
#undef BEGINNING
#undef END
#undef CONVOLVE
#undef FIR_LENGTH
#undef FUNCTION
#undef PHASE_BITS
