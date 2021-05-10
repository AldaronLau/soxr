/* SoX Resampler Library      Copyright (c) 2007-16 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details. */

/* Resample using a non-interpolated poly-phase FIR with length LEN. */
/* Input must be followed by FIR_LENGTH-1 samples. */


#define _ sum += coefs[j]*at[j], ++j;

static void FUNCTION(stage_t * p, fifo_t * output_fifo)
{
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
        float const * const __restrict coefs = (float *)COEFS + FIR_LENGTH * rem;
        CONVOLVE(FIR_LENGTH);
        output[i] = sum;
    }

    assert(i == num_out);
    fifo_read(&p->fifo, at / p->L, NULL);
    p->at.integer = at % p->L;
  }
}

#undef _
#undef CORE
#undef cc
#undef core
#undef MIDDLE
#undef CONVOLVE
#undef FIR_LENGTH
#undef FUNCTION
