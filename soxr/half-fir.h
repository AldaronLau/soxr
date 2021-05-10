/* SoX Resampler Library      Copyright (c) 2007-16 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details. */

/* Decimate by 2 using a FIR with odd length (LEN). */
/* Input must be preceded and followed by LEN >> 1 samples. */

#define COEFS ((float const *)p->coefs)

#define BEGINNING float sum = input[0] * .5f
#define ____ __ __
#define __ _ _
#define _ sum += (input[-(2*j +1)] + input[(2*j +1)]) * COEFS[j], ++j;
#define END output[i] = sum

static void FUNCTION_H(stage_t * p, fifo_t * output_fifo)
{
  float const * __restrict input = stage_read_p(p);
  int num_in = min(stage_occupancy(p), p->input_size);
  int i, num_out = (num_in + 1) >> 1;
  float * __restrict output = fifo_reserve(output_fifo, num_out);

  for (i = 0; i < num_out; ++i, input += 2) {
    int j = 0;
    BEGINNING; CONVOLVE; END;
  }
  fifo_read(&p->fifo, 2 * num_out, NULL);
}



#undef _
#undef __
#undef ____
#undef BEGINNING
#undef END
#undef COEFS
#undef CONVOLVE
#undef FUNCTION_H
