/* https://bitbucket.org/jpommier/pffft/raw/483453d8f7661058e74aa4e7cf5c27bcd7887e7a/pffft.h
 * with minor changes for libsoxr. */

#if !defined PFFT_MACROS_ONLY

/* Copyright (c) 2013  Julien Pommier ( pommier@modartt.com )

   Based on original fortran 77 code from FFTPACKv4 from NETLIB,
   authored by Dr Paul Swarztrauber of NCAR, in 1985.

   As confirmed by the NCAR fftpack software curators, the following
   FFTPACKv5 license applies to FFTPACKv4 sources. My changes are
   released under the same terms.

   FFTPACK license:

   http://www.cisl.ucar.edu/css/software/fftpack5/ftpk.html

   Copyright (c) 2004 the University Corporation for Atmospheric
   Research ("UCAR"). All rights reserved. Developed by NCAR's
   Computational and Information Systems Laboratory, UCAR,
   www.cisl.ucar.edu.

   Redistribution and use of the Software in source and binary forms,
   with or without modification, is permitted provided that the
   following conditions are met:

   - Neither the names of NCAR's Computational and Information Systems
   Laboratory, the University Corporation for Atmospheric Research,
   nor the names of its sponsors or contributors may be used to
   endorse or promote products derived from this Software without
   specific prior written permission.

   - Redistributions of source code must retain the above copyright
   notices, this list of conditions, and the disclaimer below.

   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions, and the disclaimer below in the
   documentation and/or other materials provided with the
   distribution.

   THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT
   HOLDERS BE LIABLE FOR ANY CLAIM, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE
   SOFTWARE.
*/

/*
   PFFFT : a Pretty Fast FFT.

   This is basically an adaptation of the single precision fftpack
   (v4) as found on netlib taking advantage of SIMD instruction found
   on cpus such as intel x86 (SSE1), powerpc (Altivec), and arm (NEON).

   For architectures where no SIMD instruction is available, the code
   falls back to a scalar version.

   Restrictions:

   - 1D transforms only, with 32-bit single precision.

   - supports only transforms for inputs of length N of the form
   N=(2^a)*(3^b)*(5^c), a >= 5, b >=0, c >= 0 (32, 48, 64, 96, 128,
   144, 160, etc are all acceptable lengths). Performance is best for
   128<=N<=8192.

   - all (float*) pointers in the functions below are expected to
   have an "simd-compatible" alignment, that is 16 bytes on x86 and
   powerpc CPUs.

   You can allocate such buffers with the functions
   pffft_aligned_malloc / pffft_aligned_free (or with stuff like
   posix_memalign..)

*/

#ifndef PFFFT_H
#define PFFFT_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#if PFFFT_DOUBLE
#define float double
#endif

  /* opaque struct holding internal stuff (precomputed twiddle factors)
     this struct can be shared by many threads as it contains only
     read-only data.
  */
  typedef struct PFFFT_Setup PFFFT_Setup;

  /* direction of the transform */
  typedef enum { PFFFT_FORWARD, PFFFT_BACKWARD } pffft_direction_t;

  /* type of transform */
  typedef enum { PFFFT_REAL, PFFFT_COMPLEX } pffft_transform_t;

  /*
     Perform a multiplication of the frequency components of dft_a and
     dft_b and accumulate them into dft_ab. The arrays should have
     been obtained with pffft_transform(.., PFFFT_FORWARD) and should
     *not* have been reordered with pffft_zreorder (otherwise just
     perform the operation yourself as the dft coefs are stored as
     interleaved complex numbers).

     the operation performed is: dft_ab += (dft_a * fdt_b)*scaling

     The dft_a, dft_b and dft_ab pointers may alias.
  */
  void pffft_zconvolve_accumulate(PFFFT_Setup *setup, const float *dft_a, const float *dft_b, float *dft_ab, float scaling);

#undef float

#ifdef __cplusplus
}
#endif

#endif

#endif
