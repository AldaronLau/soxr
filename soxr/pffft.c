/* https://bitbucket.org/jpommier/pffft/raw/483453d8f7661058e74aa4e7cf5c27bcd7887e7a/pffft.c
 * with minor changes for libsoxr. */

/* Copyright (c) 2013  Julien Pommier ( pommier@modartt.com )

   Based on original fortran 77 code from FFTPACKv4 from NETLIB
   (http://www.netlib.org/fftpack), authored by Dr Paul Swarztrauber
   of NCAR, in 1985.

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


   PFFFT : a Pretty Fast FFT.

   This file is largerly based on the original FFTPACK implementation, modified in
   order to take advantage of SIMD instructions of modern CPUs.
*/

/*
  ChangeLog:
  - 2011/10/02, version 1: This is the very first release of this file.
*/

#include "pffft.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

/* detect compiler flavour */
# define COMPILER_GCC

#if defined(COMPILER_GCC)
#  define ALWAYS_INLINE(return_type) inline return_type __attribute__ ((always_inline))
#  define NEVER_INLINE(return_type) return_type __attribute__ ((noinline))
#  define RESTRICT __restrict
#  define VLA_ARRAY_ON_STACK(type__, varname__, size__) type__ varname__[size__];
#endif


/*
   vector support macros: the rest of the code is independant of
   SSE/Altivec/NEON -- adding support for other platforms with 4-element
   vectors should be limited to these macros
*/


/* define PFFFT_SIMD_DISABLE if you want to use scalar code instead of simd code */
/*#define PFFFT_SIMD_DISABLE */

/*
   Altivec support macros
*/
#if !defined(PFFFT_SIMD_DISABLE) && (defined(__ppc__) || defined(__ppc64__))
typedef vector float v4sf;
#  define SIMD_SZ 4
#  define VZERO() ((vector float) vec_splat_u8(0))
#  define VMUL(a,b) vec_madd(a,b, VZERO())
#  define VADD(a,b) vec_add(a,b)
#  define VMADD(a,b,c) vec_madd(a,b,c)
#  define VSUB(a,b) vec_sub(a,b)
inline v4sf ld_ps1(const float *p) { v4sf v=vec_lde(0,p); return vec_splat(vec_perm(v, v, vec_lvsl(0, p)), 0); }
#  define LD_PS1(p) ld_ps1(&p)
#  define INTERLEAVE2(in1, in2, out1, out2) { v4sf tmp__ = vec_mergeh(in1, in2); out2 = vec_mergel(in1, in2); out1 = tmp__; }
#  define UNINTERLEAVE2(in1, in2, out1, out2) {                           \
    vector unsigned char vperm1 =  (vector unsigned char)(0,1,2,3,8,9,10,11,16,17,18,19,24,25,26,27); \
    vector unsigned char vperm2 =  (vector unsigned char)(4,5,6,7,12,13,14,15,20,21,22,23,28,29,30,31); \
    v4sf tmp__ = vec_perm(in1, in2, vperm1); out2 = vec_perm(in1, in2, vperm2); out1 = tmp__; \
  }
#  define VTRANSPOSE4(x0,x1,x2,x3) {              \
    v4sf y0 = vec_mergeh(x0, x2);               \
    v4sf y1 = vec_mergel(x0, x2);               \
    v4sf y2 = vec_mergeh(x1, x3);               \
    v4sf y3 = vec_mergel(x1, x3);               \
    x0 = vec_mergeh(y0, y2);                    \
    x1 = vec_mergel(y0, y2);                    \
    x2 = vec_mergeh(y1, y3);                    \
    x3 = vec_mergel(y1, y3);                    \
  }
#  define VSWAPHL(a,b) vec_perm(a,b, (vector unsigned char)(16,17,18,19,20,21,22,23,8,9,10,11,12,13,14,15))
#  define VALIGNED(ptr) ((((long)(ptr)) & 0xF) == 0)

/*
  SSE1 support macros
*/
#elif !defined(PFFFT_SIMD_DISABLE) && (defined(__x86_64__) || defined(_M_X64) || defined(i386) || defined(_M_IX86))

#  define SIMD_SZ 4 /* 4 floats by simd vector -- this is pretty much hardcoded in the preprocess/finalize functions anyway so you will have to work if you want to enable AVX with its 256-bit vectors. */

#if !PFFFT_DOUBLE
#include <xmmintrin.h>
typedef __m128 v4sf;
#  define VZERO() _mm_setzero_ps()
#  define VMUL(a,b) _mm_mul_ps(a,b)
#  define VADD(a,b) _mm_add_ps(a,b)
#  define VMADD(a,b,c) _mm_add_ps(_mm_mul_ps(a,b), c)
#  define VSUB(a,b) _mm_sub_ps(a,b)
#  define LD_PS1(p) _mm_set1_ps(p)
#  define INTERLEAVE2(in1, in2, out1, out2) { v4sf tmp__ = _mm_unpacklo_ps(in1, in2); out2 = _mm_unpackhi_ps(in1, in2); out1 = tmp__; }
#  define UNINTERLEAVE2(in1, in2, out1, out2) { v4sf tmp__ = _mm_shuffle_ps(in1, in2, _MM_SHUFFLE(2,0,2,0)); out2 = _mm_shuffle_ps(in1, in2, _MM_SHUFFLE(3,1,3,1)); out1 = tmp__; }
#  define VTRANSPOSE4(x0,x1,x2,x3) _MM_TRANSPOSE4_PS(x0,x1,x2,x3)
#  define VSWAPHL(a,b) _mm_shuffle_ps(b, a, _MM_SHUFFLE(3,2,1,0))
#  define VALIGNED(ptr) ((((long)(ptr)) & 0xF) == 0)

#else
#include "pffft-avx.h"
#endif

/*
  ARM NEON support macros
*/
#elif !defined(PFFFT_SIMD_DISABLE) && defined(__arm__)
#  include <arm_neon.h>
typedef float32x4_t v4sf;
#  define SIMD_SZ 4
#  define VZERO() vdupq_n_f32(0)
#  define VMUL(a,b) vmulq_f32(a,b)
#  define VADD(a,b) vaddq_f32(a,b)
#  define VMADD(a,b,c) vmlaq_f32(c,a,b)
#  define VSUB(a,b) vsubq_f32(a,b)
#  define LD_PS1(p) vld1q_dup_f32(&(p))
#  define INTERLEAVE2(in1, in2, out1, out2) { float32x4x2_t tmp__ = vzipq_f32(in1,in2); out1=tmp__.val[0]; out2=tmp__.val[1]; }
#  define UNINTERLEAVE2(in1, in2, out1, out2) { float32x4x2_t tmp__ = vuzpq_f32(in1,in2); out1=tmp__.val[0]; out2=tmp__.val[1]; }
#  define VTRANSPOSE4(x0,x1,x2,x3) {                                    \
    float32x4x2_t t0_ = vzipq_f32(x0, x2);                              \
    float32x4x2_t t1_ = vzipq_f32(x1, x3);                              \
    float32x4x2_t u0_ = vzipq_f32(t0_.val[0], t1_.val[0]);              \
    float32x4x2_t u1_ = vzipq_f32(t0_.val[1], t1_.val[1]);              \
    x0 = u0_.val[0]; x1 = u0_.val[1]; x2 = u1_.val[0]; x3 = u1_.val[1]; \
  }
/* marginally faster version */
/*#  define VTRANSPOSE4(x0,x1,x2,x3) { asm("vtrn.32 %q0, %q1;\n vtrn.32 %q2,%q3\n vswp %f0,%e2\n vswp %f1,%e3" : "+w"(x0), "+w"(x1), "+w"(x2), "+w"(x3)::); } */
#  define VSWAPHL(a,b) vcombine_f32(vget_low_f32(b), vget_high_f32(a))
#  define VALIGNED(ptr) ((((long)(ptr)) & 0x3) == 0)
#else
#  if !defined(PFFFT_SIMD_DISABLE)
#    warning "building with simd disabled !\n";
#    define PFFFT_SIMD_DISABLE /* fallback to scalar code */
#  endif
#endif

#if PFFFT_DOUBLE
#define float double
#endif

/* fallback mode for situations where SSE/Altivec are not available, use scalar mode instead */
#ifdef PFFFT_SIMD_DISABLE
typedef float v4sf;
#  define SIMD_SZ 1
#  define VZERO() 0.f
#  define VMUL(a,b) ((a)*(b))
#  define VADD(a,b) ((a)+(b))
#  define VMADD(a,b,c) ((a)*(b)+(c))
#  define VSUB(a,b) ((a)-(b))
#  define LD_PS1(p) (p)
#  define VALIGNED(ptr) ((((long)(ptr)) & 0x3) == 0)
#endif

/* shortcuts for complex multiplcations */
#define VCPLXMUL(ar,ai,br,bi) { v4sf tmp; tmp=VMUL(ar,bi); ar=VMUL(ar,br); ar=VSUB(ar,VMUL(ai,bi)); ai=VMUL(ai,br); ai=VADD(ai,tmp); }
#define VCPLXMULCONJ(ar,ai,br,bi) { v4sf tmp; tmp=VMUL(ar,bi); ar=VMUL(ar,br); ar=VADD(ar,VMUL(ai,bi)); ai=VMUL(ai,br); ai=VSUB(ai,tmp); }
#ifndef SVMUL
/* multiply a scalar with a vector */
#define SVMUL(f,v) VMUL(LD_PS1(f),v)
#endif

#if !defined PFFT_MACROS_ONLY

#if !defined(PFFFT_SIMD_DISABLE)
typedef union v4sf_union {
  v4sf  v;
  float f[4];
} v4sf_union;

#endif /*!PFFFT_SIMD_DISABLE */

struct PFFFT_Setup {
  int     N;
  int     Ncvec; /* nb of complex simd vectors (N/4 if PFFFT_COMPLEX, N/8 if PFFFT_REAL) */
  int ifac[15];
  pffft_transform_t transform;
  v4sf *data; /* allocated room for twiddle coefs */
  float *e;    /* points into 'data' , N/4*3 elements */
  float *twiddle; /* points into 'data', N/4 elements */
};

#if !defined(PFFFT_SIMD_DISABLE)

static ALWAYS_INLINE(void) pffft_real_finalize_4x4(const v4sf *in0, const v4sf *in1, const v4sf *in,
                            const v4sf *e, v4sf *out) {
  v4sf r0, i0, r1, i1, r2, i2, r3, i3;
  v4sf sr0, dr0, sr1, dr1, si0, di0, si1, di1;
  r0 = *in0; i0 = *in1;
  r1 = *in++; i1 = *in++; r2 = *in++; i2 = *in++; r3 = *in++; i3 = *in++;
  VTRANSPOSE4(r0,r1,r2,r3);
  VTRANSPOSE4(i0,i1,i2,i3);

  /*
    transformation for each column is:

    [1   1   1   1   0   0   0   0]   [r0]
    [1   0  -1   0   0  -1   0   1]   [r1]
    [1   0  -1   0   0   1   0  -1]   [r2]
    [1  -1   1  -1   0   0   0   0]   [r3]
    [0   0   0   0   1   1   1   1] * [i0]
    [0  -1   0   1  -1   0   1   0]   [i1]
    [0  -1   0   1   1   0  -1   0]   [i2]
    [0   0   0   0  -1   1  -1   1]   [i3]
  */

  /*cerr << "matrix initial, before e , REAL:\n 1: " << r0 << "\n 1: " << r1 << "\n 1: " << r2 << "\n 1: " << r3 << "\n"; */
  /*cerr << "matrix initial, before e, IMAG :\n 1: " << i0 << "\n 1: " << i1 << "\n 1: " << i2 << "\n 1: " << i3 << "\n"; */

  VCPLXMUL(r1,i1,e[0],e[1]);
  VCPLXMUL(r2,i2,e[2],e[3]);
  VCPLXMUL(r3,i3,e[4],e[5]);

  /*cerr << "matrix initial, real part:\n 1: " << r0 << "\n 1: " << r1 << "\n 1: " << r2 << "\n 1: " << r3 << "\n"; */
  /*cerr << "matrix initial, imag part:\n 1: " << i0 << "\n 1: " << i1 << "\n 1: " << i2 << "\n 1: " << i3 << "\n"; */

  sr0 = VADD(r0,r2); dr0 = VSUB(r0,r2);
  sr1 = VADD(r1,r3); dr1 = VSUB(r3,r1);
  si0 = VADD(i0,i2); di0 = VSUB(i0,i2);
  si1 = VADD(i1,i3); di1 = VSUB(i3,i1);

  r0 = VADD(sr0, sr1);
  r3 = VSUB(sr0, sr1);
  i0 = VADD(si0, si1);
  i3 = VSUB(si1, si0);
  r1 = VADD(dr0, di1);
  r2 = VSUB(dr0, di1);
  i1 = VSUB(dr1, di0);
  i2 = VADD(dr1, di0);

  *out++ = r0;
  *out++ = i0;
  *out++ = r1;
  *out++ = i1;
  *out++ = r2;
  *out++ = i2;
  *out++ = r3;
  *out++ = i3;

}

static ALWAYS_INLINE(void) pffft_real_preprocess_4x4(const v4sf *in,
                                             const v4sf *e, v4sf *out, int first) {
  v4sf r0=in[0], i0=in[1], r1=in[2], i1=in[3], r2=in[4], i2=in[5], r3=in[6], i3=in[7];
  /*
    transformation for each column is:

    [1   1   1   1   0   0   0   0]   [r0]
    [1   0   0  -1   0  -1  -1   0]   [r1]
    [1  -1  -1   1   0   0   0   0]   [r2]
    [1   0   0  -1   0   1   1   0]   [r3]
    [0   0   0   0   1  -1   1  -1] * [i0]
    [0  -1   1   0   1   0   0   1]   [i1]
    [0   0   0   0   1   1  -1  -1]   [i2]
    [0   1  -1   0   1   0   0   1]   [i3]
  */

  v4sf sr0 = VADD(r0,r3), dr0 = VSUB(r0,r3);
  v4sf sr1 = VADD(r1,r2), dr1 = VSUB(r1,r2);
  v4sf si0 = VADD(i0,i3), di0 = VSUB(i0,i3);
  v4sf si1 = VADD(i1,i2), di1 = VSUB(i1,i2);

  r0 = VADD(sr0, sr1);
  r2 = VSUB(sr0, sr1);
  r1 = VSUB(dr0, si1);
  r3 = VADD(dr0, si1);
  i0 = VSUB(di0, di1);
  i2 = VADD(di0, di1);
  i1 = VSUB(si0, dr1);
  i3 = VADD(si0, dr1);

  VCPLXMULCONJ(r1,i1,e[0],e[1]);
  VCPLXMULCONJ(r2,i2,e[2],e[3]);
  VCPLXMULCONJ(r3,i3,e[4],e[5]);

  VTRANSPOSE4(r0,r1,r2,r3);
  VTRANSPOSE4(i0,i1,i2,i3);

  if (!first) {
    *out++ = r0;
    *out++ = i0;
  }
  *out++ = r1;
  *out++ = i1;
  *out++ = r2;
  *out++ = i2;
  *out++ = r3;
  *out++ = i3;
}

#else

/* standard routine using scalar floats, without SIMD stuff. */

#define pffft_zreorder_nosimd pffft_zreorder
static
void pffft_zreorder_nosimd(PFFFT_Setup *setup, const float *in, float *out, pffft_direction_t direction) {
  int k, N = setup->N;
  if (setup->transform == PFFFT_COMPLEX) {
    for (k=0; k < 2*N; ++k) out[k] = in[k];
    return;
  }
  else if (direction == PFFFT_FORWARD) {
    float x_N = in[N-1];
    for (k=N-1; k > 1; --k) out[k] = in[k-1];
    out[0] = in[0];
    out[1] = x_N;
  } else {
    float x_N = in[1];
    for (k=1; k < N-1; ++k) out[k] = in[k+1];
    out[0] = in[0];
    out[N-1] = x_N;
  }
}

#define pffft_transform_internal_nosimd pffft_transform_internal

#endif /* defined(PFFFT_SIMD_DISABLE) */

#endif
