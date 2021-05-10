/* SoX Resampler Library      Copyright (c) 2007-16 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details. */

#if !defined PFFT_MACROS_ONLY

#include "math-wrap.h"

#if PFFFT_DOUBLE
  #include "util64s.h"
#else
  #include "util32s.h"
  #define sin(x) sinf(x)
  #define cos(x) cosf(x)
#endif

#define pffft_aligned_free    SIMD_ALIGNED_FREE
#define pffft_aligned_malloc  SIMD_ALIGNED_MALLOC
#define pffft_aligned_calloc  SIMD_ALIGNED_CALLOC

#undef inline
#define inline __inline

#endif

#include "pffft.c"
