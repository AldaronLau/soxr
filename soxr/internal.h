/* SoX Resampler Library      Copyright (c) 2007-16 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details. */

#if !defined soxr_internal_included
#define soxr_internal_included

#include <limits.h>
#include <stdbool.h>
#include <stdint.h>

#undef min
#undef max
#define min(a, b) ((a) <= (b) ? (a) : (b))
#define max(a, b) ((a) >= (b) ? (a) : (b))



#define range_limit(x, lower, upper) (min(max(x, lower), upper))
#define linear_to_dB(x) (log10(x) * 20)
#define array_length(a) (sizeof(a)/sizeof(a[0]))
#if !defined AL
#define AL(a) array_length(a)
#endif
#define iAL(a) (int)AL(a)
#define sqr(a) ((a) * (a))

#define UNUSED __attribute__ ((unused))

#endif
