/* SoX Resampler Library      Copyright (c) 2007-13 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details. */

/* Concurrent Control with "Readers" and "Writers", P.J. Courtois et al, 1971 */

#if !defined soxr_ccrw2_included
#define soxr_ccrw2_included

#include "internal.h"

typedef int ccrw2_t;
#define ccrw2_become_reader(x) (void)(x)
#define ccrw2_cease_reading(x) (void)(x)
#define ccrw2_become_writer(x) (void)(x)
#define ccrw2_cease_writing(x) (void)(x)
#define ccrw2_init(x) (void)(x)
#define ccrw2_clear(x) (void)(x)

#endif
