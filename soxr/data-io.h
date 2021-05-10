/* SoX Resampler Library      Copyright (c) 2007-13 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details. */

#if !defined soxr_data_io_included
#define soxr_data_io_included

#include "soxr.h"

void _soxr_deinterleave(
    double * * dest,
    void const * * src0,
    size_t n,
    unsigned ch);

void _soxr_deinterleave_f(
    float * * dest,
    void const * * src0,
    size_t n,
    unsigned ch);

size_t /* clips */ _soxr_interleave(
    void * * dest,
    double const * const * src,
    size_t n,
    unsigned ch,
    unsigned long * seed);

size_t /* clips */ _soxr_interleave_f(
    void * * dest,
    float const * const * src,
    size_t n,
    unsigned ch,
    unsigned long * seed);

#endif
