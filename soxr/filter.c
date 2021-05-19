/* SoX Resampler Library      Copyright (c) 2007-16 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details. */

#include "filter.h"

#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "fft4g.h"
#include "ccrw2.h"

static int * lsx_fft_br_f;
static float * lsx_fft_sc_f;
static int fft_len_f = -1;

void _soxr_init_fft_cache_f(void) {
  if (fft_len_f >= 0)
    return;
  assert(lsx_fft_br_f == NULL);
  assert(lsx_fft_sc_f == NULL);
  assert(fft_len_f == -1);
  fft_len_f = 0;
}

void _soxr_clear_fft_cache_f(void)
{
  assert(fft_len_f >= 0);
  free(lsx_fft_br_f);
  free(lsx_fft_sc_f);
  lsx_fft_sc_f = NULL;
  lsx_fft_br_f = NULL;
  fft_len_f = -1;
}

// Returns true if writer, false if not.
static bool update_fft_cache_f(int len) {
  _soxr_init_fft_cache_f();
  assert(lsx_is_power_of_2(len));
  assert(fft_len_f >= 0);
  if (len > fft_len_f) {
    if (len > fft_len_f) {
      int old_n = fft_len_f;
      int dft_br_len = 2ul + (1ul << (int)(log(len / 2 + .5) / log(2.)) / 2);
      int dft_sc_len = (unsigned long) len / 2;
      fft_len_f = len;
      
      lsx_fft_br_f = realloc(lsx_fft_br_f, dft_br_len * sizeof(*lsx_fft_br_f));
      lsx_fft_sc_f = realloc(lsx_fft_sc_f, dft_sc_len * sizeof(*lsx_fft_sc_f));
      if (!old_n) {
        lsx_fft_br_f[0] = 0;
        atexit(_soxr_clear_fft_cache_f);
      }
      return true;
    }
  }
  return false;
}

void _soxr_safe_cdft_f(int len, int type, float* d) {
    update_fft_cache_f(len);
    cdft(len, type, d, lsx_fft_br_f, lsx_fft_sc_f);
}

void _soxr_safe_rdft_f(int len, int type, float* d) {
    update_fft_cache_f(len);
    rdft(len, type, d, lsx_fft_br_f, lsx_fft_sc_f);
}

void _soxr_ordered_convolve_f(int n, void * not_used, float * a, const float * b)
{
  int i;
  a[0] *= b[0];
  a[1] *= b[1];
  for (i = 2; i < n; i += 2) {
    float tmp = a[i];
    a[i  ] = b[i  ] * tmp - b[i+1] * a[i+1];
    a[i+1] = b[i+1] * tmp + b[i  ] * a[i+1];
  }
  (void)not_used;
}

void _soxr_ordered_partial_convolve_f(int n, float * a, const float * b)
{
  int i;
  a[0] *= b[0];
  for (i = 2; i < n; i += 2) {
    float tmp = a[i];
    a[i  ] = b[i  ] * tmp - b[i+1] * a[i+1];
    a[i+1] = b[i+1] * tmp + b[i  ] * a[i+1];
  }
  a[1] = b[i] * a[i] - b[i+1] * a[i+1];
}

double _soxr_kaiser_beta(double att, double tr_bw)
{
  if (att >= 60) {
    static const double coefs[][4] = {
      {-6.784957e-10,1.02856e-05,0.1087556,-0.8988365+.001},
      {-6.897885e-10,1.027433e-05,0.10876,-0.8994658+.002},
      {-1.000683e-09,1.030092e-05,0.1087677,-0.9007898+.003},
      {-3.654474e-10,1.040631e-05,0.1087085,-0.8977766+.006},
      {8.106988e-09,6.983091e-06,0.1091387,-0.9172048+.015},
      {9.519571e-09,7.272678e-06,0.1090068,-0.9140768+.025},
      {-5.626821e-09,1.342186e-05,0.1083999,-0.9065452+.05},
      {-9.965946e-08,5.073548e-05,0.1040967,-0.7672778+.085},
      {1.604808e-07,-5.856462e-05,0.1185998,-1.34824+.1},
      {-1.511964e-07,6.363034e-05,0.1064627,-0.9876665+.18},
    };
    double realm = log(tr_bw/.0005)/log(2.);
    double const * c0 = coefs[range_limit(  (int)realm, 0, (int)array_length(coefs)-1)];
    double const * c1 = coefs[range_limit(1+(int)realm, 0, (int)array_length(coefs)-1)];
    double b0 = ((c0[0]*att + c0[1])*att + c0[2])*att + c0[3];
    double b1 = ((c1[0]*att + c1[1])*att + c1[2])*att + c1[3];
    return b0 + (b1 - b0) * (realm - (int)realm);
  }
  if (att > 50   ) return .1102 * (att - 8.7);
  if (att > 20.96) return .58417 * pow(att -20.96, .4) + .07886 * (att - 20.96);
  return 0;
}

double * _soxr_make_lpf(int num_taps, double Fc, double beta, double rho, double scale) {
  int i, m = num_taps - 1;
  double * h = malloc((size_t)num_taps * sizeof(*h));
  double mult = scale / _soxr_bessel_I_0(beta), mult1 = 1 / (.5 * m + rho);
  assert(Fc >= 0 && Fc <= 1);

  if (h) for (i = 0; i <= m / 2; ++i) {
    double z = i - .5 * m, x = z * M_PI, y = z * mult1;
    h[i] = x!=0? sin(Fc * x) / x : Fc;
    h[i] *= _soxr_bessel_I_0(beta * sqrt(1 - y * y)) * mult;
    if (m - i != i)
      h[m - i] = h[i];
  }
  return h;
}

void _soxr_kaiser_params(double att, double Fc, double tr_bw, double * beta, int * num_taps)
{
  *beta = *beta < 0? _soxr_kaiser_beta(att, tr_bw * .5 / Fc): *beta;
  att = att < 60? (att - 7.95) / (2.285 * M_PI * 2) :
    ((.0007528358-1.577737e-05**beta)**beta+.6248022)**beta+.06186902;
  *num_taps = !*num_taps? (int)ceil(att/tr_bw + 1) : *num_taps;
}

double * _soxr_design_lpf(
    double Fp,      /* End of pass-band */
    double Fs,      /* Start of stop-band */
    double Fn,      /* Nyquist freq; e.g. 0.5, 1, PI */
    double att,     /* Stop-band attenuation in dB */
    int * num_taps, /* 0: value will be estimated */
    int k,          /* >0: number of phases; <0: num_taps = 1 (mod -k) */
    double beta)    /* <0: value will be estimated */
{
  int n = *num_taps, phases = max(k, 1), modulo = max(-k, 1);
  double tr_bw, Fc, rho = phases == 1? .5 : att < 120? .63 : .75;

  Fp /= fabs(Fn), Fs /= fabs(Fn);        /* Normalise to Fn = 1 */
  tr_bw = .5 * (Fs - Fp); /* Transition band-width: 6dB to stop points */
  tr_bw /= phases, Fs /= phases;
  tr_bw = min(tr_bw, .5 * Fs);
  Fc = Fs - tr_bw;
  assert(Fc - tr_bw >= 0);
  _soxr_kaiser_params(att, Fc, tr_bw, &beta, num_taps);
  if (!n)
    *num_taps = phases > 1? *num_taps / phases * phases + phases - 1 :
      (*num_taps + modulo - 2) / modulo * modulo + 1;
  return Fn < 0? 0 : _soxr_make_lpf(*num_taps, Fc, beta, rho, (double)phases);
}

static double sinePhi(double x) {
    return ((2.0517e-07*x-1.1303e-04)*x+.023154)*x+.55924;
}

static double sinePsi(double x) {
    return ((9.0667e-08*x-5.6114e-05)*x+.013658)*x+1.0977;
}

static double sinePow(double x) {
    return log(0.5) / log(sin(x * 0.5));
}

double _soxr_f_resp(double t, double a) {
  double x;
  if (t > (a <= 160? .8 : .82)) {
    double a1 = a+15;
    double p = .00035*a+.375;
    double w = 1/(1-.597)*asin(pow((a1-10.6)/a1,1/p));
    double c = 1+asin(pow(1-a/a1,1/p))/w;
    return a1*(pow(sin((c-t)*w),p)-1);
  }
  if (t > .5) {
    x = sinePsi(a), x = pow(sin((1-t) * x), sinePow(x));
  } else {
    x = sinePhi(a), x = 1 - pow(sin(t * x), sinePow(x));
  }
  return linear_to_dB(x);
}

double _soxr_inv_f_resp(double drop, double a)
{
  double x = sinePhi(a), s;
  drop = exp(drop * (M_LN10 * 0.05));
  s = drop > .5 ? 1 - drop : drop;
  x = asin(pow(s, 1/sinePow(x))) / x;
  return drop > .5? x : 1 -x;
}
