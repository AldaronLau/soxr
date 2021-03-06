/* SoX Resampler Library      Copyright (c) 2007-16 robs@users.sourceforge.net
 * Licence for this file: LGPL v2.1                  See LICENCE for details.
 *
 * Constant-rate resampling common code. */

#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "filter.h" // _soxr_safe_rdft_f, etc.

#include "internal.h"

#include "cr.h"

static void rdft_forward(int length, float * H) {
    _soxr_safe_rdft_f(length, 1, H);
}

static void rdft_backward(int length, float * H) {
    _soxr_safe_rdft_f(length, -1, H);
}

static float * prepare_poly_fir_coefs(double const * coefs, int num_coefs,
    int num_phases, int interp_order,
    core_flags_t core_flags)
{
    int i, j, length = num_coefs * num_phases * (interp_order + 1);
    float * result = calloc(1, (size_t) length << 2);
    double fm1 = coefs[0], f1 = 0, f2 = 0;

    for (i = num_coefs - 1; i >= 0; --i) {
        for (j = num_phases - 1; j >= 0; --j) {
            double f0 = fm1, b = 0, c = 0, d = 0; /* = 0 to kill compiler warning */
            int pos = i * num_phases + j - 1;
            fm1 = pos > 0 ? coefs[pos - 1] : 0;
            switch (interp_order) {
                case 1: b = f1 - f0; break;
                case 2: b = f1 - (.5 * (f2+f0) - f1) - f0; c = .5 * (f2+f0) - f1; break;
                case 3: c=.5*(f1+fm1)-f0;d=(1/6.)*(f2-f1+fm1-f0-4*c);b=f1-f0-d-c; break;
                default: assert(!interp_order);
            }
            if ((core_flags & 3) == 0) {
                int fir_coef_num = num_coefs - 1 - i;

                if (interp_order > 2) coef(result, interp_order, num_coefs, j, 3, fir_coef_num) = (float)d;
                if (interp_order > 1) coef(result, interp_order, num_coefs, j, 2, fir_coef_num) = (float)c;
                if (interp_order > 0) coef(result, interp_order, num_coefs, j, 1, fir_coef_num) = (float)b;
                coef(result, interp_order, num_coefs, j, 0, fir_coef_num) = (float)f0;
            }
            f2 = f1;
            f1 = f0;
        }
    }
    return result;
}

static void dft_stage_fn(stage_t * p, fifo_t * output_fifo) {
    float* output;
    float* dft_out;
    int i;
    int j;
    int num_in = max(0, fifo_occupancy(&p->fifo));
    rate_shared_t const * s = p->shared;
    dft_filter_t const * f = &s->dft_filter[p->dft_filter_num];
    int const overlap = f->num_taps - 1;

    if ((p->at.ms.all >> 32) + p->L * num_in >= f->dft_length) {
        size_t const sizeof_real = sizeof(char) << 2;

        div_t divd = div(f->dft_length - overlap - (p->at.ms.all >> 32) + p->L - 1, p->L);
        float const * input = fifo_read_ptr(&p->fifo);
        fifo_read(&p->fifo, divd.quot, NULL);
        num_in -= divd.quot;

        output = fifo_reserve(output_fifo, f->dft_length);
        dft_out = output;

        if (lsx_is_power_of_2(p->L)) { /* F-domain */
            int portion = f->dft_length / p->L;
            memcpy(dft_out, input, (unsigned)portion * sizeof_real);
            rdft_forward(portion, dft_out);

            for (i = portion + 2; i < (portion << 1); i += 2) /* Mirror image. */
                dft_out[i] = dft_out[(portion << 1) - i],
                    dft_out[i+1] = -dft_out[(portion << 1) - i + 1];
            dft_out[portion] = dft_out[1];
            dft_out[portion + 1] = 0;
            dft_out[1] = dft_out[0];

            for (portion <<= 1; i < f->dft_length; i += portion, portion <<= 1) {
                memcpy((char *)dft_out + (size_t)i * sizeof_real, dft_out, (size_t)portion * sizeof_real);

                dft_out[i + 1] = 0;
            }
        } else {
            if (p->L == 1) {
                memcpy(dft_out, input, (size_t)f->dft_length * sizeof_real);
            } else {
                memset(dft_out, 0, (size_t)f->dft_length * sizeof_real);
                for (j = 0, i = (p->at.ms.all >> 32); i < f->dft_length; ++j, i += p->L) {
                    ((float *)dft_out)[i] = ((float *)input)[j];
                }
          
                p->at.ms.all = (p->at.ms.all & 0x00000000FFFFFFFF)
                    | (uint64_t)(p->L - 1 - divd.rem) << 32;
            }
            if ((p->step.ms.all >> 32) > 0)
                rdft_forward(f->dft_length, dft_out);
            else
                rdft_forward(f->dft_length, dft_out);
        }

        if ((p->step.ms.all >> 32) > 0) {
            _soxr_ordered_convolve_f(f->dft_length, NULL, dft_out, f->coefs);
            rdft_backward(f->dft_length, dft_out);
            if (0 && (p->step.ms.all >> 32) == 1)
                memcpy(output, dft_out, (size_t)f->dft_length * sizeof_real);
            if ((p->step.ms.all >> 32) != 1) {
                for (j = 0, i = p->remM; i < f->dft_length - overlap; ++j,
                      i += (p->step.ms.all >> 32))
                {
                    ((float *)output)[j] = ((float *)dft_out)[i];
                }
                p->remM = i - (f->dft_length - overlap);
                fifo_trim_by(output_fifo, f->dft_length - j);
            } else {
                fifo_trim_by(output_fifo, overlap);
            }
        } else { /* F-domain */
            int m = -(p->step.ms.all >> 32);
            _soxr_ordered_partial_convolve_f(f->dft_length >> m, dft_out, f->coefs);
            rdft_backward(f->dft_length >> m, dft_out);
            fifo_trim_by(output_fifo, (((1 << m) - 1) * f->dft_length + overlap) >>m);
        }
    }
    p->input_size = (f->dft_length - (p->at.ms.all >> 32) + p->L - 1) / p->L;
}

/* Set to 4 x nearest power of 2 or half of that */
/* if danger of causing too many cache misses. */
static int set_dft_length(int num_taps, int min, int large) {
    double d = log((double)num_taps) / log(2.);
    return 1 << range_limit((int)(d + 2.77), min, max((int)(d + 1.77), large));
}

static void dft_stage_init(
    unsigned instance, double Fp, double Fs, double Fn, double att,
    stage_t * p, int L, int M,
    unsigned min_dft_size, unsigned large_dft_size)
{
    dft_filter_t * f = &p->shared->dft_filter[instance];
    int num_taps = 0, dft_length = f->dft_length, i, offset;
    bool f_domain_m = abs(3-M) == 1 && Fs <= 1;
    size_t const sizeof_real = sizeof(char) << 2;

    if (!dft_length) {
        int k = lsx_is_power_of_2(L) && Fn == L? L << 1 : 4;
        double m, * h = _soxr_design_lpf(Fp, Fs, Fn, att, &num_taps, -k, -1.);

        f->post_peak = num_taps / 2;

        dft_length = set_dft_length(num_taps, (int)min_dft_size, (int)large_dft_size);
        f->coefs = calloc((size_t)dft_length, sizeof_real);
        offset = dft_length - num_taps + 1;
        m = (1. / dft_length) * 2 * L;
        for (i = 0; i < num_taps; ++i)
            ((float *)f->coefs)[(i + offset) & (dft_length - 1)] =(float)(h[i] * m);
        free(h);
    }

    if (!f->dft_length) {
        int Mp = f_domain_m? M : 1;
        if (Mp == 1) {
            rdft_forward(dft_length, f->coefs);
        } else {
            rdft_forward(dft_length, f->coefs);
        }
        f->num_taps = num_taps;
        f->dft_length = dft_length;
    }
    p->out_in_ratio = (double)L / M;
    p->core_flags = 0;
    p->fn = dft_stage_fn;
    p->preload = f->post_peak / L;
    p->at.ms.all = (p->at.ms.all & 0x00000000FFFFFFFF)
        | (uint64_t)(f->post_peak % L) << 32;
    p->L = L;
    p->step.ms.all = (p->at.ms.all & 0x00000000FFFFFFFF)
        | (uint64_t)(f_domain_m? -M/2 : M) << 32;
    p->dft_filter_num = instance;
    p->block_len = f->dft_length - (f->num_taps - 1);
    p->phase0 = (p->at.ms.all >> 32) / p->L;
    p->input_size = (f->dft_length - (p->at.ms.all >> 32) + p->L - 1) / p->L;
}

static struct half_fir_info const * find_half_fir(
    struct half_fir_info const * firs, size_t len, double att)
{
    size_t i;
    for (i = 0; i + 1 < len && att > firs[i].att; ++i);
    return &firs[i];
}

#include "stdio.h"

char const * resampler_init(
  resampler_t * const p,             /* Per audio channel. */
  double const io_ratio,        /* Input rate divided by output rate. */
  cr_core_t const * const core)
{
    double const Fp0 = 0.913743653267622;
    bool const upsample = io_ratio < 1.0;

    double arbM = io_ratio;
    double Fp1 = Fp0;
    double Fs1 = 1.0;
    double att = 21.0 * linear_to_dB(2.0);
    double attArb = att; /* +1: pass+stop */
    int arbL = 1;
    bool rational = false;
    stage_t * s;

    p->core = core;
    p->io_ratio = io_ratio;

    /* Determine stages: */
    int try;
    int L;
    int M;
    int x;
    double d;

    int shr = 0;
    int i = (int)(0.5 * io_ratio);
    for (;;) {
        i >>= 1;
        if (i == 0) {
            break;
        }
        arbM *= .5;
        shr += 1;
    }

    int preM = upsample || (arbM > 1.5 && arbM < 2);
    int postM = 1 + (arbM > 1 && preM);
    arbM /= postM;
    int preL = 1 + (!preM && arbM < 2) + (upsample);
    arbM *= preL;
    
    double frac = arbM - (int)arbM;

    double const epsilon
        = fabs(floor(frac * 4294967296.0 + .5) / (frac * 4294967296.0) - 1);
    rational = (frac == 0);
    for (i = 1; i <= 2048; i += 1) {
        d = frac * i;
        try = (int)(d + .5);
        if ((rational = fabs(try / d - 1) <= epsilon)) {    // No long doubles!
            if (try == i) {
                arbM = ceil(arbM);
                x = arbM > 3;
                shr += x;
                arbM /= 1 + x;
            } else {
                arbM = i * (int)arbM + try;
                arbL = i;
            }
            break;
        }
    }
    L = preL * arbL;
    M = (int)(arbM * postM);
    x = (L | M) & 1;
    L >>= !x;
    M >>= !x;
    if ((d = preL * arbL / arbM) > 4 && d != 5) {
        int postL = 4;
        i = (int)(d / 32);
        while (i != 0 && postL < 256) {
            postL <<= 1;
            i >>= 1;
        }
        arbM = arbM * postL / arbL / preL;
        arbL = 1;
    } else if (rational && (max(L, M) < 3 + 2 || L * M < 6)) {
        preL = L;
        preM = M;
        arbM = 1;
        arbL = 1;
        postM = 1;
    }
    // End determine stages.

    p->num_stages = shr + (preM * preL != 1) + (arbM * arbL != 1);

    printf("THIS = %d; %d %d %f %d\n", p->num_stages, preM, preL, arbM, arbL);

    p->stages = calloc((size_t)p->num_stages + 1, sizeof(*p->stages));

    for (i = 0; i < p->num_stages; ++i) {
        p->stages[i].num = i;
        p->stages[i].shared = &p->shared;
        p->stages[i].input_size = 8192;
    }
    p->stages[0].is_input = true;

    /* Att. budget: */
    if ((p->num_stages) > 1) {
        if ((arbM * arbL  != 1)) {
            att += linear_to_dB(2.);
            attArb = att;
            att += linear_to_dB(p->num_stages - 1.0);
        } else {
            att += linear_to_dB(p->num_stages);
        }
    }

    struct half_fir_info const * half_fir_info = find_half_fir(
        core->half_firs,
        core->half_firs_len,
        att
    );

    for (i = 0, s = p->stages; i < shr; ++i, ++s) {
        s->fn = NULL;
        s->coefs = half_fir_info->coefs;
        s->n = half_fir_info->num_coefs;
        s->pre_post = 4 * s->n;
        s->pre = s->pre_post >> 1;
        s->preload = s->pre;
    }

    if ((preM * preL != 1)) {
        printf("UYOYEFOEJFEFIJEFIKJFEIEJFKEJAEKFNMA\n");

        double Fn1 = preM ? max(preL, preM) : arbM / arbL;
        dft_stage_init(
            0,
            Fp0,
            Fs1,
            Fn1,
            att,
            s,
            preL,
            max(preM, 1),
            10,
            17
        );
        s += 1;
        Fp1 /= Fn1;
        Fs1 /= Fn1;
    }

    /* Higher quality arb stage: */
    printf("YYYYYYYYO!\n");
    poly_fir_t const * f = &core->poly_firs[6 * (upsample + !!preM) + 3 - !upsample];
    int order;
    int num_coefs = (int)f->interp.scalar;
    int phase_bits;
    int phases;
    size_t coefs_size;
    double at;
    double Fs;
    double Fn;
    double mult = upsample ? 1.0 : arbM / arbL;
    poly_fir1_t const * f1;

    if (!upsample && preM) {
        Fn = 2 * mult;
        Fs = 3 + fabs(Fs1 - 1);
    } else {
        Fn = 1;
        Fs = 2 - (Fp1 + (Fs1 - Fp1) * 0.7);
    }

    Fp1 = Fs - (Fs - Fp1) / (1 - _soxr_inv_f_resp(-.01f, attArb));

    i = !rational - 1;
    for(;;) {
        i += 1;
        f1 = &f->interp;
        assert(f1->fn);
        if (i != 0) {
            arbM /= arbL;
            arbL = 1;
            rational = false;
        }
        phase_bits = (int)ceil(f1->scalar - log(mult)/log(2.));
        phases = !rational? (1 << phase_bits) : arbL;
        if (f->interp.scalar==0) {
            int phases0 = max(phases, 19), n0 = 0;
            _soxr_design_lpf(Fp1, Fs, -Fn, attArb, &n0, phases0, f->beta);
            num_coefs = n0 / phases0 + 1, num_coefs += num_coefs & !preM;
        }
        if ((num_coefs & 1) != 0 && rational && (arbL & 1) != 0) {
            phases <<= 1;
            arbL <<= 1;
            arbM *= 2;
        }
        at = arbL * (s->phase0 = 0.5 * (num_coefs & 1));
        order = i;
        coefs_size = (size_t)(num_coefs * phases * (order+1)) * sizeof(float);
      
        if (i >= 2 || !f->interp.fn || coefs_size / 1000 <= 400) {
            break;
        }
    };

    if (!s->shared->poly_fir_coefs) {
      int num_taps = num_coefs * phases - 1;
      double * coefs = _soxr_design_lpf(
          Fp1, Fs, Fn, attArb, &num_taps, phases, f->beta);
      s->shared->poly_fir_coefs = prepare_poly_fir_coefs(
          coefs, num_coefs, phases, order, 0);
      free(coefs);
    }
    s->fn = f1->fn;
    s->pre_post = num_coefs - 1;
    s->preload = ((num_coefs - 1) >> 1) + (num_coefs - num_coefs);
    s->n = num_coefs;
    s->phase_bits = phase_bits;
    s->L = arbL;

    s->at.ms.all = (int64_t)(at * 4294967296.0 + 0.5);
    s->step.ms.all = (int64_t)(arbM * 4294967296.0 + 0.5);
    s->out_in_ratio = 4294967296.0 * arbL / (double)s->step.ms.all;
    // End higher-quality arb stage

    s = p->stages;
    for (i = 0; i < p->num_stages; i += 1) {
        printf(" KI %d\n", i);
  
        fifo_create(&s->fifo, (int)sizeof(float));
        memset(fifo_reserve(&s->fifo, s->preload), 0,
            sizeof(float) * (size_t)s->preload);

        s += 1;
    }
    fifo_create(&s->fifo, (int)sizeof(float));
    return 0;
}

static bool stage_process(stage_t * stage, bool flushing) {
  fifo_t * fifo = &stage->fifo;
  bool done = false;
  int want;
  while (!done && (want = stage->input_size - fifo_occupancy(fifo)) > 0) {
    if (stage->is_input) {
      if (flushing)
        memset(fifo_reserve(fifo, want), 0, fifo->item_size * (size_t)want);
      else done = true;
    }
    else done = stage_process(stage - 1, flushing);
  }
  stage->fn(stage, &stage[1].fifo);
  return done
    && fifo_occupancy(fifo) < stage->input_size;
}

void resampler_process(resampler_t * p, size_t olen) {
    int const n = p->flushing? min(-(int)p->samples_out, (int)olen) : (int)olen;
    stage_t * stage = &p->stages[p->num_stages];
    fifo_t * fifo = &stage->fifo;
    bool done = false;
    while (!done && fifo_occupancy(fifo) < (int)n) {
        done = stage->is_input || stage_process(stage - 1, p->flushing);
    }
}

float* resampler_input(resampler_t * p, float const * samples, size_t n) {
    p->samples_in += (int64_t)n;
    return fifo_write(&p->stages[0].fifo, (int)n, samples);
}

float const* resampler_output(resampler_t * p, float * samples, size_t * n0) {
    fifo_t * fifo = &p->stages[p->num_stages].fifo;
    int n = p->flushing? min(-(int)p->samples_out, (int)*n0) : (int)*n0;
    p->samples_out += n = min(n, fifo_occupancy(fifo));
    return fifo_read(fifo, (int)(*n0 = (size_t)n), samples);
}

void resampler_flush(resampler_t * p) {
    p->samples_out -= (int64_t)((double)p->samples_in / p->io_ratio + .5);
    p->samples_in = 0;
    p->flushing = true;
}
