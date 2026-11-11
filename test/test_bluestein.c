/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */

/* Sizes with a large prime factor (computed with Bluestein's algorithm in floating-point
 * builds) against a direct DFT, plus in-place, strided, caller-memory and real-input calls. */

#include "kiss_fft.h"
#include "kiss_fftr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if defined(FIXED_POINT) || defined(USE_SIMD)
int main(void)
{
    printf("test_bluestein: skipped (floating point only)\n");
    return 0;
}
#else

static int failures = 0;

/* SNR (dB) of kiss output y against the direct DFT of x (n inputs, nout bins) */
static double snr_vs_dft(const kiss_fft_cpx *x, const kiss_fft_cpx *y, int n, int nout, int inverse)
{
    const double pi = 3.14159265358979323846264338327950288;
    double sig = 0, err = 0;
    int k, t;
    for (k = 0; k < nout; ++k) {
        double re = 0, im = 0;
        for (t = 0; t < n; ++t) {
            double ph = (inverse ? 2 : -2) * pi * (double)(((long long)k * t) % n) / n;
            re += x[t].r * cos(ph) - x[t].i * sin(ph);
            im += x[t].r * sin(ph) + x[t].i * cos(ph);
        }
        sig += re * re + im * im;
        err += (re - y[k].r) * (re - y[k].r) + (im - y[k].i) * (im - y[k].i);
    }
    return 10 * log10(sig / err);
}

static void check(int ok, const char *what, int n)
{
    if (!ok) {
        printf("FAIL: %s, nfft=%d\n", what, n);
        ++failures;
    }
}

static void test_complex(int n)
{
    const int stride = 3;
    kiss_fft_cpx *x = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * n * stride);
    kiss_fft_cpx *y = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * n);
    kiss_fft_cpx *z = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * n);
    kiss_fft_cpx *xs = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * n * stride);
    int inverse, i;

    for (i = 0; i < n; ++i) {
        x[i].r = (kiss_fft_scalar)(rand() / (double)RAND_MAX - .5);
        x[i].i = (kiss_fft_scalar)(rand() / (double)RAND_MAX - .5);
        xs[i * stride] = x[i];
    }
    for (inverse = 0; inverse < 2; ++inverse) {
        size_t len = 0;
        void *mem;
        kiss_fft_cfg cfg = kiss_fft_alloc(n, inverse, NULL, NULL), cfg2;
        double snr;

        kiss_fft(cfg, x, y);
        snr = snr_vs_dft(x, y, n, n, inverse);
        printf("nfft=%d inverse=%d snr=%.1f dB\n", n, inverse, snr);
        check(snr > 100, "accuracy", n);

        memcpy(z, x, sizeof(kiss_fft_cpx) * n);
        kiss_fft(cfg, z, z);
        check(memcmp(y, z, sizeof(kiss_fft_cpx) * n) == 0, "in-place differs from out-of-place", n);

        kiss_fft_stride(cfg, xs, z, stride);
        check(memcmp(y, z, sizeof(kiss_fft_cpx) * n) == 0, "strided input differs", n);

        kiss_fft_alloc(n, inverse, NULL, &len);
        mem = KISS_FFT_MALLOC(len);
        cfg2 = kiss_fft_alloc(n, inverse, mem, &len);
        check(cfg2 != NULL, "alloc in caller memory", n);
        if (cfg2) {
            kiss_fft(cfg2, x, z);
            check(memcmp(y, z, sizeof(kiss_fft_cpx) * n) == 0, "caller-memory config differs", n);
        }
        KISS_FFT_FREE(mem);
        kiss_fft_free(cfg);
    }
    free(x); free(y); free(z); free(xs);
}

static void test_real(int n)
{
    kiss_fft_scalar *r = (kiss_fft_scalar*)malloc(sizeof(kiss_fft_scalar) * n);
    kiss_fft_cpx *x = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * n);
    kiss_fft_cpx *y = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * (n / 2 + 1));
    kiss_fftr_cfg cfg = kiss_fftr_alloc(n, 0, NULL, NULL);
    double snr;
    int i;

    for (i = 0; i < n; ++i) {
        r[i] = (kiss_fft_scalar)(rand() / (double)RAND_MAX - .5);
        x[i].r = r[i];
        x[i].i = 0;
    }
    kiss_fftr(cfg, r, y);
    snr = snr_vs_dft(x, y, n, n / 2 + 1, 0);
    printf("real nfft=%d snr=%.1f dB\n", n, snr);
    check(snr > 100, "real accuracy", n);
    kiss_fftr_free(cfg);
    free(r); free(x); free(y);
}

int main(void)
{
    static const int sizes[] = { 23, 31, 155, 1009, 2018, 4099, 6208 };
    unsigned i;

    srand(1);
    for (i = 0; i < sizeof(sizes) / sizeof(sizes[0]); ++i)
        test_complex(sizes[i]);
    test_real(2018);
    test_real(2 * 4099);
    if (failures)
        printf("%d failure(s)\n", failures);
    return failures != 0;
}
#endif
