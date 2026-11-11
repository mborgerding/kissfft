/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */

/* Regression tests for full-scale input to kiss_fftr / kiss_fftri.
 *
 * In fixed-point builds the real-FFT wrappers must not overflow on inputs that
 * the complex kiss_fft handles correctly (both scale by 1/N):
 *  1. kiss_fftr packs the reals into an N/2-point complex FFT, which can see
 *     samples up to sqrt(2)*full-scale
 *  2. kiss_fftri adds DC and Nyquist before the fixed-point divide by two
 *
 * The oracle is the complex kiss_fft on the same data (imag = 0, or the
 * Hermitian-extended spectrum), so the same test runs in float/double too.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kiss_fft.h"
#include "kiss_fftr.h"

#ifndef USE_SIMD

#ifdef FIXED_POINT
# if (FIXED_POINT == 32)
#  define FULLSCALE 2147483647.0
# else
#  define FULLSCALE 32767.0
# endif
# define TOL 16.0 /* LSB */
#else
# define FULLSCALE 1.0
# define TOL (sizeof(kiss_fft_scalar) == 4 ? 1e-4 : 1e-10) /* relative to unity input */
#endif

static int failures = 0;

static void fail(const char *msg, int a, int b)
{
    printf("FAIL: %s (%d,%d)\n", msg, a, b);
    ++failures;
}

/* tolerance for bins computed from n samples: fixed point is already 1/n scaled */
static double tol_for(int n)
{
#ifdef FIXED_POINT
    (void)n;
    return TOL;
#else
    return TOL * n;
#endif
}

/* kiss_fftr of real input x[0..n-1] must match kiss_fft of x + 0i in bins 0..n/2;
 * returns the number of mismatching bins */
static int fftr_matches_complex(const kiss_fft_scalar *x, int n, int *firstbin)
{
    int i, k, bad = 0;
    kiss_fft_cfg c = kiss_fft_alloc(n, 0, NULL, NULL);
    kiss_fftr_cfg r = kiss_fftr_alloc(n, 0, NULL, NULL);
    kiss_fft_cpx *in = malloc(sizeof(*in) * n);
    kiss_fft_cpx *want = malloc(sizeof(*want) * n);
    kiss_fft_cpx *got = malloc(sizeof(*got) * (n / 2 + 1));

    *firstbin = -1;
    if (!c || !r || !in || !want || !got) {
        fail("allocation failed in fftr_matches_complex", n, 0);
        bad = 1;
    } else {
        for (i = 0; i < n; ++i) {
            in[i].r = x[i];
            in[i].i = 0;
        }
        kiss_fft(c, in, want);
        kiss_fftr(r, x, got);
        for (k = 0; k <= n / 2; ++k) {
            if (fabs((double)got[k].r - (double)want[k].r) > tol_for(n) ||
                fabs((double)got[k].i - (double)want[k].i) > tol_for(n)) {
                if (*firstbin < 0)
                    *firstbin = k;
                ++bad;
            }
        }
    }
    free(in);
    free(want);
    free(got);
    kiss_fft_free(c);
    kiss_fftr_free(r);
    return bad;
}

/* kiss_fftri of the half spectrum f[0..n/2] must match the real part of the
 * inverse complex kiss_fft of its Hermitian extension; returns mismatching samples */
static int fftri_matches_complex(const kiss_fft_cpx *f, int n, int *firstidx)
{
    int i, bad = 0;
    kiss_fft_cfg c = kiss_fft_alloc(n, 1, NULL, NULL);
    kiss_fftr_cfg r = kiss_fftr_alloc(n, 1, NULL, NULL);
    kiss_fft_cpx *full = malloc(sizeof(*full) * n);
    kiss_fft_cpx *want = malloc(sizeof(*want) * n);
    kiss_fft_scalar *got = malloc(sizeof(*got) * n);

    *firstidx = -1;
    if (!c || !r || !full || !want || !got) {
        fail("allocation failed in fftri_matches_complex", n, 0);
        bad = 1;
    } else {
        for (i = 0; i <= n / 2; ++i)
            full[i] = f[i];
        for (i = n / 2 + 1; i < n; ++i) {
            full[i].r = f[n - i].r;
            full[i].i = -f[n - i].i;
        }
        kiss_fft(c, full, want);
        kiss_fftri(r, f, got);
        for (i = 0; i < n; ++i) {
            if (fabs((double)got[i] - (double)want[i].r) > tol_for(n)) {
                if (*firstidx < 0)
                    *firstidx = i;
                ++bad;
            }
        }
    }
    free(full);
    free(want);
    free(got);
    kiss_fft_free(c);
    kiss_fftr_free(r);
    return bad;
}

static unsigned long long rng_state = 88172645463325252ULL;
static unsigned rng(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (unsigned)(rng_state >> 33);
}

/* fixed-point repro: this input keeps the scaled spectrum in range, but the
 * packed half-length FFT sees samples of magnitude sqrt(2)*full-scale */
static void test_fftr_known_input(void)
{
    static const int pat[16] = {1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1, 0, -1, 1, -1};
    kiss_fft_scalar x[16];
    int i, first, bad;

    for (i = 0; i < 16; ++i)
        x[i] = (kiss_fft_scalar)(pat[i] * FULLSCALE);
    bad = fftr_matches_complex(x, 16, &first);
    if (bad) {
        printf("(%d of 9 bins differ from complex kiss_fft)\n", bad);
        fail("kiss_fftr n=16 full-scale input differs from kiss_fft, first bin", first, 0);
    }
}

static void test_fftr_random_fullscale(void)
{
    static const int sizes[] = {4, 6, 8, 10, 16, 30, 32, 64, 100, 128, 256};
    int s, t, trials = 0, bad = 0, firstn = 0, firsttrial = 0;

    for (s = 0; s < (int)(sizeof(sizes) / sizeof(sizes[0])); ++s) {
        for (t = 0; t < 20; ++t) {
            int n = sizes[s], i, first;
            kiss_fft_scalar *x = malloc(sizeof(*x) * n);
            if (x == NULL) {
                fail("malloc failed", n, 0);
                continue;
            }
            for (i = 0; i < n; ++i)
                x[i] = (kiss_fft_scalar)((rng() & 1) ? FULLSCALE : -FULLSCALE);
            ++trials;
            if (fftr_matches_complex(x, n, &first)) {
                if (!bad) {
                    firstn = n;
                    firsttrial = t;
                }
                ++bad;
            }
            free(x);
        }
    }
    if (bad) {
        printf("(%d of %d random +-full-scale inputs differ from complex kiss_fft)\n", bad, trials);
        fail("kiss_fftr random full-scale input, first n/trial", firstn, firsttrial);
    }
}

/* DC and Nyquist bins only: x[n] = (X0 + (-1)^n * Xnyq) / N after scaling.
 * Their sum exceeds full-scale even though the result is tiny. */
static void test_fftri_dc_nyquist(void)
{
    static const int sizes[] = {8, 30};
    static const double dc[] = {0.75, -0.75, 0.75, -0.75};
    static const double ny[] = {0.5, -0.5, -0.5, 0.5};
    int s, v, bad = 0, firstn = 0, firstv = 0, cases = 0;

    for (s = 0; s < (int)(sizeof(sizes) / sizeof(sizes[0])); ++s) {
        for (v = 0; v < 4; ++v) {
            int n = sizes[s], i, first;
            kiss_fft_cpx *f = malloc(sizeof(*f) * (n / 2 + 1));
            if (f == NULL) {
                fail("malloc failed", n, 0);
                continue;
            }
            for (i = 0; i <= n / 2; ++i) {
                f[i].r = 0;
                f[i].i = 0;
            }
            f[0].r = (kiss_fft_scalar)(dc[v] * FULLSCALE);
            f[n / 2].r = (kiss_fft_scalar)(ny[v] * FULLSCALE);
            ++cases;
            if (fftri_matches_complex(f, n, &first)) {
                if (!bad) {
                    firstn = n;
                    firstv = v;
                }
                ++bad;
            }
            free(f);
        }
    }
    if (bad) {
        printf("(%d of %d DC+Nyquist spectra differ from complex kiss_fft)\n", bad, cases);
        fail("kiss_fftri DC+Nyquist spectrum, first n/variant", firstn, firstv);
    }
}

int main(void)
{
    test_fftr_known_input();
    test_fftr_random_fullscale();
    test_fftri_dc_nyquist();
    printf("%s\n", failures ? "FAIL" : "PASS");
    return failures ? 1 : 0;
}

#else /* USE_SIMD */

int main(void)
{
    printf("(skipped: scalar arithmetic not available under USE_SIMD)\n");
    printf("PASS\n");
    return 0;
}

#endif
