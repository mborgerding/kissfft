/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */

/* Regression tests for the fast-convolution FIR filter in tools/kiss_fastfir.c,
 * driven the way the fastconv tool drives it (read chunks into inbuf at the
 * running offset, call kiss_fastfir, then flush with n_new == 0).
 *
 * The output is compared with a direct O(N*M) convolution in double. The
 * filter drops the start-up transient, so output[j] is the full convolution
 * at index j + (n_imp_resp - 1), and an input of N samples yields
 * max(0, N - (n_imp_resp - 1)) outputs.
 *
 *  1. flush with fewer input samples than n_imp_resp - 1 (including none at
 *     all) computes a negative size_t and memcpy()s it
 *  2. kiss_fastfir moves the unused input samples to the front of the buffer
 *     with memcpy() although source and destination partially overlap when
 *     the leftover is longer than the consumed part (e.g. nfft 1024 with
 *     1000 taps)
 *
 * Neither misbehavior is reliably visible in the output (an overlapping
 * memcpy happens to work with most C libraries), so memcpy() is redirected
 * while the library source is compiled and the calls are checked directly.
 *
 * Fast convolution is not available for fixed point (see "TODO" in
 * kiss_fastfir.c) or SIMD builds, so the test only runs for float/double.
 * It exercises the complex (non REAL_FASTFIR) variant.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "kiss_fft.h"
#include "_kiss_fft_guts.h"

#if defined(FIXED_POINT) || defined(USE_SIMD)

int main(void)
{
    printf("(skipped: fast convolution is only available for float and double)\n");
    printf("PASS\n");
    return 0;
}

#else

/* ---- memcpy checking ------------------------------------------------------ */

static int n_overlap = 0; /* memcpy calls with partially overlapping source and destination */
static int n_bogus = 0;   /* memcpy calls with an impossible (negative) size; not performed */

static void *checked_memcpy(void *dst, const void *src, size_t n)
{
    uintptr_t d = (uintptr_t)dst, s = (uintptr_t)src;

    if (n > SIZE_MAX / 2) {
        ++n_bogus;
        return dst; /* do not crash the test */
    }
    if (d != s && n > 0 && d < s + n && s < d + n) /* identical buffers: harmless no-op */
        ++n_overlap;
    return memmove(dst, src, n); /* correct result either way */
}

#define memcpy(dst, src, n) checked_memcpy(dst, src, n)

/* ---- library under test, compiled in with memcpy redirected -------------- */

/* the library itself is not built with -W; keep its sign-compare noise out of this test */
#if defined(__GNUC__) || defined(__clang__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wsign-compare"
#endif
#include "../tools/kiss_fastfir.c"
#if defined(__GNUC__) || defined(__clang__)
# pragma GCC diagnostic pop
#endif

#undef memcpy

/* ---- test ---------------------------------------------------------------- */

#define FIR_TOL (sizeof(kiss_fft_scalar) == 4 ? 1e-3 : 1e-9) /* relative to the output rms */

static int failures = 0;

static unsigned long long rng_state = 88172645463325252ULL;
static double rng(void) /* uniform in [-1, 1) */
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (double)(rng_state >> 11) / 4503599627370496.0 - 1.0;
}

/* Filter N random samples with an nh-tap random filter in nfft-point blocks
 * (0: let the library choose) and compare with direct convolution. */
static void check_case(const char *label, size_t nh, size_t nfft, size_t N)
{
    size_t nlag = nh - 1;
    size_t want_out = N > nlag ? N - nlag : 0;
    size_t i, j, k, ngood, nbuf, idx = 0, pos = 0, nout = 0, nread;
    int bad_count = 0;
    double rms = 0, err = 0;
    kffsamp_t *h = (kffsamp_t *)malloc(sizeof(kffsamp_t) * nh);
    kffsamp_t *x = (kffsamp_t *)malloc(sizeof(kffsamp_t) * (N + 1));
    kffsamp_t *out = (kffsamp_t *)malloc(sizeof(kffsamp_t) * (N + 1));
    kffsamp_t *inbuf, *outbuf;
    kiss_fastfir_cfg cfg;

    for (i = 0; i < nh; ++i) {
        h[i].r = (kiss_fft_scalar)rng();
        h[i].i = (kiss_fft_scalar)rng();
    }
    for (i = 0; i < N; ++i) {
        x[i].r = (kiss_fft_scalar)rng();
        x[i].i = (kiss_fft_scalar)rng();
    }

    cfg = kiss_fastfir_alloc(h, nh, &nfft, NULL, NULL);
    ngood = nfft - nh + 1;
    nbuf = nfft + 4 * ngood; /* the buffer size the fastconv tool uses */
    inbuf = (kffsamp_t *)malloc(sizeof(kffsamp_t) * nbuf);
    outbuf = (kffsamp_t *)malloc(sizeof(kffsamp_t) * nbuf);
    if (!h || !x || !out || !cfg || !inbuf || !outbuf) {
        printf("FAIL: %s: allocation failed\n", label);
        ++failures;
        goto done;
    }

    n_overlap = n_bogus = 0;
    do {
        size_t nw;
        nread = N - pos < nbuf - idx ? N - pos : nbuf - idx;
        memmove(inbuf + idx, x + pos, sizeof(kffsamp_t) * nread);
        pos += nread;
        nw = kiss_fastfir(cfg, inbuf, outbuf, nread, &idx);
        if (nw > nbuf || nout + nw > want_out) {
            bad_count = 1; /* negative count wrapped around, or more than the input can produce */
            break;
        }
        memmove(out + nout, outbuf, sizeof(kffsamp_t) * nw);
        nout += nw;
    } while (nread);

    if (n_bogus) {
        printf("FAIL: %s: memcpy of a negative size (nfft %lu, %lu taps, %lu samples)\n",
               label, (unsigned long)nfft, (unsigned long)nh, (unsigned long)N);
        ++failures;
    }
    if (n_overlap) {
        printf("FAIL: %s: memcpy with partially overlapping source and destination, %d call(s) (nfft %lu, %lu taps, %lu samples)\n",
               label, n_overlap, (unsigned long)nfft, (unsigned long)nh, (unsigned long)N);
        ++failures;
    }
    if (bad_count) {
        printf("FAIL: %s: kiss_fastfir returned an impossible output count (expected %lu samples in total)\n",
               label, (unsigned long)want_out);
        ++failures;
        goto done;
    }
    if (nout != want_out) {
        printf("FAIL: %s: produced %lu output samples, expected %lu\n", label,
               (unsigned long)nout, (unsigned long)want_out);
        ++failures;
        goto done;
    }

    /* direct convolution: y[t] = sum_i h[i] * x[t-i], output[j] = y[j + nlag] */
    for (j = 0; j < want_out; ++j) {
        double yr = 0, yi = 0;
        for (k = 0; k < nh; ++k) {
            const kffsamp_t *xs = &x[j + nlag - k];
            yr += (double)h[k].r * xs->r - (double)h[k].i * xs->i;
            yi += (double)h[k].r * xs->i + (double)h[k].i * xs->r;
        }
        rms += yr * yr + yi * yi;
        yr -= out[j].r;
        yi -= out[j].i;
        if (yr * yr + yi * yi > err)
            err = yr * yr + yi * yi;
    }
    if (want_out) {
        rms = sqrt(rms / want_out);
        if (sqrt(err) > FIR_TOL * (rms > 1 ? rms : 1)) {
            printf("FAIL: %s: output differs from direct convolution, max error %g (rms %g)\n",
                   label, sqrt(err), rms);
            ++failures;
        }
    }
done:
    free(h);
    free(x);
    free(out);
    free(cfg);
    free(inbuf);
    free(outbuf);
}

int main(void)
{
    (void)verbose;

    /* inputs shorter than the filter transient produce no output (and must not crash) */
    check_case("empty input, 100 taps", 100, 0, 0);
    check_case("10 samples, 100 taps", 100, 0, 10);
    check_case("99 samples (= taps-1), 100 taps", 100, 0, 99);
    check_case("100 samples, 100 taps", 100, 0, 100);

    /* ordinary streaming */
    check_case("5000 samples, 100 taps", 100, 0, 5000);
    check_case("3000 samples, 37 taps, nfft 256", 37, 256, 3000);

    /* little input consumed per call, so the leftover moves by an overlapping copy */
    check_case("20000 samples, 1000 taps, nfft 1024", 1000, 1024, 20000);

    printf("%s\n", failures ? "FAIL" : "PASS");
    return failures ? 1 : 0;
}

#endif
