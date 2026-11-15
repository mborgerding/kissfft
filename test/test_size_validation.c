/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */

/* Regression tests for size validation in the allocation functions:
 *  1. kiss_fftnd_alloc must not overflow the product of the dimensions (an int):
 *     {65536,65536} wraps to 0 and {46342,46342} to a negative number, so
 *     the library reports a bogus size and accepts a buffer that is far too small
 *  2. kiss_fft_alloc(0) and kiss_fftr_alloc(0) must fail
 *  3. kiss_fft_next_fast_size must not overflow int for n above the largest
 *     representable 2,3,5-smooth number (2125764000)
 *
 * Huge allocations are never attempted: the nd tests only use the size query
 * and a small caller-supplied buffer, and the nfft == 0 objects are never
 * transformed.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>
#include "kiss_fft.h"
#include "kiss_fftr.h"
#include "kiss_fftnd.h"

static int failures = 0;

static void fail(const char *msg, long a, long b)
{
    printf("FAIL: %s (%ld,%ld)\n", msg, a, b);
    ++failures;
}

/* ---- 1. kiss_fftnd_alloc: dimension product overflow ---------------------- */

#define SMALL_BUF_SIZE (16u * 1024 * 1024) /* above the wrapped size, far below the real one */
#define BUF_ALIGN 64

static void test_nd_overflow(int d0, int d1)
{
    int dims[2];
    size_t len = 0, buflen;
    double need; /* the tmpbuf alone needs this many bytes */
    void *raw, *buf;
    kiss_fftnd_cfg st;

    dims[0] = d0;
    dims[1] = d1;
    need = (double)d0 * (double)d1 * (double)sizeof(kiss_fft_cpx);

    /* size query: must report "too big" ((size_t)-1) or a size that really fits the data,
     * not a wrapped-around number */
    (void)kiss_fftnd_alloc(dims, 2, 0, NULL, &len);
    if (len != (size_t)-1 && !((double)len >= need && (double)len <= 2 * need))
        fail("kiss_fftnd_alloc size query reports a bogus size, dims", d0, d1);

    /* a 16MB buffer is far too small for the transform and must be refused */
    raw = malloc(SMALL_BUF_SIZE + BUF_ALIGN);
    if (raw == NULL) {
        fail("out of memory", d0, d1);
        return;
    }
    buf = (void *)(((uintptr_t)raw + BUF_ALIGN - 1) & ~(uintptr_t)(BUF_ALIGN - 1));
    buflen = SMALL_BUF_SIZE;
    st = kiss_fftnd_alloc(dims, 2, 0, buf, &buflen);
    if (st != NULL)
        fail("kiss_fftnd_alloc accepted a 16MB buffer for dims", d0, d1);
    free(raw);
}

/* ---- 2. nfft == 0 --------------------------------------------------------- */

static void test_zero_and_negative_size(void)
{
    kiss_fft_cfg c;
    kiss_fftr_cfg r;

    c = kiss_fft_alloc(0, 0, NULL, NULL);
    if (c != NULL)
        fail("kiss_fft_alloc(0) did not return NULL", 0, 0);
    kiss_fft_free(c);

    r = kiss_fftr_alloc(0, 0, NULL, NULL);
    if (r != NULL)
        fail("kiss_fftr_alloc(0) did not return NULL", 0, 0);
    kiss_fftr_free(r);

    c = kiss_fft_alloc(-1, 0, NULL, NULL);
    if (c != NULL)
        fail("kiss_fft_alloc(-1) did not return NULL", -1, 0);
    kiss_fft_free(c);
}

/* ---- 3. kiss_fft_next_fast_size ------------------------------------------ */

/* n has no prime factors other than 2, 3 and 5 */
static int is_smooth(long long n)
{
    if (n < 1)
        return 0;
    while (n % 2 == 0)
        n /= 2;
    while (n % 3 == 0)
        n /= 3;
    while (n % 5 == 0)
        n /= 5;
    return n == 1;
}

/* smallest smooth number >= n, computed in 64 bits */
static long long ref_next_fast_size(long long n)
{
    while (!is_smooth(n))
        ++n;
    return n;
}

static void test_next_fast_size(void)
{
    static const int spot[] = {1001, 1000001, 2000000000, 2100000001, 2125763999, 2125764000};
    static const int huge[] = {2125764001, 2147483000, INT_MAX};
    int n, i, bad = 0, firstn = 0;

    for (n = 1; n <= 2000; ++n) {
        if (kiss_fft_next_fast_size(n) != ref_next_fast_size(n)) {
            if (!bad)
                firstn = n;
            ++bad;
        }
    }
    if (bad) {
        printf("(%d of 2000 sizes in 1..2000 wrong)\n", bad);
        fail("kiss_fft_next_fast_size differs from reference, first n", firstn, 0);
    }

    for (i = 0; i < (int)(sizeof(spot) / sizeof(spot[0])); ++i) {
        int got = kiss_fft_next_fast_size(spot[i]);
        if (got != ref_next_fast_size(spot[i]))
            fail("kiss_fft_next_fast_size wrong, n", spot[i], got);
    }

    /* the next smooth number does not fit in an int: n itself (still a valid size) is returned */
    for (i = 0; i < (int)(sizeof(huge) / sizeof(huge[0])); ++i) {
        int got = kiss_fft_next_fast_size(huge[i]);
        if (got != huge[i])
            fail("kiss_fft_next_fast_size should return n when no larger fast size fits, n", huge[i], got);
    }
}

int main(void)
{
    test_nd_overflow(65536, 65536);
    test_nd_overflow(46342, 46342);
    test_zero_and_negative_size();
    test_next_fast_size();
    printf("%s\n", failures ? "FAIL" : "PASS");
    return failures ? 1 : 0;
}
