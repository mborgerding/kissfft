/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */

/* Regression tests for kiss_fftndr_alloc / kiss_fftndr:
 *  1. valid even NxN sizes must not be rejected by the allocation size check
 *  2. an odd last dimension (unsupported by kiss_fftr) must make alloc fail
 *  3. ndims == 1 must produce the same spectrum as kiss_fftr
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kiss_fftndr.h"
#include "kiss_fftr.h"

static int failures = 0;

static void fail(const char *msg, int a, int b)
{
    printf("FAIL: %s (%d,%d)\n", msg, a, b);
    ++failures;
}

/* allocate through the caller-supplied-buffer path; returns NULL on failure */
static kiss_fftndr_cfg alloc_ndr(const int *dims, int ndims, void **buf)
{
    size_t len = 0;
    *buf = NULL;
    kiss_fftndr_alloc(dims, ndims, 0, NULL, &len);
    if (len == 0 || len == (size_t)-1)
        return NULL;
    *buf = malloc(len);
    if (*buf == NULL)
        return NULL;
    return kiss_fftndr_alloc(dims, ndims, 0, *buf, &len);
}

static void test_valid_sizes_accepted(void)
{
    int n, first = 0, rejected = 0;
    for (n = 2; n <= 512; n += 2) {
        int dims[2];
        void *buf;
        dims[0] = n;
        dims[1] = n;
        if (alloc_ndr(dims, 2, &buf) == NULL) {
            if (!first)
                first = n;
            ++rejected;
        }
        free(buf);
    }
    if (rejected) {
        printf("(%d valid even NxN sizes in 2..512 rejected)\n", rejected);
        fail("kiss_fftndr_alloc rejected a valid square size, first", first, first);
    }
}

static void test_odd_last_dim_rejected(void)
{
    int dims[2] = {4, 5};
    void *buf;
    if (alloc_ndr(dims, 2, &buf) != NULL)
        fail("kiss_fftndr_alloc accepted an odd last dimension", dims[0], dims[1]);
    free(buf);
}

/* the inner kiss_fftnd counts dims[0..ndims-2] in an int: a product above INT_MAX
 * must be reported as too large by the size query (no allocation is attempted) */
static void test_int_overflow_rejected(void)
{
    int dims[3] = {46341, 46341, 2}; /* 46341^2 > INT_MAX */
    size_t len = 0;
    kiss_fftndr_alloc(dims, 3, 0, NULL, &len);
    if (len != (size_t)-1)
        fail("kiss_fftndr_alloc size query accepted dims whose product overflows int", dims[0], dims[1]);
}

#ifndef USE_SIMD
static void test_one_dimension_matches_fftr(void)
{
    enum { N = 8 };
    int dims[1] = {N};
    kiss_fft_scalar in[N];
    kiss_fft_cpx got[N / 2 + 1], want[N / 2 + 1];
    void *buf;
    kiss_fftndr_cfg ndr;
    kiss_fftr_cfg r;
    int k;

    for (k = 0; k < N; ++k)
        in[k] = (kiss_fft_scalar)(100 * (k + 1));

    ndr = alloc_ndr(dims, 1, &buf);
    r = kiss_fftr_alloc(N, 0, NULL, NULL);
    if (ndr == NULL || r == NULL) {
        fail("allocation failed for ndims == 1", N, 0);
    } else {
        kiss_fftndr(ndr, in, got);
        kiss_fftr(r, in, want);
        for (k = 0; k <= N / 2; ++k)
            if (got[k].r != want[k].r || got[k].i != want[k].i)
                fail("ndims == 1 output differs from kiss_fftr at bin", k, 0);
    }
    free(buf);
    kiss_fftr_free(r);
}
#endif

int main(void)
{
    test_valid_sizes_accepted();
    test_odd_last_dim_rejected();
    test_int_overflow_rejected();
#ifndef USE_SIMD
    test_one_dimension_matches_fftr();
#endif
    printf("%s\n", failures ? "FAIL" : "PASS");
    return failures ? 1 : 0;
}
