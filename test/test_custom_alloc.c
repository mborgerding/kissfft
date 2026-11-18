/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */

/* Regression test for user-overridden KISS_FFT_MALLOC / KISS_FFT_FREE.
 *
 * The library sources are compiled INTO this test (unity build) with
 * KISS_FFT_MALLOC/KISS_FFT_FREE replaced by tracking wrappers, so the test
 * does not link against the prebuilt library. In addition, plain malloc()
 * and free() are redirected while the library sources are compiled, to catch
 * code that bypasses the overridable allocator:
 *  1. every allocation made by kiss_fft_alloc, kiss_fftr_alloc,
 *     kiss_fftnd_alloc, kiss_fftndr_alloc and kfc_fft goes through
 *     KISS_FFT_MALLOC (no raw malloc)
 *  2. kiss_fft_free, kiss_fftr_free, kiss_fftndr_free and kfc_cleanup only
 *     release blocks that KISS_FFT_MALLOC returned, and do so through
 *     KISS_FFT_FREE (no raw free)
 *  3. nothing is left allocated afterwards
 *  4. with USE_SIMD, overriding the allocator must still compile (the
 *     KISS_FFT_ALIGN_CHECK / KISS_FFT_ALIGN_SIZE_UP macros must stay defined)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdint.h>
#include <limits.h>

/* the tracking allocator, defined below; must be visible before kiss_fft.h */
static void *test_malloc(size_t n);
static void test_free(void *p);
#define KISS_FFT_MALLOC(nbytes) test_malloc(nbytes)
#define KISS_FFT_FREE(ptr) test_free(ptr)

#include "kiss_fft.h" /* pull in every system header the library needs while  */
#include "_kiss_fft_guts.h" /* malloc/free are still the real functions       */

/* ---- tracking allocator -------------------------------------------------- */

#define ALLOC_ALIGN 32 /* enough for SSE (16) and AVX/LASX (32) kiss_fft_scalar */
#define MAX_BLOCKS 256

struct block {
    void *user;  /* pointer handed to the library */
    void *raw;   /* pointer obtained from the C library */
    int owned;   /* 1: returned by KISS_FFT_MALLOC, 0: library called raw malloc */
};

static struct block blocks[MAX_BLOCKS];
static int nblocks = 0;

/* counters since the last begin_scenario() */
static int n_wrapper_alloc = 0;    /* KISS_FFT_MALLOC calls */
static int n_raw_malloc = 0;       /* library called malloc() directly */
static int n_raw_free_owned = 0;   /* library called free() on a KISS_FFT_MALLOC block */
static int n_foreign_free = 0;     /* KISS_FFT_FREE on a block KISS_FFT_MALLOC did not return */
static int n_unknown_free = 0;     /* free of a pointer nobody allocated */

static int failures = 0;

static void *real_malloc(size_t n)
{
    return malloc(n);
}

static void real_free(void *p)
{
    free(p);
}

static int find_block(const void *p)
{
    int i;
    for (i = 0; i < nblocks; ++i)
        if (blocks[i].user == p)
            return i;
    return -1;
}

static void add_block(void *user, void *raw, int owned)
{
    if (nblocks == MAX_BLOCKS) {
        printf("FAIL: test allocation table full\n");
        exit(1);
    }
    blocks[nblocks].user = user;
    blocks[nblocks].raw = raw;
    blocks[nblocks].owned = owned;
    ++nblocks;
}

static void drop_block(int i)
{
    real_free(blocks[i].raw);
    blocks[i] = blocks[--nblocks];
}

static void *test_malloc(size_t n)
{
    char *raw = (char *)real_malloc(n + ALLOC_ALIGN);
    void *user;
    if (raw == NULL)
        return NULL;
    user = (void *)(((uintptr_t)raw + ALLOC_ALIGN - 1) & ~(uintptr_t)(ALLOC_ALIGN - 1));
    add_block(user, raw, 1);
    ++n_wrapper_alloc;
    return user;
}

static void test_free(void *p)
{
    int i;
    if (p == NULL)
        return;
    i = find_block(p);
    if (i < 0) {
        ++n_unknown_free;
        return; /* not ours: do not hand it to the C library */
    }
    if (!blocks[i].owned)
        ++n_foreign_free;
    drop_block(i);
}

/* what plain malloc()/free() inside the library sources turn into */
static void *raw_malloc(size_t n)
{
    void *p = real_malloc(n);
    if (p != NULL) {
        add_block(p, p, 0);
        ++n_raw_malloc;
    }
    return p;
}

static void raw_free(void *p)
{
    int i;
    if (p == NULL)
        return;
    i = find_block(p);
    if (i < 0) {
        ++n_unknown_free;
        return;
    }
    if (blocks[i].owned)
        ++n_raw_free_owned;
    drop_block(i);
}

#define malloc(nbytes) raw_malloc(nbytes)
#define free(ptr) raw_free(ptr)

/* ---- library under test, compiled in with the overrides above ----------- */

/* the library itself is not built with -W; keep its sign-compare / maybe-uninitialized noise out of this test */
#if defined(__GNUC__) || defined(__clang__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wsign-compare"
# if !defined(__clang__)
#  pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
# endif
#endif
#include "../kiss_fft.c"
#include "../kiss_fftr.c"
#include "../kiss_fftnd.c"
#include "../kiss_fftndr.c"
#include "../kfc.c"
#if defined(__GNUC__) || defined(__clang__)
# pragma GCC diagnostic pop
#endif

/* ---- test ---------------------------------------------------------------- */

static void begin_scenario(void)
{
    n_wrapper_alloc = n_raw_malloc = n_raw_free_owned = 0;
    n_foreign_free = n_unknown_free = 0;
}

/* call right after the *_alloc call under test */
static void expect_wrapper_alloc(const char *what)
{
    if (n_wrapper_alloc == 0) {
        printf("FAIL: %s did not allocate through KISS_FFT_MALLOC\n", what);
        ++failures;
    }
}

/* call after the object was released with the library's own free function */
static void end_scenario(const char *what)
{
    int i, live_owned = 0, live_raw = 0;

    for (i = 0; i < nblocks; ++i) {
        if (blocks[i].owned)
            ++live_owned;
        else
            ++live_raw;
    }
    if (n_raw_malloc) {
        printf("FAIL: %s called malloc() %d time(s) instead of KISS_FFT_MALLOC\n", what, n_raw_malloc);
        ++failures;
    }
    if (n_raw_free_owned) {
        printf("FAIL: %s called free() %d time(s) on KISS_FFT_MALLOC memory instead of KISS_FFT_FREE\n",
               what, n_raw_free_owned);
        ++failures;
    }
    if (n_foreign_free) {
        printf("FAIL: %s passed %d block(s) not allocated by KISS_FFT_MALLOC to KISS_FFT_FREE\n",
               what, n_foreign_free);
        ++failures;
    }
    if (n_unknown_free) {
        printf("FAIL: %s freed %d pointer(s) nobody allocated\n", what, n_unknown_free);
        ++failures;
    }
    if (live_owned || live_raw) {
        printf("FAIL: %s left %d KISS_FFT_MALLOC block(s) and %d malloc block(s) allocated\n",
               what, live_owned, live_raw);
        ++failures;
    }
    while (nblocks > 0) /* do not let one failing scenario leak into the next */
        drop_block(nblocks - 1);
}

/* the redirection itself must work, otherwise a PASS below means nothing */
static void test_harness(void)
{
    void *p;

    begin_scenario();
    p = malloc(8);
    if (n_raw_malloc != 1 || nblocks != 1) {
        printf("FAIL: test harness did not intercept malloc()\n");
        ++failures;
    }
    free(p);
    p = KISS_FFT_MALLOC(8);
    free(p); /* wrong way to release a KISS_FFT_MALLOC block */
    if (n_raw_free_owned != 1 || n_wrapper_alloc != 1 || nblocks != 0) {
        printf("FAIL: test harness did not intercept free()\n");
        ++failures;
    }
}

static void test_kiss_fft(void)
{
    kiss_fft_cpx in[7], out[7];
    kiss_fft_cfg cfg;

    memset(in, 0, sizeof(in));
    begin_scenario();
    cfg = kiss_fft_alloc(7, 0, NULL, NULL); /* prime size: generic butterfly uses a temp buffer */
    expect_wrapper_alloc("kiss_fft_alloc");
    if (cfg != NULL) {
        kiss_fft(cfg, in, out);
        kiss_fft(cfg, in, in); /* "in-place" also uses a temp buffer */
    }
    kiss_fft_free(cfg);
    end_scenario("kiss_fft");
}

static void test_kiss_fftr(void)
{
    kiss_fft_scalar in[16];
    kiss_fft_cpx out[9];
    kiss_fftr_cfg cfg;

    memset(in, 0, sizeof(in));
    begin_scenario();
    cfg = kiss_fftr_alloc(16, 0, NULL, NULL);
    expect_wrapper_alloc("kiss_fftr_alloc");
    if (cfg != NULL)
        kiss_fftr(cfg, in, out);
    kiss_fftr_free(cfg);
    end_scenario("kiss_fftr");
}

static void test_kiss_fftnd(void)
{
    static const int dims[2] = {4, 6};
    kiss_fft_cpx in[24], out[24];
    kiss_fftnd_cfg cfg;

    memset(in, 0, sizeof(in));
    begin_scenario();
    cfg = kiss_fftnd_alloc(dims, 2, 0, NULL, NULL);
    expect_wrapper_alloc("kiss_fftnd_alloc");
    if (cfg != NULL)
        kiss_fftnd(cfg, in, out);
    kiss_fft_free(cfg);
    end_scenario("kiss_fftnd");
}

static void test_kiss_fftndr(void)
{
    static const int dims[2] = {4, 6};
    kiss_fft_scalar in[24];
    kiss_fft_cpx out[4 * 4];
    kiss_fftndr_cfg cfg;

    memset(in, 0, sizeof(in));
    begin_scenario();
    cfg = kiss_fftndr_alloc(dims, 2, 0, NULL, NULL);
    expect_wrapper_alloc("kiss_fftndr_alloc");
    if (cfg != NULL)
        kiss_fftndr(cfg, in, out);
    kiss_fftndr_free(cfg);
    end_scenario("kiss_fftndr");
}

static void test_kfc(void)
{
    kiss_fft_cpx in[8], out[8];

    memset(in, 0, sizeof(in));
    begin_scenario();
    kfc_fft(8, in, out);
    kfc_ifft(8, in, out);
    expect_wrapper_alloc("kfc_fft");
    kfc_cleanup();
    end_scenario("kfc");
}

int main(void)
{
    test_harness();
    test_kiss_fft();
    test_kiss_fftr();
    test_kiss_fftnd();
    test_kiss_fftndr();
    test_kfc();
    printf("%s\n", failures ? "FAIL" : "PASS");
    return failures ? 1 : 0;
}
