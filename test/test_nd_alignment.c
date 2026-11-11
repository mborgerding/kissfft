/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */

/* Regression test for the memory layout of kiss_fftnd_alloc / kiss_fftndr_alloc.
 *
 * Both pack several objects into one block. Every object must start at an
 * address suitably aligned for its type, otherwise the transforms perform
 * misaligned accesses (undefined behavior; a SIGBUS/UBSan error on strict
 * platforms, invisible on x86 without sanitizers):
 *  1. kiss_fftnd: with an odd number of dimensions the int dims[] array shifts
 *     tmpbuf and the per-dimension kiss_fft_state objects by 4 bytes, which
 *     breaks 8-byte aligned types (double)
 *  2. kiss_fftndr: the nd state (which contains pointers) follows the real-FFT
 *     state, whose size need not be a multiple of the pointer size (int16_t)
 *
 * The private structs are visible because the library sources are compiled
 * into this test; the states are built in a caller-supplied, over-aligned
 * buffer, so the result depends only on the library's own layout code.
 * No transform is run.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdint.h>
#include "kiss_fft.h"
#include "_kiss_fft_guts.h"

/* the library itself is not built with -W; keep its sign-compare noise out of this test */
#if defined(__GNUC__) || defined(__clang__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wsign-compare"
#endif
#include "../kiss_fftr.c"
#include "../kiss_fftnd.c"
#include "../kiss_fftndr.c"
#if defined(__GNUC__) || defined(__clang__)
# pragma GCC diagnostic pop
#endif

/* C89-portable alignment requirement of a type: padding after a char member */
#define ALIGN_PROBE(tag, type) struct tag { char c; type x; }
#define ALIGNOF(tag) offsetof(struct tag, x)

ALIGN_PROBE(probe_cpx, kiss_fft_cpx);
ALIGN_PROBE(probe_fft_state, struct kiss_fft_state);
ALIGN_PROBE(probe_fftr_state, struct kiss_fftr_state);
ALIGN_PROBE(probe_fftnd_state, struct kiss_fftnd_state);
ALIGN_PROBE(probe_ptr, void *);
ALIGN_PROBE(probe_int, int);

#define BUF_ALIGN 64 /* more than any type here needs, so only the layout can misalign */

static int failures = 0;

static void fmt_dims(char *out, const int *dims, int ndims)
{
    int i;
    char *p = out;
    *p++ = '{';
    for (i = 0; i < ndims; ++i)
        p += sprintf(p, i ? ",%d" : "%d", dims[i]);
    *p++ = '}';
    *p = '\0';
}

static void check_aligned(const char *owner, const char *dimstr, const char *what,
                          int index, const void *p, size_t align)
{
    size_t off = (size_t)((uintptr_t)p % align);
    if (off != 0) {
        if (index >= 0)
            printf("FAIL: %s dims %s: %s[%d] at offset %lu is misaligned (needs %lu-byte alignment)\n",
                   owner, dimstr, what, index, (unsigned long)off, (unsigned long)align);
        else
            printf("FAIL: %s dims %s: %s at offset %lu is misaligned (needs %lu-byte alignment)\n",
                   owner, dimstr, what, (unsigned long)off, (unsigned long)align);
        ++failures;
    }
}

static void check_nd_state(const char *owner, const char *dimstr, kiss_fftnd_cfg nd)
{
    int k;
    check_aligned(owner, dimstr, "states array", -1, nd->states, ALIGNOF(probe_ptr));
    check_aligned(owner, dimstr, "dims array", -1, nd->dims, ALIGNOF(probe_int));
    check_aligned(owner, dimstr, "tmpbuf", -1, nd->tmpbuf, ALIGNOF(probe_cpx));
    for (k = 0; k < nd->ndims; ++k)
        check_aligned(owner, dimstr, "kiss_fft_state", k, nd->states[k], ALIGNOF(probe_fft_state));
}

/* allocate through the caller-supplied-buffer path, in a BUF_ALIGN-aligned block;
 * returns NULL (and *raw == NULL or to be freed) if the size is not accepted */
static void *alloc_in_aligned_buffer(size_t len, void **raw)
{
    *raw = malloc(len + BUF_ALIGN);
    if (*raw == NULL)
        return NULL;
    return (void *)(((uintptr_t)*raw + BUF_ALIGN - 1) & ~(uintptr_t)(BUF_ALIGN - 1));
}

static void test_fftnd(const int *dims, int ndims)
{
    char dimstr[64];
    size_t len = 0;
    void *raw, *buf;
    kiss_fftnd_cfg st;

    fmt_dims(dimstr, dims, ndims);
    (void)kiss_fftnd_alloc(dims, ndims, 0, NULL, &len);
    buf = alloc_in_aligned_buffer(len, &raw);
    if (buf == NULL) {
        printf("FAIL: out of memory for kiss_fftnd dims %s\n", dimstr);
        ++failures;
        return;
    }
    st = kiss_fftnd_alloc(dims, ndims, 0, buf, &len);
    if (st == NULL) {
        printf("FAIL: kiss_fftnd_alloc rejected dims %s\n", dimstr);
        ++failures;
    } else {
        check_nd_state("kiss_fftnd", dimstr, st);
    }
    free(raw);
}

static void test_fftndr(const int *dims, int ndims)
{
    char dimstr[64];
    size_t len = 0;
    void *raw, *buf;
    kiss_fftndr_cfg st;

    fmt_dims(dimstr, dims, ndims);
    (void)kiss_fftndr_alloc(dims, ndims, 0, NULL, &len);
    if (len == 0 || len == (size_t)-1) {
        printf("(note: kiss_fftndr_alloc size query rejected dims %s, skipped)\n", dimstr);
        return;
    }
    buf = alloc_in_aligned_buffer(len, &raw);
    if (buf == NULL) {
        printf("FAIL: out of memory for kiss_fftndr dims %s\n", dimstr);
        ++failures;
        return;
    }
    st = kiss_fftndr_alloc(dims, ndims, 0, buf, &len);
    if (st == NULL) {
        printf("(note: kiss_fftndr_alloc rejected dims %s, skipped)\n", dimstr);
    } else {
        check_aligned("kiss_fftndr", dimstr, "real-FFT state", -1, st->cfg_r, ALIGNOF(probe_fftr_state));
        check_aligned("kiss_fftndr", dimstr, "nd state", -1, st->cfg_nd, ALIGNOF(probe_fftnd_state));
        check_aligned("kiss_fftndr", dimstr, "tmpbuf", -1, st->tmpbuf, ALIGNOF(probe_cpx));
        if (ndims > 1)
            check_nd_state("kiss_fftndr's nd state", dimstr, st->cfg_nd);
    }
    free(raw);
}

int main(void)
{
    static const int d7[] = {7};
    static const int d345[] = {3, 4, 5};
    static const int d46[] = {4, 6};
    static const int d23456[] = {2, 3, 4, 5, 6};
    static const int d6[] = {6};
    static const int d346[] = {3, 4, 6};
    static const int d58[] = {5, 8};

    test_fftnd(d7, 1);
    test_fftnd(d345, 3);
    test_fftnd(d46, 2);
    test_fftnd(d23456, 5);

    test_fftndr(d46, 2);
    test_fftndr(d6, 1);
    test_fftndr(d346, 3);
    test_fftndr(d58, 2);
    test_fftndr(d23456, 5);

    printf("%s\n", failures ? "FAIL" : "PASS");
    return failures ? 1 : 0;
}
