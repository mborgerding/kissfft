/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 *
 *  Regression test for github issue #121: kiss_fft_next_fast_size(0)
 *  used to hang forever in the factorization loop.
 */
#include "kiss_fft.h"
#include <stdio.h>

int main(void)
{
    int n;
    /* n<=0 must return 1 per maintainer decision, not hang */
    if (kiss_fft_next_fast_size(0) != 1) {
        fprintf(stderr, "FAIL: kiss_fft_next_fast_size(0) != 1\n");
        return 1;
    }
    if (kiss_fft_next_fast_size(-3) != 1) {
        fprintf(stderr, "FAIL: kiss_fft_next_fast_size(-3) != 1\n");
        return 1;
    }
    /* positive inputs keep their pre-fix behavior */
    {
        int expected[9] = { 0, 1, 2, 3, 4, 5, 6, 8, 8 };
        for (n = 1; n <= 8; ++n) {
            if (kiss_fft_next_fast_size(n) != expected[n]) {
                fprintf(stderr, "FAIL: kiss_fft_next_fast_size(%d) = %d, expected %d\n",
                        n, kiss_fft_next_fast_size(n), expected[n]);
                return 1;
            }
        }
    }
    printf("kiss_fft_next_fast_size: all checks passed\n");
    return 0;
}
