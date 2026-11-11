/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */

/* Regression tests for kissfft_i32.hh:
 *  1. the header must be self-contained: it is deliberately the first include
 *     below, before <complex>, <cstdint> and <cmath> and without any
 *     "using namespace std"
 *  2. the twiddle-factor multiplication must not overflow int32_t for
 *     inputs of a few million (4000000 * 1024 > 2^31)
 *  3. the upscaled twiddle factors must be rounded, not truncated
 *     (truncation biases every output bin)
 *
 * Results are compared with a double precision DFT.
 */
#include "kissfft_i32.hh"

#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

typedef std::complex<std::int32_t> cpx;
typedef std::complex<double> dcpx;

static int failures = 0;

/* direct DFT of the integer input */
static std::vector<dcpx> dft(const std::vector<cpx> &in)
{
    const std::size_t n = in.size();
    const double pi = std::acos(-1.0);
    std::vector<dcpx> out(n);
    for (std::size_t k = 0; k < n; ++k) {
        dcpx acc = 0;
        for (std::size_t j = 0; j < n; ++j)
            acc += dcpx(in[j].real(), in[j].imag()) * std::polar(1.0, -2 * pi * double((k * j) % n) / double(n));
        out[k] = acc;
    }
    return out;
}

/* transform "in" with kissfft_i32 and require every bin within tol of the exact DFT */
static void check(const char *what, const std::vector<cpx> &in, double tol)
{
    const std::size_t n = in.size();
    kissfft_i32 fft(n, false);
    std::vector<cpx> out(n);
    const std::vector<dcpx> want = dft(in);
    std::size_t worst = 0;
    double worst_err = 0;

    fft.transform(&in[0], &out[0]);
    for (std::size_t k = 0; k < n; ++k) {
        const double err = std::abs(dcpx(out[k].real(), out[k].imag()) - want[k]);
        if (err > worst_err) {
            worst_err = err;
            worst = k;
        }
    }
    if (!(worst_err <= tol)) {
        std::printf("FAIL: %s: bin %lu is (%ld,%ld), expected (%.1f,%.1f) (error %g, tolerance %g)\n", what,
                    (unsigned long)worst, (long)out[worst].real(), (long)out[worst].imag(),
                    want[worst].real(), want[worst].imag(), worst_err, tol);
        ++failures;
    }
}

int main()
{
    /* 2. impulse of 4000000 at index 1: the exact result is 4000000 * exp(-i*pi*k/2),
     * all twiddles are exactly representable, so only overflow can spoil it */
    {
        std::vector<cpx> in(4, cpx(0, 0));
        in[1] = cpx(4000000, 0);
        check("n=4 impulse 4000000 (twiddle product overflows int32)", in, 2.0);
    }

    /* 3. ramp 1000..6000, n=6: rounded twiddles give an error below 1 (observed 0.2),
     * truncated ones are off by up to 7 (observed: bin 1 = (-3003,5190), exact (-3000,5196.2)) */
    {
        std::vector<cpx> in(6);
        for (std::size_t j = 0; j < in.size(); ++j)
            in[j] = cpx(std::int32_t(1000 * (j + 1)), 0);
        check("n=6 ramp 1000..6000 (twiddles truncated instead of rounded)", in, 3.0);
    }

    std::printf("%s\n", failures ? "FAIL" : "PASS");
    return failures ? 1 : 0;
}
