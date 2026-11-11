/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */

/* Regression test for kissfft<>::assign(nfft, inverse) with an unchanged nfft:
 * the documented postcondition is that the object is in the same state as a
 * freshly constructed kissfft(nfft, inverse). A direction flip must therefore
 * give the same transform as a new object, including for sizes with a radix-4
 * stage, whose butterfly reads the inverse flag.
 */
#include "kissfft.hh"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

typedef std::complex<double> cpx;

static int failures = 0;

static std::vector<cpx> run(const kissfft<double> &fft, const std::vector<cpx> &in)
{
    std::vector<cpx> out(in.size());
    fft.transform(&in[0], &out[0]);
    return out;
}

/* max |a[i]-b[i]| */
static double max_diff(const std::vector<cpx> &a, const std::vector<cpx> &b)
{
    double e = 0;
    for (std::size_t i = 0; i < a.size(); ++i)
        e = std::max(e, std::abs(a[i] - b[i]));
    return e;
}

static void check(const char *what, std::size_t n, const kissfft<double> &got, const kissfft<double> &want)
{
    std::vector<cpx> in(n);
    for (std::size_t i = 0; i < n; ++i)
        in[i] = cpx(std::cos(1.3 * i) + 0.25 * i, std::sin(0.7 * i * i) - 0.5);

    const double e = max_diff(run(got, in), run(want, in));
    if (!(e < 1e-9)) { // fresh objects agree to ~1e-15
        std::printf("FAIL: %s n=%lu differs from a freshly constructed object (max diff %g)\n",
                    what, (unsigned long)n, e);
        ++failures;
    }
}

int main()
{
    static const std::size_t sizes[] = {4, 8, 16, 12, 6, 5};

    for (std::size_t s = 0; s < sizeof(sizes) / sizeof(sizes[0]); ++s) {
        const std::size_t n = sizes[s];

        kissfft<double> flipped(n, false);
        flipped.assign(n, true);
        kissfft<double> fresh_inv(n, true);
        check("assign(n,true) on forward object", n, flipped, fresh_inv);

        // and the other way round, on an independent object
        kissfft<double> flipped_back(n, true);
        flipped_back.assign(n, false);
        kissfft<double> fresh_fwd(n, false);
        check("assign(n,false) on inverse object", n, flipped_back, fresh_fwd);
    }

    std::printf("%s\n", failures ? "FAIL" : "PASS");
    return failures ? 1 : 0;
}
