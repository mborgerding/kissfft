#ifndef KISSFFT_I32_CLASS_HH
#define KISSFFT_I32_CLASS_HH

#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

// TODO1: substitute complex<type> (behaviour not defined for nonfloats), should be faster
// TODO3: make unittests for all ffts (c, cpp, i32)

template <typename DType>
struct complex_s
{
    DType real;
    DType imag;
};

class kissfft_i32
{
private:

    using scalar_type = std::int32_t;
    using cpx_type    = std::complex<std::int32_t>;

    scalar_type _scale_factor;
    std::size_t _nfft;
    bool _inverse;
    std::vector<cpx_type> _twiddles;
    std::vector<std::size_t> _stageRadix;
    std::vector<std::size_t> _stageRemainder;

public:

    // scale_factor: upscale twiddle-factors otherwise they lie between 0..1 (out of range for integer) --> fixed point math
    kissfft_i32(const std::size_t nfft, const bool inverse, const double scale_factor = 1024.0)
            : _scale_factor(scalar_type(scale_factor)), _nfft(nfft), _inverse(inverse)
    {
        // fill twiddle factors
        _twiddles.resize(_nfft);
        const double phinc = (_inverse ? 2 : -2) * std::acos(-1.0) / _nfft;
        for (std::size_t i = 0; i < _nfft; ++i)
        {
            const std::complex<double> tw = scale_factor * std::exp(std::complex<double>(0, i * phinc));
            _twiddles[i] = cpx_type(scalar_type(std::lround(tw.real())), scalar_type(std::lround(tw.imag())));
        }
        //factorize
        //start factoring out 4's, then 2's, then 3,5,7,9,...
        std::size_t n = _nfft;
        std::size_t p = 4;
        do
        {
            while (n % p)
            {
                switch (p)
                {
                    case 4:
                        p = 2;
                        break;
                    case 2:
                        p = 3;
                        break;
                    default:
                        p += 2;
                        break;
                }
                if (p * p > n)  p = n;// no more factors
            }
            n /= p;
            _stageRadix.push_back(p);
            _stageRemainder.push_back(n);
        } while (n > 1);
    }

    /// Calculates the complex Discrete Fourier Transform.
    ///
    /// The size of the passed arrays must be passed in the constructor.
    /// The sum of the squares of the absolute values in the @c dst
    /// array will be @c N times the sum of the squares of the absolute
    /// values in the @c src array, where @c N is the size of the array.
    /// In other words, the l_2 norm of the resulting array will be
    /// @c sqrt(N) times as big as the l_2 norm of the input array.
    /// This is also the case when the inverse flag is set in the
    /// constructor. Hence when applying the same transform twice, but with
    /// the inverse flag changed the second time, then the result will
    /// be equal to the original input times @c N.
    void transform(const cpx_type * FSrc,
                   cpx_type * FDst,
                   const std::size_t stage = 0,
                   const std::size_t fstride = 1,
                   const std::size_t in_stride = 1) const
    {
        const std::size_t p = _stageRadix[stage];
        const std::size_t m = _stageRemainder[stage];
        cpx_type *const Fout_beg = FDst;
        cpx_type *const Fout_end = FDst + p * m;

        if (m == 1)
        {
            do
            {
                *FDst = *FSrc;
                FSrc += fstride * in_stride;
            } while (++FDst != Fout_end);
        }
        else
        {
            do
            {
                // recursive call:
                // DFT of size m*p performed by doing
                // p instances of smaller DFTs of size m,
                // each one takes a decimated version of the input
                transform(FSrc, FDst, stage + 1, fstride * p, in_stride);
                FSrc += fstride * in_stride;
            } while ((FDst += m) != Fout_end);
        }

        FDst = Fout_beg;

        // recombine the p smaller DFTs
        switch (p)
        {
            case 2:
                kf_bfly2(FDst, fstride, m);
                break;
            case 3:
                kf_bfly3(FDst, fstride, m);
                break;
            case 4:
                kf_bfly4(FDst, fstride, m);
                break;
            case 5:
                kf_bfly5(FDst, fstride, m);
                break;
            default:
                kf_bfly_generic(FDst, fstride, m, p);
                break;
        }
    }

private:

    // (a*b + c*d) / _scale_factor, computed in 64 bits: the product of a sample
    // and an upscaled twiddle factor easily overflows int32_t
    scalar_type mac(const scalar_type a, const scalar_type b, const scalar_type c, const scalar_type d) const
    {
        return scalar_type((std::int64_t(a) * b + std::int64_t(c) * d) / _scale_factor);
    }

    // (x * tw) / _scale_factor
    cpx_type mul_tw(const cpx_type &x, const cpx_type &tw) const
    {
        return cpx_type(mac(x.real(), tw.real(), x.imag(), -tw.imag()),
                        mac(x.real(), tw.imag(), x.imag(), tw.real()));
    }

    void kf_bfly2(cpx_type *const Fout, const std::size_t fstride, const std::size_t m) const
    {
        for (std::size_t k = 0; k < m; ++k)
        {
            const cpx_type t = mul_tw(Fout[m + k], _twiddles[k * fstride]);
            Fout[m + k] = Fout[k] - t;
            Fout[k] += t;
        }
    }

    void kf_bfly3(cpx_type *Fout, const std::size_t fstride, const std::size_t m) const
    {
        std::size_t k = m;
        const std::size_t m2 = 2 * m;
        const cpx_type *tw1, *tw2;
        cpx_type scratch[5];
        const cpx_type epi3 = _twiddles[fstride * m];

        tw1 = tw2 = &_twiddles[0];

        do
        {
            scratch[1] = mul_tw(Fout[m], *tw1);
            scratch[2] = mul_tw(Fout[m2], *tw2);

            scratch[3] = scratch[1] + scratch[2];
            scratch[0] = scratch[1] - scratch[2];
            tw1 += fstride;
            tw2 += fstride * 2;

            Fout[m] = Fout[0] - (scratch[3] / 2);
            scratch[0] = mul_tw(scratch[0], cpx_type(epi3.imag(), 0));

            Fout[0] += scratch[3];

            Fout[m2] = cpx_type(Fout[m].real() + scratch[0].imag(), Fout[m].imag() - scratch[0].real());

            Fout[m] += cpx_type(-scratch[0].imag(), scratch[0].real());
            ++Fout;
        } while (--k);
    }

    void kf_bfly4(cpx_type *const Fout, const std::size_t fstride, const std::size_t m) const
    {
        cpx_type scratch[7];
        const scalar_type negative_if_inverse = _inverse ? -1 : +1;

        for (std::size_t k = 0; k < m; ++k)
        {
            scratch[0] = mul_tw(Fout[k + m], _twiddles[k * fstride]);
            scratch[1] = mul_tw(Fout[k + 2 * m], _twiddles[k * fstride * 2]);
            scratch[2] = mul_tw(Fout[k + 3 * m], _twiddles[k * fstride * 3]);
            scratch[5] = Fout[k] - scratch[1];

            Fout[k] += scratch[1];
            scratch[3] = scratch[0] + scratch[2];
            scratch[4] = scratch[0] - scratch[2];
            scratch[4] = cpx_type(scratch[4].imag() * negative_if_inverse,
                                  -scratch[4].real() * negative_if_inverse);

            Fout[k + 2 * m] = Fout[k] - scratch[3];
            Fout[k] += scratch[3];
            Fout[k + m] = scratch[5] + scratch[4];
            Fout[k + 3 * m] = scratch[5] - scratch[4];
        }
    }

    void kf_bfly5(cpx_type *const Fout, const std::size_t fstride, const std::size_t m) const
    {
        cpx_type *Fout0, *Fout1, *Fout2, *Fout3, *Fout4;
        cpx_type scratch[13];
        const cpx_type ya = _twiddles[fstride * m];
        const cpx_type yb = _twiddles[fstride * 2 * m];

        Fout0 = Fout;
        Fout1 = Fout0 + m;
        Fout2 = Fout0 + 2 * m;
        Fout3 = Fout0 + 3 * m;
        Fout4 = Fout0 + 4 * m;

        for (std::size_t u = 0; u < m; ++u)
        {
            scratch[0] = *Fout0;

            scratch[1] = mul_tw(*Fout1, _twiddles[u * fstride]);
            scratch[2] = mul_tw(*Fout2, _twiddles[2 * u * fstride]);
            scratch[3] = mul_tw(*Fout3, _twiddles[3 * u * fstride]);
            scratch[4] = mul_tw(*Fout4, _twiddles[4 * u * fstride]);

            scratch[7] = scratch[1] + scratch[4];
            scratch[10] = scratch[1] - scratch[4];
            scratch[8] = scratch[2] + scratch[3];
            scratch[9] = scratch[2] - scratch[3];

            *Fout0 += scratch[7];
            *Fout0 += scratch[8];

            scratch[5] = scratch[0] + cpx_type(
                    mac(scratch[7].real(), ya.real(), scratch[8].real(), yb.real()),
                    mac(scratch[7].imag(), ya.real(), scratch[8].imag(), yb.real()));

            scratch[6] = cpx_type(
                    mac(scratch[10].imag(), ya.imag(), scratch[9].imag(), yb.imag()),
                    -mac(scratch[10].real(), ya.imag(), scratch[9].real(), yb.imag()));

            *Fout1 = scratch[5] - scratch[6];
            *Fout4 = scratch[5] + scratch[6];

            scratch[11] = scratch[0] + cpx_type(
                    mac(scratch[7].real(), yb.real(), scratch[8].real(), ya.real()),
                    mac(scratch[7].imag(), yb.real(), scratch[8].imag(), ya.real()));

            scratch[12] = cpx_type(
                    mac(scratch[9].imag(), ya.imag(), scratch[10].imag(), -yb.imag()),
                    mac(scratch[10].real(), yb.imag(), scratch[9].real(), -ya.imag()));

            *Fout2 = scratch[11] + scratch[12];
            *Fout3 = scratch[11] - scratch[12];

            ++Fout0;
            ++Fout1;
            ++Fout2;
            ++Fout3;
            ++Fout4;
        }
    }

    /* perform the butterfly for one stage of a mixed radix FFT */
    void kf_bfly_generic(cpx_type * const Fout, const std::size_t fstride, const std::size_t m, const std::size_t p) const
    {
        const cpx_type *twiddles = &_twiddles[0];
        cpx_type scratchbuf[p];

        for (std::size_t u = 0; u < m; ++u)
        {
            std::size_t k = u;
            for (std::size_t q1 = 0; q1 < p; ++q1)
            {
                scratchbuf[q1] = Fout[k];
                k += m;
            }

            k = u;
            for (std::size_t q1 = 0; q1 < p; ++q1)
            {
                std::size_t twidx = 0;
                Fout[k] = scratchbuf[0];
                for (std::size_t q = 1; q < p; ++q)
                {
                    twidx += fstride * k;
                    if (twidx >= _nfft)
                        twidx -= _nfft;
                    Fout[k] += mul_tw(scratchbuf[q], twiddles[twidx]);
                }
                k += m;
            }
        }
    }
};

#endif
