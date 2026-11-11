/*
 *  Copyright (c) 2003-2004, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */

#include "kiss_fftr.h"
#include "_kiss_fft_guts.h"

struct kiss_fftr_state{
    kiss_fft_cfg substate;
    kiss_fft_cpx * tmpbuf;
    kiss_fft_cpx * super_twiddles;
#ifdef USE_SIMD
    void * pad;
#endif
};

kiss_fftr_cfg kiss_fftr_alloc(int nfft,int inverse_fft,void * mem,size_t * lenmem)
{
	KISS_FFT_ALIGN_CHECK(mem)

    int i;
    kiss_fftr_cfg st = NULL;
    size_t subsize = 0, memneeded;

    if (nfft & 1) {
        KISS_FFT_ERROR("Real FFT optimization must be even.");
        return NULL;
    }
    nfft >>= 1;

    kiss_fft_alloc (nfft, inverse_fft, NULL, &subsize);
    memneeded = sizeof(struct kiss_fftr_state) + subsize + sizeof(kiss_fft_cpx) * ( nfft * 3 / 2);

    if (lenmem == NULL) {
        st = (kiss_fftr_cfg) KISS_FFT_MALLOC (memneeded);
    } else {
        if (*lenmem >= memneeded)
            st = (kiss_fftr_cfg) mem;
        *lenmem = memneeded;
    }
    if (!st)
        return NULL;

    st->substate = (kiss_fft_cfg) (st + 1); /*just beyond kiss_fftr_state struct */
    st->tmpbuf = (kiss_fft_cpx *) (((char *) st->substate) + subsize);
    st->super_twiddles = st->tmpbuf + nfft;
    kiss_fft_alloc(nfft, inverse_fft, st->substate, &subsize);

    for (i = 0; i < nfft/2; ++i) {
        double phase =
            -3.14159265358979323846264338327 * ((double) (i+1) / nfft + .5);
        if (inverse_fft)
            phase *= -1;
        kf_cexp (st->super_twiddles+i,phase);
    }
    return st;
}

/* (a+b)/2 and (a-b)/2: in fixed point the sum can exceed the integer range,
 * so form it in the wider product type before halving */
#ifdef FIXED_POINT
# define HALF_ADD(a,b) ((kiss_fft_scalar)(((SAMPPROD)(a) + (b)) >> 1))
# define HALF_SUB(a,b) ((kiss_fft_scalar)(((SAMPPROD)(a) - (b)) >> 1))
#else
# define HALF_ADD(a,b) HALF_OF((a) + (b))
# define HALF_SUB(a,b) HALF_OF((a) - (b))
#endif

void kiss_fftr(kiss_fftr_cfg st,const kiss_fft_scalar *timedata,kiss_fft_cpx *freqdata)
{
    /* input buffer timedata is stored row-wise */
    int k,ncfft;
    kiss_fft_cpx fpnk,fpk,f1k,f2k,tw,tdc;

    if ( st->substate->inverse) {
        KISS_FFT_ERROR("kiss fft usage error: improper alloc");
        return;/* The caller did not call the correct function */
    }

    ncfft = st->substate->nfft;

    /*perform the parallel fft of two real signals packed in real,imag*/
#ifdef FIXED_POINT
    /* Packed (real,imag) samples can reach sqrt(2)*full-scale, which overflows
     * inside kiss_fft, so halve the input first (the post-pass below would
     * otherwise do this /2). freqdata (ncfft+1 complex) serves as scratch,
     * which also keeps in-place use working. */
    for (k = 0; k < ncfft; ++k) {
        freqdata[k] = ((const kiss_fft_cpx*)timedata)[k];
        C_FIXDIV(freqdata[k],2);
    }
    kiss_fft( st->substate , freqdata, st->tmpbuf );
#else
    kiss_fft( st->substate , (const kiss_fft_cpx*)timedata, st->tmpbuf );
#endif
    /* The real part of the DC element of the frequency spectrum in st->tmpbuf
     * contains the sum of the even-numbered elements of the input time sequence
     * The imag part is the sum of the odd-numbered elements
     *
     * The sum of tdc.r and tdc.i is the sum of the input time sequence.
     *      yielding DC of input time sequence
     * The difference of tdc.r - tdc.i is the sum of the input (dot product) [1,-1,1,-1...
     *      yielding Nyquist bin of input time sequence
     */

    tdc.r = st->tmpbuf[0].r;
    tdc.i = st->tmpbuf[0].i;
    CHECK_OVERFLOW_OP(tdc.r ,+, tdc.i);
    CHECK_OVERFLOW_OP(tdc.r ,-, tdc.i);
    freqdata[0].r = tdc.r + tdc.i;
    freqdata[ncfft].r = tdc.r - tdc.i;
#ifdef USE_SIMD
#ifdef HAVE_LASX
    freqdata[0].i = (__m256)(__lasx_xvreplgr2vr_w(0));
    freqdata[ncfft].i = freqdata[0].i;
#elif defined(HAVE_LSX)
    freqdata[0].i = (__m128)(__lsx_vreplgr2vr_w(0));
    freqdata[ncfft].i = freqdata[0].i;
#else
    freqdata[ncfft].i = freqdata[0].i = _mm_set1_ps(0);
#endif
#else
    freqdata[ncfft].i = freqdata[0].i = 0;
#endif

    for ( k=1;k <= ncfft/2 ; ++k ) {
        fpk    = st->tmpbuf[k];
        fpnk.r =   st->tmpbuf[ncfft-k].r;
        fpnk.i = - st->tmpbuf[ncfft-k].i;

        C_ADD( f1k, fpk , fpnk );
        C_SUB( f2k, fpk , fpnk );
        C_MUL( tw , f2k , st->super_twiddles[k-1]);

        freqdata[k].r = HALF_ADD(f1k.r, tw.r);
        freqdata[k].i = HALF_ADD(f1k.i, tw.i);
        freqdata[ncfft-k].r = HALF_SUB(f1k.r, tw.r);
        freqdata[ncfft-k].i = HALF_SUB(tw.i, f1k.i);
    }
}

void kiss_fftri(kiss_fftr_cfg st,const kiss_fft_cpx *freqdata,kiss_fft_scalar *timedata)
{
    /* input buffer timedata is stored row-wise */
    int k, ncfft;
    kiss_fft_cpx tdc;

    if (st->substate->inverse == 0) {
        KISS_FFT_ERROR("kiss fft usage error: improper alloc");
        return;/* The caller did not call the correct function */
    }

    ncfft = st->substate->nfft;

    /* halve DC and Nyquist before combining them: their sum can exceed full-scale */
    tdc.r = freqdata[0].r;
    tdc.i = freqdata[ncfft].r;
    C_FIXDIV(tdc,2);
    st->tmpbuf[0].r = tdc.r + tdc.i;
    st->tmpbuf[0].i = tdc.r - tdc.i;

    for (k = 1; k <= ncfft / 2; ++k) {
        kiss_fft_cpx fk, fnkc, fek, fok, tmp;
        fk = freqdata[k];
        fnkc.r = freqdata[ncfft - k].r;
        fnkc.i = -freqdata[ncfft - k].i;
        C_FIXDIV( fk , 2 );
        C_FIXDIV( fnkc , 2 );

        C_ADD (fek, fk, fnkc);
        C_SUB (tmp, fk, fnkc);
        C_MUL (fok, tmp, st->super_twiddles[k-1]);
        C_ADD (st->tmpbuf[k],     fek, fok);
        C_SUB (st->tmpbuf[ncfft - k], fek, fok);
#ifdef USE_SIMD
#ifdef HAVE_LASX
        __m256 neg_one = (__m256)__lasx_xvreplgr2vr_w(0xBF800000); // -1.0f
        st->tmpbuf[ncfft - k].i = __lasx_xvfmul_s(st->tmpbuf[ncfft - k].i, neg_one);
#elif defined(HAVE_LSX)
        __m128 neg_one = (__m128)__lsx_vreplgr2vr_w(0xBF800000); // -1.0f
        st->tmpbuf[ncfft - k].i = __lsx_vfmul_s(st->tmpbuf[ncfft - k].i, neg_one);
#else
        st->tmpbuf[ncfft - k].i *= _mm_set1_ps(-1.0);
#endif
#else
        st->tmpbuf[ncfft - k].i *= -1;
#endif
    }
    kiss_fft (st->substate, st->tmpbuf, (kiss_fft_cpx *) timedata);
}
