/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */

#include <stdint.h>
#include "_kiss_fft_guts.h"
#include <string.h>

/* Runtime CPU dispatch (x86-64 GCC >= 12 / Clang >= 18 on ELF + glibc, where ifunc is available):
 * the hot routines get an extra AVX2/FMA clone that is selected at load time, and float
 * radix-4 stages use a hand-vectorized AVX2 butterfly when the CPU supports it.
 * Results stay within floating-point rounding of the portable code, but are not bit-identical
 * to it (FMA). Define KISS_FFT_NO_DISPATCH to build the portable code only. */
#if !defined(KISS_FFT_NO_DISPATCH) && !defined(USE_SIMD) && !defined(FIXED_POINT) \
    && defined(__x86_64__) && defined(__ELF__) && defined(__GLIBC__) && !defined(__UCLIBC__) \
    && !defined(__INTEL_COMPILER) \
    && ((defined(__clang__) && __clang_major__ >= 18) || (!defined(__clang__) && defined(__GNUC__) && __GNUC__ >= 12))
#  define KISS_FFT_DISPATCH 1
#  define KF_TARGET_CLONES __attribute__((target_clones("arch=x86-64-v3","default")))
#else
#  define KF_TARGET_CLONES
#endif
/* The guts header contains all the multiplication and addition macros that are defined for
 fixed or floating point complex numbers.  It also delares the kf_ internal functions.
 */

static void kf_bfly2(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        int m
        )
{
    kiss_fft_cpx * Fout2;
    kiss_fft_cpx * tw1 = st->twiddles;
    kiss_fft_cpx t;
    Fout2 = Fout + m;

    /* first iteration peeled: its twiddle is 1 (kept out of the loop so the loop stays branch-free) */
    C_FIXDIV(*Fout,2); C_FIXDIV(*Fout2,2);
    t = *Fout2;
    C_SUB( *Fout2 ,  *Fout , t );
    C_ADDTO( *Fout ,  t );

    while (--m) {
        ++Fout2;
        ++Fout;
        tw1 += fstride;
        C_FIXDIV(*Fout,2); C_FIXDIV(*Fout2,2);

        C_MUL (t,  *Fout2 , *tw1);
        C_SUB( *Fout2 ,  *Fout , t );
        C_ADDTO( *Fout ,  t );
    }
}

static inline void kf_bfly4_dir(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        const size_t m,
        const int inverse
        )
{
    kiss_fft_cpx *tw1,*tw2,*tw3;
    kiss_fft_cpx scratch[6];
    size_t k=m;
    const size_t m2=2*m;
    const size_t m3=3*m;


    tw3 = tw2 = tw1 = st->twiddles;

    do {
        C_FIXDIV(*Fout,4); C_FIXDIV(Fout[m],4); C_FIXDIV(Fout[m2],4); C_FIXDIV(Fout[m3],4);

        if (k != m) {
            C_MUL(scratch[0],Fout[m] , *tw1 );
            C_MUL(scratch[1],Fout[m2] , *tw2 );
            C_MUL(scratch[2],Fout[m3] , *tw3 );
        } else { /* first twiddles are 1 */
            scratch[0] = Fout[m]; scratch[1] = Fout[m2]; scratch[2] = Fout[m3];
        }

        C_SUB( scratch[5] , *Fout, scratch[1] );
        C_ADDTO(*Fout, scratch[1]);
        C_ADD( scratch[3] , scratch[0] , scratch[2] );
        C_SUB( scratch[4] , scratch[0] , scratch[2] );
        C_SUB( Fout[m2], *Fout, scratch[3] );
        tw1 += fstride;
        tw2 += fstride*2;
        tw3 += fstride*3;
        C_ADDTO( *Fout , scratch[3] );

        if(inverse) {
            Fout[m].r = scratch[5].r - scratch[4].i;
            Fout[m].i = scratch[5].i + scratch[4].r;
            Fout[m3].r = scratch[5].r + scratch[4].i;
            Fout[m3].i = scratch[5].i - scratch[4].r;
        }else{
            Fout[m].r = scratch[5].r + scratch[4].i;
            Fout[m].i = scratch[5].i - scratch[4].r;
            Fout[m3].r = scratch[5].r - scratch[4].i;
            Fout[m3].i = scratch[5].i + scratch[4].r;
        }
        ++Fout;
    }while(--k);
}

#ifdef KISS_FFT_DISPATCH
#include <immintrin.h>
#define KISS_FFT_AVX2 1
/* only for kiss_fft_scalar == float: 4 interleaved complex values per __m256 */
#define KF_AVX2_OK (sizeof(kiss_fft_scalar) == sizeof(float))

__attribute__((target("avx2,fma")))
static inline __m256 kf_cmul_avx2(__m256 a, __m256 w)
{
    __m256 wr = _mm256_moveldup_ps(w);          /* wr wr ... */
    __m256 wi = _mm256_movehdup_ps(w);          /* wi wi ... */
    __m256 as = _mm256_permute_ps(a, 0xB1);     /* ai ar ... */
    return _mm256_fmaddsub_ps(a, wr, _mm256_mul_ps(as, wi)); /* (ar*wr-ai*wi, ai*wr+ar*wi) */
}

/* radix-4 butterfly over k = 0..m-1 using contiguous twiddles tw1 = [W^k] and tw23 = [W^2k | W^3k] */
__attribute__((target("avx2,fma")))
static void kf_bfly4_avx2(kiss_fft_cpx * Fout, const kiss_fft_cpx * tw1, const kiss_fft_cpx * tw23,
                          const size_t m, const int inverse)
{
    float *F0 = (float*)Fout, *F1 = (float*)(Fout + m), *F2 = (float*)(Fout + 2*m), *F3 = (float*)(Fout + 3*m);
    const float *T1 = (const float*)tw1, *T2 = (const float*)tw23, *T3 = (const float*)(tw23 + m);
    /* multiplying by -i (forward) / +i (inverse) = swap re/im, then flip the sign bit of the imag / real lanes.
       The flip is an integer xor on purpose: as a float constant (+-0.0) it gets folded away under
       -ffast-math (-fno-signed-zeros), which the default build uses. */
    const __m256i flip = inverse ? _mm256_set_epi32(0, INT_MIN, 0, INT_MIN, 0, INT_MIN, 0, INT_MIN)
                                 : _mm256_set_epi32(INT_MIN, 0, INT_MIN, 0, INT_MIN, 0, INT_MIN, 0);
    const size_t mv = m & ~(size_t)3;
    size_t k;
    for (k = 0; k < mv; k += 4) {
        const size_t o = 2*k;
        __m256 a0 = _mm256_loadu_ps(F0 + o);
        __m256 s0 = kf_cmul_avx2(_mm256_loadu_ps(F1 + o), _mm256_loadu_ps(T1 + o));
        __m256 s1 = kf_cmul_avx2(_mm256_loadu_ps(F2 + o), _mm256_loadu_ps(T2 + o));
        __m256 s2 = kf_cmul_avx2(_mm256_loadu_ps(F3 + o), _mm256_loadu_ps(T3 + o));
        __m256 s5 = _mm256_sub_ps(a0, s1);
        __m256 a  = _mm256_add_ps(a0, s1);
        __m256 s3 = _mm256_add_ps(s0, s2);
        __m256 t  = _mm256_castsi256_ps(_mm256_xor_si256(_mm256_castps_si256(
                        _mm256_permute_ps(_mm256_sub_ps(s0, s2), 0xB1)), flip));
        _mm256_storeu_ps(F2 + o, _mm256_sub_ps(a, s3));
        _mm256_storeu_ps(F0 + o, _mm256_add_ps(a, s3));
        _mm256_storeu_ps(F1 + o, _mm256_add_ps(s5, t));
        _mm256_storeu_ps(F3 + o, _mm256_sub_ps(s5, t));
    }
    for (; k < m; ++k) {   /* scalar tail, same math */
        kiss_fft_cpx s[6];
        C_MUL(s[0], Fout[k+m],   tw1[k]);
        C_MUL(s[1], Fout[k+2*m], tw23[k]);
        C_MUL(s[2], Fout[k+3*m], tw23[m+k]);
        C_SUB(s[5], Fout[k], s[1]);
        C_ADDTO(Fout[k], s[1]);
        C_ADD(s[3], s[0], s[2]);
        C_SUB(s[4], s[0], s[2]);
        C_SUB(Fout[k+2*m], Fout[k], s[3]);
        C_ADDTO(Fout[k], s[3]);
        if (inverse) {
            Fout[k+m].r = s[5].r - s[4].i;   Fout[k+m].i = s[5].i + s[4].r;
            Fout[k+3*m].r = s[5].r + s[4].i; Fout[k+3*m].i = s[5].i - s[4].r;
        } else {
            Fout[k+m].r = s[5].r + s[4].i;   Fout[k+m].i = s[5].i - s[4].r;
            Fout[k+3*m].r = s[5].r - s[4].i; Fout[k+3*m].i = s[5].i + s[4].r;
        }
    }
}

static int kf_have_avx2(void)
{
    return __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
}
#endif

static void kf_bfly4(kiss_fft_cpx * Fout, const size_t fstride, const kiss_fft_cfg st, const size_t m, const int stage)
{
    (void)stage;
#ifdef KISS_FFT_AVX2
    if (KF_AVX2_OK && st->stage_tw[stage] && kf_have_avx2()) {
        const kiss_fft_cpx *t = st->twiddles + st->stage_tw[stage];
        if (stage == 0)   /* fstride 1: W^k is the start of the twiddles, the table has W^2k, W^3k */
            kf_bfly4_avx2(Fout, st->twiddles, t, m, st->inverse);
        else
            kf_bfly4_avx2(Fout, t, t + m, m, st->inverse);
        return;
    }
#endif
    if (st->inverse) kf_bfly4_dir(Fout,fstride,st,m,1);
    else             kf_bfly4_dir(Fout,fstride,st,m,0);
}

KF_TARGET_CLONES
static void kf_bfly3(
         kiss_fft_cpx * Fout,
         const size_t fstride,
         const kiss_fft_cfg st,
         size_t m
         )
{
     size_t k=m;
     const size_t m2 = 2*m;
     kiss_fft_cpx *tw1,*tw2;
     kiss_fft_cpx scratch[5];
     kiss_fft_cpx epi3;
     epi3 = st->twiddles[fstride*m];

     tw1=tw2=st->twiddles;

     do{
         C_FIXDIV(*Fout,3); C_FIXDIV(Fout[m],3); C_FIXDIV(Fout[m2],3);

         C_MUL(scratch[1],Fout[m] , *tw1);
         C_MUL(scratch[2],Fout[m2] , *tw2);

         C_ADD(scratch[3],scratch[1],scratch[2]);
         C_SUB(scratch[0],scratch[1],scratch[2]);
         tw1 += fstride;
         tw2 += fstride*2;

         Fout[m].r = Fout->r - HALF_OF(scratch[3].r);
         Fout[m].i = Fout->i - HALF_OF(scratch[3].i);

         C_MULBYSCALAR( scratch[0] , epi3.i );

         C_ADDTO(*Fout,scratch[3]);

         Fout[m2].r = Fout[m].r + scratch[0].i;
         Fout[m2].i = Fout[m].i - scratch[0].r;

         Fout[m].r -= scratch[0].i;
         Fout[m].i += scratch[0].r;

         ++Fout;
     }while(--k);
}

KF_TARGET_CLONES
static void kf_bfly5(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        int m
        )
{
    kiss_fft_cpx *Fout0,*Fout1,*Fout2,*Fout3,*Fout4;
    int u;
    kiss_fft_cpx scratch[13];
    kiss_fft_cpx * twiddles = st->twiddles;
    kiss_fft_cpx *tw;
    kiss_fft_cpx ya,yb;
    ya = twiddles[fstride*m];
    yb = twiddles[fstride*2*m];

    Fout0=Fout;
    Fout1=Fout0+m;
    Fout2=Fout0+2*m;
    Fout3=Fout0+3*m;
    Fout4=Fout0+4*m;

    tw=st->twiddles;
    for ( u=0; u<m; ++u ) {
        C_FIXDIV( *Fout0,5); C_FIXDIV( *Fout1,5); C_FIXDIV( *Fout2,5); C_FIXDIV( *Fout3,5); C_FIXDIV( *Fout4,5);
        scratch[0] = *Fout0;

        C_MUL(scratch[1] ,*Fout1, tw[u*fstride]);
        C_MUL(scratch[2] ,*Fout2, tw[2*u*fstride]);
        C_MUL(scratch[3] ,*Fout3, tw[3*u*fstride]);
        C_MUL(scratch[4] ,*Fout4, tw[4*u*fstride]);

        C_ADD( scratch[7],scratch[1],scratch[4]);
        C_SUB( scratch[10],scratch[1],scratch[4]);
        C_ADD( scratch[8],scratch[2],scratch[3]);
        C_SUB( scratch[9],scratch[2],scratch[3]);

        Fout0->r += scratch[7].r + scratch[8].r;
        Fout0->i += scratch[7].i + scratch[8].i;

        scratch[5].r = scratch[0].r + S_MUL(scratch[7].r,ya.r) + S_MUL(scratch[8].r,yb.r);
        scratch[5].i = scratch[0].i + S_MUL(scratch[7].i,ya.r) + S_MUL(scratch[8].i,yb.r);

        scratch[6].r =  S_MUL(scratch[10].i,ya.i) + S_MUL(scratch[9].i,yb.i);
        scratch[6].i = -S_MUL(scratch[10].r,ya.i) - S_MUL(scratch[9].r,yb.i);

        C_SUB(*Fout1,scratch[5],scratch[6]);
        C_ADD(*Fout4,scratch[5],scratch[6]);

        scratch[11].r = scratch[0].r + S_MUL(scratch[7].r,yb.r) + S_MUL(scratch[8].r,ya.r);
        scratch[11].i = scratch[0].i + S_MUL(scratch[7].i,yb.r) + S_MUL(scratch[8].i,ya.r);
        scratch[12].r = - S_MUL(scratch[10].i,yb.i) + S_MUL(scratch[9].i,ya.i);
        scratch[12].i = S_MUL(scratch[10].r,yb.i) - S_MUL(scratch[9].r,ya.i);

        C_ADD(*Fout2,scratch[11],scratch[12]);
        C_SUB(*Fout3,scratch[11],scratch[12]);

        ++Fout0;++Fout1;++Fout2;++Fout3;++Fout4;
    }
}

/* perform the butterfly for one stage of a mixed radix FFT */
KF_TARGET_CLONES
static void kf_bfly_generic(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        int m,
        int p
        )
{
    int u,k,q1,q;
    kiss_fft_cpx * twiddles = st->twiddles;
    kiss_fft_cpx t;
    int Norig = st->nfft;

    kiss_fft_cpx * scratch = (kiss_fft_cpx*)KISS_FFT_TMP_ALLOC(sizeof(kiss_fft_cpx)*p);
    if (scratch == NULL){
        KISS_FFT_ERROR("Memory allocation failed.");
        return;
    }

    for ( u=0; u<m; ++u ) {
        k=u;
        for ( q1=0 ; q1<p ; ++q1 ) {
            scratch[q1] = Fout[ k  ];
            C_FIXDIV(scratch[q1],p);
            k += m;
        }

        k=u;
        for ( q1=0 ; q1<p ; ++q1 ) {
            int twidx=0;
            Fout[ k ] = scratch[0];
            for (q=1;q<p;++q ) {
                twidx += fstride * k;
                if (twidx>=Norig) twidx-=Norig;
                C_MUL(t,scratch[q] , twiddles[twidx] );
                C_ADDTO( Fout[ k ] ,t);
            }
            k += m;
        }
    }
    KISS_FFT_TMP_FREE(scratch);
}

KF_TARGET_CLONES
static
void kf_work(
        kiss_fft_cpx * Fout,
        const kiss_fft_cpx * f,
        const size_t fstride,
        int in_stride,
        int * factors,
        const kiss_fft_cfg st
        )
{
    kiss_fft_cpx * Fout_beg=Fout;
    const int p=*factors++; /* the radix  */
    const int m=*factors++; /* stage's fft length/p */
    const int stage = (int)((factors - st->factors) / 2) - 1;
    const kiss_fft_cpx * Fout_end = Fout + p*m;

#ifdef _OPENMP
    // use openmp extensions at the
    // top-level (not recursive)
    if (fstride==1 && p<=5 && m!=1)
    {
        int k;

        // execute the p different work units in different threads
#       pragma omp parallel for
        for (k=0;k<p;++k)
            kf_work( Fout +k*m, f+ fstride*in_stride*k,fstride*p,in_stride,factors,st);
        // all threads have joined by this point

        switch (p) {
            case 2: kf_bfly2(Fout,fstride,st,m); break;
            case 3: kf_bfly3(Fout,fstride,st,m); break;
            case 4: kf_bfly4(Fout,fstride,st,m,stage); break;
            case 5: kf_bfly5(Fout,fstride,st,m); break;
            default: kf_bfly_generic(Fout,fstride,st,m,p); break;
        }
        return;
    }
#endif

    if (m==1) {
        /* leaf stage: gather and butterfly fused for radix 2 and 4, whose twiddles are all 1 */
        const size_t s = fstride*in_stride;
        switch (p) {
            case 2: {
                kiss_fft_cpx a=f[0], b=f[s];
                C_FIXDIV(a,2); C_FIXDIV(b,2);
                C_ADD(Fout[0],a,b); C_SUB(Fout[1],a,b);
                return;
            }
            case 4: {
                kiss_fft_cpx a=f[0], b=f[s], c=f[2*s], d=f[3*s], s0,s1,s2,s3;
                C_FIXDIV(a,4); C_FIXDIV(b,4); C_FIXDIV(c,4); C_FIXDIV(d,4);
                C_ADD(s0,a,c); C_SUB(s1,a,c); C_ADD(s2,b,d); C_SUB(s3,b,d);
                C_ADD(Fout[0],s0,s2); C_SUB(Fout[2],s0,s2);
                if (st->inverse) {
                    Fout[1].r=s1.r-s3.i; Fout[1].i=s1.i+s3.r; Fout[3].r=s1.r+s3.i; Fout[3].i=s1.i-s3.r;
                } else {
                    Fout[1].r=s1.r+s3.i; Fout[1].i=s1.i-s3.r; Fout[3].r=s1.r-s3.i; Fout[3].i=s1.i+s3.r;
                }
                return;
            }
        }
        do{
            *Fout = *f;
            f += fstride*in_stride;
        }while(++Fout != Fout_end );
    }else{
        do{
            // recursive call:
            // DFT of size m*p performed by doing
            // p instances of smaller DFTs of size m,
            // each one takes a decimated version of the input
            kf_work( Fout , f, fstride*p, in_stride, factors,st);
            f += fstride*in_stride;
        }while( (Fout += m) != Fout_end );
    }

    Fout=Fout_beg;

    // recombine the p smaller DFTs
    switch (p) {
        case 2: kf_bfly2(Fout,fstride,st,m); break;
        case 3: kf_bfly3(Fout,fstride,st,m); break;
        case 4: kf_bfly4(Fout,fstride,st,m,stage); break;
        case 5: kf_bfly5(Fout,fstride,st,m); break;
        default: kf_bfly_generic(Fout,fstride,st,m,p); break;
    }
}

/* twiddle W^i = e^(-+2*pi*i*I/nfft) */
static void kf_twiddle(kiss_fft_cpx * x, size_t i, int nfft, int inverse)
{
    const double pi=3.141592653589793238462643383279502884197169399375105820974944;
    double phase = -2*pi*(double)i / nfft;
    if (inverse)
        phase *= -1;
    kf_cexp(x, phase );
}

/*  facbuf is populated by p1,m1,p2,m2, ...
    where
    p[i] * m[i] = m[i-1]
    m0 = n                  */
static
void kf_factor(int n,int * facbuf)
{
    int p=4;
    double floor_sqrt;
    floor_sqrt = floor( sqrt((double)n) );

    /*factor out powers of 4, powers of 2, then any remaining primes */
    do {
        while (n % p) {
            switch (p) {
                case 4: p = 2; break;
                case 2: p = 3; break;
                default: p += 2; break;
            }
            if (p > floor_sqrt)
                p = n;          /* no more factors, skip to end */
        }
        n /= p;
        *facbuf++ = p;
        *facbuf++ = n;
    } while (n > 1);
}

/*
 *
 * User-callable function to allocate all necessary storage space for the fft.
 *
 * The return value is a contiguous block of memory, allocated with malloc.  As such,
 * It can be freed with free(), rather than a kiss_fft-specific function.
 * */
kiss_fft_cfg kiss_fft_alloc(int nfft,int inverse_fft,void * mem,size_t * lenmem )
{
    KISS_FFT_ALIGN_CHECK(mem)

    kiss_fft_cfg st=NULL;
    // check for overflow condition {memneeded > SIZE_MAX}.
    if (nfft >= (SIZE_MAX - 2*sizeof(struct kiss_fft_state))/sizeof(kiss_fft_cpx))
        return NULL;

    size_t extra = 0; /* per-stage contiguous radix-4 twiddles for the AVX2 butterfly */
    int ntw = nfft;   /* twiddles stored before the tables */
#ifdef KISS_FFT_AVX2
    if (KF_AVX2_OK && nfft > 0 && kf_have_avx2()) {   /* tables are only used by the AVX2 butterfly */
        int fac[2*MAXFACTORS];
        int i = 0, tables_only = 1;
        kf_factor(nfft, fac);
        do {
            const int p = fac[2*i], m = fac[2*i+1];
            if (p == 4 && m >= 4)
                extra += (i ? 3 : 2) * (size_t)m;   /* stage 0 takes W^k from the twiddles */
            else if (!(m == 1 && (p == 2 || p == 4)) && !(i == 0 && p == 2))
                tables_only = 0;   /* a stage that indexes the full twiddle array */
        } while (fac[2*i++ + 1] > 1);
        /* otherwise the only twiddles read directly are W^k, k < m, of the outermost stage */
        if (extra && tables_only)
            ntw = fac[1];
    }
#endif
    size_t memneeded = KISS_FFT_ALIGN_SIZE_UP(sizeof(struct kiss_fft_state)
        + sizeof(kiss_fft_cpx)*(ntw-1+extra)); /* twiddle factors*/

    if ( lenmem==NULL ) {
        st = ( kiss_fft_cfg)KISS_FFT_MALLOC( memneeded );
    }else{
        if (mem != NULL && *lenmem >= memneeded)
            st = (kiss_fft_cfg)mem;
        *lenmem = memneeded;
    }
    if (st) {
        int i;
        st->nfft=nfft;
        st->inverse = inverse_fft;

        for (i=0;i<ntw;++i)
            kf_twiddle(st->twiddles+i, (size_t)i, nfft, st->inverse);

        kf_factor(nfft,st->factors);
        memset(st->stage_tw, 0, sizeof(st->stage_tw));
#ifdef KISS_FFT_AVX2
        if (extra) {
            size_t off = (size_t)ntw, fstride = 1, k;
            int s = 0;
            do {
                const int p = st->factors[2*s], m = st->factors[2*s+1];
                if (p == 4 && m >= 4) {
                    kiss_fft_cpx *t = st->twiddles + off;
                    int q;
                    st->stage_tw[s] = (unsigned)off;
                    for (q = s ? 1 : 2; q <= 3; ++q, t += m, off += (size_t)m)  /* W^qk */
                        for (k = 0; k < (size_t)m; ++k)
                            kf_twiddle(t + k, q*k*fstride, nfft, st->inverse);
                }
                fstride *= (size_t)p;
            } while (st->factors[2*s++ + 1] > 1);
        }
#endif
    }
    return st;
}


void kiss_fft_stride(kiss_fft_cfg st,const kiss_fft_cpx *fin,kiss_fft_cpx *fout,int in_stride)
{
    if (fin == fout) {
        //NOTE: this is not really an in-place FFT algorithm.
        //It just performs an out-of-place FFT into a temp buffer
        if (fout == NULL){
            KISS_FFT_ERROR("fout buffer NULL.");
        return;
        }

        kiss_fft_cpx * tmpbuf = (kiss_fft_cpx*)KISS_FFT_TMP_ALLOC( sizeof(kiss_fft_cpx)*st->nfft);
        if (tmpbuf == NULL){
            KISS_FFT_ERROR("Memory allocation error.");
        return;
        }



        kf_work(tmpbuf,fin,1,in_stride, st->factors,st);
        memcpy(fout,tmpbuf,sizeof(kiss_fft_cpx)*st->nfft);
        KISS_FFT_TMP_FREE(tmpbuf);
    }else{
        kf_work( fout, fin, 1,in_stride, st->factors,st );
    }
}

void kiss_fft(kiss_fft_cfg cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout)
{
    kiss_fft_stride(cfg,fin,fout,1);
}


void kiss_fft_cleanup(void)
{
    // nothing needed any more
}

int kiss_fft_next_fast_size(int n)
{
    if (n <= 0) /* nothing to factor; avoid infinite loop for n<=0 */
        return 1;
    while(1) {
        int m=n;
        while ( (m%2) == 0 ) m/=2;
        while ( (m%3) == 0 ) m/=3;
        while ( (m%5) == 0 ) m/=5;
        if (m<=1)
            break; /* n is completely factorable by twos, threes, and fives */
        n++;
    }
    return n;
}
