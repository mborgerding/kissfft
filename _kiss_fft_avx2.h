/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */

/* AVX2/FMA butterflies, private to kiss_fft.c, which includes this file twice: with KF_AVX2_PS
 * defined for float (4 interleaved complex values per __m256) and without it for double (2 per
 * __m256d). Only the instantiation matching kiss_fft_scalar is ever called.
 *
 * Each kernel processes k = 0..m-1 for one stage, vector by vector; the last m % KF_L values go
 * through the same step with masked loads/stores. Twiddles come in contiguous runs of m:
 * tw1 = W^k and twr = W^2k | W^3k | ... (see stage_tw in _kiss_fft_guts.h).
 *
 * Sign flips are integer xors on purpose: as float constants (+-0.0) they get folded away under
 * -ffast-math (-fno-signed-zeros), which the default build uses. */

#ifdef KF_AVX2_PS
#  define KF_T float
#  define KF_V __m256
#  define KF_L 4
#  define KF_N(name) name##_ps
#  define KF_LOAD(p, full, mk) ((full) ? _mm256_loadu_ps(p) : _mm256_maskload_ps(p, mk))
#  define KF_STORE(p, v, full, mk) \
       do { if (full) _mm256_storeu_ps(p, v); else _mm256_maskstore_ps(p, mk, v); } while (0)
#  define KF_ADD _mm256_add_ps
#  define KF_SUB _mm256_sub_ps
#  define KF_MUL _mm256_mul_ps
#  define KF_SET1 _mm256_set1_ps
#  define KF_FMADDSUB _mm256_fmaddsub_ps
#  define KF_SWAP(v) _mm256_permute_ps(v, 0xB1)  /* (re,im) -> (im,re) */
#  define KF_DUPRE(v) _mm256_moveldup_ps(v)
#  define KF_DUPIM(v) _mm256_movehdup_ps(v)
#  define KF_TOI _mm256_castps_si256
#  define KF_FROMI _mm256_castsi256_ps
#  define KF_NEG_IM _mm256_set_epi32(INT_MIN, 0, INT_MIN, 0, INT_MIN, 0, INT_MIN, 0)
#  define KF_NEG_RE _mm256_set_epi32(0, INT_MIN, 0, INT_MIN, 0, INT_MIN, 0, INT_MIN)
#  define KF_TAILMASK(r) _mm256_cmpgt_epi32(_mm256_set1_epi32(2*(int)(r)), _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7))
#else
#  define KF_T double
#  define KF_V __m256d
#  define KF_L 2
#  define KF_N(name) name##_pd
#  define KF_LOAD(p, full, mk) ((full) ? _mm256_loadu_pd(p) : _mm256_maskload_pd(p, mk))
#  define KF_STORE(p, v, full, mk) \
       do { if (full) _mm256_storeu_pd(p, v); else _mm256_maskstore_pd(p, mk, v); } while (0)
#  define KF_ADD _mm256_add_pd
#  define KF_SUB _mm256_sub_pd
#  define KF_MUL _mm256_mul_pd
#  define KF_SET1 _mm256_set1_pd
#  define KF_FMADDSUB _mm256_fmaddsub_pd
#  define KF_SWAP(v) _mm256_permute_pd(v, 0x5)
#  define KF_DUPRE(v) _mm256_movedup_pd(v)
#  define KF_DUPIM(v) _mm256_permute_pd(v, 0xF)
#  define KF_TOI _mm256_castpd_si256
#  define KF_FROMI _mm256_castsi256_pd
#  define KF_NEG_IM _mm256_set_epi64x(INT64_MIN, 0, INT64_MIN, 0)
#  define KF_NEG_RE _mm256_set_epi64x(0, INT64_MIN, 0, INT64_MIN)
#  define KF_TAILMASK(r) _mm256_cmpgt_epi64(_mm256_set1_epi64x(2*(long long)(r)), _mm256_setr_epi64x(0, 1, 2, 3))
#endif

#define KF_STEP KF_AVX2_FN __attribute__((always_inline)) static inline

/* a*w for interleaved complex values */
KF_STEP KF_V KF_N(kf_cmul)(KF_V a, KF_V w)
{
    return KF_FMADDSUB(a, KF_DUPRE(w), KF_MUL(KF_SWAP(a), KF_DUPIM(w)));
}

/* (re,im) -> (im,re), then flip the sign of the lanes in mask: NEG_IM gives -i*v, NEG_RE gives i*v */
KF_STEP KF_V KF_N(kf_rot)(KF_V v, __m256i mask)
{
    return KF_FROMI(_mm256_xor_si256(KF_TOI(KF_SWAP(v)), mask));
}

KF_STEP void KF_N(kf_step2)(KF_T *F0, KF_T *F1, const KF_T *T1, const int full, const __m256i mk)
{
    KF_V a = KF_LOAD(F0, full, mk);
    KF_V t = KF_N(kf_cmul)(KF_LOAD(F1, full, mk), KF_LOAD(T1, full, mk));
    KF_STORE(F1, KF_SUB(a, t), full, mk);
    KF_STORE(F0, KF_ADD(a, t), full, mk);
}

KF_STEP void KF_N(kf_step3)(KF_T *F0, KF_T *F1, KF_T *F2, const KF_T *T1, const KF_T *T2,
                            const KF_V epi, const int full, const __m256i mk)
{
    KF_V a  = KF_LOAD(F0, full, mk);
    KF_V s1 = KF_N(kf_cmul)(KF_LOAD(F1, full, mk), KF_LOAD(T1, full, mk));
    KF_V s2 = KF_N(kf_cmul)(KF_LOAD(F2, full, mk), KF_LOAD(T2, full, mk));
    KF_V s3 = KF_ADD(s1, s2);
    KF_V b  = KF_SUB(a, KF_MUL(s3, KF_SET1((KF_T)0.5)));
    KF_V t  = KF_N(kf_rot)(KF_MUL(KF_SUB(s1, s2), epi), KF_NEG_IM);
    KF_STORE(F0, KF_ADD(a, s3), full, mk);
    KF_STORE(F1, KF_SUB(b, t), full, mk);
    KF_STORE(F2, KF_ADD(b, t), full, mk);
}

KF_STEP void KF_N(kf_step4)(KF_T *F0, KF_T *F1, KF_T *F2, KF_T *F3,
                            const KF_T *T1, const KF_T *T2, const KF_T *T3,
                            const __m256i rot, const int full, const __m256i mk)
{
    KF_V a0 = KF_LOAD(F0, full, mk);
    KF_V s0 = KF_N(kf_cmul)(KF_LOAD(F1, full, mk), KF_LOAD(T1, full, mk));
    KF_V s1 = KF_N(kf_cmul)(KF_LOAD(F2, full, mk), KF_LOAD(T2, full, mk));
    KF_V s2 = KF_N(kf_cmul)(KF_LOAD(F3, full, mk), KF_LOAD(T3, full, mk));
    KF_V s5 = KF_SUB(a0, s1);
    KF_V a  = KF_ADD(a0, s1);
    KF_V s3 = KF_ADD(s0, s2);
    KF_V t  = KF_N(kf_rot)(KF_SUB(s0, s2), rot);
    KF_STORE(F2, KF_SUB(a, s3), full, mk);
    KF_STORE(F0, KF_ADD(a, s3), full, mk);
    KF_STORE(F1, KF_ADD(s5, t), full, mk);
    KF_STORE(F3, KF_SUB(s5, t), full, mk);
}

/* ya = W^m, yb = W^2m: the 5th roots of unity used by kf_bfly5 */
KF_STEP void KF_N(kf_step5)(KF_T *F0, KF_T *F1, KF_T *F2, KF_T *F3, KF_T *F4,
                            const KF_T *T1, const KF_T *T2, const KF_T *T3, const KF_T *T4,
                            const KF_V yar, const KF_V yai, const KF_V ybr, const KF_V ybi,
                            const int full, const __m256i mk)
{
    KF_V s0 = KF_LOAD(F0, full, mk);
    KF_V s1 = KF_N(kf_cmul)(KF_LOAD(F1, full, mk), KF_LOAD(T1, full, mk));
    KF_V s2 = KF_N(kf_cmul)(KF_LOAD(F2, full, mk), KF_LOAD(T2, full, mk));
    KF_V s3 = KF_N(kf_cmul)(KF_LOAD(F3, full, mk), KF_LOAD(T3, full, mk));
    KF_V s4 = KF_N(kf_cmul)(KF_LOAD(F4, full, mk), KF_LOAD(T4, full, mk));
    KF_V s7 = KF_ADD(s1, s4), s10 = KF_SUB(s1, s4);
    KF_V s8 = KF_ADD(s2, s3), s9  = KF_SUB(s2, s3);
    KF_V s5  = KF_ADD(s0, KF_ADD(KF_MUL(s7, yar), KF_MUL(s8, ybr)));
    KF_V s6  = KF_N(kf_rot)(KF_ADD(KF_MUL(s10, yai), KF_MUL(s9, ybi)), KF_NEG_IM);
    KF_V s11 = KF_ADD(s0, KF_ADD(KF_MUL(s7, ybr), KF_MUL(s8, yar)));
    KF_V s12 = KF_N(kf_rot)(KF_SUB(KF_MUL(s10, ybi), KF_MUL(s9, yai)), KF_NEG_RE);
    KF_STORE(F0, KF_ADD(s0, KF_ADD(s7, s8)), full, mk);
    KF_STORE(F1, KF_SUB(s5, s6), full, mk);
    KF_STORE(F4, KF_ADD(s5, s6), full, mk);
    KF_STORE(F2, KF_ADD(s11, s12), full, mk);
    KF_STORE(F3, KF_SUB(s11, s12), full, mk);
}

/* run STEP over k = 0..m-1: full vectors, then one masked vector for the remainder */
#define KF_FOR_EACH_K(STEP)                                            \
    do {                                                               \
        const __m256i all_ = _mm256_set1_epi32(-1);                    \
        size_t k_;                                                     \
        for (k_ = 0; k_ + KF_L <= m; k_ += KF_L) { const size_t o = 2*k_; STEP(1, all_); } \
        if (k_ < m) { const size_t o = 2*k_; STEP(0, KF_TAILMASK(m - k_)); } \
    } while (0)

KF_AVX2_FN static void KF_N(kf_bfly2_avx2)(kiss_fft_cpx * Fout, const kiss_fft_cpx * tw1, const size_t m)
{
    KF_T *F = (KF_T*)Fout;
    const KF_T *T1 = (const KF_T*)tw1;
#define KF_S2(full, mk) KF_N(kf_step2)(F + o, F + 2*m + o, T1 + o, full, mk)
    KF_FOR_EACH_K(KF_S2);
#undef KF_S2
}

KF_AVX2_FN static void KF_N(kf_bfly3_avx2)(kiss_fft_cpx * Fout, const kiss_fft_cpx * tw1, const kiss_fft_cpx * twr,
                                           const size_t m, const kiss_fft_cpx epi3)
{
    KF_T *F = (KF_T*)Fout;
    const KF_T *T1 = (const KF_T*)tw1, *T2 = (const KF_T*)twr;
    const KF_V epi = KF_SET1((KF_T)epi3.i);
#define KF_S3(full, mk) KF_N(kf_step3)(F + o, F + 2*m + o, F + 4*m + o, T1 + o, T2 + o, epi, full, mk)
    KF_FOR_EACH_K(KF_S3);
#undef KF_S3
}

KF_AVX2_FN static void KF_N(kf_bfly4_avx2)(kiss_fft_cpx * Fout, const kiss_fft_cpx * tw1, const kiss_fft_cpx * twr,
                                           const size_t m, const int inverse)
{
    KF_T *F = (KF_T*)Fout;
    const KF_T *T1 = (const KF_T*)tw1, *T2 = (const KF_T*)twr, *T3 = (const KF_T*)(twr + m);
    const __m256i rot = inverse ? KF_NEG_RE : KF_NEG_IM;   /* times +i / -i */
#define KF_S4(full, mk) KF_N(kf_step4)(F + o, F + 2*m + o, F + 4*m + o, F + 6*m + o, \
                                       T1 + o, T2 + o, T3 + o, rot, full, mk)
    KF_FOR_EACH_K(KF_S4);
#undef KF_S4
}

KF_AVX2_FN static void KF_N(kf_bfly5_avx2)(kiss_fft_cpx * Fout, const kiss_fft_cpx * tw1, const kiss_fft_cpx * twr,
                                           const size_t m, const kiss_fft_cpx ya, const kiss_fft_cpx yb)
{
    KF_T *F = (KF_T*)Fout;
    const KF_T *T1 = (const KF_T*)tw1, *T2 = (const KF_T*)twr,
               *T3 = (const KF_T*)(twr + m), *T4 = (const KF_T*)(twr + 2*m);
    const KF_V yar = KF_SET1((KF_T)ya.r), yai = KF_SET1((KF_T)ya.i);
    const KF_V ybr = KF_SET1((KF_T)yb.r), ybi = KF_SET1((KF_T)yb.i);
#define KF_S5(full, mk) KF_N(kf_step5)(F + o, F + 2*m + o, F + 4*m + o, F + 6*m + o, F + 8*m + o, \
                                       T1 + o, T2 + o, T3 + o, T4 + o, yar, yai, ybr, ybi, full, mk)
    KF_FOR_EACH_K(KF_S5);
#undef KF_S5
}

#undef KF_FOR_EACH_K
#undef KF_STEP
#undef KF_T
#undef KF_V
#undef KF_L
#undef KF_N
#undef KF_LOAD
#undef KF_STORE
#undef KF_ADD
#undef KF_SUB
#undef KF_MUL
#undef KF_SET1
#undef KF_FMADDSUB
#undef KF_SWAP
#undef KF_DUPRE
#undef KF_DUPIM
#undef KF_TOI
#undef KF_FROMI
#undef KF_NEG_IM
#undef KF_NEG_RE
#undef KF_TAILMASK
