/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */

#ifndef KISS_FFT_H
#define KISS_FFT_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

// Define KISS_FFT_SHARED macro to properly export symbols
#ifdef KISS_FFT_SHARED
# ifdef _WIN32
#  ifdef KISS_FFT_BUILD
#   define KISS_FFT_API __declspec(dllexport)
#  else
#   define KISS_FFT_API __declspec(dllimport)
#  endif
# else
#  define KISS_FFT_API __attribute__ ((visibility ("default")))
# endif
#else
# define KISS_FFT_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 ATTENTION!
 If you would like a :
 -- a utility that will handle the caching of fft objects
 -- real-only (no imaginary time component ) FFT
 -- a multi-dimensional FFT
 -- a command-line utility to perform ffts
 -- a command-line utility to perform fast-convolution filtering

 Then see kfc.h kiss_fftr.h kiss_fftnd.h fftutil.c kiss_fastfir.c
  in the tools/ directory.
*/

/* User may override KISS_FFT_MALLOC and/or KISS_FFT_FREE. */
#ifdef USE_SIMD
# ifdef USE_SIMD_SIMDE
#  include <simde/x86/sse.h>
#  define kiss_fft_scalar simde__m128
#  define KISS_FFT_SET1_PS(x) simde_mm_set1_ps(x)
#  ifndef KISS_FFT_MALLOC
#   include <stdlib.h>
#   if defined(_WIN32)
#    define KISS_FFT_MALLOC(nbytes) _aligned_malloc(nbytes,16)
#    define KISS_FFT_FREE _aligned_free
#   elif (defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L) || \
         (defined(__cplusplus) && __cplusplus >= 201703L)
     /* C11 aligned_alloc or C++17 */
#    define KISS_FFT_MALLOC(nbytes) aligned_alloc(16, ((nbytes) + 15) & ~(size_t)15)
#    define KISS_FFT_FREE free
#   else
     /* POSIX posix_memalign */
     static inline void* kiss_fft_simde_malloc(size_t nbytes) {
         void* ptr = NULL;
         if (posix_memalign(&ptr, 16, nbytes) != 0) return NULL;
         return ptr;
     }
#    define KISS_FFT_MALLOC(nbytes) kiss_fft_simde_malloc(nbytes)
#    define KISS_FFT_FREE free
#   endif
#  endif
#  define KISS_FFT_ALIGN_CHECK(ptr)
#  define KISS_FFT_ALIGN_SIZE_UP(size) ((size + 15UL) & ~0xFUL)
# else
#  include <xmmintrin.h>
#  define kiss_fft_scalar __m128
#  define KISS_FFT_SET1_PS(x) _mm_set1_ps(x)
#  ifndef KISS_FFT_MALLOC
#   define KISS_FFT_MALLOC(nbytes) _mm_malloc(nbytes,16)
#   define KISS_FFT_ALIGN_CHECK(ptr)
#   define KISS_FFT_ALIGN_SIZE_UP(size) ((size + 15UL) & ~0xFUL)
#  endif
#  ifndef KISS_FFT_FREE
#   define KISS_FFT_FREE _mm_free
#  endif
# endif
#else
# define KISS_FFT_ALIGN_CHECK(ptr)
# define KISS_FFT_ALIGN_SIZE_UP(size) (size)
# ifndef KISS_FFT_MALLOC
#  define KISS_FFT_MALLOC malloc
# endif
# ifndef KISS_FFT_FREE
#  define KISS_FFT_FREE free
# endif
#endif


#ifdef FIXED_POINT
#include <stdint.h>
# if (FIXED_POINT == 32)
#  define kiss_fft_scalar int32_t
# else	
#  define kiss_fft_scalar int16_t
# endif
#else
# ifndef kiss_fft_scalar
/*  default is float */
#   define kiss_fft_scalar float
# endif
#endif

typedef struct {
    kiss_fft_scalar r;
    kiss_fft_scalar i;
}kiss_fft_cpx;

typedef struct kiss_fft_state* kiss_fft_cfg;

/* 
 *  kiss_fft_alloc
 *  
 *  Initialize a FFT (or IFFT) algorithm's cfg/state buffer.
 *
 *  typical usage:      kiss_fft_cfg mycfg=kiss_fft_alloc(1024,0,NULL,NULL);
 *
 *  The return value from fft_alloc is a cfg buffer used internally
 *  by the fft routine or NULL.
 *
 *  If lenmem is NULL, then kiss_fft_alloc will allocate a cfg buffer using malloc.
 *  The returned value should be free()d when done to avoid memory leaks.
 *  
 *  The state can be placed in a user supplied buffer 'mem':
 *  If lenmem is not NULL and mem is not NULL and *lenmem is large enough,
 *      then the function places the cfg in mem and the size used in *lenmem
 *      and returns mem.
 *  
 *  If lenmem is not NULL and ( mem is NULL or *lenmem is not large enough),
 *      then the function returns NULL and places the minimum cfg 
 *      buffer size in *lenmem.
 * */

kiss_fft_cfg KISS_FFT_API kiss_fft_alloc(int nfft,int inverse_fft,void * mem,size_t * lenmem);

/*
 * kiss_fft(cfg,in_out_buf)
 *
 * Perform an FFT on a complex input buffer.
 * for a forward FFT,
 * fin should be  f[0] , f[1] , ... ,f[nfft-1]
 * fout will be   F[0] , F[1] , ... ,F[nfft-1]
 * Note that each element is complex and can be accessed like
    f[k].r and f[k].i
 * */
void KISS_FFT_API kiss_fft(kiss_fft_cfg cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout);

/*
 A more generic version of the above function. It reads its input from every Nth sample.
 * */
void KISS_FFT_API kiss_fft_stride(kiss_fft_cfg cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout,int fin_stride);

/* If kiss_fft_alloc allocated a buffer, it is one contiguous 
   buffer and can be simply free()d when no longer needed*/
#define kiss_fft_free KISS_FFT_FREE

/*
 Cleans up some memory that gets managed internally. Not necessary to call, but it might clean up 
 your compiler output to call this before you exit.
*/
void KISS_FFT_API kiss_fft_cleanup(void);
	

/*
 * Returns the smallest integer k, such that k>=n and k has only "fast" factors (2,3,5)
 */
int KISS_FFT_API kiss_fft_next_fast_size(int n);

/* for real ffts, we need an even size */
#define kiss_fftr_next_fast_size_real(n) \
        (kiss_fft_next_fast_size( ((n)+1)>>1)<<1)

#ifdef __cplusplus
} 
#endif

#endif
