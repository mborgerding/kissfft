#ifndef KISS_FFT_H
#define KISS_FFT_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>


#ifdef FIXED_POINT
# define kiss_fft_scalar short
#else
# ifndef kiss_fft_scalar
#   define kiss_fft_scalar float
# endif
#endif


typedef struct {
    kiss_fft_scalar r;
    kiss_fft_scalar i;
}kiss_fft_cpx;

/* 
 *  fft_alloc
 *  
 *  Initialize a FFT (or IFFT) algorithm.
 *
 *  The return value from fft_alloc is a cfg buffer used internally
 *  by the fft routine.
 *
 *  Call free() on it when done using it to avoid memory leaks.
 * */
void* kiss_fft_alloc(int nfft,int inverse_fft); 
/* free() the state when done using it */

/*
 * kiss_fft(cfg,in_out_buf)
 *
 * Perform an FFT on a complex input buffer.
 * for a forward FFT,
 * fin should be  f[0] , f[1] , ... ,f[nfft-1]
 * fout will be   F[0] , F[1] , ... ,F[nfft-1]
 * Note that each element is complex and can be accessed like
    f[k].r and f[k].i

   Apologies to previous users of KISS FFT, this function has changed from 2 args to 3 args.
   To maintain the original behavior , use fout == fin 
 * */

void kiss_fft(const void * cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout);


/* allocate a 2-dimensional FFT
   the data should be stored rowwise,
   in other words, an array made up of row[0], then row[1], etc
 */
void * kiss_fft2d_alloc(int nrows,int ncols,int inverse_fft);
void kiss_fft2d(const void* cfg_from_alloc , const kiss_fft_cpx *fin,kiss_fft_cpx *fout );

/* Real optimized version can save about 45% cpu time vs. complex fft of a real seq.
 */
void * kiss_fftr_alloc(int nfft,int inverse_fft);
void kiss_fftr(const void * cfg,const kiss_fft_scalar *timedata,kiss_fft_cpx *freqdata);
void kiss_fftri(const void * cfg,const kiss_fft_cpx *freqdata,kiss_fft_scalar *timedata);


/* when done with the cfg for a given fft size and direction, simply free it*/
#define kiss_fft_free free

#endif
