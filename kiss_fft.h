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
 *  kiss_fft_alloc
 *  
 *  Initialize a FFT (or IFFT) algorithm's cfg/state buffer.
 *
 *  typical usage:      void * mycfg=kiss_fft_alloc(1024,0,NULL,NULL);
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
void* kiss_fft_alloc(int nfft,int inverse_fft,void * mem,size_t * lenmem); 

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
void kiss_fft(const void * cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout);


/* allocate a 2-dimensional FFT
   the data should be stored rowwise,
   in other words, an array made up of row[0], then row[1], etc
 */
void * kiss_fft2d_alloc(int nrows,int ncols,int inverse_fft,void * mem,size_t * lenmem);
void kiss_fft2d(const void* cfg_from_alloc , const kiss_fft_cpx *fin,kiss_fft_cpx *fout );

/* Real optimized version can save about 45% cpu time vs. complex fft of a real seq.
 */
void * kiss_fftr_alloc(int nfft,int inverse_fft,void * mem, size_t * lenmem);
void kiss_fftr(const void * cfg,const kiss_fft_scalar *timedata,kiss_fft_cpx *freqdata);
void kiss_fftri(const void * cfg,const kiss_fft_cpx *freqdata,kiss_fft_scalar *timedata);

/* 2d Real optimized version 
 */
void * kiss_fftr2d_alloc(int nfft,int inverse_fft,void * mem, size_t * lenmem);
void kiss_fftr2d(const void * cfg,const kiss_fft_scalar *timedata,kiss_fft_cpx *freqdata);
void kiss_fftr2di(const void * cfg,const kiss_fft_cpx *freqdata,kiss_fft_scalar *timedata);



/* when done with the cfg for a given fft size and direction, simply free it*/
#define kiss_fft_free free

#endif
