#ifndef KISS_FFT_H
#define KISS_FFT_H

#ifdef FIXED_POINT

# define kiss_fft_scalar short
# define FIXED_POINT_FRAC_BITS 15

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
 *  Initialize a radix 2 FFT (or IFFT) algorithm.
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
 * Perform an in-place FFT on a complex input buffer.
 * for a forward FFT,
 * the input should be  f[0] , f[1] , ... ,f[nfft-1]
 * the output will be   F[0] , F[1] , ... ,F[nfft-1]
 * Note that each element is complex.
 * */
void kiss_fft( const void* cfg_from_alloc , kiss_fft_cpx *f ); /* call for each buffer */

/* when done with the cfg for a given fft size and direction, simply free it*/
#define kiss_fft_free free

#endif
