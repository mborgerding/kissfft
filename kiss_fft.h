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
// free() the state when done using it

/*
 * kiss_fft
 *
 * Perform an in-place FFT on a complex input buffer.
 *
 * the input buffer should be real(f[0]) , imag(f[0]), ... ,real(f[nfft-1]) , imag(f[nfft-1])
 * the the output will be     real(F[0]) , imag(F[0]), ... ,real(F[nfft-1]) , imag(F[nfft-1])
 *
 * */
void kiss_fft( const void* cfg , kiss_fft_cpx *f ); // call for each buffer 

// when done with the cfg for a given fft size and direction, simply free it
#define kiss_fft_free free

#endif
