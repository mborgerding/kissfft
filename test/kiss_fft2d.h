#ifndef KISS_FFT2D_H
#define KISS_FFT2D_H
#include "kiss_fft.h"

/* allocate a 2-dimensional FFT
   the data should be stored rowwise,
   in other words, an array made up of row[0], then row[1], etc
 */
void * kiss_fft2d_alloc(int nrows,int ncols,int inverse_fft,void * mem,size_t * lenmem);
void kiss_fft2d(const void* cfg_from_alloc , const kiss_fft_cpx *fin,kiss_fft_cpx *fout );

#endif
