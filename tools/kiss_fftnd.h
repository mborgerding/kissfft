#ifndef KISS_FFTND_H
#define KISS_FFTND_H

#include "kiss_fft.h"

#ifdef __cplusplus
extern "C" {
#endif

void * kiss_fftnd_alloc(int *dims,int ndims,int inverse_fft,void*mem,size_t*lenmem);
void kiss_fftnd(void * cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout);

#ifdef __cplusplus
}
#endif
#endif
