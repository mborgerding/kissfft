#ifndef KISS_FTR_H
#define KISS_FTR_H

#include "kiss_fft.h"

/* Real optimized version can save about 45% cpu time vs. complex fft of a real seq.
 */
void * kiss_fftr_alloc(int nfft,int inverse_fft,void * mem, size_t * lenmem);
void kiss_fftr(const void * cfg,const kiss_fft_scalar *timedata,kiss_fft_cpx *freqdata);
void kiss_fftri(const void * cfg,const kiss_fft_cpx *freqdata,kiss_fft_scalar *timedata);

#endif
