#include <stdio.h>
#include <string.h>
#include <memory.h>
#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include "kiss_fft.h"

#define NFFT 1024
int main(int argc, char ** argv)
{
    int k;
    void * st;
    float fs=44100;

    short sampsin[2*NFFT];
    float lmag2[NFFT/2];
    float rmag2[NFFT/2];
    int peakr=0,peakl=0;
    int removedc=1;

    kiss_fft_cpx cbuf[NFFT];
    int nbufs=0;

    st = kiss_fft_alloc(NFFT,0);

    memset( lmag2 , 0 , sizeof(lmag2) );
    memset( rmag2 , 0 , sizeof(rmag2) );

    while ( fread( sampsin , sizeof(short) * 2*NFFT, 1 , stdin ) == 1 ) {
        //perform two ffts in parallel by packing the channels into the real and imaginary
        //
        
        for (k=0;k<NFFT;++k) {
            //cbuf[k].r = 0;
            //cbuf[k].i = 0;
            cbuf[k].r = sampsin[2*k];
            cbuf[k].i = sampsin[2*k+1];
            //cbuf[k].i = sampsin[2*k+1];
        }
        
        if (removedc){
            float dcr=0,dci=0;
            for (k=0;k<NFFT;++k){
                dcr += cbuf[k].r;
                dci += cbuf[k].i;
            }
            dcr /= NFFT;
            dci /= NFFT;

            for (k=0;k<NFFT;++k){
                cbuf[k].r -= dcr;
                cbuf[k].i -= dci;
            }
        }

        kiss_fft( st , cbuf );

        for (k=0;k<NFFT/2;++k) {
            int k2 = (NFFT-k)%NFFT;
            kiss_fft_cpx r,l;

            r.r = (cbuf[k].r + cbuf[k2].r) * 0.5;
            r.i = (cbuf[k].i - cbuf[k2].i) * 0.5;

            l.r = (cbuf[k].i + cbuf[k2].i) * 0.5;
            l.i = (cbuf[k].r - cbuf[k2].r) * 0.5;

            rmag2[k] += r.r * r.r + r.i * r.i;
            lmag2[k] += l.r * l.r + l.i * l.i;
        }
        ++nbufs;
    }
    free(st);

    for (k=0;k<NFFT/2;++k) {
        if (rmag2[peakr] < rmag2[k])
            peakr = k;
        if (lmag2[peakl] < lmag2[k])
            peakl = k;
    }

    printf("%d buffers\n",nbufs);

    printf("peak frequency R:%.3fdB @ %.1f Hz\n",10*log10(rmag2[peakr]) , peakr*fs/NFFT);
    printf("               L:%.3fdB @ %.1f Hz\n",10*log10(lmag2[peakl]) , peakl*fs/NFFT);

    return 0;
}
