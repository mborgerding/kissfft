#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "kiss_fft.h"
#include "kiss_fftr.h"
#include <limits.h>

#define pcpx(c)\
        fprintf(stderr,"%g + %gi\n",(double)((c)->r),(double)((c)->i) )

void fe_fft_test( int N, int bin1,int bin2)
{
    // static kiss_fft_cfg cfg = NULL;
    static kiss_fftr_cfg cfg	= NULL;
    static int lastN			= 0;
    // kiss_fft_cpx *kin = NULL;
    static kiss_fft_cpx *kout	= NULL;
    static kiss_fft_scalar *kin_r = NULL;
    static int halfN = 0;
    int flipN;

    int i;
    double f1 = (bin1*2.0/N)*M_PI;
    double f2 = (bin2*2.0/N)*M_PI;
    double sigpow=0;
    double noisepow=0;

    if (lastN != N) {
        // Free previous structures, if allocated.
        if (cfg)
            free(cfg);
        if (kin_r)
            free(kin_r);
        if (kout)
            free(kout);

        halfN = N / 2;

        // Allocate the KISS FFT configuration structure.
        // Although the memory for these structures will be freed 
        //  when the DLL unloads, to be squeaky clean it should be
        //  freed explicitly.
        // cfg = kiss_fft_alloc( N , 0, NULL, NULL);
        cfg = kiss_fftr_alloc(N , 0, NULL, NULL);

        // kin		= malloc(N * sizeof(kiss_fft_cpx));
        // kin_r		= malloc(N * sizeof(kiss_fft_cpx));
        kin_r		= malloc(N * sizeof(kiss_fft_scalar));
        kout	= malloc(N * sizeof(kiss_fft_cpx));


        lastN = N;
    }

#if FIXED_POINT==32
    long maxrange = LONG_MAX;
#else    
    long maxrange = SHRT_MAX;
#endif
    // Convert the floating point "in" array to fixed point.
    for (i = 0; i < N; i++)
    {
        kin_r[i] =  (maxrange>>1)*cos(f1*i)
                  + (maxrange>>1)*cos(f2*i);
    }

    // Make the kiss FFT call.
    kiss_fftr(cfg, kin_r, kout);

    sigpow = 0;
    //printf("[");
    for (i=0;i < (N/2+1);++i) {
        double tmpr = (double)kout[i].r / (double)maxrange;
        double tmpi = (double)kout[i].i / (double)maxrange;
        double mag2 = tmpr*tmpr + tmpi*tmpi;
        if (i!=0 && i!= N/2)
            mag2 *= 2; // all bins except DC and Nyquist have symmetric counterparts implied
        sigpow += mag2;

        // subtract out the expected result, anything left is noise
        if ( i!=bin1 && i != bin2 ) 
            noisepow += mag2;
        //printf("%d%+di,",(int)kout[i].r,(int)kout[i].i);
    }
    //printf("]\n");
    printf("TEST %d,%d,%d noise @ %fdB\n",N,bin1,bin2,10*log10(noisepow/sigpow +1e-20) );

    return(0);
}


int main(int argc,char ** argv)
{
    int nfft = 16;
    int i,j;
    for (i=0;i<nfft/2;++i) { 
        for (j=i;j<nfft/2;++j) {
            fe_fft_test(nfft,i,j);
        }
    }
    printf("sizeof(kiss_fft_scalar) = %d\n",sizeof(kiss_fft_scalar) );

    return 0;
}
