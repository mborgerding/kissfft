#include "kiss_fft.h"


void * kiss_fftr_alloc(int nfft,int inverse_fft);
void kiss_fftr(const void * cfg,const kiss_fft_scalar *fin,kiss_fft_cpx *fout);

double snr_compare( kiss_fft_cpx * vec1,kiss_fft_cpx * vec2, int n)
{
    int k;
    double sigpow,noisepow,err,snr,scale=0;
    sigpow = noisepow = .000000000000000000000000000001; 

    for (k=0;k<n;++k) {
        sigpow += vec1[k].r * vec1[k].i + 
                  vec1[k].i * vec1[k].i;
        err = vec1[k].r - vec2[k].r;
        noisepow += err * err;
        err = vec1[k].i - vec2[k].i;
        noisepow += err * err;

        if (vec1[k].r)
            scale += vec2[k].r / vec1[k].r;
    }
    snr = 10*log10( sigpow / noisepow );
    scale /= n;
    if (snr<10)
        printf( "\npoor snr, try a scaling factor %f\n" , scale );
    return snr;
}

#define RANDOM
#ifndef RANDOM
#define NFFT 8
#else
#define NFFT 120
#endif

void pcpx(const char * msg, kiss_fft_cpx * c)
{
    printf("%s: %g + %gi\n",msg,c->r,c->i);
}

int main()
{
    int i;
    kiss_fft_cpx cin[NFFT];
    kiss_fft_scalar sin[NFFT] = {0.309655,0.815653,0.768570,0.591841,0.404767,0.637617,0.007803,0.012665};
    kiss_fft_cpx cout[NFFT];
    kiss_fft_cpx sout[NFFT];
    
    const void * kiss_fft_state;
    const void * kiss_fftr_state;
    int inverse = 0;
    
    kiss_fft_state = kiss_fft_alloc(NFFT,inverse);
    kiss_fftr_state = kiss_fftr_alloc(NFFT,inverse);

    for (i=0;i<NFFT;++i) {
#ifdef RANDOM        
        sin[i] = (kiss_fft_scalar)(rand()-RAND_MAX/2);
#endif        
        cin[i].r = sin[i];
        cin[i].i = 0;
/*        printf("in[%d]",i);pcpx("",cin+i); */
    }

    kiss_fft(kiss_fft_state,cin,cout);
    kiss_fftr(kiss_fftr_state,sin,sout);

    printf( "nfft=%d, inverse=%d, snr=%g\n",
            NFFT,inverse, snr_compare(cout,sout,NFFT/2) );
    
    return 0;
}

