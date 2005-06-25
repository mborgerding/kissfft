#include "kiss_fftr.h"
#include "_kiss_fft_guts.h"
#include <sys/times.h>
#include <time.h>
#include <unistd.h>

static double cputime(void)
{
    struct tms t;
    times(&t);
    return (double)(t.tms_utime + t.tms_stime)/  sysconf(_SC_CLK_TCK) ;
}

static
double snr_compare( kiss_fft_cpx * vec1,kiss_fft_cpx * vec2, int n)
{
    int k;
    double sigpow=1e-10,noisepow=1e-10,err,snr,scale=0;

#ifdef USE_SIMD
    float *fv1 = (float*)vec1;
    float *fv2 = (float*)vec2;
    for (k=0;k<8*n;++k) {
        sigpow += *fv1 * *fv1;
        err = *fv1 - *fv2;
        noisepow += err*err;
        ++fv1;
        ++fv2;
    }
#else
    for (k=0;k<n;++k) {
        sigpow += (double)vec1[k].r * (double)vec1[k].r + 
                  (double)vec1[k].i * (double)vec1[k].i;
        err = (double)vec1[k].r - (double)vec2[k].r;
        noisepow += err * err;
        err = (double)vec1[k].i - (double)vec2[k].i;
        noisepow += err * err;

        if (vec1[k].r)
            scale +=(double) vec2[k].r / (double)vec1[k].r;
    }
#endif    
    snr = 10*log10( sigpow / noisepow );
    scale /= n;
    if (snr<10) {
        printf( "\npoor snr, try a scaling factor %f\n" , scale );
        exit(1);
    }
    return snr;
}
#define NFFT 8*3*5

#ifndef NUMFFTS
#define NUMFFTS 10000
#endif


int main(void)
{
    double ts,tfft,trfft;
    int i;
    kiss_fft_cpx cin[NFFT];
    kiss_fft_cpx cout[NFFT];
    kiss_fft_cpx sout[NFFT];
    kiss_fft_cfg  kiss_fft_state;
    kiss_fftr_cfg  kiss_fftr_state;

    kiss_fft_scalar rin[NFFT];
    
    srand(time(0));

    for (i=0;i<NFFT;++i) {
#ifdef USE_SIMD
        rin[i] = _mm_set1_ps(rand()-RAND_MAX/2);
        cin[i].i = _mm_set1_ps(0);
#else
        rin[i] = (kiss_fft_scalar)(rand()-RAND_MAX/2);
        cin[i].i = 0;
#endif        
        cin[i].r = rin[i];
    }

    kiss_fft_state = kiss_fft_alloc(NFFT,0,0,0);
    kiss_fftr_state = kiss_fftr_alloc(NFFT,0,0,0);
    kiss_fft(kiss_fft_state,cin,cout);
    kiss_fftr(kiss_fftr_state,rin,sout);
    printf( "nfft=%d, inverse=%d, snr=%g\n",
            NFFT,0, snr_compare(cout,sout,(NFFT/2)+1) );
    ts = cputime();
    for (i=0;i<NUMFFTS;++i) {
        kiss_fft(kiss_fft_state,cin,cout);
    }
    tfft = cputime() - ts;
    
    ts = cputime();
    for (i=0;i<NUMFFTS;++i) {
        kiss_fftr( kiss_fftr_state, rin, cout );
        /* kiss_fftri(kiss_fftr_state,cout,rin); */
    }
    trfft = cputime() - ts;

    printf("%d complex ffts took %gs, real took %gs\n",NUMFFTS,tfft,trfft);

    free(kiss_fft_state);
    free(kiss_fftr_state);

    kiss_fft_state = kiss_fft_alloc(NFFT,1,0,0);
    kiss_fftr_state = kiss_fftr_alloc(NFFT,1,0,0);

    kiss_fft(kiss_fft_state,cout,cin);
    kiss_fftri(kiss_fftr_state,cout,rin);

    for (i=0;i<NFFT;++i) {
        sout[i].r = rin[i];
#ifdef USE_SIMD
        sout[i].i ^= sout[i].i;
#else        
        sout[i].i = 0;
#endif        
    }
    
    printf( "nfft=%d, inverse=%d, snr=%g\n",
            NFFT,1, snr_compare(cin,cin,NFFT/2) );
    free(kiss_fft_state);
    free(kiss_fftr_state);

    return 0;
}

