#include "kiss_fft.h"
#include <sys/times.h>
#include <unistd.h>

static double cputime()
{
    struct tms t;
    times(&t);
    return (double)(t.tms_utime + t.tms_stime)/  sysconf(_SC_CLK_TCK) ;
}



double snr_compare( kiss_fft_cpx * vec1,kiss_fft_cpx * vec2, int n)
{
    int k;
    double sigpow,noisepow,err,snr,scale=0;
    sigpow = noisepow = .00000000000000000001; 

    for (k=0;k<n;++k) {
        sigpow += vec1[k].r * vec1[k].r + 
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
#define NFFT 1800
#endif

#ifndef NUMFFTS
#define NUMFFTS 10000
#endif

void pcpx(const char * msg, kiss_fft_cpx * c)
{
    printf("%s: %g + %gi\n",msg,c->r,c->i);
}

int main()
{
    double ts,tfft,trfft;
    int i;
    kiss_fft_cpx cin[NFFT];
    kiss_fft_scalar sin[NFFT] = {0.309655,0.815653,0.768570,0.591841,0.404767,0.637617,0.007803,0.012665};
    kiss_fft_cpx cout[NFFT];
    kiss_fft_cpx sout[NFFT];
    
    void * kiss_fft_state;
    void * kiss_fftr_state;
    

    for (i=0;i<NFFT;++i) {
#ifdef RANDOM        
        sin[i] = (kiss_fft_scalar)(rand()-RAND_MAX/2);
#endif        
        cin[i].r = sin[i];
        cin[i].i = 0;
/*        printf("in[%d]",i);pcpx("",cin+i); */
    }

    kiss_fft_state = kiss_fft_alloc(NFFT,0);
    kiss_fftr_state = kiss_fftr_alloc(NFFT,0);
    kiss_fft(kiss_fft_state,cin,cout);
    kiss_fftr(kiss_fftr_state,sin,sout);
    printf( "nfft=%d, inverse=%d, snr=%g\n",
            NFFT,0, snr_compare(cout,sout,(NFFT/2)+1) );
    free(kiss_fft_state);
    free(kiss_fftr_state);

    kiss_fft_state = kiss_fft_alloc(NFFT,1);
    kiss_fftr_state = kiss_fftr_alloc(NFFT,1);

    kiss_fft(kiss_fft_state,cout,cin);
    kiss_fftri(kiss_fftr_state,cout,sin);

    for (i=0;i<NFFT;++i) {
        sout[i].r = sin[i];
        sout[i].i = 0;
        /* printf("sin[%d] = %f\t",i,sin[i]); 
        printf("cin[%d]",i);pcpx("",cin+i); 
        printf("sout[%d]",i);pcpx("",sout+i); */ 
    }
    
    printf( "nfft=%d, inverse=%d, snr=%g\n",
            NFFT,1, snr_compare(cin,cin,NFFT/2) );
#ifdef RANDOM        
    ts = cputime();
    for (i=0;i<NUMFFTS;++i) {
        kiss_fft(kiss_fft_state,cin,cout);
    }
    tfft = cputime() - ts;
    
    ts = cputime();
    for (i=0;i<NUMFFTS;++i) {
        /* kiss_fftr(kiss_fftr_state,sin,cout); */
        kiss_fftri(kiss_fftr_state,cout,sin);
    }
    trfft = cputime() - ts;

    printf("%d complex ffts took %gs, real took %gs\n",NUMFFTS,tfft,trfft);
#endif        
    free(kiss_fft_state);
    free(kiss_fftr_state);

    return 0;
}

