#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <unistd.h>
#include "kiss_fft.h"

#include "pstats.h"


int main(int argc,char ** argv)
{
    int nfft=1024;
    int isinverse=0;
    int numffts=1000,i;
    kiss_fft_cpx * buf;
    kiss_fft_cpx * bufout;
    kiss_fft_cfg st;

    while (1) {
      int c = getopt (argc, argv, "n:ix:");
      if (c == -1)
        break;
      switch (c) {
      case 'n':
        nfft = atoi (optarg);
        if (nfft != kiss_fft_next_fast_size(nfft) ) {
            int ng = kiss_fft_next_fast_size(nfft);
            fprintf(stderr,"warning: %d might be a better choice for speed than %d\n",ng,nfft);
        }
        break;
      case 'x':
        numffts = atoi (optarg);
        break;
      case 'i':
        isinverse = 1;
        break;
      }
    }
#ifdef USE_SIMD        
    buf=(kiss_fft_cpx*)memalign(sizeof(kiss_fft_cpx),sizeof(kiss_fft_cpx) * nfft);
    bufout=(kiss_fft_cpx*)memalign(sizeof(kiss_fft_cpx),sizeof(kiss_fft_cpx) * nfft);

    numffts /= 4;
    fprintf(stderr,"since SIMD implementation does 4 ffts at a time, numffts is being reduced to %d\n",numffts);
#else
    buf=(kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * nfft);
    bufout=(kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * nfft);
#endif

    for (i=0;i<nfft;++i ) {
#ifdef USE_SIMD        
        buf[i].r = _mm_set_ps1((float)( rand() - RAND_MAX/2));
        buf[i].i = _mm_set_ps1((float)( rand() - RAND_MAX/2));
#else        
        buf[i].r = rand() - RAND_MAX/2;
        buf[i].i = rand() - RAND_MAX/2;
#endif        
    }

    pstats_init();

    st = kiss_fft_alloc( nfft ,isinverse ,0,0);

    for (i=0;i<numffts;++i)
        kiss_fft( st ,buf,bufout );

    free(st);

    free(buf); free(bufout);

    fprintf(stderr,"KISS\tnfft=%d\tnumffts=%d\n" ,nfft,numffts);
    pstats_report();

    kiss_fft_cleanup();

    return 0;
}

