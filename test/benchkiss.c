#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <unistd.h>
#include "kiss_fft.h"

#include "pstats.h"

#define CHK fprintf(stderr,"line %d\n" , __LINE__ )

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
        break;
      case 'x':
        numffts = atoi (optarg);
        break;
      case 'i':
        isinverse = 1;
        break;
      }
    }
    CHK;
#ifdef USE_SIMD        
    buf=(kiss_fft_cpx*)memalign(sizeof(kiss_fft_cpx),sizeof(kiss_fft_cpx) * nfft);
    bufout=(kiss_fft_cpx*)memalign(sizeof(kiss_fft_cpx),sizeof(kiss_fft_cpx) * nfft);

    numffts /= 4;
    fprintf(stderr,"since SIMD implementation does 4 ffts at a time, numffts is being reduced to %d\n",numffts);
#else
    buf=(kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * nfft);
    bufout=(kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * nfft);
#endif
    
    fprintf(stderr,"buf at %p, bufout at %p\n",buf,bufout);
    CHK;

    for (i=0;i<nfft;++i ) {
#ifdef USE_SIMD        
        buf[i].r = _mm_set_ps1((float)( rand() - RAND_MAX/2));
        buf[i].i = _mm_set_ps1((float)( rand() - RAND_MAX/2));
#else        
        buf[i].r = rand() - RAND_MAX/2;
        buf[i].i = rand() - RAND_MAX/2;
#endif        
    }
    CHK;

    pstats_init();
    CHK;

    st = kiss_fft_alloc( nfft ,isinverse ,0,0);
    CHK;

    for (i=0;i<numffts;++i)
        kiss_fft( st ,buf,bufout );
    CHK;

    free(st);

    free(buf); free(bufout);

    fprintf(stderr,"KISS\tnfft=%d\tnumffts=%d\n" ,nfft,numffts);
    pstats_report();

    return 0;
}

