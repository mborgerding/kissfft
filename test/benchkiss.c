#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/times.h>
#include <unistd.h>
#include "kiss_fft.h"

int main(int argc,char ** argv)
{
    int nfft=1024;
    int isinverse=0;
    int numffts=1000,i;
    struct tms t0,t1;
    float cputime;
    kiss_fft_cpx * buf;
    kiss_fft_cpx * bufout;
    void *st;

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
    buf=(kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * nfft);
    bufout=(kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * nfft);

    for (i=0;i<nfft;++i ) {
        buf[i].r = rand() - RAND_MAX/2;
        buf[i].i = rand() - RAND_MAX/2;
    }

    times(&t0);
    st = kiss_fft_alloc( nfft ,isinverse );

    for (i=0;i<numffts;++i)
        kiss_fft_io( st ,buf,bufout );

    free(st);
    times(&t1);

    free(buf); free(bufout);

    cputime = ( ((float)t1.tms_utime + t1.tms_stime + t1.tms_cutime + t1.tms_cstime ) - 
                ((float)t0.tms_utime + t0.tms_stime + t0.tms_cutime + t0.tms_cstime ) )
             / sysconf(_SC_CLK_TCK);

    fprintf(stderr,"KISS\t" 
            "nfft=%d\t"
            "numffts=%d\t"
            "cputime=%.3f\n" , 
            nfft,numffts,cputime
            );
    return 0;
}

