#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <getopt.h>
#include <sys/times.h>
#include <unistd.h>

int main(int argc,char ** argv)
{
    int nfft=1024;
    int isinverse=0;
    int numffts=1000,i;

    fftw_complex * in=NULL;
    fftw_complex * out=NULL;
    fftw_plan p;
    struct tms t0,t1;
    float cputime;

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

    in=fftw_malloc(sizeof(fftw_complex) * nfft);
    out=fftw_malloc(sizeof(fftw_complex) * nfft);
    for (i=0;i<nfft;++i ) {
        in[i][0] = rand() - RAND_MAX/2;
        in[i][1] = rand() - RAND_MAX/2;
    }

    times(&t0);
    if ( isinverse )
        p = fftw_plan_dft_1d(nfft, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    else    
        p = fftw_plan_dft_1d(nfft, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (i=0;i<numffts;++i)
        fftw_execute(p);

    fftw_destroy_plan(p);
    times(&t1);

    fftw_free(in); fftw_free(out);

    cputime = ( ((float)t1.tms_utime + t1.tms_stime + t1.tms_cutime + t1.tms_cstime ) - 
                ((float)t0.tms_utime + t0.tms_stime + t0.tms_cutime + t0.tms_cstime ) )
             / sysconf(_SC_CLK_TCK);

    fprintf(stderr,"fftw\t" 
            "nfft=%d\t"
            "numffts=%d\t"
            "cputime=%.3f\n" , 
            nfft,numffts,cputime
            );
    return 0;
}

