#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <getopt.h>
#include "pstats.h"

int main(int argc,char ** argv)
{
    int nfft=1024;
    int isinverse=0;
    int numffts=1000,i;

    fftw_complex * in=NULL;
    fftw_complex * out=NULL;
    fftw_plan p;

    pstats_init();

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

    if ( isinverse )
        p = fftw_plan_dft_1d(nfft, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    else    
        p = fftw_plan_dft_1d(nfft, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (i=0;i<numffts;++i)
        fftw_execute(p);

    fftw_destroy_plan(p);

    fftw_free(in); fftw_free(out);

    fprintf(stderr,"fftw\tnfft=%d\tnumffts=%d\n", nfft,numffts);
    pstats_report();

    return 0;
}

