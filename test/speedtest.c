#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <getopt.h>


int main(int argc,char ** argv)
{
    int nfft=1024;
    int isinverse=0;
    FILE *fin=stdin;
    FILE *fout=stdout;
    int times=1,i;

    fftw_complex * in=NULL;
    fftw_complex * out=NULL;
    fftw_plan p;

    while (1) {
      int c = getopt (argc, argv, "n:ix:");
      if (c == -1)
        break;
      switch (c) {
      case 'n':
        nfft = atoi (optarg);
        break;
      case 'x':
        times = atoi (optarg);
        break;
      case 'i':
        isinverse = 1;
        break;
      }
    }
    //fprintf(stderr,"sizeof(fftw_complex) = %d \n" , sizeof(fftw_complex) );
    fprintf(stderr,"sizeof(fftw_complex[0]) = %d \n" , sizeof((*in)[0]) );

    in=fftw_malloc(sizeof(fftw_complex) * nfft);
    out=fftw_malloc(sizeof(fftw_complex) * nfft);
    if ( isinverse )
        p = fftw_plan_dft_1d(nfft, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    else    
        p = fftw_plan_dft_1d(nfft, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    while ( fread( in , sizeof(fftw_complex) , nfft , fin ) > 0 ) {
        for (i=0;i<times;++i)
            fftw_execute(p);
        fwrite( out , sizeof(fftw_complex) , nfft , fout );
    }

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);

    return 0;
}

