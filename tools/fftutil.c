/*
Copyright (c) 2003, Mark Borgerding

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    * Neither the author nor the names of any contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>

#include "kiss_fft.h"

#ifndef NFFT
# define NFFT 1024
#endif



void fft_file(FILE * fin,FILE * fout,int nfft,int isinverse,int useascii,int times)
{
    int i;
    void *st;
    kiss_fft_cpx * buf;
    kiss_fft_cpx * bufout;

    buf = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * nfft );
    bufout = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx) * nfft );
    st = kiss_fft_alloc( nfft ,isinverse );

    while ( fread( buf , sizeof(kiss_fft_cpx) * nfft ,1, fin ) > 0 ) {
        for (i=0;i<times;++i)
            kiss_fft_io( st , buf ,bufout);
        if (useascii) {
            int i;
            for (i=0;i<nfft;++i) 
                fprintf(fout, "(%g,%g) ", (double)bufout[i].r,(double)bufout[i].i);
        }else{            
            fwrite( bufout , sizeof(kiss_fft_cpx) , nfft , fout );
        }
    }
    free(st);
    free(buf);
    free(bufout);
}

int main(int argc,char ** argv)
{
    int nfft=1024;
    int isinverse=0;
    FILE *fin=stdin;
    FILE *fout=stdout;
    int useascii=0;
    int times=1;

    while (1) {
        int c=getopt(argc,argv,"n:iax:");
        if (c==-1) break;
        switch (c) {
            case 'a':useascii=1;break;
            case 'n':nfft = atoi(optarg);break;
            case 'i':isinverse=1;break;
            case 'x':times=atoi(optarg);break;
        }
    }

    if ( optind < argc ) {
        if (strcmp("-",argv[optind]) !=0)
            fin = fopen(argv[optind],"rb");
        ++optind;
    }

    if ( optind < argc ) {
        if ( strcmp("-",argv[optind]) !=0 ) 
            fout = fopen(argv[optind],"wb");
        ++optind;
    }

    fft_file(fin,fout,nfft,isinverse,useascii,times);

    if (fout!=stdout) fclose(fout);
    if (fin!=stdin) fclose(fin);

    return 0;
}
