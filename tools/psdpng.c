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
#include <unistd.h>
#include <png.h>

#include "kiss_fft.h"
#include "kiss_fftr.h"

int nfft=1024;
int nfreqs=0;
int repeat=0;
int window=0;
int colors=256;
FILE * fin=NULL;
FILE * fout=NULL;

int navg=20;
float * valbuf=NULL;
int remove_dc=1;

size_t nrows=0;
float * vals=NULL;

void config(int argc,char** argv)
{
    while (1) {
        int c = getopt (argc, argv, "n:rs");
        if (c == -1)
            break;
        switch (c) {
        case 'n':
            nfft=atoi(optarg);
        case '?':
            fprintf (stderr, "usage options:\n"
                     "\t-n d: fft dimension(s) default = 1024\n"
                     "stereo 16 bit machine format real input is assumed\n"
                     );
        default:
            fprintf (stderr, "bad %c\n", c);
            exit (1);
            break;
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
    if (fin==NULL)
        fin=stdin;
    if (fout==NULL)
        fout=stdout;
}

#define CHECKNULL(p) if ( (p)==NULL ) do { fprintf(stderr,"CHECKNULL failed @ %s(%d): %s\n",__FILE__,__LINE__,#p );exit(1);} while(0)

typedef struct
{
    png_byte r;
    png_byte g;
    png_byte b;
} rgb_t;


rgb_t val2rgb(float x)
{
    const double pi = 3.14159265358979;
    rgb_t p;

    p.g = (int)(255*sin(x*pi));
    p.r = (int)(255*abs(sin(x*pi*3/2)));
    p.b = (int)(255*abs(sin(x*pi*5/2)));
    return p;
}


void cpx2pixels(rgb_t * res,const float * fbuf,size_t nfreqs)
{
    int i;
    float minval,maxval,valrange;
    minval=maxval=fbuf[0];

    for (i = 0; i < nfreqs; ++i) {
        if (fbuf[i] > maxval) maxval = fbuf[i];
        if (fbuf[i] < minval) minval = fbuf[i];
    }

    valrange = maxval-minval;

    for (i = 0; i < nfreqs; ++i)
        res[i] = val2rgb( (fbuf[i] - minval)/valrange );
}

void transform_signal()
{
    short *inbuf;
    void * cfg=NULL;
    kiss_fft_scalar *tbuf;
    kiss_fft_cpx *fbuf;
    kiss_fft_cpx *avgbuf;
    size_t i;
    int n;
    int avgctr=0;

    nfreqs=nfft/2+1;

    CHECKNULL( cfg=kiss_fftr_alloc(nfft,0,0,0) );
    CHECKNULL( inbuf=(short*)malloc(sizeof(short)*2*nfft ) );
    CHECKNULL( tbuf=(kiss_fft_scalar*)malloc(sizeof(kiss_fft_scalar)*nfft ) );
    CHECKNULL( fbuf=(kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx)*nfreqs ) );
    CHECKNULL( avgbuf=(kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx)*nfreqs ) );

    while ( ( n = fread(inbuf,sizeof(short)*2,nfft,fin) ) == nfft ) {
        /* pack the shorts */
        for (i=0;i<nfft;++i){
            tbuf[i] = inbuf[2*i] + inbuf[2*i+1];
        }

        /* do FFT */
        kiss_fftr(cfg,tbuf,fbuf);

        if (remove_dc) {
            fbuf[0].r = 0;
            fbuf[0].i = 0;
        }

        for (i=0;i<nfreqs;++i){
            avgbuf[i].r += fbuf[i].r;
            avgbuf[i].i += fbuf[i].i;
        }

        if (++avgctr == navg) {
            avgctr=0;
            ++nrows;
            vals = (float*)realloc(vals,sizeof(float)*nrows*nfreqs);
            
            for (i=0;i<nfreqs;++i) {
                float binpower = fbuf[i].r * fbuf[i].r + fbuf[i].i * fbuf[i].i;
                vals[(nrows - 1) * nfreqs + i] = 10 * log10 (binpower);
            }

            memset(avgbuf,0,sizeof(avgbuf[0])*nfreqs);
        }
    }

    free(cfg);
    free(inbuf);
    free(tbuf);
    free(fbuf);
    free(avgbuf);
}

void make_png()
{
    png_bytepp row_pointers=NULL;
    rgb_t * row_data=NULL;
    size_t i;

    png_structp png_ptr=NULL;
    png_infop info_ptr=NULL;
    
    CHECKNULL( png_ptr = png_create_write_struct (PNG_LIBPNG_VER_STRING,0,0,0) );
    CHECKNULL( info_ptr = png_create_info_struct(png_ptr) );


    png_init_io(png_ptr, fout );
    png_set_IHDR(png_ptr, info_ptr ,nfreqs,nrows,8,PNG_COLOR_TYPE_RGB,PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT );
    

    row_data = (rgb_t*)malloc(sizeof(rgb_t) * nrows * nfreqs) ;
    cpx2pixels(row_data, vals, nfreqs*nrows );

    row_pointers = realloc(row_pointers, nrows*sizeof(png_bytep));
    for (i=0;i<nrows;++i) {
        row_pointers[i] = (png_bytep)(row_data + i*nfreqs);
    }
    png_set_rows(png_ptr, info_ptr, row_pointers);


    fprintf(stderr,"creating %dx%d png\n",nfreqs,nrows);
    fprintf(stderr,"bitdepth %d \n",png_get_bit_depth(png_ptr,info_ptr ) );

    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY , NULL);

}

int main(int argc,char ** argv)
{
    config(argc,argv);

    transform_signal();

    make_png();

    if (fout!=stdout) fclose(fout);
    if (fin!=stdin) fclose(fin);
    return 0;
}
