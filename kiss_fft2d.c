/*
Copyright (c) 2003, Mark Borgerding

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    * Neither the author nor the names of any contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "_kiss_fft_guts.h"

typedef struct {
    int minus2; /*magic to signify a 2-d transform*/
    kiss_fft_state * rowst;
    kiss_fft_state * colst;
    kiss_fft_cpx * tmpbuf;
}kiss_fft2d_state;

void * kiss_fft2d_alloc(int nrows,int ncols,int inverse_fft)
{
    kiss_fft2d_state *st = NULL;
    int size1,size2,sizetmp;
    size1 = kf_allocsize(ncols);
    size2 = kf_allocsize(nrows);
    sizetmp = sizeof(kiss_fft_cpx)*(ncols > nrows ? ncols : nrows);

    st = (kiss_fft2d_state *) malloc ( sizeof(kiss_fft2d_state) + size1 + size2 + sizetmp );
    if (!st)
        return NULL;

    st->minus2 = -2;
    st->rowst = (kiss_fft_state *)(st+1); /*just beyond kiss_fft2d_state struct */
    st->colst = (kiss_fft_state *)( (char*)(st->rowst) + size1 );
    st->tmpbuf = (kiss_fft_cpx *)( (char*)(st->rowst) + size1 + size2 );
    kf_init_state (st->rowst, ncols, inverse_fft);
    kf_init_state (st->colst, nrows, inverse_fft);
    return st;
}

void kiss_fft2d(const void * cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout)
{
    /* input buffer fin is stored row-wise */
    kiss_fft2d_state *st = ( kiss_fft2d_state *)cfg;
    int row,col;
    int nrows,ncols;
    nrows = st->colst->nfft;
    ncols = st->rowst->nfft;

    /*fft each column*/
    for (col=0;col<ncols;++col) {
        for (row=0;row< nrows ;++row)
            st->tmpbuf[row] = fin[row*ncols + col];
        kiss_fft(st->colst,st->tmpbuf,st->tmpbuf);
        for (row=0;row< nrows ;++row) {
            fout[row*ncols + col] = st->tmpbuf[row];
        }
    }

    /*fft each row */
    for (row=0;row< nrows ;++row)
        kiss_fft(st->rowst , fout + row*ncols , fout + row*ncols );
}

