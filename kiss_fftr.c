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
    int minus3; /*magic to signify a 1-d real transform*/
    kiss_fft_state * substate;
    kiss_fft_cpx * tmpbuf;
    kiss_fft_cpx * super_twiddles;
}kiss_fftr_state;

void * kiss_fftr_alloc(int nfft,int inverse_fft)
{
    int i;
    kiss_fftr_state *st = NULL;
    int subsize;
    
    if (nfft&1) {
        /*fprintf(stderr,"Real FFT optimization must be even.\n"); */
        return NULL;
    }
    nfft >>= 1;

    subsize = kf_allocsize(nfft);

    st = (kiss_fftr_state *) malloc ( sizeof(kiss_fftr_state)
            + subsize + sizeof(kiss_fft_cpx) * ( nfft * 2) );
    if (!st)
        return NULL;

    st->minus3 = -3;
    st->substate = (kiss_fft_state *)(st+1); /*just beyond kiss_fftr_state struct */
    st->tmpbuf =  (kiss_fft_cpx*)(((char*)st->substate) + subsize);
    st->super_twiddles =  st->tmpbuf + nfft;
    kf_init_state (st->substate, nfft, inverse_fft);

    for (i=0;i<nfft;++i) {
        double phase = -3.14159265358979323846264338327 * i / nfft;
        st->super_twiddles[i] = kf_cexp( phase );
    }
    return st;
}

static void pcpx( kiss_fft_cpx * c)
{
    printf("%g + %gi\n",c->r,c->i);
}

void kiss_fftr(const void * cfg,const kiss_fft_scalar *fin,kiss_fft_cpx *fout)
{
    /* input buffer fin is stored row-wise */
    kiss_fftr_state *st = ( kiss_fftr_state *)cfg;
    int k,N;

    if (st->minus3 != -3 || st->substate->inverse ) {
        fprintf(stderr,"kiss fft usage error: improper alloc\n");
        exit(1);
    }

    N = st->substate->nfft;

    /*perform the parallel fft of two real signals packed in real,imag*/
    kiss_fft( st->substate , (const kiss_fft_cpx*)fin, st->tmpbuf );
 
    for (k=0;k<N;++k) {
        kiss_fft_cpx fpnk,fpk,f1k,f2k,scr[3];
        
        fpk = st->tmpbuf[k];
        /* fpk is Fp[k] */

        if (k==0)
            fpnk = st->tmpbuf[0];
        else
            fpnk = st->tmpbuf[N-k];
        fpnk.i = fpnk.i * -1;
        /* fpnk is Fp[-k].conjugate() */

        C_ADD( f1k, fpk , fpnk );
        /*printf("F1k[%d]:",k); pcpx(&f1k); */

        C_SUB( scr[0], fpk , fpnk );
        f2k.r = scr[0].i;
        f2k.i = -scr[0].r;
        /*
        printf("F2k[%d]:",k); pcpx(&f2k);
        */

        C_MUL( scr[1], f2k , st->super_twiddles[k] );

        C_ADD(fout[k],f1k,scr[1]);
        fout[k].r /= 2;
        fout[k].i /= 2;
    }
}

