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
    int minus4; /*magic to signify a N-d transform*/
    int dimprod; /* dimsum would be mighty tasty right now */
    int ndims;
    int *dims;
    void **states;
    kiss_fft_cpx * tmpbuf;
}kiss_fftnd_state;

void * kiss_fftnd_alloc(int *dims,int ndims,int inverse_fft,void*mem,size_t*lenmem)
{
    kiss_fftnd_state *st = NULL;
    int i;
    int dimprod=1;
    size_t memneeded = sizeof(kiss_fftnd_state);
    char * ptr;

    for (i=0;i<ndims;++i) {
        size_t sublen;
        kiss_fft_alloc (dims[i], inverse_fft, NULL, &sublen);
        memneeded += sublen;   /* st->states[i] */
        dimprod *= dims[i];
    }
    memneeded += sizeof(int) * ndims;/*  st->dims */
    memneeded += sizeof(void*) * ndims;/* st->states  */
    memneeded += sizeof(kiss_fft_cpx) * dimprod; /* st->tmpbuf */

    if (lenmem == NULL) {
        st = (kiss_fftnd_state *) malloc (memneeded);
    } else {
        if (*lenmem >= memneeded)
            st = (kiss_fftnd_state *) mem;
        *lenmem = memneeded;
    }
    if (!st)
        return NULL;

    st->minus4 = -4;
    st->dimprod = dimprod;
    st->ndims = ndims;
    ptr=(char*)(st+1);

    st->states = (void**)ptr;
    ptr += sizeof(void*) * ndims;

    st->dims = (int*)ptr;
    ptr += sizeof(int) * ndims;

    st->tmpbuf = (kiss_fft_cpx*)ptr;
    ptr += sizeof(kiss_fft_cpx) * dimprod;

    for (i=0;i<ndims;++i) {
        size_t len;
        st->dims[i] = dims[i];
        kiss_fft_alloc (st->dims[i], inverse_fft, NULL, &len);
        st->states[i] = ptr;
        ptr += len;
        kiss_fft_alloc (st->dims[i], inverse_fft, st->states[i], &len);
    }
    return st;
}

void kiss_fftnd(const void * cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout)
{
    int i,k;
    kiss_fftnd_state *st = ( kiss_fftnd_state *)cfg;
    const kiss_fft_cpx * bufin=fin;
    kiss_fft_cpx * bufout;

    /*arrange it so the last bufout == fout*/
    if ( st->ndims & 1 )
        bufout = fout;
    else
        bufout = st->tmpbuf;

    for ( k=0; k < st->ndims; ++k) {
        int curdim = st->dims[k];
        int stride = st->dimprod / curdim;
        for (i=0;i<stride;++i) 
            kiss_fft_stride( st->states[k], bufin+i , bufout+i*curdim,stride );

        /*toggle back and forth between the two buffers*/
        if (bufout == st->tmpbuf){
            bufout = fout;
            bufin = st->tmpbuf;
        }else{
            bufout = st->tmpbuf;
            bufin = fout;
        }
    }
}
