/*
Copyright (c) 2003-2004, Mark Borgerding

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    * Neither the author nor the names of any contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "kiss_fftndr.h"
#include "_kiss_fft_guts.h"


// copied from kiss_fftnd.c
struct kiss_fftnd_state{
    int dimprod; 
    int ndims;
    int *dims;
    kiss_fft_cfg *states; /* cfg states for each dimension */
    kiss_fft_cpx * tmpbuf; /*buffer capable of hold the entire input */
};



struct kiss_fftndr_state
{
    int dimReal;
    int dimOther;
    kiss_fftr_cfg cfg_r;
    kiss_fftnd_cfg cfg_nd;
    void * tmpbuf;
};

static int prod(const int *dims, int ndims)
{
    int x=1;
    while (ndims--) 
        x *= *dims++;
    return x;
}

#define CHK fprintf(stderr,"line=%d\t",__LINE__)
kiss_fftndr_cfg kiss_fftndr_alloc(const int *dims,int ndims,int inverse_fft,void*mem,size_t*lenmem)
{
    kiss_fftndr_cfg st = NULL;
    size_t nr=0 , nd=0,ntmp=0;
    int dimReal = dims[0];
    int dimOther = prod(dims+1,ndims-1);
    
    (void)kiss_fftr_alloc(dims[0],inverse_fft,NULL,&nr);
    (void)kiss_fftnd_alloc(dims+1,ndims-1,inverse_fft,NULL,&nd);
    ntmp =
        dimReal * sizeof(kiss_fft_scalar)   // time buffer for one pass
        + (dimReal+2) * sizeof(kiss_fft_scalar)  // freq buffer for one pass
        + dimOther*(dimReal+2) * sizeof(kiss_fft_scalar);  // large enough to hold entire input in case of in-place

    
    size_t memneeded = sizeof( struct kiss_fftndr_state ) + nr + nd + ntmp;

    if (lenmem==NULL) {
        st = (kiss_fftndr_cfg) malloc(memneeded);
    }else{
        if (*lenmem >= memneeded)
            st = (kiss_fftndr_cfg)mem;
        *lenmem = memneeded; 
    }
    if (st==NULL)
        return NULL;
    memset( st , 0 , memneeded);
    
    st->dimReal = dimReal;
    st->dimOther = dimOther;
    st->cfg_r = kiss_fftr_alloc( dims[0],inverse_fft,st+1,&nr);
    st->cfg_nd = kiss_fftnd_alloc(dims+1,ndims-1,inverse_fft, ((char*) st->cfg_r)+nr,&nd);

    // tell the N-D complex FFT to stride the input
    st->cfg_nd->dimprod *= (dims[0]/2+1);
    
    st->tmpbuf = (char*)st->cfg_nd + nd;

    return st;
}

void kiss_fftndr(kiss_fftndr_cfg st,const kiss_fft_scalar *timedata,kiss_fft_cpx *freqdata)
{
    int k1,k2;
    int dimReal = st->dimReal;
    int dimOther = st->dimOther;

    kiss_fft_scalar * t1 = (kiss_fft_scalar*)st->tmpbuf; // enough to hold dimReal scalars
    kiss_fft_cpx * t2 = (kiss_fft_cpx*)(t1+dimReal); // enough to hold the entire input (only used if in==out)

    fprintf(stderr,"dimReal=%d\n",dimReal);
    fprintf(stderr,"dimOther=%d\n",dimOther);
    fprintf(stderr,"t1=%p\n",t1);
    fprintf(stderr,"t2=%p\n",t2);

    
    // take the real data, fft it and transpose
    for (k1=0;k1<dimOther;++k1) {
        for (k2=0;k2<dimReal;++k2)
            t1[k2] = timedata[k1 + k2*dimOther]; // select the appropriate samples into a contiguous buffer

        //then fft the N real samples to create N/2+1 complex points
        kiss_fftr( st->cfg_r, t1, t2 + k1*(dimReal/2+1) );

        fprintf(stderr,"t1=fft([ " );
        for (k2=0;k2<dimReal;++k2) 
            fprintf(stderr,"%f ",t1[k2] );
        fprintf(stderr,"])\n" );
    }

        // fft the remaining dimensions just like the N-D complex case
        kiss_fftnd(st->cfg_nd, t2,freqdata);

    fprintf(stderr,"out=[ " );
    for (k1=0;k1<dimOther;++k1) {
        for (k2=0;k2<dimReal/2+1;++k2) {
            kiss_fft_cpx c = freqdata[k1*(dimReal/2+1) +k2];
            fprintf(stderr,"%f%+fi ",c.r , c.i );
        }
        fprintf(stderr,"\n" );
    }
    fprintf(stderr,"]\n" );
}

void kiss_fftndri(kiss_fftndr_cfg st,const kiss_fft_cpx *freqdata,kiss_fft_scalar *timedata)
{
    memset(timedata,0,sizeof(kiss_fft_scalar) * st->dimReal * st->dimOther);
}
