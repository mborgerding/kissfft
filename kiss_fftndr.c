/*
 *  Copyright (c) 2003-2004, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */

#include "kiss_fftndr.h"
#include "_kiss_fft_guts.h"
#define MAX(x,y) ( ( (x)<(y) )?(y):(x) )

struct kiss_fftndr_state
{
    int dimReal;
    size_t dimOther;
    kiss_fftr_cfg cfg_r;
    kiss_fftnd_cfg cfg_nd;
    void * tmpbuf;
};

kiss_fftndr_cfg kiss_fftndr_alloc(const int *dims,int ndims,int inverse_fft,void*mem,size_t*lenmem)
{
    KISS_FFT_ALIGN_CHECK(mem)

    kiss_fftndr_cfg st = NULL;
    size_t nr=0 , nd=0;
    int dimReal;
    size_t dimOther = 1, nrbins;
    size_t memneeded = 0;
    char * ptr = NULL;
    int k, dimprod = 1;

    if (dims == NULL || ndims < 1)
        goto fail;
    for (k=0; k<ndims; ++k) {
        /* The multidimensional transform uses int-sized indices. */
        if (dims[k] <= 0 || dimprod > INT_MAX / dims[k])
            goto fail;
        dimprod *= dims[k];
        if (k < ndims-1)
            dimOther = (size_t)dimprod;
    }
    dimReal = dims[ndims-1];
    nrbins = (size_t)dimReal / 2 + 1;

    (void)kiss_fftr_alloc(dimReal,inverse_fft,NULL,&nr);
    (void)kiss_fftnd_alloc(dims,ndims-1,inverse_fft,NULL,&nd);
    if (nr == 0 || nr == (size_t)-1 || nd == 0 || nd == (size_t)-1
        || !kf_alloc_size_add(&memneeded, 1, sizeof(struct kiss_fftndr_state))
        || !kf_alloc_size_add(&memneeded, 1, nr)
        || !kf_alloc_size_add(&memneeded, 1, nd)
        /* Frequency buffer for one pass, followed by the in-place input buffer. */
        || !kf_alloc_size_add(&memneeded, MAX(nrbins, dimOther), sizeof(kiss_fft_cpx))
        || !kf_alloc_size_add(&memneeded, dimOther * nrbins, sizeof(kiss_fft_cpx)))
        goto fail;

    if (lenmem==NULL) {
        ptr = (char*) malloc(memneeded);
    }else{
        if (*lenmem >= memneeded)
            ptr = (char *)mem;
        *lenmem = memneeded; 
    }
    if (ptr==NULL)
        return NULL;
    
    st = (kiss_fftndr_cfg) ptr;
    memset( st , 0 , memneeded);
    ptr += KISS_FFT_ALIGN_SIZE_UP(sizeof(struct kiss_fftndr_state));
    
    st->dimReal = dimReal;
    st->dimOther = dimOther;
    st->cfg_r = kiss_fftr_alloc( dimReal,inverse_fft,ptr,&nr);
    ptr += KISS_FFT_ALIGN_SIZE_UP(nr);
    st->cfg_nd = kiss_fftnd_alloc(dims,ndims-1,inverse_fft, ptr,&nd);
    ptr += KISS_FFT_ALIGN_SIZE_UP(nd);
    st->tmpbuf = ptr;

    return st;

fail:
    if (lenmem != NULL)
        *lenmem = (size_t)-1;
    return NULL;
}

void kiss_fftndr(kiss_fftndr_cfg st,const kiss_fft_scalar *timedata,kiss_fft_cpx *freqdata)
{
    size_t k1,k2;
    size_t dimReal = (size_t)st->dimReal;
    size_t dimOther = st->dimOther;
    int nrbins = dimReal/2+1;

    kiss_fft_cpx * tmp1 = (kiss_fft_cpx*)st->tmpbuf; 
    kiss_fft_cpx * tmp2 = tmp1 + MAX(nrbins,dimOther);

    // timedata is N0 x N1 x ... x Nk real

    // take a real chunk of data, fft it and place the output at correct intervals
    for (k1=0;k1<dimOther;++k1) {
        kiss_fftr( st->cfg_r, timedata + k1*dimReal , tmp1 ); // tmp1 now holds nrbins complex points
        for (k2=0;k2<nrbins;++k2)
           tmp2[ k2*dimOther+k1 ] = tmp1[k2];
    }

    for (k2=0;k2<nrbins;++k2) {
        kiss_fftnd(st->cfg_nd, tmp2+k2*dimOther, tmp1);  // tmp1 now holds dimOther complex points
        for (k1=0;k1<dimOther;++k1) 
            freqdata[ k1*(nrbins) + k2] = tmp1[k1];
    }
}

void kiss_fftndri(kiss_fftndr_cfg st,const kiss_fft_cpx *freqdata,kiss_fft_scalar *timedata)
{
    int k1,k2;
    int dimReal = st->dimReal;
    int dimOther = st->dimOther;
    int nrbins = dimReal/2+1;
    kiss_fft_cpx * tmp1 = (kiss_fft_cpx*)st->tmpbuf; 
    kiss_fft_cpx * tmp2 = tmp1 + MAX(nrbins,dimOther);

    for (k2=0;k2<nrbins;++k2) {
        for (k1=0;k1<dimOther;++k1) 
            tmp1[k1] = freqdata[ k1*(nrbins) + k2 ];
        kiss_fftnd(st->cfg_nd, tmp1, tmp2+k2*dimOther);
    }

    for (k1=0;k1<dimOther;++k1) {
        for (k2=0;k2<nrbins;++k2)
            tmp1[k2] = tmp2[ k2*dimOther+k1 ];
        kiss_fftri( st->cfg_r,tmp1,timedata + k1*dimReal);
    }
}
