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
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "kiss_fft.h"
/*
 * kiss_fft.h
 * defines kiss_fft_scalar as either short or a float type
 * and defines
 * typedef struct {
 *     kiss_fft_scalar r;
 *     kiss_fft_scalar i;
 * }kiss_fft_cpx;
 */

typedef struct {
    int nfft;
    int inverse;
    int *factors;
    kiss_fft_cpx * twiddles;
    kiss_fft_cpx * tmpbuf;
    kiss_fft_cpx * scratch;
}kiss_fft_state;

#ifdef FIXED_POINT
#define  C_SUB( res, a,b)\
    do {    (res).r=((a).r-b.r+1)>>1;  (res).i=((a).i-b.i+1)>>1;  }while(0)
#define C_ADDTO( res , a)\
    do {    (res).r=((res).r+(a).r+1)>>1;  (res).i=((res).i+(a).i+1)>>1;  }while(0)
    /*  We don't have to worry about overflow from multiplying by twiddle factors since they
     *  all have unity magnitude.  Still need to shift away fractional bits after adding 1/2 for
     *  rounding. */
#   define C_MUL(m,a,b) \
      do{ (m).r = ( ( (a).r*(b).r - (a).i*(b).i)  + (1<<14) ) >> 15;\
          (m).i = ( ( (a).r*(b).i + (a).i*(b).r)  + (1<<14) ) >> 15;\
      }while(0)
#else  // not FIXED_POINT

#define  C_SUB( res, a,b)\
    do {    (res).r=(a).r-(b).r;  (res).i=(a).i-(b).i;  }while(0)
#define C_ADDTO( res , a)\
    do {    (res).r += (a).r;  (res).i += (a).i;  }while(0)
#define C_MUL(m,a,b) \
    do{ (m).r = (a).r*(b).r - (a).i*(b).i;\
        (m).i = (a).r*(b).i + (a).i*(b).r; }while(0)
#endif

static
kiss_fft_cpx cexp(double phase)
{
    kiss_fft_cpx x;
#ifdef FIXED_POINT    
    x.r = (kiss_fft_scalar) ( 32767*cos(phase) );
    x.i = (kiss_fft_scalar) ( -32767*sin(phase) );
#else
    x.r = cos(phase);
    x.i = -sin(phase);
#endif
    return x;
}

static 
void fft_work(
        kiss_fft_cpx * Fout,
        const kiss_fft_cpx * f,
        int fstride,
        int * factors,
        const kiss_fft_state * st
        );

// the heart of the fft
static 
void fft2work(
        kiss_fft_cpx * Fout,
        const kiss_fft_cpx * f,
        int fstride,
        int * factors,
        const kiss_fft_state * st
        )
{
    kiss_fft_cpx * Fout2;
    int u,m;
    kiss_fft_cpx t;
    kiss_fft_cpx * twiddles = st->twiddles;

    factors++;// p==2
    m=*factors++;

    if (*factors == 2){
        fft2work(Fout,f,fstride*2,factors,st);
        fft2work(Fout+m,f+fstride,fstride*2,factors,st);
    }else if (m==1) {
        Fout[0] = f[0];
        Fout[1] = f[fstride];
    } else {
        fft_work( Fout , f, fstride*2,factors,st);
        fft_work( Fout + m, f+fstride, fstride*2,factors,st);
    }

    Fout2 = Fout + m;
    for ( u=0; u<m; ++u ) {
        C_MUL (t,  *Fout2 , *twiddles);
        twiddles += fstride;
        C_SUB( *Fout2 ,  *Fout , t );
        C_ADDTO( *Fout ,  t );
        ++Fout2;
        ++Fout;
    }
}

static 
void fft_work(
        kiss_fft_cpx * Fout,
        const kiss_fft_cpx * f,
        int fstride,
        int * factors,
        const kiss_fft_state * st
        )
{
    int m,p,q,q1,u,k;
    kiss_fft_cpx t;
    kiss_fft_cpx * scratch = st->scratch;
    kiss_fft_cpx * twiddles = st->twiddles;
    int Norig = st->nfft;

#if 1    
    switch (*factors) {
        case 2:
            fft2work(Fout,f,fstride,factors,st);
            return;
        default:
            break;
    }
#endif
    p=*factors++;
    m=*factors++;

    for (q=0;q<p;++q) {
        if (m==1) 
            Fout[q] = *f;
        else
            fft_work( Fout + m*q, f, fstride*p,factors,st);
        f+= fstride;
    }

    for ( u=0; u<m; ++u ) {
        k=u;
        for ( q1=0 ; q1<p ; ++q1 ) {
            scratch[q1] = Fout[ k  ];
#ifdef FIXED_POINT
            scratch[q1].r >>= 1;
            scratch[q1].i >>= 1;
#endif            
            k += m;
        }

        k=u;
        for ( q1=0 ; q1<p ; ++q1 ) {
            int twidx=0;
            Fout[ k ] = scratch[0];
            for (q=1;q<p;++q ) {
                twidx += fstride * k;
                if (twidx>=Norig)
                    twidx-=Norig;
                C_MUL(t,scratch[q] , twiddles[twidx] );
                Fout[ k ].r += t.r;
                Fout[ k ].i += t.i;
            }
            k += m;
        }
    }
}

/*
 *      void * kiss_fft_alloc(int nfft,int inverse_fft)
 *
 * User-callable function to allocate all necessary scratch space for the fft.
 *
 * The return value is a contiguous block of memory, allocated with malloc.  As such,
 * It can be freed with free(), rather than a kiss_fft-specific function.
 * */
void * kiss_fft_alloc(int nfft,int inverse_fft)
{
    int allocsize;
    int nstages=0;
    int i;
    kiss_fft_state * st=NULL;

    allocsize =  sizeof(kiss_fft_state)
        + sizeof(kiss_fft_cpx)*nfft // twiddle factors
        + sizeof(kiss_fft_cpx)*nfft // tmpbuf
        + sizeof(int)*nfft // factors
        + sizeof(kiss_fft_cpx)*nfft; // scratch
    
    st = ( kiss_fft_state *)malloc( allocsize );
    if (!st)
        return NULL;

    st->nfft=nfft;
    st->inverse = inverse_fft;
    st->twiddles = (kiss_fft_cpx*)(st+1); // just beyond struct
    st->tmpbuf = (kiss_fft_cpx*)(st->twiddles + nfft);//  just after twiddles
    st->scratch = (kiss_fft_cpx*)(st->tmpbuf + nfft);
    st->factors = (int*)(st->scratch + nfft); // just after tmpbuf


    for (i=0;i<nfft;++i) {
        const double pi=3.14159265358979323846264338327;
        double phase = ( 2*pi /nfft ) * i;
        if (st->inverse)
            phase *= -1;
        st->twiddles[i] = cexp( phase );
    }

    while (nfft>1) {
        const int primes[] = {2,3,5,7,11,13,17,-1};
        int p=nfft;
        i=0;
        while ( primes[i] != -1 ) {
            if ( nfft %  primes[i] == 0){
                p =  primes[i];
                break;
            }
            ++i;
        }
        st->factors[2*nstages] = p;
        nfft /= p;
        st->factors[2*nstages+1] = nfft;
        ++nstages;
    }

    return st;
}

void kiss_fft(const void * cfg,kiss_fft_cpx *f)
{
    int n;
    const kiss_fft_state * st = cfg;

    n = st->nfft;
    memcpy(st->tmpbuf,f,sizeof(kiss_fft_cpx)*n);
    fft_work( f, st->tmpbuf, 1, st->factors,st );
}
