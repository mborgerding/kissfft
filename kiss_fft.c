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

/* kiss_fft.h
   defines kiss_fft_scalar as either short or a float type
   and defines
   typedef struct { kiss_fft_scalar r; kiss_fft_scalar i; }kiss_fft_cpx; */
#include "kiss_fft.h"

typedef struct {
    int nfft;
    int inverse;
    int *factors;
    kiss_fft_cpx * twiddles;
    kiss_fft_cpx * tmpbuf;
    kiss_fft_cpx * scratch;
}kiss_fft_state;

/*
  Explanation of macros dealing with complex math:

   C_MUL(m,a,b)         : m = a*b
   C_FIXDIV( c , div )  : if a fixed point impl., c /= div. noop otherwise
   C_SUB( res, a,b)     : res = a - b 
   C_SUBFROM( res , a)  : res -= a
   C_ADDTO( res , a)    : res += a
   C_ROTADDTO(sum,c,q)  : sum += c * exp(-j*q*pi/4)
 * */
#ifdef FIXED_POINT

#   define C_MUL(m,a,b) \
      do{ (m).r = ( ( (a).r*(b).r - (a).i*(b).i)  + (1<<14) ) >> 15;\
          (m).i = ( ( (a).r*(b).i + (a).i*(b).r)  + (1<<14) ) >> 15;\
      }while(0)
#   define C_FIXDIV(c,div) \
    do{ (c).r /= div; (c).i /=div; }while(0)

#else  /* not FIXED_POINT*/

#define C_MUL(m,a,b) \
    do{ (m).r = (a).r*(b).r - (a).i*(b).i;\
        (m).i = (a).r*(b).i + (a).i*(b).r; }while(0)
#   define C_FIXDIV(c,div) /* NOOP */
#endif

#define  C_SUB( res, a,b)\
    do {    (res).r=(a).r-(b).r;  (res).i=(a).i-(b).i;  }while(0)
#define C_ADDTO( res , a)\
    do {    (res).r += (a).r;  (res).i += (a).i;  }while(0)
#define C_SUBFROM( res , a)\
    do {    (res).r -= (a).r;  (res).i -= (a).i;  }while(0)
#define C_ROTADDTO(sum,c,q) \
    do{ switch (q) {\
            case 0: C_ADDTO(sum,c); break;\
            case 1: (sum).r += (c).i; (sum).i -= (c).r; break;\
            case 2: C_SUBFROM(sum,c); break;\
            case 3: (sum).r -= (c).i; (sum).i += (c).r; break;\
        } }while(0) 

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

/* bfly2 is a optimization of bfly_generic for p==2 */
void bfly2(
        kiss_fft_cpx * Fout,
        int fstride,
        const kiss_fft_state * st,
        int m
        )
{
    kiss_fft_cpx * Fout2;
    kiss_fft_cpx * tw1 = st->twiddles;
    kiss_fft_cpx t;
    Fout2 = Fout + m;
    do{
        C_MUL (t,  *Fout2 , *tw1);
        tw1 += fstride;
        C_FIXDIV(*Fout,2); C_FIXDIV(t,2);
        C_SUB( *Fout2 ,  *Fout , t );
        C_ADDTO( *Fout ,  t );
        ++Fout2;
        ++Fout;
    }while (--m);
}

/* bfly4 is a optimization of bfly_generic for p==4 */
void bfly4(
        kiss_fft_cpx * Fout,
        int fstride,
        const kiss_fft_state * st,
        int m
        )
{
    kiss_fft_cpx *Fout1,*Fout2,*Fout3;
    kiss_fft_cpx t1,t2,t3;
    kiss_fft_cpx *tw1,*tw2,*tw3;

    Fout1 = Fout + m;
    Fout2 = Fout + 2*m;
    Fout3 = Fout + 3*m;
    tw3 = tw2 = tw1 = st->twiddles;

    do {
        C_FIXDIV(*Fout,4); C_FIXDIV(*Fout1,4); C_FIXDIV(*Fout2,4); C_FIXDIV(*Fout3,4);

        C_MUL(t1,*Fout1 , *tw1 );
        tw1 += fstride;
        C_MUL(t2,*Fout2 , *tw2 );
        tw2 += fstride*2;
        C_MUL(t3,*Fout3 , *tw3 );
        tw3 += fstride*3;

        *Fout3 = *Fout2 = *Fout1 = *Fout;

        C_ADDTO(*Fout,t1);
        C_ADDTO(*Fout,t2);
        C_ADDTO(*Fout,t3);
        C_SUBFROM(*Fout2,t3);
        C_SUBFROM(*Fout2,t1); 
        C_ADDTO(  *Fout2,t2);
        C_SUBFROM(*Fout3,t2);
        C_SUBFROM(*Fout1,t2);

        if(st->inverse) {
            C_ROTADDTO(*Fout1,t1,3);
            C_ROTADDTO(*Fout1,t3,1);
            C_ROTADDTO(*Fout3,t1,1);
            C_ROTADDTO(*Fout3,t3,3);
        }else{
            C_ROTADDTO(*Fout1,t1,1);
            C_ROTADDTO(*Fout1,t3,3);
            C_ROTADDTO(*Fout3,t1,3);
            C_ROTADDTO(*Fout3,t3,1);
        }
        ++Fout; ++Fout1; ++Fout2; ++Fout3;
    }while(--m);
}

/* perform the butterfly for one stage of a mixed radix FFT */
void bfly_generic(
        kiss_fft_cpx * Fout,
        int fstride,
        const kiss_fft_state * st,
        int m,
        int p
        )
{
    int u,k,q1,q;
    kiss_fft_cpx * scratch = st->scratch;
    kiss_fft_cpx * twiddles = st->twiddles;
    kiss_fft_cpx t;
    int Norig = st->nfft;

    for ( u=0; u<m; ++u ) {
        k=u;
        for ( q1=0 ; q1<p ; ++q1 ) {
            scratch[q1] = Fout[ k  ];
            C_FIXDIV(scratch[q1],p);
            k += m;
        }

        k=u;
        for ( q1=0 ; q1<p ; ++q1 ) {
            int twidx=0;
            Fout[ k ] = scratch[0];
            for (q=1;q<p;++q ) {
                twidx += fstride * k;
                if (twidx>=Norig) twidx-=Norig;
                C_MUL(t,scratch[q] , twiddles[twidx] );
                C_ADDTO( Fout[ k ] ,t);
            }
            k += m;
        }
    }
}

void fft_work(
        kiss_fft_cpx * Fout,
        const kiss_fft_cpx * f,
        int fstride,
        int * factors,
        const kiss_fft_state * st
        )
{
    int m,p,q;
    p=*factors++;
    m=*factors++;

    for (q=0;q<p;++q) {
        if (m==1) 
            Fout[q] = *f;
        else
            fft_work( Fout + m*q, f, fstride*p,factors,st);
        f += fstride;
    }

    switch (p) {
        case 4: bfly4(Fout,fstride,st,m); break;
        case 2: bfly2(Fout,fstride,st,m); break;
        default: bfly_generic(Fout,fstride,st,m,p); break;
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
        + sizeof(kiss_fft_cpx)*nfft /* twiddle factors*/
        + sizeof(kiss_fft_cpx)*nfft /* tmpbuf*/
        + sizeof(int)*nfft /* factors*/
        + sizeof(kiss_fft_cpx)*nfft; /* scratch*/
    
    st = ( kiss_fft_state *)malloc( allocsize );
    if (!st)
        return NULL;

    st->nfft=nfft;
    st->inverse = inverse_fft;
    st->twiddles = (kiss_fft_cpx*)(st+1); /* just beyond struct*/
    st->tmpbuf = (kiss_fft_cpx*)(st->twiddles + nfft);/*  just after twiddles*/
    st->scratch = (kiss_fft_cpx*)(st->tmpbuf + nfft);
    st->factors = (int*)(st->scratch + nfft); /* just after tmpbuf*/


    for (i=0;i<nfft;++i) {
        const double pi=3.14159265358979323846264338327;
        double phase = ( 2*pi /nfft ) * i;
        if (st->inverse)
            phase *= -1;
        st->twiddles[i] = cexp( phase );
    }

    while (nfft>1) {
        /* If you add a new radix, don't forget to put it here */
        const int primes[] = {4,2,-1};
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

/* original form of processing function, first release of KISS FFT was in-place.  This maintains API. */
void kiss_fft(const void * cfg,kiss_fft_cpx *f)
{
    const kiss_fft_state * st = cfg;
    memcpy(st->tmpbuf,f,sizeof(kiss_fft_cpx)*st->nfft);
    fft_work( f, st->tmpbuf, 1, st->factors,st );
}

/* two buffer version of above */
void kiss_fft_io(const void * cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout)
{
    const kiss_fft_state * st = cfg;
    fft_work( fout, fin, 1, st->factors,st );
}
