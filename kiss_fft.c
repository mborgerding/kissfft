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

typedef struct {
    int minus2; /*signify a 2-d transform*/
    kiss_fft_state * rowst;
    kiss_fft_state * colst;
}kiss_fft2d_state;


/*
  Explanation of macros dealing with complex math:

   C_MUL(m,a,b)         : m = a*b
   C_FIXDIV( c , div )  : if a fixed point impl., c /= div. noop otherwise
   C_SUB( res, a,b)     : res = a - b 
   C_SUBFROM( res , a)  : res -= a
   C_ADDTO( res , a)    : res += a
 * */
#ifdef FIXED_POINT
#   define S_MUL(a,b) ( ( (a)*(b) + (1<<14) )>>15 )
#   define C_MUL(m,a,b) \
      do{ (m).r = ( ( (a).r*(b).r - (a).i*(b).i)  + (1<<14) ) >> 15;\
          (m).i = ( ( (a).r*(b).i + (a).i*(b).r)  + (1<<14) ) >> 15;\
      }while(0)
#   define C_FIXDIV(c,div) \
    do{ (c).r /= div; (c).i /=div; }while(0)

#   define C_MULBYSCALAR( c, s ) \
    do{ (c).r = ( ( (c).r*(s) ) + (1<<14) ) >> 15;\
        (c).i = ( ( (c).i*(s) ) + (1<<14) ) >> 15; }while(0)

#else  /* not FIXED_POINT*/

#   define S_MUL(a,b) ( (a)*(b) )
#define C_MUL(m,a,b) \
    do{ (m).r = (a).r*(b).r - (a).i*(b).i;\
        (m).i = (a).r*(b).i + (a).i*(b).r; }while(0)
#   define C_FIXDIV(c,div) /* NOOP */
#   define C_MULBYSCALAR( c, s ) \
    do{ (c).r *= (s);\
        (c).i *= (s); }while(0)
#endif

#define  C_ADD( res, a,b)\
    do {    (res).r=(a).r+(b).r;  (res).i=(a).i+(b).i;  }while(0)
#define  C_SUB( res, a,b)\
    do {    (res).r=(a).r-(b).r;  (res).i=(a).i-(b).i;  }while(0)
#define C_ADDTO( res , a)\
    do {    (res).r += (a).r;  (res).i += (a).i;  }while(0)
#define C_SUBFROM( res , a)\
    do {    (res).r -= (a).r;  (res).i -= (a).i;  }while(0)

kiss_fft_cpx cexp(double phase) /* returns e ** (j*phase)   */
{
    kiss_fft_cpx x;
#ifdef FIXED_POINT    
    x.r = (kiss_fft_scalar) ( 32767*cos(phase) );
    x.i = (kiss_fft_scalar) ( 32767*sin(phase) );
#else
    x.r = cos(phase);
    x.i = sin(phase);
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
        C_FIXDIV(*Fout,2); C_FIXDIV(*Fout2,2);

        C_MUL (t,  *Fout2 , *tw1);
        tw1 += fstride;
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
    kiss_fft_cpx *tw1,*tw2,*tw3;
    kiss_fft_cpx scratch[6];

    Fout1 = Fout + m;
    Fout2 = Fout + 2*m;
    Fout3 = Fout + 3*m;
    tw3 = tw2 = tw1 = st->twiddles;

    do {
        C_FIXDIV(*Fout,4); C_FIXDIV(*Fout1,4); C_FIXDIV(*Fout2,4); C_FIXDIV(*Fout3,4);

        C_MUL(scratch[0],*Fout1 , *tw1 );
        C_MUL(scratch[1],*Fout2 , *tw2 );
        C_MUL(scratch[2],*Fout3 , *tw3 );

        C_SUB( scratch[5] , *Fout, scratch[1] );
        C_ADDTO(*Fout, scratch[1]);
        C_ADD( scratch[3] , scratch[0] , scratch[2] );
        C_SUB( scratch[4] , scratch[0] , scratch[2] );
        C_SUB( *Fout2, *Fout, scratch[3] );
        tw1 += fstride;
        tw2 += fstride*2;
        tw3 += fstride*3;
        C_ADDTO( *Fout , scratch[3] );

        if(st->inverse) {
            Fout1->r = scratch[5].r - scratch[4].i;
            Fout1->i = scratch[5].i + scratch[4].r;
            Fout3->r = scratch[5].r + scratch[4].i;
            Fout3->i = scratch[5].i - scratch[4].r;
        }else{
            Fout1->r = scratch[5].r + scratch[4].i;
            Fout1->i = scratch[5].i - scratch[4].r;
            Fout3->r = scratch[5].r - scratch[4].i;
            Fout3->i = scratch[5].i + scratch[4].r;
        }
        ++Fout; ++Fout1; ++Fout2; ++Fout3;
    }while(--m);
}

/* bfly3 is a optimization of bfly_generic for p==3 */
void bfly3(
         kiss_fft_cpx * Fout,
         int fstride,
         const kiss_fft_state * st,
         int m
         )
{
     kiss_fft_cpx *Fout0,*Fout1,*Fout2;
     kiss_fft_cpx *tw1,*tw2;
     kiss_fft_cpx scratch[5];
     kiss_fft_cpx epi3;
     epi3 = st->twiddles[fstride*m];

     Fout0=Fout;
     Fout1=Fout0+m;
     Fout2=Fout0+2*m;
     tw1=tw2=st->twiddles;

     do{
         C_FIXDIV(*Fout0,3); C_FIXDIV(*Fout1,3); C_FIXDIV(*Fout2,3);

         C_MUL(scratch[1],*Fout1 , *tw1);
         C_MUL(scratch[2],*Fout2 , *tw2);

         C_ADD(scratch[3],scratch[1],scratch[2]);
         C_SUB(scratch[0],scratch[1],scratch[2]);

         Fout1->r = Fout0->r - scratch[3].r/2;
         Fout1->i = Fout0->i - scratch[3].i/2;

         C_MULBYSCALAR( scratch[0] , epi3.i );

         C_ADDTO(*Fout0,scratch[3]);

         Fout2->r = Fout1->r + scratch[0].i;
         Fout2->i = Fout1->i - scratch[0].r;

         Fout1->r -= scratch[0].i;
         Fout1->i += scratch[0].r;

         ++Fout0;++Fout1;++Fout2;
         tw1 += fstride;
         tw2 += fstride*2;
     }while(--m);
}

/* bfly5 is a optimization of bfly_generic for p==5 */
void bfly5(
        kiss_fft_cpx * Fout,
        int fstride,
        const kiss_fft_state * st,
        int m
        )
{
    kiss_fft_cpx *Fout0,*Fout1,*Fout2,*Fout3,*Fout4;
    int u;
    kiss_fft_cpx scratch[13];
    kiss_fft_cpx * twiddles = st->twiddles;
    kiss_fft_cpx *tw;
    kiss_fft_cpx y1,y2;
    y1 = twiddles[fstride*m];
    y2 = twiddles[fstride*2*m];

    Fout0=Fout;
    Fout1=Fout0+m;
    Fout2=Fout0+2*m;
    Fout3=Fout0+3*m;
    Fout4=Fout0+4*m;

    tw=st->twiddles;
    for ( u=0; u<m; ++u ) {
        C_FIXDIV( *Fout0,5); C_FIXDIV( *Fout1,5); C_FIXDIV( *Fout2,5); C_FIXDIV( *Fout3,5); C_FIXDIV( *Fout4,5);
        scratch[0] = *Fout0;

        C_MUL(scratch[1] ,*Fout1, tw[u*fstride]);
        C_MUL(scratch[2] ,*Fout2, tw[2*u*fstride]);
        C_MUL(scratch[3] ,*Fout3, tw[3*u*fstride]);
        C_MUL(scratch[4] ,*Fout4, tw[4*u*fstride]);

        C_ADD( scratch[7],scratch[1],scratch[4]);
        C_SUB( scratch[10],scratch[1],scratch[4]);
        C_ADD( scratch[8],scratch[2],scratch[3]);
        C_SUB( scratch[9],scratch[2],scratch[3]);

        Fout0->r += scratch[7].r + scratch[8].r;
        Fout0->i += scratch[7].i + scratch[8].i;

        scratch[5].r = scratch[0].r + S_MUL(scratch[7].r,y1.r) + S_MUL(scratch[8].r,y2.r);
        scratch[5].i = scratch[0].i + S_MUL(scratch[7].i,y1.r) + S_MUL(scratch[8].i,y2.r);

        scratch[6].r =  S_MUL(scratch[10].i,y1.i) + S_MUL(scratch[9].i,y2.i);
        scratch[6].i = -S_MUL(scratch[10].r,y1.i) - S_MUL(scratch[9].r,y2.i);

        C_SUB(*Fout1,scratch[5],scratch[6]);
        C_ADD(*Fout4,scratch[5],scratch[6]);
 
        scratch[11].r = scratch[0].r + S_MUL(scratch[7].r,y2.r) + S_MUL(scratch[8].r,y1.r);
        scratch[11].i = scratch[0].i + S_MUL(scratch[7].i,y2.r) + S_MUL(scratch[8].i,y1.r);
        scratch[12].r = - S_MUL(scratch[10].i,y2.i) + S_MUL(scratch[9].i,y1.i);
        scratch[12].i = S_MUL(scratch[10].r,y2.i) - S_MUL(scratch[9].r,y1.i);

        C_ADD(*Fout2,scratch[11],scratch[12]);
        C_SUB(*Fout3,scratch[11],scratch[12]);

        ++Fout0;++Fout1;++Fout2;++Fout3;++Fout4;
    }
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
        case 2: bfly2(Fout,fstride,st,m); break;
        case 3: bfly3(Fout,fstride,st,m); break;
        case 4: bfly4(Fout,fstride,st,m); break;
        case 5: bfly5(Fout,fstride,st,m); break;
        default: bfly_generic(Fout,fstride,st,m,p); break;
    }
}

int allocsize(int nfft)
{
    int allocsize =  sizeof(kiss_fft_state)
        + sizeof(kiss_fft_cpx)*nfft /* twiddle factors*/
        + sizeof(kiss_fft_cpx)*nfft /* tmpbuf*/
        + sizeof(int)*nfft /* factors*/
        + sizeof(kiss_fft_cpx)*nfft; /* scratch*/
    return allocsize;   
}

void init_state(kiss_fft_state * st,int nfft,int inverse_fft)
{
    int nstages=0;
    int i;
    st->nfft=nfft;
    st->inverse = inverse_fft;
    st->twiddles = (kiss_fft_cpx*)(st+1); /* just beyond struct*/
    st->tmpbuf = (kiss_fft_cpx*)(st->twiddles + nfft);/*  just after twiddles*/
    st->scratch = (kiss_fft_cpx*)(st->tmpbuf + nfft);
    st->factors = (int*)(st->scratch + nfft);

    for (i=0;i<nfft;++i) {
        const double pi=3.14159265358979323846264338327;
        double phase = ( -2*pi /nfft ) * i;
        if (st->inverse)
            phase *= -1;
        st->twiddles[i] = cexp( phase );
    }

    while (nfft>1) {
        /* If you want a new radix, don't forget to put it here */
        const int divisors[] = {
		4,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,
		71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,
		149,151,157,163,167,173,179,181,191,193,197,199,-1};
        int p=nfft;
        i=0;
        while ( divisors[i] != -1 ) {
            if ( nfft %  divisors[i] == 0){
                p =  divisors[i];
                break;
            }
            ++i;
        }
        st->factors[2*nstages] = p;
        nfft /= p;
        st->factors[2*nstages+1] = nfft;
        ++nstages;
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
    kiss_fft_state * st=NULL;

    st = ( kiss_fft_state *)malloc( allocsize(nfft) );
    if (!st)
        return NULL;
    init_state( st ,nfft,inverse_fft );
    return st;
}

void * kiss_fft2d_alloc(int nrows,int ncols,int inverse_fft)
{
    kiss_fft2d_state *st = NULL;
    int size1,size2;
    size1 = allocsize(ncols);
    size2 = allocsize(nrows);

    st = (kiss_fft2d_state *) malloc ( sizeof(kiss_fft2d_state) + size1 + size2 );
    if (!st)
        return NULL;
    
    st->rowst = (kiss_fft_state *)(st+1); /*just beyond kiss_fft2d_state struct */
    st->colst = (kiss_fft_state *)( (char*)(st->rowst) + size1 );
    init_state (st->rowst, ncols, inverse_fft);
    init_state (st->colst, nrows, inverse_fft);
    return st;
}

void kiss_fft2d(const void * cfg,kiss_fft_cpx *f)
{
    /* 
     f is stored row-wise
     * */
    fprintf(stderr,"not yet implemented\n");
    exit(1);
}

void kiss_fft2d_io(const void * cfg,const kiss_fft_cpx * fin,kiss_fft_cpx * fout)
{
    /*just use the in-place version sinc the multi-dim needs two passes anyway*/
    kiss_fft2d_state *st = (kiss_fft2d_state *)cfg;
    memcpy(fout,fin,sizeof(kiss_fft_cpx) * st->rowst->nfft * st->colst->nfft );
    kiss_fft2d(cfg,fout);
}

/* original form of processing function, first release of KISS FFT was in-place.  This maintains API. */
void kiss_fft(const void * cfg,kiss_fft_cpx *f)
{
    const kiss_fft_state * st = cfg;
    if (st->nfft < 0) {
        kiss_fft2d(cfg,f);
    }else{
        memcpy(st->tmpbuf,f,sizeof(kiss_fft_cpx)*st->nfft);
        fft_work( f, st->tmpbuf, 1, st->factors,st );
    }
}

/* two buffer version of above */
void kiss_fft_io(const void * cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout)
{
    const kiss_fft_state * st = cfg;
    if (st->nfft < 0) {
        kiss_fft2d_io(cfg,fin,fout);
    }else{
        fft_work( fout, fin, 1, st->factors,st );
    }
}
