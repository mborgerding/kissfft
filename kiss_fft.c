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

const double pi=3.14159265358979323846264338327;

#define MAX_STAGES 20
typedef struct {
    int nfft;
    int inverse;
}kiss_fft_state;

#define C_ADD(x,a,b) \
    do{ (x).r = (a).r+(b).r;\
        (x).i = (a).i+(b).i;}while(0)
#define C_MUL(m,a,b) \
    do{ (m).r = (a).r*(b).r - (a).i*(b).i;\
        (m).i = (a).r*(b).i + (a).i*(b).r; }while(0)

static
kiss_fft_cpx cexp(double phase)
{
    kiss_fft_cpx x;
    x.r = cos(phase);
    x.i = -sin(phase);
    return x;
}

static kiss_fft_cpx cadd(kiss_fft_cpx a,kiss_fft_cpx b)
{
    kiss_fft_cpx c;
    C_ADD(c,a,b);
    return c;
}

static kiss_fft_cpx cmul(kiss_fft_cpx a,kiss_fft_cpx b)
{
    kiss_fft_cpx c;
    C_MUL(c,a,b);
    return c;
}

// the heart of the fft
static 
void fft_work(
        kiss_fft_cpx * Fout,
        const kiss_fft_cpx * f,
        int fstride,
        int n,
        int inverse,
        kiss_fft_cpx * scratch
        )
{
    int m,p=0,q,q1,u,k;
    kiss_fft_cpx t;
    double two_pi_divN;


    if (n==1) {
        *Fout = *f;
        return;
    }

    two_pi_divN = 2*pi/n;
    if (inverse)
        two_pi_divN *= -1;

    if (n%2 == 0) p=2;
    else if(n%3 == 0) p=3;
    else if(n%5 == 0) p=5;
    else if(n%7 == 0) p=7;
    else {
        fprintf(stderr,"%d is not divisible by %d\n",n,p);
        p=n;
        exit(1);
    }
    m = n/p;
    /* n = stage->n; m = stage->m; p = stage->p; */

#ifdef LOUD        
    printf("n=%d,p=%d\n",n,p);
    for (k=0;k<n;++k) {
        t=f[k*fstride];
        printf("(%.3f,%.3f) ",t.r,t.i);
    }
    printf(" <<< fin \n");
#endif

    for (q=0;q<p;++q) {
#ifdef LOUD
        for (k=0;k<m;++k) {
            t = (f+q*fstride)[fstride*p*k];
            printf("(%.3f,%.3f) ",t.r,t.i);
        }
        printf(" <<< f[fstride*k+q] \n");
#endif
        fft_work( Fout + m*q, f+q*fstride, fstride*p, m,inverse, scratch );
    }

#ifdef LOUD
    printf("twiddling n=%d,p=%d\n",n,p);
#endif
    
    for ( u=0; u<m; ++u ) {
        for ( q1=0 ; q1<p ; ++q1 ) {
            scratch[q1] = Fout[ u+q1*m  ];
#ifdef LOUD
            printf("(%.3f,%.3f) ",scratch[q1].r,scratch[q1].i);
#endif
        }
#ifdef LOUD
        printf(" <<< scratch \n");
#endif

        for ( q1=0 ; q1<p ; ++q1 ) {
            k=q1*m+u;
            Fout[ k ] = scratch[0];
            for (q=1;q<p;++q ) {
                //t = e ** ( j*2*pi*k*q/n );
                t = cexp( two_pi_divN * k * q );
                Fout[ k ] = cadd( Fout[ k ] , 
                                  cmul( scratch[q] , t ) );
            }
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
    kiss_fft_state * st=NULL;
    st = ( kiss_fft_state *)malloc( sizeof(kiss_fft_state) );
    st->nfft=nfft;
    st->inverse = inverse_fft;

    return st;
}

void kiss_fft(const void * cfg,kiss_fft_cpx *f)
{
    int i,n;
    const kiss_fft_state * st = cfg;
    kiss_fft_cpx tmpbuf[1024];
    kiss_fft_cpx scratch[1024];
    n = st->nfft;

    for (i=0;i<n;++i)
        tmpbuf[i] = f[i];

    fft_work( f, tmpbuf, 1, n, st->inverse, scratch );
}
