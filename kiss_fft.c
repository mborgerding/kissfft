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
    kiss_fft_cpx* twiddle;
    int * swap_indices;
    int nswap;
}kiss_fft_state;

#ifdef FIXED_POINT

#   define TWIDDLE_RANGE ( (1<<15) - 1 )
#   define C_ADD(x,a,b) \
     do{ (x).r = ( ( (a).r+(b).r +1) >> 1 );\
         (x).i = ( ( (a).i+(b).i +1) >> 1 ); }while(0)
#   define C_SUB(x,a,b) \
     do{ (x).r = ( ( (a).r-(b).r +1) >> 1 );\
         (x).i = ( ( (a).i-(b).i +1) >> 1 ); }while(0)
#   define C_MUL(m,a,b) \
     do{ (m).r = ( ( (a).r*(b.r) - (a).i*(b).i)  + (1<<14) ) >> 15;\
         (m).i = ( ( (a).r*(b.i) + (a).i*(b).r)  + (1<<14) ) >> 15;\
         }while(0)

#else // floating point

#define C_ADD(x,a,b) \
    do{ (x).r = (a).r+(b).r;\
        (x).i = (a).i+(b).i;}while(0)
#define C_SUB(x,a,b) \
    do{ (x).r = (a).r-(b).r;\
        (x).i = (a).i-(b).i;}while(0)
#define C_MUL(m,a,b) \
    do{ (m).r = (a).r*(b.r) - (a).i*(b).i;\
        (m).i = (a).r*(b.i) + (a).i*(b).r; }while(0)

#endif


// create the twiddle factors
static 
void make_twiddles(kiss_fft_cpx * twiddle,int nfft,int inverse_fft) 
{
    const double pi=3.14159265358979323846264338327;
    int n;
    for (n=0;n< (nfft>>1);++n) {
#ifdef FIXED_POINT
        twiddle[n].r =  (kiss_fft_scalar)(TWIDDLE_RANGE*cos(2*pi*n/nfft)+.5);
        twiddle[n].i = -(kiss_fft_scalar)(TWIDDLE_RANGE*sin(2*pi*n/nfft)+.5);
#else
        twiddle[n].r = cos(2*pi*n/nfft);
        twiddle[n].i = -sin(2*pi*n/nfft);
#endif
        if (inverse_fft)
            twiddle[n].i *= -1;
    }
}


// reverse the bits for numbers from 0 to N-1
static
int bit_reverse(int i,int N)
{
    int rev=0;
    while ( N >>= 1 ) {
        rev = rev*2 + (i&1);
        i>>=1;
    }
    return rev;
}

// make a list of index swaps that need to be done for bit-reversed addressing
static
int make_bit_reverse_indices(int *swap_indices ,int nfft)
{
    int n,nswap;
    nswap=0;
    // set up swap indices to unwrap bit-reversed addressing
    for ( n=0; n<nfft; ++n ) {
        swap_indices[nswap*2] = bit_reverse(n,nfft);
        if ( n < swap_indices[nswap*2] ) {
            swap_indices[nswap*2+1] = n;
            ++nswap;
        }
    }
    return nswap;
}

static
void undo_bit_rev( kiss_fft_cpx * f, const int * swap_indices,int nswap)
{
    int i,i0,i1;
    kiss_fft_cpx tmp;

    // unwrap bit-reverse indexing
    for ( i=0 ; i<nswap ; ++i ) {
        i0= *swap_indices++;
        i1= *swap_indices++;
        tmp = f[i0];
        f[i0] = f[i1];
        f[i1] = tmp;
    }
}

// the heart of the fft
static
void fft_work(int N,kiss_fft_cpx *f,const kiss_fft_cpx * twid,int twid_step)
{
    N >>= 1;
    {   // declare automatic variables in here to reduce stack usage
        int n;
        kiss_fft_cpx csum,cdiff;
 
        for (n=0;n<N;++n) {
            C_ADD(csum,  f[n] , f[N+n] );
            C_SUB(cdiff,  f[n] , f[N+n] );
            f[n]=csum;
            C_MUL(f[N+n], cdiff, twid[n*twid_step] );
        }
    }
    if (N==1) return;
    fft_work(N,f,twid,twid_step*2); 
    fft_work(N,f+N,twid,twid_step*2);
}


void * kiss_fft_alloc(int nfft,int inverse_fft)
{
    kiss_fft_state * st=NULL;
    // allocate one large buffer to hold the state and the buffers for twiddles and bit-rev indices
    int size = sizeof(kiss_fft_state) + (nfft>>1)*sizeof(kiss_fft_cpx) + nfft*sizeof(int);

    st = ( kiss_fft_state *)malloc(size);
    if (st) { 
        st->twiddle = (kiss_fft_cpx*)(st+1);
        st->swap_indices = (int*)( st->twiddle + (nfft>>1) );
        st->nfft = nfft;
        make_twiddles(st->twiddle,nfft,inverse_fft);
        st->nswap = make_bit_reverse_indices(st->swap_indices,nfft);
    }
    return st;
}

void kiss_fft(const void * cfg,kiss_fft_cpx *f)
{
    const kiss_fft_state * st = cfg;
    fft_work(st->nfft, f, st->twiddle , 1 );
    undo_bit_rev(f,st->swap_indices,st->nswap);
}

#ifdef FFT_UTIL

#include <stdio.h>
#include <string.h>

int main(int argc,char ** argv)
{
    kiss_fft_cpx * buf=NULL;
    void *st;
    int nfft=1024;
    int inverse_fft=0;

    fprintf(stderr,"sizeof(kiss_fft_cpx) == %d\n",sizeof(kiss_fft_cpx) );

    if (argc>1 && 0==strcmp("-i",argv[1])) {
        inverse_fft = 1;
        --argc;++argv;
    }

    if (argc>1) {
        nfft = atoi(argv[1]);
    }

    buf=(kiss_fft_cpx*)malloc( sizeof(kiss_fft_cpx) * nfft );
    st = kiss_fft_alloc( nfft ,inverse_fft );

    while ( fread( buf , sizeof(kiss_fft_cpx) , nfft , stdin ) > 0 ) {
        kiss_fft( st , buf );
        fwrite( buf , sizeof(kiss_fft_cpx) , nfft , stdout );
    }

    free(st);
    free(buf);

    return 0;
}

#endif
