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
    int nr;    // N remaining
    int radix; // radix of the stage
    kiss_fft_cpx * twiddle;
}kiss_fft_stage;



#define MAX_STAGES 20
typedef struct {
    int nstages;
    int allocsize;
    kiss_fft_stage stages[MAX_STAGES];
    int * unscramble_indices;
}kiss_fft_state;

#ifdef FIXED_POINT

#   define TWIDDLE_RANGE ( (1<<FIXED_POINT_FRAC_BITS) - 1 )

    /* The fixed point versions of C_ADD and C_SUB divide by two
     * after the add/sub.
     * This is to prevent overflow.
     * The cumulative effect is to multiply the sequence by 1/Nfft.
     * */
#   define C_ADD(x,a,b) \
     do{ (x).r = ( ( (a).r+(b).r +1) >> 1 );\
         (x).i = ( ( (a).i+(b).i +1) >> 1 ); }while(0)
#   define C_SUB(x,a,b) \
     do{ (x).r = ( ( (a).r-(b).r +1) >> 1 );\
         (x).i = ( ( (a).i-(b).i +1) >> 1 ); }while(0)

    /*  We don't have to worry about overflow from multiplying by twiddle factors since they
     *  all have unity magnitude.  Still need to shift away fractional bits after adding 1/2 for
     *  rounding.
     * */
#   define C_MUL(m,a,b) \
     do{ (m).r = ( ( (a).r*(b).r - (a).i*(b).i)  + (TWIDDLE_RANGE>>1) ) >> FIXED_POINT_FRAC_BITS;\
         (m).i = ( ( (a).r*(b).i + (a).i*(b).r)  + (TWIDDLE_RANGE>>1) ) >> FIXED_POINT_FRAC_BITS;\
         }while(0)
#else // not FIXED_POINT

#define C_ADD(x,a,b) \
    do{ (x).r = (a).r+(b).r;\
        (x).i = (a).i+(b).i;}while(0)
#define C_SUB(x,a,b) \
    do{ (x).r = (a).r-(b).r;\
        (x).i = (a).i-(b).i;}while(0)
#define C_MUL(m,a,b) \
    do{ (m).r = (a).r*(b).r - (a).i*(b).i;\
        (m).i = (a).r*(b).i + (a).i*(b).r; }while(0)

#endif


// create the twiddle factors
static 
void make_twiddles(kiss_fft_cpx * twiddle,int ntwid,int symmetry,int inverse_fft) 
{
    const double pi=3.14159265358979323846264338327;
    int n;
    int denom = ntwid*symmetry;
    
    for (n=0;n<ntwid;++n) {
        twiddle[n].r = cos(2*pi*n/denom);
        twiddle[n].i = -sin(2*pi*n/denom);
        // inverse fft uses complex conjugate twiddle factors
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
void fft_work( kiss_fft_cpx * f,const  kiss_fft_stage * stage , int nstages )
{
    int n;
    // only works for 2s
    
    {   // declare automatic variables in here to reduce stack usage
        kiss_fft_cpx csum,cdiff;
        kiss_fft_cpx *f1 = f;
        kiss_fft_cpx *f2 = f+ stage->nr;
        const kiss_fft_cpx * twid = stage->twiddle;
        for ( n = 0 ; n < stage->nr ; ++n )
        {
            C_ADD( csum,  *f1 , *f2 );
            C_SUB( cdiff, *f1 , *f2 );
            *f1++ = csum;
            C_MUL(*f2, cdiff, *twid);
            ++f2;
            ++twid;
        }
    }
    if ( nstages == 1 )
        return;
    
    for (n=0;n<stage->radix;++n)
        fft_work( f + n * stage->nr ,stage+1,nstages-1);
}
static
kiss_fft_state * decompose( kiss_fft_state * st, int * pnfft,int radix,int inverse_fft)
{
    int nfft = *pnfft;

    while (nfft>1 && (nfft%radix)==0) {
        int newsize;
        int nstages=st->nstages;
        nfft /= radix;
        fprintf(stderr,"%d ",radix);

        newsize = st->allocsize + nfft * sizeof(kiss_fft_cpx);

        st=(kiss_fft_state*)realloc(st,newsize);

        st->stages[nstages].radix=radix;
        st->stages[nstages].nr=nfft;
        st->stages[nstages].twiddle = (kiss_fft_cpx*)((char*)st + st->allocsize);

        make_twiddles( st->stages[nstages].twiddle,nfft,radix,inverse_fft );
        st->allocsize = newsize;
        st->nstages = nstages+1;
    }
    *pnfft = nfft;
    return st;
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
    int tmp;
    int i;
    int allocsize=sizeof(kiss_fft_state) + nfft*sizeof(int);

    st = ( kiss_fft_state *)malloc( allocsize );
    st->nstages=0;
    st->allocsize=allocsize;
    st->unscramble_indices = (int*)(st+1);// just past end of buffer

    for (i=0;i<nfft;++i)
        st->unscramble_indices[i] = i;

    fprintf(stderr,"factoring %d: ",nfft);

    tmp=nfft;
    st = decompose( st,&tmp,2,inverse_fft);
    fprintf(stderr,"\n");

    // TODO factor 3,5,7,11 ???
    if ( tmp > 1 ) {
        fprintf(stderr,"%d not factorable by 2,3 (%d remaining)\n",nfft,tmp);
        free(st);
        return NULL;
    }

    return st;
}

void kiss_fft(const void * cfg,kiss_fft_cpx *f)
{
    const kiss_fft_state * st = cfg;
    fft_work( f, st->stages , st->nstages );
}
