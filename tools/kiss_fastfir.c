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

#ifdef REAL_FASTFIR
#include "kiss_fftr.h"
typedef kiss_fft_scalar kffsamp_t;
#define FFT_ALLOC kiss_fftr_alloc
#define FFTFWD kiss_fftr
#define FFTINV kiss_fftri
#else
typedef kiss_fft_cpx kffsamp_t;
#define FFT_ALLOC kiss_fft_alloc
#define FFTFWD kiss_fft
#define FFTINV kiss_fft
#endif

void * kiss_fastfir_alloc(const kffsamp_t * imp_resp,size_t n_imp_resp,
        int nfft,void * mem,size_t*lenmem);

size_t kiss_fastfir(const void * cfg,
        const kffsamp_t *in, size_t nin,
        kffsamp_t *out, size_t nout);

typedef struct {
    int minus5; /*magic */
    int nfft;
    size_t n_scrap;
    void * fftcfg;
    void * ifftcfg;
    kiss_fft_cpx * fir_freq_resp;
    size_t n_freq_bins;
    size_t bufin_idx;
    kffsamp_t * bufin;
    size_t bufout_idx;
    kffsamp_t * bufout;
    kiss_fft_cpx * tmpbuf;
}kiss_fastfir_state;


void * kiss_fastfir_alloc(const kffsamp_t * imp_resp,size_t n_imp_resp,
        int nfft,void * mem,size_t*lenmem)
{
    kiss_fastfir_state *st = NULL;
    size_t len_fftcfg,len_ifftcfg;
    size_t memneeded = sizeof(kiss_fastfir_state);
    char * ptr;
    size_t i;
    float scale;
    int n_freq_bins;

    if (nfft<=0) {
        /* determine fft size as next power of two at least 2x 
         the impulse response length*/
        int i=n_imp_resp-1;
        nfft=2;
        do{
             nfft<<=1;
        }while (i>>=1);
    }

#ifdef REAL_FASTFIR
    n_freq_bins = nfft/2 + 1;
#else    
    n_freq_bins = nfft;
#endif
    /*fftcfg*/
    FFT_ALLOC (nfft, 0, NULL, &len_fftcfg);
    memneeded += len_fftcfg;  
    /*ifftcfg*/
    FFT_ALLOC (nfft, 1, NULL, &len_ifftcfg);
    memneeded += len_ifftcfg;  

    /* fir_freq_resp */
    memneeded += sizeof(kiss_fft_cpx) * n_freq_bins;
    /* bufin */
    memneeded += sizeof(kffsamp_t) * nfft;
    /* bufout */
    memneeded += sizeof(kffsamp_t) * nfft;
    /* tmpbuf */
    memneeded += sizeof(kiss_fft_cpx) * n_freq_bins;

    if (lenmem == NULL) {
        st = (kiss_fastfir_state *) malloc (memneeded);
    } else {
        if (*lenmem >= memneeded)
            st = (kiss_fastfir_state *) mem;
        *lenmem = memneeded;
    }
    if (!st)
        return NULL;

    st->minus5 = -5;
    st->nfft = nfft;
    st->n_scrap = n_imp_resp-1;
    st->n_freq_bins = n_freq_bins;
    st->bufin_idx = 0;
    st->bufout_idx = nfft;
    ptr=(char*)(st+1);

    st->fftcfg = (void*)ptr;
    ptr += len_fftcfg;

    st->ifftcfg = (void*)ptr;
    ptr += len_ifftcfg;

    st->fir_freq_resp = (kiss_fft_cpx*)ptr;
    ptr += sizeof(kiss_fft_cpx) * n_freq_bins;

    st->bufin = (kffsamp_t*)ptr;
    ptr += sizeof(kffsamp_t) * nfft;

    st->bufout = (kffsamp_t*)ptr;
    ptr += sizeof(kffsamp_t) * nfft;
    
    st->tmpbuf = (kiss_fft_cpx*)ptr;
    ptr += sizeof(kiss_fft_cpx) * n_freq_bins;

    FFT_ALLOC (nfft,0,st->fftcfg , &len_fftcfg);
    FFT_ALLOC (nfft,1,st->ifftcfg , &len_ifftcfg);

    memset(st->bufin,0,sizeof(kffsamp_t)*nfft);
    memcpy(st->bufin,imp_resp,sizeof(kffsamp_t)*n_imp_resp);
    FFTFWD(st->fftcfg,st->bufin,st->fir_freq_resp);

    scale = 1.0 / st->nfft;

    for ( i=0; i < st->n_freq_bins; ++i ) {
        st->fir_freq_resp[i].r *= scale;
        st->fir_freq_resp[i].i *= scale;
    }

    return st;
}

static
size_t write_output(kiss_fastfir_state *st,
        kffsamp_t *out,size_t * pnout,size_t zpadded)
{
    size_t nout = *pnout;
    size_t n2flush = st->nfft - st->bufout_idx;
    if (zpadded)
        n2flush -= zpadded;
    if ( nout < n2flush )
        n2flush=nout;
    memcpy(out,st->bufout + st->bufout_idx, sizeof(kffsamp_t)*n2flush );
    st->bufout_idx += n2flush;
    *pnout = nout - n2flush;
    return n2flush;
}

static void do_fastconv(kiss_fastfir_state *st)
{
    int i;
    if ( st->bufout_idx < st->nfft ) {
        fprintf(stderr,"kiss_fastfir warning: "
                " output buffer size must be >= input buffer size,"
                " %d samples lost\n",st->nfft - st->bufout_idx );
    }
    /*FFT st->bufin to st->tmpbuf*/
    FFTFWD(st->fftcfg,st->bufin,st->tmpbuf);

    /* shift tail to front of input buffer*/
    memcpy( st->bufin, 
            st->bufin + st->nfft - st->n_scrap,
            sizeof(kiss_fft_cpx)*st->n_scrap);
    /*set input idx to the next input spot*/
    st->bufin_idx = st->n_scrap;

    /* multiply the frequency response of the input signal by*/
    /* that of the fir filter*/
    for ( i=0; i<st->n_freq_bins; ++i ) {
        kiss_fft_cpx tmpsamp; 
        C_MUL(tmpsamp,st->tmpbuf[i],st->fir_freq_resp[i]);
        st->tmpbuf[i] = tmpsamp;
    }

    /* perform the inverse fft*/
    FFTINV(st->ifftcfg,st->tmpbuf,st->bufout);

    /* need to skip over junk caused by circular convolution*/
    st->bufout_idx = st->n_scrap;
}
    
size_t kiss_fastfir(const void * cfg,
        const kffsamp_t *in, size_t nin,
        kffsamp_t *out, size_t nout_avail)
{
    size_t nout_orig=nout_avail;
    kiss_fastfir_state *st = ( kiss_fastfir_state *)cfg;

    out += write_output(st,out,&nout_avail,0);

    if ( nin <= 0 ) {
        size_t zero_pad = st->nfft - st->bufin_idx;
        memset( st->bufin + st->bufin_idx, 0, zero_pad*sizeof(kiss_fft_cpx) );
        st->bufin_idx = st->nfft;
        do_fastconv(st);
        fprintf(stderr,"padded with %d zeros\n",zero_pad);
        return write_output(st,out,&nout_avail,zero_pad);
    }
    
    while (nin--) {
        /* copy the input sample to bufin*/
        st->bufin[st->bufin_idx++] = *in++;

        /* when the input buffer is full, perform fast convolution*/
        if ( st->bufin_idx == st->nfft ) {
            do_fastconv(st);
            /* write the output buffer*/
            out += write_output(st,out,&nout_avail,0);
        }
    }
    return nout_orig - nout_avail;
}

#ifdef FAST_FILT_UTIL

#define BUFLEN 1024

void do_filter(
        FILE * fin,
        FILE * fout,
        const kffsamp_t * imp_resp,
        size_t n_imp_resp,
        size_t nfft)
{
    void * cfg = kiss_fastfir_alloc(imp_resp,n_imp_resp,nfft,0,0);
    kffsamp_t inbuf[BUFLEN],outbuf[BUFLEN];
    size_t ninbuf,noutbuf;
    do{
        ninbuf = fread(inbuf,sizeof(inbuf[0]),BUFLEN,fin );
        /* when ninbuf <= 0, that signals a flush*/
        noutbuf = kiss_fastfir(cfg,inbuf,ninbuf,outbuf,BUFLEN);
        if ( fwrite(outbuf,sizeof(outbuf[0]),noutbuf,fout) != noutbuf ) {
            fprintf(stderr,"short write\n");
            exit(1);
        }
    }while(ninbuf>0);
    fclose(fout);
    free(cfg);
}

#include <unistd.h>
int main(int argc,char**argv)
{
    kffsamp_t * h;
    size_t nh,nfft=0;
    FILE *fin=stdin;
    FILE *fout=stdout;
    FILE *filtfile=NULL;
    while (1) {
        int c=getopt(argc,argv,"n:h:i:o:");
        if (c==-1) break;
        switch (c) {
            case 'n':
                nfft=atoi(optarg);
                break;
            case 'i':
                fin = fopen(optarg,"rb");
                if (fin==NULL) {
                    perror(optarg);
                    exit(1);
                }
                break;
            case 'o':
                fout = fopen(optarg,"wb");
                if (fout==NULL) {
                    perror(optarg);
                    exit(1);
                }
                break;
            case 'h':
                filtfile = fopen(optarg,"rb");
                if (filtfile==NULL) {
                    perror(optarg);
                    exit(1);
                }
                break;
            case '?':
                     fprintf(stderr,"usage options:\n"
                            "\t-i filename: input file\n"
                            "\t-o filename: output(filtered) file\n"
                            "\t-h filename: impulse response\n");
                     exit (1);
            default:fprintf(stderr,"bad %c\n",c);break;
        }
    }
    if (filtfile==NULL) {
        fprintf(stderr,"You must supply the FIR coeffs via -h\n");
        exit(1);
    }
    fseek(filtfile,0,SEEK_END);
    nh = ftell(filtfile) / sizeof(kffsamp_t);
    fprintf(stderr,"%d samples in FIR filter\n",nh);
    h = (kffsamp_t*)malloc(sizeof(kffsamp_t)*nh);
    fseek(filtfile,0,SEEK_SET);
    fread(h,sizeof(kffsamp_t),nh,filtfile);
    fclose(filtfile);
 
    do_filter( fin, fout, h,nh,nfft);

    if (fout!=stdout) fclose(fout);
    if (fin!=stdin) fclose(fin);

    return 0;
}
#endif
