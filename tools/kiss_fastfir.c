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

/*
 Some definitions that allow real or complex filtering
*/
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

static int verbose=0;

void * kiss_fastfir_alloc(const kffsamp_t * imp_resp,size_t n_imp_resp,
        size_t * nfft,void * mem,size_t*lenmem);

/* n : the size of inbuf and outbuf in samples
   return value: the number of samples completely processed
   n-retval samples should be copied to the front of the next input buffer */
size_t kff_nocopy( void *st, const kffsamp_t * inbuf, kffsamp_t * outbuf, size_t n);



typedef struct {
    int nfft;
    size_t ngood;
    void * fftcfg;
    void * ifftcfg;
    kiss_fft_cpx * fir_freq_resp;
    kiss_fft_cpx * freqbuf;
    size_t n_freq_bins;
    kffsamp_t * tmpbuf;
}kiss_fastfir_state;


void * kiss_fastfir_alloc(
        const kffsamp_t * imp_resp,size_t n_imp_resp,
        size_t *pnfft, /* if <= 0, an appropriate size will be chosen */
        void * mem,size_t*lenmem)
{
    kiss_fastfir_state *st = NULL;
    size_t len_fftcfg,len_ifftcfg;
    size_t memneeded = sizeof(kiss_fastfir_state);
    char * ptr;
    size_t i;
    int nfft=0;
    float scale;
    int n_freq_bins;
    if (pnfft)
        nfft=*pnfft;

    if (nfft<=0) {
        /* determine fft size as next power of two at least 2x 
         the impulse response length*/
        i=n_imp_resp-1;
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
    /* tmpbuf */
    memneeded += sizeof(kffsamp_t) * nfft;
    /* fir_freq_resp */
    memneeded += sizeof(kiss_fft_cpx) * n_freq_bins;
    /* freqbuf */
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

    st->nfft = nfft;
    st->ngood = nfft - n_imp_resp + 1;
    st->n_freq_bins = n_freq_bins;
    ptr=(char*)(st+1);

    st->fftcfg = (void*)ptr;
    ptr += len_fftcfg;

    st->ifftcfg = (void*)ptr;
    ptr += len_ifftcfg;

    st->tmpbuf = (kffsamp_t*)ptr;
    ptr += sizeof(kffsamp_t) * nfft;

    st->freqbuf = (kiss_fft_cpx*)ptr;
    ptr += sizeof(kiss_fft_cpx) * n_freq_bins;
    
    st->fir_freq_resp = (kiss_fft_cpx*)ptr;
    ptr += sizeof(kiss_fft_cpx) * n_freq_bins;

    FFT_ALLOC (nfft,0,st->fftcfg , &len_fftcfg);
    FFT_ALLOC (nfft,1,st->ifftcfg , &len_ifftcfg);

    memset(st->tmpbuf,0,sizeof(kffsamp_t)*nfft);
    /*zero pad in the middle to left-rotate the impulse response 
      This puts the scrap samples at the end of the inverse fft'd buffer */
    st->tmpbuf[0] = imp_resp[ n_imp_resp - 1 ];
    for (i=0;i<n_imp_resp - 1; ++i) {
        st->tmpbuf[ nfft - n_imp_resp + 1 + i ] = imp_resp[ i ];
    }

    FFTFWD(st->fftcfg,st->tmpbuf,st->fir_freq_resp);

    /* TODO: this won't work for fixed point */
    scale = 1.0 / st->nfft;

    for ( i=0; i < st->n_freq_bins; ++i ) {
        st->fir_freq_resp[i].r *= scale;
        st->fir_freq_resp[i].i *= scale;
    }
    if (pnfft)
        *pnfft = nfft;
    return st;
}

static void fastconv1buf(const kiss_fastfir_state *st,const kffsamp_t * in,kffsamp_t * out)
{
    int i;
    /* multiply the frequency response of the input signal by
     that of the fir filter*/
    FFTFWD( st->fftcfg, in , st->freqbuf );
    for ( i=0; i<st->n_freq_bins; ++i ) {
        kiss_fft_cpx tmpsamp; 
        C_MUL(tmpsamp,st->freqbuf[i],st->fir_freq_resp[i]);
        st->freqbuf[i] = tmpsamp;
    }

    /* perform the inverse fft*/
    FFTINV(st->ifftcfg,st->freqbuf,out);
}

/* n : the size of inbuf and outbuf in samples
   return value: the number of samples completely processed
   n-retval samples should be copied to the front of the next input buffer */
size_t kff_nocopy(
        void *vst,
        const kffsamp_t * inbuf, 
        kffsamp_t * outbuf,
        size_t n)
{
    kiss_fastfir_state *st=(kiss_fastfir_state *)vst;
    size_t norig=n;
    while (n >= st->nfft ) {
        fastconv1buf(st,inbuf,outbuf);
        inbuf += st->ngood;
        outbuf += st->ngood;
        n -= st->ngood;
    }
    return norig - n;
}

size_t kff_flush(void *vst,const kffsamp_t * inbuf,kffsamp_t * outbuf,size_t n)
{
    kiss_fastfir_state *st=(kiss_fastfir_state *)vst;
    n -= kff_nocopy(vst,inbuf,outbuf,n);
    if (n) {
        memset(st->tmpbuf,0,sizeof(kffsamp_t)*st->nfft );
        memcpy(st->tmpbuf,inbuf,sizeof(kffsamp_t)*n );
        fastconv1buf(st,st->tmpbuf,st->tmpbuf);
        memcpy(outbuf,st->tmpbuf,sizeof(kffsamp_t)*n);
    }
    return n;
}

#ifdef FAST_FILT_UTIL
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>

int do_mmap_filter(int fd_src,int fd_dst,const kffsamp_t * imp_resp, size_t n_imp_resp, size_t nfft )
{
    int error_code=1;
    off_t src_len=0;
    void *psrc=NULL;
    void *pdst=NULL;

    src_len = lseek(fd_src,0,SEEK_END );
    if (src_len == (off_t)-1 ) goto nommap;
    if ( ftruncate(fd_dst,src_len) ) goto nommap;
    if ( lseek(fd_src,0,SEEK_SET) == (off_t)-1 ) goto nommap;
    if ( lseek(fd_dst,0,SEEK_SET) == (off_t)-1 ) goto nommap;
    psrc = mmap(0,src_len,PROT_READ,MAP_SHARED,fd_src,0);
    if (psrc==NULL) goto nommap;
    pdst = mmap(0,src_len,PROT_WRITE,MAP_SHARED,fd_dst,0);
    if (psrc==NULL) goto nommap;
    {
        const kffsamp_t *inbuf= (const kffsamp_t *)psrc;
        kffsamp_t *outbuf= (kffsamp_t *)pdst;

        size_t noutbuf;
        size_t nsamps_in = src_len / sizeof(kffsamp_t);
        
        void * cfg = kiss_fastfir_alloc(imp_resp,n_imp_resp,&nfft,0,0);
        
        noutbuf = kff_nocopy( cfg, (kffsamp_t*)psrc, (kffsamp_t*)pdst, nsamps_in );

        kff_flush( cfg,inbuf+noutbuf,outbuf+noutbuf,nsamps_in-noutbuf);

        free(cfg);
    }
    error_code=0;
nommap:
    if (psrc) munmap(psrc,src_len);
    if (pdst) munmap(pdst,src_len);
    if (error_code)
        perror("mapping");
    return error_code;
}


void do_file_filter(
        FILE * fin,
        FILE * fout,
        const kffsamp_t * imp_resp,
        size_t n_imp_resp,
        size_t nfft )
{
    void * cfg = kiss_fastfir_alloc(imp_resp,n_imp_resp,&nfft,0,0);
    size_t max_buf=4*nfft;
    kffsamp_t *inbuf = (kffsamp_t*)malloc(max_buf*sizeof(kffsamp_t));
    kffsamp_t *outbuf = (kffsamp_t*)malloc(max_buf*sizeof(kffsamp_t));
    size_t ninbuf,noutbuf;
    size_t idx_inbuf=0;
    int zpad=0;

    do{
        ninbuf = fread(&inbuf[idx_inbuf],sizeof(inbuf[0]),max_buf-idx_inbuf,fin );
        if (ninbuf==0) {
            /* zero pad the input to the fft size */
            zpad = nfft - idx_inbuf;
            if (verbose) fprintf(stderr,"zero padding %d samples\n",zpad);

            /* We don't strictly need to do the zero padding,since the output samples
             * will get thrown away anyway.  But the zeros should help keep the fft magnitudes in
             * check and thus decrease errors.  */
            memset(&inbuf[idx_inbuf],0,sizeof(inbuf[0])*zpad);
            ninbuf = nfft;/* force an fft */
        }else{
            ninbuf += idx_inbuf;
        }

        /*perform the fast convolution filtering*/
        noutbuf = kff_nocopy( cfg, inbuf, outbuf, ninbuf );
        if (verbose) fprintf(stderr,"kff_nocopy(,,,%d) -> %d\n",ninbuf,noutbuf);

        /* move the unconsumed input samples to the front */
        idx_inbuf = ninbuf - noutbuf;
        memmove( inbuf , &inbuf[noutbuf] , sizeof(inbuf[0])*(ninbuf-noutbuf) );
        if (verbose) fprintf(stderr,"moved %d samples to front of buffer\n",idx_inbuf);

        /*whatever padding was done to get the size up to nfft, must now be ignored*/
        noutbuf -= zpad;

        if ( fwrite( outbuf, sizeof(outbuf[0]), noutbuf, fout) != noutbuf ) {
            fprintf(stderr,"short write %d \n",noutbuf);
            fprintf(stderr,"zpad= %d \n",zpad);
            exit(1);
        }
    }while ( ! zpad );
    fclose(fout);
    free(cfg);
    free(inbuf);
    free(outbuf);
}

int main(int argc,char**argv)
{
    kffsamp_t * h;
    size_t nh,nfft=0;
    FILE *fin=stdin;
    FILE *fout=stdout;
    FILE *filtfile=NULL;
    while (1) {
        int c=getopt(argc,argv,"n:h:i:o:v");
        if (c==-1) break;
        switch (c) {
            case 'v':
                verbose=1;
                break;
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
                fout = fopen(optarg,"w+b");
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
    if (verbose) fprintf(stderr,"%d samples in FIR filter\n",nh);
    h = (kffsamp_t*)malloc(sizeof(kffsamp_t)*nh);
    fseek(filtfile,0,SEEK_SET);
    fread(h,sizeof(kffsamp_t),nh,filtfile);
    fclose(filtfile);
 
#if 0
    if (do_mmap_filter( fileno(fin), fileno(fout), h,nh,nfft ) )
#endif        
        do_file_filter( fin, fout, h,nh,nfft);

    if (fout!=stdout) fclose(fout);
    if (fin!=stdin) fclose(fin);

    return 0;
}
#endif
