#include "kiss_fft.h"

typedef struct
{
    int nfft;
    int inverse;
    const void * cfg;
    void * next;
} cached_fft;

static cached_fft *cache_root=NULL;

static const void * find_cached_fft(int nfft,int inverse)
{
    size_t len;
    cached_fft * cur=cache_root;
    cached_fft * prev=NULL;
    while ( cur ) {
        if ( cur->nfft == nfft && inverse == cur->inverse )
            return cur->cfg;
        prev = cur;
        cur = (cached_fft*)prev->next;
    }

    // no cached node found, need to create a new one
    kiss_fft_alloc(nfft,inverse,0,&len);
    cur = (cached_fft*)malloc(sizeof(cached_fft) + len );
    if (cur == NULL)
        return NULL;
    cur->cfg = kiss_fft_alloc(nfft,inverse,cur+1,&len);
    cur->nfft=nfft;
    cur->inverse=inverse;
    if ( prev )
        prev->next = cur;
    else
        cache_root = cur;
    return cur->cfg;
}

void kfc_fft(int nfft, const kiss_fft_cpx * fin,kiss_fft_cpx * fout)
{
    kiss_fft( find_cached_fft(nfft,0),fin,fout );
}

void kfc_ifft(int nfft, const kiss_fft_cpx * fin,kiss_fft_cpx * fout)
{
    kiss_fft( find_cached_fft(nfft,1),fin,fout );
}

void kfc_cleanup()
{
    cached_fft * cur=cache_root;
    cached_fft * next=NULL;
    while (cur){
        next = (cached_fft*)cur->next;
        free(cur);
        cur=(cached_fft*)next;
    }
    cache_root = NULL;
}
