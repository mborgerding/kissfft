#include "kiss_fft.h"

typedef struct
{
    int nfft;
    int inverse;
    void * cfg;
    void * next;
} cached_fft;

static cached_fft *cache_root=NULL;
static int ncached=0;

static const void * find_cached_fft(int nfft,int inverse)
{
    size_t len;
    cached_fft * cur=cache_root;
    cached_fft * prev=NULL;
    while ( cur ) {
        if ( cur->nfft == nfft && inverse == cur->inverse )
            break;//found the right node
        prev = cur;
        cur = (cached_fft*)prev->next;
    }
    if (cur== NULL) {
        // no cached node found, need to create a new one
        kiss_fft_alloc(nfft,inverse,0,&len);
        cur = (cached_fft*)malloc(sizeof(cached_fft) + len );
        if (cur == NULL)
            return NULL;
        cur->cfg = cur+1;
        kiss_fft_alloc(nfft,inverse,cur->cfg,&len);
        cur->nfft=nfft;
        cur->inverse=inverse;
        cur->next = NULL;
        if ( prev )
            prev->next = cur;
        else
            cache_root = cur;
        ++ncached;
    }
    return cur->cfg;
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
    ncached=0;
    cache_root = NULL;
}
void kfc_fft(int nfft, const kiss_fft_cpx * fin,kiss_fft_cpx * fout)
{
    kiss_fft( find_cached_fft(nfft,0),fin,fout );
}

void kfc_ifft(int nfft, const kiss_fft_cpx * fin,kiss_fft_cpx * fout)
{
    kiss_fft( find_cached_fft(nfft,1),fin,fout );
}

#ifdef KFC_TEST
static void check(int nc)
{
    if (ncached != nc) {
        fprintf(stderr,"ncached should be %d,but it is %d\n",nc,ncached);
        exit(1);
    }
}

int main(int argc,char ** argv)
{
    kiss_fft_cpx buf1[1024],buf2[1024];
    memset(buf1,0,sizeof(buf1));
    check(0);
    kfc_fft(512,buf1,buf2);
    check(1);
    kfc_fft(512,buf1,buf2);
    check(1);
    kfc_ifft(512,buf1,buf2);
    check(2);
    kfc_cleanup();
    check(0);
    return 0;
}
#endif
