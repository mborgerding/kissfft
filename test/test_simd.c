#include <kiss_fftnd.h>

static void test1(void)
{
    int is_inverse = 1;
    int n[2] = {256,256};
    size_t nbytes = sizeof(kiss_fft_cpx)*n[0]*n[1];

    kiss_fft_cpx * inbuf = KISS_FFT_MALLOC(nbytes);
    kiss_fft_cpx * outbuf = KISS_FFT_MALLOC(nbytes);
    memset(inbuf,0,nbytes);
    memset(outbuf,0,nbytes);

    kiss_fftnd_cfg cfg = kiss_fftnd_alloc(n,2,is_inverse,0,0);
    kiss_fftnd(cfg,inbuf,outbuf);
    kiss_fft_free(cfg);
    KISS_FFT_FREE(inbuf);
    KISS_FFT_FREE(outbuf);
}

int main(void)
{
    test1();
    return 0;
}
