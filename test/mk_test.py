#!/usr/local/bin/python2.3
import FFT
import sys
import random
import re
j=complex(0,1)

def randvec(n,iscomplex):
    if iscomplex:
        return [
                int(random.uniform(-32768,32767) ) + j*int(random.uniform(-32768,32767) )
                for i in range(n) ]
    else:                
        return [ int(random.uniform(-32768,32767) ) for i in range(n) ]
    
def c_format(v,round=0):
    if round:
        return ','.join( [ '{%d,%d}' %(int(c.real),int(c.imag) ) for c in v ] ) 
    else:
        s= ','.join( [ '{%.60f ,%.60f }' %(c.real,c.imag) for c in v ] ) 
        return re.sub(r'\.?0+ ',' ',s)

def test_vec( v,inverse ):
    if inverse:
        tvecout = FFT.inverse_fft(v)
        tvecout = [ c * len(v) for c in tvecout ]
    else:
        tvecout = FFT.fft(v)

    s="""#define NFFT %d""" % len(v) + """
    {
        double snr;
        kiss_fft_cpx test_vec_in[NFFT] = { """  + c_format(v) + """};
        kiss_fft_cpx test_vec_out[NFFT] = {"""  + c_format( tvecout ) + """};
        kiss_fft_cpx testbuf[NFFT];
        void * cfg = kiss_fft_alloc(NFFT,%d);""" % inverse + """

        kiss_fft(cfg,test_vec_in,testbuf);

        snr = snr_compare(test_vec_out,testbuf,NFFT);
        printf("FFT n=%d, inverse=%d, snr = %g dB\\n",NFFT,""" + str(inverse) + """,snr);
        free(cfg);
    }
#undef NFFT    
"""
    return s

def compare_func():
    s="""
double snr_compare( kiss_fft_cpx * test_vec_out,kiss_fft_cpx * testbuf, int n)
{
    int k;
    double sigpow,noisepow,err,snr,scale=0;
    sigpow = noisepow = .000000000000000000000000000001; 

    for (k=0;k<n;++k) {
        sigpow += test_vec_out[k].r * test_vec_out[k].i + 
                  test_vec_out[k].i * test_vec_out[k].i;
        err = test_vec_out[k].r - testbuf[k].r;
        noisepow += err * err;
        err = test_vec_out[k].i - testbuf[k].i;
        noisepow += err * err;

        if (test_vec_out[k].r)
            scale += testbuf[k].r / test_vec_out[k].r;
    }
    snr = 10*log10( sigpow / noisepow );
    scale /= n;
    if (snr<10)
        printf( "\\npoor snr, try a scaling factor %f\\n" , scale );
    return snr;
}
"""
    return s

def main():
    fftsizes = sys.argv[1:]
    if not fftsizes:
        fftsizes = [ 1800 ]
    print '#include "kiss_fft.h"'
    print compare_func()
    print "int main() {"
    for n in fftsizes:
        v = randvec(int(n),1)
        print test_vec(v,0)
        print test_vec(v,1)
    print """
    return 0;
}
"""

if __name__ == "__main__":
    main()
