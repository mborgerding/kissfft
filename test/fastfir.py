#!/usr/local/bin/python2.3
from Numeric import *
from FFT import *

def make_random(len):
    import random
    res=[]
    for i in range(int(len)):
        r=random.uniform(-1,1)
        i=random.uniform(-1,1)
        res.append( complex(r,i) )
    return res

def slowfilter(sig,h):
    translen = len(h)-1
    return convolve(sig,h)[translen:-translen]

def nextpow2(x):
    return 2 ** math.ceil(math.log(x)/math.log(2))

def fastfilter(sig,h,nfft=None):
    if nfft is None:
        nfft = int( nextpow2( 2*len(h) ) )
    H = fft( h , nfft )
    scraplen = len(h)-1
    keeplen = nfft-scraplen
    res=[]
    isdone = 0
    lastidx = nfft
    idx0 = 0
    while not isdone:
        idx1 = idx0 + nfft
        if idx1 >= len(sig):
            idx1 = len(sig)
            lastidx = idx1-idx0
            if lastidx <= scraplen:
                break
            isdone = 1
        Fss = fft(sig[idx0:idx1],nfft)
        fm = Fss * H
        m = inverse_fft(fm)
        res.append( m[scraplen:lastidx] )
        idx0 += keeplen
    return concatenate( res )

def main():
    siglen = 1e5
    hlen = 500
    nfft = 1024*4
    print 'nfft=%d'%nfft
    # make a signal
    sig = make_random( siglen )
    # make an impulse response
    h = make_random( hlen )
    #h=[1]*2+[0]*3

    # perform MAC filtering
    yslow = slowfilter(sig,h)
    #print '<YSLOW>',yslow,'</YSLOW>'
    yfast = fastfilter(sig,h,nfft)
    #print yfast
    print 'len(yslow)=%d'%len(yslow)
    print 'len(yfast)=%d'%len(yfast)
    diff = yslow-yfast
    snr = 10*log10( vdot(yslow,yslow) / vdot(diff,diff) )
    print 'snr=%s' % snr


if __name__ == "__main__":
    main()
