#!/usr/local/bin/python2.3
import math
pi=math.pi
e=math.e
j=complex(0,1)

def fft(f):
    n=len(f)
    if n==1:
        return f

    for p in 2,3,5:
        if n%p==0:
            break
    else:
        raise Exception('%s not factorable ' % n)

    print 'n=%d,p=%d' % (n,p)
    print f,' << fin'
    m = n/p
    Fout=[]
    for q in range(p): # 0,1
        fp = f[q::p]
        print fp,'<< fp'
        Fp = fft( fp )
        Fout.extend( Fp )

    for u in range(m):
        scratch = Fout[u::m] # u to end in strides of m
        print scratch
        for q1 in range(p):
            k = q1*m + u  # indices to Fout above that became scratch
            Fout[ k ] = scratch[0] # cuz e**0==1 in loop below
            for q in range(1,p):
                t = e ** ( j*2*pi*k*q/n )
                Fout[ k ] += scratch[q] * t

    return Fout

def test(f=range(1024),ntimes=10):
    import time
    t0 = time.time()
    for i in range(ntimes):
        fft(f)
    t1 = time.time()
    print '%ss'% (t1-t0)

