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

    #print 'n=%d,p=%d' % (n,p)
    #print f,' << fin'
    m = n/p
    Fout=[]
    for q in range(p): # 0,1
        fp = f[q::p]
        #print fp,'<< fp'
        Fp = fft( fp )
        Fout.extend( Fp )

    for u in range(m):
        scratch = Fout[u::m] # u to end in strides of m
        #print scratch
        for q1 in range(p):
            k = q1*m + u  # indices to Fout above that became scratch
            Fout[ k ] = scratch[0] # cuz e**0==1 in loop below
            for q in range(1,p):
                t = e ** ( j*2*pi*k*q/n )
                Fout[ k ] += scratch[q] * t

    return Fout

def real_fft( f ):
    N = len(f) / 2

    res = f[::2]
    ims = f[1::2]

    fp = [ complex(r,i) for r,i in zip(res,ims) ]

    Fp = fft( fp )

    F = []
    for k in range(N):
        s2 = ( Fp[k] + Fp[-k] ) * .5
        d2 = ( Fp[k] - Fp[-k] ).conjugate() * .5
        F1k = complex( s2.real , d2.imag )
        F2k = complex( s2.imag , d2.real )
        F.append( F1k +  e ** ( -j*pi*k/N ) * F2k )

    F.append( complex( Fp[0].real - Fp[0].imag , 0 ) ) 
    return F

def test(f=range(1024),ntimes=10):
    import time
    t0 = time.time()
    for i in range(ntimes):
        fft(f)
    t1 = time.time()
    print '%ss'% (t1-t0)

