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

    m = n/p
    Fm=[]
    for q in range(p): # 0,1
        fp = f[q::p]
        Fp = fft( fp )
        Fm.extend( Fp )

    Fout=[0]*n
    for k in range(n):
        val = 0
        for q in range(p):
            t = e ** ( j*2*pi*k*q/n )
            val += Fm[ q*m + (k%m) ] * t
        Fout[k] = val

    return Fout

def test(f=range(256),ntimes=10):
    import time
    t0 = time.time()
    for i in range(ntimes):
        fft(f)
    t1 = time.time()
    print '%ss'% (t1-t0)

