#!/usr/local/bin/python2.3
import math
pi=math.pi
e=math.e
j=complex(0,1)

def T(n,k):
    return e ** -2*j*pi*k/n

def fft(f):
    n=len(f)
    if n==1:
        return f

    for p in 4,3,2,5:
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
            t = e ** ( -j*2*pi*k*q/n )
            val += Fm[ q*m + (k%m) ] * t
        Fout[k] = val

    return Fout

