#!/usr/local/bin/python2.3
import math

def T(n,i):
    return math.e ** complex( 0,-2*math.pi*i/n )

def fft(f):
    n=len(f)
    if n%2 == 0:
        np = n/2
        fe=[0 ] * np
        fo=[0 ] * np
        for i in range(np):
            fe[i] = f[i] + f[i+np]
            fo[i] = (f[i] - f[i+np]) * T(n,i)
        Fe=fft(fe)
        Fo=fft(fo)
        F=[0 ] * n
        F[::2] = Fe
        F[1::2] = Fo
    elif n==1:
        F=f
    else:
        raise exceptions.Exception('cannot factor %d by 2' % n)
    return F

