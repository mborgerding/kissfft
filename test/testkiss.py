#!/usr/local/bin/python2.3
import math
import sys
import os
import random
import Numeric
import FFT
import struct
import popen2

pi=math.pi
e=math.e
j=complex(0,1)

doreal=0

datatype = os.environ.get('DATATYPE','float')

util = '../tools/fft'
fmt='f'
minsnr=90
if datatype == 'double':
    util = '../tools/fft_double'
    fmt='d'
elif datatype=='short':
    util = '../tools/fft_short'
    fmt='h'
    minsnr=10

def dopack(x,cpx=1):
    x = Numeric.reshape( x, ( Numeric.size(x),) )
    if cpx:
        s = ''.join( [ struct.pack(fmt*2,c.real,c.imag) for c in x ] )
    else:
        s = ''.join( [ struct.pack(fmt,c) for c in x ] )
    return s

def dounpack(x,cpx):
    uf = fmt * ( len(x) / struct.calcsize(fmt) )
    s = struct.unpack(uf,x)
    if cpx:
        return Numeric.array(s[::2]) + Numeric.array( s[1::2] )*j
    else:
        return Numeric.array(s )

def make_random(dims=[1]):
    res = []
    for i in range(dims[0]):
        if len(dims)==1:
            r=random.uniform(-1,1)
            i=random.uniform(-1,1)
            res.append( complex(r,i) )
        else:
            res.append( make_random( dims[1:] ) )
    return Numeric.array(res)

def flatten(x):
    ntotal = Numeric.product(Numeric.shape(x))
    return Numeric.reshape(x,(ntotal,))

def randmat( ndims ):
    dims=[]
    for i in range( ndims ):
        curdim = int( random.uniform(2,4) )
        dims.append( curdim )
    return make_random(dims )

def test_fft(ndims):
    if ndims == 1:
        x=Numeric.array(make_random( [ int(random.uniform(500,1025)) ] ))
    else:
        x=randmat( ndims )

    xver = FFT.fftnd(x)
    x2=dofft(x)
    err = xver - x2
    errf = flatten(err)
    xverf = flatten(xver)
    errpow = Numeric.vdot(errf,errf)+1e-10
    sigpow = Numeric.vdot(xverf,xverf)+1e-10
    snr = 10*math.log10(abs(sigpow/errpow) )
    print 'dimensions=%s' % str( Numeric.shape(x) ),
    print 'SNR (compared to NumPy) : %.1fdB' % float(snr)

    if snr<minsnr:
        print 'xver=',xver
        print 'x2=',x2
        print 'err',err
        sys.exit(1)
 
def dofft(x):
    dims=Numeric.shape(x)
    x = flatten(x)
    iscomp = (type(x[0]) == complex)

    scale=1
    if datatype=='short':
        x = 32767 * x
        scale = len(x) / 32767.0

    cmd='%s -n ' % util
    cmd += ','.join([str(d) for d in dims])

    p = popen2.Popen3(cmd )

    p.tochild.write( dopack( x , iscomp ) )
    p.tochild.close()

    res = dounpack( p.fromchild.read() , iscomp )

    res = scale * res

    p.wait()
    return Numeric.reshape(res,dims)

def main():
    for dim in range(1,9):
        test_fft( dim )
    print 'We crossed the 8th dimension.  Buckaroo would be proud'


if __name__ == "__main__":
    main()

