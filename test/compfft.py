#!/usr/local/bin/python2.3

# use FFTPACK as a baseline
import FFT
from Numeric import *
import math
import random
import sys
import struct

pi=math.pi
e=math.e
j=complex(0,1)
lims=(-32768,32767)

def randbuf(n,cpx=1):
    res = array( [ random.uniform( lims[0],lims[1] ) for i in range(n) ] )
    if cpx:
        res = res + j*randbuf(n,0)
    return res

def main():
    from getopt import getopt
    import popen2
    opts,args = getopt( sys.argv[1:],'u:n:Rt:' )
    opts=dict(opts)

    util = opts.get('-u','./kf_float')

    try:
        n = int(opts['-n'])
        cpx = opts.get('-R') is None
        fmt=opts.get('-t','f')
    except KeyError:
        sys.stderr.write("""
        usage: compfft.py 
        -n nfft
        -u utilname : see sample_code/fftutil.c
        -R : real-optimized version
        """)

    x = randbuf(n,cpx)

    cmd = '%s -n %d ' % ( util, n )
    if cpx:
        xout = FFT.fft(x)
    else:        
        cmd += '-R '
        xout = FFT.real_fft(x)

    proc = popen2.Popen3( cmd , bufsize=len(x) )

    proc.tochild.write( dopack( x , fmt ,cpx ) )
    proc.tochild.close()
    xoutcomp = dounpack( proc.fromchild.read( ) , fmt ,1 )

    sigpow = sum( xout * conjugate(xout) )
    print len(xout)
    print len(xoutcomp)
    
    diff = xout-xoutcomp
    noisepow = sum( diff * conjugate(diff) )

    snr = 10 * math.log10(abs( sigpow / noisepow ) )
    print 'NFFT=%d,SNR = %f dB' % (n,snr)

def dopack(x,fmt,cpx):
    if cpx:
        s = ''.join( [ struct.pack('ff',c.real,c.imag) for c in x ] )
    else:
        s = ''.join( [ struct.pack('f',c) for c in x ] )
    return s 

def dounpack(x,fmt,cpx):
    uf = fmt * ( len(x) / 4 )
    s = struct.unpack(uf,x)
    if cpx:
        return array(s[::2]) + array( s[1::2] )*j
    else:    
        return array(s )

if __name__ == "__main__":
    main()
