#!/usr/local/bin/python2.3
import math
import sys
import random
import Numeric
import struct

pi=math.pi
e=math.e
j=complex(0,1)

def dopack(x,fmt='f',cpx=1):
    x = Numeric.reshape( x, ( Numeric.size(x),) )
    if cpx:
        s = ''.join( [ struct.pack(fmt*2,c.real,c.imag) for c in x ] )
    else:
        s = ''.join( [ struct.pack(fmt,c) for c in x ] )
    return s

def dounpack(x,fmt,cpx):
    uf = fmt * ( len(x) / 4 )
    s = struct.unpack(uf,x)
    if cpx:
        return Numeric.array(s[::2]) + Numeric.array( s[1::2] )*j
    else:
        return Numeric.array(s )


def main():
    #fft_func = fft
    fft_func = real_fft

    tvec = [0.309655,0.815653,0.768570,0.591841,0.404767,0.637617,0.007803,0.012665]
    Ftvec = [ complex(r,i) for r,i in zip(
                [3.548571,-0.378761,-0.061950,0.188537,-0.566981,0.188537,-0.061950,-0.378761],
                [0.000000,-1.296198,-0.848764,0.225337,0.000000,-0.225337,0.848764,1.296198] ) ]

    F = fft_func( tvec,0 )

    nerrs= 0
    for i in range(len(Ftvec)/2 + 1):
        if abs( F[i] - Ftvec[i] )> 1e-5:
            print 'F[%d]: %s != %s' % (i,F[i],Ftvec[i])
            nerrs += 1

    print '%d errors in forward fft' % nerrs
    if nerrs:
        return

    trec = fft_func( F , 1 )

    for i in range(len(trec) ):
        trec[i] /= len(trec)

    for i in range(len(tvec) ):
        if abs( trec[i] - tvec[i] )> 1e-5:
            print 't[%d]: %s != %s' % (i,tvec[i],trec[i])
            nerrs += 1

    print '%d errors in reverse fft' % nerrs


def make_random(dims=[1]):
    import Numeric 
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
    import Numeric
    ntotal = Numeric.product(Numeric.shape(x))
    return Numeric.reshape(x,(ntotal,))

def randmat( ndims ):
    dims=[]
    for i in range( ndims ):
        curdim = int( random.uniform(2,4) )
        dims.append( curdim )
    return make_random(dims )

def test_fftnd(ndims=3):
    import FFT
    import Numeric

    x=randmat( ndims )
    print 'dimensions=%s' % str( Numeric.shape(x) )
    #print 'x=%s' %str(x)
    xver = FFT.fftnd(x)
    x2=myfftnd(x)
    err = xver - x2
    errf = flatten(err)
    xverf = flatten(xver)
    errpow = Numeric.vdot(errf,errf)+1e-10
    sigpow = Numeric.vdot(xverf,xverf)+1e-10
    snr = 10*math.log10(abs(sigpow/errpow) )
    print 'SNR(compared to Python FFT module) =%sdB' % str( snr )
    if snr<80:
        print xver
        print x2
        sys.exit(1)
 
def myfftnd(x):
    import Numeric
    xf = flatten(x)
    Xf = fftndwork( xf , Numeric.shape(x) )
    return Numeric.reshape(Xf,Numeric.shape(x) )

def fftndwork(x,dims):
    import popen2

    cmd = '../tools/fft -n '
    cmd += ','.join([str(d) for d in dims])
    p = popen2.Popen3(cmd )
    p.tochild.write( dopack( x , 'f' ,1 ) )
    p.tochild.close()
    res = dounpack( p.fromchild.read() , 'f' ,1 )
    p.wait()
    return res

    #import Numeric
    #dimprod=Numeric.product( dims )
#
    #for k in range( len(dims) ):
        #cur_dim=dims[ k ]
        #stride=dimprod/cur_dim
        #next_x = [complex(0,0)]*len(x)
        #for i in range(stride):
            #next_x[i*cur_dim:(i+1)*cur_dim] = fft(x[i:(i+cur_dim)*stride:stride],0)
        #x = next_x
    #return x

if __name__ == "__main__":
    try:
        nd = int(sys.argv[1])
    except:
        nd=None
    if nd:    
        test_fftnd( nd )
    else:    
        sys.exit(0)
