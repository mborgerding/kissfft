#!/usr/local/bin/python2.3

import sys
import random
import struct

def main():
    from getopt import getopt
    opts,args = getopt(sys.argv[1:],'n:N:t:l:h:')
    opts=dict(opts)
    nbufs = int( opts.get('-n','10000') )
    nfft = int( opts.get('-N','1024') )
    type = opts.get('-t','f')
    lo = float(opts.get('-l','-32768') )
    hi = float(opts.get('-h','32767') )

    format = {'float':'f','short':'h','double':'d'}[type]

    nums = [ random.uniform(lo,hi) for i in range(nfft*2) ]
    buf = struct.pack( '%d%s' % ( len( nums ) , format ) , *nums )

    for i in range(nbufs):
        sys.stdout.write( buf )

if __name__ == '__main__':
    main()
