/*
Copyright (c) 2003, Mark Borgerding

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    * Neither the author nor the names of any contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>


#ifndef BUFSIZE
# define BUFSIZE 1024
#endif

int main(int argc, char ** argv)
{
    FILE * f1=stdin;
    FILE * f2=stdin;

    int i,n;
    TYPE buf1[BUFSIZE];
    TYPE buf2[BUFSIZE];

    double sigpower=0;
    double noisepower=0;
    double snrdb=0;
    double scale=0;
    long ntotal=0;

    if (argc>1 && strcmp( argv[1] ,"-") != 0 ) {
        f1 = fopen(argv[1],"rb");
    }
    if (argc>2 && strcmp( argv[2] ,"-") != 0 ) {
        f2 = fopen(argv[2],"rb");
    }

    // TODO LEFT OFF HERE
    while ( ( n = fread( buf1 , sizeof(TYPE) , BUFSIZE , f1 ) ) > 0 ) {
        if ( fread( buf2 , sizeof(TYPE) , BUFSIZE , f2 ) != n ) {
            fprintf(stderr,"premature end of file 2\n");
            exit(1);
        }
        for (i=0;i<n;++i) {
            double s=buf1[i];
            double n = s - buf2[i];
            sigpower += s*s;
            noisepower += n*n;
            if (s!=0) {
                ++ntotal;
                scale += buf2[i] / s;
            }
        }
    }
    
    if  ( fread( buf2 , sizeof(TYPE) , BUFSIZE , f2 ) > 0 ) {
        fprintf(stderr,"premature end of file 1\n");
        exit(1);
    }
    scale /= ntotal;
    
    if (noisepower>0)
        snrdb = 10*log10( sigpower / noisepower );
    else
        snrdb = 200;

    printf("snr = %.2f dB\n",snrdb);
    printf("average output/input = %.5e\n",scale);

    return 0;
}
