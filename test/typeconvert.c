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


#define BUFSIZE 1024
int main(int argc, char ** argv )
{
    int i,n;
    TYPE1 buf1[BUFSIZE];
    TYPE2 buf2[BUFSIZE];
    double mult=1.0;

    for (i=1;i<argc;++i){
        if (strcmp("-m",argv[i])==0){
            mult *= atof(argv[i+1] );
            ++i;
        }else if (strcmp("-d",argv[i])==0){
            mult /= atof(argv[i+1] );
            ++i;
        }
    }

    while ( ( n = fread( buf1 , sizeof(TYPE1) , BUFSIZE , stdin ) ) > 0 ) {
        for (i=0;i<n;++i)
            buf2[i] = (TYPE1)(mult * buf1[i]);
        fwrite( buf2 , sizeof(TYPE2) , n , stdout );
    }

    return 0;
}
