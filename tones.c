
#include <stdio.h>
#include <string.h>
#include <memory.h>
#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include "kiss_fft.h"

#define PI 3.14159

int main(int argc, char ** argv)
{
    int k;
    float fs=44100;

    float fr=1,fl=300;

    float th[2] = {0,0};
    float thinc[2] = {2*PI*fr/fs,2*PI*fl/fs };

    while (1){
        for (k=0;k<2;++k){
            short s;
            th[k] += thinc[k];
            if (th[k] > 2*PI){
                th[k] -= 2*PI;
            }
            s=(short)32767*cos( th[k] );
            fwrite(&s,sizeof(s),1,stdout);
        }
    }

    return 0;
}
