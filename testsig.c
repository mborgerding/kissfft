#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int usage()
{
    fprintf(stderr,"usage:testsig nsamps\n");
    exit(1);
    return 1;
}

double randphase()
{
    return (double)rand()*2*3.14159/RAND_MAX;
}

int main(int argc, char ** argv)
{
    float samps[2];
    int nsamps;

    if (argc != 2)
        return usage();
    nsamps = atoi( argv[1] );
    
    while (nsamps-- > 0) {
        samps[0]=sin( randphase() );
        samps[1]=sin( randphase() );
        fwrite(samps,sizeof(samps),1,stdout);
    }
    return 0;
}
