
all: kiss_fft_s kiss_fft_f kiss_fft_d

kiss_fft_s: kiss_fft.h kiss_fft.c
	gcc -Wall -O3 -o kiss_fft_s -DFIXED_POINT -DFFT_UTIL kiss_fft.c -lm

kiss_fft_f: kiss_fft.h kiss_fft.c
	gcc -Wall -O3 -o kiss_fft_f -Dkiss_fft_scalar=float -DFFT_UTIL kiss_fft.c -lm

kiss_fft_d: kiss_fft.h kiss_fft.c
	gcc -Wall -O3 -o kiss_fft_d -Dkiss_fft_scalar=double -DFFT_UTIL kiss_fft.c -lm

clean:
	rm -f kiss_fft_s kiss_fft_f kiss_fft_d *~ fftin.dat fftout.dat 

test: clean all
	./test.oct

tarball: clean
	tar -czf kiss_fft.tar.gz .
