SPEEDTESTFILE=/dev/shm/kissfft_speedtest
SPEEDTESTNSAMPS=1000

all: kiss_fft_s kiss_fft_f kiss_fft_d

kiss_fft_s: kiss_fft.h kiss_fft.c
	gcc -Wall -O3 -o kiss_fft_s -DFIXED_POINT -DFFT_UTIL kiss_fft.c -lm

kiss_fft_f: kiss_fft.h kiss_fft.c
	gcc -Wall -O3 -o kiss_fft_f -Dkiss_fft_scalar=float -DFFT_UTIL kiss_fft.c -lm

kiss_fft_d: kiss_fft.h kiss_fft.c
	gcc -Wall -O3 -o kiss_fft_d -Dkiss_fft_scalar=double -DFFT_UTIL kiss_fft.c -lm
	
clean:
	rm -f kiss_fft_s kiss_fft_f kiss_fft_d *~ fftin.dat fftout.dat $(SPEEDTESTFILE)

test: all
	./test.oct

speedf: kiss_fft_f  $(SPEEDTESTFILE)
	time ./kiss_fft_f < $(SPEEDTESTFILE) > /dev/null

$(SPEEDTESTFILE):
	dd if=/dev/zero bs=8192 count=$(SPEEDTESTNSAMPS) of=$(SPEEDTESTFILE)

tarball: clean
	tar -czf kiss_fft.tar.gz .
