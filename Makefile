
all: kiss_fft_s kiss_fft_f kiss_fft_d freqpeak tones testsig

kiss_fft_s: kiss_fft.h kiss_fft.c
	gcc -Wall -O3 -o kiss_fft_s -DFIXED_POINT -DFFT_UTIL kiss_fft.c -lm

kiss_fft_f: kiss_fft.h kiss_fft.c
	gcc -Wall -O3 -o kiss_fft_f -Dkiss_fft_scalar=float -DFFT_UTIL kiss_fft.c -lm

kiss_fft_d: kiss_fft.h kiss_fft.c
	gcc -Wall -O3 -o kiss_fft_d -Dkiss_fft_scalar=double -DFFT_UTIL kiss_fft.c -lm
	
freqpeak: kiss_fft.h kiss_fft.c freqpeak.c
	gcc -Wall -O3 -o freqpeak freqpeak.c kiss_fft.c -lm

testsig: testsig.c
	gcc -Wall -O3 -o testsig testsig.c -lm

tones: tones.c
	gcc -Wall -O3 -o tones tones.c -lm
	
clean:
	rm -f kiss_fft_s kiss_fft_f kiss_fft_d *~ fftin.dat fftout.dat \
	freqpeak testsig tones

test: all
	./test.oct

speedtest: /dev/shm/junk kiss_fft_f
	time ./kiss_fft_f < /dev/shm/junk > /dev/null

/dev/shm/junk:
	dd if=/dev/urandom bs=8192 count=1000 of=/dev/shm/junk

tarball: clean
	tar -czf kiss_fft.tar.gz .
