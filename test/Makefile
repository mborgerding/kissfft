ifeq "$(NFFT)" ""
 NFFT=1800
endif
ifeq "$(NUMFFTS)" ""
 NUMFFTS=10000
endif
ifeq "$(DATATYPE)" ""
 DATATYPE=float
endif

UTIL=fftutil_$(DATATYPE)
BENCH=bm_$(DATATYPE)

all: $(UTIL) $(BENCH)

ifeq  "$(DATATYPE)" "short"
 TYPEFLAGS=-DFIXED_POINT -Dkiss_fft_scalar=short
else
 TYPEFLAGS=-Dkiss_fft_scalar=$(DATATYPE)
endif

CFLAGS=-Wall -O3 -ansi -pedantic
$(UTIL): ../kiss_fft.c fftutil.c
	gcc -o $@ $(CFLAGS) -I.. $(TYPEFLAGS) ../kiss_fft.c fftutil.c -lm

$(BENCH): benchkiss.c ../kiss_fft.c
	gcc -o $@ $(CFLAGS) -I.. benchkiss.c $(TYPEFLAGS) ../kiss_fft.c -lm 

fftw: bm_fftw
	@[ -x ./bm_fftw ] && \
		./bm_fftw -x $(NUMFFTS) -n $(NFFT)

bm_fftw: benchfftw.c
	@gcc -o $@ $(CFLAGS) benchfftw.c -lm -lfftw3 -L /usr/local/lib/ \
	 || echo 'Cannot build FFTW test script'

time: all
	@./$(BENCH) -x $(NUMFFTS) -n $(NFFT) 

POW2=256 512 1024 2048
POW3=243 729 2187
POW5=25 125 625 3125

mtime: all bm_fftw
	for n in $(POW3) ;do \
		./$(BENCH) -x $(NUMFFTS) -n $$n;\
		[ "$(DATATYPE)" == "double" ] && [ -x ./bm_fftw ] && ./bm_fftw -x $(NUMFFTS) -n $$n || true ; \
	done

snr: all
	@echo "### testing SNR for $(NFFT) point $(DATATYPE) FFTs"
	@echo "testkiss( $(NFFT) , '$(DATATYPE)' );" | octave -q


ifeq "$(DATATYPE)" "double"
test: snr time fftw
else
test: snr time 
endif

clean:
	rm -f *~ fftutil_* bm_*
