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
BENCHKISS=bm_kiss_$(DATATYPE)
BENCHFFTW=bm_fftw_$(DATATYPE)

all: $(UTIL) $(BENCHKISS)

ifeq  "$(DATATYPE)" "short"
 TYPEFLAGS=-DFIXED_POINT -Dkiss_fft_scalar=short
else
 TYPEFLAGS=-Dkiss_fft_scalar=$(DATATYPE)
endif

CFLAGS=-Wall -O3 -ansi -pedantic -march=pentiumpro -ffast-math -fomit-frame-pointer
# If the above flags do not work, try the following
#CFLAGS=-Wall -O3 -ansi -pedantic

$(UTIL): ../kiss_fft.c ../kiss_fft2d.c fftutil.c
	gcc -o $@ $(CFLAGS) -I.. $(TYPEFLAGS) ../kiss_fft.c ../kiss_fft2d.c fftutil.c -lm

$(BENCHKISS): benchkiss.c ../kiss_fft.c pstats.c
	gcc -o $@ $(CFLAGS) -I.. benchkiss.c $(TYPEFLAGS) ../kiss_fft.c pstats.c -lm 

fftw: $(BENCHFFTW)
	./$(BENCHFFTW) -x $(NUMFFTS) -n $(NFFT)

$(BENCHFFTW): benchfftw.c pstats.c
	gcc -o $@ $(CFLAGS) -DDATATYPE$(DATATYPE) benchfftw.c pstats.c -lm -lfftw3f -lfftw3 -L /usr/local/lib/ 

time: all
	@./$(BENCHKISS) -x $(NUMFFTS) -n $(NFFT) 

POW2=256 512 1024 2048 4096 8192
POW3=243 729 2187
POW5=25 125 625 3125

mtime: all $(BENCHFFTW)
	for n in $(POW2) $(POW3) $(POW5) ;do \
		echo ============================;\
		./$(BENCHKISS) -x $(NUMFFTS) -n $$n;\
		[ -x ./$(BENCHFFTW) ] && ./$(BENCHFFTW) -x $(NUMFFTS) -n $$n || true ; \
	done

snr: all
	@echo "### testing SNR for $(NFFT) point $(DATATYPE) FFTs"
	@echo "testkiss( $(NFFT) , '$(DATATYPE)' );" | octave -q

test: snr time fftw

clean:
	rm -f *~ fftutil_* bm_* *.dat
