ifeq "$(NFFT)" ""
 NFFT=1800
endif
ifeq "$(NUMFFTS)" ""
 NUMFFTS=10000
endif

ifeq "$(DATATYPE)" ""
 DATATYPE=float
endif

BENCHKISS=bm_kiss_$(DATATYPE)
BENCHFFTW=bm_fftw_$(DATATYPE)
SELFTEST=st_$(DATATYPE)
TESTREAL=tr_$(DATATYPE)

ifeq  "$(DATATYPE)" "short"
 TYPEFLAGS=-DFIXED_POINT -Dkiss_fft_scalar=short
 SELFTESTSRC=selftest_short.c
else
 TYPEFLAGS=-Dkiss_fft_scalar=$(DATATYPE)
 SELFTESTSRC=selftest.c
endif

ifeq  "$(DATATYPE)" "float"
# fftw needs to be built with --enable-float to build this lib
 FFTWLIB=fftw3f
else
 FFTWLIB=fftw3
endif


all: $(BENCHKISS) $(SELFTEST) $(BENCHFFTW) $(TESTREAL)

#CFLAGS=-Wall -O3 -ansi -pedantic -march=pentiumpro -ffast-math -fomit-frame-pointer
# If the above flags do not work, try the following
CFLAGS=-Wall -O3

$(SELFTEST): ../kiss_fft.c $(SELFTESTSRC)
	$(CC) -o $@ $(CFLAGS) -I.. $(TYPEFLAGS) -lm $+

$(TESTREAL): ../kiss_fft.c kiss_fftr.c test_real.c
	$(CC) -o $@ $(CFLAGS) -I.. $(TYPEFLAGS) -lm $+

$(BENCHKISS): benchkiss.c ../kiss_fft.c pstats.c
	$(CC) -o $@ $(CFLAGS) -I.. benchkiss.c $(TYPEFLAGS) ../kiss_fft.c pstats.c -lm 

$(BENCHFFTW): benchfftw.c pstats.c
	@$(CC) -o $@ $(CFLAGS) -DDATATYPE$(DATATYPE) benchfftw.c pstats.c -lm -l$(FFTWLIB) -L/usr/local/lib/ || echo "FFTW not available for comparison"

test: all
	@echo "======SELF TEST $(DATATYPE)"
	@./$(SELFTEST)
	@echo "======REAL FFT TEST $(DATATYPE)"
	@./$(TESTREAL)
	@echo "======TIMING TEST $(DATATYPE)"
	@./$(BENCHKISS) -x $(NUMFFTS) -n $(NFFT) 
	@[ -x ./$(BENCHFFTW) ] && ./$(BENCHFFTW) -x $(NUMFFTS) -n $(NFFT) ||true

testall:
	@(export DATATYPE=double && make test )
	@(export DATATYPE=float && make test )
	@(export DATATYPE=short && make test )

selftest.c:
	./mk_test.py 10 12 14 > selftest.c
selftest_short.c:
	./mk_test.py -s 10 12 14 > selftest_short.c

clean:
	rm -f *~ bm_* st_* tr_*
