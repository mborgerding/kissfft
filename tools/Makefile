ifeq "$(NFFT)" ""
 NFFT=1800
endif
ifeq "$(NUMFFTS)" ""
 NUMFFTS=10000
endif
NROWS=30

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

all: $(BENCHKISS) $(SELFTEST) $(BENCHFFTW) $(TESTREAL)

#CFLAGS=-Wall -O3 -ansi -pedantic -march=pentiumpro -ffast-math -fomit-frame-pointer
# If the above flags do not work, try the following
CFLAGS=-Wall -O3

$(SELFTEST): ../kiss_fft.c $(SELFTESTSRC)
	$(CC) -o $@ $(CFLAGS) -I.. $(TYPEFLAGS) -lm $+

$(TESTREAL): ../kiss_fft.c ../kiss_fftr.c test_real.c
	$(CC) -o $@ $(CFLAGS) -I.. $(TYPEFLAGS) -lm $+

$(BENCHKISS): benchkiss.c ../kiss_fft.c pstats.c
	$(CC) -o $@ $(CFLAGS) -I.. benchkiss.c $(TYPEFLAGS) ../kiss_fft.c pstats.c -lm 

$(BENCHFFTW): benchfftw.c pstats.c
	@$(CC) -o $@ $(CFLAGS) -DDATATYPE$(DATATYPE) benchfftw.c pstats.c -lm -lfftw3f -lfftw3 -L /usr/local/lib/ || echo "FFTW not available for comparison"

POW2=256 512 1024 2048 4096 8192
POW3=243 729 2187
POW5=25 125 625 3125

mtime: all 
	@for n in $(POW2) $(POW3) $(POW5) ;do \
		echo ============================;\
		./$(BENCHKISS) -x $(NUMFFTS) -n $$n;\
		[ -x ./$(BENCHFFTW) ] && ./$(BENCHFFTW) -x $(NUMFFTS) -n $$n || true ; \
	done

test: all $(TESTREAL)
	@echo "======SELF TEST"
	@./$(SELFTEST)
	@echo "======REAL FFT TEST"
	@./$(TESTREAL)
	@echo "======TIMING TEST"
	@./$(BENCHKISS) -x $(NUMFFTS) -n $(NFFT) 
	@[ -x ./$(BENCHFFTW) ] && ./$(BENCHFFTW) -x $(NUMFFTS) -n $(NFFT) ||true

selftest.c:
	./mk_test.py 10 12 14 > selftest.c
selftest_short.c:
	./mk_test.py -s 10 12 14 > selftest_short.c

clean:
	rm -f *~ bm_* st_* tr_* *.dat
