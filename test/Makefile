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
TESTKFC=tkfc_$(DATATYPE)
FASTFILT=ff_$(DATATYPE)
FASTFILTREAL=ffr_$(DATATYPE)

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

SRCFILES=../kiss_fft.c ../tools/kiss_fftnd.c ../tools/kiss_fftr.c pstats.c ../tools/kfc.c

all: tools $(BENCHKISS) $(SELFTEST) $(BENCHFFTW) $(TESTREAL) $(TESTKFC)

tools:
	cd ../tools && make all

CFLAGS=-Wall -O3 -pedantic -march=pentiumpro -ffast-math -fomit-frame-pointer   -I.. -I../tools
# If the above flags do not work, try the following
#CFLAGS=-Wall -O3 -I.. -I../tools

$(SELFTEST): $(SELFTESTSRC) $(SRCFILES)
	$(CC) -o $@ $(CFLAGS) $(TYPEFLAGS) -lm $+

$(TESTKFC): $(SRCFILES) 
	$(CC) -o $@ $(CFLAGS)  -DKFC_TEST $(TYPEFLAGS) -lm $+
	
$(TESTREAL): test_real.c $(SRCFILES)
	$(CC) -o $@ $(CFLAGS) $(TYPEFLAGS) -lm $+

$(BENCHKISS): benchkiss.c $(SRCFILES)
	$(CC) -o $@ $(CFLAGS) $(TYPEFLAGS) -lm  $+

$(BENCHFFTW): benchfftw.c pstats.c
	@echo "======attempting to build FFTW benchmark"
	@$(CC) -o $@ $(CFLAGS) -DDATATYPE$(DATATYPE) $+ -lm -l$(FFTWLIB) -L/usr/local/lib/ || echo "FFTW not available for comparison"

test: all
	@./$(TESTKFC)
	@echo "======1d & 2-d complex fft self test (type= $(DATATYPE) )"
	@./$(SELFTEST)
	@echo "======real FFT (type= $(DATATYPE) )"
	@./$(TESTREAL)
	@echo "======timing test (type=$(DATATYPE))"
	@./$(BENCHKISS) -x $(NUMFFTS) -n $(NFFT) 
	@[ -x ./$(BENCHFFTW) ] && ./$(BENCHFFTW) -x $(NUMFFTS) -n $(NFFT) ||true
	@echo "======higher dimensions type=$(DATATYPE))"
	@for dim in 2 3 4 5 6 7 8 9; do ./testkiss.py $$dim ;done

selftest.c:
	./mk_test.py 10 12 14 > selftest.c
selftest_short.c:
	./mk_test.py -s 10 12 14 > selftest_short.c

clean:
	rm -f *~ bm_* st_* tr_* kf_* tkfc_* ff_* ffr_* *.pyc *.pyo *.dat
