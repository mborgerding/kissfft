
ifeq "$(NFFT)" ""
 NFFT=1024
endif

ifeq "$(NUMFFTS)" ""
 NUMFFTS=10000
endif

UTILSRC=../kiss_fft.c fftutil.c

CFLAGS=-Wall -O3 -lm -I.. -ansi -pedantic

ifeq "$(DATATYPE)" ""
  DATATYPE=float
endif
UTIL=fftutil_$(DATATYPE)


all: $(UTIL)

ifeq  "$(DATATYPE)" "short"
 TYPEFLAGS=-DFIXED_POINT -Dkiss_fft_scalar=short
else	
 TYPEFLAGS=-Dkiss_fft_scalar=$(DATATYPE)
endif

$(UTIL): $(UTILSRC)
	gcc -o $@ $(CFLAGS) $(TYPEFLAGS) $(UTILSRC)

RANDDAT=rand_$(DATATYPE)_$(NFFT).dat
$(RANDDAT):
	./rand_fft_data.py -n 1 -N $(NFFT) -t $(DATATYPE) > $(RANDDAT)

time: all $(RANDDAT)
	@echo
	@echo -n "#### timing $(NUMFFTS) x $(NFFT) point FFTs.  "; factor $(NFFT)
	@if [ -x ~/fftw/st ] && [ $(DATATYPE) == "double" ]; then \
		echo "#### FFTW FFT $(DATATYPE)"; \
		time ~/fftw/st -x $(NUMFFTS) -n $(NFFT) < $(RANDDAT) > /dev/null;\
	fi
	@echo "#### KISS FFT $(DATATYPE)"
	@time ./$(UTIL) -x $(NUMFFTS) -n $(NFFT) < $(RANDDAT) > /dev/null
	@rm $(RANDDAT)

POW2=256 512 1024 2048
POW3=243 729 2187
POW5=25 125 625 
mtime: all
	@for n in $(POW5) ;do \
		export NFFT=$$n;make time; \
		done

snr: all
	@which octave
	@echo
	@echo "### testing SNR for $(NFFT) point FFTs"
	@echo "#### $(DATATYPE)"
	@echo "testkiss( $(NFFT) , '$(DATATYPE)' );" | octave -q

test: snr time 

clean:
	rm -f *.dat *~ fftutil_*

