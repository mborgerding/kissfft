ifeq "$(DATATYPE)" ""
 DATATYPE=float
endif

FFTUTIL=fft_$(DATATYPE)
FASTFILT=fastconv_$(DATATYPE)
FASTFILTREAL=fastconvr_$(DATATYPE)

ifeq  "$(DATATYPE)" "short"
 TYPEFLAGS=-DFIXED_POINT -Dkiss_fft_scalar=short
else
 TYPEFLAGS=-Dkiss_fft_scalar=$(DATATYPE)
endif


all: $(FFTUTIL) $(FASTFILT) $(FASTFILTREAL)

CFLAGS=-Wall -O3 -pedantic -march=pentiumpro -ffast-math -fomit-frame-pointer 
# If the above flags do not work, try the following
#CFLAGS=-Wall -O3

$(FASTFILTREAL): ../kiss_fft.c kiss_fastfir.c kiss_fftr.c
	$(CC) -o $@ $(CFLAGS) -I.. $(TYPEFLAGS) -DREAL_FASTFIR -lm $+ -DFAST_FILT_UTIL

$(FASTFILT): ../kiss_fft.c kiss_fastfir.c
	$(CC) -o $@ $(CFLAGS) -I.. $(TYPEFLAGS) -lm $+ -DFAST_FILT_UTIL

$(FFTUTIL): ../kiss_fft.c fftutil.c kiss_fftnd.c kiss_fftr.c
	$(CC) -o $@ $(CFLAGS) -I.. $(TYPEFLAGS) -lm $+

clean:
	rm -f *~ fft_* fastconv_* fastconvr_*
