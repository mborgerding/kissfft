WARNINGS=-W -Wall -Wstrict-prototypes -Wmissing-prototypes -Waggregate-return \
    -Wcast-align -Wcast-qual -Wnested-externs -Wshadow -Wbad-function-cast \
    -Wwrite-strings

ifeq "$(DATATYPE)" ""
 DATATYPE=float
endif


ifeq  "$(DATATYPE)" "short"
 TYPEFLAGS=-DFIXED_POINT -Dkiss_fft_scalar=short
else
 TYPEFLAGS=-Dkiss_fft_scalar=$(DATATYPE)
endif

ifeq  "$(DATATYPE)" "float"
 FFTUTIL=fft
 FASTFILT=fastconv
 FASTFILTREAL=fastconvr
 PSDPNG=psdpng
 DUMPHDR=dumphdr
else
 FFTUTIL=fft_$(DATATYPE)
 FASTFILT=fastconv_$(DATATYPE)
 FASTFILTREAL=fastconvr_$(DATATYPE)
 PSDPNG=psdpng_$(DATATYPE)
 DUMPHDR=dumphdr_$(DATATYPE)
endif

all: $(FFTUTIL) $(FASTFILT) $(FASTFILTREAL) $(PSDPNG) $(DUMPHDR)

#CFLAGS=-Wall -O3 -pedantic -march=pentiumpro -ffast-math -fomit-frame-pointer $(WARNINGS)
# If the above flags do not work, try the following
CFLAGS=-Wall -O3 $(WARNINGS)

$(FASTFILTREAL): ../kiss_fft.c kiss_fastfir.c kiss_fftr.c
	$(CC) -o $@ $(CFLAGS) -I.. $(TYPEFLAGS) -DREAL_FASTFIR -lm $+ -DFAST_FILT_UTIL

$(FASTFILT): ../kiss_fft.c kiss_fastfir.c
	$(CC) -o $@ $(CFLAGS) -I.. $(TYPEFLAGS) -lm $+ -DFAST_FILT_UTIL

$(FFTUTIL): ../kiss_fft.c fftutil.c kiss_fftnd.c kiss_fftr.c
	$(CC) -o $@ $(CFLAGS) -I.. $(TYPEFLAGS) -lm $+

$(PSDPNG): ../kiss_fft.c psdpng.c kiss_fftr.c
	$(CC) -o $@ $(CFLAGS) -I.. $(TYPEFLAGS) -lm -lpng $+

$(DUMPHDR): ../kiss_fft.c dumphdr.c
	$(CC) -o $@ $(CFLAGS) -I.. $(TYPEFLAGS) -lm $+

clean:
	rm -f *~ fft fft_* fastconv fastconv_* fastconvr fastconvr_* psdpng psdpng_*
