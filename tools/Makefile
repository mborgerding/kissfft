
ifeq "$(NFFT)" ""
	NFFT=1024
endif

ALLUTILS=kfft kffts kfftd

ifeq "$(NUMFFTS)" ""
	NUMFFTS=10000
endif

UTILSRC=../kiss_fft.c fftutil.c

CFLAGS=-Wall -O3 -lm -I.. -ansi -pedantic

all: $(ALLUTILS)

kfft: $(UTILSRC)
	echo $(FUNCDECL)
	gcc -o $@ $(CFLAGS) $(UTILSRC)
kffts: $(UTILSRC)
	gcc -o $@ $(CFLAGS) -DFIXED_POINT $(UTILSRC)
kfftd: $(UTILSRC)
	gcc -o $@ $(CFLAGS) -Dkiss_fft_scalar=double $(UTILSRC)

time: all
	@echo 
	@echo -n "#### timing $(NUMFFTS) x $(NFFT) point FFTs.  "; factor $(NFFT)
	@echo "#### DOUBLE"
	@time -f 'Elapsed:%E user:%U sys:%S' \
	 sh -c 'dd if=/dev/zero bs=$$((16*$(NFFT))) count=$(NUMFFTS) 2>/dev/null |./kfftd>/dev/null'
	@echo "#### FLOAT"
	@time -f 'Elapsed:%E user:%U sys:%S' \
	 sh -c 'dd if=/dev/zero bs=$$((8*$(NFFT))) count=$(NUMFFTS) 2>/dev/null |./kfft>/dev/null'
	@echo "#### SHORT"
	@time -f 'Elapsed:%E user:%U sys:%S' \
	 sh -c 'dd if=/dev/zero bs=$$((4*$(NFFT))) count=$(NUMFFTS) 2>/dev/null |./kffts>/dev/null'


POW2=256 512 1024 2048
POW3=243 729 2187
mtime: all
	@for n in $(POW2) ;do \
		export NFFT=$$n;make time; \
		done

snr: all
	@which octave
	@echo
	@echo "### testing SNR for  $(NFFT) point FFTs"
	@echo "#### DOUBLE"
	@echo "testkiss($(NFFT),'double',1,1/$(NFFT));" | octave -q
	@echo "#### FLOAT"
	@echo "testkiss($(NFFT),'float',1,1/$(NFFT));" | octave -q
	@echo "#### SHORT"
	@echo "testkiss($(NFFT),'short',$(NFFT),$(NFFT));" | octave -q

test: snr time 

clean:
	rm -f $(ALLUTILS) *.dat *~

