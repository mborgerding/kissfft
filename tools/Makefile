NFFT=1024
ALLUTILS=kfft kffts kfftd
NUMFFTS=10000
UTILSRC=../kiss_fft.c fftutil.c

all: $(ALLUTILS)

kfft: $(UTILSRC)
	gcc -Wall -O3 -o $@ -lm -I.. $(UTILSRC)
kffts: $(UTILSRC)
	gcc -Wall -O3 -o $@ -lm -I.. -DFIXED_POINT $(UTILSRC)
kfftd: $(UTILSRC)
	gcc -Wall -O3 -o $@ -lm -I.. -Dkiss_fft_scalar=double $(UTILSRC)

time: all
	@echo 
	@echo "#### timing $(NUMFFTS) x $(NFFT) point FFTs"
	@echo "#### DOUBLE"
	@time -f 'Elapsed:%E user:%U sys:%S' \
	 sh -c 'dd if=/dev/zero bs=$$((16*$(NFFT))) count=$(NUMFFTS) 2>/dev/null |./kfftd>/dev/null'
	@echo "#### FLOAT"
	@time -f 'Elapsed:%E user:%U sys:%S' \
	 sh -c 'dd if=/dev/zero bs=$$((8*$(NFFT))) count=$(NUMFFTS) 2>/dev/null |./kfft>/dev/null'
	@echo "#### SHORT"
	@time -f 'Elapsed:%E user:%U sys:%S' \
	 sh -c 'dd if=/dev/zero bs=$$((4*$(NFFT))) count=$(NUMFFTS) 2>/dev/null |./kffts>/dev/null'

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

