NFFT=256
ALLUTILS=kfft kifft kffts kiffts kfftd kifftd s2f s2d sigdiff sigdiffs sigdiffd scale scales scaled 
NUMFFTS=1000
UTILSRC=../kiss_fft.c fftutil.c

all: $(ALLUTILS)

kfft: $(UTILSRC)
	gcc -Wall -O3 -o $@ -lm -I.. -DNFFT=$(NFFT) $(UTILSRC)
kifft: ../kiss_fft.h ../kiss_fft.c fftutil.c
	gcc -Wall -O3 -o $@ -lm -I.. -DNFFT=$(NFFT) -DINVERSE_FFT $(UTILSRC)
kffts: $(UTILSRC)
	gcc -Wall -O3 -o $@ -lm -I.. -DNFFT=$(NFFT) -DFIXED_POINT $(UTILSRC)
kiffts: $(UTILSRC)
	gcc -Wall -O3 -o $@ -lm -I.. -DNFFT=$(NFFT) -DFIXED_POINT -DINVERSE_FFT $(UTILSRC)
kfftd: $(UTILSRC)
	gcc -Wall -O3 -o $@ -lm -I.. -DNFFT=$(NFFT) -Dkiss_fft_scalar=double $(UTILSRC)
kifftd: ../kiss_fft.h ../kiss_fft.c fftutil.c
	gcc -Wall -O3 -o $@ -lm -I.. -DNFFT=$(NFFT) -Dkiss_fft_scalar=double -DINVERSE_FFT $(UTILSRC)

s2f: typeconvert.c
	gcc -Wall -O3 -o $@ $^ -DTYPE1=short -DTYPE2=float 
s2d: typeconvert.c
	gcc -Wall -O3 -o $@ $^ -DTYPE1=short -DTYPE2=double 

scale: typeconvert.c
	gcc -Wall -O3 -o $@ $^ -DTYPE1=float -DTYPE2=float 
scales: typeconvert.c
	gcc -Wall -O3 -o $@ $^ -DTYPE1=short -DTYPE2=short 
scaled: typeconvert.c
	gcc -Wall -O3 -o $@ $^ -DTYPE1=double -DTYPE2=double 

sigdiff: sigdiff.c
	gcc -Wall -O3 -o $@ $^ -DTYPE=float -lm
sigdiffs: sigdiff.c
	gcc -Wall -O3 -o $@ $^ -DTYPE=short -lm
sigdiffd: sigdiff.c
	gcc -Wall -O3 -o $@ $^ -DTYPE=double -lm

randshort.dat:
	dd if=/dev/urandom bs=$$((4 * $(NFFT) )) count=$(NUMFFTS) of=randshort.dat
randfloat.dat: randshort.dat s2f 
	./s2f < $< > $@
randdouble.dat: randshort.dat s2d 
	./s2d < $< > $@

data: randshort.dat randfloat.dat randdouble.dat

test: all data
	@echo "####################Testing double fft|ifft"
	@./kfftd < randdouble.dat | ./kifftd | ./scaled -d $(NFFT) | ./sigdiffd randdouble.dat -
	@echo "####################Testing float fft|ifft"
	./kfft < randfloat.dat | ./kifft | ./scale -d $(NFFT) | ./sigdiff randfloat.dat -
	@echo "####################Testing short fft|ifft"
	./kffts < randshort.dat | ./kiffts | ./scales -m $(NFFT) | ./sigdiffs randshort.dat -

clean:
	rm -f $(ALLUTILS) *.dat *~

