


testall:
	export DATATYPE=short && cd test && make test
	export DATATYPE=float && cd test && make test
	export DATATYPE=double && cd test && make test

tarball: clean
	find | grep -i -v cvs | zip kiss_fft.zip -@
	tar --exclude CVS --exclude .cvsignore --exclude kiss_fft.zip -cvzf kiss_fft.tar.gz .

clean:
	cd test && make clean
	cd tools && make clean
	rm -f kiss_fft.tar.gz *~ *.pyc kiss_fft.zip 

