message:
	@echo "Nothing to make here.  Move on down to test/ for self test stuff"
	@echo "or tools/ for real FFTs, multi-d FFTs, fast convolution filtering, cacher"

tarball: clean
	find | grep -i -v cvs | zip kiss_fft.zip -@
	tar --exclude CVS --exclude .cvsignore --exclude kiss_fft.zip -cvzf kiss_fft.tar.gz .

clean:
	cd test && make clean
	cd tools && make clean
	rm -f kiss_fft.tar.gz *~ *.pyc kiss_fft.zip 

