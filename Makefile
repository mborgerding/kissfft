message:
	@echo "Nothing to make here.  Move on down to sample_code for ... you guessed it! Sample Code!"

tarball: clean
	find | grep -i -v cvs | zip kiss_fft.zip -@
	tar --exclude CVS --exclude .cvsignore --exclude kiss_fft.zip -cvzf kiss_fft.tar.gz .

clean:
	cd sample_code && make clean
	rm -f kiss_fft.tar.gz *~ *.pyc kiss_fft.zip fft.py

