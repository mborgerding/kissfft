message:
	@echo "Nothing to make here.  Move on down to sample_code for ... you guessed it! Sample Code!"

tarball: clean
	tar --exclude CVS -czf kiss_fft.tar.gz .

clean:
	cd sample_code && make clean
	rm -f kiss_fft.tar.gz *~

