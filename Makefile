KFVER=111

DISTDIR=kiss_fft_v$(KFVER)
TARBALL=kiss_fft_v$(KFVER).tar.gz
ZIPFILE=kiss_fft_v$(KFVER).zip


testall:
	export DATATYPE=short && cd test && make test
	export DATATYPE=float && cd test && make test
	export DATATYPE=double && cd test && make test

tarball: clean
	tar --exclude CVS --exclude .cvsignore --exclude $(TARBALL) -cvzf $(TARBALL) .

clean:
	cd test && make clean
	cd tools && make clean
	rm -f kiss_fft*.tar.gz *~ *.pyc kiss_fft*.zip 
	rm -rf $(DISTDIR)

dist: tarball
	mkdir $(DISTDIR)
	cd $(DISTDIR) && tar -zxf ../$(TARBALL)
	rm $(TARBALL)
	tar -czf $(TARBALL) $(DISTDIR)
	zip -r $(ZIPFILE) $(DISTDIR)
	rm -rf $(DISTDIR)

upload: dist
	ncftpput upload.sourceforge.net incoming $(ZIPFILE) $(TARBALL)
