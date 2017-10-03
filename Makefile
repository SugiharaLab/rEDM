rtest:
	rm -rvf *.gz
	R CMD build --no-build-vignettes .
	R CMD INSTALL rEDM_0.6.0.tar.gz
	./test.r
