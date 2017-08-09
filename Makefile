rtest:
	R CMD build --no-build-vignettes .
	R CMD INSTALL rEDM_0.5.10.tar.gz
	./test.r
