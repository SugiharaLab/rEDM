rtest:
	R -e "roxygen2::roxygenise()"
	R CMD build --no-build-vignettes .
	R CMD INSTALL rEDM_0.5.10.tar.gz
	R -e "roxygen2::roxygenise()"
	./test.r
