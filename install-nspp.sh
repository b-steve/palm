#!/bin/bash
rm -rfv man
rm -fv NAMESPACE
#R --slave -e "library(Rcpp); compileAttributes()"
R --slave -e "library(roxygen2); roxygenise('.')"
R CMD build .
mkdir -p package-build
mv nspp_*.tar.gz package-build/
R CMD check package-build/nspp_*.tar.gz
R CMD INSTALL package-build/nspp_*.tar.gz
