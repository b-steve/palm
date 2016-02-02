#!/bin/bash
cd ~/GitHub/nspp
rm -rfv man
rm -fv NAMESPACE
R --slave -e "library(Rcpp); compileAttributes()"
R --slave -e "library(roxygen2); roxygenise('.')"
rm -rfv ..Rcheck/ ..pdf
rm -rfv src/*.o src/*.so src/*.rds
rm -rfv src-i386/ src-x64/
R CMD build .
mkdir -p package-build
mv nspp_*.tar.gz package-build/
R CMD check package-build/nspp_*.tar.gz
R CMD INSTALL --install-tests package-build/nspp_*.tar.gz
