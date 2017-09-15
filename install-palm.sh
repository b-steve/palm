#!/bin/bash
cd ~/GitHub/palm
rm -rfv man
rm -fv NAMESPACE
rm -fv src/*.o src/RcppExports.cpp src/*.so R/RcppExports.R
rm -rfv package-build
R --slave -e "library(roxygen2); roxygenise('.')"
R --slave -e "library(Rcpp); compileAttributes()"
rm -rfv ..Rcheck/ ..pdf
rm -rfv src/*.o src/*.so src/*.rds
rm -rfv src-i386/ src-x64/
R CMD build .
mkdir -p package-build
mv palm_*.tar.gz package-build/
R CMD check package-build/palm_*.tar.gz
R CMD INSTALL --install-tests package-build/palm_*.tar.gz
