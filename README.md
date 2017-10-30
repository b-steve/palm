# palm

[![Build Status](https://travis-ci.org/b-steve/palm.svg?branch=master)](https://travis-ci.org/b-steve/palm)
[![downloads](http://cranlogs.r-pkg.org/badges/R2palm)](http://cranlogs.r-pkg.org/badges/palm)

This package provides functions for the fitting of point process models using the Palm likelihood. Maximisation of the Palm likelihood can provide computationally efficient parameter estimation in situations where the full likelihood is intractable. This package is chiefly focussed on Neyman-Scott point processes, but can also fit void processes. Estimation via the Palm likelihood was first proposed by Tanaka, Ogata, and Stoyan (2008; Biometrical Journal) and further generalised by both Stevenson, Borchers, and Fewster (in review) and Jones-Todd (in prep).

The development of this package was motivated by the analysis of capture-recapture surveys on which individuals cannot be identified---the data from which can conceptually be seen as a clustered point process. Some of the functions in this package are specifically for the estimation of cetacean density from two-camera aerial surveys.

## Installation

The stable version of this package can be installed from CRAN:

```
install.packages("palm")
```
