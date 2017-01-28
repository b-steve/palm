# palm

This package provides functions for the fitting of point process models using the Palm likelihood. Maximisation of the Palm likelihood can provide computationally efficient parameter estimation in situations where the full likelihood is intractable. This package is chiefly focussed on Neyman-Scott point processes, but can also fit void processes. Estimation via the Palm likelihood was first proposed by Tanaka, Ogata, and Stoyan (2008; Biometrical Journal) and further generalised by both Stevenson, Borchers, and Fewster (in prep) and Jones-Todd (2017; PhD thesis, University of St Andrews).

The development of this package was motivated by the analysis of capture-recapture surveys on which individuals cannot be identified---the data from which can conceptually be seen as a clustered point process. Some of the functions in this package are specifically for the estimation of cetacean density from two-camera aerial surveys.

## Installation

The stable version of this package can be installed from CRAN:

```
install.packages("palm")
```
