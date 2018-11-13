# palm: A package to fit point process models via the Palm likelihood

First proposed by Tanaka, Ogata, and Stoyan (2008), maximisation of the Palm likelihood can provide computationally efficient parameter estimation for point process models in situations where the full likelihood is intractable. This package contains functions to fit a variety of point process models, but is chiefly concerned with Neyman-Scott point processes (NSPPs).

The development of this package was motivated by the analysis of capture-recapture surveys on which individuals cannot be identified---the data from which can conceptually be seen as a NSPP (Fewster, Stevenson, and Borchers, 2016). As such, some of the functions in this package are specifically for the estimation of cetacean density from two-camera aerial surveys; see Stevenson, Borchers, and Fewster (in press).

This package can also fit void processes, which, along with NSPPs, have been fitted to patterns of colon cancer and stroma cell locations (Jones-Todd et al., in press).

## Installation

The stable version of this package can be installed from CRAN:

```
install.packages("palm")
```

## References

Fewster, R. M., Stevenson, B. C., and Borchers, D. L. (2016) Trace-contrast methods for capture-recapture without capture histories. *Statistical Science*, **31**: 245--258.

Jones-Todd, C. M., Caie, P., Illian, J. B., Stevenson, B. C., Savage, A., Harrison, D. J., and Bown, J. L. (in press). Identifying prognostic structural features in tissue sections of colon cancer patients using point pattern analysis. *Statistics in Medicine*.

Stevenson, B. C., Borchers, D. L., and Fewster, R. M. (in press) Cluster capture-recapture to account for identification uncertainty on aerial surveys of animal populations. *Biometrics*.

Tanaka, U., Ogata, Y., and Stoyan, D. (2008) Parameter estimation and model selection for Neyman-Scott point processes. *Biometrical Journal*, **50**: 43--57.