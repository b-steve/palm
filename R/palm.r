## Main package-level documentation.

#' palm: A package to fit point process models via the Palm
#' likelihood
#'
#' First proposed by Tanaka, Ogata, and Stoyan (2008), maximisation of
#' the Palm likelihood can provide computationally efficient parameter
#' estimation for point process models in situations where the full
#' likelihood is intractable. This package contains functions to fit a
#' variety of point process models, but is chiefly concerned with
#' Neyman-Scott point processes (NSPPs).
#'
#' The development of this package was motivated by the analysis of
#' capture-recapture surveys on which individuals cannot be
#' identified---the data from which can conceptually be seen as a NSPP
#' (Fewster, Stevenson, and Borchers, 2016). As such, some of the
#' functions in this package are specifically for the estimation of
#' cetacean density from two-camera aerial surveys; see Stevenson,
#' Borchers, and Fewster (2019).
#'
#' This package can also fit void processes, which, along with NSPPs,
#' have been fitted to patterns of colon cancer and stroma cell
#' locations (Jones-Todd et al., 2019).
#'
#' The main functions of this package are summarised below.
#' 
#' @section Model fitting:
#' 
#' \itemize{
#'
#' \item The \link{fit.ns} function fits NSPPs.
#'
#' \item The \link{fit.twocamera} function estimates animal density
#'       from two-camera aerial surveys. This model is a NSPP and can
#'       be fitted using \link{fit.ns}, but it is more straightforward
#'       to use \link{fit.twocamera}.
#' 
#' \item The \link{fit.void} function fits void point processes.
#' 
#' }
#'
#' @section Variance estimation:
#' 
#' Variance estimation is achieved by parametric bootstrap. The
#' \link{boot.palm} function carries out this procedure from an object
#' generated by one of the fitting functions, above. Confidence
#' intervals and standard errors can be calculated from an object
#' returned by \link{boot.palm} using \link{confint.palm} and
#' \link{coef.palm}, respectively.
#'
#' @section Data simulation:
#' 
#' \itemize{
#'
#' \item The \link{sim.ns} function simulates data from NSPPs.
#'
#' \item The \link{sim.twocamera} function simulates detection data
#'       from two-camera aerial surveys.
#'
#' \item The \link{sim.void} function simulates data from void point
#'       processes.
#'
#' }
#'
#' @references Fewster, R. M., Stevenson, B. C., and Borchers,
#'     D. L. (2016) Trace-contrast methods for capture-recapture
#'     without capture histories. \emph{Statistical Science},
#'     \strong{31}: 245--258.
#' @references Jones-Todd, C. M., Caie, P., Illian, J. B., Stevenson,
#'     B. C., Savage, A., Harrison, D. J., and Bown, J. L. (in
#'     press). Identifying prognostic structural features in tissue
#'     sections of colon cancer patients using point pattern
#'     analysis. \emph{Statistics in Medicine}, \strong{38}:
#'     1421--1441.
#' @references Stevenson, B. C., Borchers, D. L., and Fewster,
#'     R. M. (2019) Cluster capture-recapture to account for
#'     identification uncertainty on aerial surveys of animal
#'     populations. \emph{Biometrics}, \strong{75}: 326--336.
#' @references Tanaka, U., Ogata, Y., and Stoyan, D. (2008) Parameter
#'     estimation and model selection for Neyman-Scott point
#'     processes. \emph{Biometrical Journal}, \strong{50}: 43--57.
#'
#' @docType package
#' @name palm
NULL

## Roxygen code for NAMESPACE.

#' @import methods Rcpp R6
#' @importFrom graphics lines par plot title
#' @importFrom gsl hyperg_2F1
#' @importFrom minqa bobyqa
#' @importFrom mvtnorm rmvnorm
#' @importFrom Rcpp sourceCpp
#' @importFrom spatstat crossdist
#' @importFrom stats coef dist integrate nlminb pbeta pgamma pnorm printCoefmat qnorm quantile rbinom rnorm rpois runif sd
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @useDynLib palm, .registration = TRUE
NULL

## Data documentation.

#' 1-dimensional example data
#' 
#' Simulated data from a Neyman-Scott point process, with children
#' points generated in the interval [0, 1]. The number of children
#' spawned by each parent is from a Binomial(4, 0.5) distribution.
#' 
#' @name example.1D
#' @format A matrix.
#' @usage example.1D
#' @docType data
#' @keywords datasets
NULL

#' 2-dimensional example data
#'
#' Simulated data from a Neyman-Scott point process, with children
#' points generated on the unit square. The number of children spawned
#' by each parent is from a Binomial(2, 0.5) distribution.
#'
#' @name example.2D
#' @format A matrix.
#' @usage example.2D
#' @docType data
#' @keywords datasets
NULL

#' Two-camera example data.
#'
#' Simulated data from a two-camera aerial survey.
#'
#' @name example.twocamera
#' @format A list.
#' @usage example.twocamera
#' @docType data
#' @keywords datasets
NULL

#' Two-camera porpoise data.
#'
#' Synthetic data constructed from circle-back aerial survey data; see
#' Stevenson, Borchers, and Fewster (2019).
#'
#' @name porpoise.data
#' @format A list.
#' @usage porpoise.data
#' @docType data
#' @keywords datasets
#' @references Stevenson, B. C., Borchers, D. L., and Fewster,
#'     R. M. (2019) Cluster capture-recapture to account for
#'     identification uncertainty on aerial surveys of animal
#'     populations. \emph{Biometrics}, \strong{75}: 326--336.
NULL
