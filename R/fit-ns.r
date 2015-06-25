#' Fitting a model to a Neyman-Scott point process
#'
#' Estimates parameters for a Neyman-Scott point process using the
#' approach of Tanaka (2008).
#'
#' @references Tanaka, U., Ogata, Y., and Stoyan, D. (2008) Parameter
#' estimation and model selection for Neyman-Scott point
#' processes. \emph{Biometrical Journal}, \strong{50}: 43--57.
#'
#' @return An object with information that can be extracted via other
#' utility functions.
#'
#' @param points A matrix containing locations of observed points,
#' where each row corresponds to a point and each column corresponds
#' to a dimension.
#' @param lims A matrix with two columns, corresponding to the upper
#' and lower limits of each dimension, respectively.
#' @param R Truncation distance for the difference process.
#' @param sigma.sv The start value for \code{sigma} in optimisation.
#' @param sigma.bounds The bounds of the \code{sigma} parameter in
#' optimisation.
#' @param child.dist An argument describing the distribution
#' generating the number of children per parent. It must be a list
#' with four components: (1) A component named \code{mean} that
#' provides a function returning the expectation of the distribution,
#' as an argument of a single parameter to be estimated, (2) A
#' component named \code{var} that provides a function returning the
#' variance of a distribution, as an argument of a singla parameter to
#' be estimated, (3) A component named \code{sv} providing the start
#' value for the parameter to be estimated, and (4) A component named
#' \code{bounds} providing a vector of length two that gives the
#' parameter bounds.
#' @param siblings A named list, containing the following three
#' components: (1) A component named \code{matrix}, where the jth
#' column of the ith row is \code{TRUE} if the ith and jth observed
#' points are known siblings, \code{FALSE} if they are known
#' non-siblings, and \code{NA} if it is not known whether or not they
#' are siblings, (2) A component named pT, containing a scalar
#' providing the probability that the element (i, j) in the component
#' \code{matrix} is \code{TRUE}, conditional on the ith and jth points
#' being siblings, (3) A component named pF, containing a scalar
#' providing the probability that the element (i, j) in the component
#' \code{matrix} is \code{FALSE}, conditional on the ith and jth
#' points not being siblings.
#' @param trace Logical, if \code{TRUE}, parameter values are printed
#' to the screen for each iteration of the optimisation procedure.
#'
#' @export
fit.ns <- function(points = NULL, lims = NULL, R, sigma.sv = 0.1*R,
                   sigma.bounds = c(0, R),
                   child.dist = list(mean = function(x) x, var = function(x) x,
                       sv = 5, bounds = c(1e-8, 1e8)),
                   siblings = NULL, trace = FALSE){
    ## Saving arguments.
    arg.names <- names(as.list(environment()))
    args <- vector(mode = "list", length = length(arg.names))
    names(args) <- arg.names
    for (i in arg.names){
        if (!is.null(get(i))){
            args[[i]] <- get(i)
        }
    }
    ## Errors for inconsistent dimensions.
    if (!is.matrix(points)){
        points <- matrix(points, nrow = 1)
    }
    if (!is.matrix(lims)){
        lims <- matrix(lims, nrow = 1)
    }
    ## Error for incompatible dimensions.
    error.dims(points, lims)
    n.points <- nrow(points)
    n.dims <- nrow(lims)
    ## Vectorising siblings matrix, and setting intensity function.
    if (!is.null(siblings)){
        siblings <- vectorise.siblings(siblings)
        intensity.fun <- palm.intensity.siblings
    } else {
        intensity.fun <- palm.intensity
    }
    ## Declaring function to calculate nu.
    nu.fun <- function(x, child.dist){
        ## Mean and variance for number of children.
        var.c <- child.dist$var(x)
        mean.c <- child.dist$mean(x)
        (var.c + mean.c^2)/mean.c - 1
    }
    ## Calculating survey area.
    area <- prod(apply(lims, 1, diff))
    ## Calculating distances.
    dists <- pbc_distances(points = points, lims = lims)
    ## Truncating distances.
    dists <- dists[dists <= R]
    n.dists <- length(dists)
    ## Sorting out start values.
    nu.sv <- nu.fun(child.dist$sv, child.dist)
    Dc.sv <- analytic.Dc(nu.sv, sigma.sv, n.dists, n.points, R, n.dims)
    sv <- c(Dc.sv, nu.sv, sigma.sv)
    names(sv) <- c("Dc", "nu", "sigma")
    ## Sorting out bounds.
    Dc.bounds <- c(0, Inf)
    nu.bounds <- nu.fun(child.dist$bounds, child.dist)
    if (is.nan(nu.bounds[1])) nu.bounds[1] <- 0
    lower <- c(Dc.bounds[1], nu.bounds[1], sigma.bounds[1])
    upper <- c(Dc.bounds[2], nu.bounds[2], sigma.bounds[2])
    fit <-  optimx(par = log(sv), fn = ns.nll,
                   method = "L-BFGS-B",
                   lower = log(lower),
                   upper = log(upper),
                   n.points = n.points, dists = dists, R = R,
                   d = n.dims, nu.fun = nu.fun,
                   par.names = names(sv), siblings = siblings,
                   intensity.fun = intensity.fun, trace = trace)
    ## Extracting sigma and nu estimates.
    opt.pars <- exp(coef(fit)[1, ])
    names(opt.pars) <- names(sv)
    sigma.par <- opt.pars["sigma"]
    nu.par <- opt.pars["nu"]
    ## Calculating Dc estimate.
    ##Dc.par <- analytic.Dc(nu.par, sigma.par, n.dists, n.points, R)
    Dc.par <- opt.pars["Dc"]
    ## Search bounds slightly outwith optimised bounds.
    search.bounds <- child.dist$bounds
    search.bounds[1] <- search.bounds[1] - 0.04*diff(child.dist$bounds)
    search.bounds[2] <- search.bounds[2] + 0.04*diff(child.dist$bounds)
    ## Calculating parameter of child distribution.
    child.par <- try(uniroot(function(x, child.dist, nu) nu.fun(x, child.dist) - nu,
                         interval = search.bounds, child.dist = child.dist,
                             nu = nu.par)$root, silent = FALSE)
    ## Calculating mu.
    mu.par <- child.dist$mean(child.par)
    D.par <- Dc.par/mu.par
    pars <- c(D.par, sigma.par, child.par, Dc.par, mu.par, nu.par)
    names(pars) <- c("D", "sigma", "child.par", "Dc", "mu", "nu")
    out <- list(pars = pars, args = args)
    class(out) <- "nspp"
    out
}

ns.nll <- function(pars, n.points, dists, R, d, nu.fun, par.names, siblings,
                   intensity.fun, trace){
    ## Extracting parameters.
    names(pars) <- par.names
    pars <- exp(pars)
    Dc <- pars["Dc"]
    nu <- pars["nu"]
    sigma <- pars["sigma"]
    ## Can work out Dc analytically.
    ll1 <- sum(log(n.points*intensity.fun(dists, Dc, nu, sigma, d, siblings)))
    ## Contribution from integral.
    ll2 <- n.points*(Dc*Vd(R, d) + nu*Fd(R, sigma, d))
    ll <- ll1 - ll2
    ## Printing parameter values.
    if (trace){
        cat("Dc = ", Dc, ", sigma = ", sigma, ", nu = ", nu, ", LL = ", ll, "\n", sep = "")
    }
    -ll
}

## Roxygen code for NAMESPACE.
#' @import Rcpp
#' @importFrom mvtnorm rmvnorm
#' @importFrom optimx optimx
#' @useDynLib nspp
NULL

## Data documentation.

#' 1-dimensional example data
#'
#' Simulated data, with children points generated from a Binomial(4,
#' 0.5) distribution.
#'
#' @name example.1D
#' @format A matrix.
#' @usage example.1D
#' @docType data
#' @keywords datasets
NULL

#' 2-dimensional example data
#'
#' Simulated data, with children points generated from a Binomial(2,
#' 0.5) distribution.
#'
#' @name example.2D
#' @format A matrix.
#' @usage example.2D
#' @docType data
#' @keywords datasets
NULL
