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
#' @param non.siblings An object containing information about known
#' siblings (or known non-siblings). Not yet implemented.
#' @param sv A vector (or list?) of start values for optimisation.
#'
#' @export
fit.ns <- function(points = NULL, lims = NULL, R, sigma.sv = 0.1*R,
                   sigma.bounds = c(0, R),
                   child.dist = list(mean = function(x) x, var = function(x) x,
                       sv = 5, bounds = c(1e-8, 1e8)),
                   non.siblings = NULL){
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
    if (ncol(points) != nrow(lims)){
        stop("The number of columns in 'points' and 'lims' must both equal the number of dimensions.")
    }
    n.points <- nrow(points)
    n.dims <- nrow(lims)
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
    Dc.sv <- analytic.Dc(nu.sv, sigma.sv, n.dists, n.points, R)
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
                   n.dims = n.dims, nu.fun = nu.fun,
                   par.names = names(sv))
    ## Extracting sigma and nu estimates.
    opt.pars <- exp(coef(fit)[1, ])
    names(opt.pars) <- names(sv)
    sigma.par <- opt.pars["sigma"]
    nu.par <- opt.pars["nu"]
    ## Calculating Dc estimate.
    ##Dc.par <- analytic.Dc(nu.par, sigma.par, n.dists, n.points, R)
    Dc.par <- opt.pars["Dc"]
    ## Calculating parameter of child distribution.
    child.par <- uniroot(function(x, child.dist, nu) nu.fun(x, child.dist) - nu,
                         interval = child.dist$bounds, child.dist = child.dist,
                         nu = nu.par)$root
    ## Calculating mu.
    mu.par <- child.dist$mean(child.par)
    D.par <- Dc.par/mu.par
    pars <- c(D.par, sigma.par, child.par, Dc.par, mu.par, nu.par)
    names(pars) <- c("D", "sigma", "child.par", "Dc", "mu", "nu")
    out <- list(pars = pars, args = args)
    class(out) <- "nspp"
    out
}

ns.nll <- function(pars, n.points, dists, R, n.dims, nu.fun, par.names){
    ## Extracting parameters.
    names(pars) <- par.names
    pars <- exp(pars)
    Dc <- pars["Dc"]
    nu <- pars["nu"]
    sigma <- pars["sigma"]
    ## Can work out Dc analytically.
    ##Dc <- analytic.Dc(nu, sigma, length(dists), n.points, R)
    vol <- n.points*(pi*Dc*R^2 + nu - nu*exp((-R^2)/(4*sigma^2)))
    ##if (nu > (length(dists)/n.points + nu*exp((-R^2)/(4*sigma^2)))){
    ##    ll <- NA
    ##} else {
        ll <- sum(log(n.points*ns.palm(dists, Dc, nu, sigma, n.dims))) - vol
    ##}
    ## Printing parameter values.
    cat("Dc = ", Dc, ", sigma = ", sigma, ", nu = ", nu, ", LL = ", ll, "\n", sep = "")
    -ll
}

## To generalise for n dimensions.
ns.palm <- function(r, Dc, nu, sigma, n.dims){
    Dc + nu*exp((-r^2)/(4*sigma^2))/(4*pi*sigma^2)
}

## Roxygen code for NAMESPACE.
#' @import Rcpp
#' @importFrom mvtnorm rmvnorm
#' @importFrom optimx optimx
#' @useDynLib nspp
NULL

