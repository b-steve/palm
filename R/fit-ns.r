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
#'     where each row corresponds to a point and each column
#'     corresponds to a dimension.
#' @param lims A matrix with two columns, corresponding to the upper
#'     and lower limits of each dimension, respectively.
#' @param R Truncation distance for the difference process.
#' @param disp A character string, providing the distribution for the
#'     dispersion of points around the parents. This can either be
#'     \code{"gaussian"} for Gaussian dispersion, or \code{"matern"}
#'     for uniform dispersion (giving a matern process).
#' @param sigma.sv The start value for \code{sigma} in optimisation.
#' @param sigma.bounds The bounds of the \code{sigma} parameter in
#'     optimisation.
#' @param child.dist An argument describing the distribution
#'     generating the number of children per parent. It must be a list
#'     with four components: (1) A component named \code{mean} that
#'     provides a function returning the expectation of the
#'     distribution, with its only argument being the single parameter
#'     to be estimated, (2) A component named \code{var} that provides
#'     a function returning the variance of a distribution, with its
#'     only argument being the single parameter to be estimated, (3) A
#'     component named \code{sv} providing the start value for the
#'     parameter to be estimated, and (4) A component named
#'     \code{bounds} providing a vector of length two that gives the
#'     parameter bounds.
#' @param edge.correction The method used for the correction of edge
#'     effects. Either \code{"pbc"} for periodic boundary conditions,
#'     or \code{"buffer"} for a buffer-zone correction.
#' @param siblings A named list, containing the following three
#'     components: (1) A component named \code{matrix}, where the jth
#'     element of the ith row is \code{TRUE} if the ith and jth
#'     observed points are known siblings, \code{FALSE} if they are
#'     known non-siblings, and \code{NA} if it is not known whether or
#'     not they are siblings, (2) A component named pT, containing a
#'     scalar providing the probability that the element (i, j) in the
#'     component \code{matrix} is \code{TRUE}, conditional on the ith
#'     and jth points being siblings, (3) A component named pF,
#'     containing a scalar providing the probability that the element
#'     (i, j) in the component \code{matrix} is \code{FALSE},
#'     conditional on the ith and jth points not being
#'     siblings. Defaults to an argument representing a Poisson
#'     distribution.
#' @param trace Logical, if \code{TRUE}, parameter values are printed
#'     to the screen for each iteration of the optimisation procedure.
#'
#' @examples
#'## Poisson number of children per parent.
#'fit.2D.pois <- fit.ns(example.2D, lims = rbind(c(0, 1), c(0, 1)), R = 0.5)
#'## Binomial number of children per parent.
#'fit.1D.binom <- fit.ns(example.1D, lims = cbind(0, 1), R = 0.5,
#'                       child.dist = list(mean = function(x) 4*x,
#'                           var = function(x) 4*x*(1 - x),
#'                           sv = 0.5, bounds = c(0, 1)))
#' 
#' @export
fit.ns <- function(points, lims = NULL, R, disp = "gaussian",
                   sigma.sv = 0.1*R, sigma.bounds = c(1e-10, R),
                   child.dist = list(mean = function(x) x, var = function(x) x,
                                     sv = 5, bounds = c(1e-8, 1e8)),
                   edge.correction = "pbc", siblings = NULL, trace = FALSE){
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
    error.dims(points, lims)
    n.points <- nrow(points)
    n.dims <- nrow(lims)
    ## Error for sibling matrix not matching.
    if (!is.null(siblings)){
        if (nrow(siblings$matrix) != nrow(points)){
            stop("Sibling matrix does not have a row for every detection.")
        }
    }
    ## Declaring function to calculate nu.
    protect.child.var <- function(x, child.dist, sigma){
        if (any(names(formals(child.dist$var)) == "sigma")){
            if (is.null(sigma)){
                stop("Value of sigma not specified.")
            }
            out <- child.dist$var(x, sigma)
        } else {
            out <- child.dist$var(x)
        }
        out
    }
    protect.child.mean <- function(x, child.dist, sigma){
        if (any(names(formals(child.dist$mean)) == "sigma")){
            if (is.null(sigma)){
                stop("Value of sigma not specified.")
            }
            out <- child.dist$mean(x, sigma)
        } else {
            out <- child.dist$mean(x)
        }
        out
    }
    nu.fun <- function(x, child.dist, sigma = NULL){
        ## Mean and variance for number of children.
        var.c <- protect.child.var(x, child.dist, sigma)
        mean.c <- protect.child.mean(x, child.dist, sigma)
        (var.c + mean.c^2)/mean.c - 1
    }
    ## Calculating survey area.
    area <- prod(apply(lims, 1, diff))
    ## Calculating distances.
    if (edge.correction == "pbc"){
        buffer.keep <- NULL
        dists <- pbc_distances(points = points, lims = lims)
    } else if (edge.correction == "buffer"){
        buffer.keep <- buffer_keep(points = points, lims = lims, R = R)
        dists <- as.vector(as.matrix(dist(points))[buffer.keep])
    } else {
        stop("Edge correction method not recognised.")
    }
    ## Sorting out sibling stuff.
    if (is.null(siblings)){
        v.siblings <- NULL
        intensity.fun <- palm.intensity
    } else {
        v.siblings <- vectorise.siblings(siblings, edge.correction, buffer.keep)
        intensity.fun <- palm.intensity.siblings
    }
    ## Truncating distances.
    keep <- dists <= R
    dists <- dists[keep]
    v.siblings$ns.multipliers <- v.siblings$ns.multipliers[keep]
    v.siblings$s.multipliers <- v.siblings$s.multipliers[keep]
    n.dists <- length(dists)
    ## Sorting out start values.
    nu.sv <- nu.fun(child.dist$sv, child.dist, sigma.sv)
    Dc.sv <- analytic.Dc(nu.sv, sigma.sv, n.dists, n.points, R, n.dims, disp)
    if (Dc.sv <= 0){
        Dc.sv <- n.points/area
    }
    sv <- c(Dc.sv, nu.sv, sigma.sv)
    names(sv) <- c("Dc", "nu", "sigma")
    ## Sorting out bounds.
    Dc.bounds <- c(0, Inf)
    nu.optim.fun <- function(par, child.dist, max){
        if (max){
            out <- -nu.fun(par[1], child.dist, par[2])
        } else {
            out <- nu.fun(par[1], child.dist, par[2])
        }
        if (is.nan(out)){
            out <- NA
        }
        out
    }
    ## Sorting out bounds for nu (complicated as it is affected by sigma and child.par).
    if (any(names(formals(child.dist$mean)) == "sigma") | any(names(formals(child.dist$var)) == "sigma")){
        nu.max <- -optim(c(child.dist$sv, sigma.sv), nu.optim.fun, child.dist = child.dist, max = TRUE,
                         method = "L-BFGS-B", lower = c(child.dist$bounds[1], sigma.bounds[1]),
                         upper = c(child.dist$bounds[2], sigma.bounds[2]))$value
        nu.min <- optim(c(child.dist$sv, sigma.sv), nu.optim.fun, child.dist = child.dist, max = FALSE,
                        method = "L-BFGS-B", lower = c(child.dist$bounds[1], sigma.bounds[1]),
                        upper = c(child.dist$bounds[2], sigma.bounds[2]))$value
        nu.bounds <- c(nu.min, nu.max)
    } else {
        nu.bounds <- nu.fun(child.dist$bounds, child.dist)
        if (is.nan(nu.bounds[1])) nu.bounds[1] <- 0
        nu.bounds <- sort(nu.bounds)
    }
    lower <- c(Dc.bounds[1], nu.bounds[1], sigma.bounds[1])
    upper <- c(Dc.bounds[2], nu.bounds[2], sigma.bounds[2])
    fit <-  optimx(par = log(sv), fn = ns.nll,
                   method = "L-BFGS-B",
                   lower = log(lower),
                   upper = log(upper),
                   n.points = n.points, dists = dists, R = R,
                   d = n.dims,
                   par.names = names(sv), siblings = v.siblings,
                   intensity.fun = intensity.fun, disp = disp,
                   trace = trace)
    ## Estimation from system of partial derivatives.
    ## fit.system <- nleqslv(c(sv["Dc"]/sv["nu"], sv["nu"], sv["sigma"),
    ##                       function(x, n.points, dists, R) c(dldD(x[1], x[2], x[3], n.points, dists, R),
    ##                                                         dldnu(x[1], x[2], x[3], n.points, dists, R),
    ##                                                         dldsigma(x[1], x[2], x[3], n.points, dists, R)),
    ##                       n.points = n.points, dists = dists, R = R)
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
    #search.bounds[1] <- search.bounds[1] - 0.04*diff(child.dist$bounds)
    #search.bounds[2] <- search.bounds[2] + 0.04*diff(child.dist$bounds)
    ## Calculating parameter of child distribution.
    ##child.par <- try(uniroot(function(x, child.dist, nu, sigma)
    ##    nu.fun(x, child.dist, sigma) - nu,
    ##                         interval = search.bounds, child.dist = child.dist,
    ##                         nu = nu.par, sigma = sigma.par)$root, silent = FALSE)
    child.par.optim <- function(child.par, child.dist, nu, sigma){
        (nu.fun(child.par, child.dist, sigma) - nu)^2
    }
    ## Kludgey fix for parameter hitting 0 in secondary optimisation.
    if (child.dist$bounds[1] == 0){
        child.dist$bounds[1] <- .Machine$double.xmin
    }
    child.par.optim <- optim(child.dist$sv, child.par.optim, child.dist = child.dist,
                       sigma = sigma.par, nu = nu.par, method = "L-BFGS-B",
                             lower = child.dist$bounds[1], upper = child.dist$bounds[2])
    child.par <- child.par.optim$par
    ## Calculating mu.
    mu.par <- protect.child.mean(child.par, child.dist, sigma.par)
    D.par <- Dc.par/mu.par
    pars <- c(D.par, sigma.par, child.par, Dc.par, mu.par, nu.par)
    names(pars) <- c("D", "sigma", "child.par", "Dc", "mu", "nu")
    out <- list(pars = pars, args = args)
    class(out) <- "nspp"
    out
}

fit.ns_refclass <- function(points, lims, R, disp = "gaussian",
                            child.dist = list(child.expectation = function(x) x["child.par"],
                                              child.variance = function(x) x["child.par"],
                                              child.link = log),
                            edge.correction = "pbc", siblings = NULL, trace = FALSE){
    use.thomas.class <- FALSE
    use.matern.class <- FALSE
    if (disp == "gaussian"){
        use.thomas.class <- TRUE
    } else if (disp == "uniform"){
        use.matern.class <- TRUE
    } else {
        stop("Dispersion type not recognised; use either 'gaussian' or 'uniform'.")
    }
    use.pbc.class <- FALSE
    use.buffer.class <- FALSE
    if (edge.correction == "pbc"){
        use.pbc.class <- TRUE
    } else if (edge.correction == "buffer"){
        use.buffer.class <- TRUE
    } else {
        stop("Dispersion type not recognised; use either 'pbc' or 'buffer'.")
    }
    use.sibling.class <- FALSE
    if (!is.null(siblings)){
        use.sibling.class <- TRUE
    }
    final.contains <- c("thomas"[use.thomas.class],
                        "matern"[use.matern.class],
                        "pbc"[use.pbc.class],
                        "buffer"[use.buffer.class],
                        "sibling"[use.sibling.class])
    final.class <- setRefClass("final", contains = final.contains)
    obj <- final.class$new(points = points, lims = lims, R = 0.5, child.dist = child.dist)
    obj$get.invlinks()
    obj$initialise.par.start.link()
    obj$fit()
    obj$par.fitted
}

ns.nll <- function(pars, n.points, dists, R, d, par.names, siblings,
                   intensity.fun, disp, trace){
    ## Gettind CDF of between-sibling distances.
    Fd <- get.Fd(disp)
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
        if (d == 2){
            cat("Partial derivative for D: ", dldD(Dc/nu, nu, sigma, n.points, dists, R), "\n",
                sep = "")
            cat("Partial derivative for nu: ", dldnu(Dc/nu, nu, sigma, n.points, dists, R), "\n",
                sep = "")
            cat("Partial derivative for sigma: ", dldsigma(Dc/nu, nu, sigma, n.points, dists, R),
                "\n", sep = "")
        }
    }
    -ll
}

fit.void <- function(points, lims = NULL, R, edge.correction = "pbc",
                     trace = FALSE){
    stop("Fitting of void processes not yet implemented.")
}

#' Estimation of animal density from two-plane surveys.
#'
#' Estimates animal density (amongst other parameters) from two-plane
#' aerial surveys. This conceptualises sighting locations as a
#' Neyman-Scott point pattern---estimation is carried out via
#' \code{fit.ns()}.
#'
#' @return An object of class \code{"nspp"} that can be extracted via
#' the same utility functions fit for objects created using
#' \code{fit.ns()}.
#'
#' @param points A vector (or single-column matrix) containing the
#'     distance along the transect that each detection was made.
#' @param planes A vector containing the plane ID (either \code{1} or
#'     \code{2}) that made the corresponding detection in
#'     \code{points}.
#' @param l The length of the transect flown (in km).
#' @param w The distance from the transect to which detection of
#'     individuals on the surface is certain. This is equivalent to
#'     the half-width of the detection zone.
#' @param b The distance from the transect to the edge of the area of
#'     interest. Conceptually, the distance between the transect and
#'     the furthest distance a whale could be on the passing on the
#'     first plane and plausibly move into the detection zone by the
#'     passing of the second plane.
#' @param t The lag between planes (in seconds).
#' @param C Mean dive-cycle duration (in seconds).
#' @param R Truncation distance (see \link{fit.ns}).
#' @param edge.correction The method used for the correction of edge
#'     effects. Either \code{"pbc"} for periodic boundary conditions,
#'     or \code{"buffer"} for a buffer-zone correction.
#' @param trace Logical, if \code{TRUE}, parameter values are printed
#'     to the screen for each iteration of the optimisation procedure.
#' 
#' @export
fit.twoplane <- function(points, planes = NULL, l, w, b, t, C, R,
                         edge.correction = "pbc", trace = FALSE){
    if (!is.matrix(points)){
        points <- matrix(points, nrow = length(points), ncol = 1)
    }
    ## Saving arguments.
    arg.names <- names(as.list(environment()))
    args <- vector(mode = "list", length = length(arg.names))
    names(args) <- arg.names
    for (i in arg.names){
        if (!is.null(get(i))){
            args[[i]] <- get(i)
        }
    }
    if (is.null(planes)){
        sibling.mat = NULL
    } else {
        sibling.mat <- siblings.twoplane(planes)
    }
    child.dist <- make.twoplane.child.dist(t, C, w, b)
    fit <- fit.ns(points = points, lims = rbind(c(0, l)), R, sigma.sv = b/5,
                  child.dist = child.dist, siblings = sibling.mat,
                  edge.correction = edge.correction, trace = trace)
    fit$pars <- c(fit$pars, fit$pars[1]/(2*b))
    names(fit$pars)[length(fit$pars)] <- "D.2D"
    fit$args.twoplane <- args
    class(fit) <- c("twoplane.nspp", class(fit))
    fit
}


## Roxygen code for NAMESPACE.
#' @import methods Rcpp
#' @importFrom graphics abline axis box lines par plot.new plot.window points
#' @importFrom mvtnorm rmvnorm
#' @importFrom optimx optimx
#' @importFrom stats coef dist integrate optim pgamma pnorm printCoefmat qnorm quantile rnorm rpois runif sd var
#' @importFrom utils setTxtProgressBar txtProgressBar
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
