#' Fitting a Neyman-Scott point process model
#'
#' Estimates parameters for a Neyman-Scott point process by maximising
#' the Palm likelihood. This approach was first proposed by Tanaka et
#' al. (2008) for two-dimensional Thomas processes. Further
#' generalisations were made by Stevenson, Borchers, and Fewster (in
#' press) and Jones-Todd et al. (in press).
#' 
#' The parameter \code{D} is the density of parent points, which is
#' always estimated. Possible additional parameters are
#' \itemize{
#'   \item \code{lambda}, the expected number of children generated per
#'         parent (when \code{child.dist = "pois"}).
#' 
#'   \item \code{p}, the proportion of the \code{x} possible children
#'         that are generated (when \code{child.dist = "binomx"}).
#'
#'   \item \code{kappa}, the average length of the surface phase of a
#'         diving cetacean (when \code{child.dist = "twocamera"}; see
#'         Stevenson, Borchers, and Fewster, in press).
#'
#'   \item \code{sigma}, the standard deviation of dispersion along
#'         each dimension (when \code{disp} = "gaussian").
#'
#'   \item \code{tau}, the maximum distance a child can be from its
#'         parent (when \code{disp} = "uniform").
#'
#' }
#'
#' The \code{"child.info"} argument is required when \code{child.dist}
#' is set to \code{"twocamera"}. It must be a list that comprises (i) a
#' component named \code{w}, providing the halfwidth of the detection
#' zone; (ii) a component named \code{b}, providing the halfwidth of
#' the survey area; (iii) a component named \code{l}, providing the
#' time lag between cameras (in seconds); and (iv) a component named
#' \code{tau}, providing the mean dive-cycle duration. See Stevenson,
#' Borchers, and Fewster (in press) for details.
#'
#' @references Jones-Todd, C. M., Caie, P., Illian, J. B., Stevenson,
#'     B. C., Savage, A., Harrison, D. J., and Bown, J. L. (in
#'     press). Identifying prognostic structural features in tissue
#'     sections of colon cancer patients using point pattern
#'     analysis. \emph{Statistics in Medicine}.
#' @references Stevenson, B. C., Borchers, D. L., and Fewster,
#'     R. M. (in press) Cluster capture-recapture to account for
#'     identification uncertainty on aerial surveys of animal
#'     populations. \emph{Biometrics}.
#' @references Tanaka, U., Ogata, Y., and Stoyan, D. (2008) Parameter
#'     estimation and model selection for Neyman-Scott point
#'     processes. \emph{Biometrical Journal}, \strong{50}: 43--57.
#'
#' @param points A matrix or list of matrices containing locations of
#'     observed points, where each row corresponds to a point and each
#'     column corresponds to a dimension. If a list, then the patterns
#'     are assumed to be independent and a single process is fitted to
#'     all.
#' @param lims A matrix or list of matrices with two columns,
#'     corresponding to the upper and lower limits of each dimension,
#'     respectively. If a list, then each matrix provides the limits
#'     for the corresponding pattern in \code{points}.
#' @param disp A character string indicating the distribution of
#'     children around their parents. Use \code{"gaussian"} for
#'     multivariate normal dispersion with standard deviation
#'     \code{sigma}, or \code{"uniform"} for uniform dispersion within
#'     distance \code{tau} of the parent.
#' @param R Truncation distance for the difference process.
#' @param child.dist The distribution of the number of children
#'     generated by a randomly selected parent. For a Poisson
#'     distribution, use \code{"pois"}; for a binomial distribution,
#'     use \code{"binomx"}, where \code{"x"} is replaced by the fixed
#'     value of the number of independent trials (e.g.,
#'     \code{"binom5"} for a Binomial(5, p) distribution, and
#'     \code{"binom50"} for a Binomial(50, p) distribution); and
#'     \code{"twocamera"} for a child distribution appropriate for a
#'     two-camera aerial survey.
#' @param child.info A list of further information that is required
#'     about the distribution for the number of children generated by
#'     parents. See `Details'.
#' @param sibling.list An optional list that comprises (i) a component
#'     named \code{sibling.mat}, containing a matrix such that the jth
#'     entry in the ith row is \code{TRUE} if the ith and jth points
#'     are known siblings, \code{FALSE} if they are known nonsiblings,
#'     and \code{NA} if their sibling status is not known; (ii) alpha,
#'     providing the probability that a sibling is successfully
#'     identified as a sibling; and (iii) beta, providing the
#'     probability that a nonsibling is successfully identified as a
#'     nonsibling. For multi-pattern fitting, this object must be a
#'     list of such lists, one for each pattern.
#' @param edge.correction The method used for the correction of edge
#'     effects. Either \code{"pbc"} for periodic boundary conditions,
#'     or \code{"buffer"} for a buffer-zone correction.
#' @param start A named vector of starting values for the model
#'     parameters.
#' @param bounds A list with named components. Each component should
#'     be a vector of length two, giving the upper and lower bounds
#'     for the named parameter.
#' @param trace Logical; if \code{TRUE}, parameter values are printed
#'     to the screen for each iteration of the optimisation procedure.
#' @param use.bobyqa Logial; if \code{TRUE} the \link{bobyqa} function
#'     is used for optimisation. Otherwise the \link{nlminb} function
#'     is used. Note that \link{bobyqa} seems to be less stable than
#'     \link{nlminb}, but does not require calculation of the Palm
#'     likelihood's partial derivatives.
#'
#' @inheritParams fit.ns
#' 
#' @return An R6 reference class object.
#'
#' @seealso Use \link{coef.palm} to extract estimated parameters, and
#'     \link{plot.palm} to plot the estimated Palm intensity
#'     function. Use \link{boot.palm} to run a parametric bootstrap,
#'     allowing calculation of standard errors and confidence
#'     intervals.
#'
#' @seealso See \link{sim.ns} to simulate from a Neyman-Scott point
#'     process.
#'
#' @examples
#' ## Fitting model to example data.
#' fit <- fit.ns(example.2D, lims = rbind(c(0, 1), c(0, 1)), R = 0.5)
#' ## Printing estimates.
#' coef(fit)
#' ## Plotting the estimated Palm intensity.
#' plot(fit)
#' \dontrun{
#' ## Simulating data and fitting additional models.
#' set.seed(1234)
#' ## One-dimensional Thomas process.
#' data.thomas <- sim.ns(c(D = 10, lambda = 5, sigma = 0.025), lims = rbind(c(0, 1)))
#' ## Fitting a model to these data.
#' fit.thomas <- fit.ns(data.thomas$points, lims = rbind(c(0, 1)), R = 0.5)
#' ## Three-dimensional Matern process.
#' data.matern <- sim.ns(c(D = 10, lambda = 10, tau = 0.1), disp = "uniform",
#'                       lims = rbind(c(0, 1), c(0, 2), c(0, 3)))
#' ## Fitting a model to these data.
#' fit.matern <- fit.ns(data.matern$points, lims = rbind(c(0, 1), c(0, 2), c(0, 3)),
#'                      R = 0.5, disp = "uniform")
#' }
#'
#' @export
fit.ns <- function(points, lims, R, disp = "gaussian", child.dist = "pois", child.info = NULL,
                   sibling.list = NULL, edge.correction = "pbc", start = NULL, bounds = NULL,
                   use.bobyqa = FALSE, trace = FALSE){
    classes.list <- setup.classes(fit = TRUE, family = "ns",
                                  family.info = list(child.dist = child.dist,
                                                     child.info = child.info,
                                                     disp = disp,
                                                     sibling.list = sibling.list),
                                  fit.info = list(edge.correction = edge.correction,
                                                  use.bobyqa = use.bobyqa))
    obj <- create.obj(classes = classes.list$classes, points = points, lims = lims, R = R,
                      child.list = classes.list$child.list, parent.locs = NULL,
                      sibling.list = sibling.list, trace = trace, start = start, bounds = bounds)
    obj$fit()
    obj
}

#' Simulating points from a Neyman-Scott point process
#'
#' Generates points from a Neyman-Scott point process using parameters
#' provided by the user.
#'
#' For a list of possible parameter names, see \link{fit.ns}.
#' 
#' The \code{"child.info"} argument is required when \code{child.dist}
#' is set to \code{"twocamera"}. It must be a list that comprises (i) a
#' component named \code{w}, providing the halfwidth of the detection
#' zone; (ii) a component named \code{b}, providing the halfwidth of
#' the survey area; (iii) a component named \code{l}, providing the
#' time lag between cameras (in seconds); and (iv) a component named
#' \code{tau}, providing the mean dive-cycle duration. See Stevenson,
#' Borchers, and Fewster (in press) for details.
#'
#' @references Stevenson, B. C., Borchers, D. L., and Fewster,
#'     R. M. (in press) Cluster capture-recapture to account for
#'     identification uncertainty on aerial surveys of animal
#'     populations. \emph{Biometrics}.
#'
#' @param pars A named vector containing the values of the parameters
#'     of the process that generates the points.
#' 
#' @inheritParams fit.ns
#' @param parents An optional matrix containing locations of
#'     parents. If this is provided, then the parameter \code{D} is
#'     not required in \code{pars}. If this is not provided, then
#'     parents are generated from a homogeneous Poisson point process
#'     with intensity \code{D}.
#'
#' @return A list. The first component gives the Cartesian coordinates
#'     of the generated points. The second component returns the
#'     parent locations. A third component may provide sibling
#'     information.
#'
#' @examples
#' ## Simulating from a one-dimensional Thomas process.
#' data.thomas <- sim.ns(c(D = 10, lambda = 5, sigma = 0.025), lims = rbind(c(0, 1)))
#' ## Simulating from a three-dimensional Matern process.
#' data.matern <- sim.ns(c(D = 10, lambda = 10, tau = 0.1), disp = "uniform",
#'                       lims = rbind(c(0, 1), c(0, 2), c(0, 3)))
#' 
#' @export
sim.ns <- function(pars, lims, disp = "gaussian", child.dist = "pois", parents = NULL, child.info = NULL){
    classes.list <- setup.classes(fit = FALSE, family = "ns", family.info = list(child.dist = child.dist,
                                                                                 child.info = child.info,
                                                                                 parent.locs = parents,
                                                                                 disp = disp),
                                  fit.info = NULL)
    obj <- create.obj(classes = classes.list$classes, points = NULL, lims = lims, R = NULL,
                      child.list = classes.list$child.list, parent.locs = classes.list$parent.locs,
                      sibling.list = NULL, trace = NULL, start = NULL, bounds = NULL)
    obj$simulate(pars)
}

#' Fitting a model to a void point process
#'
#' Estimates parameters for a void point process by maximising the
#' Palm likelihood. This approach was first proposed by Tanaka et
#' al. (2008) for two-dimensional Thomas processes. Generalisation to
#' d-dimensional void processes was made by Jones-Todd et al. (in press).
#'
#' Parameters to estimate are as follows:
#' \itemize{
#'   \item \code{Dc}, the baseline density of points prior to the deletion process.
#'
#'   \item \code{Dp}, the density of unobserved parents that cause voids.
#'
#'   \item \code{tau}, the radius of the deletion process centred at each parent.
#' }
#'
#' @references Jones-Todd, C. M., Caie, P., Illian, J. B., Stevenson,
#'     B. C., Savage, A., Harrison, D. J., and Bown, J. L. (in
#'     press). Identifying prognostic structural features in tissue
#'     sections of colon cancer patients using point pattern
#'     analysis. \emph{Statistics in Medicine}.
#' @references Tanaka, U., Ogata, Y., and Stoyan, D. (2008) Parameter
#'     estimation and model selection for Neyman-Scott point
#'     processes. \emph{Biometrical Journal}, \strong{50}: 43--57.
#'
#' @inheritParams fit.ns
#'
#' @return An R6 reference class object.
#' 
#' @seealso Use \link{coef.palm} to extract estimated parameters, and
#'     \link{plot.palm} to plot the estimated Palm intensity
#'     function. Use \link{boot.palm} to run a parametric bootstrap,
#'     allowing calculation of standard errors and confidence
#'     intervals.
#'
#' @seealso See \link{sim.void} to simulate from a void process.
#'
#' @examples
#' \dontrun{
#' set.seed(1234)
#' ## Simulating a two-dimensional void process.
#' void.data <- sim.void(c(Dc = 1000, Dp = 10, tau = 0.05), rbind(c(0, 1), c(0, 1)))
#' ## Fitting model.
#' fit <- fit.void(void.data$points, rbind(c(0, 1), c(0, 1)), R = 0.5)
#' }
#' 
#' @export
fit.void <- function(points, lims, R, edge.correction = "pbc", start = NULL, bounds = NULL,
                     use.bobyqa = FALSE, trace = FALSE){
    classes.list <- setup.classes(fit = TRUE, family = "void", family.info = NULL,
                                  fit.info = list(edge.correction = edge.correction,
                                                  use.bobyqa = use.bobyqa))
    obj <- create.obj(classes = classes.list$classes, points = points, lims = lims, R = R,
                      child.list = NULL, parent.locs = NULL, sibling.list = NULL,
                      trace = trace, start = start, bounds = bounds)
    obj$fit()
    obj
}

#' Simulating points from a void point process.
#'
#' Generates points from a void point process using parameters provided by the user.
#'
#' For a list of possible parameter names, see \link{fit.ns}.
#'
#' @inheritParams fit.void
#' @inheritParams sim.ns
#' @param parents An optional matrix containing locations of
#'     parents. If this is provided, then the parameter \code{D} is
#'     not required in \code{pars}. If this is not provided, then
#'     parents are generated from a homogeneous Poisson point process
#'     with intensity \code{Dp}.
#'
#' @return A list. The first component gives the Cartesian coordinates
#'     of the generated points. The second component returns the
#'     parent locations.
#'
#' @examples
#' ## Two-dimensional void process.
#' void.data <- sim.void(c(Dc = 1000, Dp = 10, tau = 0.05), rbind(c(0, 1), c(0, 1)))
#' ## Plotting the data.
#' plot(void.data$points)
#' points(void.data$parents, pch = 16, col = "red")
#'
#' @export
sim.void <- function(pars, lims, parents = NULL){
    classes.list <- setup.classes(fit = FALSE, family = "void", family.info = list(parent.locs = parents),
                                  fit.info = NULL)
    obj <- create.obj(classes = classes.list$classes, points = NULL, lims = lims, R = NULL,
                      child.list = NULL, parent.locs = classes.list$parent.locs,
                      sibling.list = NULL, trace = NULL, start = NULL, bounds = NULL)
    obj$simulate(pars)
}

setup.classes <- function(fit, family, family.info, fit.info){
    ## Initialising all classes to FALSE.
    use.fit.class <- FALSE
    use.nlminb.class <- FALSE
    use.bobyqa.class <- FALSE
    use.pbc.class <- FALSE
    use.buffer.class <- FALSE
    use.ns.class <- FALSE
    use.sibling.class <- FALSE
    use.poischild.class <- FALSE
    use.binomchild.class <- FALSE
    use.twocamerachild.class <- FALSE
    child.list <- NULL
    parent.locs <- NULL
    use.thomas.class <- FALSE
    use.matern.class <- FALSE
    use.void.class <- FALSE
    use.totaldeletion.class <- FALSE
    use.giveparent.class <- FALSE
    ## Sorting out fitting classes.
    if (fit){
        use.fit.class <- TRUE
        ## Sorting out optimisation class.
        if (fit.info$use.bobyqa){
            use.bobyqa.class <- TRUE
        } else {
            use.nlminb.class <- TRUE
        }
        ## Sorting out boundary condition class.
        if (fit.info$edge.correction == "pbc"){
            use.pbc.class <- TRUE
        } else if (fit.info$edge.correction == "buffer"){
            use.buffer.class <- TRUE
        } else {
            stop("Edge correction method not recognised; use either 'pbc' or 'buffer'.")
        }
    }
    ## Sorting out family (ns/void) class.
    if (family == "ns"){
        ## Stuff for ns class.
        use.ns.class <- TRUE
        if (!is.null(family.info$sibling.list)){
            use.sibling.class <- TRUE
        }
        ## Sorting out child distribution class.
        if (family.info$child.dist == "pois"){
            use.poischild.class <- TRUE
        } else if (substr(family.info$child.dist, 1, 5) == "binom"){
            use.binomchild.class <- TRUE
            n <- as.numeric(substr(family.info$child.dist, 6, nchar(family.info$child.dist)))
            child.list <- list(size = n)
        } else if (family.info$child.dist == "twocamera"){
            use.twocamerachild.class <- TRUE
            child.list <- list(twocamera.w = family.info$child.info$w,
                               twocamera.b = family.info$child.info$b,
                               twocamera.l = family.info$child.info$l,
                               twocamera.tau = family.info$child.info$tau)
        } else {
            stop("Only 'pois', 'binomx', or 'twocamera' can currently be used for 'child.dist'.")
        }
        ## Sorting out dispersion class.
        if (family.info$disp == "gaussian"){
            use.thomas.class <- TRUE
        } else if (family.info$disp == "uniform"){
            use.matern.class <- TRUE
        } else {
            stop("Dispersion type not recognised; use either 'gaussian' or 'uniform'.")
        }
    } else if (family == "void"){
        ## Stuff for void class.
        use.void.class <- TRUE
        use.totaldeletion.class <- TRUE
    }
    ## Sorting out parent location class.
    if (!is.null(family.info$parent.locs)){
        use.giveparent.class <- TRUE
        parent.locs <- family.info$parent.locs
    }
    classes <- c("fit"[use.fit.class],
                 "bobyqa"[use.bobyqa.class],
                 "nlminb"[use.nlminb.class],
                 "pbc"[use.pbc.class],
                 "buffer"[use.buffer.class],
                 "ns"[use.ns.class],
                 "sibling"[use.sibling.class],
                 "poischild"[use.poischild.class],
                 "binomchild"[use.binomchild.class],
                 "twocamerachild"[use.twocamerachild.class],
                 "thomas"[use.thomas.class],
                 "matern"[use.matern.class],
                 "void"[use.void.class],
                 "totaldeletion"[use.totaldeletion.class],
                 "giveparent"[use.giveparent.class])
    list(classes = classes, child.list = child.list, parent.locs = parent.locs)
}

#' Estimation of animal density from two-camera surveys.
#'
#' Estimates animal density (amongst other parameters) from two-camera
#' aerial surveys. This conceptualises sighting locations as a
#' Neyman-Scott point pattern.
#'
#' This function is simply a wrapper for \code{fit.ns}, and
#' facilitates the fitting of the model proposed by Stevenson,
#' Borchers, and Fewster (in press). This function presents the
#' parameter \code{D.2D} (two-dimensional cetacean density in
#' cetaceans per square km) rather than \code{D} for enhanced
#' interpretability.
#'
#' For further details on the cluster capture-recapture estimation
#' approach, see Fewster, Stevenson and Borchers (2016).
#'
#' @references Fewster, R. M., Stevenson, B. C., and Borchers,
#'     D. L. (2016) Trace-contrast methods for capture-recapture
#'     without capture histories. \emph{Statistical Science},
#'     \strong{31}: 245--258.
#' @references Stevenson, B. C., Borchers, D. L., and Fewster,
#'     R. M. (in press) Cluster capture-recapture to account for
#'     identification uncertainty on aerial surveys of animal
#'     populations. \emph{Biometrics}.
#'
#' @param points A vector (or single-column matrix) containing the
#'     distance along the transect that each detection was made.
#' @param cameras An optional vector containing the camera ID (either
#'     \code{1} or \code{2}) that made the corresponding detection in
#'     \code{points}.
#' @param d The length of the transect flown (in km).
#' @param w The distance from the transect to which detection of
#'     individuals on the surface is certain. This is equivalent to
#'     the half-width of the detection zone.
#' @param b The distance from the transect to the edge of the area of
#'     interest. Conceptually, the distance between the transect and
#'     the furthest distance a whale could be on the passing of the
#'     first camera and plausibly move into the detection zone by the
#'     passing of the second camera.
#' @param l The lag between cameras (in seconds).
#' @param tau Mean dive-cycle duration (in seconds).
#' @param R Truncation distance (see \link{fit.ns}).
#' @inheritParams fit.ns
#'
#' @return An R6 reference class object.
#'
#' @seealso Use \link{coef.palm} to extract estimated parameters, and
#'     \link{plot.palm} to plot the estimated Palm intensity
#'     function. Use \link{boot.palm} to run a parametric bootstrap,
#'     allowing calculation of standard errors and confidence
#'     intervals.
#'
#' @seealso See \link{sim.twocamera} to simulate sightings from a
#'     two-camera aerial survey.
#'
#' @examples
#' ## Fitting model.
#' fit <- fit.twocamera(points = example.twocamera$points, cameras = example.twocamera$cameras,
#'                      d = 500, w = 0.175, b = 0.5, l = 20, tau = 110, R = 1)
#' ## Printing estimates.
#' coef(fit)
#' ## Plotting the estimated Palm intensity.
#' plot(fit)
#' 
#' @export
fit.twocamera <- function(points, cameras = NULL, d, w, b, l, tau, R,
                         edge.correction = "pbc", start = NULL,
                         bounds = NULL, trace = FALSE){
    if (is.vector(points)){
        points <- matrix(points, ncol = 1)
    }
    if (is.null(cameras)){
        sibling.list <- NULL
    } else {
        sibling.list <- siblings.twocamera(cameras)
    }
    if (is.null(bounds)){
        bounds <- list(sigma = c(0, min(R, b/3)))
    } else if (!any(names(bounds) == "sigma")){
        bounds[["sigma"]] <- c(0, min(R, b/3))
    }
    fit.ns(points = points, lims = rbind(c(0, d)), R = R,
           child.dist = "twocamera",
           child.info = list(w = w, b = b, l = l, tau = tau),
           sibling.list = sibling.list, edge.correction = edge.correction,
           start = start, bounds = bounds, trace = trace)
}

#' Simulating data from two-camera aerial surveys.
#'
#' @param pars A vector containing elements named \code{D.2D},
#'     \code{kappa}, and \code{sigma}, providing values of animal
#'     density (animals per square km), average duration of surface
#'     phase (s), and dispersion (km).
#' @param parents An optional vector containing the parent locations
#'     for all animals within the area of interest, given in distance
#'     along the transect (in km). If this is provided, then the
#'     parameter \code{D.2D} is not required in \code{pars}. If this
#'     is not provided, then parent locations are generated from a
#'     homogeneous Poisson point process with intensity \code{D.2D}.
#' @inheritParams fit.twocamera
#'
#' @return A list. The first component gives the distance along the
#'     transect of detected individuals. The second gives the parent
#'     locations. The third identifies which parent location generated
#'     each detected individual. The fourth gives the distance from
#'     the transect centre line of the detection location. The fifth
#'     provides observed sibling information.
#'
#' @examples
#' twocamera.data <- sim.twocamera(c(D.2D = 1.3, kappa = 27, sigma = 0.02), d = 500,
#'                                 w = 0.175, b = 0.5, l = 20, tau = 110)
#'
#' @export
sim.twocamera <- function(pars, d, w, b, l, tau, parents = NULL){
    if (!is.null(parents)){
        parents <- matrix(parents, ncol = 1)
    }
    family.info <- list(child.dist = "twocamera",
                        child.info = list(w = w, b = b, l = l, tau = tau),
                        parent.locs = parents, disp = "gaussian")
    classes.list <- setup.classes(fit = FALSE, family = "ns",
                                  family.info = family.info)
    obj <- create.obj(classes = classes.list$classes, points = NULL, lims = rbind(c(0, d)),
                      R = NULL, child.list = classes.list$child.list,
                      parent.locs = classes.list$parent.locs, sibling.list = NULL,
                      trace = NULL, bounds = NULL)
    obj$simulate(pars)
}

#' Bootstrapping for fitted models
#'
#' Carries out a parametric bootstrap procedure for models fitted
#' using the \code{palm} package.
#'
#' @return The original model object containing additional information
#'     from the bootstrap procedure. These are accessed by functions
#'     such as \link{summary.palm} and \link{confint.palm}. The
#'     bootstrap parameter estimates can be found in the \code{boots}
#'     component of the returned object.
#'
#' @param fit A fitted object.
#' @param N The number of bootstrap resamples.
#' @param prog Logical, if \code{TRUE}, a progress bar is printed to
#'     the console.
#'
#' @examples
#' ## Fit model.
#' fit <- fit.ns(example.2D, lims = rbind(c(0, 1), c(0, 1)), R = 0.5)
#' ## Carry out bootstrap.
#' fit <- boot.palm(fit, N = 100)
#' ## Inspect standard errors and confidence intervals.
#' summary(fit)
#' confint(fit)
#' ## Estimates are very imprecise---these data were only used as
#' ## they can be fitted and bootstrapped quickly for example purposes.
#' 
#' @export
boot.palm <- function(fit, N, prog = TRUE){
    fit$boot(N, prog)
    fit
}



