#' Simulating Neyman-Scott point process data
#'
#' Simulates point locations from a Neyman-Scott point process.
#'
#' @return A matrix containing simulated point locations.
#'
#' @param pars A named vector (or list?) of parameter values. Required
#' parameters are \code{D}, the density of parent points, and
#' \code{sigma}, for the square root of the main diagonal of the
#' variance-covariance matrix of the multivariate normal distribution
#' of children point locations around their parent.
#' @param lims A matrix with two colums, corresponding to the upper
#' and lower limits of each dimension, respectively.
#' @param rchild A function for the generation of random
#' variables from the distribution of the number of children for each
#' parent. Parameter values for this distribution are specified in
#' \code{...}. The first argument of this function must be named
#' \code{n}, and set the number of random values to generate. Built in
#' functions (e.g., \code{rpois}) are suitable.
#' @param non.siblings An argument specifying the type of sibling
#' information to return (not sure what this should look like yet).
#' @param plot.points Logical, if \code{TRUE}, simulated parent and children
#' point locations will be plotted (only for two dimensions).
#' @param plot.empirical Logical, if \code{TRUE}, the empirical Palm
#' intensity is plotted, with the true Palm intensity from the
#' provided parameters overlain. The latter is only approximate as the
#' mean and variance of the number of children per parent are
#' calculated via simulation from rchild.
#' @param ... Further parameters for rchild.dist.
#'
#' @export
sim.ns <- function(pars = NULL, lims = rbind(c(0, 1), c(0, 1)), rchild = rpois, non.siblings = NULL, plot.points = FALSE, plot.empirical = FALSE, ...){
    ## Allowing lims to be a vector if only one dimension.
    if (!is.matrix(lims)){
        lims <- matrix(lims, nrow = 1)
    }
    ## Parameter values.
    D <- pars["D"]
    sigma <- pars["sigma"]
    ## Number of dimensions.
    n.dims <- nrow(lims)
    ## Calculating survey area.
    area <- prod(apply(lims, 1, diff))
    ## Calculating expected number of parents.
    expected.n.parents <- D*area
    ## Generating number of parents.
    n.parents <- rpois(n = 1, lambda = expected.n.parents)
    ## Error if no parents generated.
    if (n.parents == 0){
        stop("No parents generated.")
    }
    ## Generating parent location points.
    parent.locs <- matrix(0, nrow = n.parents, ncol = n.dims)
    for (i in 1:n.dims){
        parent.locs[, i] <- runif(n.parents, lims[i, 1], lims[i, 2])
    }
    ## Generating the number of children spawned by each parent.
    n.childs <- rchild(n = n.parents, ...)
    ## Total number of children.
    n.children <- sum(n.childs)
    ## Error if no children generated.
    if (n.children == 0){
        stop("No children generated.")
    }
    ## Generating children dispersion from parent.
    child.disp <- rmvnorm(n = n.children, mean = rep(0, n.dims),
                          sigma = sigma^2*diag(n.dims))
    child.locs <- matrix(rep(parent.locs, times = rep(n.childs, n.dims)), ncol = n.dims) + child.disp
    ## Adjusting for periodic boundary constraints.
    child.locs <- pbc.fix(child.locs, lims)
    if (plot.points){
        if (n.dims == 2){
            plot.new()
            plot.window(xlim = lims[1, ], ylim = lims[2, ])
            box()
            axis(1)
            axis(2)
            points(parent.locs, pch = 4, lwd = 2)
            points(child.locs)
        } else {
            warning("Plotting points only implemented for two dimensions.")
        }
        if (plot.empirical){
            warning("Both 'plot.points' and 'plot.empirical' are TRUE, the latter is being ignored.")
        }
    } else if (plot.empirical){
        empirical.palm(child.locs, lims)
        rs <- rchild(10000, ...)
        rs.mean <- mean(rs)
        rs.var <- var(rs)
        Dc <- D*rs.mean
        nu <- (rs.var + rs.mean^2)/rs.mean - 1
        analytic.palm(Dc, nu, sigma, n.dims, c(0, 1), add = TRUE)
    }
    child.locs
}

#' Simulating two-plane whale survey data
#'
#' Simulates observed whale locations and plane IDs from a two-plane
#' whale survey.
#'
#' @return A list containing observed whale locations and associated
#' plane IDs.
#'
#' @param pars A named vector of parameter values. Required parameters
#' are \code{D}, whale density, \code{sigma}, whale movement,
#' \code{p01}, the probability that a whale is on the surface when the
#' second plane flies over, given that it was submerged when the first
#' plane flew over, and \code{10}, the probability that a whale is
#' submerged with the second plane flies over, given that it was on
#' the surface when the first plane flew over.
#' @param lims The limits of the survey transect.
sim.twoplane <- function(pars, lims){
    ## Extracting parameters.
    D <- pars["D"]
    sigma <- pars["sigma"]
    p01 <- pars["p01"]
    p10 <- pars["p10"]
    p11 <- 1 - p10
    p00 <- 1 - p01
    ## Simulating number of whales.
    n.whales <- rpois(n = 1, lambda = D*diff(lims))
    ## Simulating whale locations for first flyover.
    pos.plane1 <- runif(n.whales, min = lims[1], max = lims[2])
    ## Simulating whale movement.
    movement <- sample(c(0, 1), size = n.whales, replace = TRUE)*
        sigma*sqrt(2)*sqrt(rchisq(n.whales, 1))
    ## Calculating whale locations for second flyover.
    pos.plane2 <- pos.plane1 + movement
    ## Simulating detections for first flyover.
    det.plane1 <- sample(c(TRUE, FALSE), size = n.whales, replace = TRUE,
                         prob = c(p01/(p10 + p01), p10/(p10 + p01)))
    ## Simulating detections for second flyover.
    det.plane2 <- logical(n.whales)
    det.plane2[det.plane1] <- sample(c(TRUE, FALSE), size = sum(det.plane1),
                                     replace = TRUE, prob = c(p11, p10))
    det.plane2[!det.plane1] <- sample(c(TRUE, FALSE), size = sum(!det.plane1),
                                      replace = TRUE, prob = c(p01, p00))
    ## Concatenating detection locations.
    points <- c(pos.plane1[det.plane1], pos.plane2[det.plane2])
    ## Fixing points for periodic boundary conditions.
    points <- pbc.fix(points, lims)
    points
}

    
