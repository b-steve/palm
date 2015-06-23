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
#' @param plot Logical, if \code{TRUE}, simulated parent and children
#' point locations will be plotted (only for two dimensions).
#' @param ... Further parameters for rchild.dist.
#'
#' @export
sim.ns <- function(pars = NULL, lims = rbind(c(0, 1), c(0, 1)), rchild = rpois, non.siblings = NULL, plot = FALSE, ...){
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
    child.locs <- matrix(rep(parent.locs, times = rep(n.childs, 2)), ncol = 2) + child.disp
    ## Adjusting for periodic boundary constraints.
    child.locs <- pbc.fix(child.locs, lims)
    if (plot){
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
    }
    child.locs
}
