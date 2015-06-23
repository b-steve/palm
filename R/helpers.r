## Taking points and moving those outside the limits back into the limits.
pbc.fix <- function(points, lims){
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
    lim.diffs <- apply(lims, 1, diff)
    for (i in 1:n.points){
        for (j in 1:n.dims){
            ## Checking if ith point's jth dimension is within the limits.
            in.lims <- FALSE
            while (!in.lims){
                ## If too low, increment by the distance between the limits,
                ## If too high, decrement by the distance between the limits.
                if (points[i, j] < lims[j, 1]){
                    points[i, j] <- points[i, j] + lim.diffs[j]
                } else if (points[i, j] > lims[j, 2]){
                    points[i, j] <- points[i, j] - lim.diffs[j]
                }
                in.lims <- points[i, j] >= lims[j, 1] & points[i, j] <= lims[j, 2]
            }
        }
    }
    points
}

#' Plotting the empirical Palm intensity
#'
#' Plots the Palm intensity from empirical data.
#'
#' @return A plot showing the empirical Palm intensity.
#' 
#' @inheritParams fit.ns
#' @param breaks The (approximate) number of points plotted.
#' @param add Logical, if \code{TRUE} then the line is added to a
#' plot.
#' 
#' @export
empirical.palm <- function(points, lims, breaks = NULL, add = FALSE){
    if (is.null(breaks)){
        breaks <- "Sturges"
    }
    dists <- pbc_distances(points = points, lims = lims)
    n.points <- nrow(points)
    hist.obj <- hist(dists, plot = FALSE, breaks = breaks)
    midpoints <- hist.obj$mids
    h <- hist.obj$breaks[2] - hist.obj$breaks[1]
    intensities <- numeric(length(midpoints))
    for (i in 1:length(midpoints)){
        n.interval <- sum(dists <= (midpoints[i] + 0.5*h)) -
            sum(dists <= (midpoints[i] - 0.5*h))
        area <- pi*(midpoints[i] + 0.5*h)^2 - pi*(midpoints[i] - 0.5*h)^2
        intensities[i] <- n.interval/(n.points*area)
    }
    if (!add){
        par(xaxs = "i")
        par(yaxs = "i")
        plot.new()
        plot.window(xlim = c(0, midpoints[length(midpoints)]), ylim = c(0, max(intensities)*1.04))
        box()
        axis(1)
        axis(2)
    }
    lines(midpoints, intensities)
}

#' Plotting an estimated Palm intensity function.
#'
#' Plots a fitted Palm intensity function from an object returned by
#' \link{fit.ns}().
#'
#' @param x A fitted model from \link{fit.ns}().
#' @param emp Logical, if \code{TRUE} then the empirical Palm
#' intensity is also plotted.
#' @inheritParams empirical.palm
#' @param ... Other parameters (for S3 generic compatibility).
#'
#' @method plot nspp
#'
#' @export
plot.nspp <- function(x, emp = FALSE, breaks = NULL, ...){
    Dc <- x$pars["Dc"]
    nu <- x$pars["nu"]
    sigma <- x$pars["sigma"]
    R <- x$args$R
    xx <- seq(0, R, length.out = 500)
    yy <- ns.palm(xx, Dc, nu, sigma, n.dims = 2)
    par(xaxs = "i")
    plot.new()
    plot.window(xlim = c(0, R), ylim = c(0, max(yy)))
    box()
    abline(h = 0, col = "lightgrey")
    lines(xx, yy, col = "red")
    if (emp){
        empirical.palm(x$args$points, x$args$lims,
                       breaks = breaks, add = TRUE)
    }
    axis(1)
    axis(2)
}

## Analytic value for Dc given sigma and nu.
analytic.Dc <- function(nu, sigma, n.dists, n.points, R){
    (n.dists/n.points - nu + nu*exp((-R^2)/(4*sigma^2)))/(pi*R^2)
}
