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
    error.dims(points, lims)
    if (is.null(breaks)){
        breaks <- "Sturges"
    }
    n.dims <- ncol(points)
    dists <- pbc_distances(points = points, lims = lims)
    n.points <- nrow(points)
    hist.obj <- hist(dists, plot = FALSE, breaks = breaks)
    midpoints <- hist.obj$mids
    h <- hist.obj$breaks[2] - hist.obj$breaks[1]
    intensities <- numeric(length(midpoints))
    for (i in 1:length(midpoints)){
        n.interval <- sum(dists <= (midpoints[i] + 0.5*h)) -
            sum(dists <= (midpoints[i] - 0.5*h))
        area <- Vd(midpoints[i] + 0.5*h, n.dims) -  Vd(midpoints[i] - 0.5*h, n.dims)
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
    xx <- seq(0, R, length.out = 500)[-1]
    yy <- palm.intensity(xx, Dc, nu, sigma, ncol(x$args$points))
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
