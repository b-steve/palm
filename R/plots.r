#' Plotting the empirical Palm intensity
#'
#' Plots the Palm intensity from empirical data.
#'
#' @return A plot showing the empirical Palm intensity.
#' 
#' @inheritParams fit.ns
#' @param breaks The (approximate) number of points plotted.
#' @param xlim The x-axis limits for the plot.
#' @param add Logical, if \code{TRUE} then the line is added to a
#' plot.
#' @param ... Graphical parameters (e.g., to be passed to
#' \link{par}().
#'
#' @examples
#' empirical.palm(example.2D, lims = rbind(c(0, 1), c(0, 1)), breaks = 25)
#' 
#' @export
empirical.palm <- function(points, lims, breaks = 50, xlim = NULL, add = FALSE, ...){
    error.dims(points, lims)
    n.dims <- ncol(points)
    dists <- pbc_distances(points = points, lims = lims)
    n.points <- nrow(points)
    if (is.null(xlim)){
        xlim <- 0.5*min(apply(lims, 1, diff))
    }
    if (is.null(breaks)){
        breaks <- 50
    }
    midpoints <- seq(0, max(xlim), length.out = breaks)
    midpoints <- midpoints[-length(midpoints)]
    h <- diff(midpoints[c(1, 2)])
    midpoints[1] <- midpoints[1] + h/2
    intensities <- numeric(length(midpoints))
    for (i in 1:length(midpoints)){
        halfwidth <- ifelse(i == 1, 0.25*h, 0.5*h)
        n.interval <- sum(dists <= (midpoints[i] + halfwidth)) -
            sum(dists <= (midpoints[i] - halfwidth))
        area <- Vd(midpoints[i] + halfwidth, n.dims) -  Vd(midpoints[i] - halfwidth, n.dims)
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
    lines(midpoints, intensities, ...)
}

#' Plotting an estimated Palm intensity function.
#'
#' Plots a fitted Palm intensity function from an object returned by
#' \link{fit.ns}().
#'
#' @param x A fitted model from \link{fit.ns}().
#' @param plot.empirical Logical, if \code{TRUE} then the empirical
#' Palm intensity is also plotted.
#' @param xlim The x-axis limits for the plot.
#' @param ylim The y-axis limits for the plot.
#' @inheritParams empirical.palm
#' @param ... Graphical parameters (e.g., to be passed to
#' \link{par}().
#'
#' @method plot nspp
#'
#' @examples
#' ## Fitting a model.
#' fit <- fit.ns(example.2D, lims = rbind(c(0, 1), c(0, 1)), R = 0.5)
#' ## Plotting.
#' plot(fit, plot.empirical = TRUE, breaks = 30, ylim = c(0, 200))
#' 
#' @export
plot.nspp <- function(x, plot.empirical = FALSE, breaks = NULL,
                      xlim = NULL, ylim = NULL, ...){
    Dc <- x$pars["Dc"]
    nu <- x$pars["nu"]
    sigma <- x$pars["sigma"]
    R <- x$args$R
    if (is.null(xlim)){
        xlim <- c(0, R)
    }
    analytic.palm(Dc, nu, sigma, ncol(x$args$points), xlim = xlim,
                  ylim = ylim, lty = ifelse(plot.empirical, "dashed", "solid"), ...)
    if (plot.empirical){
        empirical.palm(x$args$points, x$args$lims,
                       breaks = breaks, xlim = xlim, add = TRUE)
    }
}

## Plots the analytic Palm intensity.
analytic.palm <- function(Dc, nu, sigma, n.dims, xlim = c(0, 1), ylim = NULL, add = FALSE, ...){
    xx <- seq(xlim[1], xlim[2], length.out = 500)
    yy <- palm.intensity(xx, Dc, nu, sigma, n.dims)
    if (is.null(ylim)){
        ylim <- c(0, max(yy))
    }
    if (!add){
        par(xaxs = "i")
        plot.new()
        plot.window(xlim = xlim, ylim = ylim)
        box()
        abline(h = 0, col = "lightgrey")
        axis(1)
        axis(2)
    }
    lines(xx, yy, ...)
}
