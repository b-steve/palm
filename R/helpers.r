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

empirical.palm <- function(points, lims, plot = FALSE){
    dists <- pbc_distances(points = points, lims = lims)
    n.points <- nrow(points)
    hist.obj <- hist(dists, plot = FALSE)
    midpoints <- hist.obj$mids
    h <- hist.obj$breaks[2] - hist.obj$breaks[1]
    intensities <- numeric(length(midpoints))
    for (i in 1:length(midpoints)){
        n.interval <- sum(dists <= (midpoints[i] + 0.5*h)) -
            sum(dists <= (midpoints[i] - 0.5*h))
        area <- pi*(midpoints[i] + 0.5*h)^2 - pi*(midpoints[i] - 0.5*h)^2
        intensities[i] <- n.interval/(n.points*area)
    }
    if (plot){
        plot(midpoints, intensities, ylim = c(0, max(intensities)), type = "l")
    }
    list(intensity = intensities, distance = midpoints)
}
