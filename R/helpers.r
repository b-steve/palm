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
