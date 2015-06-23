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

## Analytic value for Dc given sigma and nu.
analytic.Dc <- function(nu, sigma, n.dists, n.points, R, d){
    (n.dists/n.points - nu*Fd(R, sigma, d))/Vd(R, d)
}

## Analytic value for nu given Dc and sigma.
analytic.nu <- function(Dc, sigma, n.dists, n.points, R, d){
    (n.dists/n.points - Dc*Vd(R, d))/Fd(R, sigma, d)
}

## Surface area of d-dimensional hypersphere with radius r.
Sd <- function(r, d){
    d*pi^(d/2)*r^(d - 1)/gamma(d/2 + 1)
}

## Volume of d-dimensional hypersphere with radius r.
Vd <- function(r, d){
    pi^(d/2)*r^d/gamma(d/2 + 1)
}

## PDF of between-sibling distances.
fd <- function(r, sigma, d){
    2^(1 - d/2)*(r/(sigma*sqrt(2)))^(d - 1)*exp(-r^2/(4*sigma^2))/(sigma*sqrt(2)*gamma(d/2))
}

## CDF of between-sibling distances.
Fd <- function(r, sigma, d){
    pgamma(r^2/(4*sigma^2), d/2)
}

palm.intensity <- function(r, Dc, nu, sigma, d){
    Dc + nu/Sd(r, d)*fd(r, sigma, d)
}
