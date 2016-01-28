## Taking points and moving those outside the limits back into the limits.
pbc.fix <- function(points, lims){
    ## Errors for inconsistent dimensions.
    if (!is.matrix(points)){
        points <- matrix(points, ncol = 1)
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
    2^(1 - d/2)*r^(d - 1)*exp(-r^2/(4*sigma^2))/((sigma*sqrt(2))^d*gamma(d/2))
}

## CDF of between-sibling distances.
Fd <- function(r, sigma, d){
    pgamma(r^2/(4*sigma^2), d/2)
}

## Note that Dc + nu/Sd(r, d)*fd(r, sigma, d) is a correct
## formulation, but the below cancels the r^(d - 1) from both Sd(r, d)
## and fd(r, sigma, d).
palm.intensity <- function(r, Dc, nu, sigma, d, siblings = NULL){
    Dc + nu/(pi^(d/2)*d/gamma(d/2 + 1))*
        2^(1 - d/2)*exp(-r^2/(4*sigma^2))/((sigma*sqrt(2))^d*gamma(d/2))
}

## Separate intensity function for known sibling information to
## optimise performance when there isn't any.
palm.intensity.siblings <- function(r, Dc, nu, sigma, d, siblings){
    ns.intensity <- Dc
    s.intensity <- nu/(pi^(d/2)*d/gamma(d/2 + 1))*
        2^(1 - d/2)*exp(-r^2/(4*sigma^2))/((sigma*sqrt(2))^d*gamma(d/2))
    siblings$ns.multipliers*ns.intensity +
        siblings$s.multipliers*s.intensity
}

## Takes matrix component of siblings and turns it into a vector that
## matches up to the distances computed by pbc_distance().
vectorise.siblings <- function(siblings){
    if (nrow(siblings$matrix) != ncol(siblings$matrix)){
        stop("Sibling matrix is not square.")
    }
    n.points <- nrow(siblings$matrix)
    n.comparisons <- n.points^2 - n.points
    vec <- numeric(n.comparisons)
    ns.multipliers <- numeric(n.comparisons)
    s.multipliers <- numeric(n.comparisons)
    k <- 1
    for (i in 1:(n.points - 1)){
        for (j in (i +1):n.points){
            vec[k] <- vec[k + 1] <- siblings$matrix[i, j]
            ## Multipliers for TRUE.
            if (!is.na(siblings$matrix[i, j]) & siblings$matrix[i, j]){
                ns.multipliers[k] <- ns.multipliers[k + 1] <- 0
                s.multipliers[k] <- s.multipliers[k + 1] <- siblings$pT
            }
            ## Multipliers for FALSE.
            if (!is.na(siblings$matrix[i, j]) & !siblings$matrix[i, j]){
                ns.multipliers[k] <- ns.multipliers[k + 1] <- siblings$pF
                s.multipliers[k] <- s.multipliers[k + 1] <- 0
            }
            ## Multipliers for NA.
            if (is.na(siblings$matrix[i, j])){
                ns.multipliers[k] <- ns.multipliers[k + 1] <- 1 - siblings$pF
                s.multipliers[k] <- s.multipliers[k + 1] <- 1 - siblings$pT
            }
            k <- k + 2
        }
    }
    list(vector = vec, ns.multipliers = ns.multipliers, s.multipliers = s.multipliers)
}

## Error function for incompatible dimensions.
error.dims <- function(points, lims){
    if (!is.matrix(points)){
        stop("Argument 'points' must be a matrix.")
    }
    if (!is.matrix(lims)){
        stop("Argument 'lims' must be a matrix.")
    }
    if (ncol(points) != nrow(lims)){
        stop("The number of columns in 'points' and 'lims' must both equal the number of dimensions.")
    }
}

## Functions below calculate partial derivatives for model parameters
## for two-dimensional processes with a Poisson distribution for the
## number of children.
dldD <- function(D, nu, sigma, n.points, dists, R){
    sum(1/(D + exp(-dists^2/(4*sigma^2))/(4*pi*sigma^2))) - n.points*pi*nu*R^2
}

dldnu <- function(D, nu, sigma, n.points, dists, R){
    length(dists)/nu - n.points*pi*D*R^2 - n.points + n.points*exp(-R^2/(4*sigma^2))
}

dldsigma <- function(D, nu, sigma, n.points, dists, R){
    sum((-n.points*nu*exp(-dists^2/(4*sigma^2))/(2*pi*sigma^3) +
             n.points*nu*dists^2*exp(-dists^2/(4*sigma^2))/(8*pi*sigma^5))/
        (n.points*D*nu + n.points*nu*exp(-dists^2/(4*sigma^2))/(4*pi*sigma^2))) +
        n.points*nu*(R^2)*exp(-R^2/(4*sigma^2))/(2*sigma^3)
}
