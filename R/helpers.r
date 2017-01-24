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
analytic.Dc <- function(nu, sigma, n.dists, n.points, R, d, disp){
    Fd <- get.Fd(disp)
    (n.dists/n.points - nu*Fd(R, sigma, d))/Vd(R, d)
}

## Analytic value for nu given Dc and sigma.
analytic.nu <- function(Dc, sigma, n.dists, n.points, R, d, disp){
    Fd <- get.Fd(disp)
    (n.dists/n.points - Dc*Vd(R, d))/Fd(R, sigma, d)
}

## Surface volume of d-dimensional hypersphere with radius r.
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

## Closure for Fd function.
get.Fd <- function(disp){
    if (disp == "gaussian"){
        out <- function(r, sigma, d){
            pgamma(r^2/(4*sigma^2), d/2)
        }
    }
    else if (disp == "matern"){
        out <- function(r, sigma, d){

        }
    }
    out
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
vectorise.siblings <- function(siblings, edge.correction, buffer.keep = NULL){
    if (nrow(siblings$matrix) != ncol(siblings$matrix)){
        stop("Sibling matrix is not square.")
    }
    n.points <- nrow(siblings$matrix)
    if (edge.correction == "pbc"){
        n.comparisons <- n.points^2 - n.points
        vec <- numeric(n.comparisons)
        ns.multipliers <- numeric(n.comparisons)
        s.multipliers <- numeric(n.comparisons)
        k <- 1
        for (i in 1:(n.points - 1)){
            for (j in (i + 1):n.points){
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
    } else if (edge.correction == "buffer"){
        warning("Buffer edge correction with siblings is untested.")
        ns.mat <- matrix(0, nrow = n.points, ncol = n.points)
        s.mat <- matrix(0, nrow = n.points, ncol = n.points)
        ## Matrices for TRUE.
        ns.mat[!is.na(siblings$matrix) & siblings$matrix] <- 0
        s.mat[!is.na(siblings$matrix) & siblings$matrix] <- siblings$pT
        ## Matrices for FALSE.
        ns.mat[!is.na(siblings$matrix) & !siblings$matrix] <- siblings$pF
        s.mat[!is.na(siblings$matrix) & !siblings$matrix] <- 0
        ## Matrices for NA.
        ns.mat[is.na(siblings$matrix)] <- 1 - siblings$pF
        s.mat[is.na(siblings$matrix)] <- 1 - siblings$pT
        ## Turning into vectors.
        vec <- as.vector(siblings$matrix[buffer.keep])
        ns.multipliers <- as.vector(ns.mat[buffer.keep])
        s.multipliers <- as.vector(s.mat[buffer.keep])
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

## Calculation of two-plane detection probabilities.
twoplane.probs <- function(t, C, w, b, S, sigma){
    #browser(expr = S > 9.9)
    ## Up/down probabilities.

    ## Marginal probabilities.
    p.up <- S/C
    p.down <- 1 - S/C
    ## Conditional probabilities.
    ## Workaround for an edge case. 
    if (S >= C){
        p.up.up <- 1
        p.down.up <- 0
        p.up.down <- 1
        p.down.down <- 0
    } else {
        p.up.up <- S/C + ((C - S)/C)*exp(-t*(1/S + 1/(C - S)))
        p.down.up <- 1 - p.up.up
        p.up.down <- p.up*p.down.up/p.down
        p.down.down <- 1 - p.up.down
    }
    ## In/out probabilities.
    
    ## Probability of being in the transect for second plane, given in
    ## transect for first plane.
    f.y2.cond.y1 <- function(y1, w, sigma){
        F.pw <- pnorm(w, mean = y1, sd = sigma*sqrt(2))
        F.nw <- pnorm(-w, mean = y1, sd = sigma*sqrt(2))
        F.pw - F.nw
    }
    f.y1.cond.in <- function(y1, w){
        ifelse(abs(y1) < w, 1/(2*w), 0)
    }
    ## Marginal probabilities.
    p.in <- w/b
    p.out <- 1 - w/b
    ## Conditional probabilities.
    p.in.in <- integrate(function(y1, w, sigma) f.y2.cond.y1(y1, w, sigma)*
                             f.y1.cond.in(y1, w),
                         lower = -b, upper = b, w = w, sigma = sigma)$value
    p.out.in <- 1 - p.in.in
    p.in.out <- p.out.in*p.in/p.out
    p.out.out <- 1 - p.in.out

    ## Detection probabilities.

    ## Some inbetweeny stuff we need to calculate them.
    p.out.or.down <- p.out + p.down - p.out*p.down
    p.in.and.out.or.down <- (p.in.out*p.out*p.down +
                                 p.in.out*p.out*p.up +
                                     p.in.in*p.in*p.down)
    p.in.out.or.down <- p.in.and.out.or.down/p.out.or.down
    p.up.and.out.or.down <- (p.up.down*p.down*p.out +
                                 p.up.up*p.up*p.out +
                                     p.up.down*p.down*p.in)
    p.up.out.or.down <- p.up.and.out.or.down/p.out.or.down

    p.in.and.up.and.out.or.down <- (p.in.out*p.out*p.up.down*p.down +
                                        p.in.out*p.out*p.up.up*p.up +
                                        p.in.in*p.in*p.up.down*p.down)
    ## Workaround for an edge case.
    if (S == 0){
        p.in.up.and.out.or.down <- 0
    } else {
        p.in.up.and.out.or.down <- p.in.and.up.and.out.or.down/p.up.and.out.or.down
    }
    ## Conditional probabilities.
    p.11 <- p.up.up*p.in.in
    p.01 <- 1 - p.11
    p.10 <- p.in.up.and.out.or.down*p.up.out.or.down
    p.00 <- 1 - p.10

    ## Outputting conditional detection probabilities.
    out <- list(p.11 = p.11, p.01 = p.01, p.10 = p.10, p.00 = p.00,
                p.up = p.up, p.down = p.down, p.in = p.in, p.out = p.out,
                p.up.down = p.up.down, p.down.up = p.down.up,
                p.in.out = p.in.out, p.out.in = p.out.in)
    out
}

## Closurey function to make two-plane child.dist.
make.twoplane.child.dist <- function(t, C, w, b){
    mean.1D <- function(S, sigma){
        probs <- twoplane.probs(t, C, w, b, S, sigma)
        2*probs$p.10/(probs$p.10 + probs$p.01)
    }
    var.1D <- function(S, sigma){
        probs <- twoplane.probs(t, C, w, b, S, sigma)
        2*probs$p.10*probs$p.01*(2 - probs$p.10 - probs$p.01)/(probs$p.10 + probs$p.01)^2
    }
    list(mean = mean.1D, var = var.1D, sv = 0.5*C, bounds = c(1e-10, 0.999*C))  
}

## The incomplete beta function.
incomplete.beta <- function(x, a, b){
    pbeta(x, a, b)*beta(a, b)
}

logit <- function(p){
    log(p/(1 - p))
}

invlogit <- function(y){
    1/(1 + exp(-y))
}
