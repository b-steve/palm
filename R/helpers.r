## Surface volume of d-dimensional hypersphere with radius r.
Sd <- function(r, d){
    d*pi^(d/2)*r^(d - 1)/gamma(d/2 + 1)
}

## Volume of d-dimensional hypersphere with radius r.
Vd <- function(r, d){
    pi^(d/2)*r^d/gamma(d/2 + 1)
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

## Calculation of two-camera detection probabilities.
twocamera.probs <- function(t, C, w, b, S, sigma){
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
    
    ## Probability of being in the transect for second camera, given in
    ## transect for first camera.
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

## Obtaining sibling matrix from camera IDs.
siblings.twocamera <- function(camera.id){
    n.points <- length(camera.id)
    out <- matrix(NA, nrow = n.points, ncol = n.points)
    for (i in 1:n.points){
        for (j in i:n.points){
            if (camera.id[i] == camera.id[j]){
                out[i, j] <- out[j, i] <- FALSE
            }
        }
    }
    list(sibling.mat = out, alpha = 0, beta = 0.5)
}
