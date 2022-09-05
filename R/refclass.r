######
## General base class.
######

base.class.R6 <- R6Class("palm",
                         public = list(
                             ## Setting fields.
                             n.patterns = NULL,
                             lims = NULL,
                             lims.list = NULL,
                             vol = NULL,
                             vol.list = NULL,
                             dim = NULL,
                             dim.list = NULL,
                             classes = NULL,
                             ## Initialisation method.
                             initialize = function(lims.list, classes, ...){
                                 self$lims.list <- lims.list
                                 self$vol.list <- lapply(self$lims.list,
                                                         function(x) prod(apply(x, 1, diff)))
                                 self$dim.list <- lapply(self$lims.list,
                                                         function(x) nrow(x))
                                 self$n.patterns <- length(lims.list)
                                 self$classes <- classes
                             },
                             ## A method to set up a new pattern.
                             setup.pattern = function(pattern){
                                 if (pattern > self$n.patterns){
                                     stop("Pattern index exceeds the number of patterns")
                                 }
                                 self$lims <- self$lims.list[[pattern]]
                                 self$vol <-  self$vol.list[[pattern]]
                                 self$dim <- self$dim.list[[pattern]]
                             },
                             ## A method to simulate multiple patterns.
                             simulate = function(pars = self$par.fitted){
                                 out <- vector(mode = "list", length = self$n.patterns)
                                 for (i in 1:self$n.patterns){
                                     self$setup.pattern(i)
                                     out[[i]] <- self$simulate.pattern(pars)
                                 }
                                 if (self$n.patterns == 1){
                                     out <- out[[1]]
                                 }
                                 out
                             },
                             ## A method to simulate a single pattern.
                             simulate.pattern = function(pars){},
                             ## A method to trim points to the observation window.
                             trim.points = function(points, output.indices = FALSE){
                                 in.window <- rep(TRUE, nrow(points))
                                 for (i in 1:self$dim){
                                     in.window <- in.window & (self$lims[i, 1] <= points[, i] & self$lims[i, 2] >= points[, i])
                                 }
                                 if (output.indices){
                                     out <- which(in.window)
                                 } else {
                                     out <- points[in.window, , drop = FALSE]
                                 }
                                 out
                             }
                         ))

######
## Template for new classes.
######

## Replace CLASSNAME with class name, then add fields and methods.
set.CLASSNAME.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("CLASSNAME.inherit", class, envir = class.env)
    R6Class("palm_CLASSNAME",
            inherit = class.env$CLASSNAME.inherit,
            public = list(

            ))
}

######
## Class for models fitted to data.
######

set.fit.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("fit.inherit", class, envir = class.env)
    R6Class("palm_fit",
            inherit = class.env$fit.inherit,
            public = list(
                points.list = NULL,
                points = NULL,
                n.points = NULL,
                n.points.list = NULL,
                contrasts = NULL,
                contrasts.list = NULL,
                n.contrasts = NULL,
                n.contrasts.list = NULL,
                contrast.pairs = NULL,
                contrast.pairs.list = NULL,
                R = NULL,
                boots = NULL,
                trace = NULL,
                conv.code = NULL,
                converged = NULL,
                n.par = NULL,
                set.start = NULL,
                set.bounds = NULL,
                par.names = NULL,
                fixed.names = NULL,
                par.names.link = NULL,
                par.start = NULL,
                par.start.link = NULL,
                par.fixed = NULL,
                par.fixed.link = NULL,
                par.links = NULL,
                par.invlinks = NULL,
                par.fitted = NULL,
                par.fitted.link = NULL,
                par.lower = NULL,
                par.lower.link = NULL,
                par.upper = NULL,
                par.upper.link = NULL,
                pi.multiplier = NULL,
                ## Initialisation method.
                initialize = function(points.list, R, trace, start, bounds, ...){
                    super$initialize(...)
                    self$points.list <- points.list
                    ## Checking list of points and list of lims have the same length.
                    if (length(self$points.list) != length(self$lims.list)){
                        stop("The list 'points' must have the same number of components as 'lims'.")
                    }
                    self$R <- R
                    ## Creating numbers of points and contrasts for each pattern.
                    self$n.points.list <- lapply(points.list, nrow)
                    self$contrasts.list <- list()
                    self$contrast.pairs.list <- list()
                    self$n.contrasts.list <- list()
                    for (i in 1:self$n.patterns){
                        self$get.contrasts(i)
                    }
                    ## Starting out with the first pattern for start values and such.
                    self$setup.pattern(1)
                    self$trace <- trace
                    self$set.start <- start
                    self$set.bounds <- bounds
                    self$get.pars()
                    self$n.par <- length(self$par.start)
                    self$par.start.link <- self$link.pars(self$par.start)
                    self$get.link.bounds()
                },
                ## Overwriting the method to set up a new pattern.
                setup.pattern = function(pattern, do.contrasts = TRUE){
                    super$setup.pattern(pattern)
                    self$points <- self$points.list[[pattern]]
                    self$n.points <- self$n.points.list[[pattern]]
                    if (do.contrasts){
                        self$contrast.pairs <- self$contrast.pairs.list[[pattern]]
                        self$contrasts <- self$contrasts.list[[pattern]]
                        self$n.contrasts <- self$n.contrasts.list[[pattern]]
                    }
                },
                ## An empty method for getting contrasts.
                get.contrasts = function(pattern){
                    self$n.contrasts.list[[pattern]] <- length(self$contrasts.list[[pattern]])
                },
                ## A method for getting the parameters across all classes.
                get.pars = function(){
                    self$fetch.pars()
                },
                ## An empty method for fetching parameters from a particular class.
                fetch.pars = function(){
                },
                ## A method for adding new parameters.
                add.pars = function(name, link.name = NULL, link, invlink = NULL,
                                    start, lower, upper){
                    self$par.names <- c(self$par.names, name)
                    ## Sorting out inverse link.
                    if (identical(link, log)){
                        full.link <- function(pars) log(pars[name])
                        ## Default name and inverse for log link.
                        if (is.null(link.name)){
                            link.name <- paste("log", name, sep = ".")
                        }
                        if (is.null(invlink)){
                            full.invlink <- function(pars) exp(pars[link.name])
                        }
                    } else if (identical(link, logit)){
                        full.link <- function(pars) link(pars[name])
                        ## Default name and inverse for logit link.
                        if (is.null(link.name)){
                            link.name <- paste("logit", name, sep = ".")
                        }
                        if (is.null(invlink)){
                            full.invlink <- function(pars) invlogit(pars[link.name])
                        }
                    } else {
                        if (is.null(link.name)){
                            stop("Please provide link name.")
                        }
                        if (is.null(invlink)){
                            stop("Please provide inverse link function.")
                        }
                        full.link <- link
                        full.invlink <- invlink
                    }
                    self$par.links <- c(self$par.links, full.link)
                    names(self$par.links) <- self$par.names
                    self$par.invlinks <- c(self$par.invlinks, full.invlink)
                    names(self$par.invlinks) <- self$par.names
                    self$par.names.link <- c(self$par.names.link, link.name)
                    ## Overwriting start parameter, if one provided.
                    if (any(name == names(self$set.start))){
                        start <- self$set.start[name]
                    }
                    if (any(name == names(self$set.bounds))){
                        lower <- self$set.bounds[[name]][1]
                        upper <- self$set.bounds[[name]][2]
                    }
                    self$par.upper <- c(self$par.upper, upper)
                    names(self$par.upper) <- self$par.names
                    self$par.lower <- c(self$par.lower, lower)
                    ## Adjustment for start values on the bounds.
                    if (start >= upper){
                        start <- lower + 0.9*(upper - lower)
                    } else if (start <= lower){
                        start <- lower + 0.9*(upper - lower)
                    }
                    self$par.start <- c(self$par.start, start)
                    names(self$par.start) <- self$par.names
                    names(self$par.lower) <- self$par.names
                },
                ## A method to set the upper and lower parameter bounds on the link scale.
                get.link.bounds = function(){
                    self$par.lower.link <- self$link.pars(self$par.lower)
                    self$par.upper.link <- self$link.pars(self$par.upper)
                },
                ## A method for converting parameters to their link scales.
                link.pars = function(pars){
                    out <- numeric(self$n.par)
                    for (i in 1:self$n.par){
                        out[i] <- self$par.links[[i]](pars)
                    }
                    names(out) <- self$par.names.link
                    out
                },
                ## A method for converting parameters to their real scale.
                invlink.pars = function(pars){
                    out <- numeric(self$n.par)
                    for (i in 1:self$n.par){
                        out[i] <- self$par.invlinks[[i]](pars)
                    }
                    names(out) <- self$par.names
                    out
                },               

                ## A method for the objective function.
                obj.fun = function(link.pars, fixed.link.pars = NULL,
                                         est.names = NULL, fixed.names = NULL){
                    combined.pars <- c(link.pars, fixed.link.pars)
                    names(combined.pars) <- c(est.names, fixed.names)
                    link.pars <- combined.pars[self$par.names.link]
                    pars <- self$invlink.pars(link.pars)
                    obj.fun.components <- numeric(self$n.patterns)
                    ## Summing over all the patterns.
                    for (i in 1:self$n.patterns){
                        self$setup.pattern(i)
                        obj.fun.components[i] <- self$neg.log.palm.likelihood(pars)
                    }
                    out <- sum(obj.fun.components)
                    ll <- -out
                    if (self$trace){
                        for (i in 1:self$n.par){
                            cat(self$par.names[i], ": ", sep = "")
                            cat(pars[i], ", ", sep = "")
                        }
                        cat("Log-lik: ", ll, "\n", sep = "")
                    }
                    out
                },
                ## An empty method for the Palm intensity.
                palm.intensity = function(r, pars){},
                ## An empty method for the sum of the log Palm intensities.
                sum.log.intensities = function(pars){
                    sum(log(self$pi.multiplier*self$palm.intensity(self$contrasts, pars)))
                },
                ## A default method for the integral in the Palm likelihood.
                palm.likelihood.integral = function(pars){
                    f <- function(r, pars){
                        self$palm.intensity(r, pars)*Sd(r, self$dim)
                    }
                    -self$pi.multiplier*integrate(f, lower = 0, upper = self$R, pars = pars)$value
                },
                ## A default method for the log of the Palm likelihood function.
                log.palm.likelihood = function(pars){
                    self$sum.log.intensities(pars) + self$palm.likelihood.integral(pars)
                },
                ## A method for the negative Palm likelihood function.
                neg.log.palm.likelihood = function(pars){
                    -self$log.palm.likelihood(pars)
                },
                ## An empty method for optimisation.
                palm.optimise = function(){},
                ## A method for non-optimisation.
                palm.not.optimise = function(){
                    self$converged <- FALSE
                    self$par.fitted.link <- rep(NA, self$n.par)
                    names(self$par.fitted.link) <- self$par.names.link
                    self$par.fitted <- rep(NA, self$n.par)
                    self$conv.code <- 5
                },
                ## A method for model fitting.
                fit = function(){
                    ## If there are contrasts, fit the model. Otherwise, return something else.
                    if (self$n.contrasts <= 1){
                        warning("Not enough observed points to compute parameter estimates.")
                        self$palm.not.optimise()
                    } else {
                        self$palm.optimise()
                    }
                },
                ## A method for bootstrapping.
                boot = function(N, prog = TRUE){
                    boots <- matrix(0, nrow = N, ncol = self$n.par)
                    ## Setting up progress bar.
                    if (prog){
                        pb <- txtProgressBar(min = 0, max = N, style = 3)
                    }
                    for (i in 1:N){
                        zero.points <- TRUE
                        while (zero.points){
                            sim.obj <- self$simulate()
                            if (self$n.patterns > 1){
                                sim.obj$points <- lapply(sim.obj, function(x) x$points)
                                sim.obj$sibling.list <- lapply(sim.obj, function(x) x$sibling.list)
                                zero.points <- sum(sapply(sim.obj$points, nrow)) == 0
                            } else {
                                zero.points <- nrow(sim.obj$points) == 0
                            }
                        }
                        ## Doesn't matter that sibling.list is non-null below, as sibling class is not passed.
                        obj.boot <- create.obj(classes = self$classes, points = sim.obj$points, lims = self$lims,
                                               R = self$R, child.list = self$child.list,
                                               sibling.list = sim.obj$sibling.list, trace = FALSE,
                                               start = self$par.fitted, bounds = self$bounds)
                        
                        obj.boot$fit()
                        if (obj.boot$converged){
                            boots[i, ] <- obj.boot$par.fitted  
                        } else {
                            boots[i, ] <- rep(NA, self$n.par)
                        }
                        ## Updating progress bar.
                        if (prog){
                            setTxtProgressBar(pb, i)
                        }
                    }
                    cat("\n")
                    self$boots <- boots
                    colnames(self$boots) <- self$par.names
                },
                ## A method for plotting.
                plot = function(xlim = NULL, ylim = NULL, show.empirical = TRUE, breaks = 50, ...){
                    if (show.empirical){
                        emp <- self$empirical(xlim, breaks)  
                    }
                    if (is.null(xlim)){
                        xlim <- self$get.xlim()
                    }
                    xx <- seq(xlim[1], xlim[2], length.out = 1000)
                    yy <- self$palm.intensity(xx, self$par.fitted)
                    if (is.null(ylim)){
                        if (show.empirical){
                            ylim <- c(0, max(c(emp$y, yy)))
                        } else {
                            ylim <- c(0, max(yy))
                        }
                    }
                    par(xaxs = "i")
                    plot(xx, yy, xlab = "", ylab = "", xlim = range(xx), ylim = ylim, type = "n", ...)
                    title(xlab = "r", ylab = "Palm intensity")
                    lines(xx, yy)
                    if (show.empirical){
                        lines(emp$x, emp$y, lty = "dashed")
                    }
                },
                ## A method for default xlim.
                get.xlim = function(){
                    c(0, self$R)
                },
                ## A method for empirical Palm intensity.
                empirical = function(xlim = NULL, breaks = 50){
                    if (is.null(xlim)){
                        xlim <- self$get.xlim()
                    }
                    dists <- pbc_distances(self$points, self$lims)
                    midpoints <- seq(0, xlim[2], length.out = breaks)
                    midpoints <- midpoints[-length(midpoints)]
                    h <- diff(midpoints[c(1, 2)])
                    midpoints[1] <- midpoints[1] + h/2
                    intensities <- numeric(length(midpoints))
                    for (i in 1:length(midpoints)){
                        halfwidth <- ifelse(i == 1, 0.25*h, 0.5*h)
                        n.interval <- sum(dists <= (midpoints[i] + halfwidth)) -
                            sum(dists <= (midpoints[i] - halfwidth))
                        area <- Vd(midpoints[i] + halfwidth, self$dim) -  Vd(midpoints[i] - halfwidth, self$dim)
                        intensities[i] <- n.interval/(self$pi.multiplier*area)
                    }
                    list(x = midpoints, y = intensities)
                }
            ))
}

######
## Class for bobyqa() optimisation.
######

set.bobyqa.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("bobyqa.inherit", class, envir = class.env)
    R6Class("palm_bobyqa",
            inherit = class.env$bobyqa.inherit,
            public = list(
                ## A method to minimise the negative Palm likelihood.
                palm.optimise = function(){
                    out <- try(bobyqa(self$par.start.link, self$obj.fun,
                                      lower = self$par.lower.link, upper = self$par.upper.link,
                                      fixed.link.pars = self$par.fixed.link,
                                      est.names = self$par.names.link,
                                      fixed.names = self$fixed.names.link,
                                      control = list(maxfun = 2000, npt = 2*self$n.par + 1)), silent = TRUE)
                    if (class(out)[1] == "try-error"){
                        self$palm.not.optimise()
                    } else {
                        if (out$ierr > 0){
                            warning("Failed convergence.")
                            self$converged <- FALSE
                        } else {
                            self$converged <- TRUE
                        }
                        self$par.fitted.link <- out$par
                        names(self$par.fitted.link) <- self$par.names.link
                        self$par.fitted <- self$invlink.pars(self$par.fitted.link)
                        self$conv.code <- out$ierr
                    }
                }
            ))
}
    
######
## Class for nlminb() optimisation.
######

set.nlminb.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("nlminb.inherit", class, envir = class.env)
    R6Class("palm_nlminb",
            inherit = class.env$nlminb.inherit,
            public = list(
                ## A method to minimise the negative Palm likelihood.
                palm.optimise = function(){
                    out <- nlminb(self$par.start.link, self$obj.fun,
                                  fixed.link.pars = self$par.fixed.link,
                                  est.names = self$par.names.link, fixed.names = self$fixed.names.link,
                                  control = list(eval.max = 2000, iter.max = 1500),
                                  lower = self$par.lower.link, upper = self$par.upper.link)
                    if (out$convergence > 1){
                        warning("Failed convergence.")
                        self$converged <- FALSE
                    } else {
                        self$converged <- TRUE
                    }
                    self$par.fitted.link <- out$par
                    names(self$par.fitted.link) <- self$par.names.link
                    self$par.fitted <- self$invlink.pars(self$par.fitted.link)
                    self$conv.code <- out$convergence
                }
            ))
}
    
######
## Class for periodic boundary conditions.
######

set.pbc.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("pbc.inherit", class, envir = class.env)
    R6Class("palm_pbc",
            inherit = class.env$pbc.inherit,
            public = list(
                ## Overwriting the method to set up a new pattern.
                setup.pattern = function(pattern, do.contrasts = TRUE){
                    super$setup.pattern(pattern, do.contrasts)
                    self$pi.multiplier <- self$n.points/2
                },
                ## A method to generate contrasts.
                get.contrasts = function(pattern){
                    self$setup.pattern(pattern, do.contrasts = FALSE)
                    ## Saving which contrast applies to which pair of observations.
                    contrast.pairs <- matrix(0, nrow = (self$n.points^2 - self$n.points)/2, ncol = 2)
                    k <- 1
                    if (self$n.points > 1){
                        for (i in 1:(self$n.points - 1)){
                            for (j in (i + 1):self$n.points){
                                contrast.pairs[k, ] <- c(i, j)
                                k <- k + 1
                            }
                        }
                        contrasts <- pbc_distances(points = self$points, lims = self$lims)
                        self$contrast.pairs.list[[pattern]] <- contrast.pairs[contrasts <= self$R, ]
                        self$contrasts.list[[pattern]] <- contrasts[contrasts <= self$R]
                    } else {
                        self$contrast.pairs.list[[pattern]] <- contrast.pairs
                        self$contrasts[[pattern]] <- numeric(0)
                    }
                    super$get.contrasts(pattern)
                }
            ))
}

######
## Class for buffer-zone boundary conditions.
######

set.buffer.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("buffer.inherit", class, envir = class.env)
    R6Class("palm_buffer",
            inherit = class.env$buffer.inherit,
            public = list(
                ## Overwriting the method to set up a new pattern.
                setup.pattern = function(pattern, do.contrasts = TRUE){
                    super$setup.pattern(pattern, do.contrasts)
                    self$pi.multiplier <- self$n.points
                },
                ## A method to generate contrasts.
                get.contrasts = function(pattern){
                    self$setup.pattern(pattern, do.contrasts = FALSE)
                    ## Getting rid of contrasts between two external points.
                    buffer.keep <- buffer_keep(points = self$points, lims = self$lims,
                                               R = self$R)
                    contrasts <- as.vector(as.matrix(dist(self$points))[buffer.keep])
                    contrast.pairs.1 <- matrix(rep(1:self$n.points, self$n.points), nrow = self$n.points)
                    contrast.pairs.2 <- matrix(rep(1:self$n.points, self$n.points), nrow = self$n.points, byrow = TRUE)
                    contrast.pairs.1 <- as.vector(contrast.pairs.1[buffer.keep])
                    contrast.pairs.2 <- as.vector(contrast.pairs.2[buffer.keep])
                    contrast.pairs <- cbind(contrast.pairs.1, contrast.pairs.2)
                    ## Now truncating to contrasts less than R.
                    self$contrasts.list[[pattern]] <- contrasts[contrasts <= self$R] 
                    self$contrast.pairs.list[[pattern]] <- contrast.pairs[contrasts <= self$R, ]
                    super$get.contrasts(pattern)
                }
            ))
}


######
## Class for Neyman-Scott processes.
######

set.ns.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("ns.inherit", class, envir = class.env)
    R6Class("palm_ns",
            inherit = class.env$ns.inherit,
            public = list(
                child.list = NULL,
                ## Adding density parameter.
                fetch.pars = function(){
                    super$fetch.pars()
                    link <- function(pars){
                        log(pars["D"]*self$child.expectation(pars))
                    }
                    ## An inverse link that uses the other unlinked paramters.
                    invlink <- function(pars, out){
                        exp(pars["log.Dc"])/self$child.expectation(out)
                    }
                    self$add.pars(name = "D", link.name = "log.Dc", link = link, invlink = invlink,
                                  start = sqrt(self$n.points/self$vol), lower = 0, upper = Inf)
                },
                ## Overwriting the method for converting parameters to their real scale.
                invlink.pars = function(pars){
                    out <- numeric(self$n.par)
                    names(out) <- self$par.names
                    which.D <- which(self$par.names == "D")
                    for (i in (1:self$n.par)[-which.D]){
                        out[i] <- self$par.invlinks[[i]](pars)
                    }
                    out[which.D] <- self$par.invlinks[[which.D]](pars, out)
                    out
                },
                ## Overwriting simulation methods.
                simulate.pattern = function(pars = self$par.fitted){
                    parent.locs <- self$get.parents(pars)
                    n.parents <- nrow(parent.locs)
                    sim.n.children <- self$simulate.n.children(n.parents, pars)
                    n.children <- sim.n.children$n.children
                    sibling.list <- sim.n.children$sibling.list
                    child.locs <- matrix(0, nrow = sum(n.children), ncol = self$dim)
                    j <- 0
                    for (i in seq_along(n.children)){
                        if (n.children[i] > 0){
                            child.locs[(j + 1):(j + n.children[i]), ] <-
                                self$simulate.children(n.children[i], parent.locs[i, ], pars)
                            j <- j + n.children[i]
                        }
                    }
                    parent.ids <- rep(seq_along(n.children), times = n.children)
                    trimmed <- self$trim.points(child.locs, output.indices = TRUE)
                    list(points = child.locs[trimmed, , drop = FALSE],
                         parents = parent.locs,
                         parent.ids = parent.ids[trimmed],
                         child.ys = sim.n.children$child.ys[trimmed],
                         sibling.list = self$trim.siblings(sibling.list, trimmed))
                },
                ## A method to trim the sibling list.
                trim.siblings = function(sibling.list, trimmed){
                    sibling.list
                },
                ## A method to get parents by simulation.
                get.parents = function(pars){
                    expected.parents <- pars["D"]*self$vol
                    n.parents <- rpois(1, expected.parents)
                    parent.locs <- matrix(0, nrow = n.parents, ncol = self$dim)
                    for (i in 1:self$dim){
                        parent.locs[, i] <- runif(n.parents, self$lims[i, 1], self$lims[i, 2])
                    }
                    parent.locs
                },
                ## Overwriting method for the Palm intensity.
                palm.intensity = function(r, pars){
                    self$sibling.pi(r, pars) + self$nonsibling.pi(pars)
                },
                ## An empty method for the expectation of the child distribution.
                child.expectation = function(pars){},
                ## An empty method for the variance of the child distribution.
                child.variance = function(pars){},
                ## A method for the expected number of siblings of a randomly chosen point.
                sibling.expectation = function(pars){
                    (self$child.variance(pars) + self$child.expectation(pars)^2)/self$child.expectation(pars) - 1
                },
                ## An empty method for the PDF of Q, the between-sibling distances.
                fq = function(r, pars){},
                ## An empty method for the CDF of Q.
                Fq = function(r, pars){},
                ## A default method for the quotient of the PDF of Q and the surface volume.
                ## Numerically unstable; better to replace in child classes.
                fq.over.s = function(r, pars){
                    self$fq(r, pars)/Sd(r, dim)
                },
                ## A method for the Palm intensity of nonsibling points.
                nonsibling.pi = function(pars){
                    pars["D"]*self$child.expectation(pars)
                },
                ## A method for the PI of sibling points.
                sibling.pi = function(r, pars){
                    self$sibling.expectation(pars)*self$fq.over.s(r, pars)
                }
            ))
}

######
## Class for known sibling relationships.
######

set.sibling.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("sibling.inherit", class, envir = class.env)
    R6Class("palm_sibling",
            inherit = class.env$sibling.inherit,
            public = list(
                siblings = NULL,
                ## Lol one day you'll get a bug because you're an
                ## idiot and you have objects called 'sibling.list'
                ## and 'siblings.list' and you'll be so mad.
                siblings.list = NULL,
                sibling.list = NULL,
                sibling.mat = NULL,
                sibling.alpha = NULL,
                sibling.beta = NULL,
                ## Initialisation method.
                initialize = function(sibling.list, ...){
                    self$sibling.list <- sibling.list
                    super$initialize(...)
                    self$siblings.list <- list()
                    for (i in 1:self$n.patterns){
                        self$get.siblings(i)
                    }
                },
                ## Overwriting the method to set up a new pattern.
                setup.pattern = function(pattern, do.contrasts = TRUE, do.siblings = TRUE){
                    super$setup.pattern(pattern, do.contrasts)
                    self$sibling.mat <- self$sibling.list[[pattern]]$sibling.mat
                    self$sibling.alpha <- self$sibling.list[[pattern]]$alpha
                    self$sibling.beta <- self$sibling.list[[pattern]]$beta
                    if (do.siblings){
                        self$siblings <- self$siblings.list[[pattern]]
                    }
                },
                ## A method to get vector of sibling relationships
                ## that matches with the contrasts.
                get.siblings = function(pattern){
                    self$setup.pattern(pattern, do.siblings = FALSE)
                    siblings <- numeric(self$n.contrasts)
                    for (i in 1:self$n.contrasts){
                        pair <- self$contrast.pairs[i, ]
                        siblings[i] <- self$sibling.mat[pair[1], pair[2]]
                    }
                    ## 0 for known nonsiblings.
                    ## 1 for known siblings.
                    ## 2 for unknown sibling status.
                    siblings <- as.numeric(siblings)
                    siblings[is.na(siblings)] <- 2
                    self$siblings.list[[pattern]] <- siblings
                },
                ## Overwriting the sum of the log intensities.
                sum.log.intensities = function(pars){
                    all.palm.intensities <- numeric(self$n.contrasts)
                    all.palm.intensities[self$siblings == 0] <-
                        self$sibling.beta*self$nonsibling.pi(pars)
                    all.palm.intensities[self$siblings == 1] <-
                        self$sibling.alpha*
                        self$sibling.pi(self$contrasts[self$siblings == 1], pars)
                    all.palm.intensities[self$siblings == 2] <-
                        (1 - self$sibling.beta)*self$nonsibling.pi(pars) +
                        (1 - self$sibling.alpha)*
                        self$sibling.pi(self$contrasts[self$siblings == 2], pars)
                    sum(log(self$n.points*all.palm.intensities))
                }
            ))
}

######
## Class for Poisson number of children.
######
set.poischild.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("poischild.inherit", class, envir = class.env)
    R6Class("palm_poischild",
            inherit = class.env$poischild.inherit,
            public = list(
                ## Adding lambda parameter.
                fetch.pars = function(){
                    super$fetch.pars()
                    self$add.pars(name = "lambda", link = log, start = sqrt(self$n.points/self$vol),
                                  lower = 0, upper = self$n.points)
                },
                ## Simulation method for the number of children per parent.
                simulate.n.children = function(n, pars){
                    list(n.children = rpois(n, pars["lambda"]))
                },
                ## A method for the expectation of the child distribution.
                child.expectation = function(pars){
                    pars["lambda"]
                },
                ## A method for the variance of the child distribution.
                child.variance = function(pars){
                    pars["lambda"]
                }
            ))
}

######
## Class for binomial number of children.
######
set.binomchild.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("binomchild.inherit", class, envir = class.env)
    R6Class("palm_binomchild",
            inherit = class.env$binomchild.inherit,
            public = list(
                binom.size = NULL,
                initialize = function(child.list, ...){
                    self$child.list <- child.list
                    self$binom.size <- child.list$size
                    super$initialize(...)
                },
                ## Adding p parameter.
                fetch.pars = function(){
                    super$fetch.pars()
                    p.start <- sqrt(self$n.points/self$vol)/self$binom.size
                    p.start <- ifelse(p.start > 1, 0.9, p.start)
                    p.start <- ifelse(p.start < 0, 0.1, p.start)
                    self$add.pars(name = "p", link = logit, start = p.start, lower = 0, upper = 1)
                },
                ## Simulation method for the number of children per parent.
                simulate.n.children = function(n, pars){
                    list(n.children = rbinom(n, self$binom.size, pars["p"]))
                },
                ## A method for the expectation of the child distribution.
                child.expectation = function(pars){
                    self$binom.size*pars["p"]
                },
                ## A method for the variance of the child distribution.
                child.variance = function(pars){
                    self$binom.size*pars["p"]*(1 - pars["p"])
                }
            ))
}

######
## Class for two-camera children distribution.
######
set.twocamerachild.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("twocamerachild.inherit", class, envir = class.env)
    R6Class("palm_twocamerachild",
            inherit = class.env$twocamerachild.inherit,
            public = list(
                ## Detection zone halfwidth.
                twocamera.w = NULL,
                ## Survey area halfwidth.
                twocamera.b = NULL,
                ## Lag between cameras.
                twocamera.l = NULL,
                ## Mean dive-cycle duration.
                twocamera.tau = NULL,
                initialize = function(child.list, ...){
                    self$child.list <- child.list
                    self$twocamera.w <- child.list$twocamera.w
                    self$twocamera.b <- child.list$twocamera.b
                    self$twocamera.l <- child.list$twocamera.l
                    self$twocamera.tau <- child.list$twocamera.tau
                    super$initialize(...)
                    self$setup.pattern(1)
                    if (self$dim != 1){
                        stop("Analysis of two-camera surveys is only implemented for one-dimensional processes.")
                    }
                },
                ## Adding kappa parameter.
                fetch.pars = function(){
                    super$fetch.pars()
                    self$add.pars(name = "kappa", link = log, start = 0.1*self$twocamera.tau, lower = 0,
                                  upper = self$twocamera.tau)
                },
                ## Overwriting base simulation method in case D.2D is provided.
                simulate.pattern = function(pars = self$par.fitted){
                    if (any(names(pars) == "D.2D")){
                        which.D <- which(names(pars) == "D.2D")
                        pars[which.D] <- 2*pars[which.D]*self$twocamera.b
                        names(pars)[which.D] <- "D"
                    }
                    super$simulate.pattern(pars)
                },
                ## Simulation method for the number of children per parent.
                simulate.n.children = function(n, pars){
                    probs <- twocamera.probs(self$twocamera.l, self$twocamera.tau, self$twocamera.w,
                                             self$twocamera.b, pars["kappa"], pars["sigma"])
                    ## Generating some y values.
                    parent.ys <- runif(n, -self$twocamera.b, self$twocamera.b)
                    child.ys <- matrix(nrow = n, ncol = 2)
                    for (i in 1:n){
                        child.ys[i, ] <- rnorm(2, parent.ys[i], pars["sigma"])
                    }
                    camera1.in <- ifelse(child.ys[, 1] < self$twocamera.w & child.ys[, 1] > -self$twocamera.w,
                                         TRUE, FALSE)
                    camera2.in <- ifelse(child.ys[, 2] < self$twocamera.w & child.ys[, 2] > -self$twocamera.w,
                                         TRUE, FALSE)
                    camera1.up <- sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(probs$p.up, probs$p.down))
                    camera2.up <- logical(n)
                    for (i in 1:n){
                        p.up <- ifelse(camera1.up[i], 1 - probs$p.down.up, probs$p.up.down)
                        camera2.up[i] <- sample(c(TRUE, FALSE), 1, prob = c(p.up, 1 - p.up))
                    }
                    camera1.det <- camera1.in & camera1.up
                    camera2.det <- camera2.in & camera2.up
                    child.ys <- t(child.ys)[t(cbind(camera1.det, camera2.det))]
                    n.children <-  camera1.det + camera2.det
                    cameras <- numeric(sum(n.children))
                    j <- 1
                    for (i in 1:n){
                        if (n.children[i] > 0){
                            cameras[j:(j + n.children[i] - 1)] <- c(1[camera1.det[i]], 2[camera2.det[i]])
                        }
                        j <- j + n.children[i]
                    }
                    sibling.list <- siblings.twocamera(cameras)
                    sibling.list$cameras <- cameras
                    list(n.children = n.children, sibling.list = sibling.list, child.ys = child.ys)
                },
                ## Overwriting the method to trim the sibling list for children outside the window.
                trim.siblings = function(sibling.list, trimmed){
                    sibling.list$sibling.mat <- sibling.list$sibling.mat[trimmed, trimmed, drop = FALSE]
                    sibling.list$cameras <- sibling.list$cameras[trimmed]
                    super$trim.siblings(sibling.list, trimmed)
                },
                ## A method for the expectation of the child distribution.
                child.expectation = function(pars){
                    probs <- twocamera.probs(self$twocamera.l, self$twocamera.tau, self$twocamera.w,
                                            self$twocamera.b, pars["kappa"], pars["sigma"])
                    2*probs$p.10/(probs$p.10 + probs$p.01)
                },
                ## A method for the variance of the child distribution.
                child.variance = function(pars){
                    probs <- twocamera.probs(self$twocamera.l, self$twocamera.tau, self$twocamera.w,
                                            self$twocamera.b, pars["kappa"], pars["sigma"])
                    2*probs$p.10*probs$p.01*(2 - probs$p.10 - probs$p.01)/(probs$p.10 + probs$p.01)^2
                }
            ))
}

######
## Class for Thomas processes.
######

set.thomas.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("thomas.inherit", class, envir = class.env)
    R6Class("palm_thomas",
            inherit = class.env$thomas.inherit,
            public = list(
                ## Adding sigma paremter.
                fetch.pars = function(){
                    super$fetch.pars()
                    self$add.pars(name = "sigma", link = log, start = 0.1*self$R, lower = 0,
                                  upper = self$R)
                },
                ## Simulation method for children locations.
                simulate.children = function(n, parent.loc, pars){
                    rmvnorm(n, parent.loc, pars["sigma"]^2*diag(self$dim))
                },
                ## Overwriting method for the integral in the Palm likelihood.
                palm.likelihood.integral = function(pars){
                    -self$pi.multiplier*(pars["D"]*self$child.expectation(pars)*Vd(self$R, self$dim) +
                                    self$sibling.expectation(pars)*self$Fq(self$R, pars))
                },
                ## Overwriting method for the PDF of Q.
                fq = function(r, pars){
                    2^(1 - self$dim/2)*r^(self$dim - 1)*exp(-r^2/(4*pars["sigma"]^2))/
                        ((pars["sigma"]*sqrt(2))^self$dim*gamma(self$dim/2))
                },
                ## Overwriting method for the CDF of Q.
                Fq = function(r, pars){
                    pgamma(r^2/(4*pars["sigma"]^2), self$dim/2)
                },
                ## Overwriting method for the quotient of the PDF of Q and the surface volume.
                fq.over.s = function(r, pars){
                    exp(-r^2/(4*pars["sigma"]^2))/((2*pars["sigma"])^self$dim*pi^(self$dim/2))
                },
                ## A method for xlim.
                get.xlim = function(){
                    c(0, min(self$R, 10*self$par.fitted["sigma"]))
                }
            ))
}

######
## Class for Matern processes.
######

set.matern.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("matern.inherit", class, envir = class.env)
    R6Class("palm_matern",
            inherit = class.env$matern.inherit,
            public = list(
                ## Adding tau parameter.
                fetch.pars = function(){
                    super$fetch.pars()
                    self$add.pars(name = "tau", link = log, start = 0.1*self$R, lower = 0, upper = self$R)
                },
                ## Simulation method for children locations via rejection.
                simulate.children = function(n, parent.loc, pars){
                    out <- matrix(0, nrow = n, ncol = self$dim)
                    for (i in 1:n){
                        keep <- FALSE
                        while (!keep){
                            proposal <- runif(self$dim, -pars["tau"], pars["tau"])
                            if (sqrt(sum(proposal^2)) <= pars["tau"]){
                                proposal <- proposal + parent.loc
                                out[i, ] <- proposal
                                keep <- TRUE
                            }
                        }
                    }
                    out
                },
                ## Overwriting method for the PDF of Q.
                fq = function(r, pars){
                    ifelse(r > 2*pars["tau"], 0,
                           2*self$dim*r^(self$dim - 1)*(pars["tau"]*hyperg_2F1(0.5, 0.5 - self$dim/2, 1.5, 1) -
                                                        r/2*hyperg_2F1(0.5, 0.5 - self$dim/2, 1.5, r^2/(4*pars["tau"]^2)))/
                           (beta(self$dim/2 + 0.5, 0.5)*pars["tau"]^(self$dim + 1)))
                },
                ## Overwriting method for the CDF of Q.
                Fq = function(r, pars){
                    alpha <- r^2/(4*pars["tau"]^2)
                    r^2/pars["tau"]^2*(1 - pbeta(alpha, 0.5, self$dim/2 + 0.5)) +
                        2^self$dim*incomplete.beta(alpha, self$dim/2 + 0.5, self$dim/2 + 0.5)/
                            beta(0.5, self$dim/2 + 0.5)
                },
                ## Overwriting method for the quotient of the PDF of Q and the surface volume.
                fq.over.s = function(r, pars){
                    ifelse(r > 2*pars["tau"], 0,
                           2*(pars["tau"]*hyperg_2F1(0.5, 0.5 - self$dim/2, 1.5, 1) -
                              r/2*hyperg_2F1(0.5, 0.5 - self$dim/2, 1.5, r^2/(4*pars["tau"]^2)))*gamma(self$dim/2 + 1)/
                           (beta(self$dim/2 + 0.5, 0.5)*pars["tau"]^(self$dim + 1)*pi^(self$dim/2)))
                },
                ## A method for xlim.
                get.xlim = function(){
                    c(0, min(self$R, 10*self$par.fitted["tau"]))
                }
            ))
}

######
## Class for void processes.
######

set.void.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("void.inherit", class, envir = class.env)
    R6Class("palm_void",
            inherit = class.env$void.inherit,
            public = list(
                ## Adding children and parent density parameters.
                fetch.pars = function(){
                    super$fetch.pars()
                    self$add.pars(name = "Dp", link = log, start = 10/self$vol, lower = 0, upper = Inf)
                    self$add.pars(name = "Dc", link = log, start = self$n.points/self$vol, lower = 0, upper = Inf)
                },
                ## Overwriting simulation method.
                simulate.pattern = function(pars = self$par.fitted){
                    ## Generating children.
                    expected.children <- pars["Dc"]*self$vol
                    n.children <- rpois(1, expected.children)
                    child.locs <- matrix(0, nrow = n.children, ncol = self$dim)
                    for (i in 1:self$dim){
                        child.locs[, i] <- runif(n.children, self$lims[i, 1], self$lims[i, 2])
                    }
                    ## Generating parents.
                    parent.locs <- self$get.parents(pars)
                    list(points = self$delete.points(child.locs, parent.locs, pars), parents = parent.locs)
                },
                ## A method to get parents by simulation.
                get.parents = function(pars){
                    expected.parents <- pars["Dp"]*self$vol
                    n.parents <- rpois(1, expected.parents)
                    parent.locs <- matrix(0, nrow = n.parents, ncol = self$dim)
                    for (i in 1:self$dim){
                        parent.locs[, i] <- runif(n.parents, self$lims[i, 1], self$lims[i, 2])
                    }
                    parent.locs
                },
                ## Overwriting method for the Palm intensity.
                palm.intensity = function(r, pars){
                    pars["Dc"]*self$prob.safe(r, pars)
                }
            ))
}

set.totaldeletion.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("totaldeletion.inherit", class, envir = class.env)
    R6Class("palm_totaldeletion",
            inherit = class.env$totaldeletion.inherit,
            public = list(
                fetch.pars = function(){
                    super$fetch.pars()
                    self$add.pars(name = "tau", link = log, start = 0.1*self$R, lower = 0, upper = self$R)
                },
                ## Probability of being safe, given distance r from a known point.
                prob.safe = function(r, pars){
                    exp(-pars["Dp"]*Vd(pars["tau"], self$dim)*(1 - pbeta(1 - (r/(2*pars["tau"])), (self$dim + 1)/2, 0.5)))
                },
                ## Function to delete children, given children and parent locations.
                delete.points = function(child.locs, parent.locs, pars){
                    dists <- euc_distances(child.locs[, 1], child.locs[, 2], parent.locs[, 1], parent.locs[, 2])
                    child.locs[apply(dists, 1, min) > pars["tau"], ]
                },
                ## A method for xlim.
                get.xlim = function(){
                    c(0, min(self$R, 10*self$par.fitted["tau"]))
                }
            ))
}

######
## Class for simulating given known parent locations.
######

set.giveparent.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("giveparent.inherit", class, envir = class.env)
    R6Class("palm_giveparent",
            inherit = class.env$giveparent.inherit,
            public = list(
                parent.locs = NULL,
                initialize = function(parent.locs, ...){
                    self$parent.locs <- parent.locs
                    super$initialize(...)
                },
                ## Overwriting method for getting parents.
                get.parents = function(pars){
                    if (!(ncol(self$parent.locs) == self$dim)){
                        stop("Incorrection dimensions of parent locations.")
                    }
                    self$parent.locs
                }
            ))
}


## Function to create R6 object with correct class hierarchy.
create.obj <- function(classes, points, lims, R, child.list, parent.locs, sibling.list, trace, start, bounds){
    class <- base.class.R6
    n.classes <- length(classes)
    class.env <- new.env()
    for (i in 1:n.classes){
        set.class <- get(paste("set", classes[i], "class", sep = "."))
        class <- set.class(class, class.env)
    }
    if (any(classes == "twocamerachild") & !any(classes == "thomas")){
        stop("Analysis of two-camera surveys is only implemented for Thomas processes.")
    }
    if (!is.null(points) & !(is.matrix(points) | is.list(points))){
        stop("The argument 'points' must be a matrix, or a list of matrices.")
    }
    if (is.matrix(points)){
        points <- list(points)
        sibling.list <- list(sibling.list)
    }
    if (is.matrix(lims)){
        lims <- list(lims)
    }
    class$new(points.list = points, lims.list = lims, R = R,
              child.list = child.list, parent.locs = parent.locs,
              sibling.list = sibling.list, trace = trace,
              classes = classes, start = start, bounds = bounds)
}

## Some objects to get around R CMD check.
super <- NULL
self <- NULL


