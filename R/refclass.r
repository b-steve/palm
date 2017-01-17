######
## General base class.
######

base.class.R6 <- R6Class("nspp_r6",
                      public = list(
                          ## Setting fields.
                          points = NULL,
                          n.points = NULL,
                          contrasts = NULL,
                          lims = NULL,
                          vol = NULL,
                          dim = NULL,
                          R = NULL,
                          child.dist = NULL,
                          boots = NULL,
                          n.par = NULL,
                          set.start = NULL,
                          par.names = NULL,
                          par.link.names = NULL,
                          par.start = NULL,
                          par.start.link = NULL,
                          par.links = NULL,
                          par.invlinks = NULL,
                          par.fitted = NULL,
                          par.fitted.link = NULL,
                          par.lower = NULL,
                          par.lower.link = NULL,
                          par.upper = NULL,
                          par.upper.link = NULL,
                          gr = NULL,
                          trace = NULL,
                          conv.code = NULL,
                          classes = NULL,
                          ## Initialisation method.
                          initialize = function(points, lims, R, trace, classes, start = NULL){
                              self$points <- points
                              self$n.points <- nrow(points)
                              self$lims <- lims
                              self$vol <- prod(apply(self$lims, 1, diff))
                              self$dim <- ncol(points)
                              self$R <- R
                              self$set.start <- start
                              self$get.contrasts()
                              self$get.pars()
                              self$get.link.names()
                              self$n.par <- length(self$par.start)
                              self$get.invlinks()
                              self$par.start.link <- self$link.pars(self$par.start)
                              self$get.link.bounds()
                              self$trace <- trace
                              self$classes <- classes
                          },
                          ## An empty method for getting contrasts.
                          get.contrasts = function(){
                          },
                          ## A method for getting the parameters across all classes.
                          get.pars = function(){
                              self$fetch.pars()
                          },
                          ## An empty method for fetching parameters from a particular class.
                          fetch.pars = function(){
                          },
                          ## A method for adding new parameters.
                          add.pars = function(names, links, start, lower, upper){
                              self$par.names <- c(self$par.names, names)
                              self$par.link.names <- self$par.names
                              self$par.links <- c(self$par.links, links)
                              names(self$par.links) <- self$par.names
                              for (i in 1:length(names)){
                                  if (any(names == names(self$set.start))){
                                      start[i] <- self$set.start[names[i]]
                                  }
                              }
                              self$par.start <- c(self$par.start, start)
                              names(self$par.start) <- self$par.names
                              self$par.upper <- c(self$par.upper, upper)
                              names(self$par.upper) <- self$par.names
                              self$par.lower <- c(self$par.lower, lower)
                              names(self$par.lower) <- self$par.names
                          },
                          ## A default method for getting the names of the link parameters.
                          get.link.names = function(){
                              self$par.link.names <- self$par.names
                          },
                          ## A method for setting inverse links.
                          get.invlinks = function(){
                              par.invlinks <- vector(mode = "list", length = self$n.par)
                              par.lower <-  numeric(self$n.par)
                              par.upper <-  numeric(self$n.par)
                              names(par.invlinks) <- self$par.names
                              for (i in 1:length(self$par.links)){
                                  if (identical(self$par.links[[i]], identity)){
                                      par.invlinks[[i]] <- identity
                                  } else if (identical(self$par.links[[i]], log)){
                                      par.invlinks[[i]] <- exp
                                  } else {
                                      stop("Link functions must be either identity or log.")
                                  }
                              }
                              self$par.invlinks <- par.invlinks
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
                                  out[i] <- self$par.links[[i]](pars[i])
                              }
                              names(out) <- self$par.link.names
                              out
                          },
                          ## A method for converting parameters to their real scale.
                          invlink.pars = function(pars){
                              out <- numeric(self$n.par)
                              for (i in 1:self$n.par){
                                  out[i] <- self$par.invlinks[[i]](pars[i])
                              }
                              names(out) <- self$par.names
                              out
                          },
                          ## A method for the empirical Palm intensity function.
                          empirical.palm = function(){
                              
                          },
                          ## An empty method for simulation.
                          simulate = function(pars){},
                          ## A method to trim points to the observation window.
                          trim.points = function(points){
                              in.window <- rep(TRUE, nrow(points))
                              for (i in 1:self$dim){
                                  in.window <- in.window & (self$lims[i, 1] <= points[, i] & self$lims[i, 2] >= points[, i])
                              }
                              points[in.window, , drop = FALSE]
                          },
                          ## An empty method for the Palm intensity.
                          palm.intensity = function(r, pars){},
                          ## A default method for the sum of the log Palm intensities.
                          sum.log.intensities = function(pars){
                              sum(log(self$n.points*self$palm.intensity(self$contrasts, pars)))
                          },
                          ## A default method for the integral in the Palm likelihood.
                          palm.likelihood.integral = function(pars){
                              f <- function(r, pars){
                                  self$palm.intensity(r, pars)*Sd(r, self$dim)
                              }
                              -self$n.points*integrate(f, lower = 0, upper = self$R, pars = pars)$value
                          },
                          ## A default method for the log of the Palm likelihood function.
                          log.palm.likelihood = function(pars){
                              out <- self$sum.log.intensities(pars) + self$palm.likelihood.integral(pars)
                              if (self$trace){
                                  for (i in 1:self$n.par){
                                      cat(self$par.names[i], ": ", sep = "")
                                      cat(pars[i], ", ", sep = "")
                                  }
                                  cat("Log-lik: ", out, "\n", sep = "")
                              }
                              out
                          },
                          ## A method for the negative Palm likelihood function.
                          neg.log.palm.likelihood = function(pars){
                              -self$log.palm.likelihood(pars)
                          },
                          ## A method for the negative Palm likelihood function with linked parameters.
                          link.neg.log.palm.likelihood = function(link.pars){
                              names(link.pars) <- names(self$par.start.link)
                              pars <- self$invlink.pars(link.pars)
                              self$neg.log.palm.likelihood(pars)
                          },
                          ## A method for model fitting.
                          fit = function(){
                              optim.obj <- nlminb(self$par.start.link, self$link.neg.log.palm.likelihood,
                                                  control = list(eval.max = 2000, iter.max = 1500),
                                                  lower = self$par.lower.link, upper = self$par.upper.link)
                              #optim.obj <- optim(self$par.start.link, self$link.neg.log.palm.likelihood,
                              #                   lower = self$par.lower.link, upper = self$par.upper.link)
                              #optim.obj <- nlm(self$link.neg.log.palm.likelihood, self$par.start.link)
                              if (optim.obj$convergence != 0){
                                  warning("Failed convergence.")
                              }
                              self$par.fitted.link <- optim.obj$par
                              names(self$par.fitted.link) <- self$par.link.names
                              self$par.fitted <- self$invlink.pars(self$par.fitted.link)
                              self$conv.code <- optim.obj$convergence
                          },
                          ## A method for bootstrapping.
                          boot = function(N, prog = TRUE){
                              boots <- matrix(0, nrow = N, ncol = self$n.par)
                              for (i in 1:N){
                                  points.boot <- self$simulate()
                                  obj.boot <- create.obj(self$classes, points.boot, self$lims,
                                                         self$R, FALSE, self$par.fitted)
                                  obj.boot$fit()
                                  if (obj.boot$conv.code != 0){
                                      boots[i, ] <- rep(NA, self$n.par)
                                  } else {
                                      boots[i, ] <- coef(obj.boot)
                                  }
                              }
                              self$boots <- boots
                          }
                      ))

######
## Template for new classes.
######

## Replace CLASSNAME with class name, then add fields and methods.
set.CLASSNAME.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("CLASSNAME.inherit", class, envir = class.env)
    R6Class("nspp_r6",
            inherit = class.env$CLASSNAME.inherit,
            public = list(

            ))
}

######
## Class for periodic boundary conditions.
######

set.pbc.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("pbc.inherit", class, envir = class.env)
    R6Class("pbc",
            inherit = class.env$pbc.inherit,
            public = list(
                ## A method to generate contrasts.
                get.contrasts = function(){
                    contrasts <- pbc_distances(points = self$points, lims = self$lims)
                    self$contrasts <- contrasts[contrasts <= self$R]
                }
            ))
}

######
## Class for Neyman-Scott processes.
######

set.ns.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("ns.inherit", class, envir = class.env)
    R6Class("ns",
            inherit = class.env$ns.inherit,
            public = list(
                par.name.child = NULL,
                initialize = function(points, lims, R, trace, classes, start = NULL){
                    self$get.par.name.child()
                    super$initialize(points, lims, R, trace, classes, start)
                    self$par.link.names[self$par.link.names %in%
                                        c("D", self$par.name.child)] <- c("Dc", "nu")
                },
                ## An empty method for getting the name of the child parameter.
                get.par.name.child = function(){},
                ## Adding density parameter.
                fetch.pars = function(){
                    super$fetch.pars()
                    self$add.pars("D", log, 1/(4*pi*(0.1*self$R)^2), 0, Inf)
                },
                ## Overwriting simulation method.
                simulate = function(pars = self$par.fitted){
                    expected.parents <- pars["D"]*self$vol
                    n.parents <- rpois(1, expected.parents)
                    parent.locs <- matrix(0, nrow = n.parents, ncol = self$dim)
                    for (i in 1:self$dim){
                        parent.locs[, i] <- runif(n.parents, self$lims[i, 1], self$lims[i, 2])
                    }
                    n.children <- self$simulate.n.children(n.parents, pars)
                    child.locs <- matrix(0, nrow = sum(n.children), ncol = self$dim)
                    j <- 0
                    for (i in 1:n.parents){
                        if (n.children[i] > 0){
                            child.locs[(j + 1):(j + n.children[i]), ] <-
                                self$simulate.children(n.children[i], parent.locs[i, ], pars)
                            j <- j + n.children[i]
                        }
                    }
                    self$trim.points(child.locs)
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
## Class for Poisson number of children.
######
set.poischild.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("poischild.inherit", class, envir = class.env)
    R6Class("nspp_r6",
            inherit = class.env$poischild.inherit,
            public = list(
                ## A method for getting the parameter name.
                get.par.name.child = function(){
                    self$par.name.child <- "lambda"
                },
                ## Adding lambda parameter.
                fetch.pars = function(){
                    super$fetch.pars()
                    self$add.pars("lambda", log, 4*pi*(0.1*self$R)^2*self$n.points/self$vol, 0, self$n.points)
                },
                ## Simulation method for the number of children per parent.
                simulate.n.children = function(n, pars){
                    rpois(n, pars["lambda"])
                },
                ## An method for the expectation of the child distribution.
                child.expectation = function(pars){
                    pars["lambda"]
                },
                ## An method for the variance of the child distribution.
                child.variance = function(pars){
                    pars["lambda"]
                }
            ))
}

######
## Class for Thomas processes.
######

set.thomas.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("thomas.inherit", class, envir = class.env)
    R6Class("nspp_r6",
            inherit = class.env$thomas.inherit,
            public = list(
                ## Adding sigma paremter.
                fetch.pars = function(){
                    super$fetch.pars()
                    self$add.pars("sigma", log, 0.1*self$R, 0, self$R)
                },
                ## Simulation method for children locations.
                simulate.children = function(n, parent.loc, pars){
                    rmvnorm(n, parent.loc, pars["sigma"]^2*diag(self$dim))
                },
                ## Overwriting method for the integral in the Palm likelihood.
                palm.likelihood.integral = function(pars){
                    -self$n.points*(pars["D"]*self$child.expectation(pars)*Vd(self$R, self$dim) +
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
                }
            ))
}

######
## Class for Matern processes.
######

set.matern.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("matern.inherit", class, envir = class.env)
    R6Class("nspp_r6",
            inherit = class.env$matern.inherit,
            public = list(
                ## Adding tau parameter.
                fetch.pars = function(){
                    super$fetch.pars()
                    self$add.pars("tau", log, 0.1*self$R, 0, self$R)
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
                }
            ))
}

######
## Class for void processes.
######

set.void.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("void.inherit", class, envir = class.env)
    R6Class("nspp_r6",
            inherit = class.env$void.inherit,
            public = list(

            ))
}


## Function to create R6 object with correct class hierarchy.
create.obj <- function(classes, points, lims, R, trace, start){
    class <- base.class.R6
    n.classes <- length(classes)
    class.env <- new.env()
    for (i in 1:n.classes){
        set.class <- get(paste("set", classes[i], "class", sep = "."))
        class <- set.class(class, class.env)
    }
    class$new(points, lims, 0.5, trace, classes, start)
}

## Some objects to get around R CMD check.
super <- NULL
self <- NULL

