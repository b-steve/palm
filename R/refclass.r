######
## General base class.
######

base.class.R6 <- R6Class("nspp-r6",
                      public = list(
                          ## Setting fields.
                          points = NULL,
                          n.points = NULL,
                          contrasts = NULL,
                          lims = NULL,
                          dim = NULL,
                          R = NULL,
                          child.dist = NULL,
                          par.names = NULL,
                          par.start = NULL,
                          par.start.link = NULL,
                          par.links = NULL,
                          par.invlinks = NULL,
                          par.fitted = NULL,
                          par.fitted.link = NULL,
                          ## Initialisation method.
                          initialize = function(points, lims, R){
                              self$points <- points
                              self$n.points <- nrow(points)
                              self$lims <- lims
                              self$dim <- ncol(points)
                              self$R <- R
                              self$get.contrasts()
                              self$get.pars()
                              self$get.invlinks()
                              self$par.start.link <- self$link.pars(self$par.start)
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
                          add.pars = function(names, links, start){
                              self$par.names <- c(self$par.names, names)
                              self$par.links <- c(self$par.links, links)
                              names(self$par.links) <- self$par.names
                              self$par.start <- c(self$par.start, start)
                              names(self$par.start) <- self$par.names
                          },
                          ## A method for setting inverse links.
                          get.invlinks = function(){
                              par.invlinks <- vector(mode = "list", length = length(self$par.links))
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
                          ## A method for converting parameters to their link scales.
                          link.pars = function(pars){
                              n.pars <- length(pars)
                              out <- numeric(n.pars)
                              for (i in 1:n.pars){
                                  out[i] <- self$par.links[[i]](pars[i])
                              }
                              names(out) <- names(pars)
                              out
                          },
                          ## A method for converting parameters to their real scale.
                          invlink.pars = function(pars){
                              n.pars <- length(pars)
                              out <- numeric(n.pars)
                              for (i in 1:n.pars){
                                  out[i] <- self$par.invlinks[[i]](pars[i])
                              }
                              names(out) <- names(pars)
                              out
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
                              self$sum.log.intensities(pars) + self$palm.likelihood.integral(pars)
                          },
                          ## A method for the negative Palm likelihood function.
                          neg.log.palm.likelihood = function(pars){
                              -self$log.palm.likelihood(pars)
                          },
                          ## A method for the negative Palm likelihood function with linked parameters.
                          link.neg.log.palm.likelihood = function(link.pars){
                              pars <- self$invlink.pars(link.pars)
                              self$neg.log.palm.likelihood(pars)
                          },
                          ## A method for model fitting.
                          fit = function(){
                              self$par.fitted.link <- optim(self$par.start.link, self$link.neg.log.palm.likelihood)$par
                              self$par.fitted <- self$invlink.pars(self$par.fitted.link)
                          }
                      ))

######
## Template for new classes.
######

## Replace CLASSNAME with class name, then add fields and methods.
set.CLASSNAME.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("CLASSNAME.inherit", class, envir = class.env)
    R6Class("nspp-r6",
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
                ## Adding density parameter.
                fetch.pars = function(){
                    super$fetch.pars()
                    vol <- prod(apply(self$lims, 1, diff))
                    self$add.pars("D", log, sqrt(self$n.points/vol))
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
    R6Class("nspp-r6",
            inherit = class.env$poischild.inherit,
            public = list(
                ## Adding lambda parameter.
                fetch.pars = function(){
                    super$fetch.pars()
                    vol <- prod(apply(self$lims, 1, diff))
                    self$add.pars("lambda", log, sqrt(self$n.points/vol))
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
    R6Class("nspp-r6",
            inherit = class.env$thomas.inherit,
            public = list(
                ## Adding sigma paremter.
                fetch.pars = function(){
                    super$fetch.pars()
                    self$add.pars("sigma", log, 0.1*self$R)
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

## Replace CLASSNAME with class name, then add fields and methods.
set.matern.class <- function(class, class.env){
    ## Saving inherited class to class.env.
    assign("matern.inherit", class, envir = class.env)
    R6Class("nspp-r6",
            inherit = class.env$matern.inherit,
            public = list(
                ## Adding tau parameter.
                fetch.pars = function(){
                    super$fetch.pars()
                    self$add.pars("tau", log, 0.1*self$R)
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

## Function to create R6 object with correct class hierarchy.
create.obj <- function(classes, points, lims, R){
    class <- base.class.R6
    n.classes <- length(classes)
    class.env <- new.env()
    for (i in 1:n.classes){
        set.class <- get(paste("set", classes[i], "class", sep = "."))
        class <- set.class(class, class.env)
    }
    class$new(points, lims, 0.5)
}

## Some objects to get around R CMD check.
super <- NULL
self <- NULL
