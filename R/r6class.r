######
## General base class.
######

base.class.R6 <- R6Class("base",
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
                              self$contrasts <- self$get.contrasts()
                              #self$par.start <- self$get.start()
                              #self$par.names <- names(self$par.start)
                              #self$par.links <- self$get.links()
                              #self$par.invlinks <- self$get.invlinks()
                              #self$par.start.link <- self$link.pars(self$par.start)
                          },
                          ## An empty method for getting contrasts.
                          get.contrasts = function(){
                          },
                          ## An empty method for setting start values.
                          get.start = function(){},
                          ## An empty method for getting link functions.
                          get.links = function(){},
                          ## A method for setting inverse links.
                          ## get.invlinks = function(){
                          ##     par.invlinks.save <- vector(mode = "list", length = length(self$par.links))
                          ##     names(par.invlinks.save) <- self$par.names
                          ##     for (i in 1:length(self$par.links)){
                          ##         if (identical(self$par.links[[i]], identity)){
                          ##             par.invlinks.save[[i]] <- identity
                          ##         } else if (identical(self$par.links[[i]], log)){
                          ##             par.invlinks.save[[i]] <- exp
                          ##         } else {
                          ##             stop("Link functions must be either identity or log.")
                          ##         }
                          ##     }
                          ##     par.invlinks.save
                          ## },
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
                              neg.log.parlm.likelihood(pars)
                          },
                          ## A method for model fitting.
                          fit = function(){
                              self$par.fitted.link <- optim(self$par.start.link, self$link.neg.log.palm.likelihood)$par
                              self$par.fitted <- self$invlink.pars(self$par.fitted.link)
                          }
                      ))

######
## Class for periodic boundary conditions.
######

use.pbc.class <- function(class, class.env){
    ## Saving inherited class to environment of creat.obj().
    assign("pbc.inherit", class, envir = class.env)
    R6Class("nspp-r6",
            inherit = class.env$pbc.inherit,
            public = list(
                ## A method to generate contrasts.
                get.contrasts = function(){
                    contrasts.save <- pbc_distances(points = self$points, lims = self$lims)
                    contrasts.save[contrasts.save <= self$R]
                }
            ))
}

create.obj <- function(classes, points, lims, R){
    class <- base.class.R6
    n.classes <- length(classes)
    class.env <- new.env()
    for (i in 1:n.classes){
        use.class <- get(paste("use", classes[i], "class", sep = "."))
        class <- use.class(class, class.env)
    }
    class$new(example.1D, rbind(c(0, 1)), 0.5)
}
