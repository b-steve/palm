## Reference class stuff.

######
## Overall class containing base stuff.
######
base.class <- setRefClass("base", fields = c("points",
                                                 "lims",
                                                 "R",
                                                 "par.names",
                                                 "par.links",
                                                 "par.invlinks",
                                                 "start.pars",
                                                 "link.start.pars"))
## A method for converting functions to their link scale.
base.class$methods(link.pars = function(pars){
    n.pars <- length(pars)
    out <- numeric(n.pars)
    for (i in 1:n.pars){
        out[i] <- par.links[[i]](pars[i])
    }
    out
})
## A method for converting functions to their real scale.
base.class$methods(invlink.pars = function(pars){
    n.pars <- length(pars)
    out <- numeric(n.pars)
    for (i in 1:n.pars){
        out[i] <- par.invlinks[[i]](pars[i])
    }
    out
}
## A method for generating link.start.pars.
base.class$methods(initialise.link.start.pars(){
    link.start.pars <<- link.pars(start.pars)
})
## An empty method for the Palm intensity.
base.class$methods(palm.intensity = function(pars){})
## An empty method for the Palm likelihood function.
base.class$methods(palm.likelihood = function(pars){})
## A method for the negative Palm likelihood function.
base.class$methods(neg.palm.likelihood = function(pars){
    -palm.likelihood(pars)
})
## A method for model fitting.
base.class$methods(fit = function(.self) optim(.self$start.pars, .self$nloglik))

######
## Class for Neyman-Scott processes.
######
ns.class <- setRefClass("ns", fields = "child.dist", contains = "base")
## Overwriting method for the Palm intensity.
ns.class$methods(palm.intensity = function(pars) sibling.pi(pars) + nonsibling.pi(pars))

######
## Class for sibling models.
######
sibling.class <- setRefClass("sibling", fields = c("sibling.mat",
                                                   "sibling.pT",
                                                   "sibling.fF"),
                             contains = "ns")

######
## Class for Thomas processes.
######
thomas.class <- setRefClass("thomas", contains = "ns")

######
## Class for Matern processes.
######
matern.class <- setRefClass("matern", contains = "ns")


