## Reference class stuff.

######
## Overall class containing base stuff.
######
base.class <- setRefClass("base", fields = c("points",
                                             "n.points",
                                             "contrasts",
                                             "lims",
                                             "dim",
                                             "R",
                                             "par.links",
                                             "par.invlinks",
                                             "start.pars",
                                             "link.start.pars"))
## Initialisation.
base.class$methods(initialize = function(...){
    n.points <<- nrow(points)
})
## A method for converting functions to their link scale.
base.class$methods(link.pars = function(pars){
    n.pars <- length(pars)
    out <- numeric(n.pars)
    for (i in 1:n.pars){
        out[i] <- par.links[[i]](pars[i])
    }
    names(out) <- names(pars)
    out
})
## A method for converting functions to their real scale.
base.class$methods(invlink.pars = function(pars){
    n.pars <- length(pars)
    out <- numeric(n.pars)
    for (i in 1:n.pars){
        out[i] <- par.invlinks[[i]](pars[i])
    }
    names(out) <- names(pars)
    out
})
## A method for generating link.start.pars.
base.class$methods(initialise.link.start.pars = function(){
    link.start.pars <<- link.pars(start.pars)
})
## An empty method for the Palm intensity.
base.class$methods(palm.intensity = function(r, pars){})
## A default method for the sum of the log Palm intensities.
base.class$methods(sum.log.intensities = function(pars){
    sum(log(n.points*palm.intensity(contrasts, pars)))
})
## A default method for the integral in the Palm likelihood.
base.class$methods(palm.likelihood.integral = function(pars){
    f <- function(r, pars){
        palm.intensity(r, pars)*Sd(r, dim)
    }
    -n.points*integrate(f, lower = 0, upper = R, pars = pars)$value
})
## A method for the Palm likelihood function.
base.class$methods(log.palm.likelihood = function(pars){
    sum.log.intensities(pars) + palm.likelihood.integral(pars)
})
## A method for the negative Palm likelihood function.
base.class$methods(neg.log.palm.likelihood = function(pars){
    -log.palm.likelihood(pars)
})
## A method for the negative Palm likelihood function with linked parameters.
base.class$methods(link.neg.log.palm.likelihood = function(link.pars){
    pars <- invlink.pars(link.pars)
    neg.log.palm.likelihood(pars)
})
## A method for model fitting.
base.class$methods(fit = function(.self){
    optim(.self$link.start.pars, .self$link.neg.log.palm.likelihood)
})

######
## Class for periodic boundary conditions.
######
pbc.class <- setRefClass("pbc", contains = "base")
## Creating contrasts.
pbc.class$methods(get.contrasts = function(){
    contrasts.save <- pbc_distances(points = points, lims = lims)
    contrasts <<- contrasts.save[contrasts.save <= R]
})

######
## Class for Neyman-Scott processes.
######
ns.class <- setRefClass("ns", contains = "base")
## Overwriting method for the Palm intensity.
ns.class$methods(palm.intensity = function(r, pars) sibling.pi(r, pars) + nonsibling.pi(pars))
## An empty method for the expectation of the child distribution (filled by fit.ns()).
ns.class$methods(child.expectation = function(pars){})
## An empty method for the variance of the child distribution (filled by fit.ns()).
ns.class$methods(child.variance = function(pars){})
## A method for the expected number of siblings from a randomly chosen point.
ns.class$methods(sibling.expectation = function(pars){
    (child.variance(pars) + child.expectation(pars)^2)/child.expectation(pars) - 1
})
## An empty method for the PDF of Q, the between-sibling distances.
ns.class$methods(fq = function(r, pars){})
## An empty method for the CDF of Q.
ns.class$methods(Fq = function(r, pars){})
## A default method for the quotient of the PDF of Q and the surface volume.
ns.class$methods(q.over.s = function(r, pars) fq(r, pars)/Sd(r, dim))
## A method for the PI of nonsibling points.
ns.class$methods(nonsibling.pi = function(pars) pars["D"]*child.expectation(pars))
## A method for the PI of sibling points.
ns.class$methods(sibling.pi = function(r, pars) sibling.expectation(pars)*q.over.s(r, pars))

######
## Class for Thomas processes.
######
thomas.class <- setRefClass("thomas", contains = "ns")
## Overwriting method for the PDF of Q.
thomas.class$methods(fq = function(r, pars){
    2^(1 - dim/2)*r^(dim - 1)*exp(-r^2/(4*pars["sigma"]^2))/((pars["sigma"]*sqrt(2))^dim*gamma(dim/2))
})
## Overwriting method for the CDF of Q.
thomas.class$methods(Fq = function(r, pars){
    pgamma(r^2/(4*pars["sigma"]^2), dim/2)
})
## Overwriting method for the quotient of the PDF of Q and the surface volume.
thomas.class$methods(q.over.s = function(r, pars){
    exp(-r^2/(4*pars["sigma"]^2))/((2*pars["sigma"])^dim*pi^(dim/2))
})

######
## Class for Matern processes.
######
matern.class <- setRefClass("matern", contains = "ns")

######
## Class for sibling models.
######
sibling.class <- setRefClass("sibling", fields = c("sibling.mat",
                                                   "sibling.pT",
                                                   "sibling.fF"),
                             contains = "ns")


