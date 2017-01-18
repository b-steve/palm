context("Testing model fits")

test_that(
    "1D fitting",
    {
        ## Testing fits with PBC edge correction.
        fit.pois.1D <- fit.ns(points = example.1D, lims = rbind(c(0, 1)), R = 0.5,
                              child.dist = list(mean = function(x) x,
                                                var = function(x) x,
                                                sv = 5, bounds = c(0, 1e8)))
        expect_that(abs(coef(fit.pois.1D)[1] - 46.530524823) < 1e-4, is_true())
        fit.pois.1D.r6 <- fit.ns_r6(example.1D, lims = rbind(c(0, 1)), R = 0.5, trace = FALSE)
        par.pois.1D.r6 <- coef(fit.pois.1D.r6)
        diff.ratio.1D <- par.pois.1D.r6[c(1, 3, 2)]/coef(fit.pois.1D)
        names(diff.ratio.1D) <- NULL
        expect_equal(diff.ratio.1D, expected = rep(1, 3), tolerance = 0.01)
        fit.bin.1D <- fit.ns(points = example.1D, lims = rbind(c(0, 1)), R = 0.5,
                             child.dist = list(mean = function(x) 4*x,
                                 var = function(x) 4*x*(1 - x),
                                 sv = 0.5, bounds = c(0, 1)))
        expect_that(abs(coef(fit.bin.1D)[1] - 34.907561984) < 1e-4, is_true())
        ## Testing fits with buffer-zone edge correction.
        fit.pois.1D.buffer <- fit.ns(points = example.1D, lims = rbind(c(0, 1)), R = 0.1,
                                     child.dist = list(mean = function(x) x,
                                                       var = function(x) x,
                                                       sv = 5, bounds = c(0, 1e8)),
                                     edge.correction = "buffer")
        expect_that(abs(coef(fit.pois.1D.buffer)[1] - 38.7933732) < 1e-4, is_true())
        ## Testing bootstrapping.
        set.seed(5432)
        fit.pois.1D.boot <- boot.ns(fit = fit.pois.1D,
                                    rchild = function(n, child.par){
                                        rpois(n, lambda = child.par)
                                    },
                                    N = 5, prog = FALSE)
        expect_that(abs(mean(fit.pois.1D.boot$boots[, 1]) - 48.78834) < 1e-4, is_true())
        ## Testing plotting.
        expect_that(plot(fit.pois.1D.boot), is_null())
    })

test_that(
    "2D fitting",
    {
        ## With Poisson response (in R6).
        fit.pois.2D <- fit.ns_r6(example.2D, lims = rbind(c(0, 1), c(0, 1)), R = 0.5)
        pars.fit.pois.2D <- coef(fit.pois.2D)
        names(pars.fit.pois.2D) <- NULL
        expect_equal(pars.fit.pois.2D,
                     c(41.1579245303134, 0.694235909141885, 0.0270400534209158),
                     tolerance = 0.01)
        ## With binomial response (in R6).             
        fit.bin.2D <- fit.ns_r6(example.2D, lims = rbind(c(0, 1), c(0, 1)), R = 0.5,
                                child.dist = "binom2")
        par.fit.bin.2D <- coef(fit.bin.2D)
        names(par.fit.bin.2D) <- NULL
        expect_equal(par.fit.bin.2D, c(20.54872790, 0.69515914, 0.02698736), tolerance = 0.001)
    })

test_that(
    "Matern fitting",
    {
        fit.matern <- fit.ns_r6(example.2D, lims = rbind(c(0, 1), c(0, 1)), R = 0.5, disp = "uniform")
        par.fit.matern <- coef(fit.matern)
        names(par.fit.matern) <- NULL
        expect_equal(par.fit.matern, c(39.163600, 0.7284017, 0.0524428), tolerance = 0.001)
    })

test_that(
    "Two-plane simulation and model fitting.",
    {
        set.seed(4321)
        D.2D <- 1.6025; sigma <- 0.025; S <- 28; l <- 500; w <- 0.175;
        b <- w + 5*sigma; t = 20; C = 110
        ## Simulating data.
        plane.data <- sim.twoplane(D = D.2D, sigma = sigma, S = S, l = l, w = w,
                                   b = b, t = t, C = C)
        points <- plane.data$points
        planes <- plane.data$planes
        expect_that(abs(plane.data$points[1, 1] - 379.5013) < 1e-4, is_true())
        ## Fitting model.
        fit <- fit.twoplane(points = points, planes = planes, l = l, w = w,
                            b = b, t = t, C = C, R = 1)
        ## coef(fit)[1]/(2*b) ## This is an estimate of D.2D.
        expect_that(abs(coef(fit, all = TRUE)[7] - 0.4958427) < 1e-4, is_true())
        ## Test using R6.

        sibling.list <- siblings.twoplane(planes)
        names(sibling.list) <- c("sibling.mat", "alpha", "beta")
        
        fit.r6 <- fit.ns_r6(points = points, lims = rbind(c(0, l)), R = 1, child.dist = "twoplane",
                            child.info = list(w = w, b = b, l = t, tau = C),
                            sibling.list = sibling.list, trace = TRUE)
        pars.old <- coef(fit)[c(1, 3, 2)]
        pars.new <- coef(fit.r6)
        pars.new[1] <- pars.new[1]/(2*b)
        names(pars.old) <- names(pars.new) <- NULL
        expect_equal(pars.old, pars.new, tolerance = 0.001)
    })
