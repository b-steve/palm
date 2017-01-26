context("Testing model fits")

test_that(
    "1D fitting",
    {
        ## Testing fits with PBC edge correction.
        fit.pois.1D <- fit.ns(example.1D, lims = rbind(c(0, 1)), R = 0.5)
        expect_equal(coef(fit.pois.1D), expected = c(D = 46.545192991,
                                                     lambda = 2.061197102,
                                                     sigma = 0.005430651),
                     tolerance = 0.01)
        fit.bin.1D <- fit.ns(points = example.1D, lims = rbind(c(0, 1)), R = 0.5,
                                   child.dist = "binom4")
        expect_equal(coef(fit.bin.1D), expected = c(D = 34.908870817,
                                                    p = 0.687066075,
                                                    sigma = 0.005430652),
                     tolerance = 0.01)
        ## Testing fit with buffer-zone edge correction.
        fit.pois.1D.buffer.r6 <- fit.ns(example.1D, lims = rbind(c(0, 1)), R = 0.1,
                                           edge.correction = "buffer")
        expect_equal(coef(fit.pois.1D.buffer.r6), expected = c(D = 38.791827160,
                                                               lambda = 1.922548805,
                                                               sigma = 0.005606419),
                     tolerance = 0.01)
        
        ## Testing bootstrapping.
        set.seed(5432)
        fit.pois.1D <- boot(fit.pois.1D, 10, FALSE)
        expect_equal(mean(fit.pois.1D$boots[, 1]), 50.47325, tolerance = 0.01)
    })

test_that(
    "2D fitting",
    {
        ## With Poisson response (in R6).
        fit.pois.2D <- fit.ns(example.2D, lims = rbind(c(0, 1), c(0, 1)), R = 0.5)
        expect_equal(coef(fit.pois.2D), expected = c(D = 41.1579245303134,
                                                     lambda = 0.694235909141885,
                                                     sigma = 0.0270400534209158),
                     tolerance = 0.01)
        ## With binomial response (in R6).             
        fit.bin.2D <- fit.ns(example.2D, lims = rbind(c(0, 1), c(0, 1)), R = 0.5,
                                child.dist = "binom2")
        par.fit.bin.2D <- coef(fit.bin.2D)
        names(par.fit.bin.2D) <- NULL
        expect_equal(par.fit.bin.2D, c(20.54872790, 0.69515914, 0.02698736), tolerance = 0.001)
    })

test_that(
    "Matern fitting",
    {
        fit.matern <- fit.ns(example.2D, lims = rbind(c(0, 1), c(0, 1)), R = 0.5, disp = "uniform")
        par.fit.matern <- coef(fit.matern)
        names(par.fit.matern) <- NULL
        expect_equal(par.fit.matern, c(38.76459565, 0.73561864, 0.05327295), tolerance = 0.001)
    })

test_that(
    "Two-plane simulation and model fitting.",
    {
        set.seed(4321)
        ## Simulation.
        twoplane.data <- sim.ns(c(D = 0.9615, kappa = 28, sigma = 0.025), lims = rbind(c(0, 500)),
                                   disp = "gaussian", child.dist = "twoplane",
                                   child.info = list(w = 0.175, b = 0.175 + 5*0.025, l = 20, tau = 110))
        expect_equal(twoplane.data$points[1, 1], expected = 399.969972, tolerance = 0.01)
        ## With plane information.
        fit <- fit.ns(points = twoplane.data$points, lims = rbind(c(0, 500)), R = 1,
                         child.dist = "twoplane",
                         child.info = list(w = 0.175, b = 0.175 + 5*0.025, l = 20, tau = 110),
                         sibling.list = twoplane.data$sibling.list)
        expect_equal(coef(fit, report.2D = FALSE), expected = c(D = 0.81322222, kappa = 26.02728291, sigma = 0.02046535),
                     tolerance = 0.01)
        ## Without plane information.
        fit <- fit.ns(points = twoplane.data$points, lims = rbind(c(0, 500)), R = 1,
                         child.dist = "twoplane",
                         child.info = list(w = 0.175, b = 0.175 + 5*0.025, l = 20, tau = 110))
        expect_equal(coef(fit, report.2D = FALSE), expected = c(D = 0.76768111, kappa = 26.91820235, sigma = 0.02054313),
                     tolerance = 0.01)
        fit.noplane <- fit.twoplane(points = twoplane.data$points, NULL, 500, 0.175, 0.175 + 5*0.025, 20, 110, 1)
        expect_equal(coef(fit.noplane, report.2D = FALSE), expected = c(D = 0.76768111, kappa = 26.91820235, sigma = 0.02054313),
                     tolerance = 0.01)
    })

