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
        ## Testing the above with bobyqa().
        fit.bobyqa.pois.1D <- fit.ns(example.1D, lims = rbind(c(0, 1)), R = 0.5, use.bobyqa = TRUE)
        expect_equal(coef(fit.bobyqa.pois.1D), expected = c(D = 46.545192991,
                                                     lambda = 2.061197102,
                                                     sigma = 0.005430651),
                     tolerance = 0.01)
        ## Testing the above with a binomial number of children.
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
        fit.pois.1D <- boot.palm(fit.pois.1D, 10, FALSE)
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
        ## With bobyqa(), as nlminb() appears more unstable for this.
        fit.bobyqa.matern <- fit.ns(example.2D, lims = rbind(c(0, 1), c(0, 1)), R = 0.5, disp = "uniform",
                             use.bobyqa = TRUE)
        par.fit.bobyqa.matern <- coef(fit.bobyqa.matern)
        names(par.fit.bobyqa.matern) <- NULL
        expect_equal(par.fit.bobyqa.matern, c(38.88561077, 0.73344422, 0.05298069), tolerance = 0.001)
    })

test_that(
    "Two-camera simulation and model fitting.",
    {
        set.seed(4321)
        ## Simulation.
        twocamera.data <- sim.ns(c(D = 0.9615, kappa = 28, sigma = 0.025), lims = rbind(c(0, 500)),
                                   disp = "gaussian", child.dist = "twocamera",
                                   child.info = list(w = 0.175, b = 0.175 + 5*0.025, l = 20, tau = 110))
        expect_equal(twocamera.data$points[1, 1], expected = 399.969972, tolerance = 0.01)
        ## With camera information.
        fit <- fit.ns(points = twocamera.data$points, lims = rbind(c(0, 500)), R = 1,
                      child.dist = "twocamera",
                      child.info = list(w = 0.175, b = 0.175 + 5*0.025, l = 20, tau = 110),
                      sibling.list = twocamera.data$sibling.list)
        expect_equal(coef(fit, report.2D = FALSE), expected = c(D = 1.74321506, kappa = 19.75875843, sigma = 0.02202638),
                     tolerance = 0.001)
        fit.camera <- fit.twocamera(points = twocamera.data$points, cameras = twocamera.data$sibling.list$cameras,
                                    d = 500, w = 0.175, b = 0.175 + 5*0.025, l = 20, tau = 110, R = 1)
        expect_equal(coef(fit.camera, report.2D = FALSE), expected = c(D = 1.74321506, kappa = 19.75875843, sigma = 0.02202638),
                     tolerance = 0.001)
        ## Without camera information.
        fit <- fit.ns(points = twocamera.data$points, lims = rbind(c(0, 500)), R = 1,
                         child.dist = "twocamera",
                         child.info = list(w = 0.175, b = 0.175 + 5*0.025, l = 20, tau = 110))
        expect_equal(coef(fit, report.2D = FALSE), expected = c(D = 1.97975828, kappa = 18.08215974, sigma = 0.02144923),
                     tolerance = 0.001)
        fit.nocamera <- fit.twocamera(points = twocamera.data$points, NULL, 500, 0.175, 0.175 + 5*0.025, 20, 110, 1)
        expect_equal(coef(fit.nocamera, report.2D = FALSE), expected = c(D = 1.97975828, kappa = 18.08215974, sigma = 0.02144923),
                     tolerance = 0.001)
    })

test_that(
    "Void process simulation and model fitting.",
    {
        set.seed(2468)
        ## Simulation.
        void.data <- sim.void(c(Dc = 500, Dp = 100, tau = 0.05), rbind(c(0, 1), c(0, 1)))
        expect_equal(void.data$points[22, 1], expected = 0.2294614180, tolerance = 1e-6)
        ## Model fitting.
        fit <- fit.void(void.data$points, lims = rbind(c(0, 1), c(0, 1)), R = 0.5,
                        start = c(Dc = 1000, Dp = 20, tau = 0.075))
        expect_equal(fit$log.palm.likelihood(fit$par.fitted),
                     expected = 166313.248625127, tolerance = 0.001)
    })

test_that(
    "Multipattern simulation and model fitting.",
    {
        set.seed(3210)
        ## Spatial domain limits.
        multi.lims <- list(rbind(c(0, 1), c(0, 1)), rbind(c(0, 1), c(0, 2)))
        ## Simulation.
        multipattern.data <- sim.ns(c(D = 5, lambda = 10, sigma = 0.025), lims = multi.lims)
        expect_equal(multipattern.data[[1]]$points[1, 1], expected = 0.3199291, tolerance = 1e-6)
        expect_equal(multipattern.data[[2]]$points[1, 1], expected = 0.4434522, tolerance = 1e-6)
        ## Model fitting.
        fit <- fit.ns(lapply(multipattern.data, function(x) x$points), lims = multi.lims, R = 0.5)
        expect_equal(coef(fit), expected = c(D = 3.4247809, lambda = 9.3783646, sigma = 0.0203175),
                     tolerance = 1e-3)
    })
