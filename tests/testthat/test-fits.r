context("Testing model fits")

test_that(
    "1D fitting",
    {
        fit.pois.1D <- fit.ns(points = example.1D, lims = rbind(c(0, 1)), R = 0.5,
                              child.dist = list(mean = function(x) x,
                                                var = function(x) x,
                                                sv = 5, bounds = c(0, 1e8)))
        expect_that(abs(coef(fit.pois.1D)[1] - 46.530524823) < 1e-4, is_true())
        fit.bin.1D <- fit.ns(points = example.1D, lims = rbind(c(0, 1)), R = 0.5,
                             child.dist = list(mean = function(x) 4*x,
                                 var = function(x) 4*x*(1 - x),
                                 sv = 0.5, bounds = c(0, 1)), trace = TRUE)
        expect_that(abs(coef(fit.bin.1D)[1] - 34.907561984) < 1e-4, is_true())
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
        fit.pois.2D <- fit.ns(points = example.2D, lims = rbind(c(0, 1), c(0, 1)), R = 0.5,
                              child.dist = list(mean = function(x) x, var = function(x) x,
                                  sv = 20, bounds = c(1e-6, nrow(example.2D))))
        expect_that(abs(coef(fit.pois.2D)[1] - 41.10442073) < 1e-4, is_true())
        fit.bin.2D <- fit.ns(points = example.2D, lims = rbind(c(0, 1), c(0, 1)), R = 0.5,
                             child.dist = list(mean = function(x) 2*x,
                                 var = function(x) 2*x*(1 - x),
                                 sv = 0.5, bounds = c(0, 1)))
        expect_that(abs(coef(fit.bin.2D)[1] - 20.54872790) < 1e-4, is_true())
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
        expect_that(abs(coef(fit)[1] - 0.4958427) < 1e-4, is_true())
    })


