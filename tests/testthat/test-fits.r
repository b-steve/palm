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
                                 sv = 0.5, bounds = c(0, 1)))
        expect_that(abs(coef(fit.bin.1D)[1] - 34.908843190) < 1e-4, is_true())
        ## Testing bootstrapping.
        set.seed(5432)
        fit.pois.1D.boot <- boot.ns(fit = fit.pois.1D,
                                    rchild = function(n, child.par){
                                        rpois(n, lambda = child.par)
                                    },
                                    N = 5, prog = FALSE)
        expect_that(abs(mean(fit.pois.1D.boot$boots[, 1]) - 48.787782395) < 1e-4, is_true())

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
                                 sv = 0.5, bounds = c(1e-6, 1)))
        expect_that(abs(coef(fit.bin.2D)[1] - 20.54872790) < 1e-4, is_true())
    })

test_that(
    "Two-plane simulation and model fitting.",
    {
        set.seed(4321)
        lims <- rbind(c(0, 100))
        pars <- c(D = 1, sigma = 0.025, p01 = 0.2, p10 = 0.1)
        ## Simulating data.
        plane.data <- sim.twoplane(pars = pars, lims = lims)
        points <- plane.data$points
        planes <- plane.data$planes
        expect_that(abs(plane.data$points[1, 1] - 75.043889) < 1e-4, is_true())
        siblings <- siblings.twoplane(planes)
        ## Fitting model.
        fit.twoplane <- fit.ns(points = points, lims = lims, R = 0.5,
                      child.dist = list(mean = function(x) 2*0.2/(x + 0.2),
                             var = function(x) 2*x*0.2*(2 - x - 0.2)/(x + 0.2)^2,
                          sv = 0.2, bounds = c(0, 1)), siblings = siblings)
        expect_that(abs(coef(fit.twoplane)[1] - 1.387786) < 1e-4, is_true())
        ## Testing bootstrapping error.
        rchild.twoplane <- function(n, child.par){
            p10 <- child.par
            p01 <- 0.2
            p11 <- 1 - p10
            p00 <- 1 - p01
            p.up <- p01/(p10 + p01)
            out <- sample(c(0, 1), size = n, replace = TRUE, prob = c(1 - p.up, p.up))
            for (i in 1:n){
                if (out[i] == 0){
                    prob <- c(p00, p01)
                } else if (out[i] == 1){
                    prob <- c(p01, p11)
                }
                out[i] <- out[i] + sample(c(0, 1), size = 1, replace = TRUE,
                                          prob = prob)
            }
        }
        expect_that(boot.ns(fit.twoplane, rchild.fun = rchild.twoplane,
                            N = 10, prog = FALSE), throws_error())
    })


