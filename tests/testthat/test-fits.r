context("Testing model fits")

test_that(
    "1D fitting",
    {
        fit.pois.1D <- fit.ns(points = example.1D, lims = rbind(c(0, 1)), R = 0.5,
                              child.dist = list(mean = function(x) x,
                                  var = function(x) x,
                                  sv = 5, bounds = c(0, 1e8)), trace = TRUE)
        expect_that(abs(coef(fit.pois.1D)[1] - 46.530524823) < 1e-4, is_true())
        fit.bin.1D <- fit.ns(points = example.1D, lims = rbind(c(0, 1)), R = 0.5,
                             child.dist = list(mean = function(x) 4*x,
                                 var = function(x) 4*x*(1 - x),
                                 sv = 0.5, bounds = c(0, 1)), trace = TRUE)
        expect_that(abs(coef(fit.bin.1D)[1] - 34.908843190) < 1e-4, is_true())
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
