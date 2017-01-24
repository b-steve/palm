#' Bootstrapping a Neymann-Scott point process model
#'
#' Carries out a bootstrap for Neymann-Scott point process models
#' fitted by \link{fit.ns}.
#'
#' The \code{rchild} function may only take a single distributional
#' parameter. If the distribution for the number of children generated
#' by each parent is Poisson, then the native \code{rpois} is
#' appropriate, as this distribution has a single parameter. For
#' distributions with two or more parameters, those other than
#' \code{child.par} must be hard-coded into \code{rchild}. For
#' example, if a Binomial(n, 2, p) is required, then \code{function(n,
#' p) rbinom(n = n, size = 2, prob = p)} would be an appropriate
#' function for \code{rchild}.
#'
#' @return The first argument, with added information from the
#' bootstrap procedure.
#'
#' @param fit A fitted object from \link{fit.ns}().
#' @param N The number of bootstrap resamples.
#' @param prog Logical, if \code{TRUE}, a progress bar is printed to
#' the console.
#' @inheritParams sim.ns
#'
#' @examples
#' \dontrun{
#' ## Fitting a model.
#' fit <- fit.ns(example.2D, lims = rbind(c(0, 1), c(0, 1)), R = 0.5)
#' ## Bootstrapping.
#' fit.boot <- boot.ns(fit, rchild = rpois, N = 100)
#' }
#' 
#' @export
boot.ns <- function(fit, rchild, N, prog = TRUE){
    ## Setting up progress bar.
    if (prog){
        pb <- txtProgressBar(min = 0, max = N, style = 3)
    }
    ## Shitty method dispatch.
    if (class(fit)[1] == "nspp"){
        args <- fit$args
        lims <- args$lims
        ## Error for fits with known (non-)siblings.
        if (!is.null(args$siblings)){
            stop("Bootstrapping not implemented for models with known (non-)siblings.")
        }
    } else if (class(fit)[1] == "twoplane.nspp"){
        args <- fit$args.twoplane
    } else {
        stop("Model class not recognised.")
    }
    pars <- fit$pars
    n.pars <- length(pars)
    boots <- matrix(0, nrow = N, ncol = n.pars)
    ## More shitty method dispatch.
    if (class(fit)[1] == "nspp"){
        for (i in 1:N){
            args$points <- sim.ns(pars = pars[c("D", "sigma", "child.par")], lims = lims,
                                  rchild = rchild)
            args$sigma.sv <- pars["sigma"]
            args$child.dist$sv <- pars["child.par"]
            args$trace <- FALSE
            fit.boot <- do.call("fit.ns", args)
            boots[i, ] <- fit.boot$pars
            ## Updating progress bar.
            if (prog){
                setTxtProgressBar(pb, i)
            }
        }
    } else if (class(fit)[1] == "twoplane.nspp"){
        for (i in 1:N){
            sim.obj <- sim.twoplane(D = pars["D.2D"], sigma = pars["sigma"], S = pars["child.par"],
                                    l = args$l, w = args$w, b = args$b, t = args$t, C = args$C)
            args$points <- sim.obj$points
            if (!is.null(args$planes)){
                args$planes <- sim.obj$planes
            }
            fit.boot <- do.call("fit.twoplane", args)
            boots[i, ] <- fit.boot$pars
            ## Updating progress bar.
            if (prog){
                setTxtProgressBar(pb, i)
            }
        }
    }
    if (prog){
        close(pb)
    }
    colnames(boots) <- names(pars)
    fit$boots <- boots
    fit$se <- apply(boots, 2, sd)
    class(fit) <- c("boot.nspp", class(fit))
    fit
}

