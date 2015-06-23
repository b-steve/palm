#' Bootstrapping a Neymann-Scott point process model
#'
#' Carries out a bootstrap for Neymann-Scott point process models
#' fitted by \link{fit.ns}().
#'
#' @return The first argument, with added information from the
#' bootstrap procedure.
#'
#' @param fit A fitted object from \link{fit.ns}().
#' @inheritParams sim.ns
#' @param N The number of bootstrap resamples.
#' 
#' @export
boot.ns <- function(fit, rchild, N, ...){
    ## Extracting information.
    pars <- fit$pars
    n.pars <- length(pars)
    lims <- fit$args$lims
    boots <- matrix(0, nrow = N, ncol = n.pars)
    for (i in 1:N){
        args$points <- sim.ns(pars = pars[c("D", "sigma")], lims = lims,
                              rchild = rchild)
        args$child.dist$sv <- pars["child.par"]
    }
}
