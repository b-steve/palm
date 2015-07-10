#' Extract Neyman-Scott point process parameter estimates
#'
#' Extracts estimated and derived parameters from a model fitted using
#' \link{fit.ns}.
#'
#' @param object A fitted model from \link{fit.ns}.
#' @param ... Other parameters (for S3 generic compatibility).
#'
#' @method coef nspp
#'
#' @export
coef.nspp <- function(object, ...){
    object$pars
}

#' Extracts Neyman-Scott point process parameter confidence intervals.
#'
#' Extracts confidence intervals for estimated and derived parameters
#' from a model fitted using \link{fit.ns}, then bootstrapped using
#' \link{boot.ns}.
#'
#' Calculation of confidence intervals is via a normal approximation,
#' whereby standard errors are calculated from the standard deviations
#' of the parameter estimates across the bootstrap
#' resamples. Bootstrap parameter estimates can be found in the
#' \code{boot} component of the model object, so alternative
#' confidence interval methods can be calculated by hand.
#'
#' @param object A fitted model from \link{fit.ns}, bootstrapped using
#' \link{boot.ns}.
#' @param parm A vector of parameter names, specifying which
#' parameters are to be given confidence intervals.
#' @param level The confidence level required.
#' @param ... Other parameters (for S3 generic compatability).
#'
#' @method confint boot.nspp
#'
#' @export
confint.boot.nspp <- function(object, parm = NULL, level = 0.95, ...){
    if (is.null(parm)){
        parm <- 1:length(object$se)
    }
    ests <- coef(object)[parm]
    ses <- object$se[parm]
    out <- cbind(ests + qnorm((1 - level)/2)*ses,
                 ests - qnorm((1 - level)/2)*ses)
    colnames(out) <- c(paste(100*(1 - level)/2, "%"),
                       paste(100*((1 + level)/2), "%"))
    out
}

#' Summarising nspp model fits
#'
#' Provides a useful summary of the model fit.
#'
#' @param object A fitted model from \link{fit.ns}.
#' @param ... Other parameters (for S3 generic compatibility).
#'
#' @method summary nspp
#'
#' @export
summary.nspp <- function(object, ...){
    coefs <- coef(object)
    ses <- object$se
    if (is.null(ses)){
        ses <- rep(NA, length(coefs))
        names(ses) <- names(coefs)
    }
    out <- list(coefs = coefs, ses = ses)
    class(out) <- c("summary.nspp", class(out))
    out
}

#' @method print summary.nspp
#'
#' @export
print.summary.nspp <- function(x, ...){
    n.coefs <- length(x$coefs)
    mat <- matrix(0, nrow = n.coefs, ncol = 2)
    rownames(mat) <- names(x$coefs)
    colnames(mat) <- c("Estimate", "Std. Error")
    mat[, 1] <- x$coefs
    mat[, 2] <- x$ses
    cat("Coefficients:", "\n")
    printCoefmat(mat, na.print = "")
}
