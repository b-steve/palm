#' Extract Neyman-Scott point process parameter estimates
#'
#' Extracts estimated and derived parameters from a model fitted using
#' \link{fit.ns}.
#'
#' @param object A fitted model from \link{fit.ns}.
#' @param all Logical, if \code{TRUE} derived parameters are also
#' included.
#' @param ... Other parameters (for S3 generic compatibility).
#'
#' @method coef nspp
#'
#' @export
coef.nspp <- function(object, all = FALSE, ...){
    if (all){
        parm <- 1:6
    } else {
        parm <- 1:3
    }
    object$pars[parm]
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
#' @param method A character string specifying the method used to
#' calculate confidence intervals. Choices are "normal", for a normal
#' approximation, and "percentile", for the percentile method.
#' @param ... Other parameters (for S3 generic compatability).
#'
#' @method confint boot.nspp
#'
#' @export
confint.boot.nspp <- function(object, parm = c("D", "sigma", "child.par"), level = 0.95, method = "percentile", ...){
    if (is.null(parm)){
        parm <- 1:length(object$se)
    }
    if (method == "normal"){
        ests <- coef(object)[parm]
        ses <- object$se[parm]
        out <- cbind(ests + qnorm((1 - level)/2)*ses,
                     ests - qnorm((1 - level)/2)*ses)
    } else if (method == "percentile"){
        out <- t(apply(object$boots[, parm, drop = FALSE], 2, quantile, probs = c((1 - level)/2, 1 - (1 - level)/2),
                       na.rm = TRUE))
    }
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
#' @inheritParams coef.nspp
#' 
#' @method summary nspp
#'
#' @export
summary.nspp <- function(object, all = FALSE, ...){
    if (all){
        parm <- 1:6
    } else {
        parm <- 1:3
    }
    coefs <- coef(object)[parm]
    ses <- object$se[parm]
    if (is.null(ses)){
        ses <- rep(NA, length(coefs))
        names(ses) <- names(coefs)[parm]
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
