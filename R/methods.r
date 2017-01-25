#' Extract parameter estimates.
#'
#' Extracts estimated parameters from an object returned by
#' \link{fit.ns_r6}.
#'
#' @param object A fitted model object returned by \link{fit.ns_r6}.
#' @param se Logical, if \code{TRUE} standard errors are presented
#'     (if available) instead of parameter estimates.
#' @param ... Other parameters (for S3 generic compatibility).
#'
#' @method coef nspp_r6
#' 
#'@export
coef.nspp_r6 <- function(object, se = FALSE, ...){
    if (se){
        if (is.null(object$boots)){
            stop("Standard errors not available as the model object has not been bootstrapped.")
        }
        out <- object$par.se
    } else {
        out <- object$par.fitted
    }
    out
}

#' Extracts Neyman-Scott point process parameter confidence intervals.
#'
#' Extracts confidence intervals for estimated and derived parameters
#' from a model fitted using \link{fit.ns}, then bootstrapped using
#' \link{boot.ns}.
#'
#' Bootstrap parameter estimates can be found in the
#' \code{boot} component of the model object, so alternative
#' confidence interval methods can be calculated by hand.
#'
#' @param object A fitted model returned by \link{fit.ns_r6},
#'     bootstrapped using \link{boot}.
#' @param parm A vector of parameter names, specifying which
#'     parameters are to be given confidence intervals. Defaults to
#'     all parameters.
#' @param level The confidence level required.
#' @param method A character string specifying the method used to
#'     calculate confidence intervals. Choices are "normal", for a
#'     normal approximation, and "percentile", for the percentile
#'     method.
#' @param ... Other parameters (for S3 generic compatability).
#'
#' @method confint nspp_r6
#'
#' @export
confint.nspp_r6 <- function(object, parm = NULL, level = 0.95, method = "percentile", ...){
    if (is.null(object$boots)){
        stop("Confidence intervals not available as the model object has not been bootstrapped.")
    }
    if(is.null(parm)){
        parm <- object$par.names
    }
    if (method == "normal"){
        ests <- coef(object)
        ses <- coef(object, se = TRUE)
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
#' @param object A fitted model returned by \link{fit.ns_r6}.
#' @param ... Other parameters (for S3 generic compatibility).
#' @inheritParams coef.nspp
#' 
#' @method summary nspp_r6
#'
#' @export
summary.nspp_r6 <- function(object, all = FALSE, ...){
    coefs <- coef(object)
    ses <- coef(object, se = TRUE)
    out <- list(coefs = coefs, ses = ses)
    class(out) <- c("summary.nspp_r6", class(out))
    out
}

#' @method print summary.nspp_r6
#' 
#' @export
print.summary.nspp_r6 <- function(x, ...){
    n.coefs <- length(x$coefs)
    mat <- matrix(0, nrow = n.coefs, ncol = 2)
    rownames(mat) <- names(x$coefs)
    colnames(mat) <- c("Estimate", "Std. Error")
    mat[, 1] <- x$coefs
    mat[, 2] <- x$ses
    cat("Coefficients:", "\n")
    printCoefmat(mat, na.print = "")
}

#' Plotting an estimated Palm intensity function.
#'
#' Plots a fitted Palm intensity function from an object returned by
#' \link{fit.ns_r6}().
#'
#' @param x A fitted model from \link{fit.ns_r6}.
#' @param ... Other parameters (for S3 generic compatibility).
#'
#' @export
plot.nspp_r6 <- function(x, ...){
    x$plot(...)
}


