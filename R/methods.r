#' Extract parameter estimates.
#'
#' Extracts estimated parameters from an object returned by
#' \link{fit.ns} or \link{fit.void}.
#'
#' @param object A fitted model object.
#' @param se Logical, if \code{TRUE} standard errors are presented
#'     (if available) instead of parameter estimates.
#' @param ... Other parameters (for S3 generic compatibility).
#'
#' @method coef palm
#' 
#'@export
coef.palm <- function(object, se = FALSE, ...){
    if (se){
        if (is.null(object$boots)){
            stop("Standard errors not available as the model object has not been bootstrapped.")
        }
        out <- apply(object$boots, 2, sd, na.rm = TRUE)
    } else {
        out <- object$par.fitted
    }
    out
}

#' Extract parameter estimates.
#'
#' Extracts estimated parameters from an object returned by
#' \link{fit.ns} with \code{"twocamera"} dispersion or
#' \link{fit.twocamera}.
#'
#' @param object A fitted model object.
#' @param se Logical, if \code{TRUE} standard errors are presented (if
#'     available) instead of parameter estimates.
#' @param report.2D Logical, if \code{TRUE}, two-dimensional density
#'     is reported.
#' @param ... Other parameters (for S3 generic compatibility).
#'
#' @method coef palm_twocamerachild
#' 
#'@export
coef.palm_twocamerachild <- function(object, se = FALSE, report.2D = TRUE, ...){
    if (se){
        if (is.null(object$boots)){
            stop("Standard errors not available as the model object has not been bootstrapped.")
        }
        boots <- get.boots(object, report.2D = report.2D)
        out <- apply(boots, 2, sd, na.rm = TRUE)
    } else {
        out <- object$par.fitted
        if (report.2D){
            which.D <- which(object$par.names == "D")
            out[which.D] <- out[which.D]/(2*object$twocamera.b)
            names(out)[which.D] <- "D.2D"
        }
    }
    out
}

#' Extracts Neyman-Scott point process parameter confidence intervals.
#'
#' Extracts confidence intervals for estimated and derived parameters
#' from a model fitted using \link{fit.ns}, then bootstrapped using
#' \link{boot}.
#'
#' Bootstrap parameter estimates can be found in the
#' \code{boot} component of the model object, so alternative
#' confidence interval methods can be calculated by hand.
#'
#' @param object A fitted model returned by \link{fit.ns},
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
#' @method confint palm
#'
#' @export
confint.palm <- function(object, parm = NULL, level = 0.95, method = "percentile", ...){
    if (is.null(object$boots)){
        stop("Confidence intervals not available as the model object has not been bootstrapped.")
    }
    if (method == "normal"){
        ests <- coef(object, ...)
        ses <- coef(object, se = TRUE, ...)
        out <- cbind(ests + qnorm((1 - level)/2)*ses,
                     ests - qnorm((1 - level)/2)*ses)
    } else if (method == "percentile"){
        boots <- get.boots(object, ...)
        out <- t(apply(boots, 2, quantile, probs = c((1 - level)/2, 1 - (1 - level)/2),
                       na.rm = TRUE))
    }
    colnames(out) <- c(paste(100*(1 - level)/2, "%"),
                       paste(100*((1 + level)/2), "%"))
    if (!is.null(parm)){
        out <- out[parm, , drop = FALSE]
    }
    out
}

## Getting bootstraps.
get.boots <- function(object, ...){
    UseMethod("get.boots")
}
## For normal fits.
#' @method get.boots palm
get.boots.palm <- function(object, ...){
    object$boots
}
## For twocamera fits.
#' @method get.boots palm_twocamerachild
get.boots.palm_twocamerachild <- function(object, report.2D = TRUE, ...){
    boots <- object$boots
    if (report.2D){
        which.D <- which(object$par.names == "D")
        boots[, which.D] <- boots[, which.D]/(2*object$twocamera.b)
        colnames(boots)[which.D] <- "D.2D"
    }
    boots
}

#' Summarising palm model fits
#'
#' Provides a useful summary of the model fit.
#'
#' @param object A fitted model returned by \link{fit.ns}.
#' @param ... Other parameters (for S3 generic compatibility).
#' @inheritParams coef.palm
#' 
#' @method summary palm
#'
#' @export
summary.palm <- function(object, ...){
    coefs <- coef(object, ...)
    ses <- coef(object, se = TRUE, ...)
    out <- list(coefs = coefs, ses = ses)
    class(out) <- c("summary.palm", class(out))
    out
}

#' @method print summary.palm
#' 
#' @export
print.summary.palm <- function(x, ...){
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
#' \link{fit.ns}().
#'
#' @param x A fitted model from \link{fit.ns}.
#' @param xlim Numeric vector giving the x-coordinate range.
#' @param ylim Numeric vector giving the y-coordinate range.
#' @param show.empirical Logical, if \code{TRUE} the empirical Palm
#'     intensity is also plotted.
#' @param breaks The number of breakpoints between cells for the
#'     empirical Palm intensity.
#' @param ... Other parameters (for S3 generic compatibility).
#'
#' @export
plot.palm <- function(x, xlim = NULL, ylim = NULL, show.empirical = TRUE, breaks = 50, ...){
    x$plot(xlim, ylim, show.empirical, breaks, ...)
}


