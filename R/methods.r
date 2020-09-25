#' Extract parameter estimates.
#'
#' Extracts estimated parameters from an object returned by the
#' fitting functions in this package, such as \link{fit.ns},
#' \link{fit.void}, and \link{fit.twocamera}.
#'
#' @param object A fitted model object.
#' @param se Logical, if \code{TRUE} standard errors are presented
#'     (if available) instead of parameter estimates.
#' @param ... Other parameters (for S3 generic compatibility).
#'
#' @method coef palm
#' 
#' @export
coef.palm <- function(object, se = FALSE, ...){
    if (se){
        if (is.null(object$boots)){
            warning("Standard errors not available as the model object has not been bootstrapped.")
            out <- rep(NA, length(object$par.fitted))
        } else {
            boots <- get.boots(object)
            out <- apply(boots, 2, sd, na.rm = TRUE)
        }
    } else {
        out <- object$par.fitted
    }
    out
}

#' @param report.2D Logical, for two-camera model fits only. If
#'     \code{TRUE}, two-dimensional animal density is reported.
#'
#' @method coef palm_twocamerachild
#'
#' @rdname coef.palm
#'
#' @export
coef.palm_twocamerachild <- function(object, se = FALSE, report.2D = TRUE, ...){
    if (se){
        if (is.null(object$boots)){
            warning("Standard errors not available as the model object has not been bootstrapped.")
            out <- rep(NA, length(object$par.fitted))
        } else {
            boots <- get.boots(object, report.2D = report.2D)
            out <- apply(boots, 2, sd, na.rm = TRUE)
        }
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
#' from a model fitted using \link{fit.ns}, \link{fit.void}, or
#' \link{fit.twocamera}, then bootstrapped using \link{boot.palm}.
#'
#' Bootstrap parameter estimates can be found in the \code{boots}
#' component of the model object, so alternative confidence interval
#' methods can be calculated by hand.
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
#' @param ... Other parameters (for S3 generic compatibility).
#'
#' @method confint palm
#'
#' @examples
#' ## Fitting model.
#' fit <- fit.ns(example.2D, lims = rbind(c(0, 1), c(0, 1)), R = 0.5)
#' ## Carrying out bootstrap.
#' fit <- boot.palm(fit, N = 100)
#' ## Calculating 95% confidence intervals.
#' confint(fit)
#' ## Estimates are very imprecise---these data were only used as
#' ## they can be fitted and bootstrapped quickly for example purposes.
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
#' \link{fit.ns}.
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
#' @examples
#' ## Fit model.
#' fit <- fit.ns(example.2D, lims = rbind(c(0, 1), c(0, 1)), R = 0.5)
#' ## Plot fitted Palm intensity.
#' plot(fit)
#' 
#' @export
plot.palm <- function(x, xlim = NULL, ylim = NULL, show.empirical = TRUE, breaks = 50, ...){
    x$plot(xlim, ylim, show.empirical, breaks, ...)
}


