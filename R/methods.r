#' Extract Neyman-Scott point process parameter estiamtes
#'
#' Extracts estimated and derived parameters from a model fitted using
#' \link{fit.ns}().
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
