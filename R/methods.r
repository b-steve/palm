#' Extract Neyman-Scott point process parameter estiamtes
#'
#' Extracts estimated and derived parameters from a model fitted using
#' \link{fit.ns}().
#'
#' @param object A fitted model from \link{fit.ns}().
#' @param ... Other parameters (for S3 generic compatibility).
#'
#' @method coef nspp
#'
#' @export
coef.nspp <- function(object, ...){
    object$pars
}
