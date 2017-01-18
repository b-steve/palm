## Functions for creating the sibling matrix.

#' Creating a sibling matrix for two-plane surveys
#'
#' Creates an object suitable for the \code{siblings} argument to
#' \link{fit.ns} for two-plane surveys of whale populations. There is
#' partial sibling information in that observations made by the same
#' plane cannot be siblings.
#'
#' @return A list, a suitable object for the \code{siblings} argument
#' to \link{fit.ns}.
#'
#' @param plane.id A vector indicating which plane collected each
#' observed location.
#'
#' @export
siblings.twoplane <- function(plane.id){
    n.points <- length(plane.id)
    out <- matrix(NA, nrow = n.points, ncol = n.points)
    for (i in 1:n.points){
        for (j in i:n.points){
            if (plane.id[i] == plane.id[j]){
                out[i, j] <- FALSE
            }
        }
    }
    list(matrix = out, pT = 0, pF = 0.5)
}

## Internal alternative for R6 function.
siblings.twoplane_r6<- function(plane.id){
    n.points <- length(plane.id)
    out <- matrix(NA, nrow = n.points, ncol = n.points)
    for (i in 1:n.points){
        for (j in i:n.points){
            if (plane.id[i] == plane.id[j]){
                out[i, j] <- FALSE
            }
        }
    }
    list(sibling.mat = out, alpha = 0, beta = 0.5)
}
