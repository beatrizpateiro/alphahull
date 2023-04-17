#' Length of the boundary of the alpha-convex hull
#' 
#' This function calculates the length of the boundary of the
#' \eqn{\alpha}-convex hull of a given sample.
#' 
#' The function \code{lengthahull} is internally called by the function
#' \code{\link{ahull}}.
#' 
#' @param ahull.arcs Output matrix of arcs returned by \code{\link{ahull}}.
#' @return \item{length}{Length of the boundary of the \eqn{\alpha}-convex
#' hull.}
#' @seealso \code{\link{ahull}}.
#' @keywords nonparametric
#' @examples
#' 
#' \dontrun{
#' # Random sample in the unit square
#' x <- matrix(runif(100), nc = 2)
#' # Value of alpha
#' alpha <- 0.2
#' # alpha-convex hull
#' ahull.obj <- ahull(x, alpha = alpha)
#' # Length of the alpha-convex hull
#' ahull.obj$length
#' }
#' 
#' @export lengthahull
lengthahull <-
function (ahull.arcs) 
{
    gamma <- 2 * ahull.arcs[, "theta"]
    length <- gamma * ahull.arcs[, "r"]
    return(length = sum(length))
}
