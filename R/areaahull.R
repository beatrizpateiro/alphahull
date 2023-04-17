#' Area of the alpha-convex hull
#'
#' This function calculates the area of the \eqn{\alpha}-convex hull of a
#' sample of points.
#'
#'
#' @param x Object of class \code{"ahull"}.
#' @param timeout A numeric specifying the maximum number of seconds the
#' expression is allowed to run before being interrupted by the timeout.
#' @return \item{area}{Area of the \eqn{\alpha}-convex hull. If the area cannot
#' be computed, the output will be NA with a warning.}
#' @seealso \code{\link{ahull}}.
#' @keywords nonparametric
#' @examples
#'
#' \dontrun{
#' # Random sample in the unit square
#' x <- matrix(runif(500), nc = 2)
#' # Value of alpha
#' alpha <- 1
#' # alpha-convex hull
#' ahull.obj <- ahull(x, alpha = alpha)
#' # Area of the alpha-convex hull
#' areaahull(ahull.obj)
#' }
#'
#' @export areaahull
areaahull <-
function (x, timeout = 5)
{
    area <- R.utils::withTimeout(try(areaahulleval(x), silent = TRUE),
        timeout = timeout)
    if (!is.numeric(area)) {
        warning("Problem in area computation (Returns NA)")
        area <- NA
    }
    if (is.numeric(area) & area < 0) {
        warning("Problem in area computation (Returns NA)")
        area <- NA
    }
    return(area)
}
