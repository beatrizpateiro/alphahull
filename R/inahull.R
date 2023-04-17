#' Determines for one or more points whether they belong to the alpha-convex
#' hull
#' 
#' This function determines for one or more points \eqn{p} whether they belong
#' to the \eqn{\alpha}-convex hull of a sample.
#' 
#' The complement of the \eqn{\alpha}-convex hull of a sample is calculated by
#' \code{\link{complement}}. The function \code{\link{inahull}} checks whether
#' each point in \eqn{p} belongs to any of the open balls or halfplanes that
#' define the complement.
#' 
#' @param ahull.obj Object of class \code{"ahull"} returned by the funcion
#' \code{\link{ahull}}.
#' @param p Numeric vector with two components describing a point in the plane
#' or two-column matrix of points.
#' @return \item{in.ahull}{A logical vector specifying whether each point in
#' \eqn{p} belongs to the \eqn{\alpha}-convex hull.}
#' @seealso \code{\link{ahull}}, \code{\link{complement}}.
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
#' # Check if the point (0.5, 0.5) belongs to the alpha-convex hull
#' inahull(ahull.obj, p = c(0.5, 0.5))
#' # Check if the points (0.5, 0.5) and (2, 2) belong to the alpha-convex hull
#' inahull(ahull.obj, p = rbind(c(0.5, 0.5), c(2, 2)))
#' }
#' 
#' @export inahull
inahull <-
function (ahull.obj, p) 
{
    p <- matrix(p, ncol = 2)
    compl <- ahull.obj$complement
    halfpl <- which(compl[, "r"] < 0)
    n.halfpl <- length(halfpl)
    ball <- which(compl[, "r"] > 0)
    n.ball <- length(ball)
    in.compl <- logical(length = dim(p)[1])
    if (n.halfpl >= 1) {
        h <- 1
        while ((h <= n.halfpl)) {
            sig = compl[halfpl[h], 3]
            a = compl[halfpl[h], 1]
            b = compl[halfpl[h], 2]
            if (sig <= -3) {
                in.compl[p[, 1] > a & sig == 3] <- TRUE
                in.compl[p[, 1] < a & sig == -4] <- TRUE
            }
            else {
                in.compl[p[, 2] > (a + b * p[, 1]) & sig == -1] <- TRUE
                in.compl[p[, 2] < (a + b * p[, 1]) & sig == -2] <- TRUE
            }
            h <- h + 1
        }
    }
    wp <- which(!in.compl)
    k <- 1
    while ((k <= n.ball)) {
        r = compl[ball[k], 3]
        c1 = compl[ball[k], 1]
        c2 = compl[ball[k], 2]
        d <- sqrt((p[wp, 1] - c1)^2 + (p[wp, 2] - c2)^2)
        in.compl[wp[d < r]] <- TRUE
        k <- k + 1
    }
    return(in.ahull = !in.compl)
}
