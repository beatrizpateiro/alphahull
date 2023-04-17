#' Random generation on a Koch snowflake curve
#'
#' This function generates ramdom points from a uniform distribution on a Koch
#' snowflake.
#'
#'
#' @param n Number of observations.
#' @param side Side length of the initial equilateral triangle of the Koch
#' snowflake curve.
#' @param niter Number of iterations in the development of the snowflake curve.
#' @return A 2-column matrix with the coordinates of generated points.
#' @seealso \code{\link{koch}}.
#' @keywords nonparametric
#' @examples
#'
#' \dontrun{
#' unifkoch <- rkoch(2000, side = 1, niter = 3)
#' plot(unifkoch, asp = TRUE)
#' }
#'
#' @export rkoch
rkoch <-
function (n, side = 3, niter = 5)
{
    vert <- koch(side, niter)
    square <- c(min(vert[, 1]), max(vert[, 1]), min(vert[, 2]),
        max(vert[, 2]))
    m = 0
    xy <- matrix(nrow = n, ncol = 2)
    while (m < n) {
        x <- stats::runif(1, square[1], square[2])
        y <- stats::runif(1, square[3], square[4])
        if (sgeostat::in.polygon(x, y, vert[, 1], vert[, 2])) {
            xy[m + 1, ] <- c(x, y)
            m <- m + 1
        }
    }
    return(xy)
}
