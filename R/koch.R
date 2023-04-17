#' Construct a Kock snowflake curve
#' 
#' This function uses recursion to construct a Kock snowflake curve.
#' 
#' The Koch snowflake is a fractal curve described by the Swedish mathematician
#' Helge von Koch in 1904. It is built by starting with an equilateral
#' triangle, removing the inner third of each side, building another
#' equilateral triangle at the location where the side was removed, and then
#' repeating the process.
#' 
#' @param side Side length of the initial equilateral triangle.
#' @param niter Number of iterations in the development of the snowflake curve.
#' @return \item{vertices}{A 2-column matrix with the coordinates of the
#' snowflake vertices.}
#' @seealso \code{\link{rkoch}}.
#' @references von Koch, H. (1904). Sur une courbe continue sans tangente,
#' obtenue par une construction geometrique elementaire. \emph{Arkiv for
#' Matematik}, 1, pp.681-704.
#' @keywords nonparametric
#' @examples
#' 
#' \dontrun{
#' # The first four iterations of a Koch snowflake 
#' # with side length of the initial equilateral triangle equal to 3.
#' vertices <- koch(side = 2, niter = 4)
#' plot(vertices[, 1], vertices[, 2], type = "l", asp = TRUE, 
#' main = "Koch snowflake", xlab = "", ylab = "", col = 4)
#' polygon(vertices[, 1], vertices[, 2] , col = 4)
#' }
#' 
#' @export koch
koch <-
function (side = 3, niter = 5) 
{
    npoints = 3 * 4^(niter - 1)
    p <- matrix(nrow = npoints, ncol = 2)
    index <- c(1, npoints/3 + 1, 2 * npoints/3 + 1)
    p[1, ] <- c(-side/2, 0)
    p[index[3], ] <- c(side/2, 0)
    p[index[2], ] <- p[1, ] + rotation(p[index[3], ] - p[1, ], 
        -pi/3)
    if (niter >= 2) {
        for (k in 2:niter) {
            npoints2 = 3 * 4^(niter - k)
            index1 <- index + npoints2/3
            index2 <- index + npoints2
            index3 <- index1 + npoints2/3
            aux <- diff(rbind(p[!is.na(p[, 1]), ], p[1, ]))/3
            p[index1, ] <- p[index, ] + aux
            p[index2, ] <- p[index, ] + 2 * aux
            aux1 <- matrix(c(cos(-pi/3), sin(-pi/3)), nrow = length(index1), 
                ncol = 2, byrow = TRUE)
            aux2 <- matrix(c(-sin(-pi/3), cos(-pi/3)), nrow = length(index1), 
                ncol = 2, byrow = TRUE)
            v <- p[index2, ] - p[index1, ]
            p[index3, ] <- p[index1, ] + cbind(apply(aux1 * v, 
                1, sum), apply(aux2 * v, 1, sum))
            index <- sort(c(index, index1, index2, index3))
        }
    }
    return(vertices = p)
}
