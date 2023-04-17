#' Intersection of two circumferences
#' 
#' This function calculates the intersection of two circumferences, given their
#' centers and radius \eqn{c1,r1} and \eqn{c2,r2}, respectively.
#' 
#' The function \code{inter} is internally called by the function
#' \code{\link{ahull}}.
#' 
#' @param c11 \emph{X}-coordinate of the center \eqn{c1}.
#' @param c12 \emph{Y}-coordinate of the center \eqn{c1}.
#' @param r1 Radius \eqn{r1}.
#' @param c21 \emph{X}-coordinate of the center \eqn{c2}.
#' @param c22 \emph{Y}-coordinate of the center \eqn{c2}.
#' @param r2 Radius \eqn{r2}.
#' @return A list with the following components: \item{n.cut}{Number of
#' intersection points (0,1,2, or Inf).} \item{v1}{If there are two
#' intersection points, \code{v1} is the numeric vector whose components are
#' the coordinates of the unitary vector that has its origin in \eqn{c1} and
#' it's perpendicular to the chord that joins the intersection points of the
#' two circumferences. Otherwise, \code{v1=(0,0)}} \item{theta1}{Angle that
#' forms \code{v1} with the radius that joins the center \eqn{c1} with an
#' intersection point. } \item{v2}{If there are two intersection points,
#' \code{v2} is the numeric vector whose components are the coordinates of the
#' unitary vector that has its origin in \eqn{c2} and it's perpendicular to the
#' chord that joins the intersection points of the two circumferences.
#' Otherwise, \code{v2=(0,0)}} \item{theta2}{Angle that forms \code{v2} with
#' the radius that joins the center \eqn{c2} with an intersection point.}
#' @keywords nonparametric
#' @export inter
inter <-
function (c11, c12, r1, c21, c22, r2) 
{
    d <- sqrt((c21 - c11)^2 + (c22 - c12)^2)
    v.x <- 0
    v.y <- 0
    theta1 <- 0
    theta2 <- 0
    if (d == 0 & r1 == r2) {
        n.cut <- Inf
    }
    else if (d > (r1 + r2)) {
        n.cut <- 0
    }
    else {
        if (d < max(r1, r2) - min(r1, r2)) {
            n.cut <- 0
        }
        else if (d == max(r1, r2) - min(r1, r2)) {
            n.cut <- 1
        }
        else if (d == r1 + r2) {
            n.cut <- 1
        }
        else {
            n.cut <- 2
            d1 <- (d^2 - r2^2 + r1^2)/(2 * d)
            d2 <- d - d1
            v.x <- (c21 - c11)/d
            v.y <- (c22 - c12)/d
            theta1 <- acos(d1/r1)
            theta2 <- acos(d2/r2)
        }
    }
    return(list(n.cut = n.cut, v1 = c(v.x, v.y), theta1 = theta1, 
        v2 = c(-v.x, -v.y), theta2 = theta2))
}
