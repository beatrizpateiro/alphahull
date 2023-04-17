#' Add an arc to a plot
#'
#' This function adds the arc of \eqn{B(c,r)} between the angles that
#' \eqn{A_\theta v} and \eqn{A_{-\theta} v} form with the axis \emph{OX}, where
#' \eqn{A_\theta v} represents the clockwise rotation of angle \eqn{\theta} of
#' the vector \eqn{v}.
#'
#'
#' @param c Center \eqn{c} of the arc.
#' @param r Radius \eqn{r} of the arc.
#' @param v Vector \eqn{v} in the plane.
#' @param theta Angle \eqn{\theta} (in radians).
#' @param \dots Arguments to be passed to methods, such as graphical parameters
#' (see \code{\link{par}}).
#' @seealso \code{\link{plot.ahull}}.
#' @keywords nonparametric
#' @examples
#'
#' \dontrun{
#' # Plot of the circumference of radius 1
#' theta <- seq(0, 2*pi, length = 100)
#' r <- 1
#' plot(r*cos(theta), r*sin(theta), type = "l")
#' # Add in red the arc between pi/4 and 3*pi/4
#' arc(c(0,0), 1, c(0,1), pi/4, col = 2, lwd = 2)
#' }
#'
#' @export arc
arc <-
function (c, r, v, theta, ...)
{
    angles <- anglesArc(v, theta)
    seqang <- seq(angles[1], angles[2], length = 100)
    x <- c[1] + r * cos(seqang)
    y <- c[2] + r * sin(seqang)
    graphics::lines(x, y, ...)
}
