#' Angles of the extremes of an arc
#' 
#' Given a vector \eqn{v} and an angle \eqn{\theta}, \code{anglesArc} returns
#' the angles that \eqn{A_\theta v} and \eqn{A_{-\theta} v} form with the axis
#' \emph{OX}, where \eqn{A_\theta v} represents the clockwise rotation of angle
#' \eqn{\theta} of the vector \eqn{v}.
#' 
#' The angle that forms the vector \eqn{v} with the axis \emph{OX} takes its
#' value in \eqn{[0,2\pi)}.
#' 
#' @param v Vector \eqn{v} in the plane.
#' @param theta Angle \eqn{\theta} (in radians).
#' @return \item{angs}{Numeric vector with two components.}
#' @keywords nonparametric
#' @examples
#' 
#' \dontrun{
#' # Let v=c(0,1) and theta=pi/4
#' # Consider the arc such that v is the internal angle bisector that 
#' # divides the angle 2*theta into two equal angles
#' # The angles that the arc forms with the OX axis are pi/4 and 3*pi/4
#' v <- c(0,1)
#' theta <- pi/4
#' anglesArc(v,theta)
#' }
#' 
#' @export anglesArc
anglesArc <-
function (v, theta) 
{
    theta.OX <- ifelse(v[2] >= 0, acos(v[1]), 2 * pi - acos(v[1]))
    angs <- c(theta.OX - theta, theta.OX + theta)
    return(angs)
}
