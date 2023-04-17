#' Clockwise rotation
#' 
#' This function calculates the clockwise rotation of angle \eqn{\theta} of a
#' given vector \eqn{v} in the plane.
#' 
#' 
#' @param v Vector \eqn{v} in the plane.
#' @param theta Angle \eqn{\theta} (in radians).
#' @return \item{v.rot}{Vector after rotation.}
#' @keywords nonparametric
#' @examples
#' 
#' \dontrun{
#' # Rotation of angle pi/4 of the vector (0,1)
#' rotation(v = c(0, 1), theta = pi/4)
#' }
#' 
#' @export rotation
rotation <-
function (v, theta) 
{
    v.rot <- numeric(2)
    v.rot[1] <- cos(theta) * v[1] + sin(theta) * v[2]
    v.rot[2] <- -sin(theta) * v[1] + cos(theta) * v[2]
    return(v.rot)
}
