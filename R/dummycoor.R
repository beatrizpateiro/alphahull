#' Semi-infinite edge of the Voronoi diagram
#'
#' This function determines fictitious coordinates for the boundless extreme of
#' a semi-infinite edge of the Voronoi diagram.
#'
#' When a triangle of the Delaunay triangulation has one of its edges (given by
#' the segment that joins the sample points with indexes \code{l1} and
#' \code{l2}) on the convex hull, the corresponding segment of the Voronoi
#' diagram is semi-infinite. The finite extreme coincides with the circumcenter
#' of the triangle and the direction of the line is given by the perpendicular
#' bisector of the edge that lies on the convex hull.
#'
#' @param tri.obj Object of class \code{"triSht"}. See
#' \code{\link[interp]{tri.mesh}} in package \pkg{interp}.
#' @param l1 Index of the sample point correponding to one vertex of a triangle
#' of Delaunay that lies on the convex hull, see Details.
#' @param l2 Index of the sample point correponding to other vertex of a
#' triangle of Delaunay that lies on the convex hull, see Details.
#' @param m Index of the circumcenter of the triangle of Delaunay with one edge
#' on the convex hull.
#' @param away Constant that determines how far away the fictitious boundless
#' extreme is located.
#' @return \item{dum}{Fictitious coordinates of the boundless extreme.}
#' @seealso \code{\link{delvor}}.
#' @keywords nonparametric
#' @export dummycoor
dummycoor <-
function (tri.obj, l1, l2, m, away)
{
    v <- l2 - l1
    v <- c(v[2], -v[1])
    norm <- sum(v^2)
    if (norm > 0) {
        v <- v/norm
    }
    mp <- (l1 + l2)/2
    eps <- 1e-05
    test <- mp + eps * v
    inconv <- interp::in.convex.hull(tri.obj, test[1], test[2])
    if (inconv) {
        dum <- m - away * v
    }
    else {
        dum <- m + away * v
    }
    return(dum)
}
