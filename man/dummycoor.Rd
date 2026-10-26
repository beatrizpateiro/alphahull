\name{dummycoor}
\alias{dummycoor}
\title{Fictitious endpoint of a semi-infinite Voronoi edge}
\description{
   Computes a fictitious endpoint in the outward direction of a
   semi-infinite Voronoi edge, for plotting purposes.
}
\usage{
dummycoor(tri.obj, l1, l2, m, away)
}
\arguments{
  \item{tri.obj}{Sample coordinates as a two-column matrix or as a list
    with numeric components \code{x} and \code{y}.}
  \item{l1}{Numeric vector of length two containing the coordinates of
    the first endpoint of a Delaunay edge on the convex hull.}
  \item{l2}{Numeric vector of length two containing the coordinates of
    the second endpoint of that Delaunay edge.}
  \item{m}{Numeric vector of length two containing the circumcentre
    of the adjacent Delaunay triangle.}
  \item{away}{Positive numeric distance used to place the fictitious
    endpoint along the outward ray.}
}
\details{
   A Delaunay edge on the convex hull corresponds to a semi-infinite
   Voronoi edge whose finite endpoint is the circumcentre \code{m}.
   The direction is perpendicular to the Delaunay edge joining
   \code{l1} and \code{l2}. The outward orientation is selected using
   the centroid of the sample coordinates in \code{tri.obj}.

   The fictitious endpoint is placed along the ray at a distance
   depending on \code{away} and the projection of the Delaunay edge
   midpoint. This ensures that the plotting endpoint extends beyond
   that projection. The returned point is not an actual Voronoi vertex.
}
\value{
   A numeric vector of length two giving the fictitious endpoint
   of the semi-infinite Voronoi ray.
}
\seealso{ \code{\link{delvor}}.}
\examples{
xy <- rbind(c(0, 0), c(1, 0), c(0, 1))
dummycoor(xy, c(0, 0), c(1, 0), c(0.5, 0.5), away = 2)
}
\keyword{nonparametric}
