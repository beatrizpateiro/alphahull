\name{delvor}
\alias{delvor}
\title{Delaunay triangulation and Voronoi diagram}
\description{
   This function returns a list containing the Delaunay triangulation
   and Voronoi diagram of a given planar sample.
}
\usage{
delvor(x, y = NULL)
}
\arguments{
   \item{x, y}{The \code{x} and \code{y} arguments provide the coordinates
   of a set of planar points. Alternatively, a single argument \code{x}
   can be provided, see Details.}
}
\details{
The coordinates are obtained using \code{\link[grDevices]{xy.coords}}.
At least three non-collinear points are required.

The Delaunay triangulation is computed with
\code{\link[geometry]{delaunayn}} from the \pkg{geometry} package,
which uses Qhull. The corresponding Voronoi diagram is reconstructed
from the triangle incidences and circumcentres. For an interior Delaunay
edge, its dual Voronoi edge joins the circumcentres of the two adjacent
triangles.

For a Delaunay edge on the convex hull, the dual Voronoi edge is
semi-infinite. The function \code{\link{dummycoor}} supplies fictitious
coordinates along the outward-pointing ray for plotting; these are not
additional Voronoi vertices. Its direction is determined using the
sample centroid.

The matrix \code{mesh} contains one row for each Delaunay edge:
\code{ind1}, \code{ind2} are its sample-point indices;
\code{x1}, \code{y1}, \code{x2}, \code{y2} are its endpoint coordinates;
\code{mx1}, \code{my1}, \code{mx2}, \code{my2} are the coordinates of
the corresponding Voronoi edge endpoints; and \code{bp1}, \code{bp2}
identify fictitious endpoints of unbounded Voronoi rays.
}
\value{
   A list of class \code{"delvor"} with the following components:
   \item{mesh}{A twelve-column matrix with one row for each distinct
   Delaunay edge and its dual Voronoi edge.}
   \item{x}{A two-column matrix with the coordinates of the sample points.}
   \item{tri.obj}{An integer matrix with three columns, containing the
   indices of the three sample points forming each Delaunay triangle.}
}
\references{
Barber, C. B., Dobkin, D. P. and Huhdanpaa, H. (1996).
The Quickhull algorithm for convex hulls.
\emph{ACM Transactions on Mathematical Software}, 22(4), pp. 469--483.
}
\seealso{\code{\link{plot.delvor}}.}
\examples{
\dontrun{
# Random sample in the unit square
x <- matrix(runif(20), ncol = 2)
# Delaunay triangulation and Voronoi diagram calculation
delvor.obj <- delvor(x)
}
}
\keyword{nonparametric}
