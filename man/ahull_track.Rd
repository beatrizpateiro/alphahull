\name{ahull_track}
\alias{ahull_track}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{alpha-convex hull calculation of tracking data}
\description{
This function approximates the \eqn{\alpha}-convex hull of tracking data and returns a list of geom_path objects of the boundary.
}
\usage{
ahull_track(x, y = NULL, alpha, nps = 10000, sc = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{The \code{x} and \code{y} arguments provide the \code{x} and \code{y} coordinates of a set of points. Alternatively, a single argument \code{x} can be provided, see Details.}
  \item{alpha}{Value of \eqn{\alpha}.}
  \item{nps}{Number of points to generate in each segment connecting two locations, see Details}
  \item{sc}{Scale factor.}
}
\details{

An attempt is made to interpret the arguments x and y in a way suitable for computing the \eqn{\alpha}-convex hull. Any reasonable way of defining the coordinates is acceptable, see \code{\link{xy.coords}}. 

Increase \code{nps} if the trajectory is not contained in the computed estimator.
}
\value{
A list of geom_path objects defining the boundary of the \eqn{\alpha}-convex
}
\references{Cholaquidis, A., Fraiman, R., Lugosi, G. and Pateiro-Lopez, B. (2014) Set estimation from reflected Brownian motion. \emph{arXiv:1411.0433}.

Wikelski, M., and Kays, R. (2014). Movebank: archive, analysis and sharing of animal movement data. World Wide Web electronic publication.
}
\examples{
\dontrun{
# Data from Movebank
# Study Name: Missouri Bison Tracking Project
# Principal Investigator: Stephen Blake
library(move)
# Login to Movebank (requires a Movebank account)
curl <- movebankLogin(username = "xxx", password = "zzz")
# Download study data (requires agreed to the license terms in the web interface)
study <- "Missouri Bison Tracking Project"
animalName <-"AHP2059"
track <- getMovebankData(study = study, animalName = animalName, login = curl)
dat <- track@data
# Plot animal locations
p <- ggplot2::ggplot(
  dat,
  ggplot2::aes(x = location_long, y = location_lat)
) +
  ggplot2::geom_path(
    colour = "red",
    linewidth = 0.3
  ) +
  ggplot2::coord_quickmap()
p
# Alpha-convex hull of the animal trajectory
ah_out <- ahull_track(x = dat[, c("location_long", "location_lat")], alpha = 0.005)
# Plot the trajectory and the alpha-convex hull
p + ah_out
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
