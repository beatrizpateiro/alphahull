\name{dw_track}
\alias{dw_track}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{RBM-sausage calculation of tracking data}
\description{
This function approximates the RBM-sausage of tracking data and returns a list of geom_path objects of the boundary.
}
\usage{
dw_track(x, y = NULL, eps, nps = 20000, sc = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{The \code{x} and \code{y} arguments provide the \code{x} and \code{y} coordinates of a set of points. Alternatively, a single argument \code{x} can be provided, see Details.}
  \item{eps}{Value of \eqn{\epsilon}.}
  \item{nps}{Number of points to generate in each segment connecting two locations, see Details.}
  \item{sc}{Scale factor.}
}
\details{

An attempt is made to interpret the arguments x and y in a way suitable for computing the RBM-sausage. Any reasonable way of defining the coordinates is acceptable, see \code{\link{xy.coords}}. 
Given a trajectory in the plane, the estimator is defined as the set of points whose distance to the trajectory is less than or equal to \eqn{\epsilon} (this estimator is analogous to the one of Devroye and Wise (1980) for i.i.d. data). Increase \code{nps} if the trajectory is not contained in the computed estimator.
}
\value{
A list of geom_path objects defining the boundary of the estimator
}
\references{Cholaquidis, A., Fraiman, R., Lugosi, G. and Pateiro-Lopez, B. (2014) Set estimation from reflected Brownian motion. \emph{arXiv:1411.0433}.

Devroye, L. and Wise, G. (1980) Detection of abnormal behaviour via nonparametric estimation of the support. \emph{SIAM J. Appl. Math.} 3, pp. 480-488.

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
# Devroye-Wise of the animal trajectory
dw_out <- dw_track(x = dat[, c("location_long", "location_lat")], eps = 0.001)
# Plot the trajectory and the Devroye-Wise
p + dw_out
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
