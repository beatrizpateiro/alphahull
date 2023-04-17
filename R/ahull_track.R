#' alpha-convex hull calculation of tracking data
#'
#' This function approximates the \eqn{\alpha}-convex hull of tracking data and
#' returns a list of geom_path objects of the boundary.
#'
#'
#' An attempt is made to interpret the arguments x and y in a way suitable for
#' computing the \eqn{\alpha}-convex hull. Any reasonable way of defining the
#' coordinates is acceptable, see \code{\link{xy.coords}}.
#'
#' Increase \code{nps} if the trajectory is not contained in the computed
#' estimator.
#'
#' @param x,y The \code{x} and \code{y} arguments provide the \code{x} and
#' \code{y} coordinates of a set of points. Alternatively, a single argument
#' \code{x} can be provided, see Details.
#' @param alpha Value of \eqn{\alpha}.
#' @param nps Number of points to generate in each segment connecting two
#' locations, see Details
#' @param sc Scale factor.
#' @return A list of geom_path objects defining the boundary of the
#' \eqn{\alpha}-convex
#' @references Cholaquidis, A., Fraiman, R., Lugosi, G. and Pateiro-Lopez, B.
#' (2014) Set estimation from reflected Brownian motion.
#' \emph{arXiv:1411.0433}.
#'
#' Wikelski, M., and Kays, R. (2014). Movebank: archive, analysis and sharing
#' of animal movement data. World Wide Web electronic publication.
#' @keywords nonparametric
#' @examples
#'
#' \dontrun{
#' library(move)
#' library(ggmap)
#' # Data from Movebank
#' # Study Name: Dunn Ranch Bison Tracking Project
#' # Principal Investigator: Stephen Blake, Randy Arndt, Doug Ladd
#' # Max Planck Institute for Ornithology Radolfzell Germany
#' study <- "Dunn Ranch Bison Tracking Project"
#' cainfo <- system.file("CurlSSL", "cacert.pem", package = "RCurl")
#' options(RCurlOptions = list(verbose = FALSE, capath = cainfo, ssl.verifypeer = FALSE))
#' # Login to movebank (first create the login object)
#' curl <- movebankLogin(username = "xxx", password = "zzz")
#' # Downloads study stored in Movebank
#' track <- getMovebankData(study = study, login = curl)
#' dat <- track@data[track@data[, "deployment_id"] == 13848432,]
#' # Map of animal locations
#' bbox <- ggmap::make_bbox(dat[,"location_long"], dat[,"location_lat"], f = 0.3)
#' map_loc <- get_map(location = bbox, source = "google", maptype = 'satellite')
#' map <- ggmap(map_loc, extent = 'panel', maprange=FALSE)
#' p <- map + geom_path(data = dat, aes(x = location_long, y = location_lat), col=2, size=0.3)
#' p
#' ah_gp <- ahull_track(x = dat[, c("location_long", "location_lat")], alpha = 0.005)
#' p + ah_gp
#' }
#'
#' @export ahull_track
#' @importFrom ggplot2 geom_path aes
ahull_track <-
function (x, y = NULL, alpha, nps = 10000, sc = 100)
{
    X <- grDevices::xy.coords(x, y)
    dat <- data.frame(location_long = X$x, location_lat = X$y)
    colnames(dat) <- c("location_long", "location_lat")
    np <- dim(dat)[1]
    seg <- spatstat.geom::psp(dat[-np, 1], dat[-np, 2], dat[-1, 1], dat[-1,
        2], window = spatstat.geom::owin(range(dat[, 1]), range(dat[, 2])))
    pseg <- spatstat.random::runifpointOnLines(nps, seg)
    xah <- c(dat[, 1], pseg$x)
    yah <- c(dat[, 2], pseg$y)
    xah <- jitter(xah) * sc
    yah <- jitter(yah) * sc
    ah <- ahull(xah, yah, alpha = alpha * sc)
    aux <- ah$arcs[ah$arcs[, 3] > 0, ]
    na <- dim(aux)[1]
    if (na >= 1) {
        dfr <- list()
        for (j in 1:na) {
            c <- aux[j, 1:2]
            r <- aux[j, 3]
            v <- aux[j, 4:5]
            theta <- aux[j, 6]
            angles <- anglesArc(v, theta)
            seqang <- seq(angles[1], angles[2], length = 100)
            x <- c[1] + r * cos(seqang)
            y <- c[2] + r * sin(seqang)
            daux <- data.frame(x = x/sc, y = y/sc)
            dfr[[j]] <- ggplot2::geom_path(data = daux, aes(x = x, y = y),
                col = 4)
        }
        return(dfr)
    }
    else {
        warning("Please, choose a larger value of alpha")
    }
}
