#' Delaunay triangulation and Voronoi diagram
#'
#' This function returns a matrix with information about the Delaunay
#' triangulation and Voronoi diagram of a given sample.
#'
#' An attempt is made to interpret the arguments x and y in a way suitable for
#' computing the Delaunay triangulation and Voronoi diagram . Any reasonable
#' way of defining the coordinates is acceptable, see \code{\link{xy.coords}}.
#'
#' The function \code{\link[interp]{tri.mesh}} from package \pkg{interp}
#' calculates the Delaunay triangulation of at least three non-collinear
#' points. Using the Delaunay triangulation, the function \code{delvor}
#' calculates the correspondig Voronoi diagram. For each edge of the Delaunay
#' triangulation there is a segment in the Voronoi diagram, given by the union
#' of the circumcenters of the two neighbour triangles that share the edge. For
#' those triangles with edges on the convex hull, the corresponding line in the
#' Voronoi diagram is a semi-infinite segment, whose boundless extreme is
#' calculated by the function \code{\link{dummycoor}}. The function
#' \code{delvor} returns the sample, the output object of class \code{"triSht"}
#' from the function \code{\link[interp]{tri.mesh}} and a matrix \code{mesh}
#' with all the necessary information of the Delaunay triangulation and Voronoi
#' diagram. Thus, for each edge of the Delaunay triangulation the output matrix
#' contains the indexes and coordinates of the sample points that form the
#' edge, the indexes and coordinates of the extremes of the corresponding
#' segment in the Voronoi diagram, and an indicator that takes the value 1 for
#' those extremes of the Voronoi diagram that represent a boundless extreme.
#'
#' @param x,y The \code{x} and \code{y} arguments provide the \code{x} and
#' \code{y} coordinates of a set of points. Alternatively, a single argument
#' \code{x} can be provided, see Details.
#' @return A list with the following components: \item{mesh}{A
#' \eqn{n.edges}-row matrix, where \eqn{n.edges} is the total number of
#' different edges of the Delaunay triangulation.} \item{x}{A 2-column matrix
#' with the coordinates of the sample points.} \item{tri.obj}{Object of class
#' \code{"tri"}. See \code{\link[interp]{tri.mesh}} in package \pkg{interp}.}
#' @seealso \code{\link{plot.delvor}}.
#' @references Renka, R. J. (1996). Algorithm 751: TRIPACK: a constrained
#' two-dimensional Delaunay triangulation package, \emph{ACM Trans. Math.
#' Softw.}, 22(1), pp.1-8.
#' @keywords nonparametric
#' @examples
#'
#' \dontrun{
#' # Random sample in the unit square
#' x <- matrix(runif(20), nc = 2)
#' # Delaunay triangulation and Voronoi diagram calculation
#' delvor.obj <- delvor(x)
#' }
#'
#' @export delvor
delvor <-
function (x, y = NULL)
{
    X <- grDevices::xy.coords(x, y)
    x <- cbind(X$x, X$y)
    if (dim(x)[1] <= 2) {
        stop("At least three non-collinear points are required")
    }
    tri.obj <- interp::tri.mesh(X)
    tri <- interp::triangles(tri.obj)
    nt <- nrow(tri)
    circenter <- matrix(nrow = nt, ncol = 2)
    colnames(circenter) <- c("circumx", "circumy")
    for (i in 1:nt){
		aux <- interp::circum(c(x[tri[i, 1], 1], x[tri[i, 2], 1], x[tri[i, 3], 1]), c(x[tri[i, 1], 2], x[tri[i, 2], 2], x[tri[i, 3], 2]))
		circenter[i, ] <- c(aux$x, aux$y)
    }
    tri.info<-cbind(tri, circenter)

    n.tri <- dim(tri.info)[1]
    n.arc <- max(tri.info[, 7:9])
    if (n.tri == 1) {
        aux1 <- cbind(matrix(tri.info[, c("arc1", "node2", "node3")],
            ncol = 3, nrow = 1), 1:n.tri, tri.info[, "tr1"])
        aux2 <- cbind(matrix(tri.info[, c("arc2", "node1", "node3")],
            ncol = 3, nrow = 1), 1:n.tri, tri.info[, "tr2"])
        aux3 <- cbind(matrix(tri.info[, c("arc3", "node1", "node2")],
            ncol = 3, nrow = 1), 1:n.tri, tri.info[, "tr3"])
    }
    else {
        aux1 <- cbind(tri.info[, c("arc1", "node2", "node3")],
            1:n.tri, tri.info[, "tr1"])
        aux2 <- cbind(tri.info[, c("arc2", "node1", "node3")],
            1:n.tri, tri.info[, "tr2"])
        aux3 <- cbind(tri.info[, c("arc3", "node1", "node2")],
            1:n.tri, tri.info[, "tr3"])
    }
    aux <- rbind(aux1, aux2, aux3)
    repeated <- duplicated(aux[, 1])
    aux <- aux[!repeated, ]
    colnames(aux) <- c("arc", "ind1", "ind2", "indm1", "indm2")
    bp1 <- (aux[, "indm1"] == 0)
    bp2 <- (aux[, "indm2"] == 0)
    is.dummy <- which(bp2)
    n.dummy <- length(is.dummy)
    circumcentres <- tri.info[, c("circumx", "circumy")]
    away <- max(diff(range(x[, 1])), diff(range(x[, 2])))
    for (i in is.dummy) {
        n.tri <- n.tri + 1
        dum <- dummycoor(tri.obj, x[aux[i, "ind1"], ], x[aux[i,
            "ind2"], ], tri.info[aux[i, "indm1"], c("circumx",
            "circumy")], away)
        circumcentres <- rbind(circumcentres, dum)
        aux[i, "indm2"] <- n.tri
    }
    mesh <- cbind(aux[, c("ind1", "ind2")], x[aux[, "ind1"],
        ], x[aux[, "ind2"], ], circumcentres[aux[, "indm1"],
        ], circumcentres[aux[, "indm2"], ], bp1, bp2)
    colnames(mesh) <- c("ind1", "ind2", "x1", "y1", "x2", "y2",
        "mx1", "my1", "mx2", "my2", "bp1", "bp2")
    delvor.obj <- list(mesh = mesh, x = x, tri.obj = tri.obj)
    class(delvor.obj) <- "delvor"
    invisible(delvor.obj)
}
