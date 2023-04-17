#' alpha-shape calculation
#'
#' This function calculates the \eqn{\alpha}-shape of a given sample for
#' \eqn{\alpha>0}.
#'
#'
#' An attempt is made to interpret the arguments x and y in a way suitable for
#' computing the \eqn{\alpha}-shape, see \code{\link{xy.coords}}.
#'
#' The \eqn{\alpha}-shape is defined for any finite number of points. However,
#' since the algorithm is based on the Delaunay triangulation, at least three
#' non-collinear points are required.
#'
#' If \code{y} is NULL and \code{x} is an object of class \code{"delvor"}, then
#' the \eqn{\alpha}-shape is computed without invoking again the function
#' \code{\link{delvor}} (it reduces the computational cost).
#'
#' The function \code{\link{ashape}} returns (among other values) the matrix
#' \code{edges}. The structure of \code{edges} is that of matrix \code{mesh}
#' returned by the function \code{\link{delvor}}. Note that the
#' \eqn{\alpha}-shape is a subgraph of the Delaunay triangulation and,
#' therefore, \code{edges} is a submatrix of \code{mesh}.
#'
#' @param x,y The \code{x} and \code{y} coordinates of a set of points.
#' Alternatively, a single argument \code{x} can be provided, see Details.
#' @param alpha Value of \eqn{\alpha}.
#' @return A list with the following components: \item{edges}{A
#' \emph{n.seg}-row matrix with the coordinates and indexes of the edges of the
#' Delaunay triangulation that form the \eqn{\alpha}-shape. The number of rows
#' \emph{n.seg} coincides with the number of segments of the
#' \eqn{\alpha}-shape. The matrix also includes information of the Voronoi
#' extremes corresponding to each segment.} \item{length}{Length of the
#' \eqn{\alpha}-shape.} \item{alpha}{Value of \eqn{\alpha}.}
#' \item{alpha.extremes}{Vector with the indexes of the sample points that are
#' \eqn{\alpha}-extremes. See Edelsbrunnner \emph{et al.} (1983).}
#' \item{delvor.obj}{Object of class \code{"delvor"} returned by the
#' \code{\link{delvor}} function.} \item{x}{A 2-column matrix with the
#' coordinates of the set of points.}
#' @seealso \code{\link{plot.ashape}}, \code{\link{delvor}}.
#' @references Edelsbrunner, H., Kirkpatrick, D.G. and Seidel, R. (1983). On
#' the shape of a set of points in the plane. \emph{IEEE Transactions on
#' Information Theory}, 29(4), pp.551-559.
#' @keywords nonparametric
#' @examples
#'
#' \dontrun{
#' # Uniform sample of size n=300 in the annulus B(c,0.5)\B(c,0.25),
#' # with c=(0.5,0.5).
#' n <- 300
#' theta<-runif(n,0,2*pi)
#' r<-sqrt(runif(n,0.25^2,0.5^2))
#' x<-cbind(0.5+r*cos(theta),0.5+r*sin(theta))
#' # Value of alpha
#' alpha <- 0.1
#' # alpha-shape
#' ashape.obj <- ashape(x, alpha = alpha)
#' # If we change the value of alpha there is no need to compute
#' # again the Delaunay triangulation and Voronoi Diagram
#' alpha <- 0.4
#' ashape.obj.new <- ashape(ashape.obj$delvor.obj, alpha = alpha)
#'
#' # Random sample  from a uniform distribution on a Koch snowflake
#' # with initial side length 1 and 3 iterations
#' x <- rkoch(2000, side = 1, niter = 3)
#' # Value of alpha
#' alpha <- 0.05
#' # alpha-shape
#' ashape.obj <- ashape(x, alpha = alpha)
#' }
#'
#' @export ashape
ashape <-
function (x, y = NULL, alpha)
{
    if (alpha < 0) {
        stop("Parameter alpha must be greater or equal to zero")
    }
    if (!inherits(x, "delvor")) {
        dd.obj <- delvor(x, y)
    }
    else {
        dd.obj <- x
    }
    xy.data <- dd.obj$x
    mesh <- dd.obj$mesh
    dm1 <- sqrt((mesh[, "x1"] - mesh[, "mx1"])^2 + (mesh[, "y1"] -
        mesh[, "my1"])^2)
    dm2 <- sqrt((mesh[, "x1"] - mesh[, "mx2"])^2 + (mesh[, "y1"] -
        mesh[, "my2"])^2)
    dm1[mesh[, "bp1"] == 1] <- Inf
    dm2[mesh[, "bp2"] == 1] <- Inf
    n <- dim(xy.data)[1]
    ind <- 1:n
    ind.on <- grDevices::chull(xy.data)
    ind.in <- ind[-ind.on]
    n.on <- length(ind.on)
    n.in <- length(ind.in)
    if (n.in > 0) {
        aux <- rbind(cbind(mesh[, c("ind1", "ind2")], dm1), cbind(mesh[,
            c("ind1", "ind2")], dm2), cbind(mesh[, c("ind2",
            "ind1")], dm1), cbind(mesh[, c("ind2", "ind1")],
            dm2))
        fc <- factor(aux[, 1])
        aux <- tapply(aux[, 3], fc, max)
        alpha.max <- cbind(aux, as.numeric(levels(fc)))
        alpha.ext <- c(ind.on, ind.in[stats::na.omit(match(alpha.max[alpha <
            alpha.max[, 1], 2], ind.in))])
    }
    else {
        alpha.ext <- ind.on
    }
    n.edges <- dim(mesh)[1]
    is.edge <- numeric()
    ind.is <- 0
    i1 <- match(mesh[, 1], alpha.ext)
    i2 <- match(mesh[, 2], alpha.ext)
    is.edge <- which(i1 & i2)
    aux <- mesh[is.edge, ]
    n.pos <- dim(aux)[1]
    pm.x <- (aux[, "x1"] + aux[, "x2"]) * 0.5
    pm.y <- (aux[, "y1"] + aux[, "y2"]) * 0.5
    dm <- sqrt((aux[, "x1"] - aux[, "x2"])^2 + (aux[, "y1"] -
        aux[, "y2"])^2) * 0.5
    betw = rep(NA, n.pos)
    for (i in 1:n.pos) {
        if (aux[i, "mx1"] == aux[i, "mx2"]) {
            if (rank(c(aux[i, "my1"], aux[i, "my2"], pm.y[i]))[3] ==
                2) {
                betw[i] = 1
            }
        }
        else {
            if (rank(c(aux[i, "mx1"], aux[i, "mx2"], pm.x[i]))[3] ==
                2) {
                betw[i] = 1
            }
        }
    }
    l.min <- apply(cbind(dm1[is.edge], dm2[is.edge], dm * betw),
        1, min, na.rm = TRUE)
    l.max <- apply(cbind(dm1[is.edge], dm2[is.edge], dm * betw),
        1, max, na.rm = TRUE)
    in.ashape <- (l.min <= alpha & alpha <= l.max)
    edges <- matrix(t(aux[in.ashape, ]), byrow = TRUE, ncol = 12)
    colnames(edges) <- colnames(aux)
    ashape.obj <- list(edges = edges, length = sum(2 * dm[in.ashape]),
        alpha = alpha, alpha.extremes = alpha.ext, delvor.obj = dd.obj,
        x = xy.data)
    class(ashape.obj) <- "ashape"
    invisible(ashape.obj)
}
