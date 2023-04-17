#' Plot the alpha-convex hull
#'
#' This function returns a plot of the \eqn{\alpha}-convex hull. If desired, it
#' also adds the Delaunay triangulation, Voronoi diagram and \eqn{\alpha}-shape
#' of the sample.
#'
#'
#' @param x Object of class \code{"ahull"}.
#' @param add Logical, if TRUE add to a current plot.
#' @param do.shape Logical, indicates if the \eqn{\alpha}-shape should also be
#' plotted.
#' @param wlines "Which lines?". I.e. should the Delaunay triangulation be
#' plotted (wlines='del'), should the Voronoi diagram be plotted
#' (wlines='vor'), should both be plotted (wlines='both'), or none
#' (wlines='none', the default)?
#' @param wpoints Logical, indicates if sample points should be plotted.
#' @param number Logical, defaulting to FALSE; if TRUE then the points plotted
#' will be labelled with their index numbers.
#' @param col The colour numbers for plotting the \eqn{\alpha}-convex hull,
#' \eqn{\alpha}-shape, data points, Delaunay triangulation, Voronoi diagram,
#' and the point numbers, in that order; defaults to c(1,1,1,1,1,1). If fewer
#' than six numbers are given, they are recycled. (If more than six numbers are
#' given, the redundant ones are ignored.)
#' @param xlim The limits on the x-axis.
#' @param ylim The limits on the y-axis.
#' @param lwd The line widths for plotting the tesselations, the
#' \eqn{\alpha}-shape, and the \eqn{\alpha}-convex hull, in that order;
#' defaults to c(1,1,2).
#' @param \dots Arguments to be passed to methods, such as graphical parameters
#' (see \code{\link{par}}).
#' @seealso \code{\link{ahull}}, \code{\link{ashape}}.
#' @keywords nonparametric
#' @examples
#'
#' \dontrun{
#' # Random sample in the unit square
#' x <- matrix(runif(100), nc = 2)
#' # Value of alpha
#' alpha <- 0.2
#' # alpha-convex hull
#' ahull.obj <- ahull(x, alpha = alpha)
#' # Plot including the alpha-convex hull in pink, alpha-shape in blue,
#' # sample points in black, voronoi diagram in green
#' # and Delaunay triangulation in red
#' plot(ahull.obj, do.shape = TRUE, wlines = "both", col = c(6, 4, 1, 2, 3))
#'
#' # Random sample  from a uniform distribution on a Koch snowflake
#' # with initial side length 1 and 3 iterations
#' x <- rkoch(2000, side = 1, niter = 3)
#' # Value of alpha
#' alpha <- 0.05
#' # Alpha-convex hull
#' ahull.obj <- ahull(x, alpha = alpha)
#' plot(ahull.obj)
#' }
#' @export
plot.ahull <-
function (x, add = FALSE, do.shape = FALSE, wlines = c("none",
    "both", "del", "vor"), wpoints = TRUE, number = FALSE, col = NULL,
    xlim = NULL, ylim = NULL, lwd = NULL, ...)
{
    wlines <- match.arg(wlines)
    if (is.null(class(x)) || !inherits(x, "ahull")) {
        cat("Argument is not of class ahull.\n")
        return(invisible())
    }
    if (is.null(col)) {
        col <- c(1, 1, 1, 1, 1, 1)
    }
    else {
        col <- rep(col, length.out = 6)
    }
    if (is.null(lwd)) {
        lwd <- c(1, 1, 2)
    }
    else {
        lwd <- rep(lwd, length.out = 3)
    }
    wlines <- match.arg(wlines)
    plot.dd <- switch(wlines, none = TRUE, both = FALSE, del = FALSE,
        vor = FALSE)
    if (do.shape) {
        plot.ashape(x$ashape.obj, add = add, wlines = wlines,
            wpoints = wpoints, number = number, col = col[2:6],
            xlim = xlim, ylim = ylim, lwd = lwd[1:2], ...)
    }
    else {
        if (plot.dd) {
            if (!add) {
                if (is.null(xlim))
                  xlim <- range(x$ashape.obj$x[, 1])
                if (is.null(ylim))
                  ylim <- range(x$ashape.obj$x[, 2])
                plot(0, 0, type = "n", xlim = xlim, ylim = ylim,
                  axes = FALSE, ...)
                graphics::axis(side = 1)
                graphics::axis(side = 2)
            }
            if (wpoints) {
              graphics::points(x$ashape.obj$x, col = col[3], ...)
            }
            if (number) {
                xoff <- 0.02 * diff(range(x$ashape.obj$x[, 1]))
                yoff <- 0.02 * diff(range(x$ashape.obj$x[, 2]))
                graphics::text(x$ashape.obj$x[, 1] + xoff, x$ashape.obj$x[,
                  2] + yoff, 1:(dim(x$ashape.obj$x)[1]), col = col[6],
                  ...)
            }
        }
        else {
            plot.delvor(x$ashape.obj$delvor.obj, add = add, wlines = wlines,
                wpoints = wpoints, number = number, col = col[3:6],
                lwd = lwd[1], xlim = xlim, ylim = ylim, ...)
        }
    }
    arcs <- which(x$arcs[, 3] > 0)
    if (length(arcs) > 0) {
        for (i in arcs) {
            arc(x$arcs[i, 1:2], x$arcs[i, 3], x$arcs[i, 4:5],
                x$arcs[i, 6], col = col[1], lwd = lwd[3])
        }
    }
    points <- which(x$arcs[, 3] == 0)
    if (length(points) > 0) {
        for (i in points) {
          graphics::points(x$arcs[i, 1], x$arcs[i, 2], col = col[1],
                pch = 19)
        }
    }
}
