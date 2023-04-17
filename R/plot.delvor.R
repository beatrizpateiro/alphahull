#' Plot the Voronoi diagram and Delaunay traingulation
#'
#' This function returns a plot of the Voronoi diagram and Delaunay
#' traingulation.
#'
#'
#' @param x An object of class \code{"delvor"} as constructed by the function
#' delvor.
#' @param add Logical, if TRUE add to a current plot.
#' @param wlines "Which lines?". I.e. should the Delaunay triangulation be
#' plotted (wlines='del'), should the Voronoi diagram be plotted
#' (wlines='vor'), or should both be plotted (wlines='both', the default)?
#' @param wpoints Logical, indicates if sample points should be plotted.
#' @param number Logical, defaulting to FALSE; if TRUE then the points plotted
#' will be labelled with their index numbers.
#' @param col The colour numbers for plotting the data points, Delaunay
#' triangulation, Voronoi diagram, and the point numbers, in that order;
#' defaults to c(1,1,1,1). If fewer than four numbers are given, they are
#' recycled. (If more than four numbers are given, the redundant ones are
#' ignored.)
#' @param xlim The limits on the x-axis.
#' @param ylim The limits on the y-axis.
#' @param \dots Arguments to be passed to methods, such as graphical parameters
#' (see \code{\link{par}}).
#' @seealso \code{\link{delvor}}.
#' @keywords nonparametric
#' @examples
#'
#' \dontrun{
#' # Random sample in the unit square
#' x <- matrix(runif(100), nc = 2)
#' # Delaunay triangulation and Voronoi diagram
#' delvor.obj <- delvor(x)
#' # Plot Voronoi diagram and Delaunay triangulation
#' plot(delvor.obj)
#' }
#' @export
plot.delvor <-
function (x, add = FALSE, wlines = c("both", "del", "vor"), wpoints = TRUE,
    number = FALSE, col = NULL, xlim = NULL, ylim = NULL, ...)
{
    wlines <- match.arg(wlines)
    if (is.null(class(x)) || !inherits(x, "delvor")) {
        cat("Argument is not of class delvor.\n")
        return(invisible())
    }
    if (is.null(col)) {
        col <- c(1, 1, 1, 1)
    }
    else {
        col <- rep(col, length.out = 4)
    }
    plot.del <- switch(wlines, both = TRUE, del = TRUE, vor = FALSE)
    plot.vor <- switch(wlines, both = TRUE, del = FALSE, vor = TRUE)
    if (!add) {
        if (is.null(xlim))
            xlim <- range(x$x[, 1])
        if (is.null(ylim))
            ylim <- range(x$x[, 2])
        plot(0, 0, type = "n", xlim = xlim, ylim = ylim, axes = FALSE,
            ...)
        graphics::axis(side = 1)
        graphics::axis(side = 2)
    }
    if (wpoints) {
      graphics::points(x$x, col = col[1], ...)
    }
    if (number) {
        xoff <- 0.02 * diff(range(x$x[, 1]))
        yoff <- 0.02 * diff(range(x$x[, 2]))
        graphics::text(x$x[, 1] + xoff, x$x[, 2] + yoff, 1:(dim(x$x)[1]),
            col = col[4], ...)
    }
    if (plot.vor) {
        n.edge <- dim(x$mesh)[1]
        for (i in 1:n.edge) {
            lty.vor = ifelse(x$mesh[i, "bp1"] == 1 | x$mesh[i,
                "bp2"] == 1, 2, 1)
            graphics::segments(x$mesh[i, "mx1"], x$mesh[i, "my1"], x$mesh[i,
                "mx2"], x$mesh[i, "my2"], lty = lty.vor, col = col[3],
                ...)
        }
    }
    if (plot.del) {
        plot(x$tri.obj, add = TRUE, do.points = FALSE, col = col[2],
            ...)
    }
}
