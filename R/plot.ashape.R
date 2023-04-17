#' Plot the alpha-shape
#'
#' This function returns a plot of the \eqn{\alpha}-shape. If desired, it also
#' adds the Delaunay triangulation and Voronoi diagram of the sample.
#'
#'
#' @param x Object of class \code{"ashape"}.
#' @param add Logical, if TRUE add to a current plot.
#' @param wlines "Which lines?". I.e. should the Delaunay triangulation be
#' plotted (wlines='del'), should the Voronoi diagram be plotted
#' (wlines='vor'), should both be plotted (wlines='both'), or none
#' (wlines='none', the default)?
#' @param wpoints Logical, indicates if sample points should be plotted.
#' @param number Logical, defaulting to FALSE; if TRUE then the points plotted
#' will be labelled with their index numbers.
#' @param col The colour numbers for plotting the \eqn{\alpha}-shape, data
#' points, Delaunay triangulation, Voronoi diagram, and the point numbers, in
#' that order; defaults to c(1,1,1,1,1). If fewer than five numbers are given,
#' they are recycled. (If more than five numbers are given, the redundant ones
#' are ignored.)
#' @param xlim The limits on the x-axis.
#' @param ylim The limits on the y-axis.
#' @param lwd The line widths for plotting the tesselations and the
#' \eqn{\alpha}-shape; defaults to c(1,2).
#' @param \dots Arguments to be passed to methods, such as graphical parameters
#' (see \code{\link{par}}).
#' @seealso \code{\link{ashape}}.
#' @keywords nonparametric
#' @examples
#'
#' \dontrun{
#' # Uniform sample of size n=300 in the annulus B(c, 0.5)\B(c, 0.25)
#' # with c=(0.5, 0.5).
#' n <- 300
#' theta<-runif(n,0,2*pi)
#' r<-sqrt(runif(n,0.25^2,0.5^2))
#' x<-cbind(0.5+r*cos(theta),0.5+r*sin(theta))
#' # Value of alpha
#' alpha <- 0.1
#' # alpha-shape
#' ashape.obj <- ashape(x, alpha = alpha)
#' # Plot alpha-shape in blue, sample points in black,
#' # and Delaunay triangulation in red
#' plot(ashape.obj, wlines= "del", col = c(4, 1, 2))
#'
#' # Random sample  from a uniform distribution on a Koch snowflake
#' # with initial side length 1 and 3 iterations
#' x <- rkoch(2000, side = 1, niter = 3)
#' # Value of alpha
#' alpha <- 0.05
#' # alpha-shape
#' ashape.obj <- ashape(x, alpha = alpha)
#' # Plot alpha-shape in blue
#' plot(ashape.obj, col = c(4, 1))
#' }
#' @export
plot.ashape <-
function (x, add = FALSE, wlines = c("none", "both", "del", "vor"),
    wpoints = TRUE, number = FALSE, col = NULL, xlim = NULL,
    ylim = NULL, lwd = NULL, ...)
{
    wlines <- match.arg(wlines)
    if (is.null(class(x)) || !inherits(x, "ashape")) {
        cat("Argument is not of class ashape.\n")
        return(invisible())
    }
    if (is.null(col)) {
        col <- c(1, 1, 1, 1, 1)
    }
    else {
        col <- rep(col, length.out = 5)
    }
    if (is.null(lwd)) {
        lwd <- 1:2
    }
    else {
        lwd <- rep(lwd, length.out = 2)
    }
    wlines <- match.arg(wlines)
    plot.dd <- switch(wlines, none = TRUE, both = FALSE, del = FALSE,
        vor = FALSE)
    if (plot.dd) {
        if (!add) {
            if (is.null(xlim))
                xlim <- range(x$x[, 1])
            if (is.null(ylim))
                ylim <- range(x$x[, 2])
            plot(0, 0, type = "n", xlim = xlim, ylim = ylim,
                axes = FALSE, ...)
            graphics::axis(side = 1)
            graphics::axis(side = 2)
        }
        if (wpoints) {
          graphics::points(x$x, col = col[2], ...)
        }
        if (number) {
            xoff <- 0.02 * diff(range(x$x[, 1]))
            yoff <- 0.02 * diff(range(x$x[, 2]))
            graphics::text(x$x[, 1] + xoff, x$x[, 2] + yoff, 1:(dim(x$x)[1]),
                col = col[5], ...)
        }
    }
    else {
        plot.delvor(x$delvor.obj, add = add, wlines = wlines,
            wpoints = wpoints, number = number, col = col[2:5],
            lwd = lwd[1], xlim = xlim, ylim = ylim, ...)
    }
    ashape <- x$edges
    n2 <- dim(ashape)[1]
    if (n2 >= 1) {
        for (i in 1:n2) {
          graphics::segments(ashape[i, "x1"], ashape[i, "y1"], ashape[i,
                "x2"], ashape[i, "y2"], col = col[1], lwd = lwd[2])
        }
    }
}
