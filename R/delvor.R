delvor <- function(x, y = NULL)
{
    X <- xy.coords(x, y)
    x <- cbind(X$x, X$y)
    if (nrow(x) < 3L || any(!is.finite(x))) {
        stop("At least three finite, non-collinear points are required")
    }
    if (anyDuplicated(as.data.frame(x))) {
        stop("Duplicate points are not supported")
    }

    # Each row contains the original point indices of a Delaunay triangle.
    # Reconstruct adjacency instead of using interp's neighbour flags.
    if (nrow(x) == 3L) {
        tri <- matrix(1:3, nrow = 1L)
    }
    else {
        tri <- geometry::delaunayn(x)
    }
    tri <- matrix(as.integer(tri), ncol = 3L)
    if (nrow(tri) == 0L) {
        stop("At least three non-collinear points are required")
    }
    if (any(tri < 1L | tri > nrow(x)) ||
        length(unique(as.vector(tri))) != nrow(x)) {
        stop("The triangulation omits or misindexes input points")
    }
    colnames(tri) <- c("node1", "node2", "node3")
    nt <- nrow(tri)

    # Circumcentres relative to the first vertex, avoiding cancellation
    # from differences of large absolute squared coordinates.
    a <- x[tri[, 1L], , drop = FALSE]
    u <- x[tri[, 2L], , drop = FALSE] - a
    v <- x[tri[, 3L], , drop = FALSE] - a
    det <- u[, 1L] * v[, 2L] - u[, 2L] * v[, 1L]
    if (any(det == 0) || any(!is.finite(det))) {
        stop("Degenerate Delaunay triangle")
    }
    ru <- rowSums(u * u) / 2
    rv <- rowSums(v * v) / 2
    circumcentres <- cbind(
        a[, 1L] + (ru * v[, 2L] - u[, 2L] * rv) / det,
        a[, 2L] + (u[, 1L] * rv - ru * v[, 1L]) / det
    )
    if (any(!is.finite(circumcentres))) {
        stop("Non-finite Delaunay circumcentre")
    }

    # Triangle incidences of each undirected Delaunay edge.
    sides <- rbind(tri[, c(1L, 2L), drop = FALSE],
                   tri[, c(2L, 3L), drop = FALSE],
                   tri[, c(3L, 1L), drop = FALSE])
    triangle.id <- rep(seq_len(nt), 3L)
    lo <- pmin(sides[, 1L], sides[, 2L])
    hi <- pmax(sides[, 1L], sides[, 2L])
    ord <- order(lo, hi)
    lo <- lo[ord]
    hi <- hi[ord]
    triangle.id <- triangle.id[ord]
    first <- c(TRUE, lo[-1L] != lo[-length(lo)] |
                         hi[-1L] != hi[-length(hi)])
    start <- which(first)
    count <- diff(c(start, length(lo) + 1L))
    if (any(count > 2L)) {
        stop("An edge has more than two incident triangles")
    }

    endpoints <- cbind(lo[start], hi[start])
    t1 <- triangle.id[start]
    t2 <- rep(NA_integer_, length(start))
    interior <- count == 2L
    t2[interior] <- triangle.id[start[interior] + 1L]
    boundary <- !interior

    # Each boundary edge must lie on a supporting line of the convex hull.
    # This also permits multiple exactly collinear sites along a hull side:
    # chull() reports only its extremes, but Qhull can subdivide that side.
    hull.points <- x[chull(x), , drop = FALSE]
    scale <- max(diff(range(x[, 1L])), diff(range(x[, 2L])))
    for (i in which(boundary)) {
        p <- x[endpoints[i, 1L], ]
        d <- x[endpoints[i, 2L], ] - p
        cross <- d[1L] * (hull.points[, 2L] - p[2L]) -
                 d[2L] * (hull.points[, 1L] - p[1L])
        tol <- 128 * .Machine$double.eps * scale * sqrt(sum(d * d))
        if (any(cross > tol) && any(cross < -tol)) {
            stop("Delaunay boundary contains a non-supporting edge")
        }
    }

    m1 <- circumcentres[t1, , drop = FALSE]
    m2 <- matrix(NA_real_, nrow(endpoints), 2L)
    m2[interior, ] <- circumcentres[t2[interior], , drop = FALSE]
    away <- max(diff(range(x[, 1L])), diff(range(x[, 2L])))
    points.obj <- list(x = X$x, y = X$y)
    for (i in which(boundary)) {
        m2[i, ] <- dummycoor(points.obj,
                            x[endpoints[i, 1L], ],
                            x[endpoints[i, 2L], ],
                            m1[i, ], away)
    }

    mesh <- cbind(endpoints,
                  x[endpoints[, 1L], , drop = FALSE],
                  x[endpoints[, 2L], , drop = FALSE],
                  m1, m2,
                  bp1 = rep(FALSE, nrow(endpoints)),
                  bp2 = boundary)
    colnames(mesh) <- c("ind1", "ind2", "x1", "y1", "x2", "y2",
                        "mx1", "my1", "mx2", "my2", "bp1", "bp2")

    # Documented change: tri.obj is now a triangle index matrix.
    delvor.obj <- list(mesh = mesh, x = x, tri.obj = tri)
    class(delvor.obj) <- "delvor"
    invisible(delvor.obj)
}
