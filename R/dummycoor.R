dummycoor <-
function (tri.obj, l1, l2, m, away) 
{
    if (is.list(tri.obj) && !is.null(tri.obj$x) && !is.null(tri.obj$y)) {
        center <- c(mean(tri.obj$x), mean(tri.obj$y))
    }
    else {
        xy <- xy.coords(tri.obj)
        center <- c(mean(xy$x), mean(xy$y))
    }
    v <- l2 - l1
    v <- c(v[2], -v[1])
    norm <- sqrt(sum(v^2))
    if (norm == 0) stop("Coincident Delaunay edge endpoints")
    v <- v/norm
    mp <- (l1 + l2)/2
    if (sum(v * (mp - center)) < 0) v <- -v

    # The marker must extend beyond the midpoint's projection whenever
    # that projection lies on the outward semi-infinite Voronoi edge.
    projection <- sum((mp - m) * v)
    dum <- m + (away + 2 * abs(projection)) * v
    return(dum)
}
