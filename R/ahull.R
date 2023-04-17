#' alpha-convex hull calculation
#' 
#' This function calculates the \eqn{\alpha}-convex hull of a given sample of
#' points in the plane for \eqn{\alpha>0}.
#' 
#' 
#' An attempt is made to interpret the arguments x and y in a way suitable for
#' computing the \eqn{\alpha}-convex hull. Any reasonable way of defining the
#' coordinates is acceptable, see \code{\link{xy.coords}}.
#' 
#' The \eqn{\alpha}-convex hull is defined for any finite number of points.
#' However, since the algorithm is based on the Delaunay triangulation, at
#' least three non-collinear points are required.
#' 
#' If \code{y} is NULL and \code{x} is an object of class \code{"delvor"}, then
#' the \eqn{\alpha}-convex hull is computed with no need to invoke again the
#' function \code{\link{delvor}} (it reduces the computational cost).
#' 
#' The complement of the \eqn{\alpha}-convex hull can be written as the union
#' of \eqn{O(n)} open balls and halfplanes, see \code{\link{complement}}.  The
#' boundary of the \eqn{\alpha}-convex hull is formed by arcs of open balls of
#' radius \eqn{\alpha} (besides possible isolated sample points). The arcs are
#' determined by the intersections of some of the balls that define the
#' complement of the \eqn{\alpha}-convex hull. The extremes of an arc are given
#' by \eqn{c+rA_\theta v} and \eqn{c+rA_{-\theta}v} where \eqn{c} and \eqn{r}
#' represent the center and radius of the arc, repectively, and \eqn{A_\theta
#' v} represents the clockwise rotation of angle \eqn{\theta} of the unitary
#' vector \eqn{v}. Joining the end points of adjacent arcs we can define
#' polygons that help us to determine the area of the estimator , see
#' \code{\link{areaahull}}.
#' 
#' @param x,y The \code{x} and \code{y} arguments provide the \code{x} and
#' \code{y} coordinates of a set of points. Alternatively, a single argument
#' \code{x} can be provided, see Details.
#' @param alpha Value of \eqn{\alpha}.
#' @return A list with the following components: \item{arcs}{For each arc in
#' the boundary of the \eqn{\alpha}-convex hull, the columns of the matrix
#' \code{arcs} store the center \eqn{c} and radius \eqn{r} of the arc, the
#' unitary vector \eqn{v}, the angle \eqn{\theta} that define the arc and the
#' indices of the end points, see Details. The coordinates of the end points of
#' the arcs are stored in \code{xahull}. For isolated points in the boundary of
#' the \eqn{\alpha}-convex hull, columns 3 to 6 of the matrix \code{arcs} are
#' equal to zero. } \item{xahull}{A 2-column matrix with the coordinates of the
#' original set of points besides possible new end points of the arcs in the
#' boundary of the \eqn{\alpha}-convex hull.} \item{length}{Length of the
#' boundary of the \eqn{\alpha}-convex hull, see \code{\link{lengthahull}}.}
#' \item{complement}{Output matrix from \code{\link{complement}}.}
#' \item{alpha}{Value of \eqn{\alpha}.} \item{ashape.obj}{Object of class
#' \code{"ashape"} returned by the function \code{\link{ashape}}.}
#' @seealso \code{\link{plot.ahull}}.
#' @references Edelsbrunner, H., Kirkpatrick, D.G. and Seidel, R. (1983). On
#' the shape of a set of points in the plane. \emph{IEEE Transactions on
#' Information Theory}, 29(4), pp.551-559.
#' 
#' Rodriguez-Casal, R. (2007). Set estimation under convexity type assumptions.
#' \emph{Annales de l'I.H.P.- Probabilites & Statistiques}, 43, pp.763-774.
#' 
#' Pateiro-Lopez, B. (2008). \emph{Set estimation under convexity type
#' restrictions}. Phd. Thesis. Universidad de Santiago de Compostela. ISBN
#' 978-84-9887-084-8.
#' @keywords nonparametric
#' @examples
#' 
#' \dontrun{
#' # Random sample in the unit square
#' x <- matrix(runif(100), nc = 2)
#' # Value of alpha
#' alpha <- 0.2
#' # Alpha-convex hull
#' ahull.obj <- ahull(x, alpha = alpha)
#' plot(ahull.obj)
#' 
#' # Uniform sample of size n=300 in the annulus B(c,0.5)\B(c,0.25), 
#' # with c=(0.5,0.5). 
#' n <- 300
#' theta<-runif(n,0,2*pi)
#' r<-sqrt(runif(n,0.25^2,0.5^2))
#' x<-cbind(0.5+r*cos(theta),0.5+r*sin(theta))
#' # Value of alpha
#' alpha <- 0.1
#' # Alpha-convex hull
#' ahull.obj <- ahull(x, alpha = alpha)
#' # The arcs defining the boundary of the alpha-convex hull are ordered
#' plot(x)
#' for (i in 1:dim(ahull.obj$arcs)[1]){
#' arc(ahull.obj$arcs[i,1:2],ahull.obj$arcs[i,3],ahull.obj$arcs[i,4:5],
#' ahull.obj$arcs[i,6],col=2)
#' Sys.sleep(0.5)
#' }
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
#' 
#' @export ahull
ahull <-
function (x, y = NULL, alpha) 
{
    ashape.obj <- ashape(x, y, alpha)
    compl <- complement(ashape.obj$delvor.obj, alpha = alpha)
    pm.x <- (compl[, "x1"] + compl[, "x2"]) * 0.5
    pm.y <- (compl[, "y1"] + compl[, "y2"]) * 0.5
    dm <- sqrt((compl[, "x1"] - compl[, "x2"])^2 + (compl[, "y1"] - 
        compl[, "y2"])^2) * 0.5
    ashape.edges <- matrix(ashape.obj$edges[, c("ind1", "ind2")], 
        ncol = 2, byrow = FALSE)
    noforget <- ashape.obj$alpha.extremes
    ind2 <- integer()
    j <- 0
    nshape <- length(ashape.edges) * 0.5
    arcs <- matrix(0, nrow = nshape, ncol = 6)
    indp <- matrix(0, nrow = nshape, ncol = 2)
    cutp <- ashape.obj$x
    if (nshape > 0) {
        for (i in 1:nshape) {
            ind <- which(ashape.edges[i, 1] == compl[, "ind1"] & 
                ashape.edges[i, 2] == compl[, "ind2"])
            if (length(ind) > 0) {
                if (!((1 <= sum((compl[ind, "ind"] == 1))) & 
                  (sum((compl[ind, "ind"] == 1)) < length(ind)))) {
                  which <- which(compl[ind, "r"] == min(compl[ind[compl[ind, 
                    "r"] > 0], "r"]))
                  j <- j + 1
                  arcs[j, ] <- c(compl[ind[which], 1], compl[ind[which], 
                    2], compl[ind[which], 3], compl[ind[which], 
                    "v.x"], compl[ind[which], "v.y"], compl[ind[which], 
                    "theta"])
                  vaux <- compl[ind[which], c("x1", "y1")] - 
                    cbind(pm.x, pm.y)[ind[which], ]
                  theta.aux <- rotation(compl[ind[which], c("v.x", 
                    "v.y")], compl[ind[which], "theta"])
                  a2 <- sum(vaux * theta.aux)
                  if (a2 > 0) {
                    indp[j, ] <- compl[ind[which], c("ind1", 
                      "ind2")]
                  }
                  else {
                    indp[j, ] <- compl[ind[which], c("ind2", 
                      "ind1")]
                  }
                }
                ind2 <- c(ind2, ind)
            }
        }
    }
    arcs.old <- arcs[arcs[, 3] > 0, ]
    colnames(arcs.old) <- c("c1", "c2", "r", "v.x", "v.y", "theta")
    arcs <- arcs.old
    indp <- indp[indp[, 1] != 0 & indp[, 2] != 0, ]
    n.arc <- dim(arcs)[1]
    watch <- 1
    j <- 1
    if (n.arc > 0) {
        while (watch <= n.arc) {
            ind.arc <- 1:n.arc
            while (j <= n.arc) {
                if (j != watch) {
                  intersection <- inter(arcs[watch, 1], arcs[watch, 
                    2], arcs[watch, 3], arcs[j, 1], arcs[j, 2], 
                    arcs[j, 3])
                  if (intersection$n.cut == 2) {
                    v.arc <- c(arcs[watch, "v.x"], arcs[watch, 
                      "v.y"])
                    if (v.arc[2] >= 0) {
                      ang.OX <- acos(v.arc[1])
                    }
                    else {
                      ang.OX <- 2 * pi - acos(v.arc[1])
                    }
                    v.arc.rot <- rotation(v.arc, ang.OX)
                    v.int <- intersection$v1
                    v.int.rot <- rotation(v.int, ang.OX)
                    if (v.int.rot[2] >= 0) {
                      ang.v.int.rot.OX <- acos(v.int.rot[1])
                      angles <- c(-arcs[watch, "theta"], arcs[watch, 
                        "theta"], ang.v.int.rot.OX - intersection$theta1, 
                        ang.v.int.rot.OX + intersection$theta1)
                      names(angles) <- c("theta1", "theta2", 
                        "beta1", "beta2")
                      order <- names(sort(angles))
                      theta1 <- ang.OX - arcs[watch, "theta"]
                      theta2 <- ang.OX + arcs[watch, "theta"]
                      beta1 <- ang.v.int.rot.OX + ang.OX - intersection$theta1
                      beta2 <- ang.v.int.rot.OX + ang.OX + intersection$theta1
                    }
                    else {
                      ang.v.int.rot.OX <- acos(v.int.rot[1])
                      angles <- c(-arcs[watch, "theta"], arcs[watch, 
                        "theta"], -ang.v.int.rot.OX - intersection$theta1, 
                        -ang.v.int.rot.OX + intersection$theta1)
                      names(angles) <- c("theta1", "theta2", 
                        "beta1", "beta2")
                      order <- names(sort(angles))
                      theta1 <- ang.OX - arcs[watch, "theta"]
                      theta2 <- ang.OX + arcs[watch, "theta"]
                      beta1 <- -ang.v.int.rot.OX + ang.OX - intersection$theta1
                      beta2 <- -ang.v.int.rot.OX + ang.OX + intersection$theta1
                    }
                    if (sum(match(indp[watch, ], indp[j, ], nomatch = 0)) > 
                      0) {
                      coinc <- indp[j, sum(match(indp[watch, 
                        ], indp[j, ], nomatch = 0))]
                      if (indp[watch, 1] == indp[j, 2]) {
                        if (all(order == c("beta1", "beta2", 
                          "theta1", "theta2"))) {
                          case <- 1
                        }
                        else if (all(order == c("theta1", "beta1", 
                          "beta2", "theta2"))) {
                          case <- 2
                        }
                        else if (all(order == c("beta1", "theta1", 
                          "beta2", "theta2"))) {
                          ang.control <- (angles["theta1"] - 
                            angles["beta1"])/2
                          if (abs(ang.control) < 1e-05) {
                            case <- 2
                          }
                          else {
                            case <- 1
                          }
                        }
                        if (case == 2) {
                          ang.middle2 <- (angles["theta2"] - 
                            angles["beta2"])/2
                          v.new2 <- rotation(c(1, 0), -arcs[watch, 
                            "theta"] + ang.middle2 - ang.OX)
                          cutp <- rbind(cutp, arcs[watch, 1:2] + 
                            arcs[watch, 3] * rotation(v.new2, 
                              ang.middle2))
                          inn <- dim(cutp)[1]
                          arcs[watch, 4:6] <- c(v.new2, ang.middle2)
                          indp[watch, 1] <- inn
                          pmaux <- (cutp[inn, ] + cutp[indp[j, 
                            1], ]) * 0.5
                          dmaux <- pmaux - cutp[indp[j, 1], ]
                          ndmaux <- sqrt(sum(dmaux^2))
                          vaux <- pmaux - arcs[j, 1:2]
                          nvaux <- sqrt(sum(vaux^2))
                          th <- atan(ndmaux/nvaux)
                          arcs[j, 4:6] <- c(vaux/nvaux, th)
                          indp[j, 2] <- inn
                        }
                      }
                      else if (indp[watch, 2] == indp[j, 1]) {
                        if (all(order == c("theta1", "theta2", 
                          "beta1", "beta2"))) {
                          case <- 1
                        }
                        else if (all(order == c("theta1", "beta1", 
                          "beta2", "theta2"))) {
                          case <- 2
                        }
                        else if (all(order == c("theta1", "beta1", 
                          "theta2", "beta2"))) {
                          ang.control <- (angles["theta2"] - 
                            angles["beta2"])/2
                          if (abs(ang.control) < 1e-05) {
                            case <- 2
                          }
                          else {
                            case <- 1
                          }
                        }
                        if (case == 2) {
                          ang.middle <- (angles["beta1"] - angles["theta1"])/2
                          v.new <- rotation(c(1, 0), arcs[watch, 
                            "theta"] - ang.middle - ang.OX)
                          cutp <- rbind(cutp, arcs[watch, 1:2] + 
                            arcs[watch, 3] * rotation(v.new, 
                              -ang.middle))
                          inn <- dim(cutp)[1]
                          arcs[watch, 4:6] <- c(v.new, ang.middle)
                          indp[watch, 2] <- inn
                          pmaux <- (cutp[inn, ] + cutp[indp[j, 
                            2], ]) * 0.5
                          dmaux <- pmaux - cutp[indp[j, 2], ]
                          ndmaux <- sqrt(sum(dmaux^2))
                          vaux <- pmaux - arcs[j, 1:2]
                          nvaux <- sqrt(sum(vaux^2))
                          th <- atan(ndmaux/nvaux)
                          arcs[j, 4:6] <- c(vaux/nvaux, th)
                          indp[j, 1] <- inn
                        }
                      }
                    }
                    else if (all(order == c("theta1", "beta1", 
                      "beta2", "theta2"))) {
                      v.arcj <- c(arcs[j, "v.x"], arcs[j, "v.y"])
                      if (v.arcj[2] >= 0) {
                        ang.OXj <- acos(v.arcj[1])
                      }
                      else {
                        ang.OXj <- 2 * pi - acos(v.arcj[1])
                      }
                      v.arc.rotj <- rotation(v.arcj, ang.OXj)
                      v.intj <- intersection$v2
                      v.int.rotj <- rotation(v.intj, ang.OXj)
                      if (v.int.rotj[2] >= 0) {
                        ang.v.int.rot.OXj <- acos(v.int.rotj[1])
                        anglesj <- c(-arcs[j, "theta"], arcs[j, 
                          "theta"], ang.v.int.rot.OXj - intersection$theta2, 
                          ang.v.int.rot.OXj + intersection$theta2)
                        names(anglesj) <- c("theta1", "theta2", 
                          "beta1", "beta2")
                        orderj <- names(sort(anglesj))
                        theta1j <- ang.OXj - arcs[j, "theta"]
                        theta2j <- ang.OXj + arcs[j, "theta"]
                        beta1j <- ang.v.int.rot.OXj + ang.OXj - 
                          intersection$theta2
                        beta2j <- ang.v.int.rot.OXj + ang.OXj + 
                          intersection$theta2
                      }
                      else {
                        ang.v.int.rot.OXj <- acos(v.int.rotj[1])
                        anglesj <- c(-arcs[j, "theta"], arcs[j, 
                          "theta"], -ang.v.int.rot.OXj - intersection$theta2, 
                          -ang.v.int.rot.OXj + intersection$theta2)
                        names(anglesj) <- c("theta1", "theta2", 
                          "beta1", "beta2")
                        orderj <- names(sort(anglesj))
                        theta1j <- ang.OXj - arcs[j, "theta"]
                        theta2j <- ang.OXj + arcs[j, "theta"]
                        beta1j <- -ang.v.int.rot.OXj + ang.OXj - 
                          intersection$theta2
                        beta2j <- -ang.v.int.rot.OXj + ang.OXj + 
                          intersection$theta2
                      }
                      if (all(orderj == c("theta1", "beta1", 
                        "beta2", "theta2"))) {
                        ang.middle <- (angles["beta1"] - angles["theta1"])/2
                        v.new <- rotation(c(1, 0), arcs[watch, 
                          "theta"] - ang.middle - ang.OX)
                        ang.middle2 <- (angles["theta2"] - angles["beta2"])/2
                        v.new2 <- rotation(c(1, 0), -arcs[watch, 
                          "theta"] + ang.middle2 - ang.OX)
                        arcs <- rbind(arcs, c(arcs[watch, 1], 
                          arcs[watch, 2], arcs[watch, 3], v.new2[1], 
                          v.new2[2], ang.middle2))
                        arcs[watch, ] <- c(arcs[watch, 1], arcs[watch, 
                          2], arcs[watch, 3], v.new[1], v.new[2], 
                          ang.middle)
                        n.arc <- n.arc + 1
                        np1 <- arcs[watch, 1:2] + arcs[watch, 
                          3] * rotation(v.new, -ang.middle)
                        np2 <- arcs[watch, 1:2] + arcs[watch, 
                          3] * rotation(v.new2, ang.middle2)
                        indold <- indp[watch, 2]
                        inn1 <- dim(cutp)[1] + 1
                        inn2 <- dim(cutp)[1] + 2
                        indp[watch, 2] <- inn1
                        indp <- rbind(indp, c(inn2, indold))
                        cutp <- rbind(cutp, np1)
                        cutp <- rbind(cutp, np2)
                        indold <- indp[j, 1]
                        indp[j, 1] <- inn1
                        indp <- rbind(indp, c(indold, inn2))
                        pmaux <- (cutp[inn1, ] + cutp[indp[j, 
                          2], ]) * 0.5
                        dmaux <- pmaux - cutp[indp[j, 2], ]
                        ndmaux <- sqrt(sum(dmaux^2))
                        vaux <- pmaux - arcs[j, 1:2]
                        nvaux <- sqrt(sum(vaux^2))
                        th <- atan(ndmaux/nvaux)
                        arcs[j, 4:6] <- c(vaux/nvaux, th)
                        pmaux <- (cutp[inn2, ] + cutp[indold, 
                          ]) * 0.5
                        dmaux <- pmaux - cutp[indold, ]
                        ndmaux <- sqrt(sum(dmaux^2))
                        vaux <- pmaux - arcs[j, 1:2]
                        nvaux <- sqrt(sum(vaux^2))
                        th <- atan(ndmaux/nvaux)
                        arcs <- rbind(arcs, c(arcs[j, 1:3], vaux/nvaux, 
                          th))
                        n.arc <- n.arc + 1
                      }
                    }
                  }
                  case <- 0
                }
                j <- j + 1
            }
            watch <- watch + 1
            j <- 1
        }
        ord.old <- 1:dim(indp)[1]
        ord.new <- numeric()
        while (length(ord.new) < length(ord.old)) {
            if (length(ord.new) == 0) {
                ord.new <- 1
            }
            else {
                ord.new <- c(ord.new, ord.old[-ord.new][1])
            }
            coinc <- match(indp[ord.new[length(ord.new)], 2], 
                indp[-ord.new, 1])
            while (!is.na(coinc)) {
                ord.new <- c(ord.new, ord.old[-ord.new][coinc])
                coinc <- match(indp[ord.new[length(ord.new)], 
                  2], indp[-ord.new, 1])
            }
        }
        indp <- indp[ord.new, ]
        ahull.arcs <- cbind(arcs[ord.new, ], indp)
        colnames(ahull.arcs) <- c("c1", "c2", "r", "v.x", "v.y", 
            "theta", "end1", "end2")
        lengthah <- lengthahull(arcs)
        addp <- noforget[is.na(match(noforget, indp))]
        num <- length(addp)
        if (num > 0) {
            mat.noforget <- cbind(matrix(ashape.obj$x[addp, 1:2], 
                ncol = 2, byrow = FALSE), rep(0, num), rep(0, 
                num), rep(0, num), rep(0, num), addp, addp)
            ahull.arcs <- rbind(ahull.arcs, mat.noforget)
        }
    }
    else {
        num <- length(noforget)
        ahull.arcs <- cbind(ashape.obj$x[noforget, ], rep(0, 
            num), rep(0, num), rep(0, num), rep(0, num), noforget, 
            noforget)
        colnames(ahull.arcs) <- c("c1", "c2", "r", "v.x", "v.y", 
            "theta", "end1", "end2")
        lengthah <- 0
    }
    ahull.obj <- list(arcs = ahull.arcs, xahull = cutp, length = lengthah, 
        complement = compl, alpha = alpha, ashape.obj = ashape.obj)
    class(ahull.obj) <- "ahull"
    invisible(ahull.obj)
}
