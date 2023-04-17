#' Ahull according to the IUCN Red List Categories and Criteria
#'
#' @param x The x and y arguments provide the x and y coordinates of a set of points. Alternatively, a single argument x can be provided.
#' @param y The x and y arguments provide the x and y coordinates of a set of points. Alternatively, a single argument x can be provided.
#' @param alpha A number.
#' @return An ahull.IUCN object.
#' @examples
#' n <- 500
#' theta<-runif(n,0,2*pi)
#' r<-sqrt(runif(n,0.25^2,0.5^2))
#' x<-cbind(0.5+r*cos(theta),0.5+r*sin(theta))
#' alpha<-2.5
#' ah.IUCN.obj<-ahull.IUCN(x,alpha=alpha)
#' @references IUCN. 2022. The IUCN Red List of Threatened Species. Version 2022-2. https://www.iucnredlist.org
#' @references Burgman, M., & Fox, J. (2003). Bias in species range estimates from minimum convex polygons: Implications for conservation and options for improved planning. Animal Conservation Forum, 6(1), 19-28. doi:10.1017/S1367943003003044
#' @export
#' @importFrom grDevices xy.coords
#' @importFrom interp arcs area triangles tri.mesh

ahull.IUCN <- function(x, y = NULL, alpha){
  X <- grDevices::xy.coords(x, y)
  x <- cbind(X$x, X$y)
  n<-nrow(x)
  if (dim(x)[1] <= 2) {
    stop("At least three non-collinear points are required")
  }

  # Algorithm Burgman and Fox (2003)
  # ------------------------------------
  ## 1) Delauney triangulation
  tri.obj <- interp::tri.mesh(X)
  tri  <- interp::triangles(tri.obj)
  arcs <- interp::arcs(tri.obj)

  ## 2) Length of lines
  lengtharcs<-sqrt(rowSums((x[arcs[,1],]-x[arcs[,2],])^2))
  avl<-mean(lengtharcs)

  ## 3) Delete all lines that are longer than a multiple (alpha) of the average line length.
  karcs<-(lengtharcs<=alpha*avl)       # arcs to keep (T/F)
  wkarcs<-which(lengtharcs<=alpha*avl) # which arcs to keep
  ktri<-apply(matrix(tri[,c("arc1","arc2","arc3")]%in%wkarcs,ncol=3),1,all) # triangles to keep (those with all 3 arcs in wkarcs)
  wktri<-which(ktri)                   # which triangles to keep

  aux<-tri[,c("tr1","tr2","tr3")]%in%wktri
  triIUCN<-tri[,c("tr1","tr2","tr3")]
  triIUCN[!aux]<-0

  aux<-tri[,c("arc1","arc2","arc3")]%in%wkarcs
  arcIUCN<-tri[,c("arc1","arc2","arc3")]
  arcIUCN[!aux]<-0

  # For each triangle in Delaunay triangulation (for each row in tri)
  # triIUCN[,1]   ---- TRUE/FALSE indicating whether the complete triangle i belongs to the ahull
  # triIUCN[,2:4] ---- indices of neighbour triangles belonging to the ahul  (0 if they do no belong to the ahull)
  # triIUCN[,5:7] ---- indices of arcs of the triangle belonging to the ahul (0 if they do no belong to the ahull)

  triIUCN<-cbind("tri.in.ah"=ktri,triIUCN,arcIUCN)

  ##(4) Calculate the area of habitat by summing the areas of all remaining triangles.
  areaIUCN<-sum(interp::area(tri.obj)[ktri]) # Area of the ahull (sum of the areas of the triangles belonging to the ahull)

  aux<-as.numeric(triIUCN[triIUCN[,1]==0,c("arc1","arc2","arc3")])
  aux<-aux[aux!=0]
  barcsIUCN<-unique(aux)  # Arcs in the boundary of the ahull from triangles that don't belong to the ahull
  barcsIUCN<-sort(union(barcsIUCN,tri[ktri,7:9][tri[ktri,4:6]==0])) # Arcs in the boundary of the ahull
  isolp<-sort(setdiff(1:n,arcs[karcs])) # Isolated points belonging to the ahull (those that do not belong to an arc in the ahull)

  tri.ah.IUCN<-matrix(tri[ktri,1:3],ncol=3)
  colnames(tri.ah.IUCN)<-c("node1","node2","node3")

  bd.ah.IUCN=matrix(arcs[barcsIUCN,],ncol=2)
  colnames(bd.ah.IUCN)<-c("from","to")

  ahullIUCN.obj <- list(tri.ah.IUCN=tri.ah.IUCN, bd.ah.IUCN=bd.ah.IUCN, ip.ah.IUCN=isolp, area=areaIUCN, tri = tri.obj, alpha = alpha, x=x)
  class(ahullIUCN.obj) <- "ahull.IUCN"
  invisible(ahullIUCN.obj)
}
