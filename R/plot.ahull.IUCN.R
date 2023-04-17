#' Very simple plot function
#' @param x The ahull.IUCN object.
#' @param wplot "Which in the plot?". I.e. should all points segments and connected components be plat (wlines='all', the default), should only connected components be plotted (wlines='comp')
#' @param ... Additional plot parameters.
#' @return None.
#' @export
#' @importFrom graphics points polygon segments
#' @examples
#' n <- 500
#' theta<-runif(n,0,2*pi)
#' r<-sqrt(runif(n,0.25^2,0.5^2))
#' x<-cbind(0.5+r*cos(theta),0.5+r*sin(theta))
#' alpha<-2.5
#' ah<-ahull.IUCN(x,alpha=alpha)
#' plot(ah)


plot.ahull.IUCN<-function(x, wplot = c("all", "comp"), ...){
  ah.IUCN.obj<-x
  x <- ah.IUCN.obj$x
  tri.ah<-ah.IUCN.obj$tri.ah.IUCN
  edges.ah<-ah.IUCN.obj$bd.ah.IUCN
  ip.ah<-ah.IUCN.obj$ip.ah.IUCN

  plot(x,main=paste("IUCN ahull for alpha =",ah.IUCN.obj$alpha,"\n Area =",round(ah.IUCN.obj$area,digits=4)),xlab="",ylab="",type="n",...)
  if(nrow(tri.ah)>0){
    for(i in 1:nrow(tri.ah)){
      graphics::polygon(x[tri.ah[i,],1],x[tri.ah[i,],2],col=3,lty=2)       # Triangles in the ahull
    }
  }
  wplot <- match.arg(wplot)
  plot.all <- switch(wplot, all = TRUE, comp=FALSE)
  if(plot.all){
    graphics::segments(x[edges.ah[,1],1],x[edges.ah[,1],2],x[edges.ah[,2],1],x[edges.ah[,2],2],col=4,lwd=3) # Boundary edges in the ahull in blue
    graphics::points(x[ip.ah,],col=2,pch=19)  # Isolated points in the ahull in red (if any)
  }
}
