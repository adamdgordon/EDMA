#' Plot Object Defined by Landmarks
#' 
#' Plot a set of landmarks.
#' @param x A matrix of landmarks where rows correspond to landmarks and columns correspond to X, Y, and Z.
#' @param pt.cols ...
#' @param pt.size ...
#' @param pt.alpha ...
#' @param label ...
#' @param label.cex ... 
#' @param segs ...
#' @param tris ...
#' @param ... Arguments to be passed to other functions.  Not currently used.
#' @examples
#' plot_lmks(guenons$rawcoords[,,guenons$species=="patas"][,,30])
#' plot_lmks(guenons$rawcoords[,,guenons$species=="patas"][,,30],
#'           segs=guenons$wireframe_segments,
#' 		     tris=guenons$triangles,
#' 		     label=TRUE)
#' @export
plot_lmks <- function(x, pt.cols="blue", pt.size=0.5, pt.alpha=0.5, label=FALSE, label.cex=0.7, segs=NULL, tris=NULL, ...) {
  rgl::plot3d(x,
              asp=1,
			  type="s",
              col=pt.cols,
			  size=pt.size,
			  alpha=pt.alpha)
  if (!is.null(segs)) rgl::segments3d(x=x[as.vector(t(segs[,1:2])),])
  if (!is.null(tris)) rgl::triangles3d(x=x[as.vector(t(tris[,1:3])),],
                                       col="black", alpha=0.2)
  if (label) rgl::text3d(x=x,
                         texts=rownames(x),
			             pos=1,
						 cex=label.cex)
}

#' Plot \code{FMmean} Object
#' 
#' Plot the \code{FMmean} object produced by the function \code{meanform}.
#' @param x A matrix of landmarks where rows correspond to landmarks and columns correspond to X, Y, and Z.
#' @param pt.cols ...
#' @param pt.size ...
#' @param pt.alpha ...
#' @param label ...
#' @param label.cex ... 
#' @param segs ...
#' @param tris ...
#' @param ... Arguments to be passed to other functions.  Not currently used.
#' @examples
#' patasmean<- meanform(guenons$rawcoords[,,guenons$species=="patas"])
#' plot(patasmean)
#' plot(patasmean,
#'      segs=guenons$wireframe_segments,
#' 		tris=guenons$triangles,
#' 		label=TRUE, pt.alpha=0.2, pt.size=0.3)
#' @export
plot.FMmean <- function(x, segs=NULL, tris=NULL,
                        pt.cols="blue",
						pt.col.scale=TRUE,
						pt.size=0.5,
						pt.size.scale=TRUE,
						pt.alpha=0.5, 
						label=FALSE, label.cex=0.7,
						label.adj=1) {
  if (pt.col.scale) pt.cols <- colorRampPalette(c("blue", "red"))(256)[round((sqrt(diag(x$SigmaKstar))-min(sqrt(diag(x$SigmaKstar))))/diff(range(sqrt(diag(x$SigmaKstar))))*255,0)+1]
  if (pt.size.scale) pt.size=pt.size*sqrt(diag(x$SigmaKstar))/min(sqrt(diag(x$SigmaKstar)))
  #if (label & pt.size.scale) label.adj=(max(sqrt(diag(x$SigmaKstar))/min(sqrt(diag(x$SigmaKstar)))))^(1/2)
  plot_lmks(x=x$LMKmean, segs=segs, tris=tris, pt.cols=pt.cols, pt.size=pt.size, pt.alpha=pt.alpha, label=label, label.cex=label.cex)
}

#' Plot \code{FM} Object
#' 
#' Plot form matrix.  Uses \code{cmdscale} to reconstruct landmarks for plotting.
#' @examples
#' plot_lmks(guenons$rawcoords[,,guenons$species=="patas"][,,1],
#'           segs=guenons$wireframe_segments,
#' 		     tris=guenons$triangles,
#' 		     label=TRUE)
#' FM1 <- calcFM(guenons$rawcoords[,,guenons$species=="patas"][,,1])
#' str(FM1)
#' plot(FM1,
#'      segs=guenons$wireframe_segments,
#' 		tris=guenons$triangles,
#' 		label=TRUE)
#' @export
plot.FM <- function(x, segs=NULL, tris=NULL, pt.cols="blue", pt.size=0.5, pt.alpha=0.5, label=FALSE, label.cex=0.7) {
  lmks <- cmdscale(x,3)
  colnames(lmks) <- c("X", "Y", "Z")
  rownames(lmks) <- attr(x, "Labels")
  warning("Form matrix converted to landmarks using function 'cmdscale' -\n  some axes may be reflected from original configuration.")
  plot_lmks(lmks, segs=segs, tris=tris, pt.cols=pt.cols, pt.size=pt.size, pt.alpha=pt.alpha, label=label, label.cex=label.cex)
}

#' Plot \code{SM} Object
#' 
#' Plot shape matrix.  Uses \code{cmdscale} to reconstruct landmarks for plotting.
#' @examples
#' plot_lmks(guenons$rawcoords[,,guenons$species=="patas"][,,1],
#'           segs=guenons$wireframe_segments,
#' 		     tris=guenons$triangles,
#' 		     label=TRUE)
#' SM1 <- calcSM(guenons$rawcoords[,,guenons$species=="patas"][,,1])
#' str(SM1)
#' plot(SM1,
#'      segs=guenons$wireframe_segments,
#' 		tris=guenons$triangles,
#' 		label=TRUE)
#' @export
plot.SM <- function(x, segs=NULL, tris=NULL, pt.cols="blue", pt.size=0.5, pt.alpha=0.5, label=FALSE, label.cex=0.7) {
  x <- x*attr(x, "size")
  lmks <- cmdscale(x,3)
  colnames(lmks) <- c("X", "Y", "Z")
  rownames(lmks) <- attr(x, "Labels")
  warning("Shape matrix rescaled and converted to landmarks using function 'cmdscale' -\n  some axes may be reflected from original configuration.")
  plot_lmks(lmks, segs=segs, tris=tris, pt.cols=pt.cols, pt.size=pt.size, pt.alpha=pt.alpha, label=label, label.cex=label.cex)
}
