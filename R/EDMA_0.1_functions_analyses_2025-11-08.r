#' Calculate Form Matrix (FM)
#' 
#' Generate the Euclidean distance form matrix (FM) for a single object.
#' @param x A matrix with landmarks as rows and coordinates (X, Y, Y) as columns.
#' @return A pairwise distance matrix (of class \code{dist} and \code{FM}) for Euclidean distance 
#'   between all pairs of landmarks in the data set.  Overall object size is calculated as the 
#'   geometric mean of all distances and is returned as the attribute \code{size} in the \code{FM} 
#'   object.
#' @examples
#' FM1 <- calcFM(guenons$rawcoords[,,1])
#' str(FM1)
#' plot(FM1)
#' @export
calcFM <- function(x) {
  d <- dist(x)
  attr(d, "size") <- geomean(as.numeric(d))
  class(d) <- c(class(d), "FM")
  return(d)
}

#' Calculate Shape Matrix (SM)
#' 
#' Generate the Euclidean distance shape matrix (SM) for a single object.
#' @param x A matrix with landmarks as rows and coordinates (X, Y, Y) as columns.
#' @return A pairwise distance matrix (of class \code{dist} and \code{SM}) for Euclidean distance 
#'   between all pairs of landmarks in the data set, scaled by size.  Overall object size is 
#'   calculated as the geometric mean of all distances and is returned as the attribute \code{size} 
#'   in the \code{SM} object.
#' @examples
#' SM1 <- calcSM(guenons$rawcoords[,,1])
#' str(SM1)
#' plot(SM1)
#' @export
calcSM <- function(x) {
  d <- calcFM(x)
  d <- d/attr(d, "size")
  class(d)[class(d)=="FM"] <- "SM"
  return(d)
}

#' Calculate Shape Difference Between Two Objects
#' 
#' Calculate the shape dissimilarity between two objects sharing a common set of landmarks using the 
#'   \emph{slr} metric developed by Gordon and Wood (2013).
#' @param x A \code{FM} or \code{SM} object.
#' @param y A \code{FM} or \code{SM} object.
#' @return A measure of shape dissimilarity, measure as the standard deviation of the natural log 
#'   of ratios for all interlandmark distances, where ratios are calculated as (specimen \code{x})/(specimen \code{y}).
#' @examples
#' shapediff(x=calcFM(guenons$rawcoords[,,1]), y=calcFM(guenons$rawcoords[,,2]))
#' @export
shapediff <- function(x, y) {
  # x and y must be Euclidean distance matrices describing individual specimens
  # shape difference calculated using slr from Gordon and Wood (2013)
  slr <- sd(log(x/y))
  return(slr)
}

#' Principal Coordinates Analysis Ordination for Multiple Objects
#' 
#' Perform a Principal Coordinates Analysis (PCoA) quantifying shape variation for a set of objects 
#'   with a common set of landmarks.  Performs \code{cmdscale()} on a pairwise distance matrix 
#'   for all objects (calculated using \code{shapediff}).
#' @param x An array of landmarks for multiple specimens, with the first dimension corresponding to 
#'   landmarks, the second corresponding to X, Y, and Z, and the third corresponding to specimens.
#' @return A list return by \code{cmdscale} containing the principal coordinates and associated 
#'   eigenvalues.
#' @examples
#' PCoA <- EDMApcoa(guenons$rawcoords[,,guenons$genus %in% c("Erythrocebus", "Miopithecus")])
#' plot(PCoA$points, pch=21, bg=guenons$genus[guenons$genus %in% c("Erythrocebus", "Miopithecus")])
#' @export
EDMApcoa <- function(x) {
  # x is an array of landmarks
  FMs <- apply(x, MARGIN=3, FUN=calcFM) # this collapses each individual distance matrix into a column for each specimen
  n <- ncol(FMs)
  slrs <- matrix(NA, n, n)
  rownames(slrs) <- colnames(FMs)
  colnames(slrs) <- colnames(FMs)
  diag(slrs) <- 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      slrs[i,j] <- shapediff(FMs[,i], FMs[,j])
      slrs[j,i] <- slrs[i,j]
    }
  }
  slrs <- as.dist(slrs)
  PCoA <- cmdscale(slrs, k=n-1, eig=T)
  k <- sum(PCoA$eig >= 0)
  if (k > n-1) k <- n-1
  PCoA$points <- PCoA$points[,1:k]
  PCoA$eig <- PCoA$eig[1:k]
  return(PCoA)
}

#' Visualize Influential Landmarks
#' 
#' Calculates a measure of spatial variability for each landmark in relation to all other landmarks 
#'   for a pair of objects.  This is (mostly) the visualization method of Cole and Richtsmeier (1998), but 
#'   it presents the log of ratios as opposed to raw ratios, and shows the ratio from shape matrices as 
#'   well as form matrices.
#' @param x An array of landmarks for multiple specimens, with the first dimension corresponding to 
#'   landmarks, the second corresponding to X, Y, and Z, and the third corresponding to specimens.
#' @param highlight An integer value specifying the number of a landmark to highlight in plots. If 
#'   provided, all logged ratios for interlandmark distances involving that landmark will be highlighted 
#'   in the plot.  Defaults to \code{NULL}.
#' @param pt.bg A color value to use in plotting logged ratios. 
#' @param hl.pch An integer specifying the plotting character to use for highlighted points when a 
#'   value is supplied to \code{highlight}. 
#' @param hl.pch A color value to use for highlighted points when a value is supplied to \code{highlight}. 
#' @examples
#' influentiallmk(A=guenons$rawcoords[,,guenons$genus=="Erythrocebus"][,,1],
#'                B=guenons$rawcoords[,,guenons$genus=="Miopithecus"][,,1],
#'                highlight=57)
#' @export
influentiallmk <- function(A, B, highlight=NULL, pt.bg="#00000020", hl.pch=4, hl.col="#FF0000FF") {
  # A and B are landmark matrices
  plt=TRUE
  FMA <- calcFM(A)
  FMB <- calcFM(B)
  SMA <- calcSM(A)
  SMB <- calcSM(B)
  lmknames <- attr(FMA, "Labels")
  nlmk <- length(lmknames)
  nspec <- 2
  colpt <- NULL
  rowpt <- NULL
  for (i in 1:nlmk) colpt <- c(colpt, rep(i, nlmk-i))
  for (i in 2:nlmk) rowpt <- c(rowpt, i:nlmk)
  segs <- data.frame(lmk1=colpt, lmk2=rowpt)
  FMs <- cbind(as.numeric(FMA), as.numeric(FMB)) # this collapses each individual distance matrix into a column for each specimen
  SMs <- cbind(as.numeric(SMA), as.numeric(SMB)) # this collapses each individual distance matrix into a column for each specimen
  rownames(FMs) <- paste0(segs$lmk1, "_", segs$lmk2)
  rownames(SMs) <- paste0(segs$lmk1, "_", segs$lmk2)
  logFMs <- log(FMs)
  logSMs <- log(SMs)
  rm(FMs, SMs)
  # calculate logged shape ratios (as difference in logged shape)
  ads <- combn(nspec, 2) # this is half of the ways to pair them - will need to duplicate final results and multiply duplicate by -1 for other half
  logFMrats <- apply(ads, 2, FUN=function(x) return(logFMs[,x[1]]-logFMs[,x[2]]))
  logSMrats <- apply(ads, 2, FUN=function(x) return(logSMs[,x[1]]-logSMs[,x[2]]))
  rm(logFMs, logSMs)
  #logSMrats <- cbind(logSMrats, -1*logSMrats) # add logged ratios for same pairs with numerator and denominator swapped
  ncomp <- ncol(logSMrats)
  logFMratvec <- rep(as.numeric(logFMrats), times=2)
  logSMratvec <- rep(as.numeric(logSMrats), times=2)
  segslong <- c(rep(segs$lmk1, times=ncomp), rep(segs$lmk2, times=ncomp))
  segslongDF <- data.frame(lmk1=rep(rep(segs$lmk1, times=ncomp), times=2), lmk2=rep(rep(segs$lmk2, times=ncomp)))
  if (plt) {
	if (!is.null(highlight)) {
	par(mfrow=c(2,1))
    plot(x=segslong[segslongDF$lmk1!=highlight & segslongDF$lmk2!=highlight],
	     y=logFMratvec[segslongDF$lmk1!=highlight & segslongDF$lmk2!=highlight],
		 pch=21, col=NULL, bg=pt.bg,
         xlim=range(segslong), ylim=range(logSMratvec), xlab="landmark", ylab="log form ratio")
    abline(h=0, col="#00000060")
	points(x=segslong[segslongDF$lmk1==highlight | segslongDF$lmk2==highlight],
	       y=logFMratvec[segslongDF$lmk1==highlight | segslongDF$lmk2==highlight],
		   pch=hl.pch, col=hl.col)
    plot(x=segslong[segslongDF$lmk1!=highlight & segslongDF$lmk2!=highlight],
	     y=logSMratvec[segslongDF$lmk1!=highlight & segslongDF$lmk2!=highlight],
		 pch=21, col=NULL, bg=pt.bg,
         xlim=range(segslong), ylim=range(logSMratvec), xlab="landmark", ylab="log shape ratio")
    abline(h=0, col="#00000060")
	points(x=segslong[segslongDF$lmk1==highlight | segslongDF$lmk2==highlight],
	       y=logSMratvec[segslongDF$lmk1==highlight | segslongDF$lmk2==highlight],
		   pch=hl.pch, col=hl.col)
	}
	else{
	par(mfrow=c(2,1))
    plot(x=segslong, y=logFMratvec, pch=21, col=NULL, bg=pt.bg,
         xlab="landmark", ylab="log form ratio")
    abline(h=0, col="#00000060")
    plot(x=segslong, y=logSMratvec, pch=21, col=NULL, bg=pt.bg,
         xlab="landmark", ylab="log shape ratio")
    abline(h=0, col="#00000060")
	}
  }
  #out <- data.frame(landmark=lmknames,
  #                  #mean_ln=tapply(logSMratvec, INDEX=segslong, FUN=mean),
  #				     sd_ln=tapply(logSMratvec, INDEX=segslong, FUN=sd))
  #return(out)
}

#' Estimate Mean Form and Landmark Covariance
#' 
#' Estimates the mean form, mean landmark set, and covariance structure of lanmarks for a set of objects  
#'   following the algorithms provided in Lele and Richstmeier (2000).
#' @param x An array of landmarks for multiple specimens, with the first dimension corresponding to 
#'   landmarks, the second corresponding to X, Y, and Z, and the third corresponding to specimens.
#' @return A data frame providing the standard deviation of logged ratios associated with each landmark. 
#' @examples
#' FMmean <- meanform(guenons$rawcoords[,,guenons$species=="patas"])
#' str(FMmean)
#' plot(FMmean)
#' @export
meanform <- function(x) {
  lmknames <- attr(calcFM(x[,,1]), "Labels")
  nlmk <- dim(x)[1]
  ndim <- dim(x)[2] 
  nspec <- dim(x)[3]
  colpt <- NULL
  rowpt <- NULL
  for (i in 1:nlmk) colpt <- c(colpt, rep(i, nlmk-i))
  for (i in 2:nlmk) rowpt <- c(rowpt, i:nlmk)
  segs <- data.frame(lmk1=colpt, lmk2=rowpt)
  # Calculate mean form
  FMs <- apply(x, MARGIN=3, FUN=calcFM) # this collapses each individual distance matrix into a column for each specimen
  rownames(FMs) <- paste0(segs$lmk1, "_", segs$lmk2)
  EMs <- FMs^2
  EMmean <- apply(EMs, 1, mean)
  EMvar <- apply(EMs, 1, function(x) return(sum((x-mean(x))^2)/length(x))) # note: uses n rather than (n-1) as denominator
  # Convert back to matrices
  EMmeanmat <- matrix(0, nlmk, nlmk)
  EMvarmat <- matrix(0, nlmk, nlmk)
  counter <- 0
  for (j in 1:(nlmk-1)) {
    for(i in (j+1):nlmk) {
      counter <- counter+1
	  EMmeanmat[i,j] <- EMmean[counter]
	  EMmeanmat[j,i] <- EMmean[counter]
	  EMvarmat[i,j] <- EMvar[counter]
	  EMvarmat[j,i] <- EMvar[counter]
    }
  }
  # Estimate mean landmark shape, then calculate mean form
  H <- matrix(-1/nlmk, nlmk, nlmk)
  diag(H) <- 1-(1/nlmk)
  Bmean <- -0.5*H%*%EMmeanmat%*%t(H)
  Beig <- eigen(Bmean)
  LMKmean <- cbind(X=Beig$vectors[,1]*sqrt(Beig$values[1]),
                   Y=Beig$vectors[,2]*sqrt(Beig$values[2]),
				   Z=Beig$vectors[,3]*sqrt(Beig$values[3]))
  rownames(LMKmean) <- lmknames
  # Build in check to see if any axis needs to be inverted to preserve original orientation (i.e., check for reflection)
  FMmean <- calcFM(LMKmean)
  # Calculate SigmaKstar
  # step 1
  xcentered <- apply(x, 3, scale, scale=F)
  xcentered <- array(xcentered, c(nlmk, ndim, nspec))
  # step 2
  SKstarvals <- apply(apply(xcentered, 3, function(x) return((x%*%t(x))-(LMKmean%*%t(LMKmean)))), 1, mean)/ndim # use to fill matrix
  SigmaKstar <- matrix(SKstarvals, nlmk, nlmk)
  # step 3
  SKstareig <- eigen(SigmaKstar)
  SKstareig$values[SKstareig$values<0] <- 0
  SigmaKstar <- SKstareig$vectors%*%diag(SKstareig$values)%*%t(SKstareig$vectors)
  out <- list(LMKmean=LMKmean, FMmean=FMmean, SigmaKstar=SigmaKstar)
  class(out) <- "FMmean"
  return(out)
}



