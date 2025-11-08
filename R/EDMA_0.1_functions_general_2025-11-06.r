#' Geometric Mean
#' 
#' Function for calculating the geometric mean of a set of positive numbers.
#' @param x A vector of positive numbers.
#' @param na.rm a logical value indicating whether NA values should be stripped before the computation proceeds.
#' @return If all values of x are positive, the geometric mean is returned as a numeric vector of length one.  If 
#'    any values are non-positive, NA is returned.
#' @examples
#' x <- c(1, 10, 100)
#' mean(x)
#' geomean(x)
#' geomean(c(-1,x))
#' geomean(c(0,x))
#' geomean(c(NA,x))
#' geomean(c(NA,x), na.rm=TRUE)
#' @export
geomean <- function(x, na.rm=FALSE) {
  if (!is.numeric(x)) {
    warning("argument is not numeric: returning NA")
	return(NA)
  }
  if(na.rm) x <- x[!is.na(x)]
  if(!na.rm & sum(is.na(x)) > 0) return(NA)
  if(sum(x[!is.na(x)] <= 0) > 0) {
    warning("All values must be positive to calculate the geometric mean.")
    return(NA)
  }
  return(exp(mean(log(x))))
}

#' Read Landmark CSV file
#' 
#' Function for reading in a landmark CSV file of XY or XYZ coordinates for a single object.
#' @param file A character providing the name of the file to read in.
#' @return Produces a data frame preserving only the X, Y, and Z coordinates, along with any 
#'    landmark names that are included as row names in the file that is read in.  The name of the
#'    file is retained as an attribute attached to the data frame.
#' @export
readlmk.csv = function(file = NULL)  { 
  res <- read.csv(file = file, skip = 1, header = T)
  rownames(res) <- res[,"Name"]
  res <- as.matrix(res[,c("X", "Y", "Z")])
  attr(res, "filename") <- file
  return(res) 
}
