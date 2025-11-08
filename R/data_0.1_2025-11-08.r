#' @title GueSDat, the free Guenon Skull Database (Cardini & Elton, 2017)
#' @description Cranial landmarks for 1,315 cercopithecine crania.
#' @format ...
#' @references Cardini A, Elton S. (2017) Is there a "Wainer's rule"? Testing which sex varies most 
#'   as an example analysis using GueSDat, the free Guenon Skull Database. 
#'   \emph{Hystrix, the Italian Journal of Mammalogy}. 28:147-156. 
#'   (\href{https://doi.org/10.4404/hystrix-28.2-12139}{https://doi.org/10.4404/hystrix-28.2-12139})
#' @examples
#' plot_lmks(guenons$rawcoords[,,1],
#'      segs=guenons$wireframe_segments,
#' 		tris=guenons$triangles,
#' 		label=TRUE)
"guenons"