#' Estimate intensity surface by using non-parametric method.
#'
#' Using Epanechnikov kernel to estimate intensity surface.
#'
#' @param pattern An object of class \code{\link[spatstat]{ppp}}. It should be a
#'  2 dimesional spatial point pattern.
#' @param win Object of class \code{\link[spatstat]{owin}}.
#' @param h Kernel bandwith.  \code{h} should be a postive number.
#' @param L Number of grid for x and y axis.
#' @param kernel Kernel used to estimate intensity surface.  Currently, only
#'  support Epanechnikov kernel.
#' @param edgecorrect Logical flag indicating whether to use edge-correction in
#'  estimating intensity surface. The default is TRUE.
#' @param truncate Locgical flag indicating whether truncation is used for
#'  \code{pattern}. The default is TRUE.
#'
#' @return An object of class \code{\link[spatstat]{im}}.
#' @export
#' @examples
# generate a point pattern
#' if (require(spatstat)){
#'   mix1 <- rnormmix(3, .01, 5, square(5))
#'   pattern1 <- rsppmix(30, mix1, square(5))
#' }
#'
#' # estimate and plot the estimated intensity surface
#' if (require(spatstat)){
#'   surf1 <- est_intensity_np(pattern1, win=square(5), h=0.05, L=100)
#'   plot(surf1)
#' }
est_intensity_np <- function(pattern, win, h, L=10, kernel=c("Epanechnikov"),
                             edgecorrect=TRUE, truncate=TRUE){
  kernel <- match.arg(kernel)
  if (truncate==TRUE) {
    pattern <- pattern[spatstat::inside.owin(pattern$x,pattern$y,win)]
  }
  x <- seq(win$xrange[1],win$xrange[2],length.out = L)
  y <- seq(win$yrange[1],win$yrange[2],length.out = L)
  lenx <- length(x)
  leny <- length(y)
  intensity <- matrix(0,lenx,leny)
  loc <- expand.grid(x,y)
  pp <- cbind(pattern$x,pattern$y)
  index <- expand.grid(1:pattern$n, 1:(lenx*leny))
  xy <- loc[index[,2],]-pp[index[,1],]
  if (edgecorrect==TRUE){
    LL <- 20
    ax <- seq(win$xrange[1],win$xrange[2], length.out = LL)
    ay <- seq(win$yrange[1],win$yrange[2], length.out = LL)
    alenx <- length(ax)
    aleny <- length(ay)
    area <- (ax[2] - ax[1]) * (ay[2] - ay[1])
    centers <- expand.grid(x_center=(ax[1:(alenx-1)] + ax[2:alenx])/2,
                           y_center=(ay[1:(aleny-1)] + ay[2:aleny])/2)
    index.mass <- expand.grid(1:((alenx-1)*(aleny-1)),1:(lenx*leny))
    dist.centers <- loc[index.mass[,2],]-centers[index.mass[,1],]
    if(kernel=="Epanechnikov"){
      diff <- rowSums((dist.centers)^2)
      mass <- ifelse(diff < 1, 4*(1-diff) / (2*pi), 0)
      edgecor <- area * aggregate(mass,list(index.mass[,2]),sum) / (h^2)
      quad <- rowSums((xy/h)^2)
      val <- ifelse(quad < 1, 4*(1 - quad) / (2*pi), 0)
      intensity <- aggregate(val, list(index[, 2]), sum) / (h^2*edgecor)
    }
  } else {
    if(kernel=="Epanechnikov"){
      quad <- rowSums((xy/h)^2)
      val <- ifelse(quad < 1, 4*(1 - quad) / (2*pi), 0)
      intensity <- aggregate(val, list(index[, 2]), sum) / h^2
    }
  }
  return(spatstat::im(matrix(intensity[,2],lenx,leny,byrow=T),x,y))
}

