est_disc_mark <- function(pp, L, burnin, Rs) {
  # test if this pp has mark
  if(missing(pp$marks)) {
    stop("This point pattern doesn't include marks.")
  }
  # get the levels of mark
  marks <- as.numeric(levels(pp$marks))
  nmarks <- length(marks)
  n <- pp$n
  # calculate neighborhood
  # Rlen <- length(Rs)
  r <- hyper[2]
  Rlen <- 10
  Rmin <- min(Rs)
  Rmax <- max(Rs)
  RR <- rep(0, 10)
  for(j in 1:10) {
    RR[j] <- Rmin + j*(Rmax - Rmin)/10
  }
  if(r !=0) {
    RR[1] <- r
    Rlen <- 1
  }

  for(rr in 1:Rlen) {
  sig0 <- hyper[1]*diag(2)
  sum1 <- 0
  MHjump1 <- 0
  MHjump2 <- 0
  countksi <- rep(0, nmarks)
  stddata <- cbind(pp$x, pp$y)
  for (i in 1:n) {
    # stddata[i, ] <-
  }
  }
}
