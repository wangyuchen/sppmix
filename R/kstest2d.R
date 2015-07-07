kstest2d <- function(x1,x2,alpha=0.05){
  if((!spatstat::is.ppp(x1)) | (!spatstat::is.ppp(x2))){
    stop("Point pattern must be of class ppp.")
  }
 n1 <- length(x1$x)
 n2 <- length(x2$x)
 pp1 <- cbind(x1$x,x1$y)
 pp2 <- cbind(x2$x,x2$y)
 KSstat=-Inf
 for (i in 1:(n1+n2)){
   if(i <= n1) {
     edge=pp1[i,]
   } else {
       edge=pp2[i-n1,]
       }
   vfcdf1=apply(ecdf2d(x1,edge), 2, sum)/n1
   vfcdf2=apply(ecdf2d(x2,edge), 2, sum)/n2
   vfdiff=abs(vfcdf2-vfcdf1)
   fksts=max(vfdiff)
   if (fksts>KSstat) {KSstat=fksts}
 }
 n <- n1*n2/(n1+n2)
 Zn <- sqrt(n) * KSstat
 Zinf <- Zn/(1-0.53 * n^(-0.9))
 pValue=2*exp(-2 * (Zinf - 0.5)^2)
 if (pValue > 0.2) {pValue = 0.2}
 H <- (pValue <= alpha)
 return(c(KSstat,pValue, H))
 if (H==1) {
   print("Reject Null, No GoF")
 } else {
     print("Cannot Reject Null, we have GoF")
   }
}

ecdf2d <- function(x,edge){
  count=cbind(as.numeric(x$x>=edge[1] & x$y>=edge[2]),as.numeric(x$x<=edge[1]
                  & x$y>=edge[2]),as.numeric(x$x<=edge[1] & x$y<=edge[2]),
                  as.numeric(x$x>=edge[1] & x$y<=edge[2]))
  return(count)
}
