simdata <-
function(nchem, nassay, seed=NULL){
  
  concs <- c(0.1291468 ,  0.3913539 ,  1.1859210   ,3.5937000 , 10.8900000 , 33.0000000 ,100.0000000,200)
  nij <- length(concs)
  if(!is.null(seed)) set.seed(seed)
  N <- nchem * nassay
    
  hs <- function(x,theta) theta[1] - (theta[1]-theta[2])/(1+(x/theta[3])^theta[4])
  
  theta <- cbind(runif(N,0,10),rnorm(N,1,0.01),rnorm(N,30,10), runif(N,1,8))
  theta[,1] <- pmax(theta[,1],theta[,2])
  theta[,3] <- pmin(theta[,3],max(concs))
  theta[,3] <- pmax(theta[,3],min(concs))
  
  out <- NULL
  for(c in 1:nchem){
    for(a in 1:nassay){
      resp <- hs(concs,theta[a +nassay*(c-1),])+rnorm(nij,0,.1)
      out <- rbind(out,cbind(c,a,concs,resp))
    }
  }
  
  colnames(out) <- c("Chemical.Number","Assay.Number","Conc","Resp")
  
  return(out)
}
