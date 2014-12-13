ZIPLL <-
function(dat, nitter=2000, burnin=1000, internal.knots=0, kappa=1, mu=c(5,1,3), sigma=c(3,0.1,3)){

#organize data
dat <- dat[which(apply(dat, 1, function(x) sum(1*is.na(x)) )==0),]
dat <- dat[order(dat[,1],dat[,2]),]
nassays <- length(unique(dat[,1]))
nchems <- length(unique(dat[,2]))
nij <- max(table(dat[,1],dat[,2]), na.rm=TRUE)

#find knots
if(is.null(internal.knots)) internal.knots <- 0
if(length(which(internal.knots==0)))
external.knots <- abs(as.numeric(range(log(dat[,3])) %*% c(-1,1))*10*nij)
knots <- c(min(internal.knots)-external.knots,internal.knots,max(internal.knots)+external.knots)
knots <- knots[order(knots)]



dat.out <- .C("FIT_ZIPLL",
          resp=as.double(dat[,4]),
          conc=as.double(dat[,3]),
          dimx=as.integer(c(nchems,nassays, nij)),
          n=as.integer(c(t(table(dat[,1],dat[,2])))),
          knots=as.double(knots),
          nknots=as.integer(length(knots)),
          BE=as.double(rep(0,nchems*nassays*nij)),
          BEsd=as.double(rep(0,nchems*nassays*nij)),
          AC50=as.double(rep(0,nchems*nassays)),
          TOP=as.double(rep(0,nchems*nassays)),
          AC50sd=as.double(rep(0,nchems*nassays)),
          TOPsd=as.double(rep(0,nchems*nassays)),
          active=as.double(rep(0,nchems*nassays)),
          kappa=as.double(kappa),
          smax=as.integer(nitter),
          burnin=as.integer(burnin),
          mu=as.double(mu),
          sigma=as.double(sigma)
)


dat <- cbind(dat,dat.out$BE,dat.out$BEsd)
colnames(dat)[5:6] <- c("Pred","Pred.SD")

labels <- expand.grid( x=unique(dat[,2]) , y=unique(dat[,1]) )
parms <- cbind(labels[,2],labels[,1],dat.out$AC50,dat.out$AC50sd,dat.out$TOP,dat.out$TOPsd,dat.out$active)
colnames(parms) <- c("chemical","assay","AC50","AC50.sd","Emax","Emax.sd","Pr.Active")

return(list(dat=dat, parms=parms))
}
