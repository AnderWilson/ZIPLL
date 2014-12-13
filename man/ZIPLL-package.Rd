\name{ZIPLL-package}
\alias{ZIPLL-package}
\docType{package}
\title{
ZIPLL
}
\description{
This package fits the ZIPLL hierarchical dose-response model.
}
\details{
\tabular{ll}{
Package: \tab ZIPLL\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-05-28\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}
This package fits the ZIPLL hierarchical dose-response model.
}
\author{
Ander Wilson

Maintainer: Ander Wilson <ander_wilson@ncsu.edu>
}
\references{
forthcoming
}
\keyword{ package }
\examples{

#simulate data
dat <- simdata(nchem=10, nassay=4, seed=1234)

#fit ZIPLL
fit <- ZIPLL(dat)

#plot a curve
dose<-fit$dat[,3]
response<-fit$dat[,4]
post.mn<-fit$dat[,5]
post.sd<-fit$dat[,6]

these<-dat[,1]==6 & dat[,2]==2
plot(dose[these],response[these],xlab="Log dose",ylab="Response", log="x")
lines(dose[these],post.mn[these])
lines(dose[these],post.mn[these]-2*post.sd[these],lty=2)
lines(dose[these],post.mn[these]+2*post.sd[these],lty=2)

}
