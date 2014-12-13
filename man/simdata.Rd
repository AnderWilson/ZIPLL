\name{simdata}
\alias{simdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Simulate dose-response data.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Simulates a simple data set.
}
\usage{
simdata(nchem, nassay, seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nchem}{
%%     ~~Describe \code{nchem} here~~
Number of chemicals.
}
  \item{nassay}{
%%     ~~Describe \code{nassay} here~~
Number of assays.
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
Random seed.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  \item{dat}{
A data set of read to use in the ZIPLL function with chemical ID, assay ID, concentration, and response.
}
}
\author{
%%  ~~who you are~~
Ander Wilson
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
