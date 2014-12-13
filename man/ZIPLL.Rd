\name{ZIPLL}
\alias{ZIPLL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Zero-Inflated Piecewise Log-Logistic
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This implements the ZIPLL hierarchical dose-response estiamation.
}
\usage{
ZIPLL(dat, nitter = 2000, burnin = 1000, internal.knots = 0, kappa = 1, mu = c(5, 1, 3), sigma = c(3, 0.1, 3))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
This is a matrix with columns (Chemical ID, Assay ID, Concentration, Response).  Each column row corresponds to a unique dose-response pair for a single chemical-assay combination.  The chemical and assay IDs should be integers.
}
  \item{nitter}{
%%     ~~Describe \code{nitter} here~~
Total number of MCMC itterations.
}
  \item{burnin}{
%%     ~~Describe \code{burnin} here~~
Number of itteration discarded for burnin.  No data from these will be saved.
}
  \item{internal.knots}{
%%     ~~Describe \code{internal.knots} here~~
The location of internal knots.  The default is one internal knot at 0.
}
  \item{kappa}{
%%     ~~Describe \code{kappa} here~~
This is the minimum difference between the Emin and Emax requried for a response to be considered active.
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
The prior mean for the Emax, Emin, and AC50. Note that the AC50 is the transformed value.
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
The prior variance for the Emax, Emin, and AC50.
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{dat}{
This returns the data with the posterior predicted value and posterior standard deviation for each design point.
}

\item{parms}{
This returns posterior mean and posterior variance for the AC50 and Emax and the probability of an active response.  For chemicals-assay combinations that are not active the posterior mean and variance is drawing from the prior.
}
}
\references{
%% ~put references to the literature/web site here ~
forthcoming
}
\author{
%%  ~~who you are~~
Ander Wilson
}

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
