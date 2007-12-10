\name{tpsgrid}
\alias{tpsgrid}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Thin-plate spline transformation grids}
\description{Thin-plate spline transformation grids from one planar object
to another. A square grid on the first object is deformed smoothly using a pair of 
thin-plate splines to a curved grid on the second object. 
}
\usage{
tpsgrid(TT, YY, xbegin=-999, ybegin=-999, xwidth=-999, opt=2, ext, ngrid, cex=1, pch, col=2)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{TT}{First object (source): (k x 2 matrix)}
  \item{YY}{Second object (target): (k x 2 matrix)}
  \item{xbegin}{lowest x value for plot: if -999 then a value is determined}
  \item{ybegin}{lowest y value for plot: if -999 then a value is determined }
  \item{xwidth}{width of plot: if -999 then a value is determined}
  \item{opt}{Option 1: (just deformed grid on YY is displayed), 
  option 2: both grids are displayed}
  \item{ext}{Amount that grid extends over object.}
  \item{ngrid}{Number of grid points: size is ngrid * (ngrid -1)}
  \item{cex}{Point size}
  \item{pch}{Point symbol}
  \item{col}{Point colour}
}
\details{
Note - only for $m=2$ planar objects
}
\value{
No returned value
}
\references{
Bookstein, F.L. (1989). Principal  warps: thin-plate
       splines and  the decomposition  of  deformations,
       IEEE Transactions on Pattern Analysis and Machine Intelligence,
   11, 567--585. 

Dryden, I.L. and Mardia, K.V. (1998) Statistical Shape 
Analysis. Wiley, Chichester. Chapter 10.
}
\author{Ian Dryden}

\seealso{procGPA}

\examples{
data(gorf.dat)
data(gorm.dat)

#TPS grid with shape change exaggerated (2x)
gorf<-procGPA(gorf.dat)
gorm<-procGPA(gorm.dat)
mag<-2
TT<-gorf$mshape
YY<-gorm$mshape
par(mfrow=c(1,2))
YY<-TT+(YY-TT)*mag
tpsgrid(TT,YY,-0.6,-0.6,1.2,2,0.1,22) 
title("TPS grid: Female mean (left) to Male mean (right)")   

}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{hplot}
