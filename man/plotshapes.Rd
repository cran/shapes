\name{plotshapes}
\alias{plotshapes}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Plot 2D shapes and configurations}            
\description{
Plots 2D shapes and configurations. Either one or two groups of observations
can be plotted on the same scale. 
}
\usage{
plotshapes(A, B = 0, joinline = c(1, 1))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{A}{k x 2 x n array, or k x 2 matrix for first group} 
  \item{B}{k x 2 x n array, or k x 2 matrix for 2nd group (can be missing)}
  \item{joinline}{A vector stating which landmarks are joined up by lines,
  e.g. joinline=c(1:n,1) will start at landmark 1, join to 2, ..., join to n, 
  then re-join to landmark 1.}
} 
\value{
Just graphical output
}
\author{Ian Dryden}            
\seealso{shapepca,tpsgrid}
\examples{
data(gorf.dat)
data(gorm.dat)
plotshapes(gorf.dat)
plotshapes(gorf.dat,gorm.dat,joinline=c(1,6,7,8,2,3,4,5,1))
}
\keyword{hplot}
\keyword{multivariate}
