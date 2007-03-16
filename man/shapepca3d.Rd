\name{shapepca3d}
\alias{shapepca3d}
\title{Plot PCs for 3D data}
\description{Plot the effect of principal components for 3D data}
\usage{shapepca3d(procreg,pcno=c(1,2,3))}
\arguments{
\item{procreg}{Procrustes analysis output, as from procGPA()}
\item{pcno}{The numbers of the PCs to be plotted}
}
\value{
None
}
\references{Dryden, I.L. and Mardia, K.V. (1998) 
Statistical Shape Analysis. Wiley, Chichester. 
Chapter 5.

}
\author{Ian Dryden}

\examples{
data(macm.dat)
mac<-procGPA(macm.dat)
shapepca3d(mac,pcno=1)
}

\keyword{multivariate}
