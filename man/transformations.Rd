\name{transformations}
\alias{transformations}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Calculate similarity transformations}
\description{Calculate similarity transformations between configurations in two arrays. 
}
\usage{
transformations(Xrotated,Xoriginal)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{Xrotated}{Input k x m x n real array of the Procrustes 
transformed configurations, 
where k is the number of points, m 
is the number of dimensions, and n is the sample size. }
  \item{Xoriginal}{Input k x m x n real array of the Procrustes
original configurations,
where k is the number of points, m
is the number of dimensions, and n is the sample size. }
}
\value{A list with components
\item{translation}{The translation parameters}
\item{scale}{The scale parameters}                             
\item{rotation}{The rotation parameters} 
}
\references{Dryden, I.L. and Mardia, K.V. (1998). Statistical
Shape Analysis, Wiley, Chichester.}
\author{Ian Dryden}
\seealso{procGPA}                                      

\examples{

#2D example : female and male Gorillas (cf. Dryden and Mardia, 1998)

data(gorf.dat)

Xorig <- gorf.dat
Xrotated <- procGPA(gorf.dat)$rotated

transformations(Xrotated,Xorig)

}


\keyword{multivariate} 

