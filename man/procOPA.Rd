\name{procOPA}
\alias{procOPA}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Ordinary Procrustes analysis} 
\description{
Ordinary Procustes analysis : the matching of one configuration to  
another using translation, rotation and (possibly) scale. Reflections 
can also be included if desired. The function matches configuration B 
onto A by least squares.}                      
\usage{
procOPA(A, B, scale = TRUE, reflect = FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{A}{k x m matrix (or complex k-vector for 2D data), of 
k landmarks in m dimensions. This is the reference figure.} 
  \item{B}{k x m matrix (or complex k-vector for 2D data). This is 
the figure which is to be transformed.} 
  \item{scale}{logical indicating if scaling is required} 
  \item{reflect}{logical indicating if reflection is allowed}
}
\value{
A list with components:          
  \item{R}{The estimated rotation matrix (may be an orthogonal matrix 
if reflection is allowed)} 
  \item{s}{The estimated scale matrix}   
  \item{Ahat}{The centred configuration A}
  \item{Bhat}{The Procrustes registered configuration B}
  \item{OSS}{The ordinary Procrustes sum of squares, which is 
$\|Ahat-Bhat\|^2$}
}
\references{Dryden, I.L. and Mardia, K.V. (1998). Statistical shape 
analysis. Wiley, Chichester.} 
\author{Ian Dryden}             
\seealso{procGPA,riemdist,tpsgrid}
\examples{
data(digit3.dat)

A<-digit3.dat[,,1]
B<-digit3.dat[,,2]
ans<-procOPA(A,B) 
plotshapes(A,B,joinline=1:13)
plotshapes(ans$Ahat,ans$Bhat,joinline=1:13)
}
\keyword{multivariate}
