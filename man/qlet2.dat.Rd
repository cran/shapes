\name{qlet2.dat}
\alias{qlet2.dat}
\non_function{}
\title{Large T2 mouse vertabrae data}
\description{T2 mouse vertebrae data - large group. 6 landmarks 
in 2 dimensions, 23 individuals
}
\usage{data(qlet2.dat)}
\format{
An array of dimension 6 x 2 x 23
}
\source{
Dryden, I.L. and Mardia, K.V. (1998). 
Statistical Shape Analysis, Wiley, Chichester. p313
}
\references{
http://www.maths.nott.ac.uk/~ild/bookdata/qlet2.dat

Data from Paul O'Higgins (Hull-York Medical School) and David Johnson (Leeds)
}
\examples{
data(qlet2.dat)
k<-dim(qlet2.dat)[1]
n<-dim(qlet2.dat)[3]
par(pty="s")
plot(qlet2.dat[,,1],xlim=c(0,250),ylim=c(0,250),type="n",xlab=" ",ylab=" ")
for (i in 1:n){
points(qlet2.dat[,,i],pch=c(1:k))
}
}
\keyword{datasets}

