\name{qset2.dat}
\alias{qset2.dat}
\non_function{}
\title{Small T2 mouse vertabrae data}
\description{T2 mouse vertebrae data - small group. 6 landmarks 
in 2 dimensions, 23 individuals
}
\usage{data(qset2.dat)}
\format{
An array of dimension 6 x 2 x 23
}
\source{
Dryden, I.L. and Mardia, K.V. (1998). 
Statistical Shape Analysis, Wiley, Chichester. p313
}
\references{
http://www.maths.nott.ac.uk/~ild/bookdata/qset2.dat

Data from Paul O'Higgins (Hull-York Medical School) and David Johnson (Leeds)
}
\examples{
data(qset2.dat)
k<-dim(qset2.dat)[1]
n<-dim(qset2.dat)[3]
par(pty="s")
plot(qset2.dat[,,1],xlim=c(0,250),ylim=c(0,250),type="n",xlab=" ",ylab=" ")
for (i in 1:n){
points(qset2.dat[,,i],pch=c(1:k))
}
}
\keyword{datasets}

