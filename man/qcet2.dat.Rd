\name{qcet2.dat}
\alias{qcet2.dat}
\non_function{}
\title{Control T2 mouse vertabrae data}
\description{T2 mouse vertebrae data - control group. 6 landmarks 
in 2 dimensions, 30 individuals
}
\usage{data(qcet2.dat)}
\format{
An array of dimension 6 x 2 x 30
}
\source{
Dryden, I.L. and Mardia, K.V. (1998). 
Statistical Shape Analysis, Wiley, Chichester. p313
}
\references{
http://www.maths.nott.ac.uk/~ild/bookdata/qcet2.dat

Data from Paul O'Higgins (Hull-York Medical School) and David Johnson (Leeds)
}
\examples{
data(qcet2.dat)
k<-dim(qcet2.dat)[1]
n<-dim(qcet2.dat)[3]
par(pty="s")
plot(qcet2.dat[,,1],xlim=c(0,250),ylim=c(0,250),type="n",xlab=" ",ylab=" ")
for (i in 1:n){
points(qcet2.dat[,,i],pch=c(1:k))
}
}
\keyword{datasets}

