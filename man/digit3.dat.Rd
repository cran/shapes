\name{digit3.dat}
\alias{digit3.dat}
\non_function{}
\title{Digit 3 data}
\description{Handwritten digit `3' data. 13 landmarks 
in 2 dimensions, 30 individuals
}
\usage{data(digit3.dat)}
\format{
An array of dimension 6 x 2 x 30
}
\source{
Dryden, I.L. and Mardia, K.V. (1998). 
Statistical Shape Analysis, Wiley, Chichester. p318
}
\references{
http://www.maths.nott.ac.uk/~ild/bookdata/digit3.dat

Data from Cath Anderson
}
\examples{
data(digit3.dat)
k<-dim(digit3.dat)[1]
n<-dim(digit3.dat)[3]
digit3.dat[,2,]<- -digit3.dat[,2,]+50
par(pty="s")
plot(digit3.dat[,,1],xlim=c(0,50),ylim=c(0,50),type="n",xlab=" ",ylab=" ")
for (i in 1:n){
lines(digit3.dat[,,i])
}
}
\keyword{datasets}

