\name{gorf.dat}
\alias{gorf.dat}
\non_function{}
\title{Female gorilla data}
\description{Female gorilla skull data. 8 landmarks 
in 2 dimensions, 30 individuals
}
\usage{data(gorf.dat)}
\format{
An array of dimension 8 x 2 x 30
}
\source{
Dryden, I.L. and Mardia, K.V. (1998). 
Statistical Shape Analysis, Wiley, Chichester. p317 
}
\references{
http://www.maths.nott.ac.uk/~ild/bookdata/gorf.dat

Data from Paul O'Higgins (Hull-York Medical School)
}
\examples{
data(gorf.dat)
n<-dim(gorf.dat)[3]
par(pty="s")
plot(gorf.dat[,,1],xlim=c(-50,250),ylim=c(-50,250),type="n",xlab=" ",ylab=" ")
for (i in 1:n){
points(gorf.dat[,,i],pch=c(1:8))
}
}
\keyword{datasets}

