\name{macm.dat}
\alias{macm.dat}
\non_function{}
\title{Male macaque data}
\description{Male macaque skull data. 7 landmarks 
in 3 dimensions, 9 individuals
}
\usage{data(macm.dat)}
\format{
An array of dimension 7 x 3 x 9
}
\source{
Dryden, I.L. and Mardia, K.V. (1998). 
Statistical Shape Analysis, Wiley, Chichester. p16
}
\references{
Data from Paul O'Higgins (Hull-York Medical School)
}
\examples{
data(macm.dat)
n<-dim(macm.dat)[3]
par(pty="s")
plot(macm.dat[,1:2,1],xlim=c(0,180),ylim=c(-50,130),type="n",xlab=" ",ylab=" ")
for (i in 1:n){
points(macm.dat[,1:2,i],pch=c(1:7))
}
testmeanshapes(macm.dat,macf.dat)
}
\keyword{datasets}

