\name{macf.dat}
\alias{macf.dat}
\non_function{}
\title{Female macaque data}
\description{Female macaque skull data. 7 landmarks 
in 3 dimensions, 9 individuals
}
\usage{data(macf.dat)}
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
data(macf.dat)
n<-dim(macf.dat)[3]
par(pty="s")
plot(macf.dat[,1:2,1],xlim=c(0,180),ylim=c(-50,130),type="n",xlab=" ",ylab=" ")
for (i in 1:n){
points(macf.dat[,1:2,i],pch=c(1:7))
}
testmeanshapes(macf.dat,macm.dat)
}
\keyword{datasets}

