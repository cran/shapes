\name{gorm.dat}
\alias{gorm.dat}
\non_function{}
\title{Male gorilla data}
\description{Male gorilla skull data. 8 landmarks 
in 2 dimensions, 29 individuals
}
\usage{data(gorm.dat)}
\format{
An array of dimension 8 x 2 x 29
}
\source{
Dryden, I.L. and Mardia, K.V. (1998). 
Statistical Shape Analysis, Wiley, Chichester. p317 
}
\references{
http://www.maths.nott.ac.uk/~ild/bookdata/gorm.dat

Data from Paul O'Higgins (Hull-York Medical School)
}
\examples{
data(gorm.dat)
n<-dim(gorm.dat)[3]
par(pty="s")
plot(gorm.dat[,,1],xlim=c(-50,250),ylim=c(-50,250),type="n",xlab=" ",ylab=" ")
for (i in 1:n){
points(gorm.dat[,,i],pch=c(1:8))
}
}
\keyword{datasets}

