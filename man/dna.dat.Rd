\name{dna.dat}
\alias{dna.dat}
\non_function{}
\title{DNA data}
\description{DNA data. 22 landmarks 
in 3 dimensions, 47 individuals
}
\usage{data(dna.dat)}
\format{
An array of dimension 22 x 3 x 47
}
\source{
}
\references{

}
\examples{
data(dna.dat)
par(mfrow=c(2,2))
k<-dim(dna.dat)[1]
n<-dim(dna.dat)[3]
par(pty="s")
plot(dna.dat[,1:2,1],xlim=c(0,50),ylim=c(0,50),type="n",xlab=" ",ylab=" ")
for (i in 1:n){
lines(dna.dat[,1:2,i])
}
plot(dna.dat[,c(1,3),1],xlim=c(0,50),ylim=c(0,50),type="n",xlab=" ",ylab=" ")
for (i in 1:n){
lines(dna.dat[,c(1,3),i])
}
plot(dna.dat[,c(2,3),1],xlim=c(0,50),ylim=c(0,50),type="n",xlab=" ",ylab=" ")
for (i in 1:n){
lines(dna.dat[,c(2,3),i])
}
}
\keyword{datasets}

