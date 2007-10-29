\name{sooty.dat}
\alias{sooty.dat}
\title{Sooty mangabey data}
\description{Sooty mangabey data skull data. 12 landmarks 
in 2 dimensions, 2 individuals (juvenile and adult)
}
\usage{data(sooty.dat)}
\format{
An array of dimension 12 x 2 x 2
}
\source{
Dryden, I.L. and Mardia, K.V. (1998). 
Statistical Shape Analysis, Wiley, Chichester. p17, 42 
}
\references{
Data from Paul O'Higgins (Hull-York Medical School)
}
\examples{
data(sooty.dat)
plotshapes(sooty.dat,joinline=c(1:12,1))
}
\keyword{datasets}

