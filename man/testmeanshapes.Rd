\name{testmeanshapes}
\alias{testmeanshapes}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Tests for mean shape difference}
\description{
  Carries out Hotelling's $T^2$ or Goodall's $F$ test to examine differences in mean shape
between two independent populations, for $m>=2$ dimensional data. 
The procedure uses complex eigenanalysis for $m=2$ and iterative 
Generalised Procrustes Analysis for $m>2$ dimensions. 
}
\usage{
testmeanshapes(A, B, Hotelling = TRUE, tol1 = 1e-05, tol2 = 1e-05)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{A}{The random sample for group 1: k x m x n1 array of data, where 
k is the number of landmarks, m is dimension and n1 is the sample size}
  \item{B}{The random sample for group 2: k x m x n2 array of data, where 
k is the number of landmarks, m is dimension and n2 is the sample size}
  \item{Hotelling}{Logical. If TRUE then carry out Hotelling's $T^2$ test, 
  if FALSE carry out Goodall's $F$ test}
    \item{tol1}{Tolerance for optimal rotation for the iterative 
algorithm ($m>2$): tolerance on the mean sum of squares between successive iterations 
(depends on scale of objects)}
  \item{tol2}{tolerance for rescale/rotation step for the iterative 
algorithm ($m>2$): tolerance on the Riemannian shape distance 
between successive mean shapes}
}

\value{
  A list with components 
   \item{F}{the F statistic}
   \item{df1 and df2}{degrees of freedom of the F statistic}
   \item{pval}{p-value for the test} 
  \item{Tsq}{the $T^2$ statistic (if Hotelling)}
  \item{T.df1 and T.df2}{degrees of freedom of the $T^2$ statistic (if Hotelling)}
  \item{Tsq.partition}{the $T^2$ statistic partitioned into 
contributions from each of the pooled principal components (if Hotelling)} 
 
  \item{F.partition}{the F statistic partitioned into 
contributions from each of the pooled principal components (if Hotelling)} 
}

\references{Dryden, I.L. and Mardia, K.V. (1998) Statistical Shape Analysis, 
Wiley, Chichester. Chapter 7. 
 
Dryden and Mardia (1993) Multivariate shape analysis. Sankhya A, 55:460-480.

Goodall, C. R. (1991). Procrustes methods in the statistical analysis 
of shape (with discussion). Journal of the Royal Statistical Society, 
Series B, 53: 285-339. 
} 




\author{Ian Dryden}

\seealso{procGPA}

\examples{

#2D example : female and male Gorillas (cf. Dryden and Mardia, 1998)

data(gorf.dat)
data(gorm.dat)

#Hotelling's Tsq test
test1<-testmeanshapes(gorf.dat,gorm.dat)

#Goodall's isotropic test
test2<-testmeanshapes(gorf.dat,gorm.dat,Hotelling=FALSE)

}




\keyword{multivariate}% at least one, from doc/KEYWORDS
