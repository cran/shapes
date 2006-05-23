\name{resampletest}
\alias{resampletest}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Tests for 2D mean shape difference, including boostrap and permutation tests}
\description{
Carries out tests to examine differences in mean shape
between two independent populations, for $m=2$ dimensional data.
Four test statistics are calculated:

lambda : the asymptotically pivotal statistic $lambda_{min}$ from Amaral et al. (2006), equ.(16)

H : Hotelling $T^2$ statistic (see Amaral et al., 2006, equ.(26), Dryden and Mardia, 1998, equ.(7.4))

J : James' statistic (see Amaral et al., 2006, equ.(27) )

G : Goodall's F statistic (see Amaral et al., 2006, equ.(28), Dryden and Mardia, 1998, equ.(7.9))

p-values are given based on resampling (either a bootstrap test or a permutation test) as well
as the usual table based p-values.

}
\usage{
resampletest(A, B, resamples = 200, permutation = FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{A}{The random sample for group 1: k x 2 x n1 array of data, where
k is the number of landmarks and n1 is the sample size. (Alternatively a k x n1 complex matrix)}
  \item{B}{The random sample for group 3: k x 2 x n2 array of data, where
k is the number of landmarks and n2 is the sample size. (Alternatively a k x n2 complex matrix)}
  \item{resamples}{Integer. The number of resampling iterations.}
    \item{permutation}{Logical. If permutation = TRUE then permutation resampling is
    carried without replacement from the pooled samples. If permutation = FALSE then bootstrap resampling is carried out
    in each group (sampling with replacement).}
}

\value{
  A list with components
   \item{lambda}{$\lambda_{min}$ statistic}
   \item{lambda.pvalue}{p-value for $\lambda_{min}$ test based on resampling}
   \item{lambda.table.pvalue}{p-value for $\lambda_{min}$ test based on the asymptotic chi-squared distribution (large n1,n2)}
   \item{H}{The Hotelling $T^2$ statistic}
   \item{H.pvalue}{p-value for the Hotelling $T^2$ test based on resampling}
   \item{H.table.pvalue}{p-value for the Hotelling $T^2$ test based on the null F distribution, assuming normality
   and equal covariance matrices}
   \item{J}{The Hotelling $T^2$ statistic}
   \item{J.pvalue}{p-value for the Hotelling $T^2$ test based on resampling}
   \item{J.table.pvalue}{p-value for the Hotelling $T^2$ test based on the null F distribution, assuming normality
   and unequal covariance matrices}
   \item{G}{The Goodall $F$ statistic}
   \item{G.pvalue}{p-value for the Goodall test based on resampling}
   \item{G.table.pvalue}{p-value for the Goodall test based on the null F distribution, assuming normality and
   equal isotropic covariance matrices)}

}

\references{Amaral, G.J.A., Dryden, I.L. and Wood, A.T.A (2006). Pivotal bootstrap methods
for $k$-sample problems in directional statistics and shape analysis. University of Nottingham,
Division of Statistics Technical report.

Dryden, I.L. and Mardia, K.V. (1998) Statistical Shape Analysis,
Wiley, Chichester. Chapter 7.

Goodall, C. R. (1991). Procrustes methods in the statistical analysis
of shape (with discussion). Journal of the Royal Statistical Society,
Series B, 53: 285-339.
}




\author{Ian Dryden}

\seealso{testmeanshapes}

\examples{

#2D example : female and male Gorillas (cf. Dryden and Mardia, 1998)

data(gorf.dat)
data(gorm.dat)

#just select 3 landmarks and the first 10 observations in each group
select<-c(1,2,3)
A<-gorf.dat[select,,1:10]
B<-gorm.dat[select,,1:10]
resampletest(A,B,resamples=100)

}




\keyword{multivariate}% at least one, from doc/KEYWORDS
