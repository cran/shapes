\name{PED}
\alias{PED}
\title{Penalized Euclidean Distance regression}
\description{Penalized Euclidean Distance Regression is a sparse regression method where 
the objective function is a Euclidean distance (square root of least squares objective function)
plus a penalty involving the geometric mean of the l_1 and l_2 norms of the parameters:

              || Y - X beta || + lambda sqrt( | beta |_1 || beta || )

where ||.|| is the Euclidean (l_2) norm, and |beta|_1 is the l_1 norm.  
The objective function is minimized using the limited memory BFGS algorithm, and then 
a threshold is used to make the estimator sparse. 
}
\usage{
PED(X, Y, method="AIC")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{X}{Input n x p matrix, where n is the number of observations, p 
is the number of covariates. }
  \item{Y}{Input n vector of responses}
  \item{method}{Method for selection of penalty parameter lambda and threshold parameter delta:
"AIC" - Akaike's information criterion
"khat" - maximizing the sparsity criterion khat
"CV" - cross-validation  
} 
}
\value{A list with components
\item{coef}{fitted regression coefficients. First parameter is the intercept, followed by 
the p parameters corresponding to the covariates} 
\item{betahat}{fitted regression coefficients corresponding to the p covariates}
\item{yhat}{Fitted Y values} 
\item{resid}{Residuals: vector Y-yhat}
\item{lambda}{Penalty parameter lambda used}
\item{delta}{Threshold delta used}
}
\references{
Vasiliu, D., Dey T. and Dryden I.L. (2017).
Penalized Euclidean Distance Regression. arXiv 1405:4578.
}
\author{Ian Dryden}
\seealso{lm}

\examples{

library(MASS)
data(Boston)
X<-as.matrix(Boston[,1:13])
Y<-as.double(Boston[,14])
out<-PED(X,Y)
print(out$coef)
}


\keyword{multivariate} 


