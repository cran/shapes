\name{procGPA}
\alias{procGPA}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalised Procrustes analysis}
\description{Generalised Procrustes analysis to register landmark configurations into optimal registration using translation, 
rotation and scaling. Reflection invariance can also be chosen, and registration without scaling is also an option. Also, obtains principal 
components, and some summary statistics. 
}
\usage{
procGPA(x, scale = TRUE, reflect = FALSE, eigen2d = TRUE, 
tol1 = 1e-05, tol2 = tol1, approxtangent = TRUE, proc.output=FALSE, 
distances=TRUE, pcaoutput=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{Input k x m x n real array, (or k x n complex 
matrix for m=2 is OK), where k is the number of points, m 
is the number of dimensions, and n is the sample size. }
  \item{scale}{Logical quantity indicating if scaling is required}
  \item{reflect}{Logical quantity indicating if reflection 
  is required}
  \item{eigen2d}{Logical quantity indicating if complex eigenanalysis 
  should be used to calculate Procrustes mean for 2D case when 
  scale=TRUE, reflect=FALSE} 
   \item{tol1}{Tolerance for optimal rotation for the iterative 
algorithm: tolerance on the mean sum of squares 
(divided by size of mean squared) between successive iterations}
  \item{tol2}{tolerance for rescale/rotation step for the iterative 
algorithm: tolerance on the mean sum of squares 
(divided by size of mean squared) between successive iterations}
\item{approxtangent}{Logical quantity indicating if the approximate tangent 
coordinates (Procrustes residuals) should be given if TRUE, or the partial 
tangent coordinates (see tan below) if FALSE}
\item{proc.output}{Logical quantity indicating if printed output during the
iterations of the Procrustes GPA algorithm should be given} 
\item{distances}{Logical quantity indicating if shape distances and sizes 
should be calculated}
\item{pcaoutput}{Logical quantity indicating if PCA should be carried out} 
}
\value{A list with components
\item{k}{no of landmarks}
\item{m}{no of dimensions (m-D dimension configurations)}
\item{n}{sample size}
\item{mshape}{Procrustes mean shape. Note this is unit size if complex 
eigenanalysis used, but on the scale of the data if iterative GPA is used.} 
\item{tan}{If approxtangent=TRUE this is the mk x n matrix of 
Procrustes approximate tangent coordinates $X_i^P$ - Xbar , where Xbar = mean($X_i^P$), which are also known as the 
Procrustes residuals. If approxtangent=FALSE this is 
the km-m x n matrix of partial Procrustes tangent coordinates
with pole given by the preshape of the Procrustes mean}
\item{rotated}{the k x m x n array of full Procrustes 
rotated data} 
\item{pcar}{the columns are eigenvectors (PCs) of the sample 
covariance Sv of tan} 
\item{pcasd}{the square roots of eigenvalues of Sv using tan (s.d.'s of PCs)}
\item{percent}{the percentage of variability explained by the PCs using
tan}
\item{size}{the centroid sizes of the configurations}
\item{scores}{standardised PC scores (each with unit variance) using 
tan} 
\item{rawscores}{raw PC scores using tan} 
\item{rho}{Kendall's Riemannian distance $\rho$ to the mean shape}
\item{rmsrho}{r.m.s. of rho} 
\item{rmsd1}{r.m.s. of full Procrustes distances to the mean shape $d_F$}
}
\references{Dryden, I.L. and Mardia, K.V. (1998). Statistical 
Shape Analysis, Wiley, Chichester. 
 
Goodall, C.R. (1991). Procrustes methods in the statistical analysis 
of shape (with discussion). Journal of the Royal Statistical Society, 
Series B, 53: 285-339. 

Gower, J.C. (1975). Generalized Procrustes analysis, Psychometrika, 40, 33--50.} 

Kent, J.T. (1994). The complex Bingham distribution and shape analysis,
Journal of the Royal Statistical Society, Series B, 56, 285-299.}

Ten Berge, J.M.F. (1977). Orthogonal Procrustes rotation for two or more 
matrices. Psychometrika, 42, 267-276. 
\author{Ian Dryden, with input from Mohammad Faghihi and Alfred Kume}
\seealso{procOPA,riemdist,shapepca,testmeanshapes}

\examples{

#2D example : female and male Gorillas (cf. Dryden and Mardia, 1998)

data(gorf.dat)
data(gorm.dat)

plotshapes(gorf.dat,gorm.dat)
n1<-dim(gorf.dat)[3]
n2<-dim(gorm.dat)[3]
k<-dim(gorf.dat)[1]
m<-dim(gorf.dat)[2]
gor.dat<-array(0,c(k,2,n1+n2))
gor.dat[,,1:n1]<-gorf.dat
gor.dat[,,(n1+1):(n1+n2)]<-gorm.dat

gor<-procGPA(gor.dat)
shapepca(gor,type="r",mag=3)
shapepca(gor,type="v",mag=3)

gor.gp<-c(rep("f",times=30),rep("m",times=29))
x<-cbind(gor$size,gor$rho,gor$scores[,1:3])
pairs(x,panel=function(x,y) text(x,y,gor.gp),
   label=c("s","rho","score 1","score 2","score 3"))


##########################################################
#3D example

data(macm.dat)
out<-procGPA(macm.dat,scale=FALSE)

par(mfrow=c(2,2))
plot(out$rawscores[,1],out$rawscores[,2],xlab="PC1",ylab="PC2")
title("PC scores")
plot(out$rawscores[,2],out$rawscores[,3],xlab="PC2",ylab="PC3")
plot(out$rawscores[,1],out$rawscores[,3],xlab="PC1",ylab="PC3")
plot(out$size,out$rho,xlab="size",ylab="rho")
title("Size versus shape distance")

}


\keyword{multivariate} 


