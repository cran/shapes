\name{procGPA}
\alias{procGPA}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generalised Procrustes analysis}
\description{Generalised Procrustes analysis to register geometrical 
  objects into optimal registration using translation, 
  rotation and scaling. Reflection invariance can also be chosen, and
registration without scaling is also an option. Also, obtains principal 
  components, and some summary statistics. 
}
\usage{
procGPA(x, scale = TRUE, reflect = FALSE, tol1 = 1e-05, tol2 = 1e-05)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{Input k x m x n real array, (or k x n complex 
matrix for m=2 is OK), where k is the number of points, m 
is the number of dimensions, and n is the sample size. }
  \item{scale}{Logical quantity indicating if scaling is required}
  \item{reflect}{Logical quantity indicating if reflection 
  is required}
  \item{tol1}{Tolerance for optimal rotation for the iterative 
algorithm: tolerance on the mean sum of squares between successive iterations 
(depends on scale of objects)}
  \item{tol2}{tolerance for rescale/rotation step for the iterative 
algorithm: tolerance on the Riemannian shape distance 
between successive mean shapes}
}
\value{A list with components
\item{k}{no of landmarks}
\item{m}{no of dimensions (m-D dimension configurations)}
\item{n}{sample size}
\item{mshape}{full Procrustes mean shape} 
\item{tan}{the real mk x n matrix of full Procrustes tangent coordinates
 $X_i^P$ - Xbar , where Xbar = mean($X_i^P$)}
\item{rotated}{the k x m x n array of real full Procrustes 
rotated data} 
\item{pcar}{the columns are eigenvectors (PCs) of the sample 
covariance Sv of tan}
\item{pcasd}{the square roots of eigenvalues of Sv (s.d.'s of PCs)}
\item{percent}{the percentage of variability explained by the PCs}
\item{size}{the centroid sizes of the configurations}
\item{scores}{standardised PC scores (each with unit variance)} 
\item{rawscores}{raw PC scores} 
\item{rho}{Kendall's Procrustean (Riemannian) distance rho 
to the mean shape}
\item{rmsrho}{r.m.s. of rho} 
\item{rmsd1}{r.m.s. of full Procrustes distances 
to the mean shape d1}
}

\references{Dryden, I.L. and Mardia, K.V. (1998). 
Statistical Shape Analysis, Wiley, Chichester. 
 
Goodall, C.R. (1991). Procrustes methods in the statistical analysis 
of shape (with discussion). Journal of the Royal Statistical Society, 
Series B, 53: 285-339. 

Gower, J.C. (1975). Generalized Procrustes analysis, Psychometrika, 40, 33--50.} 

Kent, J.T. (1994). The complex Bingham distribution and shape analysis,
Journal of the Royal Statistical Society, Series B, 56, 285-299.}

\author{Ian Dryden}
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

data(dna.dat)
subsetobs<-c(1:10)
subsetk<-c(1,3,5,7,9,11,12,14,16,18,20,22)
dna<-procGPA(dna.dat[subsetk,,subsetobs],tol1=0.1)

par(mfrow=c(2,2))
plot(dna$rawscores[,1],dna$rawscores[,2],type="l",xlab="PC1",ylab="PC2")
points(dna$rawscores[1,1],dna$rawscores[1,2])
title("PC scores - joined up through time")
plot(dna$rawscores[,2],dna$rawscores[,3],type="l",xlab="PC2",ylab="PC3")
points(dna$rawscores[1,2],dna$rawscores[1,3])
plot(dna$rawscores[,1],dna$rawscores[,3],type="l",xlab="PC1",ylab="PC3")
points(dna$rawscores[1,1],dna$rawscores[1,3])
plot(dna$size,dna$rho,xlab="size",ylab="rho",type="l")
points(dna$size[1],dna$rho[1])
title("Size versus shape distance")

}


\keyword{multivariate} 


