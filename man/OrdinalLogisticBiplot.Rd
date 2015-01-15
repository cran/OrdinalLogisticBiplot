\name{OrdinalLogisticBiplot}
\alias{OrdinalLogisticBiplot}
\title{Ordinal Logistic Biplot for ordered polytomous data}
\description{
Function that calculates the parameters of the Ordinal Logistic Biplot. 
}
\usage{
OrdinalLogisticBiplot(datanom,sFormula=NULL,numFactors=2,
method="EM",rotation="varimax",metfsco="EAP",
nnodos = 10, tol = 1e-04, maxiter = 100,
penalization = 0.1,cte=TRUE, show=FALSE,ItemCurves = FALSE,initial=1,alfa=1)
}
\arguments{
  \item{datanom}{
The data set; it can be a \emph{matrix} with integers or a \emph{data frame} with factors. All variables have to be ordinal. 
}
  \item{sFormula}{
This parameter follows the unifying interface for selecting variables from a data frame for a plot,
test or model. The most common formula is of type y ~ x1+x2+x3. It has a default value of NULL if it is not specified.
}
  \item{numFactors}{
Number of dimensions of the solution. It should be lower than the number of variables. It has a default value of 2.
}
  \item{method}{
This parameter can be: "EM" or "MIRT". Method to compute the row coordinates.
}
  \item{rotation}{
Rotation method to used with "MIRT" option in  "coordinates". No effect fror other options.
}
  \item{metfsco}{
Calculation method for the fscores with  "MIRT" option in  "coordinates". No effect for other options.
}
  \item{nnodos}{
Number of nodes for gauss quadrature in the EM algorithm.
}
  \item{tol}{
Tolerance for the EM algorithm.
}
  \item{maxiter}{
Maximum number of iterations in the EM algorithm.
}
  \item{penalization}{
Penalization for the ridge regression for each variable.
}
  \item{cte}{
Include constant in the logistic regression model. Default is TRUE.
}
  \item{show}{
Show intermediate computations. Default is FALSE.
}
\item{ItemCurves}{
Show item information curves. Default is FALSE.
}
\item{initial}{
Method used to choose the initial ability in the EM algorithm. Default value is 1.
}
\item{alfa}{
Optional parameter to calculate row and column coordinates in Simple correspondence analysis if the initial parameter is equal to 1. Default value is 1.
}
}
\details{
The general algorithm used is essentially an alternated procedure in which parameters for rows and columns are computed in alternated steps repeated until convergence. Parameters for the rows are calculated by expectation (E-step) and parameters for the columns are computed by maximization (M-step), i. e., by Ordinal Logistic Regression. 

There are several options for the computation:

1.- Using the package \pkg{mirt} to obtain the row scores, i. e. using a solution obtained from a latent trait model. The column (item) parameters should be directly used by our biplot procedure but, because of the characteristics of the package that performs a default rotation after parameter estimation, we have to reestimate the item parametes to be coherent to the scores.

2.- Using our implementation of the EM algorithm alternating expected a porteriori scores and Ridge Ordinal Logistic Regression for each variable. We use here a Cumulative link model ,that is, a logistic regression model for cumulative logits.

Equations defining the set of probability response surfaces for the cumulative probabilities are sigmoidal as in the binary case (Vicente-Villardon et al.2006) and then they share its geometry. All categories have a different constant but the same slopes, that means that the prediction direction is common to all categories and just the prediction markers are different. The representation subspace can be divided into prediction regions, for each category, delimited by parallel straight lines.
}
\value{
An object of class \code{"ordinal.logistic.biplot"}.  This has some components:

  \item{dataSet }{Data set of study with all the information about the name of
                  the levels and names of the variables and individuals}
  \item{RowCoords }{Coordinates for the rows in the reduced space}
  \item{NCats }{Number of categories of each variable from the data set}  
  \item{estimObject }{Object with all the estimated information using EM alternated algorithm or MIRT procedure}
  \item{Fitting  }{matrix with the percentage of correct clasifications and pseudo R squared valued for each variable}
  \item{coefs  }{matrix with the estimated coefficients}
  \item{thresholds  }{matrix with the estimated intercept limits}
  \item{NumFactors  }{Number of dimensions selected for the study}  
  \item{Coordinates }{Type of coordinates to calculate the row positions}
  \item{Rotation }{Type of rotation if we have chosen mirt coordinates}
  \item{Methodfscores }{Method of calculation of the fscores in mirt process}
  \item{NumNodos }{Number of nodes for the gauss quadrature in EM algorithm}
  \item{tol }{Cut point to stop the EM-algorithm}
  \item{maxiter }{Maximum number of iterations in the EM-algorithm}
  \item{penalization }{Value for the correction of the ridge regression}
  \item{cte }{Boolean value to choose if the model for each variable will have
              independent term}
  \item{show }{Boolean value to indicate if we want to see the results of our analysis}
  \item{ItemCurves }{Boolean value to specify if item information curves will be plotted}
  \item{LogLik }{Logarithm of the likelihood}  
  \item{FactorLoadingsComm }{Factor loadings and communalities}
}
\references{
Vicente-Villardon, J., Galindo, M.P & Blazquez-Zaballos, A. (2006), \emph{Logistic biplots},Multiple Correspondence Analysis and related methods pp. 491--509.

Demey, J., Vicente-Villardon, J. L., Galindo, M.P.  & Zambrano, A. (2008) \emph{Identifying Molecular Markers Associated With Classification Of Genotypes Using External Logistic Biplots}. Bioinformatics, 24(24), 2832-2838.

Baker, F.B. (1992): Item Response Theory. Parameter Estimation Techniques. Marcel Dekker. New York.

Gabriel, K. (1971), \emph{The biplot graphic display of matrices with application to principal component analysis.}, Biometrika 58(3), 453--467.

Gabriel, K. R. (1998), \emph{Generalised bilinear regression}, Biometrika 85(3),
689--700.

Gabriel, K. R. & Zamir, S. (1979), \emph{Lower rank approximation of matrices by least squares with any choice of weights}, Technometrics 21(4), 489--498.

Gower, J. & Hand, D. (1996), \emph{Biplots, Monographs on statistics and applied probability}. 54. London: Chapman and Hall., 277 pp.

Chalmers,R,P (2012). \emph{mirt: A Multidimensional Item Response Theory Package for the R Environment}. Journal of Statistical Software, 48(6), 1-29. URL http://www.jstatsoft.org/v48/i06/.
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}

\seealso{
\code{\link{OrdinalLogBiplotEM}}
}
\examples{

data(LevelSatPhd)
olbo = OrdinalLogisticBiplot(LevelSatPhd)
summary(olbo)


}
\keyword{biplot}
\keyword{logistic}
\keyword{models}
