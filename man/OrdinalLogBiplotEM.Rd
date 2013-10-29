\name{OrdinalLogBiplotEM}
\alias{OrdinalLogBiplotEM}
\title{
Alternated EM algorithm for Ordinal Logistic Biplots
}
\description{
This function computes, with an alternated algorithm, the row and column parameters of an Ordinal Logistic Biplot for ordered polytomous data. The row coordinates (E-step) are computed using multidimensional Gauss-Hermite quadratures and Expected \emph{a posteriori} (EAP) scores and parameters for each variable or items (M-step)using Ridge Ordinal Logistic Regression to solve the separation problem present when the points for different categories of a variable are completely separated on the representation plane and the usual fitting methods do not converge. The separation problem is present in almost avery data set for which the goodness of fit is high. 
}
\usage{
OrdinalLogBiplotEM(x,dim = 2, nnodos = 15, tol = 0.001, maxiter = 100,
                  penalization = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Matrix with the ordinal data. The matrix must be in numerical form. 
}
  \item{dim}{
Dimension of the solution
}
  \item{nnodos}{
Number of nodes for the multidimensional Gauss-Hermite quadrature
}
 \item{tol}{
Value to stop the process of iterations.
}
  \item{maxiter}{
Maximum number of iterations in the process of solving the regression coefficients.
}
  \item{penalization}{
Penalization used in the diagonal matrix to avoid singularities.
}
}

\value{
An object of class \code{"ordinal.logistic.biplot.EM"}.This has components:

  \item{RowCoordinates}{Coordinates for the rows or the individuals}
  \item{ColumnParameters}{List with information about the Ordinal Logistic Models calculated for each variable including: estimated parameters with thresholds,percents of correct classifications,and pseudo-Rsquared}
  \item{loadings}{factor loadings}
  \item{LogLikelihood}{Logarithm of the likelihood}
  \item{r2}{R squared coefficient}
  \item{Ncats}{Number of the categories of each variable}

}
\references{
  Bock,R. & Aitkin,M. (1981),\emph{Marginal maximum likelihood estimation of item parameters: Aplication of an EM algorithm}, Phychometrika 46(4), 443-459.

}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}
\seealso{
\code{\link{pordlogist}}
}
\examples{
    data(LevelSatPhd)
    dataSet = CheckDataSet(LevelSatPhd)
    datanom = dataSet$datanom
    olb = OrdinalLogBiplotEM(datanom,dim = 2, nnodos = 10,
          tol = 0.001, maxiter = 100, penalization = 0.2)
    olb
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{EM}
\keyword{algorithm}