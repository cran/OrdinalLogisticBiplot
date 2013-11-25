\name{EstimationRowsMIRT}
\alias{EstimationRowsMIRT}
\title{
Mirt parameter estimation for a data set.
}
\description{
This function uses mirt and fscores functions of "mirt" package to estimate the coordinates for each row and the parameters for each variable.
}
\usage{
EstimationRowsMIRT(dataFile,numFactors=2,metfsco="EAP",rotation="varimax",maxiter=100)
}
\arguments{
  \item{dataFile}{
Data set with ordinal variables.
}
  \item{numFactors}{
Number of dimensions of the solution. It should be lower than the number of variables. It has a default value of 2.
}
  \item{metfsco}{
Calculation method for the fscores.
}
  \item{rotation}{
Rotation method.
}
  \item{maxiter}{
Maximum number of iterations executed by the algorithm.
}
}
\value{
An object of class \code{"EstimationRowsMIRT"}.  This has some components:
  \item{estimRows }{Coordinates for the rows or the individuals}
  \item{dataFactor }{Data set of ordinal variables}
  \item{rotation }{Rotation method}  
  \item{metfsco }{Calculation method for the fscores}
  \item{numFactors}{Number of dimensions of the solution}
  \item{coefMirt}{List with the estimated coefficients for all the variables}
  \item{sepCoefMirt}{Object with two matrices with the estimated coefficients for the independent variables and for the thresholds}
  \item{summ}{Object with the information provided by summary function for mirt one}  
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}
\examples{
data(LevelSatPhd)
erm = EstimationRowsMIRT(LevelSatPhd)
      
}
\keyword{mirt}
