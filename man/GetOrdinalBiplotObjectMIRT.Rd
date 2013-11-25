\name{GetOrdinalBiplotObjectMIRT}
\alias{GetOrdinalBiplotObjectMIRT}
\title{
Calculation of Mirt ordinal biplot object.
}
\description{
This function uses the estimation from the MIRT method to calculate all the information needed for representing the ordinal biplot: the biplot axis for each variable with the points separating the visible categories.
}
\usage{
GetOrdinalBiplotObjectMIRT(modelMirt, planex = 1, planey = 2)
}
\arguments{
  \item{modelMirt}{
Object of class "EstimationRowsMIRT" with all the estimated parameters with Mirt method.
}
  \item{planex}{
Dimension for X axis.
}
  \item{planey}{
Dimension for Y axis.
}
}
\details{
This function should only be used if the user selects "MIRT" coordinates as the estimation method in OrdinalLogisticBiplot function.
}
\value{
An object of class \code{"CategOrdBiplot"}.  This has components:
 \item{coefMirt }{Estimated coefficients from Mirt method}
 \item{sepCoefMirt}{Separated coefficients for the independent variables and for the thresholds}
 \item{numFactors }{Number of dimensions of the solution}
 \item{rotation }{Rotation method used with "MIRT" option}
 \item{metfsco }{Calculation method for the fscores with  "MIRT" option}
 \item{planex }{Dimension for X axis.}
 \item{planey}{Dimension for Y axis.}
 \item{scores}{Coordinates for the rows in the reduced space}
 \item{matBiplot }{This item keeps a matrix with the same structure than in GetOrdinalBiplotObjectPenal function}
 \item{summaryMirt }{Summary parameters for the mirt procedure}
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}
\seealso{
\code{\link{OrdinalLogisticBiplot}},\code{\link{GetOrdinalBiplotObjectPenal}}
,\code{\link{EstimationRowsMIRT}}
}
\examples{
data(LevelSatPhd)
olboM = OrdinalLogisticBiplot(LevelSatPhd,coordinates="MIRT",metfsco="EAP",rotation="varimax")
catOrdBiplot = GetOrdinalBiplotObjectMIRT(olboM$estimObject)
}
\keyword{MirtOrdinalObject}