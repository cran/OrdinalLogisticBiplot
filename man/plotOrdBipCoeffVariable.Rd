\name{plotOrdBipCoeffVariable}
\alias{plotOrdBipCoeffVariable}
\title{
Plot an ordinal variable to an existing biplot.
}
\description{
This function plots an ordinal variable to an existing biplot or to a new one, using the estimated parameters for the variable.
}
\usage{
plotOrdBipCoeffVariable(x, coeffic, D = 1, planex = 1, planey = 2,
xi = -3.5, xu = 3.5, yi = -3.5, yu = 3.5, margin = 0,
numFactors = 2, CexVar = 0.7, ColorVar = "blue",
PchVar = 0.7, addToPlot = FALSE, showIIC = TRUE,
iicxi = -2.5, iicxu = 2.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Name of the variable the user wants to plot.
}
  \item{coeffic}{
Vector with the estimated coefficients and the thresholds in this order.
}
  \item{D}{
Value that define the graded response model. In case of Mirt estimation it should be 1.702. The default value is 1.
}
  \item{planex}{
Dimension for X axis.
}
  \item{planey}{
Dimension for Y axis.
}
 \item{xi}{
Minimum value on the x-axis.
}
  \item{xu}{
Maximum value on the x-axis.
}
  \item{yi}{
Minimum value on the y-axis.
}
  \item{yu}{
Maximum value on the y-axis.
}
  \item{margin}{
This value stablishes the space between the plotted items and the border of the window.
}
  \item{numFactors}{
Number of dimensions of the solution
}
  \item{CexVar}{
Size of the category points. It can be an array with the cex information for each variable.
}
  \item{ColorVar}{
Color for the variables. It can be an array with the color information for each variable.
}
  \item{PchVar}{
Symbol for the variables. It could be an array with the pch information for each variable.
}
  \item{addToPlot}{
Boolean parameter to decide if the user wants to add the ordinal variable representation to an existing plot.
}
  \item{showIIC}{
Boolean parameter to decide if the user wants to see the item information curves for each variable. Default value is FALSE.
}
  \item{iicxi}{
Lower limit for the X-axis when plotting item information curves.
}
  \item{iicxu}{
Upper limit for the X-axis when plotting item information curves.
}
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}
\examples{
data(LevelSatPhd)
olbo = OrdinalLogisticBiplot(LevelSatPhd,sFormula=NULL,
      numFactors=2,coordinates="EM",penalization=0.2)
nameVariable="Salary"
coeffic = olbo$VariableModels[,1]$coef
plotOrdBipCoeffVariable(nameVariable,coeffic,D = 1,
      planex = 1,planey = 2,xi=-2.5,xu=2.5,yi=-2.5,
      yu=2.5,addToPlot=FALSE,iicxi=-2.5,iicxu=2.5)
}

\keyword{plot}
