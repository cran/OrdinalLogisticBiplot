\name{plot.ordinal.logistic.biplot}
\alias{plot.ordinal.logistic.biplot}
\title{
Graphical representation of an Ordinal Logistic Biplot.
}
\description{
Plotting an Ordinal Logistic Biplot. There are parameters related to the way in which the biplot is plotted. All the posible parameters have default values.
}
\usage{
# S3 method for class 'ordinal.logistic.biplot'
\method{plot}{ordinal.logistic.biplot}(x, planex = 1, planey = 2,
AtLeastR2 = 0.01, xlimi = -1.5, xlimu = 1.5, ylimi = -1.5,
ylimu = 1.5, margin = 0, ShowAxis = TRUE, PlotVars = TRUE,
PlotInd = TRUE, LabelVar = TRUE, LabelInd = TRUE, CexInd = NULL,
CexVar = NULL, ColorInd = NULL, ColorVar = NULL, PchInd = NULL,
PchVar = NULL, showIIC = FALSE, iicxi = -1.5, iicxu = 1.5,
legendPlot = FALSE,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of the class ordinal.logistic.biplot.
}
  \item{planex}{
Dimension for X axis.
}
  \item{planey}{
Dimension for Y axis.
}
  \item{AtLeastR2}{
It establishes the cutting value to plot a variable attending to its Nagelkerke pseudo R squared value. A variable is plotted if its pseudo R squared is higher than this value.
}
  \item{xlimi}{
Minimum value on the x-axis.
}
  \item{xlimu}{
Maximum value on the x-axis.
}
  \item{ylimi}{
Minimum value on the y-axis.
}
  \item{ylimu}{
Maximum value on the y-axis.
}
  \item{margin}{
This value stablishes the space between the plotted items and the border of the window.
}
 \item{ShowAxis}{
  Should the axis be shown?
}
  \item{PlotVars}{
  Should the variables (items) be plotted?
}
  \item{PlotInd}{
Should the individuals be plotted?
}
  \item{LabelVar}{
Should the variable labels be shown?
}
  \item{LabelInd}{
Should the individual labels be shown?
}
  \item{CexInd}{
Size of the individual points. It can be an array with the cex information for each row.
}
  \item{CexVar}{
Size of the category points. It can be an array with the cex information for each variable.
}
  \item{ColorInd}{
Color of the individual points. It can be an array with the color information for each row.
}
  \item{ColorVar}{
Color for the variables. It can be an array with the color information for each variable.
}
  \item{PchInd}{
Symbol for the individuals. It can be an array with the pch information for each row.
}
  \item{PchVar}{
Symbol for the variables. It could be an array with the pch information for each variable.
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
  \item{legendPlot}{
Boolean parameter to show the legend of the plot. Default value is FALSE.
}
  \item{\dots}{
Additional parameters to plot.
}
}
\details{
The function without parameters plots the ordinal.logistic.biplot object with labels in the original data and default values for colors, symbols and sizes for points and lines. Other values of colors, symbols and sizes can be supplied. A single value applies to all the points but an array with different values can be used to improve the undestanding of the plot.-
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}
\seealso{
\code{\link{OrdinalLogisticBiplot}}
}
\examples{
data(LevelSatPhd)
olbo = OrdinalLogisticBiplot(LevelSatPhd,penalization=0.2)
plot(olbo,PlotInd=TRUE,xlimi=-1.5,xlimu=1.5,ylimi=-1.5,ylimu=1.5,
  margin = 0.2, ColorVar = c("red","green","black","blue","yellow"),
  CexVar = c(0.7),showIIC=FALSE)
}

\keyword{plot}