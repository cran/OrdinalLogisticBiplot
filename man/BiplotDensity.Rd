\name{BiplotDensity}
\alias{BiplotDensity}
\title{
Density plot of a data set with overlaid contours. 
}
\description{
This function draws for a set of points a density contour lines plot.
The densities can be calculated for the whole
set of points or for the groups defines by a nominal variable.
}
\usage{
BiplotDensity(X, y = NULL, nlevels = max(y), grouplabels = 1:nlevels,
  ncontours = 6, groupcols = 1:nlevels, img = TRUE, separate = FALSE,
  ncolors = 20, ColorType = 4, xliml = -1, xlimu = 1, yliml = -1,
  ylimu = 1, plotInd = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A matrix with the  coordinates for the plane in which the the contour lines will be plotted.
}
  \item{y}{
Categorical variable defining clusters. If NULL, the density is calculated for the whole set of points.
}
  \item{nlevels}{
Number of clusters.
}
  \item{grouplabels}{
Set of labels for the centers of each cluster. It should be a vector with "nlevels" components.
}
  \item{ncontours}{
Number of contours that will be used in the representation.
}
  \item{groupcols}{
Vector whith a set of colors for the clusters.
}
  \item{img}{
Should the densty be plotted (with different colors) together with the contour lines?. Default value is TRUE. 
}
  \item{separate}{
Should the density for each cluster be represented on a different picture?. Default value is FALSE. 
}
  \item{ncolors}{
Number of colors for the densities.
}
  \item{ColorType}{
Type of color schema for the density image. It should be a number between 1 and 5.
}
 \item{xliml}{
Minimum value on the x-axis.
}
  \item{xlimu}{
Maximum value on the x-axis.
}
  \item{yliml}{
Minimum value on the y-axis.
}
  \item{ylimu}{
Maximum value on the y-axis.
}
  \item{plotInd}{
Should the individuals be plotted?
}
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}
\examples{
data(LevelSatPhd)
olbo = OrdinalLogisticBiplot(LevelSatPhd)
x = olbo$RowCoords[, 1]
y = olbo$RowCoords[, 2]
plot(x,y, cex = 0, xlim=c(-1,1),ylim=c(-1,1))
X = olbo$RowCoords
y = as.matrix(as.numeric(LevelSatPhd[,4]))
gcols = c("midnightblue","black","red","gray87")
BiplotDensity(X,y,groupcols = gcols)
}
\keyword{Density}