\name{PlotClusterVariable}
\alias{PlotClusterVariable}
\title{
Graphical representation of clustering variables.
}
\description{
This function uses a clustering variable to represent the centers or convex hulls regions of the points that share each cluster. This variable could be the result of the application of one of the main clustering methods to the data set.
}
\usage{
PlotClusterVariable(A, Groups = ones(c(nrow(A), 1)),
        colors = NULL, chulls = TRUE, centers = TRUE, ConfidentLevel = 0.95)
}
\arguments{
  \item{A}{
A matrix with the items coordinates for the plane of study. It should have only two columns.
}
  \item{Groups}{
Clustering variable: it keeps the cluster for each observation.
}
  \item{colors}{
It is a vector used to specify the colors of the centers and convex hulls regions for each cluster. 
}
  \item{chulls}{
Should convex hulls regions be plotted?
}
  \item{centers}{
Should centers of each cluster be plotted?
}
  \item{ConfidentLevel}{
Numerical value between 0 and 1. If it's value is 0.95, five percent of the points with higher distances to the center of each cluster will not be used for the convex hulls and center calculations that will be plotted.
}
}
\author{
Julio Cesar Hernandez Sanchez, Jose Luis Vicente-Villardon

Maintainer: Julio Cesar Hernandez Sanchez <juliocesar_avila@usal.es>
}
\examples{
data(LevelSatPhd)
olbo = OrdinalLogisticBiplot(LevelSatPhd)
x = olbo$ItemsCoords[, 1]
y = olbo$ItemsCoords[, 2]
plot(x,y, cex = 0.8, pch=17, xlim=c(-2,2),ylim=c(-2,2))
GroupsF = as.factor(LevelSatPhd[,4])
PlotClusterVariable(olbo$ItemsCoords, Groups = GroupsF,
    colors = c(1,2,3,4),chulls = TRUE,centers = TRUE,ConfidentLevel=NULL)

}
\keyword{cluster}
