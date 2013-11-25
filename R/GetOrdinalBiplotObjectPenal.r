# file OrdinalLogisticBiplot/R/GetOrdinalBiplotObjectPenal.R
# copyright (C) 2012-2013 J.C. Hernandez and J.L. Vicente-Villardon
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#


GetOrdinalBiplotObjectPenal <- function(ColumNames,olb,planex=1,planey=2){

    nColsdataF <- nrow(olb$Ncats)
    numFactors = ncol(olb$RowCoordinates)

    matBiplot=0
    for(nvarColum in 1:nColsdataF){
      numcat = olb$Ncats[nvarColum,1]
      coeffic = ExtractCoefficientsPenal(olb,nvarColum,numFactors,numcat)
      D = 1
      nameVariable = ColumNames[nvarColum]
      ordBipVar = CalculateOrdinalBiplotGeneral(nameVariable,numcat,coeffic,planex,planey,numFactors,D)

      matBiplot=cbind(matBiplot,ordBipVar)
    }
    matBiplot=matBiplot[,2:(nColsdataF+1)]

    result = list()
    result$models = olb
  	result$matBiplot = matBiplot
  	class(result) <- "CategOrdBiplotPen"
  	return(result)
}
