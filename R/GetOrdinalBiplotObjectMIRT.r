# file OrdinalLogisticBiplot/R/GetOrdinalBiplotObjectMIRT.R
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

GetOrdinalBiplotObjectMIRT <- function(modelMirt,planex=1,planey=2){

    nColsdataF <- ncol(modelMirt$dataFactor)
    varStudy = modelMirt$dataFactor

    matBiplot=0
    for(nvarColum in 1:nColsdataF){
      numcat = max(varStudy[,nvarColum])
      coeffic=modelMirt$coefMirt[[nvarColum]]
      numFactors = length(coeffic)- (numcat - 1)
      coeffic = c(modelMirt$coefMirt[[nvarColum]][1:numFactors],(-1)*modelMirt$coefMirt[[nvarColum]][(numFactors+1):length(coeffic)])
      D=1.702
      nameVariable = dimnames(varStudy)[[2]][nvarColum]
      ordBip = CalculateOrdinalBiplotGeneral(nameVariable,numcat,coeffic,planex,planey,modelMirt$numFactors,D)
      matBiplot=cbind(matBiplot,ordBip)
    }
    matBiplot=matBiplot[,2:(nColsdataF+1)]

    result = list()
    result$coefMirt = modelMirt$coefMirt
    result$sepCoefMirt = modelMirt$sepCoefMirt
    result$numFactors = modelMirt$numFactors
    result$rotation = modelMirt$rotation
    result$metfsco = modelMirt$metfsco
    result$planex = planex
    result$planey = planey
  	result$scores = modelMirt$estimRows
    result$matBiplot = matBiplot
    result$summaryMirt = modelMirt$summ
  	class(result) <- "CategOrdBiplot"
  	return(result)
}