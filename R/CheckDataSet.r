# file OrdinalLogisticBiplot/R/CheckDataSet.R
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



CheckDataSet <- function(datanom){

    typeDataFrame = FALSE     #Data are matrix
    nRowInit = nrow(datanom)
    #If there are NA values, with this sentence we will drop the rows with any NA value
    datanom <- na.omit(datanom)

    nRow = nrow(datanom)
    if(nRow < nRowInit){
      print("Be careful. Some rows have been deleted because they present NA values. Check your data.")
    }

    RowNames = NULL
    ColNames = NULL
    LevelNames = list()
    contLevel = 1

    notNumeric = NULL
    if(is.data.frame(datanom)){
      typeDataFrame = TRUE
      RowNames = rownames(datanom)
      ColNames = colnames(datanom)
      for(i in 1:dim(datanom)[2]){
          if(is.factor(datanom[,i])){
              LevelNames[[contLevel]] = levels(datanom[,i])
              contLevel = contLevel + 1
              datanom[,i] = as.numeric(datanom[,i])
          }else if(is.numeric(datanom[,i])){
                    if(!is.integer(datanom[,i])){
                       print(paste("Variable ",i," will be transformed to integer because it is not",sep=""))
                       datanom[,i] = as.integer(datanom[,i])
                    }
                    LevelNames[[contLevel]] = c(1:max(datanom[,i]))
                    contLevel = contLevel + 1
                }else{
                  print(paste("Variable ",i," will be omitted because it is not nominal",sep=""))
                  notNumeric = cbind(notNumeric,i)
                }
      }
      notNumeric = as.vector(notNumeric)
      dataSet = NULL
      for(j in 1:dim(datanom)[2]){
        if(length(which(notNumeric==j))==0){
          dataSet = cbind(dataSet,datanom[,j])
        }
      }
    }else if(is.matrix(datanom)){
              RowNames = dimnames(datanom)[[1]][1:nrow(datanom)]
              ColNames = dimnames(datanom)[[2]][1:ncol(datanom)]
              datanom <- as.data.frame(datanom)
              for(i in 1:dim(datanom)[2]){
                 if(!is.null(levels(datanom[,i]))){
                   LevelNames[[contLevel]] = levels(datanom[,i])
                 }else{
                   LevelNames[[contLevel]] = c(1:max(datanom[,i]))
                 }
                 contLevel = contLevel + 1
                 datanom[,i] = as.numeric(datanom[,i])
                 if(!is.integer(datanom[,i])){
                   print(paste("Variable ",i," will be transformed to integer because it is not",sep=""))
                   datanom[,i] = as.integer(datanom[,i])
                 }
              }
              dataSet = datanom
          }else{
            stop("Data set should be a data frame or a matrix")
          }
    #In case that RowNames or ColNames are NULL we fix the name of the variables and rows
   	if (is.null(RowNames)){
  		RowNames <- rownames(datanom, do.NULL = FALSE, prefix = "I")
  		dimnames(datanom)[[1]] = RowNames
    }

   	if (is.null(ColNames)){
  		ColNames <- colnames(datanom, do.NULL = FALSE, prefix = "V")
  		dimnames(datanom)[[2]] = ColNames
    }


      numVarDef <- ncol(dataSet)
      datanomcats = apply(dataSet[,1:numVarDef], 2, function(x) nlevels(as.factor(x)))
      for(i in 1:numVarDef){
        if(max(dataSet[,i]) > datanomcats[i]){
          #print(paste("Please, it would be desirable that you recode variable ", i ,
          #      " from the data set because its maximum value exceeds the distinct values it presents,
          #      so there are some categories that are not present",sep=""))
          columValuesOrd = sort(unique(dataSet[,i]))
          print(paste("Variable ",dimnames(datanom)[[1]][i]," only take the values:",sep=""))
          print(columValuesOrd)
          ActLevelNames = NULL
          newColdataSet = dataSet[,i]
          for(j in 1:datanomcats[i]){
              newColdataSet = replace(newColdataSet,newColdataSet==columValuesOrd[j],j)
              if(typeDataFrame){
                if(length(LevelNames[[i]]) > 0){
                    ActLevelNames <- c(ActLevelNames,LevelNames[[i]][columValuesOrd[j]])
                }
              }else{
                  ActLevelNames = columValuesOrd
              }
          }
          dataSet[,i] = newColdataSet
          LevelNames[[i]] = ActLevelNames
        }
      }#end for

    dimnames(dataSet)[[1]] = RowNames
    dimnames(dataSet)[[2]] = ColNames

    data.ordinal<-list()
    data.ordinal$datanom = dataSet
    data.ordinal$RowNames = RowNames
    data.ordinal$ColumNames = ColNames
    data.ordinal$LevelNames = LevelNames       

    class(data.ordinal)='data.ordinal'
    return(data.ordinal)

}