source("Scripts/R/lib/UtilityFunctions.r")

dataSetNames <- function() {
  return(c("whole", "partial", "apical", "basal", "sholl", "parameterSpace"))
}

removeNeuronWithName <- function(dataFrames, neuronName) {
  result = list()
  for (eachPosition in seq_along(dataFrames)) {
    eachName <- names(dataFrames)[eachPosition]
    eachFrame <- dataFrames[[eachName]]
    newFrame <- eachFrame[substr(eachFrame$cellName, 1,
                                 nchar(neuronName)) != neuronName,]
    result[[eachName]] <- newFrame
  }
  return(result)
}

getApicalNeuronsFromFrame <- function(dataFrame) {
  return(getNeuronsFromFrameWithNamesThatEndInString(dataFrame, "-apical"))
}

getBasalNeuronsFromFrame <- function(dataFrame) {
  return(getNeuronsFromFrameWithNamesThatEndInString(dataFrame, "-basal"))
}

getNeuronsFromFrameWithNamesThatEndInString <- function(
    dataFrame, searchString) {
  matchingRows <- sapply(dataFrame$cellName,
                         function(x) { endsWith(as.character(x),
                                                searchString) })
  result <- dataFrame[matchingRows,]
}

# Given a data frame with an "ageCategory" column, returns a list of two
# dataFrames, one young, one old, with the names "young" and "old" respectively.
splitYoungAndOldForOneFrame <- function(dataFrame) {
  result <- list()
  result$young <- dataFrame[dataFrame$ageCategory == "YOUNG",]
  result$old <- dataFrame[dataFrame$ageCategory == "OLD",]
  return(result)
}

splitYoungAndOld <- function(dataFrames) {
  result = list()
  for (eachName in dataSetNames()) {
    eachFrame <- dataFrames[[eachName]]
    result[[eachName]] <- eachFrame
    youngName <- paste0(eachName, "Young")
    oldName <- paste0(eachName, "Old")

    result[[youngName]] <- eachFrame[eachFrame["ageCategory"] == "YOUNG", ]
    result[[oldName]] <- eachFrame[eachFrame["ageCategory"] == "OLD", ]
  }
  return(result)
}
