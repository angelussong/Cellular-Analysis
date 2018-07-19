# These are the columns in the apical and basal data sets that can be translated
# into whole-neuron measurements simply by adding them.
SummableApicalBasalColumns <- function() {
  return(c("out.0Hz.nospines", "out.100Hz.nospines", "out.200Hz.nospines",
           "out.300Hz.nospines", "out.400Hz.nospines", "out.500Hz.nospines",
           "in.0Hz.nospines", "in.100Hz.nospines", "in.200Hz.nospines",
           "in.300Hz.nospines", "in.400Hz.nospines", "in.500Hz.nospines",
           "out.0Hz.spines", "out.100Hz.spines", "out.200Hz.spines",
           "out.300Hz.spines", "out.400Hz.spines", "out.500Hz.spines",
           "in.0Hz.spines", "in.100Hz.spines", "in.200Hz.spines",
           "in.300Hz.spines", "in.400Hz.spines", "in.500Hz.spines", "volume",
           "surfaceArea", "totalLength", "numberOfSections", "spineCount",
           "spineVolume", "spineSurfaceArea"))
}

MorphologicalColumns <- function() {
  return(c("somaSurfaceArea", "simRn", "firingRate.130", "firingRate.180",
           "firingRate.230", "firingRate.280", "firingRate.330",
           "firingRate.380", "surfaceArea", "volume", "totalLength",
           "numberOfSections", "spineCount", "spineVolume", "spineSurfaceArea",
           "out.0Hz.nospines", "out.100Hz.nospines", "out.200Hz.nospines",
           "out.300Hz.nospines", "out.400Hz.nospines", "out.500Hz.nospines",
           "in.0Hz.nospines", "in.100Hz.nospines", "in.200Hz.nospines",
           "in.300Hz.nospines", "in.400Hz.nospines", "in.500Hz.nospines",
           "out.0Hz.spines", "out.100Hz.spines", "out.200Hz.spines",
           "out.300Hz.spines", "out.400Hz.spines", "out.500Hz.spines",
           "in.0Hz.spines", "in.100Hz.spines", "in.200Hz.spines",
           "in.300Hz.spines", "in.400Hz.spines", "in.500Hz.spines",
           "mBPAP.nospines", "mBPAP.spines"))
}

EphysColumns <- function() {
  return(c("epspsHz", "vrestMv", "periodMs", "inputResistance",
           "thresholdMv", "amplitudeMv", "durationMs", "riseTimeMs",
           "decayTimeMs", "spikeCount30pA", "spikeCount80pA",
           "spikeCount130pA", "spikeCount180pA", "spikeCount230pA",
           "spikeCount280pA", "spikeCount330pA"))
}

LoadEmpiricalData <- function() {
  physiologyAndMorphology <- read.csv("EmpiricalData/PhysiologyAndMorphology.csv")
  somaMeasurements <- read.csv("EmpiricalData/SomaMeasurements.csv")
  stacksSummary <- read.csv("EmpiricalData/StacksSummary.csv")
  somaDepthRemeasurement <- read.csv("EmpiricalData/somaDepthRemeasurement.csv")

  empiricalData <- merge(physiologyAndMorphology, somaMeasurements,
                         c("cellName"), all=TRUE)
  empiricalData <- merge(empiricalData, stacksSummary, c("cellName"), all=TRUE)
  empiricalData <- merge(empiricalData, somaDepthRemeasurement, c("cellName"),
                         all=TRUE)

  return(empiricalData)
}

LoadWholeNeuronNumericalResults <- function() {
  firingRates <- read.csv("NumericalResults/firingRates.csv")
  inputResistance <- read.csv("NumericalResults/inputResistance.csv")
  mbpapWholeNospines <- read.csv("NumericalResults/mbpap-whole-nospines.csv")
  mbpapWholeSpiny <- read.csv("NumericalResults/mbpap-whole-spiny.csv")

  firingRates$stim <- factor(firingRates$stim)
  firingRatesWidened <- reshape(firingRates, v.names="firingRate",
                                idvar="cellName", timevar="stim",
                                direction="wide")

  result <- merge(inputResistance, firingRatesWidened,
                  c("cellName", "parameterSet"), all=TRUE)
  result <- merge(result, mbpapWholeNospines, c("cellName", "parameterSet"),
                  all=TRUE)
  result <- merge(result, mbpapWholeSpiny, c("cellName", "parameterSet"),
                  all=TRUE)

  return(result)
}

LoadWholeNeuronData <- function(apicalData, basalData) {
  empirical <- LoadEmpiricalData()
  computedWhole <- LoadWholeNeuronNumericalResults()
  result <- merge(empirical, computedWhole, c("cellName"),
                  all=TRUE)

  summableColumns <- SummableApicalBasalColumns()
  mergeableColumns <- c("cellName", "parameterSet")
  mergeableColumns[3:(length(summableColumns) + 2)] <- summableColumns

  summableApicalData <- apicalData[, mergeableColumns]
  summableBasalData <- basalData[, mergeableColumns]
  summableApicalData$cellName <- sub("-apical", "-all", apicalData$cellName)
  summableBasalData$cellName <- sub("-basal", "-all", basalData$cellName)

  # Note that the usual "all=TRUE" isn't used below, so neurons that lack an
  # apical or basal tree are dropped.
  suppressWarnings(combinedTreeData <- merge(summableApicalData, summableBasalData,
      c("cellName", "parameterSet"), suffixes=c("", "")))

  collapsedTreeData <- as.data.frame(
      sapply(
          summableColumns,
          function(x) rowSums(
              combinedTreeData[, colnames(combinedTreeData) == x, drop=FALSE])))
  collapsedTreeData <- cbind(combinedTreeData[, c("cellName", "parameterSet")],
                             collapsedTreeData)

  result <- merge(result, collapsedTreeData, c("cellName", "parameterSet"),
                  all=TRUE)

  return(result)
}

LoadPartialNeuronNumericalResults <- function() {
  attenuationPartialNospines <- read.csv("NumericalResults/attenuation-partial-nospines.csv")
  attenuationPartialSpiny <- read.csv("NumericalResults/attenuation-partial-spiny.csv")
  mbpapPartialNospines <- read.csv("NumericalResults/mbpap-partial-nospines.csv")
  mbpapPartialSpiny <- read.csv("NumericalResults/mbpap-partial-spiny.csv")
  geometry <- read.csv("NumericalResults/geometry.csv")
  ages <- read.csv("NumericalResults/partialNeuronAges.csv")

  result <- merge(ages, attenuationPartialNospines, c("cellName"), all=TRUE)
  result <- merge(result, attenuationPartialSpiny,
                  c("cellName", "parameterSet"), all=TRUE)
  result <- merge(result, mbpapPartialNospines, c("cellName", "parameterSet"),
                  all=TRUE)
  result <- merge(result, mbpapPartialSpiny, c("cellName", "parameterSet"),
                  all=TRUE)
  result <- merge(result, geometry, c("cellName", "parameterSet"),
                  all=TRUE)

  return(result)
}

LoadShollNumericalResults <- function() {
  sholl <- read.csv("NumericalResults/sholl.csv")
  spineSholl <- read.csv("NumericalResults/spineSholl.csv")
  ages <- read.csv("NumericalResults/partialNeuronAges.csv")

  sholl$radius <- factor(sholl$radius)
  spineSholl$radius <- factor(spineSholl$radius)
  shollWidened <- reshape(sholl, idvar="cellName",
                          v.names=c("numberIntersections",
                                    "intersectionArea",
                                    "meanIntersectionDiameter"),
                          timevar="radius", direction="wide")
  spineShollWidened <- reshape(spineSholl, idvar="cellName",
                               v.names=c("spineCount","sectionLength",
                                         "spineDensity"),
                               timevar="radius", direction="wide")
  result <- merge(ages, shollWidened, c("cellName"), all=TRUE)
  result <- merge(result, spineShollWidened,
                  c("cellName", "parameterSet"), all=TRUE)

  return(result)
}

LoadParameterSpaceNumericalResults <- function() {
  parameterSpace <- read.csv("NumericalResults/parameterSpace.csv")
}

