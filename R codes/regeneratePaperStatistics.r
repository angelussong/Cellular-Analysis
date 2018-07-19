#!/usr/bin/Rscript

# This script assumes it's being read from the top-level paper directory (that
# is, the one containing the Scripts/ directory.

# Clear any pre-existing environment
rm(list=ls())

cat("\n")
options(width=150)

source("Scripts/R/lib/CorrelationFunctions.r")
source("Scripts/R/lib/Filters.r")
source("Scripts/R/lib/LoadData.r")
source("Scripts/R/lib/Tests.r")

data <- list()
data$partial <- LoadPartialNeuronNumericalResults()
data$apical <- getApicalNeuronsFromFrame(data$partial)
data$basal <- getBasalNeuronsFromFrame(data$partial)
data$sholl <- LoadShollNumericalResults()
data$whole <- LoadWholeNeuronData(data$apical, data$basal)

#   May3x is excluded from the data set because, although electrophysiological
# measurements were taken, the morphology was insufficiently clear to be
# reconstructed accurately.
#   The others were excluded from the original analysis for assorted reasons;
# generally poor reconstruction or electrophysiological measurements.
data <- removeNeuronWithName(data, "May3x")
data <- removeNeuronWithName(data, "Apr6f")
data <- removeNeuronWithName(data, "Dec15d")
data <- removeNeuronWithName(data, "May4_2000")

data <- splitYoungAndOld(data)
data$parameterSpace <- LoadParameterSpaceNumericalResults()

# Verify that the soma depths are statistically different.
testSomaDepths(data$wholeYoung, data$wholeOld)

# Remove neurons such that the p-value of the difference in soma depths is about 0.5.
data <- removeNeuronWithName(data, "Aug3b")
data <- removeNeuronWithName(data, "May3c")
data <- removeNeuronWithName(data, "May3g")
data <- removeNeuronWithName(data, "May3j")
data <- removeNeuronWithName(data, "May31e")

cat("Young neurons: ")
print(data$wholeYoung$cellName)
cat("Old neurons: ")
print(data$wholeOld$cellName)
cat("\n\n")

# Verify that the soma depths are no longer statistically different.
cat("After discarding outlying neurons until p is roughly 0.5:\n")
testSomaDepths(data$wholeYoung, data$wholeOld)

testCongruentDataFrames(data$wholeYoung, data$wholeOld,
                        "NumericalResults/wholeCellTTest.csv")
testCongruentDataFrames(data$apicalYoung, data$apicalOld,
                        "NumericalResults/apicalCellTTest.csv")
testCongruentDataFrames(data$basalYoung, data$basalOld,
                        "NumericalResults/basalCellTTest.csv")

#---  The following can probably be moved into the Correlations function.

Correlations(data$whole, 'r-cortests-whole-all.csv',
             'correlation-graphs-whole-all.pdf')
Correlations(data$wholeYoung, 'r-cortests-whole-young.csv',
             'correlation-graphs-whole-young.pdf')
Correlations(data$wholeOld, 'r-cortests-whole-old.csv',
             'correlation-graphs-whole-old.pdf')
