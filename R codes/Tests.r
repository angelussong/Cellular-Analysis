testSomaDepths <- function(wholeYoung, wholeOld) {
  young <- wholeYoung$distanceToPialSurface.CoskrenFromStacks
  old <- wholeOld$distanceToPialSurface.CoskrenFromStacks

  stats <- t.test(young, old)
  youngMean <- mean(young)
  oldMean <- mean(old)
  pValue <- stats["p.value"][[1]]
  confidenceIntervalLeft = stats[["conf.int"]][1]
  confidenceIntervalRight = stats[["conf.int"]][2]

  # sink(output)
  cat("Step 1: Comparing depth of young and old neurons\n")
  cat("Mean depth of young neurons: ", youngMean, "\n")
  cat("Mean depth of old neurons: ", oldMean, "\n")
  cat("p-value (means differ): ", pValue, "\n")
  cat("\n")
  # sink()
}

# Given two data frames that consist of same-named and same-typed columns, this
# function considers each data fram a separate population, and prints t-test
# data for each numeric column.
testCongruentDataFrames <- function(wholeYoung, wholeOld, outputPath) {
  numericColumns <- sapply(wholeYoung, is.numeric)
  wholeYoungNumeric <- wholeYoung[, numericColumns]
  wholeOldNumeric <- wholeOld[, numericColumns]

  youngMeans <- lapply(wholeYoungNumeric, mean, na.rm=TRUE)
  oldMeans <- lapply(wholeOldNumeric, mean, na.rm=TRUE)
  youngStandardDeviations <- sapply(wholeYoungNumeric, sd, na.rm=TRUE)
  oldStandardDeviations <- sapply(wholeOldNumeric, sd, na.rm=TRUE)
  youngStandardErrors <- lapply(wholeYoungNumeric, stderr)
  oldStandardErrors <- lapply(wholeOldNumeric, stderr)

  stats <- data.frame(youngMeans)
  stats <- rbind(stats, youngStandardDeviations)
  stats <- rbind(stats, youngStandardErrors)
  stats <- rbind(stats, oldMeans)
  stats <- rbind(stats, oldStandardDeviations)
  stats <- rbind(stats, oldStandardErrors)

  testResults <- lapply(
      names(wholeYoungNumeric),
      function(x) {
        t.test(wholeYoungNumeric[[x]], wholeOldNumeric[[x]])
      })
  pValues <- lapply(testResults, function(x) { x["p.value"][[1]] })
  confidenceIntervalLeftValues <- lapply(testResults, function(x) { x[["conf.int"]][[1]] })
  confidenceIntervalRightValues <- lapply(testResults, function(x) { x[["conf.int"]][[2]] })
  stats <- rbind(stats, pValues)
  stats <- rbind(stats, confidenceIntervalLeftValues)
  stats <- rbind(stats, confidenceIntervalRightValues)

  testResults <- lapply(
      names(wholeYoungNumeric),
      function(x) {
        wilcox.test(wholeYoungNumeric[[x]], wholeOldNumeric[[x]], exact=FALSE)
      })
  pValues <- lapply(testResults, function(x) { x["p.value"][[1]] })
  stats <- rbind(stats, pValues)

  stat = c("Young mean", "Young std dev", "Young std err", "Old mean",
           "Old std dev", "Old std err", "t.test p-value",
           "conf. int. left", "conf. int. right", "Wilcox p-value")
  stats <- cbind(stat, stats)

  outputStream <- file(outputPath, "w")
  sink(outputStream)
  write.table(stats, sep=",", row.names=FALSE, quote=FALSE)
  sink()
  close(outputStream)
}
