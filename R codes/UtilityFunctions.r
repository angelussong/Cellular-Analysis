stderr <- function(x) {
  return(sd(x, na.rm=TRUE) / sqrt(length(x[!is.na(x)])))
}

endsWith <- function(testString, subString) {
  testLength = nchar(testString)
  subLength = nchar(subString)
  return(
      substr(testString, (testLength - subLength) + 1, testLength) == subString)
}
