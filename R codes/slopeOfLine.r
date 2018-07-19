#!/usr/bin/Rscript

args <- commandArgs(TRUE)
argsAsNumbers <- as.numeric(args)

current <- argsAsNumbers[c(TRUE, FALSE)]
voltage <- argsAsNumbers[c(FALSE, TRUE)]

voltageByCurrent <- data.frame(current, voltage)

model <- lm(voltage ~ current, data=voltageByCurrent)

cat(coefficients(model)["current"])
cat("\n")
