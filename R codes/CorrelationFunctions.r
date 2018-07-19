source("Scripts/R/lib/LoadData.r")

Correlations <- function(wholeData, output_file, pdf_file) {
  all_morpho_data <- wholeData[, MorphologicalColumns()]
  all_ephys_data <- wholeData[, EphysColumns()]

  # One point type per row, with each row's type corresponding to its age group.
  point_types <- numeric(length(wholeData$ageCategory))
  point_types[wholeData$ageCategory == "OLD"] <- 21
  point_types[wholeData$ageCategory == "YOUNG"] <- 19

  correlations_output_file <- paste("NumericalResults/", output_file,
                                    "_correlations", sep="")
  output_file <- paste("NumericalResults/", output_file, sep="")
  pdf_file <- paste("GeneratedGraphs/", pdf_file, sep="")

  legend_pos <- c(600, 40)
  legend_args <- c('old', 'young')
  legend_symbols <- c(21, 19)

  colors <- c("black")

  morpho_column_names = colnames(all_morpho_data)
  ephys_column_names = colnames(all_ephys_data)

  p_values <- numeric()
  correlations <- numeric()
  pdf(file=pdf_file)
  par(fin = c(5, 5), fig = c(0, 1, 0, 1))
  first_plot <- TRUE
  for (each_morpho_name in morpho_column_names) {
    for (each_ephys_name in ephys_column_names) {
      morpho_data <- all_morpho_data[[each_morpho_name]]
      ephys_data <- all_ephys_data[[each_ephys_name]]
      if ((length(morpho_data[morpho_data != 0 && !is.na(morpho_data)]) > 0) &&
          (length(ephys_data[ephys_data != 0 && !is.na(ephys_data)]) > 0)) {
        p_value <- cor.test(morpho_data, ephys_data)$p.value
        p_values <- c(p_values, p_value)
        correlation <- cor(morpho_data, ephys_data, use="complete.obs")
        correlations <- c(correlations, correlation)
        # Consider plotting the data regardless of the p-value
        if (!is.na(p_value) && p_value < 0.05) {
          bestfit <- line(morpho_data, ephys_data)
          fitCoefficients <- coef(bestfit)
          plot(morpho_data, ephys_data, xlab=each_morpho_name,
               ylab=each_ephys_name, bty="l", pch=point_types, col=colors)
          if (!is.na(fitCoefficients[1]) && !is.na(fitCoefficients[2])) {
            abline(fitCoefficients)
          } else {
            cat ("Could not plot best-fit line for ", each_morpho_name, " vs ",
                 each_ephys_name, ".\n",
                 sep="")
          }
          if (first_plot) {
            first_plot <- FALSE
            # This should be generalized, but for now it'll do.
            if (0 != length(legend_args)) {
               legend(legend_pos[1], legend_pos[2], legend=legend_args,
                      pch=legend_symbols)
            }
          }
        }
      } else {
        cat("Could not compute correlation for ", each_morpho_name, " vs ",
            each_ephys_name, ".  Insufficient data.\n", sep="")
        p_values <- c(p_values, NA)
        correlations <- c(correlations, NA)
      }
    }
  }
  dev.off()
  dim(p_values) <- c(length(ephys_column_names), length(morpho_column_names))
  dim(correlations) <- c(length(ephys_column_names), length(morpho_column_names))
  # Morphological measures are rows, electrophysiology measures are columns
  p_values <- t(p_values)
  correlations <- t(correlations)
  rownames(p_values) <- MorphologicalColumns()
  colnames(p_values) <- EphysColumns()
  write.csv(p_values, row.names=TRUE, quote=FALSE,
            file=output_file)
  write.csv(correlations, row.names=TRUE, quote=FALSE,
            file=correlations_output_file)
}
