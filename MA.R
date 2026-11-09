# MA plots for soil metabolic features
# Place this script and the three tab-delimited text files in the same directory.
# Command-line execution uses the script directory automatically.
# In an interactive session, set the working directory to the file directory before running.
# Only base R functions are required; R 4.2 or later is recommended.
# Three plots are written to the MA_plots directory. Input data are not modified.

arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", arguments, value = TRUE)
if (length(file_argument) > 0L) {
  script_directory <- dirname(normalizePath(sub("^--file=", "", file_argument[[1L]]), mustWork = TRUE))
} else {
  script_directory <- getwd()
}

output_directory <- file.path(script_directory, "MA_plots")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

comparison_settings <- list(
  list(
    input_file = "KvsA.txt",
    output_file = "MA_K_vs_A.png",
    plot_title = "T. japonica rhizosphere soil vs bulk soil",
    rhizosphere_mean = "mean_T_japonica"
  ),
  list(
    input_file = "KvsB.txt",
    output_file = "MA_K_vs_B.png",
    plot_title = "E. quadrinervis rhizosphere soil vs bulk soil",
    rhizosphere_mean = "mean_E_quadrinervis"
  ),
  list(
    input_file = "KvsD.txt",
    output_file = "MA_K_vs_D.png",
    plot_title = "H. contortus rhizosphere soil vs bulk soil",
    rhizosphere_mean = "mean_H_contortus"
  )
)

common_fields <- c(
  "feature_id",
  "mean_bulk_soil",
  "p_value",
  "log2_fold_change",
  "significance"
)
significance_levels <- c("upregulated", "not_significant", "downregulated")
significance_labels <- c(
  upregulated = "Upregulated",
  not_significant = "Not significant",
  downregulated = "Downregulated"
)
point_colors <- c(
  upregulated = "#B31B21",
  not_significant = "#A9A9A9",
  downregulated = "#266B69"
)

for (settings in comparison_settings) {
  input_path <- file.path(script_directory, settings$input_file)
  if (!file.exists(input_path)) {
    stop(paste("Missing input file:", settings$input_file))
  }

  data <- read.delim(
    input_path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    fileEncoding = "UTF-8",
    stringsAsFactors = FALSE
  )

  required_fields <- unique(c(common_fields, settings$rhizosphere_mean))
  missing_fields <- setdiff(required_fields, names(data))
  if (length(missing_fields) > 0L) {
    stop(paste(settings$input_file, "is missing fields:", paste(missing_fields, collapse = ", ")))
  }
  if (anyDuplicated(data[["feature_id"]])) {
    stop(paste(settings$input_file, "contains duplicate feature identifiers"))
  }

  bulk_mean <- data[["mean_bulk_soil"]]
  rhizosphere_mean <- data[[settings$rhizosphere_mean]]
  y_coordinate <- data[["log2_fold_change"]]
  significance <- data[["significance"]]

  if (any(!is.finite(bulk_mean)) || any(!is.finite(rhizosphere_mean)) || any(!is.finite(y_coordinate))) {
    stop(paste(settings$input_file, "contains missing or non-finite numeric values"))
  }
  if (any(bulk_mean <= 0) || any(rhizosphere_mean <= 0)) {
    stop(paste(settings$input_file, "contains non-positive mean signal intensities"))
  }
  if (any(is.na(significance)) || any(!significance %in% significance_levels)) {
    stop(paste(settings$input_file, "contains unrecognized significance categories"))
  }

  x_coordinate <- (log2(bulk_mean) + log2(rhizosphere_mean)) / 2
  output_path <- file.path(output_directory, settings$output_file)

  png(output_path, width = 2400, height = 1800, res = 300)
  par(mar = c(5, 6, 4, 2) + 0.1, las = 1)
  plot(
    x = x_coordinate,
    y = y_coordinate,
    pch = 17,
    cex = 0.7,
    col = point_colors[["not_significant"]],
    main = settings$plot_title,
    xlab = "Mean signal intensity (log2)",
    ylab = "Fold change (log2)"
  )
  points(
    x = x_coordinate[significance == "downregulated"],
    y = y_coordinate[significance == "downregulated"],
    pch = 17,
    cex = 0.7,
    col = point_colors[["downregulated"]]
  )
  points(
    x = x_coordinate[significance == "upregulated"],
    y = y_coordinate[significance == "upregulated"],
    pch = 17,
    cex = 0.7,
    col = point_colors[["upregulated"]]
  )
  abline(h = 0, lty = 1, lwd = 2, col = "#000000")
  abline(h = c(-1, 1), lty = 2, lwd = 2, col = "#000000")
  legend(
    "topright",
    inset = 0.01,
    title = "Significance",
    legend = significance_labels[significance_levels],
    pch = 17,
    col = point_colors[significance_levels],
    pt.cex = 1.5,
    bty = "n"
  )
  dev.off()

  message(paste("Created:", output_path))
}
