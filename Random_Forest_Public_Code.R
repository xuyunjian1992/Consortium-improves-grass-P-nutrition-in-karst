# Random Forest Analysis Public Code
# Purpose: Reproduce random forest analyses and variable-importance plots for
# commonly enriched bacteria and fungi in relation to aboveground biomass.
# Data file: Random_Forest_Public_Data.xlsx
# Usage: Place this script and the data workbook in the same directory,
# set that directory as the working directory, and run the entire script.
# The random seed, number of trees, number of permutations, and significance
# calculations are retained from the source script.

# Load required packages ---------------------------------------------------
library(randomForest)
library(ggplot2)
library(openxlsx)

# File names and analysis parameters ---------------------------------------
data_file <- "Random_Forest_Public_Data.xlsx"
output_directory <- "analysis_results"
response_variable <- "AGB_g_cm_2"
random_seed <- 123
training_fraction <- 0.7
number_of_trees <- 500
number_of_permutations <- 100

if (!file.exists(data_file)) {
  stop("Data file not found. Place this script and Random_Forest_Public_Data.xlsx in the working directory.")
}
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

# Column ranges select the populated data area and exclude blank margins.
sheet_settings <- list(
  list(name = "Bacterial_Data", columns = 2:20, label = "Bacterial", expected_columns = 18),
  list(name = "Fungal_Data", columns = 2:6, label = "Fungal", expected_columns = 4)
)

# Analyze bacterial and fungal datasets sequentially ------------------------
for (current_setting in sheet_settings) {
  source_data <- read.xlsx(
    xlsxFile = data_file,
    sheet = current_setting$name,
    startRow = 4,
    cols = current_setting$columns,
    check.names = FALSE
  )

  # Retain complete rows. The public sample identifier is not a predictor.
  complete_data <- source_data[complete.cases(source_data), , drop = FALSE]
  removed_rows <- nrow(source_data) - nrow(complete_data)
  if (removed_rows > 0) {
    cat(current_setting$name, ": removed", removed_rows, "rows with missing values.\n")
  }
  analysis_data <- complete_data[, -1, drop = FALSE]

  # Validate the data structure, response variable, and numeric columns.
  if (ncol(analysis_data) != current_setting$expected_columns) {
    stop(current_setting$name, " has an unexpected number of columns. Use the accompanying public workbook.")
  }
  if (!response_variable %in% names(analysis_data)) {
    stop("Response variable not found: ", response_variable)
  }
  if (anyDuplicated(names(analysis_data))) {
    stop(current_setting$name, " contains duplicate variable names.")
  }
  if (nrow(analysis_data) < 4) {
    stop(current_setting$name, " has too few complete samples for training and testing.")
  }
  numeric_check <- vapply(analysis_data, is.numeric, logical(1))
  if (!all(numeric_check)) {
    stop(current_setting$name, " contains non-numeric variables.")
  }

  # Split the data into training and test sets.
  set.seed(random_seed)
  training_rows <- sample(
    seq_len(nrow(analysis_data)),
    size = floor(training_fraction * nrow(analysis_data)),
    replace = FALSE
  )
  training_data <- analysis_data[training_rows, , drop = FALSE]
  test_data <- analysis_data[-training_rows, , drop = FALSE]
  model_formula <- as.formula(paste(response_variable, "~ ."))

  # Retain five candidate predictors where available; use all predictors for fungi.
  predictor_count <- ncol(analysis_data) - 1
  candidate_predictors <- min(5, predictor_count)

  # Train the source model and evaluate test-set performance.
  source_model <- randomForest(
    model_formula,
    data = training_data,
    ntree = number_of_trees,
    mtry = candidate_predictors,
    importance = TRUE,
    na.action = na.omit
  )
  source_predictions <- predict(source_model, newdata = test_data)
  observed_values <- test_data[[response_variable]]
  mean_squared_error <- mean((observed_values - source_predictions)^2)
  r_squared <- cor(observed_values, source_predictions)^2
  cat(
    current_setting$name,
    ": MSE =", mean_squared_error,
    "; R_squared =", r_squared, "\n"
  )

  # Extract importance based on node purity, consistent with the source script.
  observed_importance <- importance(source_model, type = 2, scale = TRUE)[, 1]
  feature_count <- length(observed_importance)
  permuted_r_squared <- numeric(number_of_permutations)
  permuted_importance <- matrix(nrow = feature_count, ncol = number_of_permutations)
  rownames(permuted_importance) <- names(observed_importance)

  # Permute the response variable to assess model and variable significance.
  for (permutation_index in seq_len(number_of_permutations)) {
    permuted_training_data <- training_data
    permuted_training_data[[response_variable]] <- sample(permuted_training_data[[response_variable]])
    permuted_model <- randomForest(
      model_formula,
      data = permuted_training_data,
      ntree = number_of_trees,
      mtry = candidate_predictors,
      importance = TRUE,
      na.action = na.omit
    )
    permuted_predictions <- predict(permuted_model, newdata = test_data)
    permuted_r_squared[permutation_index] <- cor(observed_values, permuted_predictions)^2
    permuted_importance[, permutation_index] <- importance(
      permuted_model, type = 2, scale = TRUE
    )[, 1]
  }

  model_p_value <- mean(permuted_r_squared >= r_squared)
  variable_p_values <- vapply(
    seq_len(feature_count),
    function(variable_index) {
      mean(permuted_importance[variable_index, ] >= observed_importance[variable_index])
    },
    numeric(1)
  )
  cat(current_setting$name, ": model p-value =", model_p_value, "\n")

  # Sort and normalize variable importance using the source calculation order.
  importance_table <- data.frame(
    Microorganism = names(observed_importance),
    Importance = as.numeric(observed_importance),
    P_value = variable_p_values,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  importance_table <- importance_table[order(-importance_table$Importance), , drop = FALSE]
  total_importance <- sum(importance_table$Importance)
  if (total_importance <= 0) {
    stop(current_setting$name, " has a non-positive total importance and cannot be normalized.")
  }
  importance_table$Normalized_importance_percent <- round(
    importance_table$Importance / total_importance * 100, 2
  )
  importance_table$Importance <- NULL
  importance_table <- importance_table[, c("Microorganism", "Normalized_importance_percent", "P_value")]

  performance_table <- data.frame(
    Data_type = current_setting$label,
    Sample_count = nrow(analysis_data),
    Training_sample_count = nrow(training_data),
    Test_sample_count = nrow(test_data),
    Number_of_trees = number_of_trees,
    Candidate_predictors = candidate_predictors,
    Permutation_count = number_of_permutations,
    MSE = mean_squared_error,
    R_squared = r_squared,
    Model_p_value = model_p_value,
    check.names = FALSE
  )

  # Save numerical results for verification and reuse.
  result_file <- file.path(
    output_directory,
    paste0(current_setting$label, "_random_forest_results.xlsx")
  )
  write.xlsx(
    list("Variable_Importance" = importance_table, "Model_Performance" = performance_table),
    file = result_file,
    rowNames = FALSE,
    overwrite = TRUE
  )

  # Plot variable importance. An asterisk denotes a p-value below 0.05.
  importance_plot <- ggplot(
    importance_table,
    aes(x = reorder(Microorganism, Normalized_importance_percent), y = Normalized_importance_percent)
  ) +
    geom_col(fill = "#4682B4", width = 0.7) +
    geom_text(
      aes(label = ifelse(P_value < 0.05, "*", "")),
      size = 6,
      color = "#D73027",
      vjust = 0.5,
      hjust = -0.2
    ) +
    coord_flip() +
    labs(
      title = paste(current_setting$label, "random forest variable importance"),
      x = "Core microorganism",
      y = "Normalized importance (%)"
    ) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5))

  # Print explicitly so that the graphics device receives the complete plot.
  figure_file <- file.path(
    output_directory,
    paste0(current_setting$label, "_random_forest_variable_importance.pdf")
  )
  cairo_pdf(figure_file, width = 8, height = 6)
  print(importance_plot)
  dev.off()

  cat("Created:", result_file, "and", figure_file, "\n")
}

cat("All analyses completed. Results are in the analysis_results directory.\n")
