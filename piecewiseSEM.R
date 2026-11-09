# Piecewise structural equation modeling: public analysis script
# Data file: piecewiseSEM.xlsx
# Prepared: 2026-11-08
# The script preserves the five composite variables and five structural paths
# used in the original analysis.

# 1. Prepare the computing environment
# Place this script and the workbook in the same directory.
# When the script is run directly, the workbook is located automatically.
# In an interactive session, set the working directory to that directory first.

required_packages <- c("piecewiseSEM", "readxl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Required packages are missing. Please install: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(piecewiseSEM)
library(readxl)

cat("piecewiseSEM version:", as.character(packageVersion("piecewiseSEM")), "\n")
cat("readxl version:", as.character(packageVersion("readxl")), "\n")

# 2. Locate and import the public dataset
script_arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", script_arguments, value = TRUE)

if (length(file_argument) > 0L) {
  script_path <- sub("^--file=", "", file_argument[[1L]])
  data_path <- file.path(dirname(normalizePath(script_path)), "piecewiseSEM.xlsx")
} else {
  data_path <- "piecewiseSEM.xlsx"
}

if (!file.exists(data_path)) {
  stop(
    "The workbook was not found. Place both public files in the same directory ",
    "and set that directory as the working directory before running again."
  )
}

# The selected range contains one English header row and sixteen data rows.
sem_data <- read_excel(
  data_path,
  sheet = "Raw_data",
  range = cell_limits(c(4, 2), c(20, 20))
)
sem_data <- as.data.frame(sem_data, check.names = FALSE)

# 3. Validate the data structure and completeness
expected_fields <- c(
  "Sample_ID",
  "Palmitoleic_acid_relative_value",
  "Linoleic_acid_relative_value",
  "Trehalose_relative_value",
  "Terrimonas_percent",
  "Unclassified_Alphaproteobacteria_percent",
  "Claroideoglomus_percent",
  "Glomus_percent",
  "Unranked_Claroideoglomeraceae_percent",
  "Unranked_Glomeraceae_percent",
  "Soil_pH",
  "Ammonium_N_mg_kg",
  "Available_P_mg_kg",
  "Total_N_g_kg",
  "Total_P_g_kg",
  "Urease_U_g",
  "Acid_phosphatase_U_g",
  "Phytase_U_g",
  "Aboveground_biomass_g_cm2"
)

if (!identical(names(sem_data), expected_fields)) {
  stop(
    "The workbook fields do not match the expected fields. ",
    "Use the public workbook supplied with this script without reordering columns."
  )
}
if (nrow(sem_data) != 16L) {
  stop("The dataset must contain sixteen samples. Check the selected data range.")
}
if (anyDuplicated(sem_data$Sample_ID)) {
  stop("Duplicate sample identifiers were found. Check Sample_ID in the workbook.")
}
if (any(!complete.cases(sem_data))) {
  stop(
    "Missing values were found. Check the source records; ",
    "this script does not impute missing values."
  )
}

analysis_fields <- expected_fields[-1L]
numeric_check <- vapply(sem_data[analysis_fields], is.numeric, logical(1))
if (any(!numeric_check)) {
  stop(
    "All analytical variables must be numeric. Check: ",
    paste(analysis_fields[!numeric_check], collapse = ", ")
  )
}

cat(
  "Data validation passed: sixteen samples, eighteen analytical variables, ",
  "and no blank cells.\n"
)

# 4. Define the composite-variable construction function
# Aboveground biomass is first regressed on all component variables.
# The composite score is the sum of each component multiplied by its
# unstandardized regression slope; the intercept is not added.
construct_composite <- function(data, predictors, composite_name) {
  cat("\n========================================\n")
  cat("Constructing", composite_name, "\n")

  initial_formula <- reformulate(
    predictors,
    response = "Aboveground_biomass_g_cm2"
  )
  initial_model <- lm(initial_formula, data = data)

  cat("\nInitial linear model summary:\n")
  print(summary(initial_model))
  cat("\nStandardized coefficients for the initial model:\n")
  print(coefs(initial_model, standardize = "scale"))

  regression_slopes <- coef(initial_model)[predictors]
  if (anyNA(regression_slopes)) {
    stop(
      composite_name,
      " has regression slopes that cannot be estimated. ",
      "Check collinearity among its components."
    )
  }

  data[[composite_name]] <- as.vector(
    as.matrix(data[predictors]) %*% regression_slopes
  )

  composite_formula <- reformulate(
    composite_name,
    response = "Aboveground_biomass_g_cm2"
  )
  composite_model <- lm(composite_formula, data = data)
  cat("\nModel summary for aboveground biomass predicted by the composite:\n")
  print(summary(composite_model))
  cat("\nStandardized coefficients for the composite model:\n")
  print(coefs(composite_model, standardize = "scale"))

  return(data)
}

# 5. Construct the five composite variables in the original order
sem_data <- construct_composite(
  sem_data,
  c("Terrimonas_percent", "Unclassified_Alphaproteobacteria_percent"),
  "Composite_1"
)

sem_data <- construct_composite(
  sem_data,
  c(
    "Claroideoglomus_percent",
    "Glomus_percent",
    "Unranked_Claroideoglomeraceae_percent",
    "Unranked_Glomeraceae_percent"
  ),
  "Composite_2"
)

sem_data <- construct_composite(
  sem_data,
  c(
    "Soil_pH",
    "Ammonium_N_mg_kg",
    "Available_P_mg_kg",
    "Total_N_g_kg",
    "Total_P_g_kg"
  ),
  "Composite_3"
)

sem_data <- construct_composite(
  sem_data,
  c("Urease_U_g", "Acid_phosphatase_U_g", "Phytase_U_g"),
  "Composite_4"
)

sem_data <- construct_composite(
  sem_data,
  c(
    "Palmitoleic_acid_relative_value",
    "Linoleic_acid_relative_value",
    "Trehalose_relative_value"
  ),
  "Composite_5"
)

# 6. Fit the piecewise structural equation model
# No paths are added to or removed from the original model.
sem_model <- psem(
  lm(
    Aboveground_biomass_g_cm2 ~ Composite_1 + Composite_2 + Composite_4 + Composite_3,
    data = sem_data
  ),
  lm(Composite_1 ~ Composite_5 + Composite_2, data = sem_data),
  lm(Composite_2 ~ Composite_5, data = sem_data),
  lm(Composite_3 ~ Composite_1 + Composite_2 + Composite_5, data = sem_data),
  lm(
    Composite_4 ~ Composite_1 + Composite_2 + Composite_5 + Composite_3,
    data = sem_data
  )
)

model_summary <- summary(sem_model, .progressBar = FALSE)
cat("\n========================================\n")
cat("Piecewise structural equation model summary:\n")
print(model_summary)

# 7. Extract model fit statistics
aic_value <- if (is.list(model_summary$AIC)) {
  model_summary$AIC[[1L]]
} else {
  as.numeric(model_summary$AIC)
}

# Calculate BIC directly when supported by the installed package version.
# Otherwise, use AIC with a sample-size penalty; no approximate parameter
# count is substituted.
calculate_bic <- function(model, data) {
  tryCatch(
    {
      result <- BIC(model)
      if (is.list(result)) result[[1L]] else as.numeric(result)
    },
    error = function(error) {
      cat("Direct BIC calculation failed; recalculating with a sample-size penalty.\n")
      result <- AIC(model, k = log(nrow(data)))
      if (is.list(result)) result[[1L]] else as.numeric(result)
    }
  )
}

bic_value <- calculate_bic(sem_model, sem_data)

cat("\nModel fit statistics:\n")
cat("AIC:", aic_value, "\n")
cat("BIC:", bic_value, "\n")
cat(
  "Chi-square:", model_summary$Chisq$Chisq,
  "; df:", model_summary$Chisq$df,
  "; p-value:", model_summary$Chisq$P.Value,
  "\n"
)
cat(
  "Fisher's C:", model_summary$Cstat$Fisher.C,
  "; df:", model_summary$Cstat$df,
  "; p-value:", model_summary$Cstat$P.Value,
  "\n"
)

# 8. Record the computing environment for reproducibility
cat("\nSession information:\n")
print(sessionInfo())
