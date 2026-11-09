# Reproducible microbial network analysis for public data release
# Place this script and both workbooks in the same directory.
# Set that directory as the working directory before running the script.

# Load the required packages
library(psych)
library(igraph)
library(readxl)

# 1. Define input files and analysis parameters
abundance_file <- "network.xlsx"
taxonomy_file <- "tax_table.xlsx"
abundance_sheet <- "bacteria-ASVs"
taxonomy_sheet <- "bacteria-phylum"
sample_columns <- c("A1", "A2", "A3", "A4")
id_field <- "ASV_ID"
retain_count <- 400
correlation_threshold <- 0.6
adjusted_p_threshold <- 0.05

# For fungi, use fungi-ASVs and fungi-phylum.
# For another sample group, replace the four sample identifiers above.

# 2. Check input files
required_files <- c(abundance_file, taxonomy_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(paste("Missing input files:", paste(missing_files, collapse = ", ")))
}

# 3. Import and validate abundance data
raw_abundance <- data.frame(
  read_excel(abundance_file, sheet = abundance_sheet, .name_repair = "minimal"),
  check.names = FALSE
)
required_fields <- c(id_field, sample_columns)
missing_fields <- setdiff(required_fields, names(raw_abundance))
if (length(missing_fields) > 0) {
  stop(paste("Missing abundance fields:", paste(missing_fields, collapse = ", ")))
}
if (anyDuplicated(raw_abundance[[id_field]])) {
  stop("Duplicate ASV identifiers were found in the abundance table.")
}

abundance <- raw_abundance[sample_columns]
abundance[] <- lapply(abundance, function(current_column) as.numeric(current_column))
if (anyNA(abundance)) {
  stop("The abundance data contain missing or non-numeric values.")
}
rownames(abundance) <- raw_abundance[[id_field]]

# 4. Retain dominant ASVs by total abundance across the selected samples
total_abundance <- rowSums(abundance)
abundance_order <- order(-total_abundance)
actual_retain_count <- min(retain_count, nrow(abundance))
top_abundance <- abundance[
  abundance_order[seq_len(actual_retain_count)],
  ,
  drop = FALSE
]

# 5. Calculate Spearman correlations with false discovery rate adjustment
correlation_results <- corr.test(
  t(top_abundance),
  method = "spearman",
  adjust = "fdr",
  alpha = 0.05
)
correlation_matrix <- correlation_results$r
adjusted_p_matrix <- correlation_results$p

# 6. Retain edges that meet significance and correlation-strength criteria
correlation_matrix[
  adjusted_p_matrix > adjusted_p_threshold |
    abs(correlation_matrix) < correlation_threshold
] <- 0
diag(correlation_matrix) <- 0

# 7. Build an undirected weighted network and remove isolated nodes
network <- graph_from_adjacency_matrix(
  correlation_matrix,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)
network <- simplify(network)
isolated_nodes <- which(degree(network) == 0)
network <- delete_vertices(network, isolated_nodes)

# 8. Add signed correlations, absolute weights, edge directions, and phyla
E(network)$correlation <- E(network)$weight
E(network)$weight <- abs(E(network)$correlation)
E(network)$type <- ifelse(E(network)$correlation > 0, "positive", "negative")

taxonomy <- data.frame(
  read_excel(taxonomy_file, sheet = taxonomy_sheet, .name_repair = "minimal"),
  check.names = FALSE
)
if (!all(c(id_field, "phylum") %in% names(taxonomy))) {
  stop("The taxonomy table must contain ASV identifiers and phylum annotations.")
}
if (anyDuplicated(taxonomy[[id_field]])) {
  stop("Duplicate ASV identifiers were found in the taxonomy table.")
}
rownames(taxonomy) <- taxonomy[[id_field]]
missing_taxonomy <- setdiff(V(network)$name, rownames(taxonomy))
if (length(missing_taxonomy) > 0) {
  stop(paste("Nodes without taxonomy:", paste(missing_taxonomy, collapse = ", ")))
}
V(network)$phylum <- taxonomy[V(network)$name, "phylum"]

# 9. Prepare node and edge tables
node_data <- data.frame(
  id = V(network)$name,
  phylum = V(network)$phylum,
  degree = degree(network),
  betweenness = betweenness(network),
  row.names = NULL,
  check.names = FALSE
)

edge_endpoints <- ends(network, es = E(network))
edge_data <- data.frame(
  source = edge_endpoints[, 1],
  target = edge_endpoints[, 2],
  weight = E(network)$weight,
  correlation = E(network)$correlation,
  type = E(network)$type,
  row.names = NULL,
  check.names = FALSE
)

# 10. Export reusable node and edge tables
write.csv(
  node_data,
  "network_nodes.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  edge_data,
  "network_edges.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# The script creates network_nodes.csv and network_edges.csv in the working directory.
