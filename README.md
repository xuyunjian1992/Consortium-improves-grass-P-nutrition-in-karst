# Data and code for selected analyses of rhizosphere bacterial–AMF associations under extreme phosphorus limitation

This repository provides analysis-ready data and R scripts associated with the manuscript:

**A rhizosphere bacterial–AMF consortium improves phosphorus nutrition and growth of a native grass under extreme phosphorus limitation in karst rocky desertification**

**Authors:** Yunjian Xu, Weihao Cai, Zhe Chen*, and Xingwu Duan*  
**Affiliation:** Yunnan University, Kunming, China  
**Manuscript year:** 2026

## Study overview

Karst rocky desertification is characterized by extremely low soil phosphorus availability, which constrains vegetation restoration. The associated study investigated three dominant grasses—*Themeda japonica*, *Eulalia quadrinervis*, and *Heteropogon contortus*—in a karst rocky desertification area in Yunnan, Southwest China, using bulk soil as a control. The study integrated soil physicochemical and enzyme assays, microbial amplicon sequencing, untargeted metabolomics, microbial network analysis, random forest analysis, piecewise structural equation modeling, and inoculation experiments.

The materials released here support selected analyses and figure-related workflows, including:

- MA plots for three rhizosphere-versus-bulk-soil metabolic-feature comparisons;
- bacterial and fungal co-occurrence network construction and export of node/edge tables;
- exploratory random forest models relating commonly enriched microbial genera to aboveground biomass; and
- piecewise structural equation modeling of relationships among core metabolites, microbial taxa, soil properties, enzyme activities, and aboveground biomass.

> **Scope:** These files contain processed, analysis-ready data rather than raw sequencing reads or raw mass-spectrometry files. They do not constitute code and data for every analysis or figure in the manuscript.

## Repository contents

| File | Description |
|---|---|
| `MA.R` | Generates MA plots for three rhizosphere-versus-bulk-soil metabolic-feature comparisons using base R. |
| `all_analytical_data（MA）.xlsx` | Source workbook for the MA plots. Contains documentation, a data dictionary, and the worksheets `K_vs_A`, `K_vs_B`, and `K_vs_D`. |
| `network.R` | Constructs an undirected weighted microbial network and exports reusable node and edge tables. |
| `network.xlsx` | Bacterial and fungal ASV abundance tables for 16 samples. Contains `bacteria-ASVs` and `fungi-ASVs` worksheets. |
| `Microbial Network Analysis Taxonomic Data.xlsx` | Phylum-level annotations and complete taxonomic/abundance records for bacterial and fungal ASVs. |
| `Random_Forest_Public_Code.R` | Runs bacterial and fungal random forest analyses, response-permutation tests, and variable-importance plotting. |
| `Random_Forest_Public_Data.xlsx` | Twelve-sample random forest dataset containing aboveground biomass, 17 commonly enriched bacterial genera, and three commonly enriched fungal genera. |
| `piecewiseSEM.R` | Constructs five composite variables, fits the piecewise structural equation model, and reports model-fit statistics. |
| `piecewiseSEM.xlsx` | Sixteen-sample SEM dataset with 18 analytical variables, documentation, a data dictionary, and model specifications. |

Each workbook also contains its own documentation and/or data-dictionary worksheets. These embedded descriptions should be consulted for variable definitions, units, and coding rules.

## Sample identifiers

The microbial abundance workbooks use the following group and replicate identifiers:

| Identifier | Soil group |
|---|---|
| `A1`–`A4` | *Themeda japonica* rhizosphere soil |
| `B1`–`B4` | *Eulalia quadrinervis* rhizosphere soil |
| `D1`–`D4` | *Heteropogon contortus* rhizosphere soil |
| `K1`–`K4` | Unplanted bulk-soil control |

For the MA comparisons, `KvsA`, `KvsB`, and `KvsD` denote *T. japonica*, *E. quadrinervis*, and *H. contortus* rhizosphere soil relative to bulk soil, respectively.

The random forest and piecewise SEM workbooks use sequential public identifiers such as `Sample_01`. These identifiers preserve source row order only; they are not experimental group labels and are not used as predictors.

## Software requirements

The original manuscript analyses used **R 4.4.1**, which is recommended for reproducing the workflows in this repository. The MA script itself requires only base R and recommends R 4.2 or later.

Install the required R packages with:

```r
install.packages(c(
  "readxl",
  "psych",
  "igraph",
  "randomForest",
  "ggplot2",
  "openxlsx",
  "piecewiseSEM"
))
```

Package versions are not locked in this release. Numerical results and plotting details may vary slightly across R and package versions. The SEM script prints its package versions and full `sessionInfo()` output; recording `sessionInfo()` is also recommended for the other workflows.

## General usage

1. Download or clone the repository.
2. Keep each R script in the same directory as its required input files.
3. Set that directory as the R working directory, or run the scripts with `Rscript` from that directory.
4. Do not reorder workbook columns or rename worksheets unless the corresponding script is updated accordingly.

The input workbooks are never modified by the scripts.

## 1. Metabolic-feature MA plots

### Input preparation

`MA.R` reads three UTF-8, tab-delimited text files, whereas the public source data are supplied in `all_analytical_data（MA）.xlsx`. Before running the script, export the populated analytical range from each worksheet as follows:

| Workbook worksheet | Required text filename | Comparison |
|---|---|---|
| `K_vs_A` | `KvsA.txt` | *T. japonica* rhizosphere vs bulk soil |
| `K_vs_B` | `KvsB.txt` | *E. quadrinervis* rhizosphere vs bulk soil |
| `K_vs_D` | `KvsD.txt` | *H. contortus* rhizosphere vs bulk soil |

In each worksheet, the analytical table begins at **row 5, columns B–G**. Export the header and all data rows as tab-delimited UTF-8 text without adding blank leading rows or columns. Each resulting file must contain these six fields:

- `feature_id`
- `mean_bulk_soil`
- the comparison-specific rhizosphere mean (`mean_T_japonica`, `mean_E_quadrinervis`, or `mean_H_contortus`)
- `p_value`
- `log2_fold_change`
- `significance`

The permitted significance categories are `upregulated`, `not_significant`, and `downregulated`.

### Run

```bash
Rscript MA.R
```

No additional R packages are required. The script validates required fields, identifiers, numeric values, and significance categories.

### Output

Three 300-dpi PNG files are written to `MA_plots/`:

- `MA_K_vs_A.png`
- `MA_K_vs_B.png`
- `MA_K_vs_D.png`

The x-axis is mean signal intensity on the log2 scale, and the y-axis is the log2 fold change for rhizosphere soil relative to bulk soil.

## 2. Bacterial and fungal co-occurrence networks

### Input preparation

The supplied taxonomy workbook is named `Microbial Network Analysis Taxonomic Data.xlsx`, while `network.R` expects `tax_table.xlsx`. Before running the unmodified script, rename or copy the taxonomy workbook to:

```text
tax_table.xlsx
```

Alternatively, update the `taxonomy_file` assignment in `network.R` to the supplied filename.

### Default analysis

```bash
Rscript network.R
```

By default, the script analyzes bacterial samples `A1`–`A4` using:

- abundance worksheet: `bacteria-ASVs`;
- taxonomy worksheet: `bacteria-phylum`;
- dominant ASVs retained: the top 400 by total abundance across selected samples;
- correlation method: Spearman correlation;
- multiple-testing correction: false discovery rate;
- edge-retention thresholds: absolute correlation ≥ 0.6 and adjusted *P* ≤ 0.05.

Isolated nodes are removed. Node tables include phylum, degree, and betweenness; edge tables include source, target, absolute weight, signed correlation, and association type.

### Other microbial groups and sample sets

Edit the following assignments near the beginning of `network.R`:

- For fungi, use `abundance_sheet <- "fungi-ASVs"` and `taxonomy_sheet <- "fungi-phylum"`.
- For another soil group, replace `sample_columns` with `B1`–`B4`, `D1`–`D4`, or `K1`–`K4` as appropriate.

Run different group analyses in separate working directories or rename the outputs between runs, because the default output filenames are reused.

### Output

The script writes the following files to the working directory:

- `network_nodes.csv`
- `network_edges.csv`

These tables can be imported into network-visualization software such as Gephi.

## 3. Random forest analysis

### Run

Keep `Random_Forest_Public_Code.R` and `Random_Forest_Public_Data.xlsx` in the working directory, then run:

```bash
Rscript Random_Forest_Public_Code.R
```

The script analyzes the `Bacterial_Data` and `Fungal_Data` worksheets sequentially. Aboveground biomass (`AGB_g_cm_2`) is the response variable; the public sample identifier is excluded from the predictors.

The fixed analytical settings are:

- random seed: `123`;
- training fraction: `0.7`;
- number of trees: `500`;
- response permutations: `100`;
- candidate predictors at each split: five where available, otherwise all predictors.

Rows with missing values are removed before model fitting. Variable importance is based on node purity, and significance is assessed by response-permutation tests.

### Output

Results are written to `analysis_results/`:

- `Bacterial_random_forest_results.xlsx`
- `Bacterial_random_forest_variable_importance.pdf`
- `Fungal_random_forest_results.xlsx`
- `Fungal_random_forest_variable_importance.pdf`

Each results workbook contains variable-importance and model-performance worksheets. An asterisk in a PDF indicates a variable-level permutation *P* value below 0.05.

> **Interpretation:** Consistent with the manuscript, these random forest models are exploratory and are used to prioritize candidate taxa. They should not be interpreted as establishing robust predictive models of aboveground biomass.

## 4. Piecewise structural equation modeling

### Run

Keep `piecewiseSEM.R` and `piecewiseSEM.xlsx` in the same directory, then run:

```bash
Rscript piecewiseSEM.R
```

The script reads `Raw_data` over the fixed range `B4:T20`, containing one header row and 16 samples. It validates all field names, sample identifiers, data completeness, and numeric variable types. Missing values are not imputed.

Five composite variables are constructed in the original analysis order:

1. **Composite_1:** *Terrimonas* and unclassified Alphaproteobacteria;
2. **Composite_2:** *Claroideoglomus*, *Glomus*, unranked Claroideoglomeraceae, and unranked Glomeraceae;
3. **Composite_3:** soil pH, ammonium N, available P, total N, and total P;
4. **Composite_4:** urease, acid phosphatase, and phytase activities;
5. **Composite_5:** palmitoleic acid, linoleic acid, and trehalose relative values.

For each composite, aboveground biomass is first regressed on the component variables. The composite score is the sum of each component multiplied by its unstandardized regression slope, without adding the intercept. The script then fits the five component models of the piecewise SEM without adding or removing structural paths.

### Output

Results are printed to the R console rather than saved as figure files. Output includes:

- initial and composite-model summaries;
- standardized coefficients;
- piecewise SEM summary;
- AIC and BIC;
- chi-square and Fisher's C statistics; and
- R session information.

These numerical results support the piecewise SEM reported in the manuscript; the script does not itself render the final SEM diagram.

## Reproducibility and data-use notes

- Use the supplied workbooks without changing worksheet names, column order, or the populated ranges read by the scripts.
- The random forest workflow uses a fixed seed, but exact results can still depend on R and package versions.
- The network workflow overwrites `network_nodes.csv` and `network_edges.csv` if they already exist in the working directory.
- The MA and SEM scripts stop with informative errors when required inputs or fields are invalid.
- Relative signal intensities and relative quantitative values should not be interpreted as absolute concentrations.
- Correlation-based networks and observational associations do not, by themselves, establish ecological interactions or causality.

## Citation

If you use these data or code, please cite the associated manuscript and this repository. The manuscript citation is:

> Xu, Y., Cai, W., Chen, Z., & Duan, X. (2026). *A rhizosphere bacterial–AMF consortium improves phosphorus nutrition and growth of a native grass under extreme phosphorus limitation in karst rocky desertification*. Associated manuscript.

Journal, volume, pages, and DOI should be added to this section after publication. Please also include the repository URL and the specific release or commit used in any published analysis.

## License

No license was supplied with the source materials used to prepare this README. A `LICENSE` file should be added before public release to specify reuse permissions for the code and data. This README does not itself grant additional reuse rights.

## Contact

For questions about the data or code, please open an issue in this repository or contact the corresponding authors, **Zhe Chen** and **Xingwu Duan**, using the contact information provided with the associated manuscript.
