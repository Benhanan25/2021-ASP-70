# =====================================================================
# COMPLETE AMR BIOINFORMATICS ANALYSIS: ESCHERICHIA COLI
# ---------------------------------------------------------------------
# One integrated R script combining:
#   1. Package setup and AMRnet data download/import
#   2. Data cleaning and missing-value summary
#   3. Antibiotic-category marker frequencies and graph
#   4. Resistance-gene frequencies
#   5. Top-20-gene heatmaps (top-marker samples and random samples)
#   6. Country-wise AMR summary and final country graph
#   7. Chi-square test, Cramer's V and expected-frequency check
#   8. Standardised residual heatmap and Bonferroni correction
#
# Run this complete file from the project folder in RStudio.
# The script creates data/, results/tables/ and results/plots/.
#
# Important:
# - This analyses the AMRnet tabular E. coli export.
# - It does not download FASTA genomes or run AMRFinderPlus.
# - AMR marker presence is not the same as laboratory-confirmed AST.
# =====================================================================


# =====================================================================
# STEP 1: INSTALL AND LOAD REQUIRED PACKAGES
# =====================================================================

required_packages <- c(
  "tidyverse",
  "janitor",
  "pheatmap",
  "effectsize"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(tidyverse)
library(janitor)
library(pheatmap)
library(effectsize)


# =====================================================================
# STEP 2: CREATE PROJECT FOLDERS
# =====================================================================

dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/plots", showWarnings = FALSE, recursive = TRUE)


# =====================================================================
# STEP 3: DOWNLOAD AND IMPORT THE AMRnet E. COLI DATASET
# =====================================================================

data_url <- paste0(
  "https://amrnet.s3.amazonaws.com/amrnet-latest/",
  "amrnetdb-Escherichia_coli.csv.gz"
)

data_file <- "data/amrnet_Escherichia_coli.csv.gz"

if (!file.exists(data_file)) {
  message("Downloading the AMRnet E. coli dataset...")
  download.file(
    url = data_url,
    destfile = data_file,
    mode = "wb",
    method = "libcurl"
  )
}

message("Importing AMRnet dataset...")

amr <- readr::read_csv(
  data_file,
  show_col_types = FALSE,
  progress = TRUE
) %>%
  janitor::clean_names()

cat("\nDataset dimensions:\n")
print(dim(amr))

cat("\nDataset column names:\n")
print(names(amr))


# =====================================================================
# STEP 4: DATA QUALITY CHECKS AND CLEANING
# =====================================================================

# Calculate missing values for every column.
missing_summary <- tibble(
  column = names(amr),
  missing_count = vapply(amr, function(x) sum(is.na(x)), numeric(1)),
  missing_percent = vapply(
    amr,
    function(x) round(mean(is.na(x)) * 100, 2),
    numeric(1)
  )
)

write_csv(
  missing_summary,
  "results/tables/missing_values_summary.csv"
)

# Remove exact duplicate rows only.
amr_clean <- amr %>%
  distinct() %>%
  clean_names()

cat("\nRows before duplicate removal:", nrow(amr), "\n")
cat("Rows after duplicate removal:", nrow(amr_clean), "\n")
cat("Columns after cleaning:", ncol(amr_clean), "\n")

# Save the cleaned data as a compressed CSV to avoid creating
# another very large uncompressed duplicate file.
readr::write_csv(
  amr_clean,
  "results/tables/amrnet_ecoli_clean.csv"
)


# =====================================================================
# STEP 5: CHECK REQUIRED DATASET COLUMNS
# =====================================================================

required_columns <- c(
  "sample_id",
  "country_only",
  "aminoglycoside",
  "bl_inhibitor",
  "carbapenemase",
  "colistin",
  "esbl",
  "fosfomycin",
  "macrolide",
  "penicllin",
  "phenicol",
  "quinolone",
  "sulfonamide",
  "tetracycline",
  "trimethoprim",
  "other_classes"
)

missing_columns <- setdiff(required_columns, names(amr_clean))

if (length(missing_columns) > 0) {
  stop(
    paste(
      "The dataset is missing required columns:",
      paste(missing_columns, collapse = ", "),
      "\nCheck the current AMRnet export and update the column names."
    )
  )
}

if (all(is.na(amr_clean$sample_id))) {
  stop("The sample_id column contains no usable sample identifiers.")
}


# =====================================================================
# STEP 6: CREATE THE ANTIBIOTIC CATEGORY DATASET
# =====================================================================

amr_categories <- amr_clean %>%
  select(
    sample_id,
    aminoglycoside,
    bl_inhibitor,
    carbapenemase,
    colistin,
    esbl,
    fosfomycin,
    macrolide,
    penicllin,
    phenicol,
    quinolone,
    sulfonamide,
    tetracycline,
    trimethoprim,
    other_classes
  )

cat("\nAntibiotic category dataset dimensions:\n")
print(dim(amr_categories))

print(head(amr_categories))


# =====================================================================
# STEP 7: CALCULATE RESISTANCE MARKER FREQUENCIES
# =====================================================================

# Count distinct samples with a reported marker in each category.
# ND and "-" are treated as no reported marker.

category_columns <- setdiff(
  names(amr_categories),
  "sample_id"
)

total_unique_samples <- n_distinct(
  amr_categories$sample_id[!is.na(amr_categories$sample_id)]
)

category_frequency <- purrr::map_dfr(
  category_columns,
  function(col) {

    values <- as.character(amr_categories[[col]])

    tibble(
      Antibiotic_Category = col,
      Samples_With_Markers = n_distinct(
        amr_categories$sample_id[
          !is.na(amr_categories$sample_id) &
            !is.na(values) &
            trimws(values) != "" &
            !trimws(values) %in% c("ND", "-")
        ]
      )
    )
  }
) %>%
  mutate(
    Percentage = round(
      Samples_With_Markers / total_unique_samples * 100,
      2
    )
  ) %>%
  arrange(desc(Samples_With_Markers))

cat("\nAMR category frequency:\n")
print(category_frequency)

write_csv(
  category_frequency,
  "results/tables/category_frequency.csv"
)


# =====================================================================
# STEP 8: CREATE THE ANTIBIOTIC CATEGORY GRAPH
# =====================================================================

category_plot <- ggplot(
  category_frequency,
  aes(
    x = reorder(Antibiotic_Category, Samples_With_Markers),
    y = Samples_With_Markers
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "AMR Marker Distribution in E. coli",
    x = "Antibiotic Category",
    y = "Number of Samples"
  ) +
  theme_minimal()

print(category_plot)

ggsave(
  filename = "results/plots/amr_category_distribution.png",
  plot = category_plot,
  width = 10,
  height = 7,
  dpi = 300
)


# =====================================================================
# STEP 9: IDENTIFY THE MOST FREQUENT RESISTANCE GENES
# =====================================================================

# Use the 13 standard antibiotic categories and exclude other_classes,
# because other_classes has a different mixed format.

gene_data <- amr_categories %>%
  select(-other_classes) %>%
  pivot_longer(
    cols = -sample_id,
    names_to = "Antibiotic_Category",
    values_to = "Gene"
  ) %>%
  mutate(
    Gene = as.character(Gene),
    Gene = trimws(Gene)
  ) %>%
  filter(
    !is.na(Gene),
    Gene != "",
    !Gene %in% c("ND", "-")
  ) %>%
  separate_rows(Gene, sep = ",") %>%
  mutate(Gene = trimws(Gene)) %>%
  filter(Gene != "") %>%
  distinct(sample_id, Antibiotic_Category, Gene)

gene_frequency <- gene_data %>%
  count(Gene, sort = TRUE, name = "Sample_Frequency")

top_genes <- gene_frequency %>%
  slice_head(n = 20)

cat("\nTop 20 resistance-associated genes:\n")
print(top_genes)

write_csv(
  gene_frequency,
  "results/tables/gene_frequency.csv"
)

# This avoids counting the same marker multiple times for the same
# sample within the same antibiotic category.


# =====================================================================
# STEP 10: CREATE A HEATMAP USING THE TOP 20 GENES
# =====================================================================

top_20_genes <- gene_frequency %>%
  slice_max(
    order_by = Sample_Frequency,
    n = 20,
    with_ties = FALSE
  ) %>%
  pull(Gene)

# Remove missing and empty sample IDs.
clean_gene_data <- gene_data %>%
  filter(
    !is.na(sample_id),
    trimws(as.character(sample_id)) != "",
    Gene %in% top_20_genes
  ) %>%
  distinct(sample_id, Gene)

# Select 50 samples containing the largest number of top-20 markers.
selected_samples <- clean_gene_data %>%
  count(sample_id, name = "Total_Markers") %>%
  arrange(desc(Total_Markers)) %>%
  slice_head(n = 50) %>%
  pull(sample_id)

heatmap_data <- clean_gene_data %>%
  filter(sample_id %in% selected_samples) %>%
  mutate(Present = 1)

heatmap_matrix <- heatmap_data %>%
  select(sample_id, Gene, Present) %>%
  pivot_wider(
    names_from = Gene,
    values_from = Present,
    values_fill = 0
  ) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

storage.mode(heatmap_matrix) <- "numeric"

cat("\nTop-marker heatmap dimensions:\n")
print(dim(heatmap_matrix))

if (nrow(heatmap_matrix) >= 2 && ncol(heatmap_matrix) >= 2) {
  pheatmap(
    heatmap_matrix,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = "AMR Gene Distribution in E. coli",
    filename = "results/plots/amr_gene_heatmap.png",
    width = 12,
    height = 10
  )
} else {
  message("Not enough samples or genes to create the top-marker heatmap.")
}


# =====================================================================
# STEP 11: CREATE AN IMPROVED HEATMAP USING RANDOM SAMPLES
# =====================================================================

# Reproducible random selection.
set.seed(123)

valid_samples <- amr_categories %>%
  filter(
    !is.na(sample_id),
    trimws(as.character(sample_id)) != ""
  ) %>%
  distinct(sample_id) %>%
  pull(sample_id)

selected_samples_random <- sample(
  valid_samples,
  size = min(50, length(valid_samples))
)

random_heatmap_data <- gene_data %>%
  filter(
    sample_id %in% selected_samples_random,
    Gene %in% top_20_genes
  ) %>%
  distinct(sample_id, Gene) %>%
  mutate(Present = 1)

# Create all sample-gene combinations so that absent genes are zero.
random_matrix <- tidyr::expand_grid(
  sample_id = selected_samples_random,
  Gene = top_20_genes
) %>%
  left_join(
    random_heatmap_data,
    by = c("sample_id", "Gene")
  ) %>%
  mutate(Present = tidyr::replace_na(Present, 0)) %>%
  pivot_wider(
    names_from = Gene,
    values_from = Present
  ) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

storage.mode(random_matrix) <- "numeric"

if (nrow(random_matrix) >= 2 && ncol(random_matrix) >= 2) {
  pheatmap(
    random_matrix,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = paste(
      "AMR Marker Distribution in Randomly Selected",
      "E. coli Samples"
    ),
    filename = "results/plots/random_amr_gene_heatmap.png",
    width = 12,
    height = 10
  )
} else {
  message("Not enough data to create the random-sample heatmap.")
}


# =====================================================================
# STEP 12: COUNTRY-WISE AMR ANALYSIS
# =====================================================================

# Identify samples with reported AMR markers and country information.

sample_country <- amr_clean %>%
  filter(
    !is.na(sample_id),
    !is.na(country_only),
    trimws(as.character(country_only)) != ""
  ) %>%
  mutate(
    has_amr_marker = if_any(
      all_of(category_columns),
      ~ {
        value <- as.character(.)
        !is.na(value) &
          trimws(value) != "" &
          !trimws(value) %in% c("ND", "-")
      }
    )
  ) %>%
  group_by(sample_id, country_only) %>%
  summarise(
    has_amr_marker = any(has_amr_marker),
    .groups = "drop"
  )

# Calculate country-wise totals and marker percentages.
country_amr <- sample_country %>%
  group_by(country_only) %>%
  summarise(
    Total_Samples = n_distinct(sample_id),
    Samples_With_Markers = sum(has_amr_marker),
    .groups = "drop"
  ) %>%
  mutate(
    Percentage = round(
      Samples_With_Markers / Total_Samples * 100,
      2
    )
  ) %>%
  arrange(desc(Total_Samples))

cat("\nTop 15 countries by sample count:\n")
print(head(country_amr, 15))

write_csv(
  country_amr,
  "results/tables/country_amr_analysis.csv"
)


# =====================================================================
# STEP 13: CREATE THE FINAL COUNTRY-WISE AMR GRAPH
# =====================================================================

# Use the 10 countries with the largest sample counts.
country_plot_data <- country_amr %>%
  slice_max(
    order_by = Total_Samples,
    n = 10,
    with_ties = FALSE
  )

country_plot <- ggplot(
  country_plot_data,
  aes(
    x = reorder(country_only, Percentage),
    y = Percentage
  )
) +
  geom_col() +
  geom_text(
    aes(
      label = paste0(
        Percentage,
        "% (n=",
        Total_Samples,
        ")"
      )
    ),
    hjust = -0.1,
    size = 3
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 115),
    breaks = seq(0, 100, 20)
  ) +
  labs(
    title = "Country-wise Distribution of AMR Markers",
    x = "Country",
    y = "Samples with Reported Markers (%)"
  ) +
  theme_minimal()

print(country_plot)

ggsave(
  filename = "results/plots/country_amr_final_v2.png",
  plot = country_plot,
  width = 11,
  height = 7,
  dpi = 300
)


# =====================================================================
# STEP 14: CHI-SQUARE STATISTICAL ANALYSIS
# =====================================================================

# Select the 10 countries with the largest sample counts.
top_countries <- country_amr %>%
  slice_max(
    order_by = Total_Samples,
    n = 10,
    with_ties = FALSE
  ) %>%
  pull(country_only)

# Filter sample-level country data.
chi_data <- sample_country %>%
  filter(country_only %in% top_countries)

# Create the country-by-AMR-status contingency table.
chi_table <- table(
  chi_data$country_only,
  chi_data$has_amr_marker
)

cat("\nChi-square contingency table:\n")
print(chi_table)

# Perform Pearson's chi-square test.
chi_result <- chisq.test(chi_table)

cat("\nPearson chi-square test:\n")
print(chi_result)

# Save chi-square result.
capture.output(
  chi_result,
  file = "results/tables/chi_square_results.txt"
)


# =====================================================================
# STEP 15: CALCULATE CRAMER'S V
# =====================================================================

cramers_v_result <- effectsize::cramers_v(chi_table)

cat("\nCramer's V:\n")
print(cramers_v_result)

capture.output(
  cramers_v_result,
  file = "results/tables/cramers_v_results.txt"
)


# =====================================================================
# STEP 16: CHECK EXPECTED FREQUENCIES
# =====================================================================

cat("\nExpected frequencies:\n")
print(chi_result$expected)

number_expected_below_5 <- sum(
  chi_result$expected < 5
)

cat(
  "\nNumber of expected cell frequencies below 5:",
  number_expected_below_5,
  "\n"
)

write.csv(
  as.data.frame.matrix(chi_result$expected),
  "results/tables/chi_square_expected_frequencies.csv"
)

# If expected cell frequencies are too small, interpret the
# Pearson chi-square approximation cautiously. Consider an
# exact or Monte Carlo test where appropriate.


# =====================================================================
# STEP 17: CALCULATE STANDARDISED RESIDUALS
# =====================================================================

standardised_residuals <- round(
  chi_result$stdres,
  2
)

cat("\nStandardised residuals:\n")
print(standardised_residuals)

residual_data <- as.data.frame(
  as.table(chi_result$stdres)
)

colnames(residual_data) <- c(
  "Country",
  "AMR_Status",
  "Residual"
)


# =====================================================================
# STEP 18: CREATE THE CHI-SQUARE RESIDUAL PLOT
# =====================================================================

residual_plot <- ggplot(
  residual_data,
  aes(
    x = Country,
    y = AMR_Status,
    fill = Residual
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = round(Residual, 2)),
    size = 3
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  labs(
    title = "Standardised Residuals of Chi-square Analysis",
    x = "Country",
    y = "AMR Marker Status",
    fill = "Residual"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

print(residual_plot)

ggsave(
  filename = "results/plots/chi_square_residual_plot.png",
  plot = residual_plot,
  width = 11,
  height = 6,
  dpi = 300
)


# =====================================================================
# STEP 19: MULTIPLE-TESTING CORRECTION
# =====================================================================

# Calculate two-sided p-values for standardised residuals.
residual_data <- residual_data %>%
  mutate(
    p_value = 2 * pnorm(-abs(Residual)),
    adjusted_p = p.adjust(
      p_value,
      method = "bonferroni"
    )
  )

cat("\nResidual analysis with Bonferroni-adjusted p-values:\n")
print(residual_data)

write_csv(
  residual_data,
  "results/tables/chi_square_residual_analysis.csv"
)


# =====================================================================
# STEP 20: SAVE A FINAL ANALYSIS SUMMARY
# =====================================================================

summary_lines <- c(
  "AMRnet E. coli analysis completed.",
  paste("Rows in imported dataset:", nrow(amr)),
  paste("Rows after exact duplicate removal:", nrow(amr_clean)),
  paste("Unique sample IDs:", total_unique_samples),
  paste("Number of AMR categories:", nrow(category_frequency)),
  paste("Number of resistance gene labels:", nrow(gene_frequency)),
  paste("Number of country labels:", nrow(country_amr)),
  paste("Countries included in chi-square test:", length(top_countries)),
  paste("Chi-square statistic:", round(unname(chi_result$statistic), 4)),
  paste("Chi-square degrees of freedom:", unname(chi_result$parameter)),
  paste("Chi-square p-value:", format.pval(chi_result$p.value)),
  paste(
    "Number of expected cells below 5:",
    number_expected_below_5
  ),
  "",
  "Interpretation note:",
  "Country percentages describe records in this AMRnet dataset,",
  "not national population prevalence.",
  "Genotypic markers do not independently confirm phenotypic resistance."
)

writeLines(
  summary_lines,
  "results/tables/final_analysis_summary.txt"
)

cat("\n")
cat("============================================================\n")
cat("AMR BIOINFORMATICS ANALYSIS COMPLETED\n")
cat("============================================================\n")
cat("Check the results/tables/ and results/plots/ folders.\n")
