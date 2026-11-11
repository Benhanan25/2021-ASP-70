# ============================================================
# ICA 3: Antimicrobial Resistance Analysis in E. coli
# Data source: AMRnet public E. coli CSV snapshot
# This dataset contains curated sample-level AMR determinants
# and metadata; it is NOT raw FASTA genome sequence data.
# ============================================================

# 1. Install and load packages
packages <- c("tidyverse", "janitor", "readr", "ggplot2", "pheatmap")
new <- packages[!packages %in% rownames(installed.packages())]
if (length(new)) install.packages(new)

library(tidyverse)
library(janitor)
library(readr)
library(ggplot2)
library(pheatmap)

# 2. Create project folders
dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# 3. Download the real E. coli dataset from AMRnet
data_url <- paste0(
  "https://amrnet.s3.amazonaws.com/amrnet-latest/",
  "amrnetdb-Escherichia_coli.csv.gz"
)
data_file <- "data/amrnet_Escherichia_coli.csv.gz"

if (!file.exists(data_file)) {
  download.file(data_url, data_file, mode = "wb", method = "libcurl")
}

# 4. Import data and standardize column names
amr <- read_csv(data_file, show_col_types = FALSE) |>
  clean_names()

# Save a copy of the imported data as a normal CSV
write_csv(amr, "results/tables/amrnet_ecoli_imported.csv")

cat("Rows:", nrow(amr), "\n")
cat("Columns:", ncol(amr), "\n")
print(names(amr))

# 5. Basic data quality checks
missing_summary <- tibble(
  column = names(amr),
  missing_count = sapply(amr, function(x) sum(is.na(x))),
  missing_percent = round(
    sapply(amr, function(x) mean(is.na(x)) * 100), 2
  )
)
write_csv(missing_summary, "results/tables/missing_values_summary.csv")

# Remove exact duplicate rows only
amr_clean <- distinct(amr)
write_csv(amr_clean, "results/tables/amrnet_ecoli_clean.csv")

# 6. Identify useful columns robustly
# AMRnet exports can evolve; inspect names above if a column differs.
find_col <- function(data, candidates) {
  found <- intersect(candidates, names(data))
  if (!length(found)) return(NA_character_)
  found[1]
}

organism_col <- find_col(amr_clean, c("organism", "species"))
country_col  <- find_col(amr_clean, c("country", "iso_country_code"))
year_col     <- find_col(amr_clean, c("collection_year", "year"))
gene_col     <- find_col(amr_clean, c(
  "amr_determinant", "determinant", "gene", "gene_symbol",
  "amr_element_symbol", "resistance_gene"
))
drug_col     <- find_col(amr_clean, c(
  "antibiotic", "antibiotic_name", "drug", "antimicrobial"
))
class_col    <- find_col(amr_clean, c(
  "class", "antibiotic_class", "amr_class", "subclass"
))
phenotype_col <- find_col(amr_clean, c(
  "resistance_phenotype", "phenotype", "sir", "ast_phenotype"
))
sample_col <- find_col(amr_clean, c(
  "biosample_id", "sample_id", "isolate", "assembly_id"
))

cat("\nDetected columns:\n")
print(list(
  organism = organism_col, country = country_col, year = year_col,
  gene = gene_col, antibiotic = drug_col, class = class_col,
  phenotype = phenotype_col, sample = sample_col
))

# 7. Organism distribution
if (!is.na(organism_col)) {
  organism_summary <- amr_clean |>
    count(.data[[organism_col]], sort = TRUE, name = "records")
  write_csv(organism_summary, "results/tables/organism_summary.csv")

  p <- ggplot(organism_summary |> slice_head(n = 15),
              aes(x = reorder(.data[[organism_col]], records),
                  y = records)) +
    geom_col() + coord_flip() + theme_minimal() +
    labs(title = "AMR records by organism",
         x = "Organism", y = "Number of records")
  ggsave("results/plots/organism_distribution.png", p,
         width = 9, height = 6, dpi = 300)
}

# 8. Gene/determinant frequency
if (!is.na(gene_col)) {
  gene_summary <- amr_clean |>
    filter(!is.na(.data[[gene_col]]), .data[[gene_col]] != "") |>
    count(.data[[gene_col]], sort = TRUE, name = "records")
  write_csv(gene_summary, "results/tables/gene_frequency.csv")

  p <- ggplot(gene_summary |> slice_head(n = 20),
              aes(x = reorder(.data[[gene_col]], records),
                  y = records)) +
    geom_col() + coord_flip() + theme_minimal() +
    labs(title = "Most frequently reported AMR determinants",
         x = "AMR determinant", y = "Number of records")
  ggsave("results/plots/top_amr_determinants.png", p,
         width = 9, height = 7, dpi = 300)
} else {
  message("No obvious gene column detected. Inspect names(amr_clean).")
}

# 9. Antibiotic / class summaries
if (!is.na(drug_col)) {
  drug_summary <- amr_clean |>
    filter(!is.na(.data[[drug_col]]), .data[[drug_col]] != "") |>
    count(.data[[drug_col]], sort = TRUE, name = "records")
  write_csv(drug_summary, "results/tables/antibiotic_summary.csv")
}

if (!is.na(class_col)) {
  class_summary <- amr_clean |>
    filter(!is.na(.data[[class_col]]), .data[[class_col]] != "") |>
    count(.data[[class_col]], sort = TRUE, name = "records")
  write_csv(class_summary, "results/tables/antibiotic_class_summary.csv")

  p <- ggplot(class_summary,
              aes(x = reorder(.data[[class_col]], records),
                  y = records)) +
    geom_col() + coord_flip() + theme_minimal() +
    labs(title = "AMR records by antibiotic class",
         x = "Antibiotic class", y = "Number of records")
  ggsave("results/plots/antibiotic_class_distribution.png", p,
         width = 9, height = 7, dpi = 300)
}

# 10. Phenotype distribution (if phenotype data are present)
if (!is.na(phenotype_col)) {
  phenotype_summary <- amr_clean |>
    filter(!is.na(.data[[phenotype_col]])) |>
    count(.data[[phenotype_col]], sort = TRUE, name = "records")
  write_csv(phenotype_summary, "results/tables/phenotype_summary.csv")

  p <- ggplot(phenotype_summary,
              aes(x = reorder(.data[[phenotype_col]], records),
                  y = records)) +
    geom_col() + coord_flip() + theme_minimal() +
    labs(title = "Reported antimicrobial susceptibility phenotypes",
         x = "Phenotype", y = "Number of records")
  ggsave("results/plots/phenotype_distribution.png", p,
         width = 8, height = 5, dpi = 300)
}

# 11. Collection-year trend (if available)
if (!is.na(year_col)) {
  year_summary <- amr_clean |>
    mutate(year_value = suppressWarnings(as.integer(.data[[year_col]]))) |>
    filter(!is.na(year_value)) |>
    count(year_value, name = "records") |>
    arrange(year_value)
  write_csv(year_summary, "results/tables/year_summary.csv")

  if (nrow(year_summary) > 1) {
    p <- ggplot(year_summary, aes(x = year_value, y = records)) +
      geom_line() + geom_point() + theme_minimal() +
      labs(title = "AMR records by collection year",
           x = "Collection year", y = "Number of records")
    ggsave("results/plots/records_by_year.png", p,
           width = 9, height = 5, dpi = 300)
  }
}

# 12. Sample-level determinant matrix and heatmap
if (!is.na(sample_col) && !is.na(gene_col)) {
  gene_matrix_data <- amr_clean |>
    filter(!is.na(.data[[sample_col]]),
           !is.na(.data[[gene_col]]),
           .data[[gene_col]] != "") |>
    distinct(sample = .data[[sample_col]],
             determinant = .data[[gene_col]]) |>
    mutate(present = 1L) |>
    pivot_wider(names_from = determinant,
                values_from = present,
                values_fill = 0L)

  write_csv(gene_matrix_data, "results/tables/determinant_presence_matrix.csv")

  if (nrow(gene_matrix_data) >= 2 && ncol(gene_matrix_data) >= 3) {
    mat <- gene_matrix_data |>
      column_to_rownames("sample") |>
      as.matrix()

    # Keep up to 40 most common determinants for a readable plot
    if (ncol(mat) > 40) {
      keep_cols <- names(sort(colSums(mat), decreasing = TRUE))[1:40]
      mat <- mat[, keep_cols, drop = FALSE]
    }

    pheatmap(
      mat,
      cluster_rows = nrow(mat) > 1,
      cluster_cols = ncol(mat) > 1,
      main = "AMR determinant presence across samples",
      filename = "results/plots/amr_determinant_heatmap.png",
      width = 12, height = 9
    )
  }
}

# 13. Optional phenotype vs genotype analysis
# This is an association analysis, not proof of causation.
if (!is.na(phenotype_col) && !is.na(gene_col) &&
    !is.na(sample_col) && !is.na(drug_col)) {

  association_data <- amr_clean |>
    filter(!is.na(.data[[phenotype_col]]),
           !is.na(.data[[gene_col]]),
           !is.na(.data[[drug_col]])) |>
    transmute(
      sample = .data[[sample_col]],
      antibiotic = .data[[drug_col]],
      phenotype = .data[[phenotype_col]],
      determinant = .data[[gene_col]]
    ) |>
    distinct()

  if (nrow(association_data) > 0) {
    association_summary <- association_data |>
      count(antibiotic, phenotype, determinant, name = "records")
    write_csv(association_summary,
              "results/tables/phenotype_genotype_association.csv")
  }
}

# 14. Save reproducibility information
writeLines(capture.output(sessionInfo()),
           "results/session_info.txt")

cat("\nAnalysis script finished.\n")
cat("Review results/tables and results/plots.\n")
cat("Dataset source: AMRnet public E. coli snapshot.\n")
cat("Interpretation must account for sampling and database biases.\n")
