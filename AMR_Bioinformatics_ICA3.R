############################################################
# BIOINFORMATICS ICA 3
# Title:
# Bioinformatics-Based Identification and Analysis of
# Antimicrobial Resistance Genes in Bacterial Genomes Using R
#
# Workflow:
# NCBI genome download -> sequence QC -> AMRFinderPlus ->
# statistical analysis -> plots -> optional ML classification
#
# IMPORTANT:
# 1. Install NCBI Datasets CLI and AMRFinderPlus separately.
# 2. Fill data/metadata.csv with real NCBI assembly accessions.
# 3. For ML, phenotype must come from a real antimicrobial
#    susceptibility test (AST), not from AMR gene predictions.
############################################################

############################
# 0. PROJECT FOLDER SETUP
############################

project_dir <- getwd()

folders <- c(
  "data", "data/downloads", "results",
  "results/amrfinder", "results/plots", "results/tables",
  "logs"
)

for (folder in folders) {
  dir.create(file.path(project_dir, folder),
             recursive = TRUE, showWarnings = FALSE)
}

############################
# 1. INSTALL / LOAD PACKAGES
############################

required_packages <- c(
  "Biostrings", "ggplot2", "dplyr", "tidyr",
  "readr", "stringr", "tibble", "pheatmap",
  "caret", "randomForest", "e1071"
)

# Biostrings and pheatmap are available through Bioconductor.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c(
  "ggplot2", "dplyr", "tidyr", "readr", "stringr",
  "tibble", "caret", "randomForest", "e1071"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

bioc_packages <- c("Biostrings", "pheatmap")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(Biostrings)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(pheatmap)
})

set.seed(123)

############################
# 2. CREATE METADATA TEMPLATE
############################

metadata_file <- file.path("data", "metadata.csv")

if (!file.exists(metadata_file)) {
  template <- data.frame(
    sample_id = c("isolate_01", "isolate_02", "isolate_03"),
    assembly_accession = c(
      "REPLACE_WITH_NCBI_ASSEMBLY_ACCESSION_1",
      "REPLACE_WITH_NCBI_ASSEMBLY_ACCESSION_2",
      "REPLACE_WITH_NCBI_ASSEMBLY_ACCESSION_3"
    ),
    organism = c(
      "Escherichia coli",
      "Escherichia coli",
      "Klebsiella pneumoniae"
    ),
    antibiotic = c(NA, NA, NA),
    phenotype = c(NA, NA, NA)
  )

  write.csv(template, metadata_file, row.names = FALSE)

  stop(
    "A template was created at data/metadata.csv. ",
    "Replace the example values with real NCBI assembly accessions. ",
    "Then run this script again."
  )
}

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

required_columns <- c(
  "sample_id", "assembly_accession", "organism"
)

missing_columns <- setdiff(required_columns, names(metadata))
if (length(missing_columns) > 0) {
  stop("Missing metadata columns: ",
       paste(missing_columns, collapse = ", "))
}

if (anyDuplicated(metadata$sample_id)) {
  stop("Each sample_id must be unique.")
}

if (any(grepl("^REPLACE_|^$", metadata$assembly_accession))) {
  stop("Replace all placeholder assembly accessions in data/metadata.csv.")
}

if (!"antibiotic" %in% names(metadata)) {
  metadata$antibiotic <- NA_character_
}
if (!"phenotype" %in% names(metadata)) {
  metadata$phenotype <- NA_character_
}

metadata$sample_id <- as.character(metadata$sample_id)
metadata$assembly_accession <- as.character(metadata$assembly_accession)
metadata$organism <- as.character(metadata$organism)

write.csv(metadata, "results/tables/metadata_used.csv",
          row.names = FALSE)

############################
# 3. CHECK EXTERNAL PROGRAMS
############################

find_program <- function(program) {
  path <- Sys.which(program)
  if (identical(path, "")) {
    stop(
      "'", program, "' was not found in PATH. ",
      "Install it and restart R/RStudio."
    )
  }
  path
}

datasets_exe <- find_program("datasets")
amrfinder_exe <- find_program("amrfinder")

message("NCBI Datasets CLI: ", datasets_exe)
message("AMRFinderPlus: ", amrfinder_exe)

############################
# 4. DOWNLOAD GENOMES FROM NCBI
############################

download_genome <- function(accession, sample_id) {

  sample_dir <- file.path("data", "downloads", sample_id)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

  zip_file <- file.path(sample_dir, "genome_dataset.zip")

  if (!file.exists(zip_file)) {
    args <- c(
      "download", "genome", "accession", accession,
      "--include", "genome",
      "--filename", zip_file
    )

    status <- system2(
      datasets_exe,
      args = args,
      stdout = TRUE,
      stderr = TRUE
    )

    exit_status <- attr(status, "status")
    if (is.null(exit_status)) exit_status <- 0

    if (exit_status != 0 || !file.exists(zip_file)) {
      writeLines(status, file.path("logs",
                                   paste0(sample_id, "_download_error.txt")))
      stop("NCBI genome download failed for ", accession,
           ". Check logs and the accession.")
    }
  }

  extracted_dir <- file.path(sample_dir, "extracted")
  dir.create(extracted_dir, recursive = TRUE, showWarnings = FALSE)

  utils::unzip(zip_file, exdir = extracted_dir)

  fasta_files <- list.files(
    extracted_dir,
    pattern = "\\.(fna|fasta|fa)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(fasta_files) == 0) {
    stop("No genome FASTA file found for sample ", sample_id)
  }

  fasta_files[1]
}

message("Downloading genome assemblies from NCBI...")
metadata$fasta_file <- mapply(
  download_genome,
  metadata$assembly_accession,
  metadata$sample_id,
  USE.NAMES = FALSE
)

write.csv(metadata, "results/tables/metadata_with_fasta_paths.csv",
          row.names = FALSE)

############################
# 5. DNA SEQUENCE QUALITY CONTROL
############################

calculate_sequence_qc <- function(fasta_path, sample_id, organism) {

  dna <- readDNAStringSet(fasta_path)

  if (length(dna) == 0) {
    stop("Empty FASTA file: ", fasta_path)
  }

  seq_lengths <- width(dna)
  total_bases <- sum(seq_lengths)

  # GC content is calculated across all bases in the assembly.
  gc_count <- sum(letterFrequency(dna, letters = c("G", "C")))
  valid_base_count <- sum(letterFrequency(
    dna, letters = c("A", "C", "G", "T")
  ))

  gc_percent <- if (valid_base_count > 0) {
    100 * gc_count / valid_base_count
  } else {
    NA_real_
  }

  tibble(
    sample_id = sample_id,
    organism = organism,
    number_of_contigs = length(dna),
    total_sequence_length_bp = total_bases,
    mean_contig_length_bp = mean(seq_lengths),
    min_contig_length_bp = min(seq_lengths),
    max_contig_length_bp = max(seq_lengths),
    gc_percent = gc_percent
  )
}

qc_table <- bind_rows(lapply(seq_len(nrow(metadata)), function(i) {
  calculate_sequence_qc(
    metadata$fasta_file[i],
    metadata$sample_id[i],
    metadata$organism[i]
  )
}))

write.csv(qc_table, "results/tables/genome_quality_control.csv",
          row.names = FALSE)

p_gc <- ggplot(qc_table,
               aes(x = reorder(sample_id, gc_percent),
                   y = gc_percent,
                   fill = organism)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "GC Content of Bacterial Genome Assemblies",
    x = "Bacterial isolate",
    y = "GC content (%)",
    fill = "Organism"
  ) +
  theme_minimal()

ggsave("results/plots/gc_content.png", p_gc,
       width = 9, height = 6, dpi = 300)

p_contigs <- ggplot(qc_table,
                    aes(x = reorder(sample_id, number_of_contigs),
                        y = number_of_contigs,
                        fill = organism)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Number of Contigs per Genome Assembly",
    x = "Bacterial isolate",
    y = "Number of contigs",
    fill = "Organism"
  ) +
  theme_minimal()

ggsave("results/plots/contig_counts.png", p_contigs,
       width = 9, height = 6, dpi = 300)

############################
# 6. RUN NCBI AMRFINDERPLUS
############################

run_amrfinder <- function(fasta_path, sample_id, organism) {

  output_file <- file.path(
    "results", "amrfinder", paste0(sample_id, "_amrfinder.tsv")
  )

  args <- c(
    "--nucleotide", fasta_path,
    "--output", output_file
  )

  # Optional organism argument helps AMRFinderPlus apply
  # organism-specific rules when the organism name is accepted.
  if (!is.na(organism) && nzchar(organism)) {
    # AMRFinderPlus organism groups are commonly specified at genus level.
    organism_group <- str_extract(organism, "^[^ ]+")
    args <- c(args, "--organism", organism_group)
  }

  result <- system2(
    amrfinder_exe,
    args = args,
    stdout = TRUE,
    stderr = TRUE
  )

  exit_status <- attr(result, "status")
  if (is.null(exit_status)) exit_status <- 0

  if (exit_status != 0 || !file.exists(output_file)) {
    writeLines(
      result,
      file.path("logs", paste0(sample_id, "_amrfinder_error.txt"))
    )
    stop("AMRFinderPlus failed for sample ", sample_id,
         ". Check logs and AMRFinderPlus installation.")
  }

  output_file
}

message("Running AMRFinderPlus on genome assemblies...")

metadata$amrfinder_file <- mapply(
  run_amrfinder,
  metadata$fasta_file,
  metadata$sample_id,
  metadata$organism,
  USE.NAMES = FALSE
)

############################
# 7. IMPORT AND COMBINE AMRFINDER RESULTS
############################

read_amr_result <- function(tsv_file, sample_id, organism) {

  x <- read.delim(
    tsv_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )

  if (nrow(x) == 0) return(NULL)

  x$sample_id <- sample_id
  x$organism_from_metadata <- organism

  # AMRFinderPlus output normally includes "Element symbol",
  # "Class", and "Subclass" columns. Keep original column names.
  x
}

amr_results <- bind_rows(lapply(seq_len(nrow(metadata)), function(i) {
  read_amr_result(
    metadata$amrfinder_file[i],
    metadata$sample_id[i],
    metadata$organism[i]
  )
}))

if (nrow(amr_results) == 0) {
  warning(
    "AMRFinderPlus did not report any AMR-associated hits. ",
    "This does not prove the isolates are susceptible."
  )
}

write.csv(
  amr_results,
  "results/tables/all_amrfinder_results.csv",
  row.names = FALSE
)

############################
# 8. STANDARDIZE AMR RESULT COLUMNS
############################

# AMRFinderPlus versions may vary. Find columns by common names.
find_column <- function(data, possible_names) {
  matched <- possible_names[possible_names %in% names(data)]
  if (length(matched) == 0) return(NA_character_)
  matched[1]
}

symbol_col <- find_column(
  amr_results,
  c("Element symbol", "Gene symbol", "Element_symbol")
)

class_col <- find_column(
  amr_results,
  c("Class", "Element class")
)

subclass_col <- find_column(
  amr_results,
  c("Subclass", "Element subclass")
)

if (!is.na(symbol_col)) {
  amr_results$gene_symbol <- as.character(amr_results[[symbol_col]])
} else {
  amr_results$gene_symbol <- NA_character_
}

if (!is.na(class_col)) {
  amr_results$amr_class <- as.character(amr_results[[class_col]])
} else {
  amr_results$amr_class <- "Not available"
}

if (!is.na(subclass_col)) {
  amr_results$amr_subclass <- as.character(amr_results[[subclass_col]])
} else {
  amr_results$amr_subclass <- "Not available"
}

############################
# 9. DESCRIPTIVE STATISTICS
############################

gene_summary <- amr_results %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  count(gene_symbol, sort = TRUE, name = "number_of_hits")

write.csv(
  gene_summary,
  "results/tables/resistance_gene_frequency.csv",
  row.names = FALSE
)

class_summary <- amr_results %>%
  count(amr_class, sort = TRUE, name = "number_of_hits")

write.csv(
  class_summary,
  "results/tables/resistance_class_frequency.csv",
  row.names = FALSE
)

sample_summary <- amr_results %>%
  count(sample_id, name = "number_of_AMR_hits") %>%
  right_join(metadata %>% select(sample_id, organism), by = "sample_id") %>%
  mutate(number_of_AMR_hits = replace_na(number_of_AMR_hits, 0L))

write.csv(
  sample_summary,
  "results/tables/amr_hits_per_sample.csv",
  row.names = FALSE
)

############################
# 10. VISUALIZATION
############################

if (nrow(gene_summary) > 0) {

  top_genes <- gene_summary %>%
    slice_max(order_by = number_of_hits, n = 20, with_ties = FALSE)

  p_genes <- ggplot(
    top_genes,
    aes(x = reorder(gene_symbol, number_of_hits),
        y = number_of_hits)
  ) +
    geom_col() +
    coord_flip() +
    labs(
      title = "Most Frequently Detected AMR-Associated Gene Symbols",
      x = "Gene symbol",
      y = "Number of detections"
    ) +
    theme_minimal()

  ggsave("results/plots/top_resistance_genes.png",
         p_genes, width = 9, height = 7, dpi = 300)
}

if (nrow(class_summary) > 0) {

  p_classes <- ggplot(
    class_summary,
    aes(x = reorder(amr_class, number_of_hits),
        y = number_of_hits)
  ) +
    geom_col() +
    coord_flip() +
    labs(
      title = "AMR-Associated Hits by Antimicrobial Class",
      x = "Antimicrobial class",
      y = "Number of detections"
    ) +
    theme_minimal()

  ggsave("results/plots/amr_class_distribution.png",
         p_classes, width = 9, height = 7, dpi = 300)
}

p_sample <- ggplot(
  sample_summary,
  aes(x = reorder(sample_id, number_of_AMR_hits),
      y = number_of_AMR_hits,
      fill = organism)
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "AMR-Associated Hits per Bacterial Isolate",
    x = "Bacterial isolate",
    y = "Number of AMR-associated hits",
    fill = "Organism"
  ) +
  theme_minimal()

ggsave("results/plots/amr_hits_per_sample.png",
       p_sample, width = 9, height = 6, dpi = 300)

############################
# 11. GENE PRESENCE / ABSENCE MATRIX + HEATMAP
############################

if (nrow(amr_results) > 0 &&
    any(!is.na(amr_results$gene_symbol))) {

  gene_presence <- amr_results %>%
    filter(!is.na(gene_symbol), gene_symbol != "") %>%
    distinct(sample_id, gene_symbol) %>%
    mutate(present = 1L) %>%
    pivot_wider(
      names_from = gene_symbol,
      values_from = present,
      values_fill = 0L
    )

  gene_presence_df <- as.data.frame(gene_presence)
  rownames(gene_presence_df) <- gene_presence_df$sample_id
  gene_presence_df$sample_id <- NULL

  write.csv(
    cbind(sample_id = rownames(gene_presence_df), gene_presence_df),
    "results/tables/gene_presence_absence_matrix.csv",
    row.names = FALSE
  )

  if (nrow(gene_presence_df) >= 2 &&
      ncol(gene_presence_df) >= 2) {

    pheatmap(
      as.matrix(gene_presence_df),
      cluster_rows = nrow(gene_presence_df) > 1,
      cluster_cols = ncol(gene_presence_df) > 1,
      main = "AMR-Associated Gene Presence / Absence",
      filename = "results/plots/amr_gene_heatmap.png",
      width = 10,
      height = 8
    )
  }
}

############################
# 12. OPTIONAL STATISTICAL COMPARISON
############################

# Compare AMR hit counts across organisms.
# Use this only if there are enough independent isolates
# in each group. A simple descriptive summary is always saved.

organism_summary <- sample_summary %>%
  group_by(organism) %>%
  summarise(
    number_of_isolates = n(),
    mean_AMR_hits = mean(number_of_AMR_hits),
    median_AMR_hits = median(number_of_AMR_hits),
    sd_AMR_hits = sd(number_of_AMR_hits),
    .groups = "drop"
  )

write.csv(
  organism_summary,
  "results/tables/amr_summary_by_organism.csv",
  row.names = FALSE
)

if (n_distinct(sample_summary$organism) == 2 &&
    all(table(sample_summary$organism) >= 2)) {

  wilcox_result <- wilcox.test(
    number_of_AMR_hits ~ organism,
    data = sample_summary,
    exact = FALSE
  )

  capture.output(
    wilcox_result,
    file = "results/tables/wilcoxon_test.txt"
  )
}

############################
# 13. OPTIONAL MACHINE LEARNING
############################

# ML is attempted only when metadata contains:
# - antibiotic: one target antibiotic per isolate
# - phenotype: laboratory AST result, e.g. Resistant/Susceptible
#
# Do not use AMRFinderPlus calls as the ground-truth phenotype.
# The model uses detected gene presence/absence as predictors.

ml_possible <- all(c("antibiotic", "phenotype") %in% names(metadata)) &&
  any(!is.na(metadata$antibiotic)) &&
  any(!is.na(metadata$phenotype))

if (ml_possible) {

  ml_meta <- metadata %>%
    filter(
      !is.na(antibiotic),
      !is.na(phenotype),
      phenotype != ""
    )

  # Select one antibiotic with the most labelled isolates.
  target_antibiotic <- ml_meta %>%
    count(antibiotic, sort = TRUE) %>%
    slice(1) %>%
    pull(antibiotic)

  ml_meta <- ml_meta %>%
    filter(antibiotic == target_antibiotic) %>%
    mutate(
      phenotype = factor(
        str_to_title(as.character(phenotype))
      )
    )

  if (nrow(ml_meta) >= 10 &&
      n_distinct(ml_meta$phenotype) == 2 &&
      all(table(ml_meta$phenotype) >= 5)) {

    if (exists("gene_presence")) {

      ml_data <- ml_meta %>%
        select(sample_id, phenotype) %>%
        inner_join(gene_presence, by = "sample_id")

      if (nrow(ml_data) >= 10) {

        x <- ml_data %>%
          select(-sample_id, -phenotype)

        # Remove genes with no variation.
        keep <- vapply(x, function(z) length(unique(z)) > 1,
                       logical(1))
        x <- x[, keep, drop = FALSE]

        if (ncol(x) >= 1) {

          ml_data_final <- data.frame(
            phenotype = droplevels(ml_data$phenotype),
            x,
            check.names = FALSE
          )

          train_index <- caret::createDataPartition(
            ml_data_final$phenotype,
            p = 0.8,
            list = FALSE
          )

          train_data <- ml_data_final[train_index, , drop = FALSE]
          test_data <- ml_data_final[-train_index, , drop = FALSE]

          if (nrow(test_data) > 0 &&
              n_distinct(train_data$phenotype) == 2) {

            rf_model <- randomForest::randomForest(
              phenotype ~ .,
              data = train_data,
              ntree = 500,
              importance = TRUE
            )

            predictions <- predict(rf_model, newdata = test_data)

            cm <- caret::confusionMatrix(
              predictions,
              test_data$phenotype
            )

            capture.output(
              cm,
              file = "results/tables/random_forest_confusion_matrix.txt"
            )

            saveRDS(
              rf_model,
              "results/tables/random_forest_model.rds"
            )

            importance_df <- as.data.frame(
              randomForest::importance(rf_model)
            ) %>%
              rownames_to_column("gene")

            write.csv(
              importance_df,
              "results/tables/random_forest_gene_importance.csv",
              row.names = FALSE
            )

            message(
              "Optional ML completed for antibiotic: ",
              target_antibiotic
            )

          } else {
            message("ML skipped: training/test split lacks both classes.")
          }

        } else {
          message("ML skipped: no variable gene predictors.")
        }

      } else {
        message("ML skipped: insufficient matched labelled isolates.")
      }

    } else {
      message("ML skipped: no AMR gene presence/absence matrix.")
    }

  } else {
    message(
      "ML skipped: need at least 10 labelled isolates, ",
      "two phenotype classes, and at least 5 samples per class."
    )
  }

} else {
  message(
    "ML skipped. Add real laboratory AST labels to metadata.csv ",
    "to perform supervised classification."
  )
}

############################
# 14. SAVE SESSION INFORMATION
############################

writeLines(
  capture.output(sessionInfo()),
  "results/sessionInfo.txt"
)

############################
# 15. FINAL SUMMARY
############################

cat("\n========================================\n")
cat("BIOINFORMATICS ANALYSIS COMPLETED\n")
cat("========================================\n")
cat("Samples analysed:", nrow(metadata), "\n")
cat("AMRFinderPlus hits:", nrow(amr_results), "\n")
cat("Unique gene symbols:",
    dplyr::n_distinct(amr_results$gene_symbol, na.rm = TRUE), "\n")
cat("\nCheck the results/ folder for tables and plots.\n")
cat("Remember: gene detection is not the same as confirmed\n")
cat("antibiotic resistance from laboratory susceptibility testing.\n")
cat("========================================\n")
