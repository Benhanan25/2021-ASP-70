# AMR Bioinformatics: *Escherichia coli* Antimicrobial Resistance Analysis

## 1. Project overview

This project analyses antimicrobial-resistance (AMR) information for *Escherichia coli* using the public AMRnet tabular data export. The main workflow uses R to inspect the dataset, clean duplicate records, summarize reported AMR categories and determinants, visualize frequently reported genes, examine country-level patterns, and perform an exploratory chi-square analysis.

The project package contains:

- The compressed AMRnet *E. coli* CSV data snapshot used for the saved analysis.
- R scripts for the tabular AMRnet analysis and a separate proposed genome-analysis workflow.
- A consolidated R script that combines the staged analysis into one end-to-end workflow.
- Result tables and figures.
- An editable R Markdown report and rendered HTML/PDF reports.

> **Important scope note:** The main analysis is based on AMRnet's tabular export. It is not a raw FASTA genome-analysis workflow. The separate `AMR_Bioinformatics_ICA3.R` script describes a genome/NCBI/AMRFinderPlus workflow, but corresponding execution outputs are not included in this package. Do not report that workflow as completed unless it has actually been run and its outputs have been checked.

## 2. Research aim

The practical aim is to describe the distribution of reported AMR determinants among *E. coli* records in the selected AMRnet data extract. The analysis focuses on:

1. Data structure and data quality.
2. Missing-value patterns and exact duplicate records.
3. Reported AMR marker frequencies across antibiotic categories.
4. Frequently reported AMR gene/determinant labels.
5. Presence/absence patterns for selected genes across samples.
6. Country-level sample counts and reported-marker proportions.
7. Exploratory association between country and reported AMR-marker status.

The outputs are descriptive and exploratory. They should not be interpreted as clinical recommendations or as national surveillance prevalence estimates.

## 3. Data source and data format

### Data source

The consolidated script downloads the AMRnet *E. coli* export from the following data endpoint when the local compressed data file is absent:

`https://amrnet.s3.amazonaws.com/amrnet-latest/amrnetdb-Escherichia_coli.csv.gz`

The project package also includes a compressed data snapshot at:

`data/amrnet_Escherichia_coli.csv.gz`

The endpoint uses a `latest` path. Therefore, a new download at a later date may contain updated records or schema changes and may not reproduce the saved results exactly. For reproducibility, use the packaged snapshot when available and record its download/access date in any new report.

### Data format

The data are supplied as a compressed CSV (`.csv.gz`) containing sample-level AMR determinants and metadata. The main scripts standardize column names with `janitor::clean_names()`.

The consolidated workflow expects the following fields after name cleaning:

- `sample_id`
- `country_only`
- `aminoglycoside`
- `bl_inhibitor`
- `carbapenemase`
- `colistin`
- `esbl`
- `fosfomycin`
- `macrolide`
- `penicllin`
- `phenicol`
- `quinolone`
- `sulfonamide`
- `tetracycline`
- `trimethoprim`
- `other_classes`

The script checks for these fields and stops with an explanatory message if one or more are missing. AMRnet may update its export structure, so inspect `names(amr)` and revise the column mapping if the current export differs.

## 4. Software and R packages

### Recommended environment

- R installed on the computer.
- RStudio Desktop is recommended for opening the project and running the scripts.
- Internet access is required only if the data file is not already present and must be downloaded.
- Sufficient free disk space is recommended because the source dataset and generated intermediate files may be large.

No specific R version is asserted by the source scripts. Use a current R release compatible with the packages below.

### Packages used by the consolidated analysis

The consolidated script installs missing packages and loads:

- `tidyverse` — data manipulation, reshaping, reading/writing CSV files, and plotting.
- `janitor` — standardizing column names.
- `pheatmap` — heatmap generation.
- `effectsize` — calculation of Cramér's V.

The separate `AMRnet_Ecoli_Analysis.R` script also uses functions available through `tidyverse`, `janitor`, `readr`, `ggplot2`, and `pheatmap`.

The separate `AMR_Bioinformatics_ICA3.R` genome workflow has additional dependencies and external-tool requirements. Review that script before attempting to run it; its optional machine-learning section refers to packages such as `caret` and `randomForest`, and the genome workflow requires the external sequence-analysis tools and input files described in the script.

## 5. Project folder structure

The cleaned project package is organized as follows:

```text
AMR_Bioinformatics/
├── AMR_Bioinformatics.Rproj
├── README.md
├── AMRnet_Ecoli_Analysis.R
├── AMR_Bioinformatics_Complete_Analysis.R
├── AMR_Bioinformatics_ICA3.R
├── AMR_Bioinformatics_Final_Report.Rmd
├── AMR_Bioinformatics_Final_Report.html
├── AMR_Bioinformatics_Final_Report.pdf
├── data/
│   └── amrnet_Escherichia_coli.csv.gz
└── results/
    ├── session_info.txt
    ├── plots/
    │   ├── category_frequency_from_table.png
    │   ├── chi_square_residual_plot.png
    │   ├── country_amr_final_v2.png
    │   ├── country_distribution.png
    │   ├── random_amr_gene_heatmap.png
    │   ├── top20_country_counts_from_table.png
    │   ├── top20_country_percent_from_table.png
    │   ├── top20_gene_frequency_from_table.png
    │   ├── amr_category_distribution.png
    │   └── amr_gene_heatmap.png
    └── tables/
        ├── category_frequency.csv
        ├── chi_square_residual_analysis.csv
        ├── chi_square_results.txt
        ├── cramers_v_results.txt
        ├── country_amr_analysis.csv
        ├── gene_frequency.csv
        └── missing_values_summary.csv
```

Some files may be generated or updated when a script is run. The folder listing above describes the main contents of the cleaned package; the consolidated script can create additional tables and figures.

## 6. Main scripts and their purposes

### A. `AMR_Bioinformatics_Complete_Analysis.R` — consolidated workflow

This is the single end-to-end script that combines the staged analysis code into one file. Run this script when you want to reproduce the main tabular analysis from data import through statistical summaries.

It performs the following steps:

1. Installs (if needed) and loads the required R packages.
2. Creates `data/`, `results/tables/`, and `results/plots/`.
3. Downloads the AMRnet compressed CSV if it is not already in `data/`.
4. Imports the CSV and standardizes column names.
5. Calculates missing-value counts and percentages for each column.
6. Removes exact duplicate rows.
7. Checks that required sample, country, and AMR-category columns exist.
8. Builds an antibiotic-category dataset.
9. Counts samples with reported markers in each category and calculates percentages.
10. Creates an antibiotic-category frequency plot.
11. Extracts resistance-gene/determinant labels, excludes placeholder/non-gene category fields as defined in the script, and summarizes the most frequent labels.
12. Builds a top-20-gene presence/absence heatmap using samples containing those markers.
13. Selects up to 50 samples using a fixed random seed (`set.seed(123)`) and creates a second heatmap.
14. Calculates country-level sample totals and the number/percentage with reported markers.
15. Creates the final country-level percentage plot for the ten countries with the largest sample counts.
16. Creates a country-by-marker-status contingency table for those ten countries.
17. Performs Pearson's chi-square test.
18. Calculates Cramér's V.
19. Saves expected cell frequencies for checking the chi-square approximation.
20. Calculates standardized residuals, creates a residual plot, applies Bonferroni adjustment to residual-based p-values, and writes a final text summary.

The script writes its outputs under `results/tables/` and `results/plots/`.

### B. `AMRnet_Ecoli_Analysis.R` — flexible AMRnet analysis script

This is the main analysis script originally included in the cleaned package. It imports the AMRnet table and uses flexible column detection for several fields (for example, organism, country, year, gene, antibiotic, class, phenotype, and sample identifiers). It can produce descriptive summaries and visualizations when the corresponding fields are present.

This script is useful for exploring a wider range of fields in the AMRnet export. It is separate from the consolidated script, so its outputs and assumptions may not be identical.

### C. `AMR_Bioinformatics_ICA3.R` — separate genome workflow

This script describes a different workflow based on genome/sample inputs and AMRFinderPlus-related outputs. It includes steps for metadata and AMR result processing, sample-level summaries, gene presence/absence, optional statistical comparisons, and an optional random-forest workflow when suitable laboratory phenotype labels are available.

**This workflow is not the source of the saved AMRnet tabular analysis results in this package.** The corresponding raw genome files, metadata, AMRFinderPlus result files, and matching run outputs are not included in the cleaned package. Treat it as a separate workflow/proposed analysis until it is run and validated independently.

## 7. How to run the consolidated analysis

### Step 1 — Extract the package

Extract `AMR_Bioinformatics_Cleaned_Project_and_Final_Report.zip` to a folder on your computer.

### Step 2 — Open the RStudio project

Open `AMR_Bioinformatics.Rproj` in RStudio. This helps keep relative file paths aligned with the project folder.

### Step 3 — Check the data file

Confirm that this file exists:

```text
data/amrnet_Escherichia_coli.csv.gz
```

If it exists, the consolidated script uses the local file. If it does not exist, the script attempts to download the latest AMRnet export from the endpoint described above.

### Step 4 — Open the consolidated script

Open:

```text
AMR_Bioinformatics_Complete_Analysis.R
```

### Step 5 — Run the script

Run the complete script from top to bottom in RStudio, or use the R console from the project directory:

```r
source("AMR_Bioinformatics_Complete_Analysis.R")
```

Do not run it from an unrelated working directory. The script uses relative paths such as `data/` and `results/`.

### Step 6 — Review generated results

After execution, inspect:

- `results/tables/` for CSV and TXT outputs.
- `results/plots/` for PNG visualizations.
- The R console for dataset dimensions, progress messages, and statistical output.

The script may create a large cleaned CSV table. Check available disk space before running it on a machine with limited storage.

## 8. Outputs and what they mean

### Data-quality output

**`missing_values_summary.csv`**

Contains one row per input column, with:

- `column` — column name.
- `missing_count` — number of missing values.
- `missing_percent` — percentage of rows missing that value.

Missingness is field-specific. A high missing percentage for a metadata field does not automatically mean the AMR marker data are unusable; it means that field has limited coverage in the extract.

### Category-level output

**`category_frequency.csv`**

Summarizes the number of unique samples with a reported marker in each antibiotic category and the percentage relative to the unique sample denominator used by the script.

**`amr_category_distribution.png`**

Visual representation of category-level marker frequencies.

The denominator and marker rules matter. The script treats missing values, blank strings, `ND`, and `-` as no reported marker for the category calculation. This is a data-processing rule, not proof that the isolate is biologically susceptible.

### Gene/determinant output

**`gene_frequency.csv`**

Contains resistance-gene/determinant labels and their reported frequencies. A determinant label may represent a gene or a mutation-associated marker, depending on the AMRnet field content.

**`amr_gene_heatmap.png`**

Shows the selected top genes across samples in a presence/absence matrix.

**`random_amr_gene_heatmap.png`**

Shows the selected top genes across up to 50 randomly selected samples. The fixed random seed makes the selection reproducible for the same input data and sample ordering.

Heatmaps visualize reported marker calls. A zero in the matrix means that the marker was not present in the selected AMRnet record representation; it should not automatically be interpreted as a confirmed susceptible laboratory result.

### Country-level output

**`country_amr_analysis.csv`**

Contains country-level sample counts, samples with reported AMR markers, and the percentage of samples with reported markers, as calculated from the available records.

**`country_amr_final_v2.png`**

Displays the reported-marker percentage for the ten countries with the largest sample counts.

Country counts and percentages describe the records available in this AMRnet extract. They are affected by data submission, sampling, surveillance coverage, and metadata completeness. They are **not** estimates of national population prevalence and should not be used to rank countries by the true burden of AMR.

### Statistical output

**`chi_square_results.txt`**

Contains the Pearson chi-square test result for the country-by-marker-status contingency table among the ten countries selected by sample count.

**`cramers_v_results.txt`**

Contains Cramér's V, an effect-size measure for the association in the contingency table.

**`chi_square_expected_frequencies.csv`**

Contains the expected cell frequencies used to assess the chi-square approximation.

**`chi_square_residual_analysis.csv`**

Contains standardized residuals and Bonferroni-adjusted p-values for the country/status cells.

**`chi_square_residual_plot.png`**

Visualizes the standardized residuals. Positive residuals indicate observed counts above the expected counts under independence; negative residuals indicate observed counts below expected counts.

A small p-value does not by itself establish a meaningful or causal relationship. Interpret the test alongside the effect size, expected-frequency diagnostics, sample composition, and data limitations. If expected counts are too small, the Pearson chi-square approximation may be unreliable; consider an appropriate exact or Monte Carlo procedure.

## 9. Interpretation of the saved analysis

The saved output tables indicate that the AMRnet extract contains many reported AMR determinant records and substantial variation in category and gene frequencies. The saved country analysis also shows that the number of records differs considerably between country labels.

Use the actual CSV/TXT result tables and the final report when quoting numerical findings. Do not treat record counts as independent patients unless the data documentation establishes that interpretation. Multiple records may relate to the same isolate or sample depending on the structure of the source export.

The statistical outputs are exploratory associations within the selected extract. They do not demonstrate that country causes AMR-marker presence.

## 10. Important scientific and methodological limitations

1. **Genotype is not phenotype.** Detection of an AMR gene or mutation does not, by itself, prove laboratory-confirmed resistance to a specific antibiotic.
2. **AMRnet extract is not raw sequence data.** The consolidated analysis works with tabular AMRnet records and marker fields. It does not perform read mapping, genome assembly, or AMRFinderPlus analysis.
3. **Country representation is uneven.** Some countries contribute far more records than others. Percentages from this extract are not national prevalence estimates.
4. **Missing metadata can affect subgroup analyses.** Country, city, disease, and other fields may have substantial missingness. Report denominators clearly.
5. **Exact duplicate removal is limited.** The consolidated script removes exact duplicate rows; it does not infer or merge biological duplicates with different metadata.
6. **AMRnet may change over time.** The `latest` download endpoint can change. Use the packaged snapshot for the saved analysis and note the date of any later download.
7. **Chi-square assumptions must be checked.** Review expected cell frequencies and consider the sampling structure. A very large dataset can produce statistically significant results for small differences.
8. **Multiple testing matters.** The residual analysis applies a Bonferroni correction to the residual-based tests. Report adjusted rather than unadjusted results when discussing those cell-wise comparisons.
9. **No causal conclusion is supported.** These are descriptive and association analyses.
10. **Genome workflow is separate.** Do not combine the ICA 3 genome workflow with the AMRnet results unless both workflows have been run on clearly documented, compatible data.

## 11. Reproducing the report

The package includes:

- `AMR_Bioinformatics_Final_Report.Rmd` — editable R Markdown source.
- `AMR_Bioinformatics_Final_Report.html` — rendered HTML report.
- `AMR_Bioinformatics_Final_Report.pdf` — rendered PDF report.

The R Markdown report is designed to use the saved analysis outputs rather than automatically rerun every step of the consolidated analysis. To refresh the report after changing the data or rerunning analyses, first regenerate the required tables and plots, then open the `.Rmd` file in RStudio and render it. Check that all referenced output paths exist.

## 12. Troubleshooting

### Error: required columns are missing

The AMRnet export may have changed. Inspect:

```r
names(amr)
```

Then compare the returned names with the required column list in Step 5 of `AMR_Bioinformatics_Complete_Analysis.R`. Update the field mapping only after confirming the meaning of the replacement fields.

### Error: cannot download the dataset

Check the internet connection and try accessing the AMRnet data endpoint in a browser. If you already have the compressed CSV snapshot, place it at:

```text
data/amrnet_Escherichia_coli.csv.gz
```

The script will use the local copy and skip downloading.

### Error: package is not installed

The consolidated script attempts to install its required packages. If installation fails, install the package manually in R:

```r
install.packages(c("tidyverse", "janitor", "pheatmap", "effectsize"))
```

### Output folders are empty

Confirm that the script was run from the project root (the folder containing the `.Rproj` file). Run:

```r
getwd()
```

If necessary, open the `.Rproj` file and rerun the script.

### Heatmap or plot is not generated

Some plotting steps require enough valid samples or marker columns. The script includes checks for insufficient data. Review the console messages and inspect the relevant input table.

## 13. Suggested citation and acknowledgement

When using this project in coursework or a report:

- Cite AMRnet as the source of the underlying AMR data according to the data provider's recommended citation and terms.
- State the date on which the data snapshot was accessed or downloaded.
- Cite the R packages and external tools actually used.
- Distinguish the AMRnet tabular analysis from the separate genome/AMRFinderPlus workflow.
- Describe the analysis as exploratory unless the study design and data provenance support stronger conclusions.

Do not invent a formal AMRnet citation. Obtain the current recommended citation from the data provider before final submission.

## 14. Project status

This package contains saved tabular AMRnet analysis outputs and a rendered report. The consolidated script is supplied as a reproducible, integrated workflow. The separate genome/AMRFinderPlus script is included as a distinct workflow, but its execution is not demonstrated by the saved outputs in this package.

---

**End of README**
