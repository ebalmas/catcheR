# catcheR

## Clonality And Treatment Controlled sHrna Effect findeR
iPS2-seq updated pipeline for our new ligation based technology
updated version from https://github.com/alessandro-bertero/catcheR

End-to-end R package for analysing shRNA plasmid library sequencing data.

## Installation

```r
# Install dependencies
install.packages(c("ggplot2", "dplyr", "tidyr", "scales",
                   "forcats", "patchwork", "data.table"))

# Install devtools if needed
install.packages("devtools")

# Install catcheR from the local folder
devtools::install("path/to/catcheR")
```

## Quick start

```r
library(catcheR)

# Step 1 — extract clones from FASTQ (calls Python internally)
catcheR_step2QC_extraction(
  fastq      = "your_file.fastq",
  barcodes   = "rc_barcodes_genes.csv",
  output_dir = "results/",
  DIs        = 300
)

# Step 2 — single-library QC
catcheR_step2QC_plasmidQC(
  results_dir      = "results/",
  DIs              = 300,
  transfect_clones = 100,
  transfect_cells  = 2000000,
  nucleofect_cells = 2000000
)

# Step 3 — compare two libraries
catcheR_step2QC_combinedQC(
  lib1_dir = "CATCHER1/results/",
  lib2_dir = "CATCHER2/results/",
  out_dir  = "combined_results/",
  DIs      = 300,
  transfect_clones = 100,
  transfect_cells  = 2000000,
  nucleofect_cells = 2000000
)

# Step 4 — publication plots
catcheR_step2QC_plots("results/", DIs = 300)
```

## Requirements

- R >= 4.0.0
- Python 3.6+ (for `run_extraction()`)
- R packages: ggplot2, dplyr, tidyr, scales, forcats, patchwork, data.table

## Functions

| Function | Description |
|---|---|
| `catcheR_step2QC_extraction()` | Calls bundled Python script to extract UMI/BC/UCI from FASTQ |
| `rcatcheR_step2QC_plasmidQC()` | Single-library QC + birthday-problem analysis |
| `run_combinedQC()` | Two-library comparison + go/no-go report |
| `catcheR_step2QC_plots()` | Publication-ready distribution and QC plots |
