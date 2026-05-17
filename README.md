# catcheR

**C**lone **A**nalysis of **T**ranscriptomic **C**onjugated **H**airpin **E**lement **R**eads

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
run_extraction(
  fastq      = "your_file.fastq",
  barcodes   = "rc_barcodes_genes.csv",
  output_dir = "results/",
  DIs        = 300
)

# Step 2 — single-library QC
run_plasmid_QC(
  results_dir      = "results/",
  DIs              = 300,
  transfect_clones = 100,
  transfect_cells  = 2000000,
  nucleofect_cells = 2000000
)

# Step 3 — compare two libraries
run_combined_QC(
  lib1_dir = "CATCHER1/results/",
  lib2_dir = "CATCHER2/results/",
  out_dir  = "combined_results/",
  DIs      = 300,
  transfect_clones = 100,
  transfect_cells  = 2000000,
  nucleofect_cells = 2000000
)

# Step 4 — publication plots
run_plasmid_plots("results/", DIs = 300)
```

## Requirements

- R >= 4.0.0
- Python 3.6+ (for `run_extraction()`)
- R packages: ggplot2, dplyr, tidyr, scales, forcats, patchwork, data.table

## Functions

| Function | Description |
|---|---|
| `run_extraction()` | Calls bundled Python script to extract UMI/BC/UCI from FASTQ |
| `run_plasmid_QC()` | Single-library QC + birthday-problem analysis |
| `run_combined_QC()` | Two-library comparison + go/no-go report |
| `run_plasmid_plots()` | Publication-ready distribution and QC plots |
