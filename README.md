# SEBULA

<!-- badges: start -->
<!-- badges: end -->

SEBULA is an R package for detecting doublets in single-cell ATAC-seq data using fragment overlap patterns and central matching statistics.

---

## 📦 Pipeline Overview

The workflow consists of four distinct, sequential steps:

1. **Barcode Generation**
2. **Fragment Overlap Calculation**
3. **Overlap Processing and Doublet Detection**
4. **Doublet Classification**

Each step generates outputs necessary for subsequent analysis, culminating in the precise classification of each cell barcode as either a "Singlet" or a "Doublet."

---

### 1. Barcode Generation

Extract valid cell barcodes from raw and filtered feature barcode matrices provided in HDF5 format.

**Input Files:**
- `filtered_feature_bc_matrix.h5`
- `raw_feature_bc_matrix.h5`

**Output File:**
- `singlecell.csv` (A CSV file listing each barcode along with a binary flag indicating whether the barcode passed quality control.)

---

### 2. Fragment Overlap Calculation

Compute overlaps between ATAC-seq fragments and specified genomic regions (e.g., autosomes).

**Input Files:**
- `atac_fragments.tsv.gz` (ATAC-seq fragment file)
- `singlecell.csv` (from previous step)
- `human_autosomes.txt` (bundled with the package)

**Output Files:**
- `Overlaps.txt`
- `OverlapSummary.txt`

---

### 3. Overlap Processing and Doublet Detection

Process overlap data to filter repetitive genomic regions (optional), summarize overlaps per cell, and generate the final count table (`colsum.csv`).

**Input Files:**
- `Overlaps.txt` and `OverlapSummary.txt` (from previous step)
- `blacklist_repeats_segdups_rmsk_hg38.bed` (optional repetitive region filter)

**Output File:**
- `colsum.csv` (summarized overlap counts per barcode, columns: `barcode`, `obs`)

---


### 4. Doublet Classification

Utilize R-based statistical modeling (`doublet_cm` function) to classify cells as doublets or singlets. This method applies Box-Cox transformations, estimates null distributions via central matching, calculates False Discovery Rates (FDR), and classifies cells based on a chosen threshold.

**Input File:**
- `colsum.csv` (output from previous step)

**R Script:**
```r
library(SEBULA)
library(readr)

filtered_h5 <- "/hits/Johann_Novaseq/fastq/8193-JF/10x_analysis_8193-JF/Sample_8193-JF-3/filtered_feature_bc_matrix.h5"
raw_h5 <- "/hits/Johann_Novaseq/fastq/8193-JF/10x_analysis_8193-JF/Sample_8193-JF-3/raw_feature_bc_matrix.h5"
atac_fragments <- "/hits/Johann_Novaseq/fastq/8193-JF/10x_analysis_8193-JF/Sample_8193-JF-3/atac_fragments.tsv.gz"
output_directory <- "/home/yuntian/Doublet_Detection/New_approach/atac_doublet_sharing_DIG_GSE200417/package_test/8193-JF-3"
rfilter_file = "/home/yuntian/Doublet_Detection/New_approach/atac_doublet_sharing_DIG_GSE200417/package_test/blacklist_repeats_segdups_rmsk_hg38.bed"

colsum_path <- run_preprocessing_pipeline(
  filtered_h5 = filtered_h5,
  raw_h5 = raw_h5,
  atac_fragments = atac_fragments,
  output_directory = output_directory,
  rfilter_file = rfilter_file
)

res <- run_detection(colsum_path)
```

---

## Output and Interpretation

Upon completion of the pipeline, you will receive:
- A detailed CSV file (`res[["result"]]`) classifying each barcode as "Doublet" or "Singlet" along with statistical metrics.
- Console outputs summarizing the total number of doublets and their proportion within your dataset.

---

## Dependencies

- **Python dependencies:**
  - `h5py`, `numpy`, `pandas`

- **R dependencies:**
  - `dplyr`, `readr`, `ggplot2`, `MASS`, `splines`

Ensure these dependencies are installed in your Python and R environments.

---

## Installation and Setup

Install the R package directly from GitHub:

```R
install.packages("devtools")
devtools::install_github("Yuntian0716/SEBULA")
```

---


## License

This project is licensed under the MIT License.



