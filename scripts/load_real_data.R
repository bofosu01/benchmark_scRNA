# Load GSE158055 monocyte subset from the core batches and save matrix + metadata
# Uses scan() streaming to bypass R's 2^31-1 integer limit on the 2.3B non-zero mtx file

library(Matrix)
library(data.table)
library(readxl)

if (basename(getwd()) == "scripts") setwd("..")

counts_path        <- "dataset/GSE158055_covid19_counts.mtx.gz"
barcodes_path      <- "dataset/GSE158055_covid19_barcodes.tsv.gz"
features_path      <- "dataset/GSE158055_covid19_features.tsv.gz"
annotation_path    <- "dataset/GSE158055_cell_annotation.csv"
metadata_path      <- "dataset/GSE158055_sample_metadata.xlsx"
patient_group_path <- "dataset_characterization/samples_per_patient.csv"

counts_out <- "dataset/counts_mono_core.rds"
meta_out   <- "dataset/meta_mono_core.rds"

# Read barcodes and features via gzfile (fread fails on gzip in this system)
con_bc <- gzfile(barcodes_path, "r")
barcodes <- readLines(con_bc)
close(con_bc)

con_ft <- gzfile(features_path, "r")
features <- readLines(con_ft)
close(con_ft)

cat("Barcodes:", length(barcodes), " Features:", length(features), "\n")

# Read cell annotation: cell -> sample -> cell type
annotation <- fread(annotation_path)
setnames(annotation, c("cellName", "sampleID", "majorType"),
                   c("cell", "sample", "cell_type"))
annotation <- annotation[, .(cell, sample, cell_type)]

# Read sample metadata: sample -> donor -> batch
sample_meta <- as.data.table(read_excel(metadata_path))
setnames(sample_meta, c("Sample name", "Patients", "characteristics:  Datasets"),
                       c("sample", "donor", "batch"))
sample_meta <- sample_meta[, .(sample, donor, batch)]

# Read patient -> disease group
patient_group <- fread(patient_group_path)
setnames(patient_group, c("Patients", "disease_group"), c("donor", "group"))
patient_group <- patient_group[, .(donor, group)]

annotation <- annotation[sample_meta, on = "sample"]
annotation <- annotation[patient_group, on = "donor"]

core_batches <- c("Batch02", "Batch03", "Batch09", "Batch15")
annotation_mono <- annotation[cell_type == "Mono" & batch %in% core_batches]

keep_barcodes <- intersect(barcodes, annotation_mono$cell)
keep_col_idx  <- which(barcodes %in% keep_barcodes)
new_col_idx   <- match(barcodes[keep_col_idx], keep_barcodes)

cat("Keeping", length(keep_col_idx), "of", length(barcodes), "cells\n")

# Read mtx header
con <- gzfile(counts_path, "r")
invisible(readLines(con, n = 2))
dims <- as.numeric(strsplit(readLines(con, n = 1), " ")[[1]])
close(con)
n_genes <- dims[1]; n_cells <- dims[2]; nnz <- dims[3]
cat("mtx dims:", n_genes, "x", n_cells, ", nnz =", nnz, "\n")

col_map <- integer(n_cells)
col_map[keep_col_idx] <- new_col_idx

kept_i <- numeric(0); kept_j <- numeric(0); kept_x <- numeric(0)

chunk_size <- 1e6
con <- gzfile(counts_path, "r")
invisible(readLines(con, n = 3))

repeat {
  vals <- scan(con, what = numeric(), n = chunk_size * 3,
               quiet = TRUE, skipNul = TRUE)
  if (length(vals) == 0) break
  m <- matrix(vals, ncol = 3, byrow = TRUE)
  rm(vals)

  rows_keep <- m[, 2] %in% keep_col_idx
  if (any(rows_keep)) {
    sub <- m[rows_keep, , drop = FALSE]
    kept_i <- c(kept_i, sub[, 1])
    kept_j <- c(kept_j, col_map[sub[, 2]])
    kept_x <- c(kept_x, sub[, 3])
    rm(sub)
  }
  rm(m)
  gc()
  cat("  chunk done, total kept:", length(kept_x), "\n")
}
close(con)

cat("Kept", length(kept_x), "non-zero entries\n")

counts_mono <- sparseMatrix(
  i = as.integer(kept_i),
  j = as.integer(kept_j),
  x = kept_x,
  dims = c(n_genes, length(keep_barcodes))
)
rm(kept_i, kept_j, kept_x); gc()

rownames(counts_mono) <- features
colnames(counts_mono) <- keep_barcodes

meta_mono <- as.data.frame(annotation_mono)
meta_mono <- meta_mono[meta_mono$cell %in% colnames(counts_mono), ]
rownames(meta_mono) <- meta_mono$cell
meta_mono$sample_id <- meta_mono$sample

cat("Cells:", ncol(counts_mono), "\n")
cat("Donors:", length(unique(meta_mono$donor)), "\n")
print(table(meta_mono$group, useNA = "ifany"))
cps <- table(meta_mono$sample_id)
cat("Cells per sample: min =", min(cps), "median =", median(cps), "max =", max(cps), "\n")

saveRDS(counts_mono, counts_out)
saveRDS(meta_mono,   meta_out)
cat("Saved", counts_out, "and", meta_out, "\n")