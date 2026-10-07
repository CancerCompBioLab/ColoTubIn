library(here)
here()
source(here("environment", "requirements.R"))


# Add TCGA complete gene matrix instead of limiting to signature genes for external validation analysis 
########################################

# Read previous subset containing only signature genes (e.g., from ours)
train_all <- readRDS(here("data", "data","cox", "tcga", "signatures", "ours", "train_all.rds"))

# Add all genes from vst transformed FPKM
counts_norm_coad_patients <- readRDS(here("data","extdata", "tcga-coad", "counts_norm_coad_patients.rds"))

# Keep only the first 9 columns
train_all_reduced <- train_all[, 1:9]

# Step 1. Identify common rows
common_rows <- intersect(
  rownames(train_all_reduced),
  rownames(counts_norm_coad_patients)
)

# Step 2. Subset both dataframes to those rows
train_sub <- train_all_reduced[common_rows, , drop = FALSE]
counts_sub <- counts_norm_coad_patients[common_rows, , drop = FALSE]

# Step 3. Ensure identical row order
counts_sub <- counts_sub[rownames(train_sub), , drop = FALSE]

# Step 4. Merge columns side-by-side
merged_df <- cbind(train_sub, counts_sub)


# Define the base directory and all signature subdirectories
base_dir <- here("data","data", "cox", "tcga_new")
signatures <- c("", "coloGuideEx", "coloGuidePro", "coloPrint", "mda114", "ours", "zhang", "ICR")

# Loop through, create directories if needed, and save
for (sig in signatures) {
  # If sig is empty, it saves to base_dir; otherwise it adds "signatures/sig_name"
  dir_path <- if (sig == "") base_dir else file.path(base_dir, "signatures", sig)
  
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  saveRDS(merged_df, file = file.path(dir_path, "train_all.rds"))
}
