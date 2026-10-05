#####################################
# Figure 3c 
####################################

library(here)
here()
source(here("environment", "requirements.R"))


# --- Make a Heatmap ---

# Load the saved matrix
res<- readRDS(here("data", "data","signatures_covariates", "signature_covariate_heatmap_matrix.rds"))

mat_log10 <- -log10(res)

col_fun <- colorRamp2(
  seq(0, max(mat_log10, na.rm = TRUE), length.out = 100), 
  viridis(100)
)

final_heatmap <- Heatmap(
  mat_log10,
  name = "-log10(p)",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_side = "left",
  column_names_rot = 45,
  col = col_fun 
)

final_heatmap

output_dir <- here("results", "figures")

# --- 5. Save the Heatmap ---
# Ensure the output directory exists
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# pdf(file.path(output_dir, "Heatmap_Figure_3c.pdf"), width = 10, height = 8)
# draw(final_heatmap)
# dev.off()
