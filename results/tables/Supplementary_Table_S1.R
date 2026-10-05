
########################################
# Supplementary Table S1
#######################################

library(here)
here()
source(here("environment", "requirements.R"))

# Load the saved matrix
res<- readRDS(here("data", "data", "signature_covariate_heatmap_matrix.rds"))

mat_log10 <- -log10(res)

# ---- Make a table -----

# --- Define Desired Order ---
desired_cols <- c("mda114", "coloGuideEx", "coloPrint", "ICR", "zhang", "coloGuidePro", "ColoTubIn")
desired_rows <- c("Consensus_clusters", "MSI", "Molecular_subtype", "braf", "apc", "stage", "tp53", "smad4", "kras")

# Subset/reorder matrix safely (keeping only columns/rows that exist in the matrix)
existing_cols <- intersect(desired_cols, colnames(mat_log10))
existing_rows <- intersect(desired_rows, rownames(mat_log10))

mat_log10 <- mat_log10[existing_rows, existing_cols, drop = FALSE]

tables_dir <- here("results", "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Convert matrix to a clean data frame with rownames as a column
df_table <- as.data.frame(mat_log10)
df_table$Signature <- rownames(df_table)
df_table <- df_table[, c("Signature", setdiff(names(df_table), "Signature"))] # Put Signature first

# Open PDF device, draw table, and close device
pdf(file.path(tables_dir, "Supplementary_Table_S1.pdf"), width = 9, height = 6)

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(1, 10), c("null", "null")))))


# Table grob
popViewport()
pushViewport(viewport(layout.pos.row = 2))
grid.draw(tableGrob(df_table, rows = NULL, theme = ttheme_minimal(base_size = 10)))

dev.off()

