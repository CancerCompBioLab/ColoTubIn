library(here)
here()
source(here("environment", "requirements.R"))

# 1. Paths & Setup
# ----------------
input_path_counts <- here("data", "data", "pp","geo")
input_path_train  <- here("data", "data","cox","tcga", "signatures")
output_folder     <- here("results", "figures")

target_signature <- "ours_new"  # Internal ID
cohorts <- c("GSE17536", "GSE17537", "GSE29621", "GSE39582","EMTAB12862", "ws4_spinal")
cvrt    <- c("vital_status", "SurvTime", "gender_male", "age_at_diagnosis", "Stage_II", "Stage_III", "Stage_IV", "MSI.Status_MSI.L.MSS", "Molecular_Subtype_noCIN")

# 2. Function to Generate PCA
# -------------------------------------------------
generate_signature_pca <- function(sig_name) {
  
  # Load TCGA and External Cohorts
  train <- readRDS(file.path(input_path_train, sig_name, "train_all.rds"))
  all_cohort_genes <- list(tcga = setdiff(colnames(train), cvrt))
  
  ext_cohorts_counts <- list()
  for (cohort in cohorts) {
    counts <- readRDS(file.path(input_path_counts, paste0(cohort, "_counts.rds"))) %>% t()
    all_cohort_genes[[cohort]] <- colnames(counts)
    ext_cohorts_counts[[cohort]] <- as.data.frame(counts)
  }
  
  common_genes <- Reduce(intersect, all_cohort_genes)
  
  # Prepare for ComBat
  train_df <- train %>% dplyr::select(dplyr::all_of(common_genes)) %>% mutate(cohort = "tcga")
  mat_list <- list(train_df)
  for (cohort in cohorts) {
    ext_data <- ext_cohorts_counts[[cohort]] %>% dplyr::select(dplyr::all_of(common_genes)) %>% mutate(cohort = cohort)
    mat_list <- append(mat_list, list(ext_data))
  }
  
  full_mat <- do.call(rbind.data.frame, mat_list)
  meta <- data.frame(cohort = full_mat$cohort, row.names = rownames(full_mat))
  genes_raw <- full_mat %>% dplyr::select(-cohort)
  
  # Run ComBat
  modcombat <- model.matrix(~1, data = meta)
  combat_res <- ComBat(dat = t(genes_raw), batch = meta$cohort, mod = modcombat, 
                       ref.batch = "tcga", par.prior = TRUE) %>% t() %>% as.data.frame()
  
  # Perform PCA
  pca_raw <- prcomp(genes_raw, scale. = TRUE)
  pca_cor <- prcomp(combat_res, scale. = TRUE)
  
  # Calculate variance explained for the labels
  var_raw <- round(pca_raw$sdev^2 / sum(pca_raw$sdev^2) * 100, 1)
  var_cor <- round(pca_cor$sdev^2 / sum(pca_cor$sdev^2) * 100, 1)
  
  # Prepare Plotting Data
  pca_raw_df <- as.data.frame(pca_raw$x) %>% dplyr::select(PC1, PC2) %>%
    mutate(cohort = meta$cohort, Type = "Raw Data")
  
  pca_cor_df <- as.data.frame(pca_cor$x) %>% dplyr::select(PC1, PC2) %>%
    mutate(cohort = meta$cohort, Type = "ComBat Corrected")
  
  #Rename Spinal cohort + tcga first in legend
  plot_data <- rbind(pca_raw_df, pca_cor_df) %>%
    mutate(cohort = ifelse(cohort == "ws4_spinal", "SPINAL", cohort))%>%
    mutate(cohort = factor(cohort, levels = c("tcga", "EMTAB12862", "GSE17536", "GSE17537", "GSE29621", "GSE39582", "SPINAL")))
  
  plot_data$Type <- factor(plot_data$Type, levels = c("Raw Data", "ComBat Corrected"))
  
  # Create the label dataframe for geom_text
  text_labels <- data.frame(
    Type = factor(c("Raw Data", "ComBat Corrected"), levels = c("Raw Data", "ComBat Corrected")),
    label = c(
      sprintf("PC1: %.1f%%\nPC2: %.1f%%", var_raw[1], var_raw[2]),
      sprintf("PC1: %.1f%%\nPC2: %.1f%%", var_cor[1], var_cor[2])
    )
  )
  
  # 3. Plotting with EXACT Requested Styling
  # ----------------------------------------
  p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = cohort)) +
    geom_point(alpha = 0.7, size = 3) +
    scale_color_viridis_d(option = "viridis") +  # Original color scheme
    stat_ellipse(type = "t", level = 0.95, size = 0.5) +
    facet_wrap(~Type, scales = "free") +
    geom_text(
      data = text_labels,
      mapping = aes(label = label),
      x = -Inf, y = Inf, hjust = -0.05, vjust = 1.05, 
      size = 3.5, inherit.aes = FALSE, color = "black"
    ) +
    labs(
      title = NULL,  # Removed headline
      color = "Cohort",
      x = "PC1",
      y = "PC2"
    ) +
    theme_bw() +
    theme(
      legend.position = "top", 
      legend.direction = "horizontal",
      legend.box = "horizontal",
      strip.background = element_rect(fill = "white"), 
      strip.text = element_text(face = "bold")
    )
  
  return(p)
}


# 4. Save
# -------
final_plot <- generate_signature_pca(target_signature)

final_plot

# Save
# ggsave(
#   filename = file.path(output_folder, "PCA_plots_Suppl_Figure_1.pdf"),
#   plot = final_plot, 
#   width = 10, 
#   height = 6
# )
