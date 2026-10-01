# Uno Index with survC1 with pertubation
# ========================================================
# Forest Plot with Random-Effects Meta-Analysis
# ========================================================

# --- Preparation of figure ---

# Get sample and event counts from data table
library(here)
here()

source(here("environment", "requirements.R"))

# --- CONFIGURATION ---
external_cohorts_inputpath <- here("data", "processed_data", "geo")
output_dir <- here("results", "tables")

cohorts <- c("GSE17536", "GSE17537", "GSE29621", "GSE39582", "EMTAB12862", "ws4_spinal", "TCGA_test")

cvrt <- c(
  "vital_status", "SurvTime", "gender_male", "age_at_diagnosis",
  "Stage_II", "Stage_III", "Stage_IV", "MSI.Status_MSI.L.MSS"
)

# --- EXTRACTION FUNCTION FOR ALL COHORTS (INCLUDING TCGA_test) ---
extract_clinical_data <- function(cohort) {
  file_path <- here("data", "processed_data", "geo", paste0(cohort, "_clinical.rds"))
  clin_data <- readRDS(file_path)
  
  names(clin_data) <- make.names(names(clin_data))
  cols_to_keep <- intersect(cvrt, names(clin_data))
  
  clin_subset <- clin_data[, cols_to_keep, drop = FALSE]
  clin_subset$Cohort <- cohort
  
  return(clin_subset)
}

# --- MERGE COHORTS ---
merged_clinical_data <- bind_rows(lapply(cohorts, extract_clinical_data))

# --- DATA PREPROCESSING & 1825-DAY CENSORING ---
if("age_at_diagnosis" %in% names(merged_clinical_data)) {
  merged_clinical_data$age_at_diagnosis <- ifelse(
    merged_clinical_data$Cohort == "GSE29621",
    NA,
    merged_clinical_data$age_at_diagnosis
  )
}

# Convert binary gender_male (0/1) to Female/Male
if("gender_male" %in% names(merged_clinical_data)) {
  merged_clinical_data$gender <- factor(
    merged_clinical_data$gender_male,
    levels = c(0, 1),
    labels = c("Female", "Male")
  )
}

# Standardize vital status safely (handles both text and numeric/factor 0/1 codes)
merged_clinical_data <- merged_clinical_data %>%
  mutate(
    vital_status_char = as.character(vital_status),
    vital_status_num = case_when(
      vital_status_char %in% c("Alive", "0", 0) ~ 0,
      vital_status_char %in% c("Dead", "1", 1) ~ 1,
      TRUE ~ suppressWarnings(as.numeric(vital_status_char))
    ),
    vital_status_censored = case_when(
      SurvTime > 1825 ~ 0,
      TRUE ~ vital_status_num
    ),
    SurvTime = pmin(SurvTime, 1825, na.rm = TRUE)
  ) %>%
  mutate(
    vital_status = factor(vital_status_censored, levels = c(0, 1), labels = c("Alive", "Dead"))
  ) %>%
  dplyr::select(-vital_status_char, -vital_status_num, -vital_status_censored) # <-- Must have dplyr:: here

merged_clinical_data <- merged_clinical_data %>%
  mutate(stage = case_when(
    Stage_IV == 1 ~ "Stage 4",
    Stage_III == 1 ~ "Stage 3",
    Stage_II == 1 ~ "Stage 2",
    Stage_IV == 0 & Stage_III == 0 & Stage_II == 0 ~ "Stage 1",
    TRUE ~ NA_character_
  ))

# Map MSI status 0/1 to words (MSS / MSI)
msi_col <- intersect(c("MSI.Status_MSI.L.MSS", "MSI.Status_MSI-L.MSS"), names(merged_clinical_data))
if(length(msi_col) > 0) 
  
  
  # Extract table text from it
  
  # Compute samples and events per cohort dynamically from merged_clinical_data
  cohort_stats <- merged_clinical_data %>%
  group_by(Cohort) %>%
  summarise(
    samples = as.character(n()),
    events = as.character(sum(vital_status == "Dead", na.rm = TRUE)),
    .groups = "drop"
  )

# Define the exact order of cohorts and their corresponding display names
cohort_order <- c("GSE17536", "GSE17537", "GSE29621", "GSE39582", "TCGA_test", "EMTAB12862", "ws4_spinal")
display_names <- c("GSE17536", "GSE17537", "GSE29621", "GSE39582", "TCGA_test", "EMTAB12862", "SPINAL")

# Match and order stats according to the specified cohort order
stats_ordered <- cohort_stats[match(cohort_order, cohort_stats$Cohort), ]

# Construct tabletext dynamically (matching your hardcoded structure with NA for Summary)
tabletext <- cbind(
  c("Study ID", display_names, "Summary"),
  c("#Samples", stats_ordered$samples, NA),
  c("#Events", stats_ordered$events, NA)
)

# View tabletext matrix
print(tabletext)

################ Plot Forestplot ########################

# --- 1. Load Data ---
inputpath <- here("data", "processed_data", "uno_index_1825_days.rds")


data <- readRDS(inputpath) %>%
  dplyr::mutate(
    signature = ifelse(signature == "ours", "ColoTubIn", signature),
    Subset = ifelse(Subset == "ws4_spinal", "SPINAL", Subset)
  )

# --- 2. Data Preparation ---
genes_cvrts_pos <- grep(" + cvrts", data$Study, fixed = TRUE)
cvrts_pos        <- grep("^cvrts$", data$Study)

genes_cvrts <- data[genes_cvrts_pos, ]

cvrts <- data[cvrts_pos, ] %>%
  dplyr::group_by(Subset) %>%
  dplyr::slice(1) %>%
  dplyr::mutate(signature = "cvrts")

genes_cvrts <- rbind.data.frame(genes_cvrts, cvrts)
sig <- data[-c(genes_cvrts_pos, cvrts_pos), ]

# Define models and cohort order
all_model_names <- c("ColoTubIn","coloGuideEx","coloGuidePro","coloPrint",
                     "mda114","zhang","ICR","cvrts")

subset_order <- c("GSE17536","GSE17537","GSE29621","GSE39582",
                  "TCGA_test","EMTAB12862","SPINAL")

# --- 3. UPDATED HELPER: Random-Effects Meta-Analysis ---
create_model_df <- function(filtered_data, model_name) {
  
  # Filter and order the cohort data
  df_model <- filtered_data %>%
    dplyr::filter(signature == model_name) %>%
    dplyr::arrange(match(Subset, subset_order)) %>%
    dplyr::select(mean = C, lower = ci_l, upper = ci_u)
  
  # Calculate Standard Error (SE) from 95% Confidence Intervals
  # SE = (Upper - Lower) / 3.92
  df_model$se <- (df_model$upper - df_model$lower) / 3.92
  
  # Perform Random-Effects Meta-Analysis (REML method is most common)
  # This weights studies by their precision (inverse variance)
  res <- tryCatch({
    metafor::rma(yi = df_model$mean, sei = df_model$se, method = "REML")
  }, error = function(e) return(NULL))
  
  if (!is.null(res)) {
    pooled_mean <- as.numeric(res$beta)
    pooled_ci_l <- res$ci.lb
    pooled_ci_u <- res$ci.ub
  } else {
    # Fallback to simple mean if meta-analysis fails (e.g., too many NAs)
    pooled_mean <- mean(df_model$mean, na.rm = TRUE)
    pooled_ci_l <- mean(df_model$lower, na.rm = TRUE)
    pooled_ci_u <- mean(df_model$upper, na.rm = TRUE)
  }
  
  summary_row <- data.frame(
    mean  = pooled_mean,
    lower = pooled_ci_l,
    upper = pooled_ci_u
  )
  
  header_row <- data.frame(mean = NA, lower = NA, upper = NA)
  
  # Bind together: Header -> Data -> Pooled Summary
  final_df <- rbind(header_row, df_model[,c("mean","lower","upper")], summary_row)
  colnames(final_df) <- c("mean","lower","upper")
  final_df
}

# --- 4. Process Dataframes ---
filtered_cov <- genes_cvrts %>%
  dplyr::filter(signature %in% all_model_names, Subset %in% subset_order)

C_cov_list <- lapply(all_model_names, function(m) create_model_df(filtered_cov, m))
names(C_cov_list) <- all_model_names

filtered_genes <- sig %>%
  dplyr::filter(signature %in% all_model_names, Subset %in% subset_order)

G_list <- lapply(all_model_names, function(m) create_model_df(filtered_genes, m))
names(G_list) <- all_model_names

# --- 5. Formatting & Layout Settings --- # Already defined above
# tabletext <- cbind(
#   c("Study ID", "GSE17536","GSE17537","GSE29621","GSE39582",
#     "TCGA_test","EMTAB12862","SPINAL", "Summary"),
#   c("#Samples", "177","55","65","185", "33","782","106", NA),
#   c("#Events", "67","20","23","44", "11","245","36", NA)
# )
tabletext_right <- cbind(rep("", nrow(tabletext)), rep("", nrow(tabletext)))
w_equal <- rep(1, nrow(C_cov_list[[1]]))
box_size_vector <- c(0.12, rep(0.2, 7), 0.12)
vid_cols <- viridis(length(all_model_names))

output_dir <- here("results", "figures")

# --- 6. Generate Plot ---
pdf(file.path(output_dir,"Forestplot_Figure_3a.pdf"),  width = 14, height = 6)
grid.newpage()

# LEFT PANEL (With Covariates)
pushViewport(viewport(x = 0.30, y = 0.5, width = 0.65, height = 0.7, just = "center", name = "left_panel"))

forestplot::forestplot(
  tabletext,
  graph.pos = 4,
  hrzl_lines = list("2" = gpar(lwd = 1, col = "black"), "9" = gpar(lwd = 1, col = "black")),
  mean  = sapply(C_cov_list, `[[`, "mean"),
  lower = sapply(C_cov_list, `[[`, "lower"),
  upper = sapply(C_cov_list, `[[`, "upper"),
  col = fpColors(box = vid_cols,lines = vid_cols, summary = vid_cols),
  clip  = c(0.3, 1),
  is.summary  = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
  zero  = 0.5,
  xticks  = c(0.3, 0.5, 0.7, 0.9),
  graphwidth  = unit(4, "cm"),
  lineheight  = unit(1, "cm"),
  boxsize  = box_size_vector,
  weights  = w_equal,
  txt_gp  = fpTxtGp(ticks = gpar(cex = 0.7), xlab = gpar(cex = 1.3), lab = gpar(cex = 1)),
  new_page  = FALSE
)

grid.text("With covariates", x = unit(0.68, "npc"), y = unit(1.025, "npc") - unit(2.2, "lines"), gp = gpar(cex = 1.1, fontface = "bold"))
upViewport()

# RIGHT PANEL (Only Genes)
pushViewport(viewport(x = 0.545, y = 0.5, width = 0.65, height = 0.7, just = "center", name = "right_panel"))

G_models_subset <- all_model_names[all_model_names != "cvrts"]
vid_cols_genes <- vid_cols[1:length(G_models_subset)]

forestplot::forestplot(
  tabletext_right,
  graph.pos = 2,
  hrzl_lines = list("2" = gpar(lwd = 1), "9" = gpar(lwd = 1)),
  mean  = sapply(G_list[G_models_subset], `[[`, "mean"),
  lower = sapply(G_list[G_models_subset], `[[`, "lower"),
  upper = sapply(G_list[G_models_subset], `[[`, "upper"),
  col = fpColors(box = vid_cols_genes,lines = vid_cols,summary = vid_cols_genes),
  clip  = c(0.3, 1),
  is.summary  = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
  zero  = 0.5,
  xticks  = c(0.3, 0.5, 0.7, 0.9),
  graphwidth  = unit(4, "cm"),
  lineheight  = unit(1, "cm"),
  boxsize  = box_size_vector,
  weights  = w_equal,
  txt_gp  = fpTxtGp(ticks = gpar(cex = 0.7), xlab = gpar(cex = 1.3), lab = gpar(cex = 1)),
  new_page  = FALSE
)

grid.text("Only genes", x = unit(0.86, "npc"), y = unit(1.025, "npc") - unit(2.2, "lines"), gp = gpar(cex = 1.1, fontface = "bold"))
upViewport()

# Shared Legend and Axis
grid.text("Uno C-Index", x = 0.47, y = 0.11, gp = gpar(cex = 1.1))

legend_labels <- c("ColoTubIn","ColoGuideEx","ColoGuidePro","ColoPrint","MDA114","Zhang","ICR","Covariates")
legend_grob   <- legendGrob(labels = legend_labels, pch = 15, 
                            gp = gpar(col = vid_cols, cex = 1.0), 
                            nrow = 1, byrow = TRUE, hgap = unit(0.3, "cm"))

grid.draw(editGrob(legend_grob, vp = viewport(x = 0.41, y = 0.94, just = "center")))

dev.off()


###########################################################################################################
###########################################################################################################

############## MSI ########################

########### Uno Index ##############

# --- 1. Load Data ---
inputpath <- here("data", "processed_data", "uno_index_MSI_stratified_final_2026.rds")
full_data <- readRDS(inputpath) %>%
  mutate(Cohort = ifelse(tolower(Cohort) == "ws4_spinal", "SPINAL", Cohort))

# --- 2. Configuration & Simple Filtering ---
subset_order <- c("GSE17536","GSE17537","GSE29621","GSE39582","EMTAB12862","SPINAL")

# Just filter for the cohorts you want; leave MSI_Status alone!
full_data <- full_data %>%
  filter(Cohort %in% subset_order)

msi_levels <- c("Combined", "MSS", "MSI")
vid_cols <- viridis(3, end = 0.8, option = "viridis")

# ==========================================
# 3. Helper Functions
# ==========================================
# This ensures exactly 7 rows (1 header + 6 data)
create_msi_df <- function(all_data, model_name, msi_status) {
  df <- all_data %>%
    dplyr::filter(Model == model_name, MSI_Status == msi_status)
  
  df_plot <- data.frame(Cohort = subset_order) %>%
    left_join(df, by = "Cohort") %>%
    dplyr::select(mean = C, lower = ci_l, upper = ci_u)
  
  rbind(data.frame(mean=NA, lower=NA, upper=NA), df_plot)
}

# Helper to bind list of dataframes into a 7x3 matrix
prepare_matrix <- function(lst, col_name) {
  m <- do.call(cbind, lapply(lst, `[[`, col_name))
  return(as.matrix(m))
}

# ==========================================
# 4. Prepare Lists & Table Text
# ==========================================
LEFT_list  <- lapply(msi_levels, function(m) create_msi_df(full_data, "ours + cvrts", m))
RIGHT_list <- lapply(msi_levels, function(m) create_msi_df(full_data, "ours", m))
names(LEFT_list) <- names(RIGHT_list) <- msi_levels

# Calculate Sample Info labels (Fixed to only count MSI/MSS rows once)
sample_info_raw <- full_data %>%
  filter(MSI_Status %in% c("MSI", "MSS")) %>%
  group_by(Cohort) %>%
  summarise(
    Total = sum(N[Model == "ours"], na.rm = TRUE),
    MSS_N = max(N[MSI_Status == "MSS" & Model == "ours"], default = 0, na.rm = TRUE),
    MSI_N = max(N[MSI_Status == "MSI" & Model == "ours"], default = 0, na.rm = TRUE),
    .groups = 'drop'
  )

# Force labels to 6 rows even if a cohort has missing values
sample_info <- data.frame(Cohort = subset_order) %>%
  left_join(sample_info_raw, by = "Cohort") %>%
  mutate(across(where(is.numeric), ~replace_na(., 0)))

# tabletext has exactly 7 rows
tabletext <- cbind(
  c("Cohort ID", subset_order),
  c("#N All", sample_info$Total),
  c("#N MSS", sample_info$MSS_N),
  c("#N MSI", sample_info$MSI_N)
)

tabletext_right <- cbind(rep("", nrow(tabletext)))

# ==========================================
# 5. Plotting Configuration
# ==========================================
box_size_vector <- c(0.12, rep(0.2, 6))
is_summary_vec  <- c(TRUE, rep(FALSE, 6))

output_dir <- here("results", "figures")
pdf(file.path(output_dir,"Forestplot_Figure_3b.pdf"), width = 14, height = 6)

grid.newpage()

# --- LEFT PANEL: ours + cvrts ---
pushViewport(viewport(x = 0.28, y = 0.5, width = 0.65, height = 0.9, just = "center"))
forestplot::forestplot(
  tabletext,
  graph.pos = 5,
  hrzl_lines = list("2" = gpar(lwd = 1)), 
  mean  = prepare_matrix(LEFT_list, "mean"),
  lower = prepare_matrix(LEFT_list, "lower"),
  upper = prepare_matrix(LEFT_list, "upper"),
  col = fpColors(box = vid_cols, lines = vid_cols, summary = vid_cols),
  clip = c(0.3, 1.0), zero = 0.5,
  is.summary = is_summary_vec,
  xticks = c(0.3, 0.5, 0.7, 0.9),
  graphwidth = unit(4, "cm"),
  lineheight = unit(0.9, "cm"),
  boxsize = box_size_vector,
  txt_gp = fpTxtGp(ticks = gpar(cex = 0.7), lab = gpar(cex = 0.75)),
  new_page = FALSE
)
grid.text("With covariates", x = unit(0.69, "npc"), y = unit(0.725, "npc"), gp = gpar(fontface = "bold", cex = 1.0))
upViewport()

# --- RIGHT PANEL: ours only ---
pushViewport(viewport(x = 0.53, y = 0.5, width = 0.55, height = 0.9, just = "center"))
forestplot::forestplot(
  tabletext_right,
  graph.pos = 2,
  hrzl_lines = list("2" = gpar(lwd = 1)), 
  mean  = prepare_matrix(RIGHT_list, "mean"),
  lower = prepare_matrix(RIGHT_list, "lower"),
  upper = prepare_matrix(RIGHT_list, "upper"),
  col = fpColors(box = vid_cols, lines = vid_cols, summary = vid_cols),
  clip = c(0.3, 1.0), zero = 0.5,
  is.summary = is_summary_vec,
  xticks = c(0.3, 0.5, 0.7, 0.9),
  graphwidth = unit(4, "cm"),
  lineheight = unit(0.9, "cm"),
  boxsize = box_size_vector,
  txt_gp = fpTxtGp(ticks = gpar(cex = 0.7)),
  new_page = FALSE
)
grid.text("Only genes", x = unit(0.86, "npc"), y = unit(0.725, "npc"), gp = gpar(fontface = "bold", cex = 1.0))
upViewport()

# --- Shared Elements ---
legend_grob <- legendGrob(
  labels = c("All", "MSS", "MSI"), 
  pch = 15, 
  gp = gpar(col = vid_cols, cex = 1.0),
  nrow = 1
)
grid.draw(editGrob(legend_grob, vp = viewport(x = 0.35, y = 0.85)))
grid.text("Uno C-Index", x = 0.46, y = 0.21, gp = gpar(cex = 1.1))

dev.off()
