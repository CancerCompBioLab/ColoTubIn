#### MSI 

library(here)
here()
source(here("environment", "requirements.R"))

###################### With Uno Index ######################

# ==========================================
# 1. Configuration & Paths
# ==========================================
inf_iter <- 1000  # Set to 10 for a quick test, 1000 for final results
tau_fixed <- 1825  # 5-year truncation (in days)

# Inputpaths
external_cohorts_inputpath <- here::here("data", "data", "pp", "geo")
cox_models_inputpath <- here::here("data", "data", "cox", "tcga_new", "signatures")
signatures_annot_inputpath <- here::here("data", "extdata", "signatures", "signatures_annot_w_EMTAB_ws4_spinal_tcga_test_ICR_signature.rds")

# Output Paths
output_dir <- here::here("data", "data", "cox", "external")
rds_output <- file.path(output_dir, "uno_index_MSI_stratified_final_2026.rds")

cvrt <- c(
  "vital_status", "SurvTime", "gender_male", "age_at_diagnosis",
  "Stage_II", "Stage_III", "Stage_IV", "MSI.Status_MSI.L.MSS", "Molecular_Subtype_noCIN"
)

selected_signature <- "ours"
cohorts <- c("GSE17536", "GSE17537", "GSE29621", "GSE39582", "TCGA_test", "EMTAB12862", "ws4_spinal")

# ==========================================
# 2. Helper Functions
# ==========================================
get_inf_stats <- function(time, status, score, tau, n_iter) {
  res <- tryCatch({
    mydata <- matrix(c(as.numeric(time), as.numeric(status), as.numeric(score)), ncol = 3)
    dimnames(mydata) <- NULL
    
    # Use 'itr' per your package version requirement
    out <- survC1::Inf.Cval(mydata, tau, itr = n_iter)
    
    # Handle naming differences in package versions (Dhat vs D)
    d_val <- if (!is.null(out$Dhat)) out$Dhat else if (!is.null(out$D)) out$D else NA
    
    list(D = d_val, L = out$low95, U = out$upp95)
  }, error = function(e) {
    message(paste("  !! Inf.Cval Internal Error:", e$message))
    return(list(D = NA, L = NA, U = NA))
  })
  return(res)
}

# ==========================================
# 3. Refined Prediction Function
# ==========================================
external_prediction_stratified <- function(signature, cohorts) {
  
  # --- 3a. Setup & Loading ---
  train <- readRDS(file.path(cox_models_inputpath, signature, "train_all.rds"))
  signatures_annot <- readRDS(signatures_annot_inputpath)
  signature_genes_tcga_pool <- setdiff(colnames(train), cvrt)
  
  all_cohort_genes <- list()
  all_cohort_genes[["tcga"]] <- signature_genes_tcga_pool
  
  ext_cohorts_counts <- list()
  ext_cohorts_clin   <- list()
  cox_models          <- list() 
  true_sig_genes_by_cohort <- list() 
  
  for (cohort in cohorts) {
    print(paste0("Loading data for ", signature, " in ", cohort))
    
    test_clinical <- readRDS(file.path(external_cohorts_inputpath, paste0(cohort, "_clinical.rds")))
    names(test_clinical) <- make.names(names(test_clinical))
    test_counts_all <- readRDS(file.path(external_cohorts_inputpath, paste0(cohort, "_counts.rds"))) %>% t()
    
    sig <- signatures_annot %>% filter(cohort_name == cohort, signature_name == signature)
    
    # ID mapping logic
    id_map <- list(
      EMTAB12862 = "ensembl_gene_id", ws4_spinal = "ensembl_gene_id",
      GSE17536 = "ensembl_gene_id", GSE17537 = "ensembl_gene_id", 
      GSE29621 = "ensembl_gene_id", GSE39582 = "ensembl_gene_id", 
      TCGA_test  = "ensembl_gene_id"
    )
    id_type <- id_map[[cohort]]; if (is.null(id_type)) id_type <- "affy_hg_u133_plus_2"
    
    all_cohort_genes[[cohort]] <- colnames(test_counts_all)
    sig_ids_true <- sig[[id_type]][!is.na(sig[[id_type]]) & sig[[id_type]] %in% colnames(test_counts_all)]
    true_sig_genes_by_cohort[[cohort]] <- sig_ids_true
    
    common_sig_genes <- intersect(signature_genes_tcga_pool, sig_ids_true)
    common_cvrts     <- intersect(cvrt, names(test_clinical))
    
    # Local Model Training
    train_subset <- train[, c(common_cvrts, common_sig_genes)]
    cox_models[[cohort]] <- list(
      all       = coxph(Surv(SurvTime, vital_status) ~ ., data = train_subset, x = TRUE),
      cvrts     = coxph(Surv(SurvTime, vital_status) ~ ., data = train_subset[, common_cvrts], x = TRUE),
      signature = coxph(Surv(SurvTime, vital_status) ~ ., data = train_subset[, c("vital_status", "SurvTime", common_sig_genes)], x = TRUE)
    )
    
    ext_cohorts_counts[[cohort]] <- as.data.frame(test_counts_all)
    ext_cohorts_clin[[cohort]]   <- test_clinical
  }
  
  # --- 3b. Batch Correction (ComBat) ---
  all_common_genes <- Reduce(intersect, all_cohort_genes)
  cohorts_to_correct <- setdiff(cohorts, "TCGA_test")
  
  # Prepare ComBat Matrix
  train_exp <- train %>% dplyr::select(all_of(all_common_genes)) %>% mutate(cohort = "tcga")
  mat_list <- list(train_exp)
  for (cohort in cohorts_to_correct) {
    ext_data <- ext_cohorts_counts[[cohort]] %>% dplyr::select(all_of(all_common_genes)) %>% mutate(cohort = cohort)
    mat_list <- append(mat_list, list(ext_data))
  }
  
  mat_to_correct <- do.call(rbind.data.frame, mat_list)
  meta_correct   <- data.frame(cohort = mat_to_correct$cohort, row.names = rownames(mat_to_correct))
  
  # Adjust everything EXCEPT TCGA_test to match the TCGA training baseline
  combat_corrected <- ComBat(dat = t(mat_to_correct %>% dplyr::select(-cohort)), 
                             batch = meta_correct$cohort,
                             mod = model.matrix(~1, data = meta_correct), 
                             ref.batch = "tcga",
                             par.prior = TRUE) %>% 
    t() %>% as.data.frame() %>% mutate(cohort = meta_correct$cohort)
  
  # Merge TCGA_test back in (keep it raw as it is already on the reference platform)
  if ("TCGA_test" %in% cohorts) {
    tcga_test_unc <- ext_cohorts_counts[["TCGA_test"]] %>% 
      dplyr::select(all_of(all_common_genes)) %>% 
      mutate(cohort = "TCGA_test")
    combat_results <- rbind(combat_corrected, tcga_test_unc)
  } else {
    combat_results <- combat_corrected
  }
  
  # --- 3c. Stratified & Combined Uno Index Calculation ---
  final_performance <- list()
  
  for (cohort in cohorts) {
    print(paste0(">>> Processing Uno Index for ", cohort))
    
    sig_genes_final <- intersect(true_sig_genes_by_cohort[[cohort]], all_common_genes)
    cohort_exp      <- combat_results %>% filter(cohort == !!cohort) %>% dplyr::select(all_of(sig_genes_final))
    cohort_clin     <- ext_cohorts_clin[[cohort]]
    
    # Merge clinical and corrected expression
    full_data <- cbind.data.frame(cohort_clin, cohort_exp) %>%
      mutate(vital_status = ifelse(vital_status %in% c("Alive", "0", FALSE, 0), 0, 1)) %>%
      filter(SurvTime > 0)
    
    # Variance Diagnostic
    vars <- apply(full_data[, sig_genes_final, drop=FALSE], 2, var, na.rm=TRUE)
    cat(sprintf("  Diagnosis: MaxTime=%.1f | Genes with 0 variance: %d/%d\n", 
                max(full_data$SurvTime, na.rm=TRUE), sum(vars == 0, na.rm=TRUE), length(sig_genes_final)))
    
    full_data$MSI_Group <- ifelse(as.numeric(as.character(full_data$MSI.Status_MSI.L.MSS)) == 1, "MSS", "MSI")
    
    # Loop through Combined, MSI, and MSS
    for (m_status in c("Combined", "MSI", "MSS")) {
      
      if (m_status == "Combined") {
        sub_data <- full_data
      } else {
        sub_data <- full_data %>% filter(MSI_Group == m_status)
      }
      
      # Survival time must reach tau_fixed for Uno's index to be valid
      if(nrow(sub_data) > 10 && sum(sub_data$vital_status) > 1 && max(sub_data$SurvTime) >= tau_fixed) {
        
        lp_all  <- predict(cox_models[[cohort]]$all, sub_data)
        lp_cvrt <- predict(cox_models[[cohort]]$cvrts, sub_data)
        lp_sig  <- predict(cox_models[[cohort]]$signature, sub_data)
        
        stats_all <- get_inf_stats(sub_data$SurvTime, sub_data$vital_status, lp_all, tau_fixed, inf_iter)
        stats_cv  <- get_inf_stats(sub_data$SurvTime, sub_data$vital_status, lp_cvrt, tau_fixed, inf_iter)
        stats_sig <- get_inf_stats(sub_data$SurvTime, sub_data$vital_status, lp_sig, tau_fixed, inf_iter)
        
        if (!is.na(stats_all$D)) {
          final_performance[[paste0(cohort, "_", m_status)]] <- data.frame(
            Cohort = cohort,
            MSI_Status = m_status,
            Model = c(paste0(signature, " + cvrts"), "cvrts", signature),
            C = c(stats_all$D, stats_cv$D, stats_sig$D),
            ci_l = c(stats_all$L, stats_cv$L, stats_sig$L),
            ci_u = c(stats_all$U, stats_cv$U, stats_sig$U),
            N = nrow(sub_data)
          )
        }
      } else {
        message(paste0("  Skipping group: ", m_status, " (Insufficient data or max follow-up < tau)"))
      }
    }
  }
  
  return(data.table::rbindlist(final_performance))
}

# ==========================================
# 4. Execution & Saving
# ==========================================
performance_results <- external_prediction_stratified(selected_signature, cohorts)

# Final Save
saveRDS(performance_results, file = rds_output)

cat("\n==========================================\n")
cat("Success! Uno Index results saved to:\n", rds_output, "\n")
print(head(performance_results))
