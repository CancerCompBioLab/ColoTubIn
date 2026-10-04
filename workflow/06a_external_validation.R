# SurvC1 with confidence interval 

library(here)
here()
source(here("environment", "requirements.R"))

# --- CONFIGURATION ---
# Set to 10 for a quick test, 1000 for final results
inf_iter <- 1000

# Inputpaths
external_cohorts_inputpath <- here::here("data", "data", "pp", "geo")
cox_models_inputpath <- here::here("data", "data", "cox", "tcga_new", "signatures")
signatures_annot_inputpath <- here::here("data","extdata", "signatures", "signatures_annot_w_EMTAB_ws4_spinal_tcga_test_ICR_signature.rds")

# Outputpath
outputpath <- here::here("data", "data", "cox", "external", "uno_index_1825_days.rds")


# Arguments
cvrt <- c(
  "vital_status", "SurvTime", "gender_male", "age_at_diagnosis",
  "Stage_II", "Stage_III", "Stage_IV", "MSI.Status_MSI.L.MSS", "Molecular_Subtype_noCIN"
)

signatures <- c("coloGuideEx", "coloGuidePro", "coloPrint", "mda114", "ours", "zhang", "ICR")

cohorts <- c("GSE17536", "GSE17537", "GSE29621", "GSE39582", "TCGA_test", "EMTAB12862", "ws4_spinal")


# --- MAIN PREDICTION FUNCTION ---
external_prediction <- function(signature, cohorts) {
  
  # 1. Load Training Data
  train <- readRDS(file.path(cox_models_inputpath, signature, "train_all.rds"))
  signatures_annot <- readRDS(signatures_annot_inputpath)
  
  # Base time horizon
  max_train_time <- 1825
  
  signature_genes_tcga_pool <- setdiff(colnames(train), cvrt)
  all_cohort_genes <- list()
  all_cohort_genes[["tcga"]] <- signature_genes_tcga_pool
  
  ext_cohorts_counts <- list()
  ext_cohorts_clin   <- list()
  cox_models         <- list()
  true_sig_genes     <- list() 
  
  # 2. Load External Cohorts and Train Models
  for (cohort in cohorts) {
    print(paste0("Loading data for ", signature, " in ", cohort))
    
    test_clinical <- readRDS(file.path(external_cohorts_inputpath, paste0(cohort, "_clinical.rds")))
    names(test_clinical) <- make.names(names(test_clinical))
    test_counts_all <- readRDS(file.path(external_cohorts_inputpath, paste0(cohort, "_counts.rds"))) %>% t()
    
    sig_map <- signatures_annot %>% filter(cohort_name == cohort, signature_name == signature)
    
    id_map <- list(
      EMTAB12862 = "ensembl_gene_id", 
      TCGA_test  = "ensembl_gene_id",
      ws4_spinal = "ensembl_gene_id", 
      GSE17536   = "ensembl_gene_id",
      GSE17537   = "ensembl_gene_id", 
      GSE29621   = "ensembl_gene_id", 
      GSE39582   = "ensembl_gene_id"
    )
    id_type <- id_map[[cohort]]; if (is.null(id_type)) id_type <- "aff_hg_u133_plus_2"
    
    all_cohort_genes[[cohort]] <- colnames(test_counts_all)
    sig_ids <- sig_map[[id_type]]
    valid <- !is.na(sig_ids) & sig_ids %in% colnames(test_counts_all)
    sig_ids_true <- sig_ids[valid]
    true_sig_genes[[cohort]] <- sig_ids_true
    
    test_counts_sig <- test_counts_all[, sig_ids_true, drop = FALSE]
    common_sig_genes <- intersect(signature_genes_tcga_pool, colnames(test_counts_sig))
    common_cvrts     <- intersect(cvrt, names(test_clinical))
    
    train_subset <- train[, c(common_cvrts, common_sig_genes)]
    
    cox_models[[cohort]] <- list(
      all = coxph(Surv(SurvTime, vital_status) ~ ., data = train_subset, x = TRUE),
      cvrts = coxph(Surv(SurvTime, vital_status) ~ ., data = train_subset[, common_cvrts], x = TRUE),
      signature = coxph(Surv(SurvTime, vital_status) ~ ., data = train_subset[, c("vital_status", "SurvTime", common_sig_genes)], x = TRUE)
    )
    
    ext_cohorts_counts[[cohort]] <- test_counts_all %>% as.data.frame()
    ext_cohorts_clin[[cohort]]   <- test_clinical 
  }
  
  # 3. ComBat Batch Correction
  all_common_genes <- Reduce(intersect, all_cohort_genes)
  cohorts_to_correct <- setdiff(cohorts, "TCGA_test")
  
  train_exp <- train %>% dplyr::select(all_of(all_common_genes)) %>% mutate(cohort = "tcga")
  mat_list <- list(train_exp)
  
  for (cohort in cohorts_to_correct) {
    ext_data <- ext_cohorts_counts[[cohort]] %>% dplyr::select(all_of(all_common_genes)) %>% mutate(cohort = cohort)
    mat_list <- append(mat_list, list(ext_data))
  }
  
  mat_full <- do.call(rbind.data.frame, mat_list)
  meta     <- data.frame(cohort = mat_full$cohort, row.names = rownames(mat_full))
  
  combat_corrected <- ComBat(dat = t(mat_full %>% dplyr::select(-cohort)), batch = meta$cohort, 
                             mod = model.matrix(~1, data = meta), ref.batch = "tcga", par.prior = TRUE) %>% 
    t() %>% as.data.frame() %>% mutate(cohort = meta$cohort)
  
  if ("TCGA_test" %in% cohorts) {
    tcga_test_unc <- ext_cohorts_counts[["TCGA_test"]] %>% dplyr::select(all_of(all_common_genes)) %>% mutate(cohort = "TCGA_test")
    final_combat_df <- rbind(combat_corrected, tcga_test_unc)
  } else { 
    final_combat_df <- combat_corrected 
  }
  
  # 4. Final Dataset Construction
  external_cohorts <- list()
  for (cohort in cohorts) {
    sig_genes_final <- intersect(true_sig_genes[[cohort]], all_common_genes)
    cohort_exp      <- final_combat_df %>% filter(cohort == !!cohort) %>% dplyr::select(all_of(sig_genes_final))
    
    external_cohorts[[cohort]] <- cbind.data.frame(ext_cohorts_clin[[cohort]], cohort_exp) %>% 
      mutate(vital_status = ifelse(vital_status == "Alive", 0, 1))
    
    # --- MOVE PRINT INSIDE HERE ---
    cat("Cohort:", cohort, "| Max Survival Time:", max(external_cohorts[[cohort]]$SurvTime, na.rm = TRUE), "\n")
  }
  
  # 5. Uno Index Inference
  uno_index_res <- list()
  for (cohort in cohorts) {
    
    print(paste0(">>> Processing Uno Index for ", cohort))
    
    sig_genes_present <- intersect(true_sig_genes[[cohort]], colnames(external_cohorts[[cohort]]))
    
    cat(paste0("\n[DIAGNOSIS] Cohort: ", cohort, " | Signature: ", signature, "\n"))
    cat(paste0("  - Expected Genes: ", length(true_sig_genes[[cohort]]), "\n"))
    cat(paste0("  - Found in Data:  ", length(sig_genes_present), "\n"))
    
    if (length(sig_genes_present) == 0) {
      cat("  !!! CRITICAL: No signature genes found in this cohort's data.\n")
    } else {
      # Check if the genes found actually vary across patients
      vars <- apply(external_cohorts[[cohort]][, sig_genes_present, drop=FALSE], 2, var, na.rm=TRUE)
      cat(paste0("  - Genes with 0 variance: ", sum(vars == 0, na.rm=TRUE), " out of ", length(sig_genes_present), "\n"))
    }
    # ---------------------------------------------
    
    lp_all <- predict(cox_models[[cohort]]$all, external_cohorts[[cohort]])
    lp_cv  <- predict(cox_models[[cohort]]$cvrts, external_cohorts[[cohort]])
    lp_sig <- predict(cox_models[[cohort]]$signature, external_cohorts[[cohort]])
    
    # Updated Helper with 'itr' and stability checks
    # Helper for Inf.Cval processing
    get_inf_stats <- function(time, status, score, tau, n_iter) {
      
      res <- tryCatch({
        mydata <- matrix(c(as.numeric(time), as.numeric(status), as.numeric(score)), ncol = 3)
        out <- survC1::Inf.Cval(mydata, tau, itr = n_iter)
        
        # Handle naming differences in package versions
        d_val <- if (!is.null(out$Dhat)) out$Dhat else if (!is.null(out$D)) out$D else NA
        
        list(D = d_val, L = out$low95, U = out$upp95)
      }, error = function(e) {
        message(paste("!! Inf.Cval Internal Error:", e$message))
        list(D = 0, L = NA, U = NA)
      })
      
      return(res)
    }
    
    stats_all <- get_inf_stats(external_cohorts[[cohort]]$SurvTime, external_cohorts[[cohort]]$vital_status, lp_all, max_train_time, inf_iter)
    stats_cv  <- get_inf_stats(external_cohorts[[cohort]]$SurvTime, external_cohorts[[cohort]]$vital_status, lp_cv,  max_train_time, inf_iter)
    stats_sig <- get_inf_stats(external_cohorts[[cohort]]$SurvTime, external_cohorts[[cohort]]$vital_status, lp_sig, max_train_time, inf_iter)
    
    uno_index_res[[cohort]] <- data.frame(
      Study = c(paste0(signature, " + cvrts"), "cvrts", signature),
      C     = c(stats_all$D, stats_cv$D, stats_sig$D),
      ci_l  = c(stats_all$L, stats_cv$L, stats_sig$L),
      ci_u  = c(stats_all$U, stats_cv$U, stats_sig$U),
      Subset = cohort,
      signature = signature
    )
  }
  return(data.table::rbindlist(uno_index_res))
}

# --- EXECUTION ---
res_list <- lapply(signatures, function(x) external_prediction(x, cohorts))
final_results <- data.table::rbindlist(res_list)

# Save
saveRDS(final_results, file = outputpath)
print(paste0("Success! Results with Confidence Intervals saved to ", outputpath))
