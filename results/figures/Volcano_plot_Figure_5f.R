
setwd("/mnt/CCBdata/users/lea/Project_ColoTubIn/Github")
library(here)
here()
source(here("environment", "requirements.R"))

# Inputpaths
tcga_train_counts_inputpath <- "/CCBdata/users/lea/YACCS_Linares/From_Linares/data/data_partitions/counts_train.rds"
tcga_train_meta_inputpath <- "/CCBdata/users/lea/YACCS_Linares/From_Linares/data/data_partitions/metadata_train.rds"
signatures_coad_inputpath <- "/CCBdata/users/lea/YACCS_Linares/From_Linares/extdata/signatures/signatures_coad.rds"
ccle_gex_inputpath <- "/CCBdata/users/lea/YACCS_Linares/From_Linares/extdata/repurposing/CCLE_expression.csv"
ccle_info_inputpath <- "/CCBdata/users/lea/YACCS_Linares/From_Linares/extdata/repurposing/sample_info.csv"

# Outputpath
outputpath <- "/CCBdata/users/lea/YACCS_Linares/From_Linares/data/repurposing/signatures_scores_new.rds"

# Load data
tcga_train <- readRDS(tcga_train_counts_inputpath)
tcga_train_meta <- readRDS(tcga_train_meta_inputpath)

signatures_coad <- readRDS(signatures_coad_inputpath)

#### Added ######

# Assuming the list is named signatures_coad and the zhang entry is already mapped:

# 1. Access the mapped vector for zhang
zhang_ensembl_named <- signatures_coad$zhang

# 2. Use unname() to remove the Gene Symbol names from the vector elements
signatures_coad$zhang <- unname(zhang_ensembl_named)

# 3. Print the result to verify the change
signatures_coad$zhang

# Add more signatures:

# ICR 

icr_ensembl_ids <- c(
  "ENSG00000271503", "ENSG00000120217", "ENSG00000153563", "ENSG00000172116", 
  "ENSG00000163599", "ENSG00000169245", "ENSG00000138755", "ENSG00000049768",
  "ENSG00000115523", "ENSG00000145649", "ENSG00000100453", "ENSG00000100450", 
  "ENSG00000131203", "ENSG00000111537", "ENSG00000113302", "ENSG00000125347",
  "ENSG00000188389", "ENSG00000180644", "ENSG00000115415", "ENSG00000073861"
)

# 2. Assign the vector to the 'ICR' element in your main list
# Note: If 'ICR' already exists, this line will overwrite it. 
signatures_coad[["ICR"]] <- icr_ensembl_ids

##################################

ccle_info <- read.delim2(
  ccle_info_inputpath,
  header = TRUE,
  row.names = 1,
  sep = ","
)

ccle_gex <- data.table::fread(
  ccle_gex_inputpath,
  header = TRUE,
  sep = ",") %>%
  tibble::column_to_rownames("V1")

# select coad cell lines
ccle_info <-
  ccle_info %>%
  filter(
    primary_disease == "Colon/Colorectal Cancer"
  )
ccle_gex <- ccle_gex[match(rownames(ccle_info), rownames(ccle_gex)), ]
ccle_gex <- ccle_gex[complete.cases(ccle_gex), ]

colnames(ccle_gex) <- sapply(strsplit(colnames(ccle_gex), " (", fixed = TRUE), "[", 1)

# 
# # Loop
i <- 4
res <- list()
for (i in seq_along(signatures_coad)) {
  
  signature_name <- names(signatures_coad)[i]
  signature_symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = signatures_coad[[i]],
    column = c("SYMBOL"),
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  signature_symbol <- make.names(signature_symbol)
  names(signature_symbol) <- signatures_coad[[i]]
  
  
  diff_genes <- setdiff(signature_symbol, colnames(ccle_gex))
  if (length(diff_genes) > 0) {
    
    to_remove <- which(signature_symbol %in% diff_genes)
    signature_symbol <- signature_symbol[-to_remove]
    signatures_coad[[i]] <- signatures_coad[[i]][-to_remove]
    
  }
  
  ccle_gex_f <- ccle_gex[, signature_symbol]
  ccle_gex_f$cohort <- "ccle"
  
  tcga_train_f <- tcga_train[, signatures_coad[[i]]]
  colnames(tcga_train_f) <- signature_symbol
  tcga_train_f <- scale(tcga_train_f, scale = FALSE) %>% as.data.frame()
  tcga_train_f$cohort <- "tcga"
  
  # Correct by combat
  df <- rbind.data.frame(tcga_train_f, ccle_gex_f)
  mat <- df %>% dplyr::select(-c(cohort))
  meta <- df %>% dplyr::select(c(cohort))
  
  modcombat <- model.matrix(~ 1, data = meta)
  combat <- ComBat(
    dat = t(mat),
    batch = meta$cohort,
    mod = modcombat,
    ref.batch = "tcga"
  )
  
  combat <- cbind.data.frame(t(combat), cohort = meta$cohort)
  
  train <- combat[which(combat$cohort == "tcga"), ]
  train <- subset(train, select = -c(cohort))
  test <- combat[which(combat$cohort == "ccle"), ]
  test <- subset(test, select = -c(cohort))
  
  # Train cox model
  train <- cbind.data.frame(
    tcga_train_meta %>% 
      dplyr::select(c(vital_status, SurvTime)) %>%
      mutate(
        vital_status = ifelse(vital_status == "Dead", TRUE, FALSE)
      ),
    train
  )
  
  cox <- coxph(Surv(SurvTime, vital_status) ~ ., data = train, x = TRUE)
  pred <- predict(cox, test)
  
  res[[i]] <- data.frame(
    cell_lines = names(pred),
    score = pred,
    signature = signature_name
  )
  
  print(paste0("Predicted using ", signature_name))
  
}

# Bind results
res_df <- data.table::rbindlist(res)

# Save
saveRDS(res_df, file = outputpath)

###############################



# Run Drug enrichment analysis
# ---------
#source("requirements.R")

# Inputpath
ccle_inputpath <- "/CCBdata/users/lea/YACCS_Linares/From_Linares/extdata/repurposing/sample_info.csv"
prism_inputpath <- "/CCBdata/users/lea/YACCS_Linares/From_Linares/extdata/repurposing/PRISM_19Q4_secondary-screen-dose-response-curve-parameters.csv"
scores_inputpath <- "/CCBdata/users/lea/YACCS_Linares/From_Linares/data/repurposing/signatures_scores_new.rds"

# Outputpath
outputpath <- "/CCBdata/users/lea/YACCS_Linares/From_Linares/data/repurposing/fgsea_by_signature_new.rds"

# Load data
scores <- readRDS(scores_inputpath) %>% as.data.frame()

# ccle_samlpe <- read.table(
#     ccle_inputpath,
#     header = TRUE,
#     sep = ",",
#     quote = '"',
#     fill = TRUE,
#     stringsAsFactors = FALSE,
#     row.names = NULL
#     )

prism2 <- read.table(
  prism_inputpath,
  header = TRUE,
  sep = ",",
  quote = '"',
  fill = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

prism2 <- prism2[, c("broad_id", "depmap_id", "auc", "moa", "name", "screen_id")]

prism2 <- reshape2::dcast(
  dat = prism2,
  formula = screen_id + name + broad_id + moa ~ depmap_id,
  fun.aggregate = sum,
  value.var = "auc"
)

prism2[prism2 == 0] <- NA

# remove duplicates
prism2 <-
  prism2 %>%
  # group_by(name) %>%
  # arrange(factor(screen_id, levels = c("MTS010", "MTS006", "MTS005", "HTS002"))) %>%
  # dplyr::slice(1) %>%
  # ungroup() %>%
  dplyr::select(-c(screen_id, name))

prism2_m <- as.matrix(prism2[, -c(1, 2)])
rownames(prism2_m) <- prism2$broad_id
prism2_m <- t(prism2_m)


# Load scores from each signature
# =================================
signature <- unique(scores$signature)
fgseaRes <- list()
# i <- 2
for (i in seq_along(signature)) {
  
  print(paste0("Running ", signature[i]))
  scores_f <- scores[which(scores$signature == signature[i]), ]
  signature_score <- scores_f$score
  names(signature_score) <- scores_f$cell_lines
  
  z1 <- intersect(names(signature_score), rownames(prism2_m))
  cor_z1 <- apply(prism2_m[z1, ], 2, function(x) {
    cor(x, signature_score[z1], use = "complete.obs", method = "sp")
  })
  
  # Establishment of drug type list
  # ===========================
  drug_type1 <- split(prism2$broad_id[which(!is.na(prism2$moa))],prism2$moa[which(!is.na(prism2$moa))])
  drug_types <- strsplit(prism2$moa, split = ";")
  names(drug_types) <- prism2$broad_id
  ctrp_cmp_ext = data.frame(broad_id = names(unlist(drug_types)), moa=unlist(drug_types),stringsAsFactors = F)
  drug_type2 = split(ctrp_cmp_ext$broad_id[which(!is.na(ctrp_cmp_ext$moa))],
                     ctrp_cmp_ext$moa[which(!is.na(ctrp_cmp_ext$moa))])
  drug_types = strsplit(ctrp_cmp_ext$moa, split = ", ")
  names(drug_types) = ctrp_cmp_ext$broad_id
  ctrp_cmp_ext2 = data.frame(broad_id = names(unlist(drug_types)), moa=unlist(drug_types),stringsAsFactors = F)
  drug_type3 = split(ctrp_cmp_ext2$broad_id[which(!is.na(ctrp_cmp_ext2$moa))],
                     ctrp_cmp_ext2$moa[which(!is.na(ctrp_cmp_ext2$moa))])
  
  
  # Run GSEA analysis
  # ===========================
  set.seed(1993)
  fgseaRes[[i]] <- fgsea::fgsea(
    pathways = drug_type3,
    stats = sort(cor_z1), sampleSize = 5,
    minSize = 1,
    maxSize = Inf, 
    nproc = 3
  )
  
  fgseaRes[[i]]$log10adjpval <- -log10(fgseaRes[[i]]$padj)
}
names(fgseaRes) <- signature

saveRDS(fgseaRes, file = outputpath)

###################

# Volcano plots figure 5f

# Load required libraries
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggrepel)

## --- Configuration ---
fgsea_inputpath <- "/mnt/CCBdata/users/lea/YACCS_Linares/From_Linares/data/repurposing/fgsea_by_signature_new.rds"

# Custom colors
COLOR_SIGNIFICANT <- "#23B880"      # Green
COLOR_NON_SIGNIFICANT <- "#482677FF" # Dark Purple

## --- 1. Load and Prepare Data ---
fgseaRes <- readRDS(fgsea_inputpath)

# Rename "ours" to "ColoTubIn"
name_to_change_index <- which(names(fgseaRes) == "ours")
if (length(name_to_change_index) > 0) {
  names(fgseaRes)[name_to_change_index] <- "ColoTubIn"
}

# Move "ColoTubIn" to the end of the list
colotubin_idx <- which(names(fgseaRes) == "ColoTubIn")
if (length(colotubin_idx) > 0) {
  other_indices <- setdiff(seq_along(fgseaRes), colotubin_idx)
  fgseaRes <- fgseaRes[c(other_indices, colotubin_idx)]
}

## --- 2. Generate Plots ---

pp <- list()
legend_grob <- NULL # We will store the extracted legend here

for (j in seq_along(fgseaRes)) {
  current_name <- names(fgseaRes)[j]
  
  # Prepare data: Identify significance and create labels
  df <- fgseaRes[[j]] %>%
    mutate(
      is_significant = padj < 0.05, 
      label_pathway = ifelse(is_significant, pathway, NA)
    )
  
  plot_j <- ggplot(
    df,
    aes(
      x = NES, 
      y = log10adjpval, 
      size = size,          
      color = is_significant,
      alpha = is_significant 
    )
  ) +
    geom_point() +
    geom_text_repel(
      aes(label = label_pathway),
      size = 3, 
      color = "black", 
      max.overlaps = 20, 
      segment.alpha = 0.5, 
      segment.colour = "gray", 
      min.segment.length = 0,
      show.legend = FALSE 
    ) +
    scale_color_manual(
      values = c("TRUE" = COLOR_SIGNIFICANT, "FALSE" = COLOR_NON_SIGNIFICANT),
      labels = c("TRUE" = "Significant", "FALSE" = "Non-Significant"),
      name = "Significance"
    ) +
    scale_size_continuous(
      range = c(1, 6), 
      name = "Gene Set Size"
    ) +
    scale_alpha_manual(
      values = c("TRUE" = 1, "FALSE" = 0.4), 
      guide = "none"
    ) +
    ylim(c(0, 6)) +
    xlim(c(-2.5, 2.5)) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      # Keep the legend temporarily so we can extract it
      legend.position = "right"
    ) +
    labs(
      title = current_name,
      y = "-log10(adj. p-val)",
      x = "NES"
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 4))
    )
  
  # Extract the legend from the first plot
  if (is.null(legend_grob)) {
    legend_grob <- ggpubr::get_legend(plot_j)
  }
  
  # Now remove the legend from the plot before adding it to the list
  plot_j <- plot_j + theme(legend.position = "none")
  
  pp[[j]] <- plot_j
}

# --- NEW: Append the extracted legend as its own panel at the end ---
# as_ggplot converts the legend grob into a plot object ggarrange can use
pp[[length(pp) + 1]] <- ggpubr::as_ggplot(legend_grob)

## --- 3. Final Arrangement and Output ---

num_plots <- length(pp) # This now includes the legend panel
NCOL <- 3 
NROW <- ceiling(num_plots / 3)

output_dir <- "/mnt/CCBdata/users/lea/YACCS_Linares/Final_Figures"

pdf(file.path(output_dir, "Volcano_Plot_Figure_6_updated_order_legend.pdf"), width = 14, height = 8)

# Arrange the plots
plot_arranged <- ggpubr::ggarrange(
  plotlist = pp, 
  ncol = NCOL,   
  nrow = NROW,   
  common.legend = FALSE # Kept false because the legend is treated as a normal plot panel now
)

# Print and save
print(plot_arranged)
dev.off()