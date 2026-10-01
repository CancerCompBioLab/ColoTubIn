# Preparation 
library(here)
here()

source(here("environment", "requirements.R"))

# Preparation 

# Calculate pathway activity scores for TCGA-COAD patients and merge with clinical metadata
# =====================

# --- 1. Load Input Data using here() ---
tcga_inputpath <- here("data", "processed_data", "counts_norm_coad_patients.rds")
signatures_inputpath <- here("data", "processed_data", "signatures_coad.rds")
metadata_inputpath <- here("data", "processed_data", "metadata_coad_patients.rds")

# Load
gex <- readRDS(tcga_inputpath)
signatures <- readRDS(signatures_inputpath)
metadata <- readRDS(metadata_inputpath)

# ICR 
icr_ensembl_ids <- c(
  "ENSG00000271503", "ENSG00000120217", "ENSG00000153563", "ENSG00000172116",  
  "ENSG00000163599", "ENSG00000169245", "ENSG00000138755", "ENSG00000049768",
  "ENSG00000115523", "ENSG00000145649", "ENSG00000100453", "ENSG00000100450",  
  "ENSG00000131203", "ENSG00000111537", "ENSG00000113302", "ENSG00000125347",
  "ENSG00000188389", "ENSG00000180644", "ENSG00000115415", "ENSG00000073861"
)

# Assign the vector to the 'ICR' element in your main list
signatures[["ICR"]] <- icr_ensembl_ids

# Colon Cancer Cluster Subtypes (Guinney, Nat Med)
# =====================
m = mapIds(org.Hs.eg.db,
           keys = colnames(gex),
           column = 'ENTREZID',
           keytype = 'ENSEMBL')
colnames(gex) = m
gex_s = scale(gex)
clusters = CMSclassifier::classifyCMS(t(gex_s), method = 'SSP')[[3]]
metadata[["clusters"]] = clusters[["SSP.nearestCMS"]]

# Biomarkers Status
# =====================
query <- GDCquery(
  project      = "TCGA-COAD",
  data.category= "Simple Nucleotide Variation",
  data.type    = "Masked Somatic Mutation",
  access       = "open"
)

GDCdownload(query)
maf <- GDCprepare(query)

genes <- c("BRAF","KRAS","TP53","APC","SMAD4")
mutMat <- maf %>%
  transmute(
    sample = Tumor_Sample_Barcode,
    patient = substr(Tumor_Sample_Barcode, 1, 12),
    gene = Hugo_Symbol
  ) %>%
  filter(gene %in% genes) %>%
  dplyr::count(patient, gene, name = "n_mut") %>%         
  pivot_wider(names_from = gene, values_from = n_mut, values_fill = 0) %>%
  as.data.frame()

metadata <- metadata %>% left_join(mutMat, by = c("patient" = "patient"))
metadata[["BRAF"]] = ifelse(metadata[["BRAF"]] > 0, T, F)
metadata[["KRAS"]] = ifelse(metadata[["KRAS"]] > 0, T, F)
metadata[["TP53"]] = ifelse(metadata[["TP53"]] > 0, T, F)
metadata[["APC"]] = ifelse(metadata[["APC"]] > 0, T, F)
metadata[["SMAD4"]] = ifelse(metadata[["SMAD4"]] > 0, T, F)

# Calculate pathway activity
method = 'plage'

# Transform back as Signature contains ENSEMBL IDs
m = mapIds(org.Hs.eg.db,
           keys = colnames(gex),
           column = 'ENSEMBL',
           keytype = 'ENTREZID')
colnames(gex) = m

expr <- t(gex)

gsva.es <- gsva(
  expr = expr,
  gset.idx.list = signatures,
  method = "plage",
  verbose = FALSE
)

gsva <- t(gsva.es) |> as.data.frame()

# Convert GSVA to data frame and add sample ID column
gsva_df <- gsva %>%
  as.data.frame() %>%
  tibble::rownames_to_column("barcode")

metadata <- as.data.frame(metadata)

# Merge metadata with GSVA
merged_data <- metadata %>%
  dplyr::left_join(gsva_df, by = "barcode") %>%
  dplyr::rename(
    ColoTubIn = ours,
    kras = KRAS,
    braf = BRAF,
    tp53 = TP53,
    apc = APC,
    smad4 = SMAD4
  ) %>%
  mutate(   
    coloGuideEx = coloGuideEx * -1,
    coloGuidePro = coloGuidePro * -1,
    mda114 = mda114 * -1
  )

# Melt by GSVA columns (pathways) and retain metadata variables
gsva_columns <- colnames(gsva)

melted_data <- merged_data %>%
  pivot_longer(
    cols = c("coloGuideEx", "coloGuidePro", "coloPrint", "mda114", "ColoTubIn", "zhang","ICR"),
    names_to = "signature",
    values_to = "value"
  ) %>%
  mutate(
    MSI.Status = case_when(
      `MSI Status_MSI-L/MSS` == 1 ~ "MSI-L/MSS",
      `MSI Status_MSI-L/MSS` == 0 ~ "MSI-H",
      TRUE ~ NA_character_
    ),
    gender = ifelse(gender_male == 1, "male", "female"),
    Molecular_Subtype = ifelse(Molecular_Subtype_noCIN == 1, "noCIN", "CIN"),
    Stage = case_when(
      Stage_II == 1 ~ "II",
      Stage_III == 1 ~ "III",
      Stage_IV == 1 ~ "IV",
      TRUE ~ "I"
    )
  ) %>%
  dplyr::select(
    c(barcode, patient, vital_status, gender, age_at_diagnosis, Molecular_Subtype,
      Stage, MSI.Status, SurvTime, clusters, braf, kras, tp53, apc, smad4, signature, value
    )
  )

# --- 2. Save Processed Output using here() ---
saveRDS(melted_data, file = here("data", "processed_data", "TCGA-COAD_clinvars_sigscores.rds"))


# --- 3. Correlation with covariates ---
data <- readRDS(here("data", "processed_data", "TCGA-COAD_clinvars_sigscores.rds"))

sig <- as.character(unique(data[["signature"]]))
vars <- c('Molecular_Subtype', 'Stage', 'MSI.Status', 'clusters', 'kras', 'braf', 'tp53', 'apc', 'smad4')
res <- list()

for (i in seq_along(sig)) {
  x <- data[which(data[["signature"]] == sig[i]), ]
  res[[i]] <- data.frame(
    signatures = sig[i],
    Molecular_subtype = t.test(x[which(x[,vars[1]] == 'CIN'),][["value"]], x[which(x[,vars[1]] == 'noCIN'),][["value"]])$p.value,
    Stage = anova(lm(x[["value"]] ~ as.numeric(as.factor(x[["Stage"]]))))[[ "Pr(>F)" ]][1],
    MSI = t.test(x[which(x[,vars[3]] == 'MSI-H'),][["value"]], x[which(x[,vars[3]] == 'MSI-L/MSS'),][["value"]])$p.value,
    Consensus_clusters = kruskal.test(x[["value"]], x[["clusters"]])$p.value,
    kras = t.test(x[which(x[,vars[5]] == F),][["value"]], x[which(x[,vars[5]] == T),][["value"]])$p.value,
    braf = t.test(x[which(x[,vars[6]] == F),][["value"]], x[which(x[,vars[6]] == T),][["value"]])$p.value,
    tp53 = t.test(x[which(x[,vars[7]] == F),][["value"]], x[which(x[,vars[7]] == T),][["value"]])$p.value,
    apc = t.test(x[which(x[,vars[8]] == F),][["value"]], x[which(x[,vars[8]] == T),][["value"]])$p.value,
    smad4 = t.test(x[which(x[,vars[9]] == F),][["value"]], x[which(x[,vars[9]] == T),][["value"]])$p.value
  )
}

res <- as.data.frame(data.table::rbindlist(res))
rownames(res) <- res[["signatures"]]
res <- res[, -1]
res <- as.data.frame(t(res))


# --- 4. Make a Heatmap ---
library(ComplexHeatmap)
library(circlize)
library(viridis)

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

# # Open the PDF device, draw the heatmap, and close the device
# pdf(file.path(output_dir, "Heatmap_Figure_3c.pdf"), width = 10, height = 8)
# draw(final_heatmap)
# dev.off()