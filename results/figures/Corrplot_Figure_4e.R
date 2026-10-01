# Plot correlation plot of coad signatures and cell cycle signatures

library(here)
here()

source(here("environment", "requirements.R"))

# Input paths
cc_signatures_inputpath <- here("data", "processed_data", "cell_cycle.rds")
tcga_inputpath <- here("data", "processed_data", "counts_norm_coad_patients.rds")
coad_signatures_inputpath <- here("data", "processed_data", "literatureSignatures_new.rds")

# Load data
train <- readRDS(tcga_inputpath)
cc_signatures <- readRDS(cc_signatures_inputpath)
coad_signatures <- readRDS(coad_signatures_inputpath)


############# Added
# Convert ENSEMBL to SYMBOL
genes <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = colnames(train),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

colnames(train)<-genes

cc_signatures <- lapply(cc_signatures, function(x) names(x))


# Define signatures
signatures <- c(cc_signatures, coad_signatures)

# Calculate signature´s activity
gsva <- gsva(
  t(train),
  signatures,
  method = "plage",
  kcdf = "Gaussian")
gsva <- as.data.frame(t(gsva))

signatures_ordered <- c(
  "G0.EarlyG1",
  "G1Phase",
  "G1S.transition",
  "S.phase",
  "G2.phase",
  "G2M.transition",
  "M.prophase",
  "M.prometaphase",
  "M.metaphase.anaphase",
  "M.telophase.cytokinesis",
  "ours",
  "MDA114",
  "zhang",
  "ColoGuideEx",
  "ColoPrint",
  "ColoGuidePro",
  "ICR"
)

#gsva <- gsva[, match(signatures_ordered, names(gsva))]

valid <- signatures_ordered[signatures_ordered %in% names(gsva)]

gsva <- gsva[, valid]


#Added
names(gsva)[grep("ours", names(gsva))] <- "ColoTubIn"

output_dir <- here("results", "figures")

# pdf(file.path(output_dir,"Corrplot_Figure_4e.pdf"), width = 6, height = 6)
# # Plotting
corrplot::corrplot(
  cor(gsva),
  #method = "number",
  type = "lower",
  col = viridis(10),
  tl.col = "black")

dev.off()