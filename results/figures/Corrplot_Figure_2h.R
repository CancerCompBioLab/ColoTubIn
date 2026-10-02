library(here)
here()
source(here("environment", "requirements.R"))

# Inputpaths
signatures_inputpath <- here("data", "processed_data", "literatureSignatures_new.rds")
counts_inputpath <- here("data", "processed_data", "counts_norm_coad_patients.rds")


# Load data
coad_signatures <- readRDS(signatures_inputpath)
counts <- readRDS(counts_inputpath) %>% t()

# Arguments
method <- "plage"

############# Added
# Convert ENSEMBL to SYMBOL
genes <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = rownames(counts),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

rownames(counts)<-genes

# Convert to
signatures_coad_symbol<-coad_signatures
# signatures_coad_symbol <- list()
# for (i in seq_along(coad_signatures)) {
#   signatures_coad_symbol[[i]] <- AnnotationDbi::mapIds(
#     org.Hs.eg.db,
#     keys = coad_signatures[[i]],
#     column = "ENSEMBL",
#     keytype = "SYMBOL",
#     multiVals = "first"
#   )
# }
#names(signatures_coad_symbol) <- names(coad_signatures)

# for (s in seq_along(signatures_coad_symbol)) {
#   signatures_coad_symbol[[s]] <- names(signatures_coad_symbol[[s]])
# }

# Run GSVA
gsva <- gsva(
  counts,
  signatures_coad_symbol,
  method = method,
  kcdf = "Gaussian"
)
gsva <- as.data.frame(t(gsva))
names(gsva)[grep("ours", names(gsva))] <- "ColoTubIn"

# order by SigCheck significance
s <- readRDS(here("data", "processed_data", "sigcheck_literature", "scKnownLiterature_ours.rds"))

sort <- names(sort(c(ours = s$survivalPval, s$survivalPvalsKnown))) # This sorting doesnt work like this 

sort<-unique(sort)
sort [sort  == "ours"] <- "ColoTubIn"
gsva <- gsva[, match(sort, colnames(gsva))]


output_dir <- here("results", "figures")
pdf(file.path(output_dir, "Corrplot_Figure_2h.pdf"), width = 6, height = 6)

# Plotting
corrplot(
  cor(gsva),
  method = "circle",
  type = "lower",
  col = viridis(10),
  tl.col = "black"
)

dev.off()
