library(here)
here()
source(here("environment", "requirements.R"))


# Figure 2 main script

# Define central data and output directories using here
sigcheck_dir <- here("data", "processed_data", "sigcheck")
output_dir <- here("results", "figures")

# ColoTubIn
# ===
random <- readRDS(file.path(sigcheck_dir, 'sigCheckRandom_ours.rds'))
features <- readRDS(file.path(sigcheck_dir, 'sigCheckPermutedFeatures_ours.rds'))
survival <- readRDS(file.path(sigcheck_dir, 'sigCheckPermutedSurvival_ours.rds'))
msigdb <- readRDS(file.path(sigcheck_dir, 'scKnownMSigDB_ours.rds'))
literature <- readRDS(file.path(sigcheck_dir, 'scKnownLiterature_ours.rds'))

pdf(file.path(output_dir, "SigCheck_Figure_2d.pdf"), width = 6, height = 6)
sigCheckPlot(random, title = 'Random signatures', nolegend = TRUE)
dev.off()

pdf(file.path(output_dir, "SigCheck_Figure_2e.pdf"), width = 6, height = 6)
sigCheckPlot(features, title = 'Permuted features', nolegend = TRUE)
dev.off()

pdf(file.path(output_dir, "SigCheck_Figure_2f.pdf"), width = 6, height = 6)
sigCheckPlot(survival, title = 'Permuted survival', nolegend = TRUE)
dev.off()

pdf(file.path(output_dir, "SigCheck_Figure_2g.pdf"), width = 6, height = 6)
sigCheckPlot(msigdb, title = 'MSigDB', nolegend = TRUE)
dev.off()

literature <- readRDS(file.path(sigcheck_dir, 'scKnownLiterature_ours.rds'))

pdf(file.path(output_dir, "SigCheck_Figure_2i.pdf"), width = 6, height = 6)
sigCheckPlot(literature, title = 'Literature signatures', nolegend = TRUE)
dev.off()



# Supplementary Material Figure 3

# --- Updated loop to fit in page and save securely ---

pdf(file.path(output_dir, "Signature_Plots_Suppl_Figure_3.pdf"), 
    width = 8.27, height = 11.69)

sigs <- c('ColoPrint', 'ColoGuideEx', 'ColoGuidePro', 'MDA114', 'zhang', 'ICR')
types <- c('Random', 'PermutedFeatures', 'PermutedSurvival')
titles <- c('Random', 'Features', 'Survival')

# List files from the local relative sigcheck directory
actual_files <- list.files(path = sigcheck_dir, pattern = "\\.rds$")

# mar = c(bottom, left, top, right)
par(mfrow = c(6, 3), 
    mar = c(5, 3, 3, 1), 
    oma = c(1, 1, 1, 1), 
    mgp = c(2, 0.6, 0))

for (s in sigs) {
  for (i in 1:3) {
    target_name <- paste0('sigCheck', types[i], '_', s, '.rds')
    match_idx <- which(tolower(actual_files) == tolower(target_name))
    
    if (length(match_idx) > 0) {
      # Safely load using file.path with sigcheck_dir
      data_obj <- readRDS(file.path(sigcheck_dir, actual_files[match_idx[1]]))
      
      sigCheckPlot(data_obj, 
                   title = paste(s, "-", titles[i]), 
                   nolegend = TRUE)
      
    } else {
      plot.new()
      text(0.5, 0.5, paste("Missing:\n", s), cex = 0.6)
    }
  }
}

dev.off()
