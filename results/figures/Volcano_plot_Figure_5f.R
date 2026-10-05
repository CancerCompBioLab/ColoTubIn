####################################
# Figure 5f
####################################

library(here)
here()
source(here("environment", "requirements.R"))


# Volcano plots 

## --- Configuration ---
fgsea_inputpath <- here("data", "data", "repurposing","fgsea_by_signature_new.rds")

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

## --- Generate Plots ---

pp <- list()
legend_grob <- NULL 

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
  
  plot_j <- plot_j + theme(legend.position = "none")
  
  pp[[j]] <- plot_j
}

pp[[length(pp) + 1]] <- ggpubr::as_ggplot(legend_grob)

## --- Final Arrangement and Output ---

num_plots <- length(pp)
NCOL <- 3 
NROW <- ceiling(num_plots / 3)

output_dir <- here("results", "figures")

pdf(file.path(output_dir, "Volcano_Plot_Figure_5f.pdf"), width = 14, height = 8)

# Arrange the plots
plot_arranged <- ggpubr::ggarrange(
  plotlist = pp, 
  ncol = NCOL,   
  nrow = NROW,   
  common.legend = FALSE 
)

# Print and save
print(plot_arranged)
dev.off()
