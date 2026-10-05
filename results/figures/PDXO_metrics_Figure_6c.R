####################################
# Figure 6c
####################################

library(here)
here()
source(here("environment", "requirements.R"))


# Load the saved results
final_results_rep <- readRDS(here("data", "data", "PDXO", "PDXO_metrics.rds"))


# Define your classifications
sensitive_models <- c("T224", "T245", "T251")
resistant_models <- c("CTAX23", "CTAX26", "T138")

# Reshape and Summarize data to get Mean and Standard Deviation
plot_summary <- final_results_rep %>%
  pivot_longer(
    cols = c(AUC, Emax), 
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  group_by(Organoid, Metric) %>%
  summarize(
    Mean_Value = mean(Value, na.rm = TRUE),
    SD_Value = sd(Value, na.rm = TRUE), 
    .groups = "drop"
  ) %>%
  # 2. Add the Colotubin Classification
  mutate(
    Colotubin_Status = case_when(
      Organoid %in% sensitive_models ~ "Sensitive",
      Organoid %in% resistant_models ~ "Resistant",
      TRUE ~ "Other"
    )
  )

# Prepare and sort the data
plot_summary <- plot_summary %>%
  # Set the status order: Resistant first, then Sensitive
  mutate(Colotubin_Status = factor(Colotubin_Status, levels = c("Resistant", "Sensitive", "Other"))) %>%
  
  # Sort the dataframe by this status
  arrange(Colotubin_Status) %>%
  
  # Lock in the X-axis Organoid order based on the sorted dataframe
  mutate(Organoid = factor(Organoid, levels = unique(Organoid)))

# 2. Generate the plot
ggplot(plot_summary, aes(x = Organoid, y = Mean_Value, fill = Colotubin_Status)) +
  geom_col(color = "black") +
  
  # Add Error Bars (Mean ± SD)
  geom_errorbar(aes(ymin = Mean_Value - SD_Value, ymax = Mean_Value + SD_Value),
                width = 0.2, position = position_dodge(0.9)) +
  
  facet_wrap(~ Metric, scales = "free", ncol = 2) +
  
  # Legend colors and breaks
  scale_fill_manual(
    breaks = c("Resistant", "Sensitive", "Other"), 
    values = c(
      "Resistant" = "#FDE725FF",  # Viridis bright yellow
      "Sensitive" = "#7AD151FF",  # Viridis teal/green
      "Other"     = "#440154FF"   # Viridis dark purple
    )
  ) +
  
  theme_bw() +
  theme(
    legend.position = "right", 
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(
    x = "Organoid Line (Sorted by Colotubin Status)",
    y = "Parameter Value",
    fill = "Colotubin\nSignature"
  )
