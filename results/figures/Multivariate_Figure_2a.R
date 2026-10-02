# Code to reproduce Figure 2a
# Forest plot of Hazard ratio: multivariate model
# --------------------------
library(here)
here()
source(here("environment", "requirements.R"))


# Inputpath
cox_inputpath <- here("data", "processed_data", "cox", "cox_cvrts_yaccs.rds")
train_inputpath <- here("data", "processed_data", "cox", "train_all.rds")

# Load data
cox <- readRDS(cox_inputpath)
train <- readRDS(train_inputpath)

# Plot
ggforest(
  model = cox,
  data = train,
  cpositions = c(0.01, NA, NA),
  fontsize = 0.8
)