library(here)
here()
source(here("environment", "requirements.R"))

# Load data 
TCGA_test_clinical <- readRDS(here("data","data", "data_partitions", "metadata_test.rds"))
cts <- readRDS(here("data","data", "data_partitions", "counts_test.rds"))

output_dir <- here("data","data", "pp", "geo")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Prepare data
################### Counts ###########################################################

cts<-t(cts)
saveRDS(cts, file = file.path(output_dir, "TCGA_test_counts.rds"))

################### Clinical ###########################################################

TCGA_test_clinical <- TCGA_test_clinical %>%
  dplyr::select(-c(1:2)) %>%
  dplyr::mutate(
    across(c(age_at_diagnosis, SurvTime), as.numeric),
    across(c(gender_male, Stage_II, Stage_III, Stage_IV, Molecular_Subtype_noCIN), as.integer)
  ) %>%
  dplyr::rename(`MSI.Status_MSI-L/MSS` = MSI.Status_MSI.L.MSS)

#Reorder
TCGA_test_clinical<- TCGA_test_clinical[, c("age_at_diagnosis", "vital_status", "SurvTime", "gender_male","Molecular_Subtype_noCIN", "Stage_II", "Stage_III", "Stage_IV", "MSI.Status_MSI-L/MSS")]

saveRDS(TCGA_test_clinical, file = file.path(output_dir, "TCGA_test_clinical.rds"))

