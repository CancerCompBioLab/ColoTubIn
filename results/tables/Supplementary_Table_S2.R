####################################
# Supplementary Table S2
####################################


library(here)
here()
source(here("environment", "requirements.R"))


# --- CONFIGURATION ---
external_cohorts_inputpath <- here("data", "data","pp", "geo")
output_dir <- here("results", "tables")

# if (!dir.exists(output_dir)) {
#   dir.create(output_dir, recursive = TRUE)
# }

cohorts <- c("GSE17536", "GSE17537", "GSE29621", "GSE39582", "EMTAB12862", "ws4_spinal")

cvrt <- c(
  "vital_status", "SurvTime", "gender_male", "age_at_diagnosis",
  "Stage_II", "Stage_III", "Stage_IV", "MSI.Status_MSI.L.MSS"
)

# --- EXTRACTION FUNCTION FOR GEO COHORTS ---
extract_clinical_data <- function(cohort) {
  file_path <- file.path(external_cohorts_inputpath, paste0(cohort, "_clinical.rds"))
  clin_data <- readRDS(file_path)
  
  names(clin_data) <- make.names(names(clin_data))
  cols_to_keep <- intersect(cvrt, names(clin_data))
  
  clin_subset <- clin_data[, cols_to_keep, drop = FALSE]
  # Keep the file prefix as ws4_spinal, but map the resulting cohort label to Spinal
  clin_subset$Cohort <- ifelse(cohort == "ws4_spinal", "SPINAL", cohort)
  
  return(clin_subset)
}

# --- LOAD TCGA CLINICAL COHORT ---
tcga_path <- here("data", "extdata", "tcga-coad","metadata_coad_patients.rds")
tcga_data <- readRDS(tcga_path)

names(tcga_data) <- make.names(names(tcga_data))
tcga_cols_to_keep <- intersect(cvrt, names(tcga_data))
tcga_subset <- tcga_data[, tcga_cols_to_keep, drop = FALSE]
tcga_subset$Cohort <- "TCGA"

# --- MERGE COHORTS ---
geo_merged <- bind_rows(lapply(cohorts, extract_clinical_data))
merged_clinical_data <- bind_rows(geo_merged, tcga_subset)

# --- DATA PREPROCESSING & 1825-DAY CENSORING ---
if("age_at_diagnosis" %in% names(merged_clinical_data)) {
  merged_clinical_data$age_at_diagnosis <- ifelse(
    merged_clinical_data$Cohort == "GSE29621", 
    NA, 
    merged_clinical_data$age_at_diagnosis
  )
}

# Convert binary gender_male (0/1) to Female/Male
if("gender_male" %in% names(merged_clinical_data)) {
  merged_clinical_data$gender <- factor(
    merged_clinical_data$gender_male,
    levels = c(0, 1),
    labels = c("Female", "Male")
  )
}

# Standardise vital status to numeric 0/1, apply 1825-day cutoff, and factorise
merged_clinical_data <- merged_clinical_data %>%
  mutate(
    vital_status_num = case_when(
      is.character(vital_status) | is.factor(vital_status) ~ ifelse(vital_status == "Alive", 0, 1),
      TRUE ~ as.numeric(vital_status)
    ),
    vital_status_censored = case_when(
      SurvTime > 1825 ~ 0, 
      TRUE ~ vital_status_num
    ),
    SurvTime = pmin(SurvTime, 1825, na.rm = TRUE)
  ) %>%
  mutate(
    vital_status = factor(vital_status_censored, levels = c(0, 1), labels = c("Alive", "Dead"))
  )

merged_clinical_data <- merged_clinical_data %>%
  mutate(stage = case_when(
    Stage_IV == 1 ~ "Stage 4",
    Stage_III == 1 ~ "Stage 3",
    Stage_II == 1 ~ "Stage 2",
    Stage_IV == 0 & Stage_III == 0 & Stage_II == 0 ~ "Stage 1",
    TRUE ~ NA_character_
  ))

# Map MSI status 0/1 to words (MSS / MSI)
msi_col <- intersect(c("MSI.Status_MSI.L.MSS", "MSI.Status_MSI-L.MSS"), names(merged_clinical_data))
if(length(msi_col) > 0) {
  merged_clinical_data$MSI_Status <- factor(
    as.numeric(merged_clinical_data[[msi_col[1]]]),
    levels = c(0, 1),
    labels = c("MSI", "MSS")
  )
}

# Set "TCGA" as the first factor level
other_cohorts <- sort(setdiff(unique(merged_clinical_data$Cohort), "TCGA"))
merged_clinical_data$Cohort <- factor(merged_clinical_data$Cohort, levels = c("TCGA", other_cohorts))

# --- CREATE SUMMARY TABLE ---
table_data <- merged_clinical_data %>%
  select(any_of(c("Cohort", "age_at_diagnosis", "SurvTime", "vital_status", "gender", "stage", "MSI_Status")))

clinical_table <- table_data %>%
  tbl_summary(
    by = Cohort,
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "ifany",
    missing_text = "Missing",
    label = list(
      age_at_diagnosis ~ "Age at Diagnosis",
      SurvTime ~ "Survival Time (Days - 5yr Capped)",
      vital_status ~ "Vital Status (at 5 yrs)",
      gender ~ "Gender",
      stage ~ "Tumor Stage",
      MSI_Status ~ "MSI Status"
    )
  ) %>%
  add_overall() %>%
  modify_header(label = "**Variable**") %>%
  bold_labels()

# Print table
clinical_table

# --- SAVE TABLE ---
output_file <- file.path(output_dir, "clinical_characteristics_1825days_all_cohorts.docx")
gt::gtsave(as_gt(clinical_table), file = output_file)

message(paste("Table successfully saved to:", output_file))
