library(here)
here()
source(here("environment", "requirements.R"))

# Load data 

# Clinical data
load(here("data", "extdata", "cohorts", "EMTAB_12862_clinical.RData"))

# Counts
cts <- read.table(here("data", "extdata", "cohorts", "CRC.SW.mRNA.TPM.txt.gz"), sep = "\t", header = TRUE, row.names = 1)
# Prepare data

################### Counts ###########################################################

#Remove Rectum
colData<-colData[colData$Tumour.Site=="Right Colon"|colData$Tumour.Site=="Left Colon",]

#Remove columns in cts which are not rows in colData
keep<-rownames(colData)
cts<-cts[,colnames(cts) %in% keep]

#Remove columns in colData which are not rows in cts
keep<-colnames(cts)
colData<-colData[rownames(colData) %in% keep,]

#Some checks
# Check order of both datasets
all(rownames(colData) %in% colnames(cts))
all(rownames(colData) == colnames(cts))

cts<-log2(cts+1)

saveRDS(cts, file = here("data", "pp", "geo", "EMTAB12862_counts.rds"))

############################ Clinical ################################################

EMTAB_12862_clinical<-colData

# Make new variables / rename 

EMTAB_12862_clinical$age.at.diagnosis = EMTAB_12862_clinical$Age.at.diagnosis
EMTAB_12862_clinical$vital_status<-EMTAB_12862_clinical$Vital.Status
EMTAB_12862_clinical$SurvTime<-EMTAB_12862_clinical$Overall.survival.days
EMTAB_12862_clinical$gender_male=0
EMTAB_12862_clinical$gender_male[EMTAB_12862_clinical$Sex=="Male"] = 1


EMTAB_12862_clinical$Stage_II = 0
EMTAB_12862_clinical$Stage_II [EMTAB_12862_clinical$Tumour.Stage=="Stage II"] = 1

EMTAB_12862_clinical$Stage_III = 0
EMTAB_12862_clinical$Stage_III [EMTAB_12862_clinical$Tumour.Stage=="Stage III"] = 1

EMTAB_12862_clinical$Stage_IV = 0
EMTAB_12862_clinical$Stage_IV[EMTAB_12862_clinical$Tumour.Stage=="Stage IV"] = 1

EMTAB_12862_clinical$`MSI.Status_MSI-L/MSS`<-0
EMTAB_12862_clinical$`MSI.Status_MSI-L/MSS`[EMTAB_12862_clinical$MSI.Status=="MSS"] = 1

# Molecular_Subtype_noCIN not provided in EMTAB12862

EMTAB_12862_clinical$CMS_1<-0
EMTAB_12862_clinical$CMS_1[EMTAB_12862_clinical$CMS.Tumour=="CMS1"] = 1

EMTAB_12862_clinical$CMS_2<-0
EMTAB_12862_clinical$CMS_2[EMTAB_12862_clinical$CMS.Tumour=="CMS2"] = 1

EMTAB_12862_clinical$CMS_3<-0
EMTAB_12862_clinical$CMS_3[EMTAB_12862_clinical$CMS.Tumour=="CMS3"] = 1

EMTAB_12862_clinical$CMS_4<-0
EMTAB_12862_clinical$CMS_4[EMTAB_12862_clinical$CMS.Tumour=="CMS4"] = 1

EMTAB_12862_clinical$BRAFmt<-0
EMTAB_12862_clinical$BRAFmt[EMTAB_12862_clinical$mutation=="p.Val640Glu"] = 1


############################

#Filter final dataset
EMTAB_12862_clinical<-dplyr::select(EMTAB_12862_clinical,contains("age.at.diagnosis"),contains("vital_status"),contains("SurvTime"),contains("gender_male"),contains("Stage_II"),contains("Stage_III"),contains("Stage_IV"),contains("MSI.Status_MSI-L/MSS"),contains("CMS_1"),contains("CMS_2"),contains("CMS_3"),contains("CMS_4"),contains("BRAFmt"))


#Remove first column 
EMTAB_12862_clinical<-EMTAB_12862_clinical[,-1]

EMTAB_12862_clinical$age.at.diagnosis<-as.numeric(EMTAB_12862_clinical$age.at.diagnosis)
EMTAB_12862_clinical$SurvTime<-as.numeric(EMTAB_12862_clinical$SurvTime)
EMTAB_12862_clinical$gender_male<-as.integer(EMTAB_12862_clinical$gender_male)
EMTAB_12862_clinical$Stage_II<-as.integer(EMTAB_12862_clinical$Stage_II)
EMTAB_12862_clinical$Stage_III<-as.integer(EMTAB_12862_clinical$Stage_III)
EMTAB_12862_clinical$Stage_IV<-as.integer(EMTAB_12862_clinical$Stage_IV)
EMTAB_12862_clinical$`MSI.Status_MSI-L/MSS`<-as.integer(EMTAB_12862_clinical$`MSI.Status_MSI-L/MSS`)
EMTAB_12862_clinical$CMS_1<-as.integer(EMTAB_12862_clinical$CMS_1)
EMTAB_12862_clinical$CMS_2<-as.integer(EMTAB_12862_clinical$CMS_2)
EMTAB_12862_clinical$CMS_3<-as.integer(EMTAB_12862_clinical$CMS_3)
EMTAB_12862_clinical$CMS_4<-as.integer(EMTAB_12862_clinical$CMS_4)
EMTAB_12862_clinical$BRAFmt<-as.integer(EMTAB_12862_clinical$BRAFmt)
EMTAB_12862_clinical$age_at_diagnosis<-EMTAB_12862_clinical$age.at.diagnosis
EMTAB_12862_clinical<-EMTAB_12862_clinical %>% dplyr::select(-age.at.diagnosis)

#Reorder
EMTAB_12862_clinical<- EMTAB_12862_clinical[, c("age_at_diagnosis", "vital_status", "SurvTime", "gender_male", "Stage_II", "Stage_III", "Stage_IV", "MSI.Status_MSI-L/MSS","CMS_1","CMS_2","CMS_3","CMS_4","BRAFmt")]

#Some checks
# Check order of both datasets
all(rownames(colData) %in% colnames(cts))
all(rownames(colData) == colnames(cts))

# Remove row names
rownames(EMTAB_12862_clinical) <- NULL

# EMTAB_12862_clinical$vital_status_new<-NULL
saveRDS(EMTAB_12862_clinical, file = here("data", "pp", "geo", "EMTAB12862_clinical.rds"))
