############################################################################################
#GSE39582
############################################################################################

# In the course of the project the gene expression data was updated with brainarray instead of jetset 
# Therefore in the following the previous scripts workflow/preprocess/GSEXXX.R are updated with brainarray


## Gene Expression ##

### Download the tar files and save them in data/extdata/cohorts
#Due to space limitations on Github wasnt saved here

# Define your target directory for this dataset
dataset_dir <- here("data", "extdata", "cohorts","GSE39582")

# # Create the directory if it doesn't exist yet
# dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)

# Download all supplementary files (this automatically fetches GSE39582_RAW.tar)
# baseDir places the GSE39582 subfolder inside your specified path
getGEOSuppFiles("GSE39582", baseDir = here("data", "extdata", "cohorts","GSE39582"))

# Point to the downloaded tar file
tar_file <- file.path(dataset_dir, "GSE39582_RAW.tar")

# Untar the CEL files into the directory
untar(tar_file, exdir = dataset_dir)

###########################################

# Processed affymetrix arrays with the custom CDF from BrainArray

#Download CDF file from Brainarray, see details in OneNote
#Then install it

# Install from local file
# install.packages("path/hgu133plus2hsensgcdf_25.0.0.tar.gz",
#                  repos = NULL, type = "source")

# GSE39582 is based on the Affymetrix Human Genome U133 Plus 2.0 Array (HG-U133_Plus_2), which has a pre-built CDF package available via Bioconductor.

ls("package:hgu133plus2hsensgcdf")[1:5]  # List a few probe sets from the package

# cdfname	-used to specify the name of an alternative cdf package. If set to NULL, then the usual cdf package based on Affymetrix's mappings will be used.

# Read all CEL files in the directory
cel_data<-ReadAffy(cdfname = "hgu133plus2hsensgcdf",celfile.path = dataset_dir)

#Check
cel_data@cdfName

# Summarize the data
summary(cel_data)

#Normalize the data (RMA method is common) 

eset <- rma(cel_data)

#5.Inspect the data

# View expression values
exprs(eset)[1:5, 1:5]

# Sample names
sampleNames(eset)

# Probe set names
featureNames(eset)


GSE39582<- exprs(eset)

#The first 62 rows in your expression matrix — 
#like "AFFX-BioB-3_at" to "AFFX-TrpnX-M_at" — are Affymetrix control probes, not actual biological gene measurements.

# Keep only rows with names starting with "ENSG"
GSE39582 <- GSE39582[grepl("^ENSG", rownames(GSE39582)), ]

# Remove everything after the dot
colnames(GSE39582) <- sub("\\..*", "", colnames(GSE39582))

#Same for rownames 
rownames(GSE39582) <- sub("_.*", "", rownames(GSE39582))

GSE39582<-as.data.frame(GSE39582)

colnames(GSE39582) <- sub("_.*", "", colnames(GSE39582))


#### Clinical ######



#################### Stage in GSE39582 ############################

GSE39582_clinical <- read_rds(here("data", "pp", "geo", "GSE39582_clinical.rds"))

GSE39582_clinical$Stage_II<-"0"

GSE39582_clinical <- GSE39582_clinical %>%
  mutate(Stage_II = if_else(Stage_III == 0 & Stage_IV == 0, 1, as.numeric(Stage_II)))

#Reorder
GSE39582_clinical<- GSE39582_clinical[, c("age_at_diagnosis", "vital_status", "SurvTime", "gender_male", "Stage_II", "Stage_III", "Stage_IV", "MSI.Status_MSI-L/MSS")]

saveRDS(GSE39582_clinical, file = here("data", "pp", "geo_2026", "GSE39582_clinical.rds"))

######################################

#Get rownames for clinical data
GSE39582_counts <- readRDS(here("data", "pp", "geo_final", "GSE39582_counts.rds"))
rownames(GSE39582_clinical)<-colnames(GSE39582_counts)

#Remove columns 
keep<-rownames(GSE39582_clinical)
GSE39582<-GSE39582[,colnames(GSE39582) %in% keep]

keep<-colnames(GSE39582)
GSE39582_clinical<-GSE39582_clinical[rownames(GSE39582_clinical) %in% keep,]

# Check order of both datasets
all(rownames(GSE39582_clinical) %in% colnames(GSE39582))
all(rownames(GSE39582_clinical) == colnames(GSE39582))

#Order predicted_gene alphabetically
GSE39582_clinical<-GSE39582_clinical[order(row.names(GSE39582_clinical)), ]
GSE39582<-GSE39582[,order(colnames(GSE39582))]

# Check order of both datasets
all(rownames(GSE39582_clinical) %in% colnames(GSE39582))
all(rownames(GSE39582_clinical) == colnames(GSE39582))

rownames(GSE39582_clinical)<-NULL

#Save
saveRDS(GSE39582, file = here("data", "pp", "geo", "GSE39582_counts.rds"))
saveRDS(GSE39582_clinical, file = here("data", "pp", "geo", "GSE39582_clinical.rds"))

############################################################################################
#GSE17537
############################################################################################
## Gene Expression ##
### Download the tar files and save them in data/extdata/cohorts
#Due to space limitations on Github wasnt saved here

# Define your target directory for this dataset
dataset_dir <- here("data", "extdata", "cohorts","GSE17537")

# # Create the directory if it doesn't exist yet
# dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)

# Download all supplementary files (this automatically fetches GSE17537_RAW.tar)
# baseDir places the GSE17537 subfolder inside your specified path
getGEOSuppFiles("GSE17537", baseDir = here("data", "extdata", "cohorts","GSE17537"))

# Point to the downloaded tar file
tar_file <- file.path(dataset_dir, "GSE17537_RAW.tar")

# Untar the CEL files into the directory
untar(tar_file, exdir = dataset_dir)


###########################################

# Processed affymetrix arrays with the custom CDF from BrainArray

#Download CDF file from Brainarray, see details in OneNote
#Then install it


# Install from local file
# install.packages("/mnt/CCBdata/users/lea/BRAF_Postdoc_Project/Datasets/Microarray/hgu133plus2hsensgcdf_25.0.0.tar.gz", 
#                  repos = NULL, type = "source")

library(hgu133plus2hsensgcdf)


# GSE17537 is based on the Affymetrix Human Genome U133 Plus 2.0 Array (HG-U133_Plus_2), which has a pre-built CDF package available via Bioconductor.
# 1. Load the package and check if it loads without error
# 2. Check if the CDF environment is available
# 3. Try using it in a ReadAffy call (if you have CEL files)

ls("package:hgu133plus2hsensgcdf")[1:5]  # List a few probe sets from the package

# cdfname	-used to specify the name of an alternative cdf package. If set to NULL, then the usual cdf package based on Affymetrix's mappings will be used.

# Read all CEL files in the directory
cel_data<-ReadAffy(cdfname = "hgu133plus2hsensgcdf",celfile.path =dataset_dir))

#Check

cel_data@cdfName

# Summarize the data
summary(cel_data)

#4.Normalize the data (RMA method is common) 

eset <- rma(cel_data)

#5.Inspect the data

# View expression values
exprs(eset)[1:5, 1:5]

# Sample names
sampleNames(eset)

# Probe set names
featureNames(eset)


GSE17537<- exprs(eset)

#The first 62 rows in your expression matrix — 
#like "AFFX-BioB-3_at" to "AFFX-TrpnX-M_at" — are Affymetrix control probes, not actual biological gene measurements.


# Keep only rows with names starting with "ENSG"
GSE17537 <- GSE17537[grepl("^ENSG", rownames(GSE17537)), ]

# Remove everything after the dot
colnames(GSE17537) <- sub("\\..*", "", colnames(GSE17537))

#Same for rownames 
rownames(GSE17537) <- sub("_.*", "", rownames(GSE17537))

GSE17537<-as.data.frame(GSE17537)

colnames(GSE17537) <- sub("_.*", "", colnames(GSE17537))


#load clinical data

GSE17537_clinical <- readRDS(here("data", "pp", "geo_final_ENSEMBL", "GSE17537_clinical.rds"))
GSE17537_counts <- readRDS(here("data", "pp", "geo_final", "GSE17537_counts.rds"))
rownames(GSE17537_clinical)<-colnames(GSE17537_counts)


#Remove columns 
keep<-rownames(GSE17537_clinical)
GSE17537<-GSE17537[,colnames(GSE17537) %in% keep]

keep<-colnames(GSE17537)
GSE17537_clinical<-GSE17537_clinical[rownames(GSE17537_clinical) %in% keep,]

# Check order of both datasets
all(rownames(GSE17537_clinical) %in% colnames(GSE17537))
all(rownames(GSE17537_clinical) == colnames(GSE17537))

#Order predicted_gene alphabetically
GSE17537_clinical<-GSE17537_clinical[order(row.names(GSE17537_clinical)), ]
GSE17537<-GSE17537[,order(colnames(GSE17537))]

# Check order of both datasets
all(rownames(GSE17537_clinical) %in% colnames(GSE17537))
all(rownames(GSE17537_clinical) == colnames(GSE17537))

rownames(GSE17537_clinical)<-NULL

#Save
saveRDS(GSE17537, file = here("data", "pp", "geo", "GSE17537_counts.rds"))
saveRDS(GSE17537_clinical, file = here("data", "pp", "geo", "GSE17537_clinical.rds"))



############################################################################################
#GSE17536
############################################################################################
## Gene Expression ##

### Download the tar files and save them in data/extdata/cohorts
#Due to space limitations on Github wasnt saved here

# Define your target directory for this dataset
dataset_dir <- here("data", "extdata", "cohorts","GSE17536")

# # Create the directory if it doesn't exist yet
# dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)

# Download all supplementary files (this automatically fetches GSE17536_RAW.tar)
# baseDir places the GSE17536 subfolder inside your specified path
getGEOSuppFiles("GSE17536", baseDir = here("data", "extdata", "cohorts","GSE17536"))

# Point to the downloaded tar file
tar_file <- file.path(dataset_dir, "GSE17536_RAW.tar")

# Untar the CEL files into the directory
untar(tar_file, exdir = dataset_dir)


###########################################

# Processed affymetrix arrays with the custom CDF from BrainArray

#Download CDF file from Brainarray, see details in OneNote
#Then install it


# Install from local file
# install.packages("/mnt/CCBdata/users/lea/BRAF_Postdoc_Project/Datasets/Microarray/hgu133plus2hsensgcdf_25.0.0.tar.gz", 
#                  repos = NULL, type = "source")

library(hgu133plus2hsensgcdf)


# GSE17536 is based on the Affymetrix Human Genome U133 Plus 2.0 Array (HG-U133_Plus_2), which has a pre-built CDF package available via Bioconductor.
# 1. Load the package and check if it loads without error
# 2. Check if the CDF environment is available
# 3. Try using it in a ReadAffy call (if you have CEL files)

ls("package:hgu133plus2hsensgcdf")[1:5]  # List a few probe sets from the package

# cdfname	-used to specify the name of an alternative cdf package. If set to NULL, then the usual cdf package based on Affymetrix's mappings will be used.

# Read all CEL files in the directory
cel_data<-ReadAffy(cdfname = "hgu133plus2hsensgcdf",celfile.path =dataset_dir)

#Check

cel_data@cdfName

# Summarize the data
summary(cel_data)

#4.Normalize the data (RMA method is common) 

eset <- rma(cel_data)

#5.Inspect the data

# View expression values
exprs(eset)[1:5, 1:5]

# Sample names
sampleNames(eset)

# Probe set names
featureNames(eset)


GSE17536<- exprs(eset)

#The first 62 rows in your expression matrix — 
#like "AFFX-BioB-3_at" to "AFFX-TrpnX-M_at" — are Affymetrix control probes, not actual biological gene measurements.


# Keep only rows with names starting with "ENSG"
GSE17536 <- GSE17536[grepl("^ENSG", rownames(GSE17536)), ]

# Remove everything after the dot
colnames(GSE17536) <- sub("\\..*", "", colnames(GSE17536))

#Same for rownames 
rownames(GSE17536) <- sub("_.*", "", rownames(GSE17536))

GSE17536<-as.data.frame(GSE17536)

colnames(GSE17536) <- sub("_.*", "", colnames(GSE17536))


#load clinical data
GSE17536_clinical <- readRDS(here("data", "pp", "geo_final_ENSEMBL", "GSE17536_clinical.rds"))
GSE17536_counts <- readRDS(here("data", "pp", "geo_final", "GSE17536_counts.rds"))
rownames(GSE17536_clinical)<-colnames(GSE17536_counts)


#Remove columns 
keep<-rownames(GSE17536_clinical)
GSE17536<-GSE17536[,colnames(GSE17536) %in% keep]

keep<-colnames(GSE17536)
GSE17536_clinical<-GSE17536_clinical[rownames(GSE17536_clinical) %in% keep,]

# Check order of both datasets
all(rownames(GSE17536_clinical) %in% colnames(GSE17536))
all(rownames(GSE17536_clinical) == colnames(GSE17536))

#Order predicted_gene alphabetically
GSE17536_clinical<-GSE17536_clinical[order(row.names(GSE17536_clinical)), ]
GSE17536<-GSE17536[,order(colnames(GSE17536))]

# Check order of both datasets
all(rownames(GSE17536_clinical) %in% colnames(GSE17536))
all(rownames(GSE17536_clinical) == colnames(GSE17536))

rownames(GSE17536_clinical)<-NULL

#Save
saveRDS(GSE17536, file = here("data", "pp", "geo_2026", "GSE17536_counts.rds"))
saveRDS(GSE17536_clinical, file = here("data", "pp", "geo_2026", "GSE17536_clinical.rds"))

############################################################################################
#GSE29621
############################################################################################
## Gene Expression ##

### Download the tar files and save them in data/extdata/cohorts
#Due to space limitations on Github wasnt saved here

# Define your target directory for this dataset
dataset_dir <- here("data", "extdata", "cohorts","GSE29621")

# # Create the directory if it doesn't exist yet
# dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)

# Download all supplementary files (this automatically fetches GSE29621_RAW.tar)
# baseDir places the GSE29621 subfolder inside your specified path
getGEOSuppFiles("GSE29621", baseDir = here("data", "extdata", "cohorts","GSE29621"))

# Point to the downloaded tar file
tar_file <- file.path(dataset_dir, "GSE29621_RAW.tar")

# Untar the CEL files into the directory
untar(tar_file, exdir = dataset_dir)


###########################################

# Processed affymetrix arrays with the custom CDF from BrainArray

#Download CDF file from Brainarray, see details in OneNote
#Then install it


# Install from local file
# install.packages("/mnt/CCBdata/users/lea/BRAF_Postdoc_Project/Datasets/Microarray/hgu133plus2hsensgcdf_25.0.0.tar.gz", 
#                  repos = NULL, type = "source")

library(hgu133plus2hsensgcdf)


# GSE29621 is based on the Affymetrix Human Genome U133 Plus 2.0 Array (HG-U133_Plus_2), which has a pre-built CDF package available via Bioconductor.
# 1. Load the package and check if it loads without error
# 2. Check if the CDF environment is available
# 3. Try using it in a ReadAffy call (if you have CEL files)

ls("package:hgu133plus2hsensgcdf")[1:5]  # List a few probe sets from the package

# cdfname	-used to specify the name of an alternative cdf package. If set to NULL, then the usual cdf package based on Affymetrix's mappings will be used.

# Read all CEL files in the directory
cel_data<-ReadAffy(cdfname = "hgu133plus2hsensgcdf",celfile.path =dataset_dir)

#Check

cel_data@cdfName

# Summarize the data
summary(cel_data)

#4.Normalize the data (RMA method is common) 

eset <- rma(cel_data)

#5.Inspect the data

# View expression values
exprs(eset)[1:5, 1:5]

# Sample names
sampleNames(eset)

# Probe set names
featureNames(eset)


GSE29621<- exprs(eset)

#The first 62 rows in your expression matrix — 
#like "AFFX-BioB-3_at" to "AFFX-TrpnX-M_at" — are Affymetrix control probes, not actual biological gene measurements.


# Keep only rows with names starting with "ENSG"
GSE29621 <- GSE29621[grepl("^ENSG", rownames(GSE29621)), ]

# Remove everything after the dot
colnames(GSE29621) <- sub("\\..*", "", colnames(GSE29621))

#Same for rownames 
rownames(GSE29621) <- sub("_.*", "", rownames(GSE29621))

GSE29621<-as.data.frame(GSE29621)

colnames(GSE29621) <- sub("_.*", "", colnames(GSE29621))


#load clinical data

GSE29621_clinical <- readRDS(here("data", "pp", "geo_final_ENSEMBL", "GSE29621_clinical.rds"))
GSE29621_counts <- readRDS(here("data", "pp", "geo_final", "GSE29621_counts.rds"))
rownames(GSE29621_clinical)<-colnames(GSE29621_counts)


#Remove columns 
keep<-rownames(GSE29621_clinical)
GSE29621<-GSE29621[,colnames(GSE29621) %in% keep]

keep<-colnames(GSE29621)
GSE29621_clinical<-GSE29621_clinical[rownames(GSE29621_clinical) %in% keep,]

# Check order of both datasets
all(rownames(GSE29621_clinical) %in% colnames(GSE29621))
all(rownames(GSE29621_clinical) == colnames(GSE29621))

#Order predicted_gene alphabetically
GSE29621_clinical<-GSE29621_clinical[order(row.names(GSE29621_clinical)), ]
GSE29621<-GSE29621[,order(colnames(GSE29621))]

# Check order of both datasets
all(rownames(GSE29621_clinical) %in% colnames(GSE29621))
all(rownames(GSE29621_clinical) == colnames(GSE29621))

rownames(GSE29621_clinical)<-NULL

#Save
saveRDS(GSE29621, file = here("data", "pp", "geo_2026", "GSE29621_counts.rds"))
saveRDS(GSE29621_clinical, file = here("data", "pp", "geo_2026", "GSE29621_clinical.rds"))


