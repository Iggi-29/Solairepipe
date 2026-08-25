#############################
#### Prepare config file ####
#############################

#### WD and libraries
library(yaml)
library(readr)

### Set the run
run <- "Run1Dx" # OPCIO en el futur

### Directories and raw data - FIX
sample_info_dir <- "/imppc/labs/eclab/ijarne/0_Clinica/999_proves_pipeline/Solairepipe/config" # OPCIO en el futur
raw_data_dir <- "/imppc/labs/eclab/ijarne/0_Clinica/999_proves_pipeline/Raw_data" # OPCIO en el futur
results_base_path <- "/imppc/labs/eclab/ijarne/0_Clinica/999_proves_pipeline/Results" # OPCIO en el futur
singulatity_images_base_path <- "/imppc/labs/eclab/Resources/Singularity_sif/Alicorn" # OPCIO en el futur
genome <- "/imppc/labs/eclab/ijarne/0_Clinica/999_proves_pipeline/Resources/data/GRCh38.primary_assembly.genome.fa.gz" # OPCIO en el futur
genome_gtf <- "/imppc/labs/eclab/ijarne/0_Clinica/999_proves_pipeline/Resources/data/gencode.v50.primary_assembly.annotation.gtf.gz" # OPCIO en el futur

### Directories and raw data - DYNAMIC
genome_unc <- gsub(pattern = "\\.gz$", replacement = "", x = genome)
genome_gtf_unc <- gsub(pattern = "\\.gz$", replacement = "", x = genome_gtf)
raw_data_dir <- paste0(raw_data_dir,"/",run)


### Import the basic yaml path
general_yaml_path <- "/imppc/labs/eclab/ijarne/0_Clinica/999_proves_pipeline/Solairepipe/config/general_config.yaml" # OPCIO en el futur
general_yaml <- yaml::read_yaml(general_yaml_path)

### Import the run data
sample_info <- list.files(path = sample_info_dir, pattern = "\\.tsv$", full.names = TRUE)
sample_info <- sample_info[grep(pattern = paste0(run,".tsv$"), x = sample_info)]
sample_info <- readr::read_tsv(file = sample_info, col_types = "c")

### Transform the sample data into a list
raw_files <- list.files(path = raw_data_dir, pattern = "\\.fastq\\.gz", full.names = TRUE, recursive = TRUE)

## Add Read files
sample_info$R1 <- ""
sample_info$R2 <- ""
for (i in 1:nrow(sample_info)) {
  sample_now <- sample_info$sample_name[i]
  reads <- grep(x = raw_files, pattern = sample_now, value = TRUE)
  sample_info$R1[i] <- grep(x = reads, pattern = "_1\\.fastq\\.gz$", value = TRUE)
  sample_info$R2[i] <- grep(x = reads, pattern = "_2\\.fastq\\.gz$", value = TRUE)
  }

### Do the yaml transformation
sample_info_yaml <- lapply(seq_len(nrow(sample_info)), function(i) {
  list(
    R1  = sample_info$R1[i],
    R2  = sample_info$R2[i],
    sex = sample_info$sex[i],
    Library = sample_info$Library[i],
    genes = sample_info$genes[i],
    HPOs = sample_info$HPOs[i],
    MUC1 = sample_info$MUC1[i]
  )
})

names(sample_info_yaml) <- sample_info$sample_name

### Add the information of the yaml file
## Paths
general_yaml$results_base_path <- paste0(results_base_path,"/",run)
general_yaml$raw_data_base_path <- raw_data_dir
general_yaml$logs_and_bmk_base_path <- paste0(paste0(results_base_path,"/",run,"/Run"))

## Ressources
general_yaml$genome <- genome
general_yaml$gtf <- genome_gtf
general_yaml$genome_unc <- genome_unc
general_yaml$gtf_unc <- genome_gtf_unc

## Images
general_yaml$images <- sapply(X = general_yaml$images, FUN = function(x) {
  full_path_name <- paste0(x = c(singulatity_images_base_path,"/",x), collapse = "")
  return(full_path_name)
}, USE.NAMES = FALSE)

## Generate the final yaml information
final_yaml <- general_yaml 
final_yaml$samples <-sample_info_yaml 


yaml::write_yaml(
  x = final_yaml,
  file = paste0("/imppc/labs/eclab/ijarne/0_Clinica/999_proves_pipeline/Solairepipe/config/",run,".yaml"))




