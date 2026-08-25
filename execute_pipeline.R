#############################
#### Prepare config file ####
#############################

#### WD and libraries
library(yaml)
library(readr)

### Directories and raw data - FIX
results_base_path <- "/imppc/labs/eclab/ijarne/0_Clinica/999_proves_pipeline/Results"
yaml_base_path <- "/imppc/labs/eclab/ijarne/0_Clinica/999_proves_pipeline/Solairepipe/config"
base_path <- "/imppc/labs/eclab/ijarne/0_Clinica/999_proves_pipeline/Solairepipe"

### Set the run
run <- "Run1Dx" # OPCIO en el futur
yaml_run <- paste0(yaml_base_path,"/",run,".yaml")
yaml_run <- yaml::read_yaml(file = yaml_run)

### Results folder to generate
results_to_generate <- yaml_run$samples
results_folders <- c()
for (i in 1:length(yaml_run$samples)) {
  result_now <- yaml_run$samples[[i]]$Library
  results_folders <- c(results_folders,result_now)
  }
results_folders <- unique(results_folders)

### Generate the Results folder and prepare the analysis
for (i in 1:length(results_folders)) {
  results_folders_now <- results_folders[i]
  results_folders_now <- paste0(results_base_path,"/",results_folders_now,"/",run)
  cmd1 <- paste0("mkdir -p ",results_folders_now,"/Pipeline")
  system(cmd1)
  
  ### Do the copies - folders
  folders_to_copy <- list.dirs(path = base_path, recursive = TRUE, full.names = TRUE)
  folders_to_copy <- folders_to_copy[-grep(pattern = "\\.git\\/", x = folders_to_copy)]
  folders_to_copy <- folders_to_copy[-grep(pattern = "\\.git", x = folders_to_copy)]
  folders_to_copy <- folders_to_copy[normalizePath(folders_to_copy) != normalizePath(base_path)]
  for (f in 1:length(x = folders_to_copy)) {
    folders_to_copy_now <- folders_to_copy[f]
    if (!grepl(pattern = "config$", x = folders_to_copy_now)) {
      cmd2 <- paste0("cp -r ",folders_to_copy_now," ",results_folders_now,"/Pipeline")
      system(cmd2)
    } else if (grepl(pattern = "config$", x = folders_to_copy_now)) {
      all_configs <- list.files(path = folders_to_copy_now, pattern = "\\.yaml$", full.names = TRUE)
      config_now <- grep(x = all_configs, pattern = paste0(run,"\\.yaml$"), value = TRUE)
      cmd3 <- paste0("mkdir -p ",results_folders_now,"/Pipeline/config")
      system(cmd3)
      cmd4 <- paste0("cp ",config_now," ",results_folders_now,"/Pipeline/config")
      system(cmd4)}
  }
  ### Do the copies - files
  files_to_copy_now <- list.files(path = base_path, full.names = TRUE, recursive = FALSE, include.dirs = FALSE, 
                                  pattern = ".sh|LICENSE|README.md|Snakefile|Execute_the_code")
  for (filee in 1:length(files_to_copy_now)) {
    cmd5 <- paste0("cp ",files_to_copy_now[filee]," ",results_folders_now)
    system(cmd5)
  }
  
  }


# for (i in 1:length(folders_to_copy))



