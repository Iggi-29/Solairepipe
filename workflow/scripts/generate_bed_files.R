#########################################################
### Generate bed file from config sample infomrmation ###
#########################################################

#### libraries 
library(dplyr)
library(readr)
library(GenomicRanges)
library(IRanges)
library(yaml)

#### sample information
# sample_now <- "25074" ## OJO - params
# sample_now_data_save_folder <- "/imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/results/25074/QC/Coverage" ## OJO - params
sample_now <- snakemake@params[["sample"]]
sample_now_data_save_folder <- snakemake@params[["saving_place"]]

#### Generate folder for results
cmd_generate <- paste0(x = c("mkdir -p ",sample_now_data_save_folder), collapse = "")
system(cmd_generate)

#### Get the files
### importation
# config_samples <- yaml::read_yaml(file = "/imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/pipeline/config/config.yaml") ## OJO - params
# hpo_data <- readr::read_tsv("/imppc/labs/eclab/Resources/AnotationDB/HPOs/06_05_2025/genes_to_phenotype.txt") ## OJO - params
# kapa_bed_data <- readr::read_tsv("/imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/bed/intersected_beds/Exome_hg38_gencode_Chr.bed",  ## OJO - params
#                             col_names = c("chr","start","end","id","score","strand"))
# kapa_annotation <- readr::read_tsv("/imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/bed/Gencode_Exome_Annotation.tsv") ### OJO - params
## naming
config_samples = snakemake@params[["config_samples"]]
hpo_data = snakemake@params[["hpo_bbdd"]]
kapa_bed_data = snakemake@params[["kapa_bed_data"]]
kapa_annotation = snakemake@params[["kapa_annotation"]]
## import
config_samples <- yaml::read_yaml(file = config_samples)
hpo_data <- readr::read_tsv(file = hpo_data)
kapa_bed_data <- readr::read_tsv(file = kapa_bed_data,
                                col_names = c("chr","start","end","id","score","strand"))
kapa_annotation <- readr::read_tsv(file = kapa_annotation)

### format the data
## sample infomration
config_samples <- config_samples[["samples"]]
config_samples <- config_samples[[sample_now]]
# genes
genes_to_work_on <- config_samples$genes
genes_to_work_on <- strsplit(x = genes_to_work_on, split = "\\, ")[[1]]
# hpos
hpos_to_work_on <- config_samples$HPOs
hpos_to_work_on <- strsplit(x = hpos_to_work_on, split = "\\, ")[[1]]

cat("OK1 \n")
#### Do the analysis
### Hpos data filtering
hpos_not_in_hpos_data <- c()
hpo_data <- hpo_data %>% 
  dplyr::filter(hpo_id %in% hpos_to_work_on)

hpos_not_in_data <- unique(hpo_data$hpo_id)
hpos_not_in_hpos_data <- hpos_to_work_on[which(!hpos_to_work_on %in% hpos_not_in_data)]

### Do the filter
hpo_data <- hpo_data %>% 
  dplyr::select(gene_symbol) %>% 
  dplyr::distinct() %>% 
  dplyr::pull(var = gene_symbol)

#### Do the analysis
### Prepare to annotate genes not in the analysis
genes_not_in_bed_genes <- c()

genes_to_look_at <- c(hpo_data, genes_to_work_on)
genes_to_look_at <- unique(genes_to_look_at)

### Do the filter on kapa annotation
kapa_annotation <- kapa_annotation %>% 
  dplyr::filter(gene_name %in% genes_to_look_at)

### Genes to worry
genes_not_in_bed_genes <- genes_to_look_at[which(!genes_to_look_at %in% kapa_annotation$gene_name)]

### Do the filter on the bed file
kapa_bed_data <- kapa_bed_data %>% 
  dplyr::filter(id %in% kapa_annotation$id)

### Mark the non-found hpos and genes 
do_the_claim <- c(
  "HPOs, not in the hpo file:",
  hpos_not_in_hpos_data,
  "Genes not in the annotation files:",
  genes_not_in_bed_genes
  )


#### Write the files
## Not in the analysis
# writeLines(
#   text = do_the_claim, 
#   con = "/imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/results/25074/QC/Coverage/Look_at_this_genes.txt")  ## OJO - output

## On target bed file
# readr::write_tsv(
#   x = kapa_bed_data, 
#   file = "/imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/results/25074/QC/Coverage/Bed_file_HPOs.bed",
#   col_names = FALSE) ## OJO - output

writeLines(
  text = do_the_claim, 
  con = snakemake@output[["not_in_the_analysis"]])  ## OJO - output

# On target bed file
readr::write_tsv(
  x = kapa_bed_data, 
  file = snakemake@output[["bed_of_hpos"]],
  col_names = FALSE) ## OJO - output




