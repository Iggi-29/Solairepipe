### Basic snakefile For Exome Analysis
import pandas as pd
from shutil import copyfile
import os, glob, re, sys

### config file importation
configfile: "config/config.yaml"

### type of panel
type_of_panel = "Exome"

## get variables
raw_data_base_path = config["raw_data_base_path"]
results_base_path = config["results_base_path"]
logs_and_bmk_base_path = config["logs_and_bmk_base_path"]

## get the sample names
samples = [s for s in config["samples"]]

## get the bed files
panel_files = config["bed_files"][type_of_panel]
on_target_bed_file = panel_files["on_target_bed_file"]
cnv_bed_file = panel_files["cnv_bed_file"]
annotated_panel_regions = panel_files["annotated_panel_regions"]

## reference files
genome = config["genome"]
genome_unc = config["genome_unc"]

### get the rules
include: "workflow/rules/qc_analysis.smk"
include: "workflow/rules/bam_generation.smk"
include: "workflow/rules/failed_regions.smk"
include: "workflow/rules/variant_calling_varscan.smk"
include: "workflow/rules/variant_calling_GATK4.smk"

rule all:
    input:
        rules.trigger_qc_analysis.input,
        rules.trigger_bam_generation.input,
        rules.trigger_vcf_varscan_generation.input,
        rules.trigger_failed_regions.input,
        rules.trigger_vcf_gatk_generation.input


# samtools view -@ 20 -b -S /imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/results/23081/alignments/23081.sam > /imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/results/23081/alignments/23081.bam 2>> /imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/results/logs/log/23081_samtools_view.log
# samtools sort -m 2G -@ 20 /imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/results/23081/alignments/23081.bam -o /imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/results/23081/alignments/23081_sorted.bam 2>> /imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/results/logs/log/23081_samtools_sort.log
# samtools index /imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/results/23081/alignments/23081_sorted.bam -o /imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/results/23081/alignments/23081_sorted.bam.bai 2>> /imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/results/logs/log/23081_samtools_index.log
# samtools flagstat /imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/results/23081/alignments/23081_sorted.bam > /imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/results/23081/alignments/23081_sorted.bam.stats 2>> /imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/results/logs/log/23081_samtools_flagstat.log
# samtools mpileup -B -f /imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/genome/GRCh38.primary_assembly.genome.fa.gz /imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/results/23081/alignments/23081_sorted.bam > /imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/results/23081/alignments/23081.mpileup 2>> /imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/results/logs/log/23081_samtools_pileup.log
