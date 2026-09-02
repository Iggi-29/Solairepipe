### Basic snakefile for HC analysis
import pandas as pd
from shutil import copyfile
import os, glob, re, sys

### config file importation
configfile: "config/config.yaml"

### type of panel
type_of_panel = "HC"

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
include: "workflow/rules/variant_calling.smk"

rule all:
    input:
        rules.trigger_qc_analysis.input,
        rules.trigger_bam_generation.input,
        rules.trigger_vcf_generation.input
