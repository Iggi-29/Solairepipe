######################################
## Failed regions report generation ##
######################################

rule failed_regions_bed_file_generation:
    input:
        R1 = lambda wc: config["samples"][wc.sample]["R1"],
        R2 = lambda wc: config["samples"][wc.sample]["R2"]
    params:
        script_place = "/imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/pipeline/workflow/scripts", ## MILLORAR AIXO !!
        sample = "{sample}",
        saving_place = results_base_path + "/{sample}/QC/Coverage",
        hpo_bbdd = "/imppc/labs/eclab/Resources/AnotationDB/HPOs/06_05_2025/genes_to_phenotype.txt",
        kapa_bed_data = "/imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/bed/intersected_beds/Exome_hg38_gencode_Chr.bed",
        kapa_annotation = "/imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/bed/Gencode_Exome_Annotation.tsv",
        config_samples = "/imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/pipeline/config/config.yaml" ## MILLORAR AIXO !!
    output:
        not_in_the_analysis = results_base_path + "/{sample}/QC/Coverage/Look_at_this_genes.txt",
        bed_of_hpos = results_base_path + "/{sample}/QC/Coverage/Bed_file_HPOs.bed"
    conda:
        "/imppc/labs/eclab/ijarne/0_Recerca/9_mostres_marc_exoma/pipeline/envs/solar_pipe.yaml", ## MILLORAR AIXO !!
    benchmark:
        logs_and_bmk_base_path + "/benchmarks/{sample}_generate_bed_hpo_file.bmk"
    log:
        logs_and_bmk_base_path + "/log/{sample}_generate_bed_hpo_file.log"
    threads: 1
    script: "{params.script_place}/generate_bed_files.R"


rule trigger_failed_regions:
    input:
        expand(results_base_path + "/{sample}/QC/Coverage/Look_at_this_genes.txt",
        results_base_path = results_base_path, sample = samples),
        expand(results_base_path + "/{sample}/QC/Coverage/Bed_file_HPOs.bed",
        results_base_path = results_base_path, sample =samples)

        
        