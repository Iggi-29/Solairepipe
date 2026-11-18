################################
## Variant calling with GATK4 ##
################################

rule mark_duplicates:
    input: 
        on_target_bam_general = results_base_path + "/{sample}/alignments/{sample}_ontarget_hpos_sorted.bam",
        sorted_bam_idx_general = results_base_path + "/{sample}/alignments/{sample}_ontarget_hpos_sorted.bam.bai", ## general
    output:
        marked_dup_file = results_base_path + "/{sample}/alignments/gatk_data/{sample}.markedDup.bam",
    params:
        metrics = results_base_path + "/{sample}/alignments/gatk_data/{sample}.markedDup.metrics.txt",
        gatk_base_path = results_base_path + "/{sample}/alignments/gatk_data",
        tmp_dir = results_base_path + "/{sample}/alignments/gatk_data/tmp"
    container:
        "/imppc/labs/eclab/Resources/Singularity_sif/Images_Recerca/gatk_latest.sif"    
    benchmark:
        logs_and_bmk_base_path + "/benchmarks/{sample}_markdups.bmk"
    log:
        logs_and_bmk_base_path + "/log/{sample}_markdups.log"
    threads: 
        10
    shell:"""
        ## Make the directory
        mkdir -p {params.gatk_base_path}
        mkdir -p {params.tmp_dir}
        
        ## Mark the duplicates
        gatk --java-options "-Xmx8g" MarkDuplicates \
        --INPUT {input.on_target_bam_general} \
        --OUTPUT {output.marked_dup_file} \
        --METRICS_FILE {params.metrics} \
        --REMOVE_DUPLICATES false \
        --REMOVE_SEQUENCING_DUPLICATES false \
        --DUPLICATE_SCORING_STRATEGY SUM_OF_BASE_QUALITIES \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        --TAGGING_POLICY All \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY STRICT \
        --TMP_DIR {params.tmp_dir} \
        > {log} 2>&1
        """

### OJO SEQUENCE_DICTIONARY <- A POSAR DINS DEL CONFIG!
### OJO REFERENCE_SEQUENCE <- A POSAR DINS DEL CONFIG!
rule coverege_assessing:
    input: 
        marked_dup_file = results_base_path + "/{sample}/alignments/gatk_data/{sample}.markedDup.bam",
        bed_file_hpos = results_base_path + "/{sample}/QC/Coverage/Bed_file_HPOs.bed"
    output:
        hs_analysis = results_base_path + "/{sample}/QC/Coverage/gatk_data/Hs_metrics.txt",
        per_target_cov = results_base_path + "/{sample}/QC/Coverage/gatk_data/per_target_cov.txt",
        targets_list = results_base_path + "/{sample}/QC/Coverage/gatk_data/targets.interval_list"
    params:
        gatk_base_path = results_base_path + "/{sample}/QC/Coverage/gatk_data",
    container:
        "/imppc/labs/eclab/Resources/Singularity_sif/Images_Recerca/gatk_latest.sif"    
    benchmark:
        logs_and_bmk_base_path + "/benchmarks/{sample}_hs_metrics_collect.bmk"
    log:
        hs_metrics = logs_and_bmk_base_path + "/log/{sample}_hs_metrics_collect.log",
        bedto_interval = logs_and_bmk_base_path + "/log/{sample}_bed_to_interval.log"
    threads: 
        10
    shell:"""
        ## Make the directory
        mkdir -p {params.gatk_base_path}
        
        ## Bed file to interval list
        gatk --java-options "-Xmx8g" BedToIntervalList \
        --INPUT {input.bed_file_hpos} \
        --OUTPUT {output.targets_list} \
        --SEQUENCE_DICTIONARY /imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/genome/GRCh38.primary_assembly.genome_PAR_masked.dict \
        --SORT true \
        --UNIQUE true \
        > {log.bedto_interval} 2>&1


        ## Coverage measure
        gatk --java-options "-Xmx8g" CollectHsMetrics \
        --INPUT {input.marked_dup_file} \
        --OUTPUT {output.hs_analysis} \
        --REFERENCE_SEQUENCE /imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/genome/GRCh38.primary_assembly.genome_PAR_masked.fa \
        --TARGET_INTERVALS {output.targets_list} \
        --BAIT_INTERVALS {output.targets_list} \
        --METRIC_ACCUMULATION_LEVEL ALL_READS \
        --MINIMUM_MAPPING_QUALITY 20 \
        --MINIMUM_BASE_QUALITY 20 \
        --CLIP_OVERLAPPING_READS true \
        --PER_TARGET_COVERAGE {output.per_target_cov} \
        --VALIDATION_STRINGENCY STRICT \
        > {log.hs_metrics} 2>&1
        """

### OJO REFERENCE_SEQUENCE <- A POSAR DINS DEL CONFIG!
### OJO sites1 <- A POSAR DINS DEL CONFIG!
### OJO sites2 <- A POSAR DINS DEL CONFIG!
### OJO sites3 <- A POSAR DINS DEL CONFIG!
rule base_recallibrator_model:
    input: 
        marked_dup_file = results_base_path + "/{sample}/alignments/gatk_data/{sample}.markedDup.bam",
    output:
        recall_table = results_base_path + "/{sample}/alignments/gatk_data/{sample}.recallModel.table",
    params:
        reference = "/imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/genome/GRCh38.primary_assembly.genome_PAR_masked.fa",
        sites1 = "/imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/genome/basre_recalibrator/dbsnp/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
        sites2 = "/imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/genome/basre_recalibrator/indels/Homo_sapiens_assembly38.known_indels.vcf.gz",
        sites3 = "/imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/genome/basre_recalibrator/mills_and_1000G/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"        
    container:
        "/imppc/labs/eclab/Resources/Singularity_sif/Images_Recerca/gatk_latest.sif"    
    benchmark:
        logs_and_bmk_base_path + "/benchmarks/{sample}_recallModel.bmk"
    log:
        logs_and_bmk_base_path + "/log/{sample}_recallModel.log"
    threads: 
        10
    shell:"""
        gatk --java-options "-Xmx8g" BaseRecalibrator \
        --reference {params.reference} \
        --input {input.marked_dup_file} \
        --known-sites {params.sites1} \
        --known-sites {params.sites2} \
        --known-sites {params.sites3} \
        --output {output.recall_table} \
        > {log} 2>&1
        """

rule base_recallibrator_apply_model:
    input: 
        marked_dup_file = results_base_path + "/{sample}/alignments/gatk_data/{sample}.markedDup.bam",
        recall_table = results_base_path + "/{sample}/alignments/gatk_data/{sample}.recallModel.table"
    output:
        recall_bam = results_base_path + "/{sample}/alignments/gatk_data/{sample}.recallModel.bam",
    params:
        reference = "/imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/genome/GRCh38.primary_assembly.genome_PAR_masked.fa",
    container:
        "/imppc/labs/eclab/Resources/Singularity_sif/Images_Recerca/gatk_latest.sif"    
    benchmark:
        logs_and_bmk_base_path + "/benchmarks/{sample}_recallModel_apply.bmk"
    log:
        logs_and_bmk_base_path + "/log/{sample}_recallModel_apply.log"
    threads: 
        10
    shell:"""
        gatk --java-options "-Xmx8g" ApplyBQSR \
        --reference {params.reference} \
        --input {input.marked_dup_file} \
        --bqsr-recal-file {input.recall_table} \
        --output {output.recall_bam} \
        --create-output-bam-index true \
        --preserve-qscores-less-than 6 \
        --quantize-quals 0 \
        > {log} 2>&1
        """

rule base_recallibrator_model_check:
    input: 
        recall_bam = results_base_path + "/{sample}/alignments/gatk_data/{sample}.recallModel.bam",
        recall_table = results_base_path + "/{sample}/alignments/gatk_data/{sample}.recallModel.table"
    output:
        recall_table = results_base_path + "/{sample}/alignments/gatk_data/{sample}.recallModel_after_bsqr.table",
        plots_of_recallibration = results_base_path + "/{sample}/QC/Base_recallibration/gatk_data/recall_data_plot.pdf"
    params:
        reference = "/imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/genome/GRCh38.primary_assembly.genome_PAR_masked.fa",
        sites1 = "/imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/genome/basre_recalibrator/dbsnp/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
        sites2 = "/imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/genome/basre_recalibrator/indels/Homo_sapiens_assembly38.known_indels.vcf.gz",
        sites3 = "/imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/genome/basre_recalibrator/mills_and_1000G/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        plot_place = results_base_path + "/{sample}/QC/Base_recallibration/gatk_data"      
    container:
        "/imppc/labs/eclab/Resources/Singularity_sif/Images_Recerca/gatk_latest.sif"    
    benchmark:
        logs_and_bmk_base_path + "/benchmarks/{sample}_recallModel_check.bmk"
    log:
        log1 = logs_and_bmk_base_path + "/log/{sample}_recallModel_check.log",
        log2 = logs_and_bmk_base_path + "/log/{sample}_recallModel_plot.log"
    threads: 
        10
    shell:"""
        gatk --java-options "-Xmx8g" BaseRecalibrator \
        --reference {params.reference} \
        --input {input.recall_bam} \
        --known-sites {params.sites1} \
        --known-sites {params.sites2} \
        --known-sites {params.sites3} \
        --output {output.recall_table} \
        > {log.log1} 2>&1
        
        mkdir -p {params.plot_place}

        gatk --java-options "-Xmx4g" AnalyzeCovariates \
        -before {input.recall_table} \
        -after {output.recall_table} \
        -plots {output.plots_of_recallibration} \
        > {log.log2} 2>&1
        """


import shlex


# Replace with the actual path to your GRCh38 XY ploidy BED.
XY_PLOIDY_BED = "/imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/genome/ploidy_bed/GRCh38.XY.ploidy.bed"


def haplotypecaller_sample_sex(wildcards):
    sample = str(wildcards.sample)
    sex = str(
        config["samples"][sample].get("sex", "")
    ).strip().upper()

    if sex not in ("F", "M"):
        raise ValueError(
            f"Sample {sample}: expected sex 'F' or 'M', "
            f"but found {sex!r}."
        )

    return sex


def haplotypecaller_ploidy_input(wildcards):
    if haplotypecaller_sample_sex(wildcards) == "M":
        return [XY_PLOIDY_BED]

    return []


def haplotypecaller_sex_args(wildcards, input):
    sex = haplotypecaller_sample_sex(wildcards)

    if sex == "F":
        return (
            "--sample-ploidy 2 "
            "--exclude-intervals chrY "
            "--exclude-intervals chrM"
        )

    return (
        "--sample-ploidy 2 "
        f"--ploidy-regions {shlex.quote(str(input.ploidy[0]))} "
        "--exclude-intervals chrY:10001-2781479 "
        "--exclude-intervals chrY:56887903-57217415 "
        "--exclude-intervals chrM"
    )


rule haplotype_caller:
    input:
        recall_bam = results_base_path + "/{sample}/alignments/gatk_data/{sample}.recallModel.bam",
        targets_list = results_base_path + "/{sample}/QC/Coverage/gatk_data/targets.interval_list",
        ploidy = haplotypecaller_ploidy_input
    output:
        gvcf = results_base_path + "/{sample}/variants/gatk_data/{sample}.g.vcf.gz",
        tbi = results_base_path + "/{sample}/variants/gatk_data/{sample}.g.vcf.gz.tbi"
    params:
        reference = "/imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/genome/GRCh38.primary_assembly.genome_PAR_masked.fa",
        gatk_base_path = results_base_path + "/{sample}/alignments/gatk_data",
        tmp_dir = results_base_path + "/{sample}/alignments/gatk_data/tmp/haplotype_caller",
        sex_args = haplotypecaller_sex_args
    container:
        "/imppc/labs/eclab/Resources/Singularity_sif/Images_Recerca/gatk_latest.sif"
    benchmark:
        logs_and_bmk_base_path + "/benchmarks/{sample}_haplotype_caller.bmk"
    log:
        logs_and_bmk_base_path + "/log/{sample}_haplotype_caller.log"
    threads:
        10
    shell:"""
        mkdir -p {params.gatk_base_path} {params.tmp_dir}

        gatk --java-options "-Xmx8g" HaplotypeCaller \
            --reference {params.reference} \
            --input {input.recall_bam} \
            --output {output.gvcf} \
            --min-base-quality-score 10 \
            --emit-ref-confidence GVCF \
            --intervals {input.targets_list} \
            --interval-padding 100 \
            {params.sex_args} \
            --native-pair-hmm-threads {threads} \
            --pcr-indel-model CONSERVATIVE \
            --create-output-variant-index true \
            --tmp-dir {params.tmp_dir} \
            > {log} 2>&1
        """

rule genotype_gvcfs:
    input:
        gvcf = results_base_path + "/{sample}/variants/gatk_data/{sample}.g.vcf.gz",
        tbi = results_base_path + "/{sample}/variants/gatk_data/{sample}.g.vcf.gz.tbi"
    output:
        vcf= results_base_path + "/{sample}/variants/gatk_data/{sample}.raw.vcf.gz",
        tbi = results_base_path + "/{sample}/variants/gatk_data/{sample}.raw.vcf.gz.tbi"
    params:
        reference = "/imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/genome/GRCh38.primary_assembly.genome_PAR_masked.fa",
        gatk_base_path = results_base_path + "/{sample}/alignments/gatk_data",
        tmp_dir = results_base_path + "/{sample}/alignments/gatk_data/tmp/genotype_gvcfs"
    container:
        "/imppc/labs/eclab/Resources/Singularity_sif/Images_Recerca/gatk_latest.sif"
    benchmark:
        logs_and_bmk_base_path + "/benchmarks/{sample}_genotype_gvcfs.bmk"
    log:
        logs_and_bmk_base_path + "/log/{sample}_genotype_gvcfs.log"
    threads:
        1
    shell:"""
        mkdir -p {params.gatk_base_path} {params.tmp_dir}

        gatk --java-options "-Xmx8g" GenotypeGVCFs \
            --reference {params.reference} \
            --variant {input.gvcf} \
            --output {output.vcf} \
            --include-non-variant-sites false \
            --create-output-variant-index true \
            --tmp-dir {params.tmp_dir} \
            > {log} 2>&1
        """

rule filter_variant_sites:
    input: 
        vcf = results_base_path + "/{sample}/variants/gatk_data/{sample}.raw.vcf.gz",
        tbi = results_base_path + "/{sample}/variants/gatk_data/{sample}.raw.vcf.gz.tbi"
    output:
        vcf = results_base_path + "/{sample}/variants/gatk_data/{sample}.site_filtered.vcf.gz",
        tbi = results_base_path + "/{sample}/variants/gatk_data/{sample}.site_filtered.vcf.gz.tbi"
    params:
        reference = "/imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/genome/GRCh38.primary_assembly.genome_PAR_masked.fa",
        gatk_base_path = results_base_path + "/{sample}/alignments/gatk_data",
        tmp_dir = results_base_path + "/{sample}/alignments/gatk_data/tmp/filter_variant_sites"
    container:
        "/imppc/labs/eclab/Resources/Singularity_sif/Images_Recerca/gatk_latest.sif"
    benchmark:
        logs_and_bmk_base_path + "/benchmarks/{sample}_filter_variant_sites.bmk"
    log:
        logs_and_bmk_base_path + "/log/{sample}_filter_variant_sites.log"
    threads:
        1
    shell:"""
        mkdir -p {params.gatk_base_path} {params.tmp_dir}

        gatk --java-options "-Xmx8g" VariantFiltration \
            --reference {params.reference} \
            --variant {input.vcf} \
            --output {output.vcf} \
            --filter-name "QD2" \
            --filter-expression 'QD < 2.0' \
            --filter-name "QUAL30" \
            --filter-expression 'QUAL < 30.0' \
            --filter-name "SNP_MQ40" \
            --filter-expression 'vc.isSNP() && MQ < 40.0' \
            --filter-name "SNP_MQRankSum12_5" \
            --filter-expression 'vc.isSNP() && MQRankSum < -12.5' \
            --filter-name "SNP_ReadPosRankSum8" \
            --filter-expression 'vc.isSNP() && ReadPosRankSum < -8.0' \
            --filter-name "INDEL_FS200" \
            --filter-expression '(vc.isIndel() || vc.isMixed()) && FS > 200.0' \
            --filter-name "INDEL_ReadPosRankSum20" \
            --filter-expression '(vc.isIndel() || vc.isMixed()) && ReadPosRankSum < -20.0' \
            --missing-values-evaluate-as-failing false \
            --create-output-variant-index true \
            --tmp-dir {params.tmp_dir} \
            > {log} 2>&1
        """
            # --filter-name "SNP_FS60" \
            # --filter-expression 'vc.isSNP() && FS > 60.0' \
            # --filter-name "SNP_SOR3" \
            # --filter-expression 'vc.isSNP() && SOR > 3.0' \

rule filter_genotype_confidence:
    input:
        vcf = results_base_path + "/{sample}/variants/gatk_data/{sample}.site_filtered.vcf.gz",
        tbi = results_base_path + "/{sample}/variants/gatk_data/{sample}.site_filtered.vcf.gz.tbi"
    output:
        vcf = results_base_path + "/{sample}/variants/gatk_data/{sample}.genotype_filtered.vcf.gz",
        tbi = results_base_path + "/{sample}/variants/gatk_data/{sample}.genotype_filtered.vcf.gz.tbi"
    params:
        reference = "/imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/genome/GRCh38.primary_assembly.genome_PAR_masked.fa",
        gatk_base_path = results_base_path + "/{sample}/alignments/gatk_data",
        tmp_dir = results_base_path + "/{sample}/alignments/gatk_data/tmp/filter_genotype_confidence",
        min_dp=4,
        min_gq=20
    container:
        "/imppc/labs/eclab/Resources/Singularity_sif/Images_Recerca/gatk_latest.sif"
    benchmark:
        logs_and_bmk_base_path + "/benchmarks/{sample}_filter_genotype_confidence.bmk"
    log:
        logs_and_bmk_base_path + "/log/{sample}_filter_genotype_confidence.log"
    threads:
        1
    shell:"""
        mkdir -p {params.gatk_base_path:q} {params.tmp_dir:q}

        gatk --java-options "-Xmx8g" VariantFiltration \
            --reference {params.reference:q} \
            --variant {input.vcf:q} \
            --output {output.vcf:q} \
            --genotype-filter-name "LowDP" \
            --genotype-filter-expression 'DP < {params.min_dp}' \
            --genotype-filter-name "LowGQ" \
            --genotype-filter-expression 'GQ < {params.min_gq}' \
            --missing-values-evaluate-as-failing true \
            --invalidate-previous-filters false \
            --set-filtered-genotype-to-no-call false \
            --create-output-variant-index true \
            --tmp-dir {params.tmp_dir:q} \
            > {log:q} 2>&1
        """

rule finalize_variants:
    input:
        vcf = results_base_path + "/{sample}/variants/gatk_data/{sample}.genotype_filtered.vcf.gz",
        tbi = results_base_path + "/{sample}/variants/gatk_data/{sample}.genotype_filtered.vcf.gz.tbi"
    output:
        vcf = results_base_path + "/{sample}/variants/gatk_data/{sample}.final.vcf.gz",
        tbi = results_base_path + "/{sample}/variants/gatk_data/{sample}.final.vcf.gz.tbi"
    params:
        reference = "/imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/genome/GRCh38.primary_assembly.genome_PAR_masked.fa",
        tmp_dir = results_base_path + "/{sample}/alignments/gatk_data/tmp/finalize_variants"
    container:
        "/imppc/labs/eclab/Resources/Singularity_sif/Images_Recerca/gatk_latest.sif"
    benchmark:
        logs_and_bmk_base_path + "/benchmarks/{sample}_finalize_variants.bmk"
    log:
        finalize1 = logs_and_bmk_base_path + "/log/{sample}_finalize_variants_f1.log",
        finalize2 = logs_and_bmk_base_path + "/log/{sample}_finalize_variants_f2.log"
    threads:
        1
    shell:"""
        mkdir -p {params.tmp_dir}

        workdir=$(mktemp -d {params.tmp_dir}/run.XXXXXX)
        trap 'rm -rf "$workdir"' EXIT

        # First pass: mask genotypes that failed Step 11.
        gatk --java-options "-Xmx8g" SelectVariants \
            --reference {params.reference} \
            --variant {input.vcf} \
            --output "$workdir/masked.vcf.gz" \
            --set-filtered-gt-to-nocall true \
            --create-output-variant-index true \
            --tmp-dir "$workdir" \
            > {log.finalize1} 2>&1

        # Second pass: retain passing sites with called alternate genotypes.
        gatk --java-options "-Xmx8g" SelectVariants \
            --reference {params.reference} \
            --variant "$workdir/masked.vcf.gz" \
            --output {output.vcf} \
            --exclude-filtered true \
            --exclude-non-variants true \
            --max-nocall-number 0 \
            --select 'vc.getGenotype(0).isCalled() && !vc.getGenotype(0).isHomRef() && !vc.getGenotype(0).isMixed()' \
            --create-output-variant-index true \
            --tmp-dir "$workdir" \
            >> {log.finalize2} 2>&1
        """

### Tigger
rule trigger_vcf_gatk_generation:
    input:
        expand(results_base_path + "/{sample}/alignments/gatk_data/{sample}.markedDup.bam",
        results_base_path = results_base_path, sample = samples),
        expand(results_base_path + "/{sample}/QC/Coverage/gatk_data/Hs_metrics.txt",
        results_base_path = results_base_path, sample = samples),
        expand(results_base_path + "/{sample}/QC/Coverage/gatk_data/per_target_cov.txt",
        results_base_path = results_base_path, sample = samples),
        expand(results_base_path + "/{sample}/QC/Coverage/gatk_data/targets.interval_list",
        results_base_path = results_base_path, sample = samples),
        expand(results_base_path + "/{sample}/alignments/gatk_data/{sample}.recallModel.table",
        results_base_path = results_base_path, sample = samples),
        expand(results_base_path + "/{sample}/alignments/gatk_data/{sample}.recallModel.bam",
        results_base_path = results_base_path, sample = samples),
        expand(results_base_path + "/{sample}/alignments/gatk_data/{sample}.recallModel_after_bsqr.table",
        results_base_path = results_base_path, sample = samples),
        expand(results_base_path + "/{sample}/QC/Base_recallibration/gatk_data/recall_data_plot.pdf",
        results_base_path = results_base_path, sample = samples),
        expand(results_base_path + "/{sample}/variants/gatk_data/{sample}.g.vcf.gz",
        results_base_path = results_base_path, sample = samples),
        expand(results_base_path + "/{sample}/variants/gatk_data/{sample}.g.vcf.gz.tbi",
        results_base_path = results_base_path, sample = samples),
        expand(results_base_path + "/{sample}/variants/gatk_data/{sample}.raw.vcf.gz",
        results_base_path = results_base_path, sample = samples),
        expand(results_base_path + "/{sample}/variants/gatk_data/{sample}.raw.vcf.gz.tbi",
        results_base_path = results_base_path, sample = samples),
        expand(results_base_path + "/{sample}/variants/gatk_data/{sample}.site_filtered.vcf.gz",
        results_base_path = results_base_path, sample = samples),
        expand(results_base_path + "/{sample}/variants/gatk_data/{sample}.site_filtered.vcf.gz.tbi",
        results_base_path = results_base_path, sample = samples),
        expand(results_base_path + "/{sample}/variants/gatk_data/{sample}.genotype_filtered.vcf.gz",
        results_base_path = results_base_path, sample = samples),
        expand(results_base_path + "/{sample}/variants/gatk_data/{sample}.genotype_filtered.vcf.gz.tbi",
        results_base_path = results_base_path, sample = samples),
        expand(results_base_path + "/{sample}/variants/gatk_data/{sample}.final.vcf.gz",
        results_base_path = results_base_path, sample = samples),
        expand(results_base_path + "/{sample}/variants/gatk_data/{sample}.final.vcf.gz.tbi",
        results_base_path = results_base_path, sample = samples),


 


        
        