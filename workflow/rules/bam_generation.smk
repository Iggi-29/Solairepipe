#############################################
## bam generation and alignemtn processing ##
#############################################

# rule trimming_trimmomatic:
#     input: 
#         R1 = lambda wc: config["samples"][wc.sample]["R1"],
#         R2 = lambda wc: config["samples"][wc.sample]["R2"]
#     params:
#         avqual = 15,
#         out_dir = results_base_path + "/{sample}/trimmed_reads/"
#     output:
#         trimmedR1 = results_base_path + "/{sample}/trimmed_reads/{sample}_FP.fastq.gz",
#         trimmedR2 = results_base_path + "/{sample}/trimmed_reads/{sample}_RP.fastq.gz",
#         trimmed_rubish_R1 = temp(results_base_path + "/{sample}/trimmed_reads/{sample}_FU.fastq.gz"),
#         trimmed_rubish_R2 = temp(results_base_path + "/{sample}/trimmed_reads/{sample}_RU.fastq.gz"),
#     container:
#         config["images"]["trimmomatic_singularity"]
#     benchmark:
#         logs_and_bmk_base_path + "/benchmarks/{sample}_trimming.bmk"
#     log:
#         logs_and_bmk_base_path + "/log/{sample}_trimming.log"
#     threads: 1
#     shell:"""
#         mkdir -p {params.out_dir}
#         java -jar /Trimmomatic-0.38/trimmomatic-0.38.jar PE \
#         {input.R1} {input.R2} {output.trimmedR1} {output.trimmed_rubish_R1} {output.trimmedR2} {output.trimmed_rubish_R2} \
#         AVGQUAL:{params.avqual} 2> {log}    
#         """
rule trimming_trimmomatic:
    input: 
        R1 = lambda wc: config["samples"][wc.sample]["R1"],
        R2 = lambda wc: config["samples"][wc.sample]["R2"]
    params:
        avqual = 15,
        out_dir = results_base_path + "/{sample}/trimmed_reads/"
    output:
        trimmedR1 = results_base_path + "/{sample}/trimmed_reads/{sample}_FP.fastq.gz",
        trimmedR2 = results_base_path + "/{sample}/trimmed_reads/{sample}_RP.fastq.gz",
        trimmed_rubish_R1 = temp(results_base_path + "/{sample}/trimmed_reads/{sample}_FU.fastq.gz"),
        trimmed_rubish_R2 = temp(results_base_path + "/{sample}/trimmed_reads/{sample}_RU.fastq.gz"),
    container:
        config["images"]["trimmomatic_singularity"]
    benchmark:
        logs_and_bmk_base_path + "/benchmarks/{sample}_trimming.bmk"
    log:
        logs_and_bmk_base_path + "/log/{sample}_trimming.log"
    threads: 1
    shell:"""
        mkdir -p {params.out_dir:q}

        tmpdir=$(mktemp -d "{params.out_dir}/subset.XXXXXX")
        trap 'rm -rf "$tmpdir"' EXIT

        gzip -cd {input.R1:q} | awk 'NR <= 4000000' > "$tmpdir/R1.fastq"
        gzip -cd {input.R2:q} | awk 'NR <= 4000000' > "$tmpdir/R2.fastq"

        java -jar /Trimmomatic-0.38/trimmomatic-0.38.jar PE \
            -threads {threads} \
            "$tmpdir/R1.fastq" "$tmpdir/R2.fastq" \
            {output.trimmedR1:q} {output.trimmed_rubish_R1:q} \
            {output.trimmedR2:q} {output.trimmed_rubish_R2:q} \
            AVGQUAL:{params.avqual} \
            2> {log:q}        
        """

rule bwa_alignment:
    input:
        trimmedR1 = results_base_path + "/{sample}/trimmed_reads/{sample}_FP.fastq.gz",
        trimmedR2 = results_base_path + "/{sample}/trimmed_reads/{sample}_RP.fastq.gz",
    params:
        out_dir = results_base_path + "/{sample}/alignments/",
        sample = "{sample}",
        genome = genome
    output:
        sam = temp(results_base_path + "/{sample}/alignments/{sample}.sam")
    container: config["images"]["bwa_singularity"]
    benchmark:
        logs_and_bmk_base_path + "/benchmarks/{sample}_alignment.bmk"
    log:
        logs_and_bmk_base_path + "/log/{sample}_alignment.log"
    threads: 10
    shell:"""
        mkdir -p {params.out_dir}
        bwa mem -t {threads} -R '@RG\\tID:{params.sample}\\tSM:{params.sample}\\tPL:ILLUMINA' \
        {params.genome} {input.trimmedR1} {input.trimmedR2} > {output.sam} 2>> {log}
        """

rule samtools:
    input:
        sam = results_base_path + "/{sample}/alignments/{sample}.sam"
    params:
        out_dir = results_base_path + "/{sample}/alignments/",
        genome_unc = genome_unc
    output:
        bam = temp(results_base_path + "/{sample}/alignments/{sample}.bam"),
        sorted_bam = results_base_path + "/{sample}/alignments/{sample}_sorted.bam",
        sorted_bam_idx = results_base_path + "/{sample}/alignments/{sample}_sorted.bam.bai",
        sorted_bam_stats = results_base_path + "/{sample}/alignments/{sample}_sorted.bam.stats",
    container: 
        config["images"]["samtools_singularity"]
    benchmark:
        logs_and_bmk_base_path + "/benchmarks/{sample}_samtools.bmk"
    log:
        samtools_view = logs_and_bmk_base_path + "/log/{sample}_samtools_view.log",
        samtools_sort = logs_and_bmk_base_path + "/log/{sample}_samtools_sort.log",
        samtools_index = logs_and_bmk_base_path + "/log/{sample}_samtools_index.log",
        samtools_flagstat = logs_and_bmk_base_path + "/log/{sample}_samtools_flagstat.log",
        samtools_pileup = logs_and_bmk_base_path + "/log/{sample}_samtools_pileup.log"
    threads: 10
    shell:"""
        samtools view -@ {threads} -b -S {input.sam} > {output.bam} 2>> {log.samtools_view}
        samtools sort -m 2G -@ {threads} {output.bam} -o {output.sorted_bam} 2>> {log.samtools_sort}
        samtools index {output.sorted_bam} -o {output.sorted_bam_idx} 2>> {log.samtools_index}
        samtools flagstat {output.sorted_bam} > {output.sorted_bam_stats} 2>> {log.samtools_flagstat}
        """

### On target parts
rule on_target_bams:
    input:
        sorted_bam = results_base_path + "/{sample}/alignments/{sample}_sorted.bam", ## anterior
        Bed_file_HPOs = results_base_path + "/{sample}/QC/Coverage/Bed_file_HPOs.bed", ## bed file to check
    params:
        bed_file = "/imppc/labs/eclab/Resources/Genomes/Hsapiens/bwa_hg38/bed/intersected_beds/Exome_hg38_gencode_Chr.bed"
    output:
        on_target_bam_general = results_base_path + "/{sample}/alignments/{sample}_ontarget_general_sorted.bam",
        on_target_bam_hpos = results_base_path + "/{sample}/alignments/{sample}_ontarget_hpos_sorted.bam",
    container:
        config["images"]["bedtools_singularity"]
    benchmark:
        logs_and_bmk_base_path + "/benchmarks/{sample}_on_target_bams.bmk"
    log:
        on_target_bam_general = logs_and_bmk_base_path + "/log/{sample}_on_target_bams_general.log",
        on_target_bam_hpos = logs_and_bmk_base_path + "/log/{sample}_on_target_bams_hpos.log"
    threads: 4
    shell:"""
    ## on.target with the raw bed file
    bedtools intersect -abam {input.sorted_bam} -b {params.bed_file} > {output.on_target_bam_general} 2>> {log.on_target_bam_general}
    ## on.target with the hpo sample file
    bedtools intersect -abam {input.sorted_bam} -b {input.Bed_file_HPOs} > {output.on_target_bam_hpos} 2>> {log.on_target_bam_hpos}
    """

rule samtools_ontarget:
    input:
        on_target_bam_general = results_base_path + "/{sample}/alignments/{sample}_ontarget_general_sorted.bam",
        on_target_bam_hpos = results_base_path + "/{sample}/alignments/{sample}_ontarget_hpos_sorted.bam",
    output:
        sorted_bam_idx_general = results_base_path + "/{sample}/alignments/{sample}_ontarget_general_sorted.bam.bai", ## general
        sorted_bam_stats_general = results_base_path + "/{sample}/alignments/{sample}_ontarget_general_sorted.stats", ## general
        sorted_bam_idx_hpos = results_base_path + "/{sample}/alignments/{sample}_ontarget_hpos_sorted.bam.bai", ## hpos
        sorted_bam_stats_hpos = results_base_path + "/{sample}/alignments/{sample}_ontarget_hpos_sorted.bam.stats", ## hpos
    container:
        config["images"]["samtools_singularity"]
    benchmark:
        logs_and_bmk_base_path + "/benchmarks/{sample}_samtools_ontarget.bmk"
    log:
        sorted_bam_idx_general = logs_and_bmk_base_path + "/log/{sample}_sorted_bam_idx_general.log",
        sorted_bam_stats_general = logs_and_bmk_base_path + "/log/{sample}_sorted_bam_stats_general.log",
        sorted_bam_idx_hpos = logs_and_bmk_base_path + "/log/{sample}_sorted_bam_idx_hpos.log",
        sorted_bam_stats_hpos = logs_and_bmk_base_path + "/log/{sample}_sorted_bam_stats_hpos.log"
    threads: 2
    shell:"""
    ## general bam
    samtools index {input.on_target_bam_general} -o {output.sorted_bam_idx_general} 2>> {log.sorted_bam_idx_general}
    samtools flagstat {input.on_target_bam_general} > {output.sorted_bam_stats_general} 2>> {log.sorted_bam_stats_general}
    ## hpo bam
    samtools index {input.on_target_bam_hpos} -o {output.sorted_bam_idx_hpos} 2>> {log.sorted_bam_idx_hpos}
    samtools flagstat {input.on_target_bam_hpos} > {output.sorted_bam_stats_hpos} 2>> {log.sorted_bam_stats_hpos}
    """

rule trigger_bam_generation:
    input:
        ## TRIMMING
        # expand(results_base_path + "/{sample}/trimmed_reads/{sample}_FP.fastq.gz",
        # results_base_path = results_base_path, sample = samples),
        # expand(results_base_path + "/{sample}/trimmed_reads/{sample}_RP.fastq.gz",
        # results_base_path = results_base_path, sample = samples)
        ## SAMTOOLS DESPRES DE ALIGNMENT AMB BWA
        # expand(results_base_path + "/{sample}/alignments/{sample}_sorted.bam",
        # results_base_path = results_base_path, sample = samples),
        # expand(results_base_path + "/{sample}/alignments/{sample}_sorted.bam.bai",
        # results_base_path = results_base_path, sample = samples),
        # expand(results_base_path + "/{sample}/alignments/{sample}_sorted.bam.stats",
        # results_base_path = results_base_path, sample = samples),
        ## ON TARGET BAMS
        # expand(results_base_path + "/{sample}/alignments/{sample}_ontarget_general_sorted.bam",
        # results_base_path = results_base_path, sample = samples),
        # expand(results_base_path + "/{sample}/alignments/{sample}_ontarget_hpos_sorted.bam",
        # results_base_path = results_base_path, sample = samples),
        ## ULTIMA RULE
        expand(results_base_path + "/{sample}/alignments/{sample}_ontarget_general_sorted.bam.bai",
        results_base_path = results_base_path, sample = samples), 
        expand(results_base_path + "/{sample}/alignments/{sample}_ontarget_general_sorted.stats",
        results_base_path = results_base_path, sample = samples), 
        expand(results_base_path + "/{sample}/alignments/{sample}_ontarget_hpos_sorted.bam.bai",
        results_base_path = results_base_path, sample = samples), 
        expand(results_base_path + "/{sample}/alignments/{sample}_ontarget_hpos_sorted.bam.stats",
        results_base_path = results_base_path, sample = samples), 


        
        