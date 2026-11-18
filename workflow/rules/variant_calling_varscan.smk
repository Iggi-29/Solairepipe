###################################
## Variant calling with Varscan2 ##
###################################

rule generate_pileup_file:
    input: 
        on_target_bam_general = results_base_path + "/{sample}/alignments/{sample}_ontarget_hpos_sorted.bam",
        sorted_bam_idx_general = results_base_path + "/{sample}/alignments/{sample}_ontarget_hpos_sorted.bam.bai", ## general
    params:
        genome = genome_unc
    output:
        mpileup_file = temp(results_base_path + "/{sample}/alignments/{sample}.mpileup")
    container:
        config["images"]["samtools_singularity"]
    benchmark:
        logs_and_bmk_base_path + "/benchmarks/{sample}_mpileup_generation.bmk"
    log:
        logs_and_bmk_base_path + "/log/{sample}_mpileup_generation.log"
    threads: 
        1
    shell:"""
        samtools mpileup -B -f {params.genome} {input.on_target_bam_general} 2>> {log} > {output.mpileup_file}
        """

rule variant_calling:
    input:
        mpileup_file = results_base_path + "/{sample}/alignments/{sample}.mpileup"
    params:
        genome = genome_unc
    output:
        vcf_snps = temp(results_base_path + "/{sample}/variants/{sample}_SNV_RELAXED.vcf"), 
        vcf_indels = temp(results_base_path + "/{sample}/variants/{sample}_INDEL_RELAXED.vcf")
    container:
        config["images"]["varscan_singularity"]
    benchmark:
        logs_and_bmk_base_path + "/benchmarks/{sample}_varscan_step1.bmk"
    log:
        log_snps = logs_and_bmk_base_path + "/log/{sample}_varscan_step1_SNV.log",
        log_inde = logs_and_bmk_base_path + "/log/{sample}_varscan_step1_INDEL.log"
    threads: 
        2
    shell:"""
        ## SNPs
        java -jar /opt/varscan/VarScan.jar mpileup2snp \
        {input.mpileup_file} \
        --min-coverage 4 --min-reads2 2 \
        --min-avg-qual 10 --min-var-freq 0.01 \
        --min-freq-for-hom 0.9 --p-value 0.1 \
        --strand-filter 0 -output-vcf > {output.vcf_snps} 2>> {log.log_snps}
        ## Indels
        java -jar /opt/varscan/VarScan.jar mpileup2indel \
        {input.mpileup_file} \
        --min-coverage 4 --min-reads2 2 \
        --min-avg-qual 10 --min-var-freq 0.01 \
        --min-freq-for-hom 0.9 --p-value 0.1 \
        --strand-filter 0 -output-vcf > {output.vcf_indels} 2>> {log.log_inde}
        """

rule join_and_sort_vcf:
    input:
        vcf_snps = results_base_path + "/{sample}/variants/{sample}_SNV_RELAXED.vcf",
        vcf_indels = results_base_path + "/{sample}/variants/{sample}_INDEL_RELAXED.vcf"
    output:
        vcf_gz_snp = results_base_path + "/{sample}/variants/{sample}_SNV_RELAXED.vcf.gz",
        vcf_gz_indel = results_base_path + "/{sample}/variants/{sample}_INDEL_RELAXED.vcf.gz"
    container:
        config["images"]["histlisb_singularity"]
    threads:
        1
    benchmark:
        logs_and_bmk_base_path + "/benchmarks/{sample}_prepare_sorting.bmk"
    log:
        logs_and_bmk_base_path + "/log/{sample}_prepare_sorting.log"
    shell:"""
        /htslib-1.16/bgzip -c {input.vcf_snps} > {output.vcf_gz_snp} 2>>{log}
        /htslib-1.16/tabix -fp vcf {output.vcf_gz_snp} 2>>{log}
        /htslib-1.16/bgzip -c {input.vcf_indels} > {output.vcf_gz_indel} 2>>{log}
        /htslib-1.16/tabix -fp vcf {output.vcf_gz_indel} 2>>{log}
        """

rule finish_varscan_vcalling:
    input:
        vcf_gz_snp = results_base_path + "/{sample}/variants/{sample}_SNV_RELAXED.vcf.gz",
        vcf_gz_indel= results_base_path + "/{sample}/variants/{sample}_INDEL_RELAXED.vcf.gz"
    output:
        vcf_out_gz = results_base_path + "/{sample}/variants/{sample}_variants_RELAXED.vcf",
        tmp_txt=temp(results_base_path + "/Run/{sample}.txt")
    params:
        sample = "{sample}"
    container:
        config["images"]["bcftools_singularity"]
    threads:
        2
    benchmark:
        logs_and_bmk_base_path + "/benchmarks/{sample}_vcf_sorting.bmk"
    log:
        logs_and_bmk_base_path + "/log/{sample}_vcf_sorting.log"
    shell:"""
        echo -e 'Sample1\t'{params.sample} > {output.tmp_txt}
        bcftools concat -a {input.vcf_gz_snp} {input.vcf_gz_indel} | bcftools sort | bcftools reheader -s {output.tmp_txt} -o {output.vcf_out_gz} 2>> {log}
        """

# ### Tigger
rule trigger_vcf_varscan_generation:
    input:
        expand(results_base_path + "/{sample}/variants/{sample}_variants_RELAXED.vcf",
        results_base_path = results_base_path, sample = samples),

 


        
        