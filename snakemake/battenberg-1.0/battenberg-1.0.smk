#!/usr/bin/env snakemake

##### ATTRIBUTION #####

# Original Author:  Shaghayegh Soudi
# Contributors:    NA 


configfile:
    "/oak/stanford/groups/emoding/analysis/shaghayegh/snakemake-pipelines/battenberg-1.0/config/config.yaml"

DIRS,SAMPLES = glob_wildcards(config['data']+"{dir}/{sample}.sorted.deduped.readgroup.bam")


rule all:
    input:
        expand("results/mapped/{dir}/{sample}.sorted.bam", zip, dir=DIRS, sample=SAMPLES)


# Symlinks the input files into the module results directory (under '00-inputs/')
rule battenberg_input_bam:
    input:
        bam = config['data']+"{dir}/{sample}.sorted.deduped.readgroup.bam"
        bai = config['data']+"{dir}/{sample}.sorted.deduped.readgroup.bai"       
    output:
        bam = "00-input/bam/{dir}/{sample}.sorted.deduped.readgroup.bam",
        bai = "00-input/bam/{dir}/{sample}.sorted.deduped.readgroup.bai"  
    shell:
        """
        ln -s {input.bam} {output.bam}
        ln -s {input.bai} {output.bai}
        """


rule infer_patient_sex:
    input: 
        normal_bam = rules.battenberg_input_bam.output,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output: 
        sex_result ="00-input/infer_sex/{dir}/{sample}.sex"
    resources:
        **CFG["resources"]["infer_sex"]
    log:
        stderr = CFG["logs"]["infer_sex"] + "{seq_type}--{genome_build}/{normal_id}_infer_sex_stderr.log"
    conda:
        "envs/samtools.yaml"
    threads: 8
    shell:
        op.as_one_line(""" 
        PATH={BATTENBERG_SCRIPT_PATH}:$PATH;
        echo "running {rule} for {wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr} ;
        calc_sex_status.sh {input.normal_bam} {input.fasta} {wildcards.normal_id} > {output.sex_result} 2>> 
{log.stderr} &&
        echo "DONE running {rule} for {wildcards.normal_id} on $(hostname) at $(date)" >> {log.stderr} 
        """)


rule _run_battenberg:
    input:
        tumour_bam = rules.battenberg_input_bam.output",
        normal_bam = rules.battenberg_input_bam.output",
        sex_result = rules.infer_patient_sex.output,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        impute_info = str(rules._battenberg_get_reference.output.impute_info)

    output:
        refit=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_refit_suggestion.txt",
        sub=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_subclones.txt",
        ac=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_alleleCounts.tab"),
        mb=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantBAF.tab"),
        mlrg=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantLogR_gcCorrected.tab"),
        mlr=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantLogR.tab"),
        nlr=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_normalLogR.tab"),
        nb=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_normalBAF.tab"),
        cp=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_cellularity_ploidy.txt"
    log:
        stdout = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_battenberg.stdout.log",
        stderr = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_battenberg.stderr.log"
    params:
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        script = CFG["inputs"]["battenberg_script"],
        out_dir = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}",
        ref = CFG["dirs"]["inputs"] + "reference/{genome_build}"
    conda:
        CFG["conda_envs"]["battenberg"]
    resources:
        **CFG["resources"]["battenberg"]
    threads:
        CFG["threads"]["battenberg"]
    shell:
       op.as_one_line("""
        if [[ $(head -c 4 {params.fasta}) == ">chr" ]]; then chr_prefixed='true'; else chr_prefixed='false'; fi;
        echo "$chr_prefixed"
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stdout};
        sex=$(cut -f 4 {input.sex_result}| tail -n 1); 
        echo "setting sex as $sex";
        Rscript {params.script} -t {wildcards.tumour_id} 
        -n {wildcards.normal_id} --tb $(readlink -f {input.tumour_bam}) --nb $(readlink -f {input.normal_bam}) -f {input.fasta} --reference $(readlink -f {params.ref})
        -o {params.out_dir} --chr_prefixed_genome $chr_prefixed --sex $sex --cpu {threads} >> {log.stdout} 2>> {log.stderr} &&  
        echo "DONE {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" >> {log.stdout}; 
        """)

 
# Convert the subclones.txt (best fit) to igv-friendly SEG files. 
rule _battenberg_to_igv_seg:
input:
    sub = rules._run_battenberg.output.sub,
    cnv2igv = ancient(CFG["inputs"]["cnv2igv"])
output:
    seg = CFG["dirs"]["battenberg"] + 
"{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_subclones.igv.seg"
log:
    stderr = CFG["logs"]["battenberg"] + 
"{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_seg2igv.stderr.log"
threads: 1
group: "battenberg_post_process"
shell:
    op.as_one_line("""
    echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > 
{log.stderr};
    python {input.cnv2igv} --mode battenberg --sample {wildcards.tumour_id} 
    {input.sub} > {output.seg} 2>> {log.stderr}
    """)



# Fill subclones.txt with empty regions for compatibility with downstream tools
rule _battenberg_fill_subclones:
    input:
        sub = str(rules._run_battenberg.output.sub)
    output:
        sub = CFG["dirs"]["fill_regions"] + 
"{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_subclones.filled.txt"
    log:
        stderr = CFG["logs"]["fill_regions"] + 
"{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_fill_subclones.stderr.log"
    threads: 1
    group: "battenberg_post_process"
    params:
        path = config["lcr-modules"]["_shared"]["lcr-scripts"] + "fill_segments/1.0/",
        script = "fill_segments.sh",
        arm_file = lambda w: "src/chromArm.hg38.bed" if "38" in str({w.genome_build}) else 
"src/chromArm.grch37.bed",
        blacklist_file = lambda w: "src/blacklisted.hg38.bed" if "38" in str({w.genome_build}) else 
"src/blacklisted.grch37.bed"
    conda:
        CFG["conda_envs"]["bedtools"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > 
{log.stderr};
        bash {params.path}{params.script}
        {params.path}{params.arm_file}
        {input.sub}
        {params.path}{params.blacklist_file}
        {output.sub}
        {wildcards.tumour_id}
        subclones
        2>> {log.stderr}
        """)
        
        

 
#due to the large number of files (several per chromosome) that are not explicit outputs, do some glob-based 
cleaning in the output directory
rule _battenberg_cleanup:
input:
    rules._battenberg_to_igv_seg.output.seg
output:
    complete = CFG["dirs"]["battenberg"] + 
"{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_cleanup_complete.txt"
group: "battenberg_post_process"
shell:
    op.as_one_line("""
    d=$(dirname {output});
    rm -f $d/*impute_input* &&
    rm -f $d/*alleleFrequencies* &&
    rm -f $d/*aplotype* &&
    rm -f $d/*BAFsegmented* && 
    touch {output.complete}
    """)
    
    
    def _battenberg_get_chain(wildcards):
        if "38" in str({wildcards.genome_build}):
            return reference_files("genomes/{genome_build}/chains/grch38/hg38ToHg19.over.chain")
        else:
            return reference_files("genomes/{genome_build}/chains/grch37/hg19ToHg38.over.chain")   
