
#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Shaghayegh Soudi
# Contributors:     



configfile:
    "PATH/TO/config.yaml"

SAMPLES, = glob_wildcards(config['data']+"/{id}_bam")



# Symlinks the input files into the module results directory (under '00-inputs/')
rule _battenberg_input_bam:
    input:
        bam = expand("{sample_id}-{group}.sorted.deduped.readgroup.bam", sample_id=SAMPLES, group=)
        bai = expand("{sample_id}-{group}.sorted.deduped.readgroup.bai", sample_id=SAMPLES, group=)

    output:
        bam = "{sample_id}/bam/{sample_id}.bam",
        bai = "{sample_id}/bam/{sample_id}.bam.bai"
    params:
        output = config['data']+"/Alignment"    
    shell:
        """
        mkdir -p {wildcards.sample_id}/bam
        ln -s {input.bam} {output.bam}
        ln -s {input.bai} {output.bai}
        """


##########################
#########################
rule _infer_patient_sex:
    input: 
        normal_bam = expand("{sample_id}.final.bam", sample=SAMPLES),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output: sex_result = CFG["dirs"]["infer_sex"] + "{seq_type}--{genome_build}/{normal_id}.sex"
    resources:
        **CFG["resources"]["infer_sex"]
    log:
        stderr = CFG["logs"]["infer_sex"] + "{seq_type}--{genome_build}/{normal_id}_infer_sex_stderr.log"
    conda:
        CFG["conda_envs"]["samtools"]
    group: "setup_run"
    threads: 8
    shell:
        op.as_one_line(""" 
        PATH={BATTENBERG_SCRIPT_PATH}:$PATH;
        echo "running {rule} for {wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr} ;
        calc_sex_status.sh {input.normal_bam} {input.fasta} {wildcards.normal_id} > {output.sex_result} 2>> 
{log.stderr} &&
        echo "DONE running {rule} for {wildcards.normal_id} on $(hostname) at $(date)" >> {log.stderr} 
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
