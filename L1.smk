####################################################
###yanzeqin
###2020/10/29
#######################################################

#configfile: "config.yaml"
SAMPLES = {}
with open(config["params"]["samples"], 'rb') as sinfo:
    for line in sinfo:
        parts = line.decode('utf-8').split()
        sample = parts[0]
        SAMPLES[sample] = [parts[1], parts[2]]

print(SAMPLES)
rule all:
    input:
        expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["fastp"],"{sample}_URMAP1.fq.gz"),sample=SAMPLES.keys())
		#expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["checkv"],"quality_summary.tsv"),sample=SAMPLES.keys()),
		#expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["blastn2NT"],"{sample}.CatRNAVirus_nt.blastn"),sample=SAMPLES.keys()),
		#expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT"],"{sample}.predicted_proteins.faa.selected.list"),sample=SAMPLES.keys()),
		#expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["cox1_findhost"],"{sample}.centrifuge.KSout.Ssum"),sample=SAMPLES.keys()),
		#expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["blastp2RdRp"],"{sample}.finish.combine"),sample=SAMPLES.keys()),
		#expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["bowtie2_depth"],"{sample}_virus_depth"),sample=SAMPLES.keys()),
		#expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["closespecies"],"{sample}.contigspecies_info"),sample=SAMPLES.keys()),
		#expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT"],"{sample}.virus.contig.info"),sample=SAMPLES.keys())


##########################################01.QC-unmap-fastp-soapnuke-prinseq-rmrRNA############################################
##########################################去除raw reads中的rRNA（减少后续处理的数据量###########################################
rule rmrRNA_URMAP:
    input:
        raw_fq1 = lambda wildcards: SAMPLES[wildcards.sample][0],
        raw_fq2 = lambda wildcards: SAMPLES[wildcards.sample][1]
    output:
        temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["fastp"],"{sample}_URMAP1.fq.gz")),
        temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["fastp"],"{sample}_URMAP2.fq.gz")),
        os.path.join(config["output"]["relative"], "{sample}", config["output"]["stats"], "raw_1_seqkit_stats.txt"),
        os.path.join(config["output"]["relative"], "{sample}", config["output"]["stats"], "raw_2_seqkit_stats.txt")
    shell:
        '''
        /hwfssz5/ST_INFECTION/GlobalDatabase/user/fengqikai/software/URMAP/bin/urmap \
        -map2 {input.raw_fq1} -reverse {input.raw_fq2} \
        -ufi /hwfssz5/ST_INFECTION/GlobalDatabase/share/public_database/rRNA/URMAP_index/rRNA.dna.8.ufi \
        -samout  /dev/stdout |/hwfssz5/ST_INFECTION/GlobalDatabase/user/fengqikai/software/.conda/envs/Trinity-2.11.0/bin/samtools fastq -f 13 \
        -1 {output[0]} -2 {output[1]} --threads 2
        seqkit stats -j 2 -a {input.raw_fq1} > {output[2]}
        seqkit stats -j 2 -a {input.raw_fq2} > {output[3]}
        '''