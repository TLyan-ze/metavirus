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
    input:expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["checkv"],"quality_summary.tsv"),sample=SAMPLES.keys()),\
	expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT"],"{sample}.virus.contig.info"),sample=SAMPLES.keys()),\
	expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["blastn2NT"],"{sample}.CatRNAVirus_nt.blastn"),sample=SAMPLES.keys()),\
	expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT"],"{sample}.predicted_proteins.faa.selected.list"),sample=SAMPLES.keys()),\
	expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["cox1_findhost"],"{sample}.centrifuge.KSout.Psum"),sample=SAMPLES.keys()),\
	expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["bowtie2_depth"],"{sample}_virus_depth"),sample=SAMPLES.keys()),\
	expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["closespecies"],"{sample}.contigspecies_info"),sample=SAMPLES.keys()),\
	expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["hs_blastn"],"{sample}.genome_hsblastn2viralContigs.blastn"),sample=SAMPLES.keys())



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

######################################################fastp质控#############################################
rule QC_fastp:
	input:
		URMAP_rmrRNA_fq1 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["fastp"], "{sample}_URMAP1.fq.gz"),
		URMAP_rmrRNA_fq2 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["fastp"], "{sample}_URMAP2.fq.gz")
	output:
		fastp_fq1 = temp(os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_1.fq.gz")),
		fastp_fq2 = temp(os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_2.fq.gz")),
		fastp_html = os.path.join(config["output"]["relative"],"{sample}", config["output"]["fastp"], "{sample}.report.html"),
		fastp_json = os.path.join(config["output"]["relative"],"{sample}", config["output"]["fastp"], "{sample}.report.json"),
		URMAP_rmrRNA_fq1_stat = protected(os.path.join(config["output"]["relative"], "{sample}", config["output"]["stats"], "URMAP_1_seqkit_stats.txt")),
		URMAP_rmrRNA_fq2_stat = protected(os.path.join(config["output"]["relative"], "{sample}", config["output"]["stats"], "URMAP_2_seqkit_stats.txt")),
		fastp_fq1_stat = protected(os.path.join(config["output"]["relative"], "{sample}", config["output"]["stats"], "fastp_1_seqkit_stats.txt")),
		fastp_fq2_stat = protected(os.path.join(config["output"]["relative"], "{sample}", config["output"]["stats"], "fastp_2_seqkit_stats.txt"))
	params:
		fastp = config["softwares"]["fastp"],
		fastp_threads = config["fastp"]["threads"],
		f_adapter = config["fastp"]["f_adapter"],
		r_adapter = config["fastp"]["r_adapter"],
		quality = config["fastp"]["q"],
		un_qualified = config["fastp"]["u"],
		length_required = config["fastp"]["length_required"]		
	log:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "logs/{sample}.fastp.logs") 
	benchmark:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "benchmarks/{sample}.fastp.benchmark.txt")
	shell:
		'''
		{params.fastp} -w {params.fastp_threads} \
		-i {input.URMAP_rmrRNA_fq1} -I {input.URMAP_rmrRNA_fq2} \
		--adapter_sequence {params.f_adapter} \
		--adapter_sequence_r2 {params.r_adapter} \
		--detect_adapter_for_pe \
		-q {params.quality} \
		--length_required {params.length_required} \
		-n 2 -y -c -p \
		--disable_trim_poly_g \
		-o {output.fastp_fq1} -O {output.fastp_fq2} \
		-h {output.fastp_html} -j {output.fastp_json}
		seqkit stats -j 2 -a  {input.URMAP_rmrRNA_fq1} > {output.URMAP_rmrRNA_fq1_stat}
		seqkit stats -j 2 -a  {input.URMAP_rmrRNA_fq1} > {output.URMAP_rmrRNA_fq2_stat}
		seqkit stats -j 2 -a  {output.fastp_fq1} > {output.fastp_fq1_stat}
		seqkit stats -j 2 -a  {output.fastp_fq2} > {output.fastp_fq2_stat}
		'''

################################################
rule QC_SOAPnuke_V2:
	input:
		fastp_fq1 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["fastp"],"{sample}.fastp_1.fq.gz"),
		fastp_fq2 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["fastp"],"{sample}.fastp_2.fq.gz")
	output:
		soapnuke_fq1 = temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["soapnuke"],"{sample}_1.fq.gz")),
		soapnuke_fq2 = temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["soapnuke"],"{sample}_2.fq.gz")),
		soapnuke_fq1_stat = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"soapnuke_1_seqkit_stats.txt"),
		soapnuke_fq2_stat = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"soapnuke_2_seqkit_stats.txt")
	log:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["soapnuke"], "logs/{sample}.soapnuke.logs") 
	benchmark:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["soapnuke"], "benchmarks/{sample}.soapnuke.benchmark.txt")
	params:
		soapnuke_V2 = config["softwares"]["soapnuke_V2"],
		f_adapter = config["soapnuke_V2"]["f_adapter"],
		r_adapter = config["soapnuke_V2"]["r_adapter"],
		outdir= os.path.join(config["output"]["relative"],"{sample}",config["output"]["soapnuke"]),
		cg="{sample}_1.fq.gz",
		dg="{sample}_2.fq.gz",
		c = config["soapnuke_V2"]["c"],
		duplications = os.path.join(config["output"]["relative"], "{sample}", config["output"]["soapnuke"], "dupReads.*.gz")
	shell:
		"""
		{params.soapnuke_V2} filter -T 4 \
		-1 {input.fastp_fq1} \
		-2 {input.fastp_fq2} \
		-f {params.f_adapter} \
		-r {params.r_adapter} \
		-C {params.cg} -D {params.dg} \
		-c {params.c} \
		-l 20 \
		-q 0.2 \
		-n 0.02 \
		-4 50 \
		-o {params.outdir} 
		seqkit stats -j 2 -a  {output.soapnuke_fq1} > {output.soapnuke_fq1_stat}
		seqkit stats -j 2 -a  {output.soapnuke_fq2} > {output.soapnuke_fq2_stat}
		rm -f {params.duplications}
		"""


rule QC_prinseq:
	input:
		soapnuke_fq1 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["soapnuke"],"{sample}_1.fq.gz"),
		soapnuke_fq2 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["soapnuke"],"{sample}_2.fq.gz")
	output:
		prinseq_fq1 = temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["prinseq"],"{sample}_good_out_R1.fastq.gz")),
		prinseq_fq2 = temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["prinseq"],"{sample}_good_out_R2.fastq.gz")),
		prinseq_fq1_stat = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"prinseq_1_seqkit_stats.txt"),
		prinseq_fq2_stat = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"prinseq_2_seqkit_stats.txt")
	log:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["prinseq"], "logs/{sample}.prinseq.logs") 
	benchmark:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["prinseq"], "benchmarks/{sample}.prinseq.benchmark.txt")
	params:
		prinseq = config["softwares"]["prinseq"]
	shell:
		'''
		{params.prinseq} \
		-threads 4 \
		-fastq {input.soapnuke_fq1} -fastq2 {input.soapnuke_fq2} \
		-lc_entropy=0.5 \
		-lc_dust=0.5 \
		-out_gz \
		-out_good {output.prinseq_fq1} \
		-out_good2 {output.prinseq_fq2} \
		-out_single /dev/null -out_single2 /dev/null -out_bad /dev/null -out_bad2 /dev/null 
		seqkit stats -j 2 -a  {output.prinseq_fq1} > {output.prinseq_fq1_stat}
		seqkit stats -j 2 -a  {output.prinseq_fq1} > {output.prinseq_fq2_stat}
		'''



rule QC_SortMeRNA_V4:
	input:
		prinseq_fq1 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["prinseq"],"{sample}_good_out_R1.fastq.gz"),
		prinseq_fq2 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["prinseq"],"{sample}_good_out_R2.fastq.gz")
	output:
		sortmerna_fq1 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rmrRNA_fwd.fq.gz"),
		sortmerna_fq2 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rmrRNA_rev.fq.gz"),
		sortmerna_fq1_stat = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"sortmerna_1_seqkit_stats.txt"),
		sortmerna_fq2_stat = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"sortmerna_2_seqkit_stats.txt"),
		sortmerna_rRNA_fq1 = temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rRNA_fwd.fq.gz")),
		sortmerna_rRNA_fq2 = temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rRNA_rev.fq.gz")),
	log:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["rmrRNA"], "logs/{sample}.sortmerna.logs") 
	benchmark:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["rmrRNA"], "benchmarks/{sample}.sortmerna.benchmark.txt")
	params:
		sortmerna_V4 = config["softwares"]["sortmerna_V4"],
		rRNA = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rRNA"),
		rmrRNA = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rmrRNA"),
		kvdb_dir = os.path.join(config["output"]["relative"], "{sample}", config["output"]["rmrRNA"], "kvdb"),
		readb_dir = os.path.join(config["output"]["relative"], "{sample}", config["output"]["rmrRNA"], "readb"),
		outdir = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"]),
		rfam_5p8s = config["database"]["SortMeRNA"]["rfam_5p8s"],
		rfam_5s = config["database"]["SortMeRNA"]["rfam_5s"],
		silva_arc_16s = config["database"]["SortMeRNA"]["silva_arc_16s"],
		silva_arc_23s = config["database"]["SortMeRNA"]["silva_arc_23s"],
		silva_bac_16s = config["database"]["SortMeRNA"]["silva_bac_16s"],
		silva_bac_23s = config["database"]["SortMeRNA"]["silva_bac_23s"],
		silva_euk_18s = config["database"]["SortMeRNA"]["silva_euk_18s"],
		silva_euk_28s = config["database"]["SortMeRNA"]["silva_euk_28s"],
		indexdir = config["database"]["SortMeRNA"]["indexdir"]
	shell:
		"""
		if [ -d {params.kvdb_dir} ];then
			rm -fr {params.kvdb_dir}
		fi
		if [ -d {params.readb_dir} ];then
			rm -fr {params.readb_dir}
		fi

		{params.sortmerna_V4} \
		--ref {params.rfam_5p8s} \
		--ref {params.rfam_5s} \
		--ref {params.silva_arc_16s} \
		--ref {params.silva_arc_23s} \
		--ref {params.silva_bac_16s} \
		--ref {params.silva_bac_23s} \
		--ref {params.silva_euk_18s} \
		--ref {params.silva_euk_28s} \
		--reads {input.prinseq_fq1} --reads {input.prinseq_fq2} \
		--workdir {params.outdir} --idx-dir {params.indexdir} \
		--fastx --aligned {params.rRNA} --other {params.rmrRNA}  \
		--no-best --num_alignments 1 --paired_out  --out2 --zip-out --threads 4
		seqkit stats -j 2 -a  {output.sortmerna_fq1} > {output.sortmerna_fq1_stat}
		seqkit stats -j 2 -a  {output.sortmerna_fq2} > {output.sortmerna_fq2_stat}
		rm -fr {params.kvdb_dir} {params.readb_dir}
		"""

rule QC_Check:
	input:
		sortmerna_fq1_stat = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"sortmerna_1_seqkit_stats.txt"),
		sortmerna_fq2_stat = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"sortmerna_2_seqkit_stats.txt")
	output:
		check_QC = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"check_QC.txt"),
		rmrRNA_done = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"rmrRNA_done.txt"),
	shell:
		"""
		seqkit=/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/bin/seqkit

		pid_1=$(cat {input.sortmerna_fq1_stat}|wc -l)
		pid_2=$(cat {input.sortmerna_fq2_stat}|wc -l)

		if [ -f {output.check_QC} ];
		then
			rm -f {output.check_QC}
		fi

		touch {output.rmrRNA_done}
		#判断一端是否完整，如果完整，将1输出到output1，如果不完整就删除output1，输出到output0
		if [ $pid_1 -eq 2 ];
		then
			echo -e "1" >> {output.rmrRNA_done}
		else
			echo -e "1" >> {output.check_QC}
			rm -f {output.rmrRNA_done}
		fi

		if [ $pid_2 -eq 2 ];
		then
			echo -e "2" >> {output.rmrRNA_done}
		else
			echo -e "2" >> {output.check_QC}
			rm -f {output.rmrRNA_done}
		fi	

		if  [[ $pid_1 -eq 2 ]] $$ [[ $pid_2 -eq 2 ]];
		then
			touch {output.rmrRNA_done}
		else
			rm -f {output.rmrRNA_done}
		fi
		touch {output.check_QC}
		echo end at:`date`
		"""


##############################################02.assemble-megahit/trinity/metaspades##############################################
rule assemble_megahit:
	input:
		sortmerna_fq1 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rmrRNA_fwd.fq.gz"),
		sortmerna_fq2 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rmrRNA_rev.fq.gz")
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit"],"{sample}.contigs.fa"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"megahit_seqkit_stats.txt"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"megahit_seqkit_fx2tab.txt")
	log:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["megahit"], "logs/{sample}.megahit.logs") 
	benchmark:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["megahit"], "benchmarks/{sample}.megahit.benchmark.txt")
	params:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit"]),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit"],"intermediate_contigs"),
		"{sample}"
	shell:
		'''
		/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/bin/megahit \
		-1 {input.sortmerna_fq1} -2 {input.sortmerna_fq2} -o {params[0]} --out-prefix {params[2]} -f
		seqkit stats -j 2 -a  {output[0]} > {output[1]}
		seqkit fx2tab -j 2 -l -n -i {output[0]} > {output[2]}
		rm -rf {params[1]}
		'''

###############################################鉴定宿主来源###################################################
rule hs_blastn:
	input:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit"],"{sample}.contigs.fa")
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["hs_blastn"],"{sample}.genome_hsblastn2viralContigs.blastn")
	log:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["hs_blastn"], "logs/{sample}.hs_blastn.logs") 
	benchmark:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["hs_blastn"], "benchmarks/{sample}.hs_blastn.benchmark.txt")
	params:
		dbtemp=os.path.join(config["output"]["relative"],"{sample}",config["output"]["hs_blastn"],"db-temp")
	shell:
		'''
		mkdir -p {params.dbtemp}
		/hwfssz5/ST_INFECTION/GlobalDatabase/user/fengqikai/software/HS-BLASTN/queries-master/Linux-amd64/bin/hs-blastn -num_threads 4 -evalue 1e-5 \
		-db_dir {params.dbtemp} \
		-keep_db -block_size 500 -outfmt 6 /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/project/wuhan/host/ref_list \
		{input[0]} > {output[0]}
		rm -rf {params.dbtemp}
		'''



rule cox1_findhost:
	input:
		sortmerna_fq1 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rmrRNA_fwd.fq.gz"),
		sortmerna_fq2 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rmrRNA_rev.fq.gz")
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["cox1_findhost"],"{sample}.centrifuge.KSout.Psum")
	log:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["cox1_findhost"], "logs/{sample}.cox1_findhost.logs") 
	benchmark:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["cox1_findhost"], "benchmarks/{sample}.cox1_findhost.benchmark.txt")
	params:
		lupath = config["output"]["relative"],
		luname = "{sample}",
		luname2 = "{sample}.centrifuge.KSout",
		cox1 =  os.path.join(config["output"]["relative"],"{sample}",config["output"]["cox1_findhost"])
	shell:
		'''
		sh /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/cox1_centrifuge.sh {params.lupath} {params.luname} {params.cox1} {input.sortmerna_fq1} {input.sortmerna_fq2}
		#得到每个级别最高的
		sh /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/cox1_top_abundant_taxa.sh {params.cox1} {params.luname} {params.cox1}
		#得到蝙蝠
		bash /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/cox1_searchbylineages.sh -t Chiroptera -s {params.luname2} -l S -i {params.cox1} -o {params.cox1}
		#得到老鼠
		bash /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/cox1_searchbylineages.sh -t Rodentia -s {params.luname2} -l S -i {params.cox1} -o {params.cox1}
		#得到丰度比较高的，在种级别
		bash /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/cox1_summary.sh -s {params.luname2} -l S -i {params.cox1} -o {params.cox1}
		#在脊椎动物
		bash /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/cox1_searchbylineages.sh -t Vertebrata -s {params.luname2} -l S -i {params.cox1} -o {params.cox1}
		#在门级别
		bash /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/cox1_summary.sh -s {params.luname2} -l P -i {params.cox1} -o {params.cox1}
		'''

###################################################机器学习识别病毒########################################

################################################初筛，优化资源消耗########################################
rule cdhitest:
	input:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit"],"{sample}.contigs.fa")
	output:
		temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["cdhitest"],"{sample}.unique.contigs.fasta"))
	benchmark:
		os.path.join(config["output"]["relative"],"{sample}/benchmark/down/{sample}.cdhitest.benchmark.txt")
	params:
		l2="{sample}"
	log:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["cdhitest"], "logs/{sample}.cdhitest.logs") 
	benchmark:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["cdhitest"], "benchmarks/{sample}.cdhitest.benchmark.txt")
	shell:
		'''
		/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/envs.multi-user/cd-hit-4.8.1/bin/cd-hit-est \
		-i {input} -o {output} -M 10240 -T 5 -c 0.95 -aS 0.9 -g 1  2> {log}
		'''



rule Blastx_Virus_Refseq_NR_IMGVR:
	input:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit"],"{sample}.contigs.fa")
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["BlastxVirusNR"],"{sample}.blast"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["BlastxVirusNR"],"{sample}.id.txt"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["BlastxVirusNR"],"{sample}.blast.fasta")
	log:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["BlastxVirusNR"], "logs/{sample}.blastx.logs") 
	benchmark:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["BlastxVirusNR"], "benchmarks/{sample}.blastx.benchmark.txt")
	shell:
		'''
		/hwfssz5/ST_INFECTION/GlobalDatabase/user/fengqikai/software/diamond blastx \
		-d /hwfssz5/ST_INFECTION/GlobalDatabase/share/public_database/tblastx_screen_database/Refseq_nr_IMGVR.virus.faa \
		-q {input[0]} -k 1 -f 6 -e 0.00001 -o {output[0]}  2> {log}
		
		cat {output[0]} | awk '{{print $1}}' > {output[1]}
		
		seqkit grep -f {output[1]} {input[0]} > {output[2]}
		'''



##############################################contigs注释##############################################
rule CAT_nr_refseq:
	input:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["BlastxVirusNR"],"{sample}.blast.fasta")
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT"],"{sample}.contig2classification.txt"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT"],"{sample}.ORF2LCA.txt"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT"],"{sample}.predicted_proteins.faa"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT"],"{sample}.contig.tax"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT"],"{sample}.contig.official.tax"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT"],"{sample}.ORF.tax"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT"],"{sample}.summary.txt"),
		temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT"],"{sample}.alignment.diamond"))
	log:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["CAT"], "logs/{sample}.CAT.logs") 
	benchmark:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["CAT"], "benchmarks/{sample}.CAT.benchmark.txt")
	params:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT"],"{sample}")
	threads:4
	shell:
		'''
		/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_pack/CAT contigs \
		-n 4 --force --block_size 5.0 --index_chunks 1 -c {input} -o {params[0]} --no_stars \
		-d /hwfssz5/ST_INFECTION/GlobalDatabase/share/Database/nr_database/CAT_nr_append \
		-t /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_tax/ \
		--path_to_prodigal /hwfssz5/ST_INFECTION/GlobalDatabase/user/fengqikai/software/.conda/envs/prokka-1.14.6/bin/prodigal \
		--path_to_diamond /hwfssz5/ST_INFECTION/GlobalDatabase/user/fengqikai/software/diamond

		/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_pack/CAT add_names -i {output[0]} -o {output[4]} \
		-t /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_tax/ --only_official
		/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_pack/CAT add_names -i {output[0]} -o {output[3]} \
		-t /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_tax/
		/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_pack/CAT add_names -i {output[1]} -o {output[5]} \
		-t /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_tax/ --only_official

		/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_pack/CAT summarise \
		-c {input[0]} -i {output[4]} -o {output[6]}  2> {log}
		'''


rule RNA_RdRp_Virus:
	input:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["BlastxVirusNR"],"{sample}.blast.fasta"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT"],"{sample}.contig.tax")
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["blastp2RdRp"],"{sample}.CatRNAVirusContigs.fasta"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["blastp2RdRp"],"{sample}.CatRNAVirusContigs.id"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT"],"{sample}.virus.contig.tax"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT"],"{sample}.virus.contig.info"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT"],"{sample}.CAT_virus_failed"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT"],"{sample}.CAT_RNA_virus_failed")
	log:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["CAT"], "logs/{sample}.CAT_virus.logs") 
	benchmark:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["CAT"], "benchmarks/{sample}.CAT_virus.benchmark.txt")
	priority: 10
	shell:
		'''
		seqkit=/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/bin/seqkit
		#获取CAT的信息
		echo start CAT_extract_visual at:`date`
		grep 'Viruses (superkingdom)' {input[1]} |awk -v OFS='\t' -v FS='\t' '{{print $0}}'|sort|uniq > {output[2]}

		touch {output[2]}

		pid_1=$(cat {output[2]}|wc -l)
		if [ $pid_1 -eq 0 ];
		then
			echo -e "virus" >> {output[4]}
		else
			python3 /hwfssz5/ST_INFECTION/GlobalDatabase/share/wangyifei_handover/genome_type/CAT_info_simplify.py \
			/hwfssz5/ST_INFECTION/GlobalDatabase/share/wangyifei_handover/genome_type/ICTV_genometype.index {output[2]} {output[3]}
		fi

		echo end CAT_extract_visual at:`date`
		#获取CAT注释为RNA病毒的contigs
		grep -v DNA |sort|uniq > {output[1]}
		touch {output[1]}

		pid_2=$(cat {output[1]}|wc -l)
		if [ $pid_2 -eq 0 ];
		then
			echo -e "RNA_virus" >> {output[5]}
		else
			$seqkit grep -f {output[1]} {input[0]} > {output[0]}
		fi

		touch {output[2]}
		touch {output[0]}
		touch {output[1]}
		touch {output[5]}
		touch {output[4]}
		echo end RNA_RdRp.fasta at:`date`
		'''


rule RNA_Virus_Checkv:
	input:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["BlastxVirusNR"],"{sample}.blast.fasta")
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["checkv"],"quality_summary.tsv") 
	log:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["checkv"], "logs/{sample}.checkv.logs") 
	benchmark:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["checkv"], "benchmarks/{sample}.checkv.benchmark.txt") 
	params:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["checkv"])
	priority: 10
	shell:
		'''
		/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/envs.multi-user/CheckV/0.8.1/bin/checkv end_to_end {input[0]} {params[0]} \
		-d /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/envs.multi-user/CheckV/0.8.1/update-2021-07-03-checkv-db-v1.0
		'''


########################################rdrp比对###########################################################################


rule Blastp_rdrp:
	input:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT"],"{sample}.predicted_proteins.faa")
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["blastp2RdRp"],"{sample}.CAT_prodigal_proteins_blastp2RdRp"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT"],"{sample}.predicted_proteins.faa.selected.list")
	log:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["blastp2RdRp"], "logs/{sample}.blastp.logs") 
	benchmark:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["blastp2RdRp"], "benchmarks/{sample}.blastp.benchmark.txt")
	params:
		diamond = config["softwares"]["diamond"],
		RdRp_index = config["database"]["RdRp_index"],
		threads = config["diamond_blastp2RdRp"]["threads"],
		evalue = config["diamond_blastp2RdRp"]["evalue"],
		max_target_seqs = config["diamond_blastp2RdRp"]["max_target_seqs"],
		block_size = config["diamond_blastp2RdRp"]["block_size"],
		index_chunks = config["diamond_blastp2RdRp"]["index_chunks"],
		outfmt = config["diamond_blastp2RdRp"]["outfmt"],
		fasta = "{sample}.predicted_proteins.faa",
		indir = os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT"])
	shell:
		'''
		{params.diamond} blastp \
		--query {input[0]} \
		--db {params.RdRp_index} \
		--sensitive --max-target-seqs {params.max_target_seqs} \
		--evalue {params.evalue} --threads {params.threads} \
		--block-size {params.block_size} \
		--index-chunks {params.index_chunks} \
		--out {output[0]} \
		--outfmt {params.outfmt}
		##这个是psi比较
		sh /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/domainblast.Orthornavirae.sh \
		{params.indir} {params.indir} {params.fasta}
		'''


######################################################NT去假阳#################################################
rule Blastn_nt:
	input:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["BlastxVirusNR"],"{sample}.blast.fasta")
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["blastn2NT"],"{sample}.CatRNAVirus_nt.blastn")
	log:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["blastn2NT"], "logs/{sample}.blastn.logs") 
	benchmark:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["blastn2NT"], "benchmarks/{sample}.blastn.benchmark.txt")
	params:
		blastn = config["softwares"]["blastn"],
		seqkit = config["softwares"]["seqkit"],
		threads = config["blastx2NR_ViralContigs_blastn2NT"]["threads"],
		NT_db = config["database"]["NT_db"],
		BLASTDB_NT = config["database"]["NT_db_liqian"],
		evalue = config["blastx2NR_ViralContigs_blastn2NT"]["evalue"],
		outfmt = '"' + config["blastx2NR_ViralContigs_blastn2NT"]["outfmt"] + '"',
		max_target_seqs = config["blastx2NR_ViralContigs_blastn2NT"]["max_target_seqs"]
	benchmark:
		os.path.join(config["output"]["relative"],"{sample}/benchmark/down/{sample}.blastn_nt.benchmark.txt")
	shell:
		'''
		export BLASTDB={params.BLASTDB_NT}
		{params.blastn} -num_threads {params.threads} -query {input[0]} \
		-db {params.NT_db} -out {output[0]} -evalue {params.evalue} \
		-max_target_seqs {params.max_target_seqs} -outfmt {params.outfmt}
		'''


#############################################近缘物种共线性分析#####################################################################
rule Tblastx_closespecies:
	input:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["BlastxVirusNR"],"{sample}.blast.fasta")
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["closespecies"],"{sample}.contigspecies_info")
	log:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["closespecies"], "logs/{sample}.tblastx.logs") 
	benchmark:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["closespecies"], "benchmarks/{sample}.tblastx.benchmark.txt")
	params:
		lupath = config["output"]["relative"],
		luname = "{sample}",
		closes =  os.path.join(config["output"]["relative"],"{sample}",config["output"]["closespecies"])
	shell:
		'''
		sh /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/get_close.species_news.sh {input[0]}  {params.luname} {params.closes}
		'''

########################################mapping back&&&coverage################################################################
rule Bowtie2_RNAVirus_coverage:
	input:
		fasta = os.path.join(config["output"]["relative"],"{sample}",config["output"]["BlastxVirusNR"],"{sample}.blast.fasta"),
		sortmerna_fq1 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rmrRNA_fwd.fq.gz"),
		sortmerna_fq2 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rmrRNA_rev.fq.gz")
	output:
		temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["BlastxVirusNR"],"{sample}_virus.sam")),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["bowtie2_depth"],"{sample}_virus_sofa.txt"),
		temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["bowtie2_depth"],"{sample}_virus_sorted.bam")),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["bowtie2_depth"],"{sample}_virus_bam.stat"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["bowtie2_depth"],"{sample}_virus_sort.bam"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["bowtie2_depth"],"{sample}_virus_depth")
	log:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["bowtie2_depth"], "logs/{sample}.bowtie2.logs") 
	benchmark:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["bowtie2_depth"], "benchmarks/{sample}.bowtie2.benchmark.txt")
	params:
		index= os.path.join(config["output"]["relative"],"{sample}",config["output"]["BlastxVirusNR"],"{sample}_spades")
	shell:
		'''
		samtools=/hwfssz5/ST_INFECTION/GlobalDatabase/user/fengqikai/software/.conda/envs/Trinity-2.11.0/bin/samtools
		bowtie2-build {input.fasta}  {params.index}
		bowtie2 -p 8 -x {params.index} -1 {input.sortmerna_fq1} -2 {input.sortmerna_fq2} -S {output[0]} 
		/ldfssz1/ST_INFECTION/P18Z10200N0164_Resistance_PN/User/zhaohailong/software/soap.coverage -cvg -sam -refsingle {input.fasta} -i {output[0]} -o {output[1]}
		$samtools view -@ 3 -bS {output[0]} -F 4 | $samtools sort -@ 3 -n -o {output[2]}
		$samtools sort -o {output[4]} {output[2]}
		$samtools stats -@ 3 -c 10,100000,1 -i 100000 {output[4]} > {output[3]}
		$samtools depth {output[4]} > {output[5]}
		'''