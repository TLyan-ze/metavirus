#####################################################
###yanzeqin
###2020/10/29
#######################################################

#configfile: "config.yaml"

SAMPLES = {}
with open(config['params']['samples'], 'rb') as f:
    for line in f:
        parts = line.decode('utf-8').split()
        key = parts[0]
        SAMPLES[key] = [[parts[1], parts[2]]

rule all:
	input:
		expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_CatRNAVirus_checkv"],"quality_summary.tsv"),sample=SAMPLES.keys()),
		expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatRNAVirus_nt.blastn"),sample=SAMPLES.keys()),
		expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"{sample}_virus_sofa.txt"),sample=SAMPLES.keys()),
		expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["cox1_findhost"],"{sample}.centrifuge.KSout.Ssum"),sample=SAMPLES.keys()),
		expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_rdrp"],"{sample}.finish.combine"),sample=SAMPLES.keys()),
		#expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit"],"final.virfinder"),sample=SAMPLES.keys()),
		expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["closespecies"],"{sample}.contigspecies_info"),sample=SAMPLES.keys()),
		expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.blast"),sample=SAMPLES.keys())



##########################################01.QC-unmap-fastp-soapnuke-prinseq-rmrRNA############################################
##########################################去除raw reads中的rRNA（减少后续处理的数据量###########################################
rule rmrRNA_URMAP:
	input:
		raw_fq1 = lambda wildcards: Samples[wildcards.sample][0],
		raw_fq2 = lambda wildcards: Samples[wildcards.sample][1]
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
		URMAP_rmrRNA_fq1 = temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["fastp"], "{sample}_URMAP1.fq.gz")),
		URMAP_rmrRNA_fq2 = temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["fastp"], "{sample}_URMAP2.fq.gz"))
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
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["soapnuke"],"{sample}.logs")
	benchmark:
		os.path.join(config["output"]["relative"],"{sample}/benchmark/QC/{sample}.soapnuke.benchmark.txt")
	params:
		soapnuke_V2 = config["softwares"]["soapnuke_V2"],
		f_adapter = config["soapnuke_V2"]["f_adapter"],
		r_adapter = config["soapnuke_V2"]["r_adapter"],
		configfile = config["soapnuke_V2"]["configfile"],
		rg= os.path.join(config["output"]["relative"],"{sample}",config["output"]["soapnuke"]),
		cg="{sample}_1.fq.gz",
		dg="{sample}_2.fq.gz",
		c = config["soapnuke_V2"]["c"],
		duplications = os.path.join(config["output"]["outdir"], "{sample}", config["output"]["soapnuke"], "dupReads.*.gz")
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
		-o {params.rg} 
		seqkit stats -j 2 -a  {output.soapnuke_fq1} > {output.soapnuke_fq1_stat}
		seqkit stats -j 2 -a  {output.soapnuke_fq2} > {output.soapnuke_fq1_stat}
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
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["prinseq"],"{sample}.logs")
	benchmark:
		os.path.join(config["output"]["relative"],"{sample}/benchmark/QC/{sample}.prinseq.benchmark.txt")
	params:
		prinseq = config["softwares"]["prinseq"],
		workdir = os.path.join(config["output"]["relative"],"{sample}",config["output"]["prinseq"],"{sample}")
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
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"rmrRNA.log")
	benchmark:
		os.path.join(config["output"]["relative"],"{sample}/benchmark/QC/{sample}.rmrRNA1.benchmark.txt")
	params:
		sortmerna_V4 = config["softwares"]["sortmerna_V4"],
		rRNA = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rRNA"),
		rmrRNA = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rmrRNA"),
		workdir = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"]),
		kvdb_dir = os.path.join(config["output"]["outdir"], "{sample}", config["output"]["rmrRNA"], "kvdb"),
		readb_dir = os.path.join(config["output"]["outdir"], "{sample}", config["output"]["rmrRNA], "readb"),
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
		--workdir {params.workdir} --idx-dir {params.indexdir} \
		--fastx --aligned {params.rRNA} --other {params.rmrRNA}  \
		--no-best --num_alignments 1 --paired_out  --out2 --zip-out --threads 4
		seqkit stats -j 2 -a  {input.prinseq_fq1} > {output.sortmerna_fq1_stat}
		seqkit stats -j 2 -a  {input.prinseq_fq2} > {output.sortmerna_fq2_stat}
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


#################################################kraken&&暂时不做###################################################################
rule kraken2:
	input:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_1_rmRNA.fastq"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_2_rmRNA.fastq")
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["kraken2"],"{sample}.kraken2.report"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["kraken2"],"{sample}.kraken2")
	log:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["kraken2"],"kraken2.log")
	benchmark:
		os.path.join(config["output"]["relative"],"{sample}/benchmark/kraken/{sample}.kraken.benchmark.txt")
	shell:
		"""
		kraken2 --db /hwfssz5/ST_INFECTION/GlobalDatabase/share/Database/K2DB/ncbik2db \
		--threads 5 --report {output[0]}  \
		--paired {input[0]} {input[1]} > {output[1]}
		"""


rule bracken:
	input:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["kraken2"],"{sample}.kraken2.report")
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["bracken"],"{sample}.Family.bracken"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["bracken"],"{sample}.Genus.bracken"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["bracken"],"{sample}.Species.bracken"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["bracken"],"{sample}.Family.bracken.report"),
                os.path.join(config["output"]["relative"],"{sample}",config["output"]["bracken"],"{sample}.Genus.bracken.report"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["bracken"],"{sample}.Species.bracken.report")
	log:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["bracken"],"{sample}.log")
	benchmark:
		os.path.join(config["output"]["relative"],"{sample}/benchmark/bracken/{sample}.bracken.benchmark.txt")
	shell:
		'''
		bracken -d /hwfssz5/ST_INFECTION/GlobalDatabase/share/Database/K2DB/ncbik2db -i {input} -o {output[0]} -w {output[3]} -l F
		bracken -d /hwfssz5/ST_INFECTION/GlobalDatabase/share/Database/K2DB/ncbik2db -i {input} -o {output[1]} -w {output[4]} -l G
		bracken -d /hwfssz5/ST_INFECTION/GlobalDatabase/share/Database/K2DB/ncbik2db -i {input} -o {output[2]} -w {output[5]} -l S 
		'''

##############################################02.assemble-megahit/trinity/metaspades##############################################
rule Assembly_megahit:
	input:
		sortmerna_fq1 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rmrRNA_fwd.fq.gz"),
		sortmerna_fq2 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rmrRNA_rev.fq.gz"),
		done = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"rmrRNA_done.txt")
	output:
		fasta = os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit"],"final.contigs.fa"),
		fasta_stat = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"megahit_seqkit_stats.txt"),
		fasta_fx2tab = os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"megahit_seqkit_fx2tab.txt")
	benchmark:
		os.path.join(config["output"]["relative"],"{sample}/benchmark/assemble/{sample}.megahit.benchmark.txt")
	params:
		megahit = config["softwares"]["megahit"],
		workdir = os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit"]),
		interdir = os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit"],"intermediate_contigs")
	shell:
		'''
		{params.megahit} \
		-1 {input.sortmerna_fq1} -2 {input.sortmerna_fq2} -o {params.workdir} -f
		seqkit stats -j 2 -a  {output.fasta} > {output.fasta_stat}
		seqkit fx2tab -j 2 -l -n -i {output.fasta} > {output.fasta_fx2tab}
		rm -rf {params.interdir}
		'''




###############################################鉴定宿主来源####################################################
rule cox1_findhost:
	input:
		sortmerna_fq1 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rmrRNA_fwd.fq.gz"),
		sortmerna_fq2 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rmrRNA_rev.fq.gz"),
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["cox1_findhost"],"{sample}.centrifuge.KSout.Ssum")
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
		#在脊椎动物
		bash /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/cox1_searchbylineages.sh -t Vertebrata -s {params.luname2} -l S -i {params.cox1} -o {params.cox1}
		#得到丰度比较高的，在种级别
		bash /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/cox1_summary.sh -s {params.luname2} -l S -i {params.cox1} -o {params.cox1}
		'''


###################################################机器学习识别病毒########################################
rule kmers_virfinder:
	input:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit"],"final.contigs.fa")
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit"],"final.virfinder")
	params:
		l1 = config["output"]["relative"],
		l2 = "{sample}"
	shell:
		'''
		sh /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/virfinder_results_new_1.sh {output[0]}
		'''


################################################初筛，优化资源消耗######################################
rule Pre_cdhitest:
	input:
		megahit_fa = os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit"],"final.contigs.fa")
	output:
		cdhitest = temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["cdhitest"],"{sample}.unique.contigs.fasta"))
	params:
		l2="{sample}"
	log:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["cdhitest"],"{sample}.log")
	shell:
		'''
		/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/envs.multi-user/cd-hit-4.8.1/bin/cd-hit-est \
		-i {input.megahit_fa} -o {output.cdhitest} -M 10240 -T 5 -c 0.95 -aS 0.9 -g 1  2> {log}
		'''


rule Pre_Blastx_VirusProtein:
	input:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["cdhitest"],"{sample}.unique.contigs.fasta")
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["blastx"],"{sample}.blast"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["blastx"],"{sample}.id.txt"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["blastx"],"{sample}.blast.fasta")
	log:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["blastx"],"blastx.log")
	shell:
		'''
		/hwfssz5/ST_INFECTION/GlobalDatabase/user/fengqikai/software/diamond blastx \
		-d /hwfssz5/ST_INFECTION/GlobalDatabase/share/public_database/tblastx_screen_database/Refseq_nr_IMGVR.virus.faa \
		-q {input[0]} -k 1 -f 6 -e 0.00001 -o {output[0]}  2> {log}
		
		cat {output[0]} | awk '{{print $1}}' > {output[1]}
		
		seqkit grep -f {output[1]} {input[0]} > {output[2]}
		'''



##############################################contigs注释##############################################
rule Annotation_CAT_RefseqNR:
	input:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["blastx"],"{sample}.blast.fasta")
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.contig2classification.txt"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.ORF2LCA.txt"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.predicted_proteins.faa"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.contig.tax"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.contig.official.tax"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.ORF.tax"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.summary.txt"),
		temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.alignment.diamond"))
	benchmark:
		os.path.join(config["output"]["relative"],"{sample}/benchmark/down/{sample}.CAT_megahit1.benchmark.txt")
	params:
		CAT = config["softwares"]["CAT"],
		prodigal = config["softwares"]["prodigal"],
		diamond = config["softwares"]["diamond"],
		workdir = os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}"),
		db = config["database"]["cat_db"],
	log:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"CAT.log")
	threads:4
	shell:
		'''
		{params.CAT} contigs \
		-n 4 --force --block_size 5.0 --index_chunks 1 -c {input} \
		-o {params.workdir} --no_stars \
		-d {params.db} \
		-t /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_tax/ \
		--path_to_prodigal {params.prodigal} \
		--path_to_diamond {params.prodigal}

		{params.CAT} add_names -i {output[0]} -o {output[4]} \
		-t /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_tax/ --only_official
		{params.CAT} add_names -i {output[0]} -o {output[3]} \
		-t /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_tax/
		{params.CAT} add_names -i {output[1]} -o {output[5]} \
		-t /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_tax/ --only_official

		/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/CAT/CAT_pack/CAT summarise \
		-c {input[0]} -i {output[4]} -o {output[6]} 
		'''


rule Annotation_get_RNAVirus:
	input:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["blastx"],"{sample}.blast.fasta"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.contig.tax"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.blast")
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatRNARdRpContigs.fasta"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatRNARdRpContigs.id"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatRNAVirusContigs.fasta"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatRNAVirusContigs.id"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatContigsProtein.blastp2rdrp.fasta"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatContigsProtein.blastp2rdrp.id"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.virus.contig.tax"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.virus.contig.info"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.CAT_virus_failed"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.CAT_RNA_virus_failed")
	shell:
		'''
		seqkit=/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/bin/seqkit
		#获取CAT的信息
		echo start CAT_extract_visual at:`date`
		grep 'Viruses (superkingdom)' {input[1]} |awk -v OFS='\t' -v FS='\t' '{{print $0}}'|sort|uniq > {output[6]}

		touch {output[6]}

		pid_1=$(cat {output[6]}|wc -l)
		if [ $pid_1 -eq 0 ];
		then
			echo -e "virus" >> {output[8]}
		else
			python3 /hwfssz5/ST_INFECTION/GlobalDatabase/share/wangyifei_handover/genome_type/CAT_info_simplify.py \
			/hwfssz5/ST_INFECTION/GlobalDatabase/share/wangyifei_handover/genome_type/ICTV_genometype.index {output[6]} {output[7]}
		fi


		echo end CAT_extract_visual at:`date`
		#获取比对到RdRp的contigs
		echo start RNA_RdRp.fasta at:`date`
		awk '{{print $1}}' {input[2]}|awk -F _ '{{print $1"_"$2}}'|sort|uniq > {output[5]}
		$seqkit grep -f {output[5]} {input[0]} > {output[4]}

		#获取CAT注释为RNA病毒的contigs
		awk '$4 ~ /RNA/ {{print $1}}' {output[7]}|sort|uniq > {output[3]}
		touch {output[3]}

		pid_2=$(cat {output[3]}|wc -l)
		if [ $pid_2 -eq 0 ];
		then
			echo -e "RNA_virus" >> {output[9]}
		else
			$seqkit grep -f {output[3]} {input[0]} > {output[2]}
		fi

		#对CAT注释为RNA病毒的contigs和比对到RdRp的contigs取交集
		sort {output[5]} {output[3]} |uniq -d > {output[1]}
		#$seqkit grep -f {output[1]} {input[0]} > {output[0]}

		touch {output[2]}
		touch {output[0]}
		touch {output[1]}
		touch {output[5]}
		touch {output[4]}
		touch {output[8]}
		touch {output[9]}
		echo end RNA_RdRp.fasta at:`date`
		'''


rule RNA_Virus_Checkv:
	input:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatRNAVirusContigs.fasta")
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_CatRNAVirus_checkv"],"quality_summary.tsv")  
	params:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_CatRNAVirus_checkv"])
	shell:
		'''
		/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/envs.multi-user/checkv/bin/checkv end_to_end {input[0]} {params[0]} \
		-d /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/envs.multi-user/checkv/checkv-db-v0.6
		'''


########################################rdrp比对###########################################################################


rule RNAVirus_Blastp_VirusRdRp:
	input:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.predicted_proteins.faa")
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.blast")
	benchmark:
		os.path.join(config["output"]["relative"],"{sample}/benchmark/down/{sample}.rdrp_blastp.benchmark.txt")
	shell:
		'''
		/hwfssz5/ST_INFECTION/GlobalDatabase/user/fengqikai/software/diamond blastp \
		-d /hwfssz5/ST_INFECTION/GlobalDatabase/share/public_database/rdrp_AAseq/rdrp.mBioANDrefseq \
		-q {input} -k 1 -f 6 -e 0.00001 -o {output[0]}


		'''


######################################################NT去假阳#################################################
rule RNAVirus_Blastn_nt:
	input:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatRNAVirusContigs.fasta")
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatRNAVirus_nt.blastn"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatRNAVirus_nt.blastn2")
	benchmark:
		os.path.join(config["output"]["relative"],"{sample}/benchmark/down/{sample}.blastn_nt.benchmark.txt")
	shell:
		'''
		/hwfssz5/ST_INFECTION/GlobalDatabase/user/liqian6/tools/ncbi-blast-2.10.1+/bin/blastn \
		-db /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/liqian6/blastdb_20201224/nt  \
		-query {input[0]}  -out {output[0]} -num_threads 32 \
		-outfmt "6 qacc sacc qlen slen pident evalue bitscore mismatch staxids sscinames scomnames sblastnames sskingdoms" \
		-max_target_seqs 1

		export BLASTDB={params.BLASTDB_NT}
		{params.blastn} -num_threads {params.threads} -query {input[0]} \
		-db {params.NT_db} -out {output[0]} -evalue {params.evalue} \
		-max_target_seqs {params.max_target_seqs} -outfmt {params.outfmt}
		'''


#############################################近缘物种共线性分析#####################################################################
rule Pre_Tblastx_Virusgenome:
	input:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["blastx"],"{sample}.blast.fasta")
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["closespecies"],"{sample}.contigspecies_info")
	params:
		lupath = config["output"]["relative"],
		luname = "{sample}",
		closes =  os.path.join(config["output"]["relative"],"{sample}",config["output"]["closespecies"])
	shell:
		'''
		sh /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/get_close.species_news.sh {params.lupath}  {params.luname} {params.closes}
		'''


##################################################建树##########################################################
rule RNAVirus_iqtree:
	input:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.blast"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["CAT_megahit"],"{sample}.predicted_proteins.faa")
	output:
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_rdrp"],"{sample}.finish.combine"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.blast.filt")
	params:
		rdrp = os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_rdrp"]),
		tree = os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_tree"]),
		visual = os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_visual"]),
		l2="{sample}"
	shell:
		'''
		perl /hwfssz5/ST_INFECTION/GlobalDatabase/user/zhaohailong/02.diamond.blast.filt.snakemake.pl \
		/hwfssz5/ST_INFECTION/GlobalDatabase/share/public_database/rdrp_AAseq/rdrp.mBioANDrefseq.faa.len {input[0]}
		touch {output[1]}
		perl /hwfssz5/ST_INFECTION/GlobalDatabase/user/zhaohailong/03.rdrp.family.snakemake.pl {output[1]} {input[1]} {params.rdrp}	
		touch {output[0]}
		'''

#ls {params.rdrp}/*.fa | perl -ne 'chomp;$_=~/\/(.*)\.(.*)\.fa/;`cat $_ /hwfssz5/ST_INFECTION/GlobalDatabase/share/public_database/rdrp_AAseq/rdrp.mBioANDrefseq.groupByFamily/$2.seqID.fa >$_.combine`;'
########################################mapping back&&&coverage################################################################
rule RNAVirus_Bowtie2_coverage:
	input:
		fasta = os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}.CatRNAVirusContigs.fasta"),
		sortmerna_fq1 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rmrRNA_fwd.fq.gz"),
		sortmerna_fq2 = os.path.join(config["output"]["relative"],"{sample}",config["output"]["rmrRNA"],"{sample}_rmrRNA_rev.fq.gz"),
	output:
		temp(os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}_virus.sam")),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"{sample}_virus_sofa.txt"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"{sample}_virus_sorted.bam"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"{sample}_virus_bam.stat"),
		os.path.join(config["output"]["relative"],"{sample}",config["output"]["stats"],"{sample}_virus_sort.bam")
	benchmark:
		os.path.join(config["output"]["relative"],"{sample}/benchmark/mapping/{sample}.bowtie2.benchmark.txt")
	params:
		samtools = config["softwares"]["samtools"],
		index = os.path.join(config["output"]["relative"],"{sample}",config["output"]["megahit_blastp"],"{sample}_spades")
	shell:
		'''
		bowtie2-build {input.fasta}  {params.index}
		bowtie2 -p 8 -x {params.index} -1 {input.sortmerna_fq1} -2 {input.sortmerna_fq2} -S {output[0]} 
		/ldfssz1/ST_INFECTION/P18Z10200N0164_Resistance_PN/User/zhaohailong/software/soap.coverage -cvg -sam -refsingle {input.fasta} -i {output[0]} -o {output[1]}
		{params.samtools} view -@ 3 -bS {output[0]} -F 4 | $samtools sort -@ 3 -n -o {output[2]}
		{params.samtools} sort -o {output[4]} {output[2]}
		{params.samtools} stats -@ 3 -c 10,100000,1 -i 100000 {output[2]} > {output[3]}
		'''

