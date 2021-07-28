usr=$1
list_dir=/jdfssz1/ST_HEALTH/P20Z10200N0206/fengqikai/zhongshan/scripts/zhongda.smk/Latest/Report/list
sdirs=$list_dir/sample.dir.list
blastp2RdRp_hits_dir_list=$list_dir/blastp2RdRp_hits_dir.list
contigs_info_stat=$list_dir/../contigsInfo/${usr}.contigs_info.stat.txt
echo -e "sample_id\tdir_name\tcontig_id\tcontig_length\tvirus_type\tcontig_Order\tcontig_Family\tcontig_Genus\tcontig_species\tcontig_rdrp_acc\trdrp_length\tcontig_rdrp_identity\tcontig_rdrp_Evalue\tcontig_rdrp_bitscore\trdrp_blastp_family" > $contigs_info_stat

if [ -f $blastp2RdRp_hits_dir_list ];then
	rm -f $blastp2RdRp_hits_dir_list
fi

cat $sdirs|while read dir
do
	sample=`basename $dir`
	blastp2RdRp_dir=$dir/07.ViralContigs_Annotation/CAT/blastp2RdRp
	blastp2RdRp_out=$blastp2RdRp_dir/${sample}.CAT_prodigal_proteins.blastp2RdRp.blastp
	if [ -s $blastp2RdRp_out ];then
		echo $dir >> $blastp2RdRp_hits_dir_list
	fi
done

if [ $usr == 'fengqikai' ];then
	grep fengqikai $blastp2RdRp_hits_dir_list > $list_dir/fengqikai.blastp2RdRp_hits_dir.list
	blastp2RdRp_hits_dir_use_list=$list_dir/fengqikai.blastp2RdRp_hits_dir.list
elif [ $usr == 'lijiandong' ];then
	grep lijiandong $blastp2RdRp_hits_dir_list > $list_dir/lijiandong.blastp2RdRp_hits_dir.list
	blastp2RdRp_hits_dir_use_list=$list_dir/lijiandong.blastp2RdRp_hits_dir.list
elif [ $usr == 'wangrong' ];then
	grep wangrong $blastp2RdRp_hits_dir_list > $list_dir/wangrong.blastp2RdRp_hits_dir.list
	blastp2RdRp_hits_dir_use_list=$list_dir/wangrong.blastp2RdRp_hits_dir.list
elif [ $usr == 'xuzheng' ];then
	grep xuzheng $blastp2RdRp_hits_dir_list > $list_dir/xuzheng.blastp2RdRp_hits_dir.list
	blastp2RdRp_hits_dir_use_list=$list_dir/xuzheng.blastp2RdRp_hits_dir.list
elif [ $usr == 'zhaohailong' ];then
	grep zhaohailong $blastp2RdRp_hits_dir_list > $list_dir/zhaohailong.blastp2RdRp_hits_dir.list
	blastp2RdRp_hits_dir_use_list=$list_dir/zhaohailong.blastp2RdRp_hits_dir.list
else
	echo invalid usr name: $usr !
fi

cat $blastp2RdRp_hits_dir_use_list|while read dir
do
	sample=`basename $dir`
	#contigs长度
	contigs_length=$dir/04.Assembly/${sample}.contigs.length
	#选出CAT注释出的pure viral contigs中潜在的非逆转录的RNA病毒contigs
	CAT_dir=$dir/07.ViralContigs_Annotation/CAT
	CAT_viralContigs_info=$CAT_dir/${sample}.AllVirus.info
	nonRT_RNAViral_contigs_info=$CAT_dir/${sample}.potential.nonRTRNAVirus.info
	awk '($4 ~ /RNA/ && $4 !~ /RT/) || ($4 ~ /Nan/ && $5 ~ /Riboviria/) {print $0}' $CAT_viralContigs_info > $nonRT_RNAViral_contigs_info
	#CheckV结果
	ChechkV_dir=$dir/06.Check_ViralContig_Completeness/CheckV/quality_summary.tsv
	CheckV_out=$ChechkV_dir/quality_summary.tsv
	#比对RdRp的结果取besthit
	blastp2RdRp_dir=$dir/07.ViralContigs_Annotation/CAT/blastp2RdRp
	blastp2RdRp_out=$blastp2RdRp_dir/${sample}.CAT_prodigal_proteins.blastp2RdRp.blastp
	blastp2RdRp_id=$blastp2RdRp_dir/${sample}.CAT_prodigal_proteins.blastp2RdRp.id
	awk '{print $1}' $blastp2RdRp_out|sort|uniq > $blastp2RdRp_id
	blastp2RdRp_out_besthit=$blastp2RdRp_dir/${sample}.CAT_prodigal_proteins.blastp2RdRp.besthit.blastp
	echo -e "qseqid\tqlen\tsseqid\tslen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tlength\tpident\tmismatch\tgaps\tstitle\tqcovhsp\tscovhsp" > $blastp2RdRp_out_besthit
	cat $blastp2RdRp_id|while read id
	do
		grep -w $id $blastp2RdRp_out|head -1 >> $blastp2RdRp_out_besthit
	done

	##输出CAT注释的潜在的非逆转录RNA病毒contigs的信息
	cat $nonRT_RNAViral_contigs_info|while read Contig Virus_gene_count All_gene_count Type Realm Kingdom Phylum Class Order Family Genus Species
	do
		sample_id=`echo $sample|awk -F '_' '{print $1}'`
		contig_id=$Contig
		virus_type=$Type
		contig_Order=$Order
		contig_Family=$Family 
		contig_Genus=$Genus
		contig_species=$Species
		contig_length=`grep -w $Contig $contigs_length|awk '{print $2}'`

		checkv_symbol=`grep -w $Contig $CheckV_out|wc -l`
		if [ $checkv_symbol -gt 0 ];then #如果CheckV的结果中包含该contigs
			checkv_quality=`grep -w $Contig $CheckV_out|awk '{print $8}'`
		else #如果CheckV的结果中不包含该contigs
			checkv_quality='NA'
		fi	

		blastp2rdrp_symbol=`grep ${Contig}_* $blastp2RdRp_out_besthit|wc -l`
		if [ $blastp2rdrp_symbol -gt 0 ];then
			contig_rdrp_identity=`grep ${Contig}_* $blastp2RdRp_out_besthit|head -1|awk -v FS='\t' '{print $12}'`
			rdrp_length=`grep ${Contig}_* $blastp2RdRp_out_besthit|head -1|awk -v FS='\t' '{print $4}'`
			contig_rdrp_Evalue=`grep ${Contig}_* $blastp2RdRp_out_besthit|head -1|awk -v FS='\t' '{print $9}'`
			contig_rdrp_bitscore=`grep ${Contig}_* $blastp2RdRp_out_besthit|head -1|awk -v FS='\t' '{print $10}'`
			contig_rdrp_acc=`grep ${Contig}_* $blastp2RdRp_out_besthit|head -1|awk -v FS='\t' '{print $3}'`
			rdrp_blastp_family=`grep ${Contig}_* $blastp2RdRp_out_besthit|head -1|awk -v FS='\t' '{print $15}'|awk -v FS='|' '{print $NF}'`
		else
			contig_rdrp_identity='NA'
			rdrp_length='NA'
			contig_rdrp_Evalue='NA'
			contig_rdrp_bitscore='NA'
			contig_rdrp_acc='NA'
			rdrp_blastp_family='NA'
		fi
		##输出
		echo -e "$sample_id\t$sample\t$contig_id\t$contig_length\t$virus_type\t$contig_Order\t$contig_Family\t$contig_Genus\t$contig_species\t$contig_rdrp_acc\t$rdrp_length\t$contig_rdrp_identity\t$contig_rdrp_Evalue\t$contig_rdrp_bitscore\t$rdrp_blastp_family" >> $contigs_info_stat
	done
done
