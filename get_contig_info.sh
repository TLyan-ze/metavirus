month=`date|awk '{print $2}'`
date=`date|awk '{print $3}'`
year=`date|awk '{print $6}'`
date=${year}${month}${date}

seqkit=/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/bin/seqkit
python=/hwfssz5/ST_INFECTION/GlobalDatabase/user/fengqikai/software/.conda/envs/Trinity-2.11.0/bin/python
cat_contig_info=/ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/get_contig_info.py
Host_info=/ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/shell/scripts/VMR_180521_MSL36.xlsx
rdrp_family=/ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/shell/Lastest/Viruses.fullnamelineage.dmp
#比对结果所在的路径
#样本名对应的文件夹名
basedir=$1
samplesbatchsdate=$2
out_info=$basedir/$samplesbatchsdate/${samplesbatchsdate}_nr_contig_info.${date}
#获取CAT注释为RNA病毒的contigs的id和信息
catsmp_file=$basedir/$samplesbatchsdate/04.AssembleFindVirus_BaseBlast/Cdhiest-BlastxVirusNR-CAT/${samplesbatchsdate}.virus.contig.info
Cat_RNAVirusContigs_info=$basedir/$samplesbatchsdate/03.stats/${samplesbatchsdate}.nr_CatRNAVirusContigs_info
Cat_RNAVirusContigs_id=$basedir/$samplesbatchsdate/03.stats/${samplesbatchsdate}.nr_CatRNAVirusContigs.id
grep -v DNA $catsmp_file|awk '{print $1}' |sort|uniq > $Cat_RNAVirusContigs_id
grep -v DNA $catsmp_file|sort|uniq > $Cat_RNAVirusContigs_info
#CheckV结果
CheckV_out=$basedir/$samplesbatchsdate/05.Check_ViralContigs/Cdhiest-BlastxVirusNR-checkv/quality_summary.tsv
#contig长度
contigs_length=$basedir/$samplesbatchsdate/03.stats/megahit_seqkit_fx2tab.txt
#得到最佳的比对rdrp的几个
blastp2RdRp_out=$basedir/$samplesbatchsdate/04.AssembleFindVirus_BaseBlast/Cdhiest-BlastxVirusNR-CAT-blastp2RdRp/${samplesbatchsdate}.CAT_prodigal_proteins_blastp2RdRp
blastp2RdRp_id=$basedir/$samplesbatchsdate/04.AssembleFindVirus_BaseBlast/Cdhiest-BlastxVirusNR-CAT-blastp2RdRp/${samplesbatchsdate}.CAT_prodigal_proteins_blastp2RdRp.id
awk '{print $1}' $blastp2RdRp_out|sort|uniq > $blastp2RdRp_id
blastp2RdRp_out_besthit=$basedir/$samplesbatchsdate/04.AssembleFindVirus_BaseBlast/Cdhiest-BlastxVirusNR-CAT-blastp2RdRp/${samplesbatchsdate}.CAT_prodigal_proteins.blastp2RdRp.besthit.blastp
blastp2RdRp_out_match=$basedir/$samplesbatchsdate/04.AssembleFindVirus_BaseBlast/Cdhiest-BlastxVirusNR-CAT-blastp2RdRp/${samplesbatchsdate}.CAT_prodigal_proteins.blastp2RdRp.match.blastp
echo -e "qseqid\tqlen\tsseqid\tslen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tlength\tpident\tmismatch\tgaps\tstitle\tqcovhsp\tscovhsp" > $blastp2RdRp_out_besthit
echo -e "qseqid\tqlen\tsseqid\tslen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tlength\tpident\tmismatch\tgaps\tstitle\tqcovhsp\tscovhsp" > $blastp2RdRp_out_match
cat $blastp2RdRp_id|while read id
do
    grep -w $id $blastp2RdRp_out|sort -k10r|head -1 >> $blastp2RdRp_out_besthit
    grep -w $id $blastp2RdRp_out|grep -v unassigned |sort -k10r|head -1 >> $blastp2RdRp_out_match
done
#比对nt库的结果
blastn_info=$basedir/$samplesbatchsdate/04.AssembleFindVirus_BaseBlast/Cdhiest-BlastxVirusNR-CAT-blastn2NT/${samplesbatchsdate}.CatRNAVirus_nt.blastn
blastn_id=$basedir/$samplesbatchsdate/04.AssembleFindVirus_BaseBlast/Cdhiest-BlastxVirusNR-CAT-blastn2NT/${samplesbatchsdate}.CatRNAVirus_nt.blastn.id
awk '{print $1}' $blastn_info|sort|uniq > $blastn_id
blastn_out=$basedir/$samplesbatchsdate/04.AssembleFindVirus_BaseBlast/Cdhiest-BlastxVirusNR-CAT-blastn2NT/${samplesbatchsdate}.CatRNAVirus_nt_best.blastn
echo -e "qacc\tsacc\tqlen\tslen\tpident\tevalue\tbitscore\tmismatch\tstaxids\tsscinames\tscomnames\tsblastnames\tsskingdoms" > $blastn_out
cat $blastn_id|while read id
do
    grep -w $id $blastn_info|sort -k7r|head -1 >> $blastn_out
done
#覆盖度结果
sofa_out=$basedir/$samplesbatchsdate/05.Check_ViralContigs/bowtie2_contig_quality/${samplesbatchsdate}_virus_sofa.txt
Cat_sofa_info=$basedir/$samplesbatchsdate/05.Check_ViralContigs/bowtie2_contig_quality/${samplesbatchsdate}_virus_coverage.txt
awk '$1 ~ /k/ {print $1,$3,$4}' $sofa_out|awk -F : '{print $1,$3,$4}'|awk '{print $1"\t"$2"\t"$4}'> $Cat_sofa_info
#共线性的结果
gong_info=$basedir/$samplesbatchsdate/04.AssembleFindVirus_BaseBlast/Cdhiest-BlastxVirusNR-tblastxVirusGenome/${samplesbatchsdate}.contigspecies_info
gong_coverage=$basedir/$samplesbatchsdate/04.AssembleFindVirus_BaseBlast/Cdhiest-BlastxVirusNR-tblastxVirusGenome/${samplesbatchsdate}.gong_coverage_txt
#psi-blast的结果
psi_blast=$basedir/$samplesbatchsdate/04.AssembleFindVirus_BaseBlast/Cdhiest-BlastxVirusNR-CAT/${samplesbatchsdate}.predicted_proteins.faa.cdds.psiblast
rps_blast=$basedir/$samplesbatchsdate/04.AssembleFindVirus_BaseBlast/Cdhiest-BlastxVirusNR-CAT/${samplesbatchsdate}.predicted_proteins.faa.cdds.rpsblast
psi_rps_blast=$basedir/$samplesbatchsdate/04.AssembleFindVirus_BaseBlast/Cdhiest-BlastxVirusNR-CAT/${samplesbatchsdate}.predicted_proteins.faa.cdds.he
cat $psi_blas > $psi_rps_blast
cat $rps_blast >> $psi_rps_blast
psi_blast_id=$basedir/$samplesbatchsdate/04.AssembleFindVirus_BaseBlast/Cdhiest-BlastxVirusNR-CAT/${samplesbatchsdate}.predicted_proteins.faa.cdds.psiblast.id
awk '{print $1}' $psi_rps_blast|sort|uniq > $psi_blast_id
psiblast_out_besthit=$basedir/$samplesbatchsdate/04.AssembleFindVirus_BaseBlast/Cdhiest-BlastxVirusNR-CAT/${samplesbatchsdate}.predicted_proteins.faa.cdds.best.psiblast
echo -e "qacc\tsacc\tqlen\tslen\tpident\tevalue\tbitscore\tmismatch\tstaxids\tsscinames\tstitle" > $psiblast_out_besthit
cat $psi_blast_id|while read id
do
    grep -w $id $psi_rps_blast|sort -k7r|head -1 >> $psiblast_out_besthit
done


#整合CAT，长度和checkv和RdRp和blastn的结果
#添加头文件
echo -e "sequencing_sample_id,sample_id,contig_id,contig_length,virus_type,contig_Order,contig_Family,contig_Genus,contig_species,checkv_quality,contig_rdrp_acc,rdrp_length,contig_rdrp_identity,contig_rdrp_Evalue,contig_rdrp_bitscore,rdrp_blastp_family,rdrp_match_family,psi_sacc,psi_cds,contig_psi_identity,contig_psi_Evalue,contig_psi_bitscore,psi_stitle,contig_nt_indentity,nt_length,contig_nt_bitscore,contig_nt_acc,contig_nt_Subject_Scientific_Name,contig_nt_Super_Kingdom,reads_to_contigs_coverage,reads_to_contigs_depth,tblastx_RNA_database_with_RTRNA,tblastx_RNA_database_with_RTRNA_length,match_virus_coverage_to_contigs" > $out_info
cat $Cat_RNAVirusContigs_info|while read Contig Virus_gene_count All_gene_count Type Realm Kingdom Phylum Class Order Family Genus Species
do
    sample_id=`echo $samplesbatchsdate|awk -F '_' '{print $1}'`
    contig_id=$Contig
    virus_type=$Type
    contig_Order=$Order
    contig_Family=$Family 
    contig_Genus=$Genus
    contig_species=$Species
    contig_length=`grep -w $Contig $contigs_length|awk '{print $2}'`

    checkv_symbol=`grep -w $Contig $CheckV_out|wc -l`
    if [ $checkv_symbol -gt 0 ];then #如果CheckV的结果中包含该contigs
        checkv_quality=`grep -w $Contig $CheckV_out|awk -v FS='\t' '{print $8}'`
    else #如果CheckV的结果中不包含该contigs
        checkv_quality='NO'
    fi	

    blastp2rdrp_symbol=`grep ${Contig}_* $blastp2RdRp_out_besthit|wc -l`
    if [ $blastp2rdrp_symbol -gt 0 ];then
        contig_rdrp_identity=`grep ${Contig}_* $blastp2RdRp_out_besthit|head -1|awk -v FS='\t' '{print $12}'`
        rdrp_length=`grep ${Contig}_* $blastp2RdRp_out_besthit|head -1|awk -v FS='\t' '{print $4}'`
        contig_rdrp_Evalue=`grep ${Contig}_* $blastp2RdRp_out_besthit|head -1|awk -v FS='\t' '{print $9}'`
        contig_rdrp_bitscore=`grep ${Contig}_* $blastp2RdRp_out_besthit|head -1|awk -v FS='\t' '{print $10}'`
        contig_rdrp_acc=`grep ${Contig}_* $blastp2RdRp_out_besthit|head -1|awk -v FS='\t' '{print $3}'`
        rdrp_blastp_family=`grep ${Contig}_* $blastp2RdRp_out_besthit|head -1|awk -v FS='\t' '{print $15}'|awk -v FS='|' '{print $NF}'`
        blastp_match_symbol=`grep ${Contig}_* $blastp2RdRp_out_match|wc -l`
        if [ $blastp_match_symbol -gt 0 ];then
            rdrp_match_family=`grep ${Contig}_* $blastp2RdRp_out_match|head -1|awk -v FS='\t' '{print $15}'|awk -v FS='|' '{print $NF}'`
        else
            rdrp_match_family='NO'
        fi
        #rdrp_id=`grep ${Contig}_* $blastp2RdRp_out_besthit|head -1|awk -v FS='\t' '{print $15}'|awk -v FS='|' '{print $2}'`
        #rdrp_blastp_order=`grep -w $rdrp_id $rdrp_family |awk -v FS='|' '{print $3}'|awk -v FS=';' '{print $6}' `
    else
        contig_rdrp_identity='NO'
        rdrp_length='NO'
        contig_rdrp_Evalue='NO'
        contig_rdrp_bitscore='NO'
        contig_rdrp_acc='NO'
        rdrp_blastp_family='NO'
        rdrp_match_family='NO'
    fi

    psiblastp_symbol=`grep ${Contig}_* $psiblast_out_besthit|wc -l`
    if [ $psiblastp_symbol -gt 0 ];then
        contig_psi_identity=`grep ${Contig}_* $psiblast_out_besthit|head -1|awk -v FS='\t' '{print $5}'`
        psi_sacc=`grep ${Contig}_* $psiblast_out_besthit|head -1|awk -v FS='\t' '{print $2}'`
        contig_psi_Evalue=`grep ${Contig}_* $psiblast_out_besthit|head -1|awk -v FS='\t' '{print $6}'`
        contig_psi_bitscore=`grep ${Contig}_* $psiblast_out_besthit|head -1|awk -v FS='\t' '{print $7}'`
        contig_psi_sscinames=`grep ${Contig}_* $psiblast_out_besthit|head -1|awk -v FS='\t' '{print $10}'`
        psi_cds=`grep ${Contig}_* $psiblast_out_besthit|head -1|awk -v FS='\t' '{print $11}'|awk -v FS=',' '{print $1}'`
        psi_stitle=`grep ${Contig}_* $psiblast_out_besthit|head -1|awk -v FS='\t' '{print $11}'|awk -v FS=',' '{print $2}'`
    else
        contig_psi_identity='NO'
        psi_sacc='NO'
        contig_psi_Evalue='NO'
        contig_psi_bitscore='NO'
        contig_psi_sscinames='NO'
        psi_cds='NO'
        psi_stitle='NO'
    fi

    blastn2nt_symbol=`grep ${Contig}_* $blastn_out|wc -l`
    if [ $blastn2nt_symbol -gt 0 ];then #如果blatn的结果中包含该contigs
        contig_nt_indentity=`grep -w $Contig $blastn_out|awk -v FS='\t' '{print $5}'`
        nt_length=`grep -w $Contig $blastn_out|awk -v FS='\t' '{print $4}'`
        contig_nt_bitscore=`grep -w $Contig $blastn_out|awk -v FS='\t' '{print $7}'`
        contig_nt_acc=`grep -w $Contig $blastn_out|awk -v FS='\t' '{print $2}'`
        contig_nt_Subject_Scientific_Name=`grep -w $Contig $blastn_out|awk -v FS='\t' '{print $10}'`
        contig_nt_Super_Kingdom=`grep -w $Contig $blastn_out|awk -v FS='\t' '{print $13}'`
    else #如果blastn的结果中不包含该contigs
        contig_nt_indentity='NO'
        nt_length='NO'
        contig_nt_bitscore='NO'
        contig_nt_acc='NO'
        contig_nt_Subject_Scientific_Name='NO'
        contig_nt_Super_Kingdom='NO'
    fi	

    soft_symbol=`grep -w $Contig $Cat_sofa_info|wc -l`
    if [ $soft_symbol -gt 0 ];then #如果soft.coverage的结果中包含该contigs
        reads_to_contigs_coverage=`grep -w $Contig $Cat_sofa_info|awk -v FS='\t' '{print $2}'`
        reads_to_contigs_depth=`grep -w $Contig $Cat_sofa_info|awk -v FS='\t' '{print $3}'`
    else #如果soft.coverage的结果中不包含该contigs
        reads_to_contigs_coverage='NO'
        reads_to_contigs_depth='NO'
    fi	
    gong_info_symbol=`grep -w $Contig $gong_info|wc -l`
    if [ $gong_info_symbol -gt 0 ];then #如果结果中包含该contigs
        tblastx_RNA_database_with_RTRNA=`grep -w $Contig $gong_info|awk -v FS='\t'  '{print $2}'`
        tblastx_RNA_database_with_RTRNA_length=`grep -w $Contig $gong_info|awk -v FS='\t' '{print $3}'`
    else #如果结果中不包含该contigs
        tblastx_RNA_database_with_RTRNA='NO'
        tblastx_RNA_database_with_RTRNA_length='NO'
    fi
    gong_coverage_symbol=`grep -w $Contig $gong_coverage|wc -l`
    if [ $gong_coverage_symbol -gt 0 ];then #如果结果中包含该contigs
        match_virus_coverage_to_contigs=`grep -w $Contig $gong_coverage|awk -v FS='\t' '{print $2}'`
    else #如果结果中不包含该contigs
        match_virus_coverage_to_contigs='NO'
    fi  	  
    ##输出
    echo -e "$samplesbatchsdate,$sample_id,$contig_id,$contig_length,$virus_type,$contig_Order,$contig_Family,$contig_Genus,$contig_species,$checkv_quality,$contig_rdrp_acc,$rdrp_length,$contig_rdrp_identity,$contig_rdrp_Evalue,$contig_rdrp_bitscore,$rdrp_blastp_family,$rdrp_match_family,$psi_sacc,$psi_cds,$contig_psi_identity,$contig_psi_Evalue,$contig_psi_bitscore,$psi_stitle,$contig_nt_indentity,$nt_length,$contig_nt_bitscore,$contig_nt_acc,$contig_nt_Subject_Scientific_Name,$contig_nt_Super_Kingdom,$reads_to_contigs_coverage,$reads_to_contigs_depth,$tblastx_RNA_database_with_RTRNA,$tblastx_RNA_database_with_RTRNA_length,$match_virus_coverage_to_contigs" >> $out_info
done


echo end $samplesbatchsdate at `date`
